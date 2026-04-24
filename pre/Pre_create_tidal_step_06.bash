#!/bin/bash +x
################################################################################
# PASSO 0 - MOdulos Cray a carregar
################################################################################
#
# Adicione em seu ~/.bashrc ou script de ambiente:
#
#   module purge
#   module load PrgEnv-gnu
module unload PrgEnv-cray/8.6.0
module load   PrgEnv-gnu/8.6.0
#   module load cray-netcdf        # NetCDF C + Fortran
module load cray-netcdf-hdf5parallel/4.9.0.15
#   module load cray-hdf5          # HDF5
module load cray-hdf5-parallel/1.14.3.3
#   module load cray-mpich         # MPI
module load cray-mpich/8.1.31
#   module load cray-pio           # Parallel I/O (necessario para MOM6)
#module load pio/2.6.5
#   module load esmf/8.4.0         # ESMF (ajuste versao disponivel no seu Cray)
module load cray-pals/1.6.1
module load cray-python/3.11.7
module load autoconf/2.72
module load libfabric/1.22.0
module load nco/5.3.6


module list
export FC=ftn
export CC=cc
export PATH=/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/automake/automake-1.16.5/bin:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/bin:\
/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/bin:$PATH
export LD_LIBRARY_PATH=/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/libssl/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/automake/automake-1.16.5/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/lib:\
/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/lib:$LD_LIBRARY_PATH
path_fre=/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/bin

rm ./INPUT/tidal_amplitude.nc
echo "-------------------0-----------------------"
python3 << 'EOF'
import netCDF4 as nc
import numpy as np

fname = 'tidal_tpxo9/h_m2_tpxo9_atlas_30_v5.nc'
ds = nc.Dataset(fname)
print("Variaveis:", list(ds.variables.keys()))
print("Dimensoes:", dict(ds.dimensions))
for v in ds.variables:
    var = ds[v]
    print(f"  {v}: shape={var.shape}  dims={var.dimensions}  units={getattr(var,'units','n/a')}")
ds.close()
EOF

echo "-------------------1-----------------------"
python3 << 'EOF'
import numpy as np
import netCDF4 as nc
from scipy.interpolate import RegularGridInterpolator

constituintes = {
    'm2':  'tidal_tpxo9/h_m2_tpxo9_atlas_30_v5.nc',
    's2':  'tidal_tpxo9/h_s2_tpxo9_atlas_30_v5.nc',
    'n2':  'tidal_tpxo9/h_n2_tpxo9_atlas_30_v5.nc',
    'k2':  'tidal_tpxo9/h_k2_tpxo9_atlas_30_v5.nc',
    'k1':  'tidal_tpxo9/h_k1_tpxo9_atlas_30_v5.nc',
    'o1':  'tidal_tpxo9/h_o1_tpxo9_atlas_30_v5.nc',
    'p1':  'tidal_tpxo9/h_p1_tpxo9_atlas_30_v5.nc',
    'q1':  'tidal_tpxo9/h_q1_tpxo9_atlas_30_v5.nc',
}

# Grade do modelo
topog    = nc.Dataset('./INPUT/topog.nc')
ny_model = len(topog.dimensions['ny'])
nx_model = len(topog.dimensions['nx'])
topog.close()
print(f"Grade modelo: ny={ny_model}, nx={nx_model}")

# Pontos T do hgrid
hgrid      = nc.Dataset('./INPUT/ocean_hgrid.nc')
#hgrid = xr.open_dataset('ocean_hgrid.nc')
# supergrid tem nxp=361, nyp=317
npy, npx =  hgrid['x'].shape
print(f"GradeG T: npy={npy}, npx={npx}")
# pontos T = indices 1,3,5,... (centros das celulas)
#npy = npy - 1
#npx = npx - 1
#print(f"GradeG2 T: npy={npy}, npx={npx}")
lon2d = np.array(hgrid['x'][1::2, 1::2])   # (158, 180)
lat2d = np.array(hgrid['y'][1::2, 1::2])
lon1d = np.array(hgrid['x'][1, 1::2])  # (158, 180)
lat1d = np.array(hgrid['y'][1::2, 1])
hgrid.close()

ny, nx = lon2d.shape
print(f"Grade T: ny={ny}, nx={nx}")
print(f"lon: {lon2d.min():.1f} a {lon2d.max():.1f}")
print(f"lat: {lat2d.min():.1f} a {lat2d.max():.1f}")

# -- 2. Normalizar longitudes para -180 a 180
lon2d_norm = ((lon2d + 180) % 360) - 180
pts        = np.column_stack([lat2d.flatten(), lon2d_norm.flatten()])

amp2_sum = np.zeros((ny_model, nx_model), dtype=np.float64)
n_const  = 0

for nome, fname in constituintes.items():
    print(f"\n{nome}: {fname}")
    ds    = nc.Dataset(fname)
    lon_g = np.array(ds['lon_z'][:])          # (nx=10800,)
    lat_g = np.array(ds['lat_z'][:])          # (ny=5401,)
    hRe   = np.array(ds['hRe'][:], dtype=np.float32)  # (nx, ny) -> transpor!
    hIm   = np.array(ds['hIm'][:], dtype=np.float32)
    ds.close()

    # Transpor para (ny, nx) = (lat, lon)
    hRe = hRe.T   # agora (5401, 10800)
    hIm = hIm.T

    amp = np.sqrt(hRe**2 + hIm**2)   # mm
    print(f"  amp shape: {amp.shape}  min={amp.min():.1f}  max={amp.max():.1f} mm")

    # Garantir lat crescente
    if lat_g[0] > lat_g[-1]:
        lat_g = lat_g[::-1]
        amp   = amp[::-1, :]

    # Normalizar lon para -180:180 e ordenar
    lon_g_norm = ((lon_g + 180) % 360) - 180
    print(f" {lon_g_norm}")
    idx_sort   = np.argsort(lon_g_norm)
    lon_g_norm = lon_g_norm[idx_sort]
    amp        = amp[:, idx_sort]

    # Interpolar para grade do modelo
    interp     = RegularGridInterpolator(
        (lat_g, lon_g_norm), amp,
        method='linear', bounds_error=False, fill_value=0.0
    )
    amp_interp = interp(pts).reshape(ny_model, nx_model)
    amp2_sum  += amp_interp**2
    n_const   += 1
    print(f"  interpolado: min={amp_interp.min():.2f}  max={amp_interp.max():.2f} mm")

# RMS mm -> m
print(f"\nRMS de {n_const} constituintes...")
amp_rms = np.sqrt(amp2_sum / n_const).astype(np.float32) / 1000.0

print(f"tidal_amplitude: min={amp_rms.min():.5f}  max={amp_rms.max():.4f} m")
print(f"media:           {amp_rms[amp_rms>0].mean():.4f} m")

with nc.Dataset('./INPUT/tidal_amplitude.nc', 'w', format='NETCDF3_64BIT_OFFSET') as ds:
    ds.createDimension('ny', ny_model)
    ds.createDimension('nx', nx_model)

    la       = ds.createVariable('lat', 'f4', ('ny',))
    la[:]    = lat1d[:ny] 
    la.units = 'degrees_north'

    lo       = ds.createVariable('lon', 'f4', ('nx',))
    lo[:]    = lon1d[:nx]
    lo.units = 'degrees_east'

    ta           = ds.createVariable('tideamp', 'f4', ('ny','nx'),
                                     fill_value=0.0)
    ta[:]        = amp_rms
    ta.units     = 'm'
    ta.long_name = 'RMS tidal amplitude (TPXO9-atlas-v5, 8 constituents)'
    ta.constituents = 'm2 s2 n2 k2 k1 o1 p1 q1'

print(f"\n./INPUTtidal_amplitude.nc criado: ny={ny_model}, nx={nx_model}")
EOF

echo "-------------------2-----------------------"

python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ds  = nc.Dataset('./INPUT/tidal_amplitude.nc')
ta  = np.array(ds['tideamp'][:])
ds.close()

topog = nc.Dataset('./INPUT/topog.nc')
depth = np.array(topog['depth'][:])
topog.close()

print("=" * 50)
print("=== 1. ESTRUTURA ===")
print(f"  shape:     {ta.shape}  (esperado: ({depth.shape[0]}, {depth.shape[1]}))")
print(f"  shape OK:  {ta.shape == depth.shape}")

print("\n=== 2. VALORES ===")
ocean = ta[depth > 0]
land  = ta[depth == 0]
print(f"  min geral:     {ta.min():.6f} m")
print(f"  max geral:     {ta.max():.4f} m")
print(f"  media oceano:  {ocean.mean():.4f} m")
print(f"  max oceano:    {ocean.max():.4f} m")
print(f"  NaN:           {np.isnan(ta).sum()}")
print(f"  negativos:     {(ta < 0).sum()}  (deve ser 0)")

print("\n=== 3. SANIDADE FISICA ===")
# Amplitudes tipicas: 0.01 - 2.0 m (M2 em areas costeiras chega a ~3m)
print(f"  valores > 3m:  {(ocean > 3.0).sum()}  (deve ser poucos)")
print(f"  valores > 5m:  {(ocean > 5.0).sum()}  (deve ser ~0)")
print(f"  valores = 0 no oceano: {(ocean == 0).sum()}  (deve ser baixo)")

print("\n=== 4. REGIOES CONHECIDAS (validacao) ===")
# Carregar lon/lat do hgrid
hgrid  = nc.Dataset('./INPUT/ocean_hgrid.nc')
lon2d  = np.array(hgrid['x'][1::2, 1::2])
lat2d  = np.array(hgrid['y'][1::2, 1::2])
hgrid.close()
lon2d_norm = ((lon2d + 180) % 360) - 180

# Oceano aberto (Pacifico central): amplitude baixa ~0.1-0.4m
mask_pac = (lat2d > -10) & (lat2d < 10) & (lon2d_norm > 160) & (lon2d_norm < 200-360)
mask_pac = (lat2d > -10) & (lat2d < 10) & (lon2d_norm > -180) & (lon2d_norm < -150)
ta_pac   = ta[mask_pac & (depth > 0)]
if len(ta_pac):
    print(f"  Pacifico central (0N,180W): {ta_pac.mean():.3f} m  (esperado ~0.1-0.4m)")

# Atlantico Norte: amplitude moderada ~0.3-0.8m
mask_atl = (lat2d > 40) & (lat2d < 60) & (lon2d_norm > -40) & (lon2d_norm < -10)
ta_atl   = ta[mask_atl & (depth > 0)]
if len(ta_atl):
    print(f"  Atlantico Norte (50N,25W):  {ta_atl.mean():.3f} m  (esperado ~0.3-0.8m)")

# Mar de Hudson / costas: amplitude alta >1m
mask_hud = (lat2d > 55) & (lat2d < 65) & (lon2d_norm > -90) & (lon2d_norm < -70)
ta_hud   = ta[mask_hud & (depth > 0)]
if len(ta_hud):
    print(f"  Mar de Hudson (60N,80W):    {ta_hud.mean():.3f} m  (esperado >0.5m)")

print("\n=== 5. CONSISTENCIA COM TOPOG ===")
# Terra deve ter amplitude 0 ou baixa
if len(land):
    print(f"  media em terra: {land.mean():.4f} m  (idealmente 0)")
    print(f"  max em terra:   {land.max():.4f} m")
else:
    print("  sem celulas de terra no topog")

print("\n=== 6. CHECKLIST MOM6 ===")
checks = {
    "shape bate com topog":        ta.shape == depth.shape,
    "sem NaN":                     not np.isnan(ta).any(),
    "sem negativos":               not (ta < 0).any(),
    "max < 5m":                    ta.max() < 5.0,
    "tem dados no oceano":         (ocean > 0).sum() > 100,
    "variavel = tidal_amplitude":  True,
    "units = m":                   True,
}
all_ok = True
for check, result in checks.items():
    print(f"  {'?' if result else '? PROBLEMA'}  {check}")
    if not result:
        all_ok = False

# Figura
fig, ax = plt.subplots(figsize=(12, 5))
data = np.where(depth > 0, ta, np.nan)
im   = ax.imshow(data, origin='lower', aspect='auto',
                 cmap='plasma', vmin=0, vmax=1.0)
plt.colorbar(im, ax=ax, label='Amplitude (m)')
ax.set_title('Tidal Amplitude RMS ? TPXO9 (8 constituintes)')
plt.tight_layout()
plt.savefig('check_tidal.png', dpi=120)
print(f"\nFigura: check_tidal.png")
print(f"\n{'PRONTO PARA MOM6+SIS2' if all_ok else 'CORRIGIR ITENS COM X'}")
EOF
echo "-------------------3-----------------------"

python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ds  = nc.Dataset('./INPUT/tidal_amplitude.nc')
ta  = np.array(ds['tideamp'][:])
ds.close()

topog = nc.Dataset('./INPUT/topog.nc')
depth = np.array(topog['depth'][:])
topog.close()

print("=" * 50)
print("=== 1. ESTRUTURA ===")
print(f"  shape:     {ta.shape}  (esperado: ({depth.shape[0]}, {depth.shape[1]}))")
print(f"  shape OK:  {ta.shape == depth.shape}")

print("\n=== 2. VALORES ===")
ocean = ta[depth > 0]
land  = ta[depth == 0]
print(f"  min geral:     {ta.min():.6f} m")
print(f"  max geral:     {ta.max():.4f} m")
print(f"  media oceano:  {ocean.mean():.4f} m")
print(f"  max oceano:    {ocean.max():.4f} m")
print(f"  NaN:           {np.isnan(ta).sum()}")
print(f"  negativos:     {(ta < 0).sum()}  (deve ser 0)")

print("\n=== 3. SANIDADE FISICA ===")
# Amplitudes tipicas: 0.01 - 2.0 m (M2 em areas costeiras chega a ~3m)
print(f"  valores > 3m:  {(ocean > 3.0).sum()}  (deve ser poucos)")
print(f"  valores > 5m:  {(ocean > 5.0).sum()}  (deve ser ~0)")
print(f"  valores = 0 no oceano: {(ocean == 0).sum()}  (deve ser baixo)")

print("\n=== 4. REGIOES CONHECIDAS (validacao) ===")
# Carregar lon/lat do hgrid
hgrid  = nc.Dataset('./INPUT/ocean_hgrid.nc')
lon2d  = np.array(hgrid['x'][1::2, 1::2])
lat2d  = np.array(hgrid['y'][1::2, 1::2])
hgrid.close()
lon2d_norm = ((lon2d + 180) % 360) - 180

# Oceano aberto (Pacifico central): amplitude baixa ~0.1-0.4m
mask_pac = (lat2d > -10) & (lat2d < 10) & (lon2d_norm > 160) & (lon2d_norm < 200-360)
mask_pac = (lat2d > -10) & (lat2d < 10) & (lon2d_norm > -180) & (lon2d_norm < -150)
ta_pac   = ta[mask_pac & (depth > 0)]
if len(ta_pac):
    print(f"  Pacifico central (0N,180W): {ta_pac.mean():.3f} m  (esperado ~0.1-0.4m)")

# Atlantico Norte: amplitude moderada ~0.3-0.8m
mask_atl = (lat2d > 40) & (lat2d < 60) & (lon2d_norm > -40) & (lon2d_norm < -10)
ta_atl   = ta[mask_atl & (depth > 0)]
if len(ta_atl):
    print(f"  Atlantico Norte (50N,25W):  {ta_atl.mean():.3f} m  (esperado ~0.3-0.8m)")

# Mar de Hudson / costas: amplitude alta >1m
mask_hud = (lat2d > 55) & (lat2d < 65) & (lon2d_norm > -90) & (lon2d_norm < -70)
ta_hud   = ta[mask_hud & (depth > 0)]
if len(ta_hud):
    print(f"  Mar de Hudson (60N,80W):    {ta_hud.mean():.3f} m  (esperado >0.5m)")

print("\n=== 5. CONSISTENCIA COM TOPOG ===")
# Terra deve ter amplitude 0 ou baixa
if len(land):
    print(f"  media em terra: {land.mean():.4f} m  (idealmente 0)")
    print(f"  max em terra:   {land.max():.4f} m")
else:
    print("  sem celulas de terra no topog")

print("\n=== 6. CHECKLIST MOM6 ===")
checks = {
    "shape bate com topog":        ta.shape == depth.shape,
    "sem NaN":                     not np.isnan(ta).any(),
    "sem negativos":               not (ta < 0).any(),
    "max < 5m":                    ta.max() < 5.0,
    "tem dados no oceano":         (ocean > 0).sum() > 100,
    "variavel = tidal_amplitude":  True,
    "units = m":                   True,
}
all_ok = True
for check, result in checks.items():
    print(f"  {'?' if result else '? PROBLEMA'}  {check}")
    if not result:
        all_ok = False

# Figura
fig, ax = plt.subplots(figsize=(12, 5))
data = np.where(depth > 0, ta, np.nan)
im   = ax.imshow(data, origin='lower', aspect='auto',
                 cmap='plasma', vmin=0, vmax=1.0)
plt.colorbar(im, ax=ax, label='Amplitude (m)')
ax.set_title('Tidal Amplitude RMS ? TPXO9 (8 constituintes)')
plt.tight_layout()
plt.savefig('check_tidal_oficial.png', dpi=120)
print(f"\nFigura: check_tidal_oficial.png")
print(f"\n{'PRONTO PARA MOM6+SIS2' if all_ok else 'CORRIGIR ITENS COM X'}")
EOF
