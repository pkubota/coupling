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
rm INPUT/seawifs_1998-2006_smoothed_2X.nc
echo "-------------------0 check input -----------------------"

python3 -c "
import h5py
with h5py.File('SEAWIFS/SEASTAR_SEAWIFS_GAC.19980101_20100131.L3b.MC.CHL.nc','r') as f:
    def show(name, obj):
        print(name, type(obj))
    f.visititems(show)
"

echo "-------------------1 check estrutura do arquivo-----------------------"
python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import h5py
import os
from scipy.ndimage import uniform_filter

# Arquivos em ordem Jan-Dez
files = [
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980101_20100131.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980201_20100228.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980301_20100331.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980401_20100430.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980501_20100531.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980601_20100630.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980701_20100731.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980801_20100831.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19970901_20100930.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19971001_20101031.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19971101_20101130.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19971201_20101231.L3b.MC.CHL.nc',
]

# Ver campos do BinIndex e BinList
print("=== Estrutura do arquivo ===")
with h5py.File(files[0], 'r') as f:
    grp = f['level-3_binned_data']
    print("BinIndex dtype:", grp['BinIndex'].dtype)
    print("BinIndex shape:", grp['BinIndex'].shape)
    print("BinList dtype:",  grp['BinList'].dtype)
    print("BinList shape:",  grp['BinList'].shape)
    print("chlor_a dtype:",  grp['chlor_a'].dtype)
    print("chlor_a shape:",  grp['chlor_a'].shape)
    print("BinIndex fields:", grp['BinIndex'].dtype.names)
    print("BinList fields:",  grp['BinList'].dtype.names)
    print("chlor_a fields:",  grp['chlor_a'].dtype.names)
    # Total de rows (NUMROWS)
    print("NUMROWS (BinIndex length):", len(grp['BinIndex']))

EOF


echo "-------------------2 processa dados brutos -----------------------"

python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import h5py
import os
from scipy.ndimage import uniform_filter

files = [
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980101_20100131.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980201_20100228.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980301_20100331.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980401_20100430.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980501_20100531.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980601_20100630.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980701_20100731.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19980801_20100831.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19970901_20100930.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19971001_20101031.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19971101_20101130.L3b.MC.CHL.nc',
    'SEAWIFS/SEASTAR_SEAWIFS_GAC.19971201_20101231.L3b.MC.CHL.nc',
]

# Grade saida = shape do topog
topog    = nc.Dataset('INPUT/topog.nc')
ny_model = len(topog.dimensions['ny'])   # 158
nx_model = len(topog.dimensions['nx'])   # 180
topog.close()
print(f"Grade modelo: {ny_model} x {nx_model}")

# Pre-calcular lookup: bin_num -> (row, col) VETORIZADO
NUMROWS = 2160
latbin  = ((np.arange(NUMROWS) + 0.5) / NUMROWS * 180.0) - 90.0
numbin  = np.maximum(1, np.round(np.cos(np.radians(latbin)) * 2 * NUMROWS).astype(np.int64))
cumbin  = np.concatenate([[0], np.cumsum(numbin)])
total_bins = int(cumbin[-1])
print(f"Total bins na grade: {total_bins:,}")

# Pre-calcular lat/lon para TODOS os bins de uma vez
print("Pre-calculando lat/lon de todos os bins...")
# row de cada bin
all_rows = np.zeros(total_bins, dtype=np.int32)
for r in range(NUMROWS):
    all_rows[cumbin[r]:cumbin[r+1]] = r

# col dentro da row
all_cols = np.arange(total_bins, dtype=np.int64)
for r in range(NUMROWS):
    all_cols[cumbin[r]:cumbin[r+1]] -= cumbin[r]

lat_all = latbin[all_rows]
lon_all = (all_cols + 0.5) / numbin[all_rows] * 360.0 - 180.0

# Indices na grade modelo para todos os bins
j_all = np.clip(np.round((lat_all + 90.0) / 180.0 * (ny_model - 1)).astype(int), 0, ny_model-1)
i_all = np.clip(np.round((lon_all + 180.0) / 360.0 * (nx_model - 1)).astype(int), 0, nx_model-1)
idx_flat = j_all * nx_model + i_all   # indice 1D na grade modelo
print(f"Lookup pre-calculado: {len(idx_flat):,} bins")

def read_l3b_fast(fname):
    """Leitura vetorizada usando lookup pre-calculado."""
    with h5py.File(fname, 'r') as f:
        grp     = f['level-3_binned_data']
        blist   = grp['BinList'][:]
        chl_raw = grp['chlor_a'][:]

    bin_num = blist['bin_num'].astype(np.int64) - 1   # zero-based
    weights = blist['weights']
    valid   = weights > 0
    chl_val = np.where(valid, chl_raw['sum'] / np.where(valid, weights, 1), np.nan)

    # Filtrar invalidos
    good    = np.isfinite(chl_val) & (chl_val > 0) & (chl_val < 100)
    bin_num = bin_num[good]
    chl_val = chl_val[good]

    print(f"  Bins validos: {good.sum():,} / {len(good):,}")

    # Mapear bin_num -> indice na grade modelo
    # Verificar range
    in_range = (bin_num >= 0) & (bin_num < total_bins)
    bin_num  = bin_num[in_range]
    chl_val  = chl_val[in_range]

    grid_idx = idx_flat[bin_num]   # indice 1D na grade modelo

    # Acumular na grade (media usando np.bincount)
    chl_sum  = np.bincount(grid_idx, weights=chl_val,  minlength=ny_model*nx_model)
    chl_cnt  = np.bincount(grid_idx, weights=np.ones(len(chl_val)), minlength=ny_model*nx_model)

    chl_grid = np.where(chl_cnt > 0, chl_sum / chl_cnt, 0.0)
    chl_grid = chl_grid.reshape(ny_model, nx_model).astype(np.float32)

    n_filled = (chl_grid > 0).sum()
    print(f"  Celulas preenchidas: {n_filled:,} / {ny_model*nx_model:,} "
          f"({100*n_filled/(ny_model*nx_model):.1f}%)")
    return chl_grid

# Processar 12 meses
chl_smooth = np.zeros((12, ny_model, nx_model), dtype=np.float32)
meses = ['Jan','Fev','Mar','Abr','Mai','Jun',
         'Jul','Ago','Set','Out','Nov','Dez']

for m, fname in enumerate(files):
    print(f"\nMes {m+1:02d} ({meses[m]}): {os.path.basename(fname)}")
    if not os.path.exists(fname):
        print(f"  FALTANDO!")
        continue

    chl_grid     = read_l3b_fast(fname)
    data         = np.where(chl_grid > 0, chl_grid, 0.0)
    chl_smooth[m] = uniform_filter(data, size=2)
    v = chl_smooth[m][chl_smooth[m] > 0]
    if len(v):
        print(f"  Apos smooth: min={v.min():.4f}  max={v.max():.4f}  "
              f"media={v.mean():.4f} mg/m3")

# Salvar
print("\nSalvando...")
lat_out = np.linspace(-90,  90,  ny_model)
lon_out = np.linspace(-180, 180, nx_model)

with nc.Dataset('./INPUT/seawifs_1998-2006_smoothed_2X.nc', 'w',
                format='NETCDF3_64BIT_OFFSET') as ds:
    ds.createDimension('time', 12)
    ds.createDimension('lat',  ny_model)
    ds.createDimension('lon',  nx_model)

    t        = ds.createVariable('time', 'f4', ('time',))
    t[:]     = np.arange(1, 13, dtype=np.float32)
    t.units  = 'months since 1998-01-01'

    la        = ds.createVariable('lat', 'f4', ('lat',))
    la[:]     = lat_out.astype(np.float32)
    la.units  = 'degrees_north'

    lo        = ds.createVariable('lon', 'f4', ('lon',))
    lo[:]     = lon_out.astype(np.float32)
    lo.units  = 'degrees_east'

    c           = ds.createVariable('CHL_A', 'f4',
                                    ('time','lat','lon'), fill_value=0.0)
    c[:]        = chl_smooth
    c.units     = 'mg m-3'
    c.long_name = 'Chlorophyll-a (SeaWiFS climatology 1998-2006, smoothed 2X)'

print(f"\nCriado: ./INPUT/seawifs_1998-2006_smoothed_2X.nc  shape=(12,{ny_model},{nx_model})")
EOF

echo "-------------------3 -check data-----------------------"


python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ds  = nc.Dataset('./INPUT/seawifs_1998-2006_smoothed_2X.nc')
chl = ds['CHL_A'][:]
lat = ds['lat'][:]
lon = ds['lon'][:]
ds.close()

print("=== ESTRUTURA ===")
print(f"  shape:    {chl.shape}  (esperado: 12 x ny x nx)")
print(f"  lat:      {lat[0]:.1f} a {lat[-1]:.1f}")
print(f"  lon:      {lon[0]:.1f} a {lon[-1]:.1f}")

print("\n=== VALORES GLOBAIS ===")
validos = chl[chl > 0]
print(f"  celulas > 0:  {len(validos)} / {chl.size}")
print(f"  min (>0):     {validos.min():.4f} mg/m3")
print(f"  max:          {chl.max():.4f} mg/m3")
print(f"  media (>0):   {validos.mean():.4f} mg/m3")
print(f"  NaN:          {np.isnan(chl).sum()}")

# Valores esperados tipicos SeaWiFS:
# - Oceano aberto: 0.01 - 0.5 mg/m3
# - Costas/upwelling: 1 - 10 mg/m3
# - Maximo absoluto: ~60 mg/m3
print("\n=== SANIDADE DOS VALORES ===")
print(f"  valores < 0:       {(chl < 0).sum()}  (deve ser 0)")
print(f"  valores > 60:      {(chl > 60).sum()} (deve ser ~0)")
print(f"  valores 0.01-10:   {((chl>=0.01)&(chl<=10)).sum()} (maioria dos oceanos)")

print("\n=== POR MES ===")
meses = ['Jan','Fev','Mar','Abr','Mai','Jun',
         'Jul','Ago','Set','Out','Nov','Dez']
for m in range(12):
    v = chl[m][chl[m] > 0]
    if len(v) > 0:
        print(f"  {meses[m]}: min={v.min():.4f}  max={v.max():.4f}  "
              f"media={v.mean():.4f}  celulas={len(v)}")
    else:
        print(f"  {meses[m]}: SEM DADOS!")

# Verificar consistencia sazonal (NH mais alta no verao boreal)
print("\n=== CONSISTENCIA SAZONAL ===")
# Hemisferio Norte (lat > 30): maior chl em Mar-Mai
nh = np.where(lat > 30)[0]
sh = np.where(lat < -30)[0]
chl_nh = [chl[m][nh, :][chl[m][nh, :] > 0].mean()
          if (chl[m][nh, :] > 0).any() else 0 for m in range(12)]
chl_sh = [chl[m][sh, :][chl[m][sh, :] > 0].mean()
          if (chl[m][sh, :] > 0).any() else 0 for m in range(12)]
print(f"  NH max mes: {meses[np.argmax(chl_nh)]} (esperado: Mar-Mai)")
print(f"  SH max mes: {meses[np.argmax(chl_sh)]} (esperado: Set-Nov)")

# Salvar figuras
fig, axes = plt.subplots(3, 4, figsize=(16, 10))
axes = axes.flatten()
for m in range(12):
    im = axes[m].imshow(
        np.where(chl[m] > 0, np.log10(chl[m] + 1e-6), np.nan),
        origin='lower', aspect='auto',
        cmap='viridis', vmin=-2, vmax=1
    )
    axes[m].set_title(meses[m])
    axes[m].axis('off')
plt.suptitle('SeaWiFS CHL log10(mg/m3) ? 12 meses', fontsize=14)
plt.colorbar(im, ax=axes, label='log10(mg/m3)', shrink=0.6)
plt.savefig('check_seawifs.png', dpi=100, bbox_inches='tight')
print("\nFigura salva: check_seawifs.png")
print("\n=== RESULTADO ===")
ok = (
    chl.shape[0] == 12 and
    len(validos) > 0 and
    validos.min() > 0 and
    chl.max() < 100 and
    np.isnan(chl).sum() == 0
)
print("  Arquivo OK!" if ok else "  PROBLEMA DETECTADO ? verificar itens acima")
EOF

echo "-------------------4-----------------------"

python3 << 'EOF'
import netCDF4 as nc
import numpy as np

ds = nc.Dataset('./INPUT/seawifs_1998-2006_smoothed_2X.nc', 'a')
t = ds['time']
print("time units atual:", t.units)
print("time values:", t[:])

# FMS precisa de unidades absolutas com calendar
t.units    = 'days since 0001-01-01 00:00:00'
t.calendar = 'julian'
t.axis     = 'T'

# 12 meses: dia 15 de cada mes do ano 1
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
t[:] = np.array(days, dtype=np.float32)

ds.close()
print("time corrigido: days=", days)
EOF

echo "-------------------5 check calendar-----------------------"

python3 << 'EOF'
import netCDF4 as nc
import numpy as np

ds  = nc.Dataset('./INPUT/seawifs_1998-2006_smoothed_2X.nc')
chl = ds['CHL_A']
t   = ds['time']
lat = ds['lat']
lon = ds['lon']

print("=== CHECKLIST MOM6+SIS2 ===")

checks = {
    'shape (12,ny,nx)':      chl.shape[0] == 12,
    'sem NaN':               not np.isnan(chl[:]).any(),
    'valores > 0':           (chl[:] > 0).any(),
    'time tem 12 valores':   len(t[:]) == 12,
    'time tem units':        hasattr(t, 'units'),
    'time tem calendar':     hasattr(t, 'calendar'),
    'lat -90 a 90':          float(lat[0]) <= -89 and float(lat[-1]) >= 89,
    'lon cobre -180 a 180':  float(lon[0]) <= -179 and float(lon[-1]) >= 179,
    'variavel CHL_A':      'CHL_A' in ds.variables,
    'units mg/m3':           'mg' in getattr(chl, 'units', ''),
}

all_ok = True
for check, result in checks.items():
    status = "?" if result else "? PROBLEMA"
    print(f"  {status}  {check}")
    if not result:
        all_ok = False

ds.close()
print(f"\n{'PRONTO PARA MOM6+SIS2' if all_ok else 'CORRIGIR ITENS COM X ANTES DE USAR'}")
EOF



python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import shutil

fname = './INPUT/seawifs_1998-2006_smoothed_2X.nc'
shutil.copy(fname, fname+'.bak2')

# Ler dados do arquivo atual
ds_in   = nc.Dataset(fname+'.bak2')
lat     = np.array(ds_in['lat'][:])
lon     = np.array(ds_in['lon'][:])
data_var = next((v for v in ds_in.variables
                 if v not in ['time','lat','lon']), None)
data    = np.array(ds_in[data_var][:])
ds_in.close()

print(f"data_var={data_var}  shape={data.shape}")

# Dias do meio de cada mes
days_mid = np.array([15.5, 45.5, 74.5, 105.0, 135.5, 166.0,
                     196.5, 227.5, 258.0, 288.5, 319.0, 349.5])

# Salvar com time como dimensao UNLIMITED
with nc.Dataset(fname, 'w', format='NETCDF3_64BIT_OFFSET') as ds:

    # CRITICO: time deve ser UNLIMITED (sem tamanho fixo)
    ds.createDimension('time', None)   # None = UNLIMITED
    ds.createDimension('lat',  len(lat))
    ds.createDimension('lon',  len(lon))

    # Tempo
    t          = ds.createVariable('time', 'f8', ('time',))
    t[:]       = days_mid
    t.units    = 'days since 0001-01-01 00:00:00'
    t.calendar = 'julian'
    t.axis     = 'T'
    t.modulo   = ' '   # climatologia ciclica

    # Lat
    la         = ds.createVariable('lat', 'f4', ('lat',))
    la[:]      = lat
    la.units   = 'degrees_north'
    la.axis    = 'Y'

    # Lon
    lo         = ds.createVariable('lon', 'f4', ('lon',))
    lo[:]      = lon
    lo.units   = 'degrees_east'
    lo.axis    = 'X'

    # Dados
    v              = ds.createVariable(data_var, 'f4',
                     ('time','lat','lon'), fill_value=1e20)
    v[:]           = data
    v.units        = 'mg m-3'
    v.long_name    = 'Chlorophyll-a concentration'
    v.missing_value = np.float32(1e20)

    ds.title   = 'SeaWiFS 1998-2006 climatologia mensal para MOM6'
    ds.history = 'time=UNLIMITED + modulo para interpolacao ciclica FMS'

# Verificar
ds2 = nc.Dataset(fname)
print(f"\nVerificacao:")
print(f"  dims: { {k:len(v) for k,v in ds2.dimensions.items()} }")
ulim = [k for k,v in ds2.dimensions.items() if v.isunlimited()]
print(f"  UNLIMITED: {ulim}  (deve ser ['time'])")
print(f"  time: {np.array(ds2['time'][:])}")
print(f"  modulo: '{getattr(ds2['time'],'modulo','NAO DEFINIDO')}'")
print(f"  {data_var}: shape={ds2[data_var].shape}  "
      f"min={np.nanmin(ds2[data_var][:]):.3f}  "
      f"max={np.nanmax(ds2[data_var][:]):.3f}")
ds2.close()
EOF
rm ./INPUT/seawifs_1998-2006_smoothed_2X.nc.bak2
