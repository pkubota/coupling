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


#5. Difusividade - KH_background_2d.nc

rm ./INPUT/KH_background_2d.nc

python3 << 'EOF6'
import numpy as np
import netCDF4 as nc
import xarray as xr
# Grade do modelo
topog    = nc.Dataset('./INPUT/topog.nc')
ny_model = len(topog.dimensions['ny'])
nx_model = len(topog.dimensions['nx'])
topog.close()
print(f"Grade modelo: ny={ny_model}, nx={nx_model}")

grid = xr.open_dataset('./INPUT/ocean_hgrid.nc')
#hgrid = xr.open_dataset('ocean_hgrid.nc')
# supergrid tem nxp=361, nyp=317
# Pegar pontos T (every other point)
npy, npx =  grid['x'].shape

lon = grid['x'][1::2, 1::2]
lat = grid['y'][1::2, 1::2]

print(f"GradeG T: npy={npy}, npx={npx}")

lon2d = np.array(grid['x'][1::2, 1::2])   # (158, 180)
lat2d = np.array(grid['y'][1::2, 1::2])
lon1d = np.array(grid['x'][1, 1::2])  # (158, 180)
lat1d = np.array(grid['y'][1::2, 1])
grid.close()


ny, nx = lon2d.shape
print(f"Grade T: ny={ny}, nx={nx}")
print(f"lon: {lon2d.min():.1f} a {lon2d.max():.1f}")
print(f"lat: {lat2d.min():.1f} a {lat2d.max():.1f}")

# -- 2. Normalizar longitudes para -180 a 180
lon2d_norm = ((lon2d + 180) % 360) - 180
pts        = np.column_stack([lat2d.flatten(), lon2d_norm.flatten()])

#amp2_sum = np.zeros((ny_model, nx_model), dtype=np.float64)
#n_const  = 0

# KH background tipico: 100-1000 m*m/s, maior nos tropicos
KH = 200.0 * np.ones_like(lat.values)
KH[np.abs(lat.values) < 30] = 600.0  # maior nos tropicos

with nc.Dataset('./INPUT/KH_background_2d.nc', 'w',
                format='NETCDF3_64BIT_OFFSET') as ds:
    # CRITICO: time deve ser UNLIMITED (sem tamanho fixo)
    ds.createDimension('time', None)
    ds.createDimension('ny', lat.shape[0])
    ds.createDimension('nx', lat.shape[1])

    t        = ds.createVariable('time', 'f4', ('time',))
    t[:]     = np.arange(0, 1, dtype=np.float32)
    t.units    = 'days since 0001-01-01 00:00:00'
    
    la       = ds.createVariable('lat', 'f4', ('ny',))
    la[:]    = lat1d[:ny] 
    la.units = 'degrees_north'

    lo       = ds.createVariable('lon', 'f4', ('nx',))
    lo[:]    = lon1d[:nx]
    lo.units = 'degrees_east'

    kh = ds.createVariable('kh', 'f4', ('time','ny','nx'))
    kh[:] = KH
    kh.units = 'm2 s-1'
    kh.long_name = 'Background horizontal diffusivity'

EOF6

echo "-------------------3 plot-----------------------"

python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

kh_file  = nc.Dataset('./INPUT/KH_background_2d.nc')
kh       = kh_file['kh'][0,:,:]
kh_file.close()

print("=== Dimensoes ===")
print(f"  shape: {kh.shape}")
print(f"  esperado: (158, 180) ou (316, 360)")

print("\n=== Valores ===")
print(f"  min:  {kh.min():.1f} m2/s")
print(f"  max:  {kh.max():.1f} m2/s")
print(f"  mean: {kh.mean():.1f} m2/s")
print(f"  NaN:  {np.isnan(kh).sum()}")

print("\n=== Consistencia com topog ===")
topog = nc.Dataset('./INPUT/topog.nc')
depth = topog['depth'][:]
topog.close()
print(f"  topog shape: {depth.shape}")
if kh.shape == depth.shape:
    print("  shapes BATEM -")
    # KH deve ser 0 onde eh terra
    kh_em_terra = kh[depth == 0]
    print(f"  KH em terra (deve ser 0): min={kh.min():.1f} max={kh.max():.1f}")
    kh_no_oceano = kh[depth > 0]
    print(f"  KH no oceano: min={kh.min():.1f} max={kh.max():.1f}")
else:
    print(f"  ERRO: shapes nao batem! kh={kh.shape} topog={depth.shape}")

# Salvar figura para inspecao visual
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].imshow(kh, origin='lower', aspect='auto', cmap='viridis')
axes[0].set_title('KH_background_2d')
plt.colorbar(axes[0].images[0], ax=axes[0], label='m2/s')
axes[1].imshow(depth, origin='lower', aspect='auto', cmap='Blues')
axes[1].set_title('topog depth')
plt.colorbar(axes[1].images[0], ax=axes[1], label='m')
plt.tight_layout()
plt.savefig('check_KH.png', dpi=100)
print("\n  Figura salva: check_KH.png")
EOF
#O erro indica que kh[depth == 0] resultou em array vazio -
# ou seja, nao ha celulas de terra no topog! 
# O topog está todo oceano. Verificar:


