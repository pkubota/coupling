#!/bin/bash +x
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


################################################################################
# Ver onde estao os 29 pontos problematicos
################################################################################


python3 << 'EOF'
import numpy as np
import netCDF4 as nc

hgrid  = nc.Dataset('./ocean_hgrid.nc')
lon_T  = np.array(hgrid['x'][1::2, 1::2])
lat_T  = np.array(hgrid['y'][1::2, 1::2])
hgrid.close()

topog  = nc.Dataset('./topog.nc')
depth  = np.array(topog['depth'][:])
topog.close()

ds_t   = nc.Dataset('INPUT/mom6_cmems_mod_glo_phy_20260101.nc')
temp   = np.array(ds_t['ptemp'][0, 0])
ds_t.close()

ds_s   = nc.Dataset('INPUT/mom6_cmems_mod_glo_phy_20260101.nc')
salt   = np.array(ds_s['salt'][0, 0])
ds_s.close()
ocean  = depth > 0
bad    = ocean & ((temp > 1e10) | (salt > 1e10))
print(f"Pontos problematicos: {bad.sum()}")
j_b, i_b = np.where(ocean)
#for k in range(len(j_b)):
#    j, i = j_b[k], i_b[k]
#    print(f"{temp[j,i]:.1f}")

j_b, i_b = np.where(bad)
for k in range(len(j_b)):
    j, i = j_b[k], i_b[k]
    print(f"  j={j:3d} i={i:3d} lon={lon_T[j,i]:8.2f} "
          f"lat={lat_T[j,i]:8.3f} depth={depth[j,i]:8.2f}m "
          f"T={temp[j,i]:.1f} S={salt[j,i]:.1f}")

print(f"\nDepths dos pontos ruins:")
depths_bad = depth[bad]
#34 line
#print(f"  min={depths_bad.min():.2f}  max={depths_bad.max():.2f}  "
#      f"media={depths_bad.mean():.2f}")
#print(f"\nSugestao MINIMUM_DEPTH > {depths_bad.max():.1f}m")
EOF
exit
# Corrigir MOM_override -- remover duplicatas e corrigir nomes
python3 << 'EOF'
# Parametros corretos
correct = {
    'TEMP_Z_INIT_FILE':      '"mom6_cmems_mod_glo_phy_20260101.nc"',
    'SALT_Z_INIT_FILE':      '"mom6_cmems_mod_glo_phy_20260101.nc"',
    'Z_INIT_FILE_PTEMP_VAR': '"ptemp"',
    'Z_INIT_FILE_SALT_VAR':  '"salt"',
    'MINIMUM_DEPTH':         '20.0',
    'MASKING_DEPTH':         '19.9',
}
# Parametros a remover (incorretos ou desnecessarios)
remove = ['TIDEAMP_VARNAME', 'TOPO_VARNAME', 'ROUGHNESS_VARNAME']

with open('../MOM_override', 'r') as f:
    lines = f.readlines()

# Filtrar: remover todos os parametros que vamos reescrever
all_keys = list(correct.keys()) + remove
clean = []
seen  = set()
for line in lines:
    key = line.split('=')[0].strip()
    if key in all_keys:
        continue  # remover todas as ocorrencias
    clean.append(line)

# Adicionar parametros corretos uma vez
clean.append('\n! Arquivos de inicializacao T/S\n')
for k, v in correct.items():
    clean.append(f'{k} = {v}\n')

with open('MOM_override', 'w') as f:
    f.writelines(clean)

print("=== MOM_override final ===")
with open('MOM_override') as f:
    content = f.read()
print(content)
EOF
