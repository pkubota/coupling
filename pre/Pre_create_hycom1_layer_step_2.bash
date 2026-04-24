!/bin/bash +x
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
echo    " =================================================="
echo    " ------step 01 -> create_ hycom1------"
echo    " =================================================="
rm ./INPUT/hycom1_75_800m_mod.nc
rm ./INPUT/hycom1_75_800m.nc
rm ./INPUT/layer_coord_vgrid.nc
rm ./INPUT/hycom1_75_800m_vgrid.nc
rm ./INPUT/hycom1_75_800m_mod_vgrid.nc

python3 << 'EOF'
import numpy as np
import netCDF4 as nc

dz = np.array([
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2.01, 2.01, 2.01, 2.02, 2.03, 2.05, 2.07, 2.09, 2.13, 2.18,
    2.24, 2.32, 2.42, 2.56, 2.73, 2.95, 3.23, 3.57, 3.99, 4.52,
    5.16, 5.95, 6.91, 8.06, 9.46, 11.13, 13.13, 15.5, 18.32, 21.65,
    25.56, 30.16, 35.54, 41.81, 49.1, 57.56, 67.34, 78.63, 91.61,
    106.51, 123.58, 143.08, 165.32, 190.61, 219.33, 251.88, 288.68,
    330.23, 377.05, 429.71, 488.84, 555.13, 629.33, 712.25, 804.76
])

sigma2_iface = np.array([
    1010, 1014.3034, 1017.8088, 1020.843, 1023.5566, 1025.813,
    1027.0275, 1027.9114, 1028.6422, 1029.2795, 1029.852, 1030.3762,
    1030.8626, 1031.3183, 1031.7486, 1032.1572, 1032.5471, 1032.9207,
    1033.2798, 1033.6261, 1033.9608, 1034.2519, 1034.4817, 1034.6774,
    1034.8508, 1035.0082, 1035.1533, 1035.2886, 1035.4159, 1035.5364,
    1035.6511, 1035.7608, 1035.8661, 1035.9675, 1036.0645, 1036.1554,
    1036.2411, 1036.3223, 1036.3998, 1036.4739, 1036.5451, 1036.6137,
    1036.68, 1036.7441, 1036.8062, 1036.8526, 1036.8874, 1036.9164,
    1036.9418, 1036.9647, 1036.9857, 1037.0052, 1037.0236, 1037.0409,
    1037.0574, 1037.0738, 1037.0902, 1037.1066, 1037.123, 1037.1394,
    1037.1558, 1037.1722, 1037.1887, 1037.206, 1037.2241, 1037.2435,
    1037.2642, 1037.2866, 1037.3112, 1037.3389, 1037.3713, 1037.4118,
    1037.475, 1037.6332, 1037.8104, 1038
])

# Layer = media entre interfaces consecutivas (75 valores)
layer = 0.5*(sigma2_iface[:-1] + sigma2_iface[1:])

print(f"dz:     len={len(dz)}")
print(f"sigma2: len={len(sigma2_iface)} (interfaces)")
print(f"Layer:  len={len(layer)}  (camadas)")
print(f"Layer[0]={layer[0]:.4f}  Layer[-1]={layer[-1]:.4f}")

for fname in ['INPUT/hycom1_75_800m_mod.nc',
             'INPUT/hycom1_75_800m.nc']:
    with nc.Dataset(fname, 'w', format='NETCDF3_64BIT_OFFSET') as ds:
        ds.createDimension('layers',     75)
        ds.createDimension('interfaces', 76)

        # dz ? espessura das camadas
        vdz           = ds.createVariable('dz', 'f8', ('layers',))
        vdz[:]        = dz
        vdz.long_name = 'z* coordinate level thickness'
        vdz.units     = 'm'

        # sigma2 nas interfaces (76 valores)
        vs2           = ds.createVariable('sigma2', 'f8', ('interfaces',))
        vs2[:]        = sigma2_iface
        vs2.long_name = 'Interface target potential density referenced to 2000 dbars'
        vs2.units     = 'kg/m3'

        # Layer = densidade no centro de cada camada (75 valores)
        vlay          = ds.createVariable('Layer', 'f8', ('layers',))
        vlay[:]       = layer
        vlay.long_name = 'Layer target potential density referenced to 2000 dbars'
        vlay.units    = 'kg/m3'

    ds2 = nc.Dataset(fname)
    print(f"\nCriado: {fname}")
    print(f"  dims: { {k:len(v) for k,v in ds2.dimensions.items()} }")
    print(f"  vars: {list(ds2.variables.keys())}")
    print(f"  Layer:  len={len(ds2['Layer'][:])}   "
          f"[{ds2['Layer'][0]:.4f}...{ds2['Layer'][-1]:.4f}]")
    print(f"  sigma2: len={len(ds2['sigma2'][:])}  "
          f"[{ds2['sigma2'][0]:.4f}...{ds2['sigma2'][-1]:.4f}]")
    print(f"  dz:     len={len(ds2['dz'][:])}   sum={ds2['dz'][:].sum():.1f}m")
    ds2.close()

print("\nArquivos prontos!")
EOF


echo    " ==============================================================="
echo    " ------step 01 -> ceate layer.nc do arquivo vgrid_75_2m.nc------"
echo    " ==============================================================="

python3 << 'EOF'
import numpy as np
import netCDF4 as nc

# Ler vgrid
vgrid   = nc.Dataset('./INPUT/vgrid_75_2m.nc')
zeta    = np.array(vgrid['zeta'][:])
vgrid.close()
z_iface = zeta[0::2]
z_mid   = 0.5*(z_iface[:-1] + z_iface[1:])
dz      = -(z_iface[:-1] - z_iface[1:])
nk      = len(z_mid)
print(f"nk={nk}")

# Coordenada sigma-2
sigma2 = np.zeros(nk)
for k in range(nk):
    z = z_mid[k]
    if z <= 200.0:
        sigma2[k] = 28.0 + (z/200.0)*5.0
    elif z <= 1000.0:
        sigma2[k] = 33.0 + ((z-200.0)/800.0)*3.0
    else:
        sigma2[k] = 36.0 + ((z-1000.0)/5500.0)*1.5
sigma2 = np.clip(sigma2, 28.0, 37.8)

# Criar os dois arquivos com variavel 'Layer' E 'sigma2'
for fname in ['./INPUT/hycom1_75_800m_vgrid.nc',
              './INPUT/hycom1_75_800m_mod_vgrid.nc']:
    with nc.Dataset(fname, 'w', format='NETCDF3_64BIT_OFFSET') as ds:
        ds.createDimension('Layer',     nk)
        ds.createDimension('Interface', nk+1)

        # Variavel 'Layer' -- nome exato que o FMS procura
        lay          = ds.createVariable('Layer', 'f8', ('Layer',))
        lay[:]       = sigma2
        lay.units    = 'kg m-3'
        lay.long_name = 'Target potential density (sigma-2)'

        # Variavel 'sigma2' -- nome usado pelo ALE_COORDINATE_CONFIG
        s2           = ds.createVariable('sigma2', 'f8', ('Layer',))
        s2[:]        = sigma2
        s2.units     = 'kg m-3'
        s2.long_name = 'Potential density referenced to 2000 dbar'

        # Interface
        iface_vals   = np.concatenate([
            [sigma2[0] - 0.5*(sigma2[1]-sigma2[0])],
            0.5*(sigma2[:-1] + sigma2[1:]),
            [sigma2[-1] + 0.5*(sigma2[-1]-sigma2[-2])]
        ])
        iface        = ds.createVariable('Interface', 'f8', ('Interface',))
        iface[:]     = iface_vals
        iface.units  = 'kg m-3'

        ds.title     = f'HYCOM1 {nk}-layer coordinate for MOM6 OM_1deg'

    # Verificar
    ds2 = nc.Dataset(fname)
    print(f"Criado: {fname}")
    print(f"  dims: { {k:len(v) for k,v in ds2.dimensions.items()} }")
    print(f"  vars: {list(ds2.variables.keys())}")
    print(f"  Layer[0]={np.array(ds2['Layer'][0]):.4f}  "
          f"  Layer[-1]={np.array(ds2['Layer'][-1]):.4f}")
    ds2.close()

# Criar os dois arquivos com variavel 'Layer' E 'sigma2'
for fname in ['./INPUT/layer_coord_vgrid.nc']:
    with nc.Dataset(fname, 'w', format='NETCDF3_64BIT_OFFSET') as ds:
        ds.createDimension('Layer',     nk)

        # Variavel 'Layer' -- nome exato que o FMS procura
        lay          = ds.createVariable('Layer', 'f8', ('Layer',))
        lay[:]       = sigma2
        lay.units    = 'kg m-3'
        lay.long_name = 'Target potential density (sigma-2)'

        # Variavel 'sigma2' -- nome usado pelo ALE_COORDINATE_CONFIG
        s2           = ds.createVariable('sigma2', 'f8', ('Layer',))
        s2[:]        = sigma2
        s2.units     = 'kg m-3'
        s2.long_name = 'Potential density referenced to 2000 dbar'

        # Variavel 'dz' -- nome usado pelo ALE_COORDINATE_CONFIG
        dz2           = ds.createVariable('dz', 'f8', ('Layer',))
        dz2[:]        = dz
        dz2.units     = 'm'
        dz2.long_name = 'dz'

        ds.title     = f'HYCOM1 {nk}-layer coordinate for MOM6 OM_1deg'

    # Verificar
    ds2 = nc.Dataset(fname)
    print(f"Criado: {fname}")
    print(f"  dims: { {k:len(v) for k,v in ds2.dimensions.items()} }")
    print(f"  vars: {list(ds2.variables.keys())}")
    print(f"  Layer[0]={np.array(ds2['Layer'][0]):.4f}  "
          f"  Layer[-1]={np.array(ds2['Layer'][-1]):.4f}")
    ds2.close()



EOF


echo    " ==============================================================="
echo    " ------step 02 -> check  NAN no arquivo hycom1_75_800m*.nc------"
echo    " ==============================================================="



# O NaN em reproducing_sum(_3d) ocorre ao calcular massa total
# Isso significa que h (espessura das camadas) tem NaN
# Isso vem da inicializacao ALE/HYCOM

# Ver se o problema e no arquivo hycom
python3 << 'EOF'
import numpy as np
import netCDF4 as nc

# Ver o arquivo hycom
ds = nc.Dataset('./INPUT/hycom1_75_800m_mod.nc')
print(f"dims: { {k:len(v) for k,v in ds.dimensions.items()} }")
for v in ds.variables:
    arr = np.array(ds[v][:])
    print(f"  {v}: shape={arr.shape} min={arr.min():.4f} max={arr.max():.4f} "
          f"NaN={np.isnan(arr).sum()} monotonic={np.all(np.diff(arr)>0)}")
ds.close()

# Verificar: sigma2 deve ser ESTRITAMENTE crescente
ds2 = nc.Dataset('./INPUT/hycom1_75_800m_mod.nc')
sig = np.array(ds2['sigma2'][:])
dif = np.diff(sig)
print(f"\nsigma2 interfaces: {len(sig)} valores")
print(f"diferencas: min={dif.min():.6f} max={dif.max():.6f}")
print(f"monotonica crescente: {np.all(dif>0)}")
if not np.all(dif>0):
    idx = np.where(dif<=0)[0]
    print(f"Problemas em indices: {idx}")
    for i in idx:
        print(f"  sigma2[{i}]={sig[i]:.6f} sigma2[{i+1}]={sig[i+1]:.6f}")
ds2.close()
EOF

echo    " =================================================="
echo    " Verificar ./hycom1_75_800m_mod.nc ./vgrid_75_2m.nc"
echo    " =================================================="
#
python3 << 'EOF'
import numpy as np
import netCDF4 as nc

# Verificar arquivo atual
ds  = nc.Dataset('./INPUT/hycom1_75_800m_mod.nc')
sig = np.array(ds['sigma2'][:])
print(f"sigma2 atual: len={len(sig)}  (deve ser 75)")
print(f"Layer atual:  len={len(ds['Layer'][:])}")
ds.close()

# Verificar vgrid
vgrid   = nc.Dataset('./INPUT/vgrid_75_2m.nc')
zeta    = np.array(vgrid['zeta'][:])
vgrid.close()
nzv     = len(zeta)
nk      = (nzv - 1) // 2
z_iface = zeta[0::2]
z_mid   = 0.5*(z_iface[:-1] + z_iface[1:])
print(f"\nvgrid: nzv={nzv}  nk={nk}")
print(f"z_mid: {z_mid[0]:.4f} a {z_mid[-1]:.1f} m")

# Verificar DIAG_COORD_DEF_Z que usa 'dz' do vgrid
dz = np.diff(z_iface)
print(f"dz: len={len(dz)}  (deve ser {nk})")
print(f"dz[0]={dz[0]:.4f}m  dz[-1]={dz[-1]:.1f}m")
EOF
#

exit



































































exit
#
echo    " =================================================="
echo    " recriar ./hycom1_75_800m_mod.nc ./vgrid_75_2m.nc"
echo    " =================================================="
exit
#

python3 << 'EOF'
import numpy as np
import netCDF4 as nc

NK = 75

# Ler vgrid
vgrid   = nc.Dataset('./INPUT/vgrid_75_2m.nc')
zeta    = np.array(vgrid['zeta'][:])
vgrid.close()
z_iface = zeta[0::2]              # 76 interfaces
z_mid   = 0.5*(z_iface[:-1] + z_iface[1:])  # 75 centros
dz      = np.diff(z_iface)        # 75 espessuras
assert len(z_mid) == NK, f"z_mid len={len(z_mid)} != {NK}"
print(f"nk={NK}  z_mid[0]={z_mid[0]:.4f}  z_mid[-1]={z_mid[-1]:.1f}")

# Coordenada sigma-2 com exatamente NK=75 valores
sigma2 = np.zeros(NK)
for k in range(NK):
    z = z_mid[k]
    if z <= 200.0:
        sigma2[k] = 28.0 + (z/200.0)*5.0
    elif z <= 1000.0:
        sigma2[k] = 33.0 + ((z-200.0)/800.0)*3.0
    else:
        sigma2[k] = 36.0 + ((z-1000.0)/5500.0)*1.5
sigma2 = np.clip(sigma2, 28.0, 37.8)
assert len(sigma2) == NK
print(f"sigma2: len={len(sigma2)}  [{sigma2[0]:.4f} ... {sigma2[-1]:.4f}]")

# Interfaces (NK+1 = 76 valores)
iface = np.zeros(NK+1)
iface[0]    = sigma2[0] - 0.5*(sigma2[1]-sigma2[0])
iface[1:-1] = 0.5*(sigma2[:-1] + sigma2[1:])
iface[-1]   = sigma2[-1] + 0.5*(sigma2[-1]-sigma2[-2])
assert len(iface) == NK+1

# Criar os dois arquivos com exatamente NK camadas
for fname in ['./INPUT/hycom1_75_800m_mod.nc',
             './INPUT/hycom1_75_800m.nc']:
    with nc.Dataset(fname, 'w', format='NETCDF3_64BIT_OFFSET') as ds:
        ds.createDimension('Layer',     NK)
        ds.createDimension('Interface', NK+1)

        # 'Layer' -- lido por COORD_FILE
        lay          = ds.createVariable('Layer', 'f8', ('Layer',))
        lay[:]       = sigma2
        lay.units    = 'kg m-3'
        lay.long_name = 'Target potential density (sigma-2)'

        # 'sigma2' -- lido por ALE_COORDINATE_CONFIG
        s2           = ds.createVariable('sigma2', 'f8', ('Layer',))
        s2[:]        = sigma2
        s2.units     = 'kg m-3'
        s2.long_name = 'Potential density referenced to 2000 dbar'

        # Interface
        si           = ds.createVariable('Interface', 'f8', ('Interface',))
        si[:]        = iface
        si.units     = 'kg m-3'
        si.long_name = 'Interface potential density (sigma-2)'

        ds.nk        = NK
        ds.title     = f'HYCOM1 {NK}-layer sigma2 coordinate for MOM6 OM_1deg'

    # Verificar imediatamente
    ds2 = nc.Dataset(fname)
    lay_check = np.array(ds2['Layer'][:])
    sig_check = np.array(ds2['sigma2'][:])
    print(f"\nCriado: {fname}")
    print(f"  dims: { {k:len(v) for k,v in ds2.dimensions.items()} }")
    print(f"  Layer:  len={len(lay_check)}  OK={len(lay_check)==NK}")
    print(f"  sigma2: len={len(sig_check)}  OK={len(sig_check)==NK}")
    ds2.close()

# Criar layer_coord.nc tambem
fname3 = './INPUT/layer_coord.nc'
with nc.Dataset(fname3, 'w', format='NETCDF3_64BIT_OFFSET') as ds:
    ds.createDimension('Layer', NK)

    lay      = ds.createVariable('Layer',  'f8', ('Layer',))
    lay[:]   = sigma2
    lay.units = 'kg m-3'

    s2       = ds.createVariable('sigma2', 'f8', ('Layer',))
    s2[:]    = sigma2
    s2.units = 'kg m-3'

    # dz -- usado por DIAG_COORD_DEF_Z = "FILE:vgrid_75_2m.nc,dz"
    dzv      = ds.createVariable('dz', 'f8', ('Layer',))
    dzv[:]   = dz
    dzv.units = 'm'

print(f"\nCriado: {fname3}")
print(f"  Layer: len={NK}  dz: len={len(dz)}")
print("\nTodos os arquivos criados com NK=75!")
EOF


echo    " =================================================="
echo    " Forcar recriacao sem symlink ./hycom1_75_800m_mod.nc ./vgrid_75_2m.nc"
echo    " =================================================="



# Forcar recriacao sem symlink
rm -f INPUT/hycom1_75_800m_mod.nc

python3 << 'EOF'
import numpy as np
import netCDF4 as nc

NK = 75

vgrid   = nc.Dataset('./INPUT/vgrid_75_2m.nc')
zeta    = np.array(vgrid['zeta'][:])
vgrid.close()
z_iface = zeta[0::2]
z_mid   = 0.5*(z_iface[:-1] + z_iface[1:])
assert len(z_mid) == NK, f"ERRO: len(z_mid)={len(z_mid)} != {NK}"

sigma2 = np.zeros(NK)
for k in range(NK):
    z = z_mid[k]
    if z <= 200.0:   sigma2[k] = 28.0 + (z/200.0)*5.0
    elif z <= 1000.0: sigma2[k] = 33.0 + ((z-200.0)/800.0)*3.0
    else:             sigma2[k] = 36.0 + ((z-1000.0)/5500.0)*1.5
sigma2 = np.clip(sigma2, 28.0, 37.8)

fname = './INPUT/hycom1_75_800m_mod.nc'
with nc.Dataset(fname, 'w', format='NETCDF3_64BIT_OFFSET') as ds:
    ds.createDimension('Layer',     NK)
    ds.createDimension('Interface', NK+1)

    lay       = ds.createVariable('Layer',  'f8', ('Layer',))
    lay[:]    = sigma2
    lay.units = 'kg m-3'

    s2        = ds.createVariable('sigma2', 'f8', ('Layer',))
    s2[:]     = sigma2
    s2.units  = 'kg m-3'
    s2.long_name = 'Potential density referenced to 2000 dbar'

    iface_v   = np.zeros(NK+1)
    iface_v[0]    = sigma2[0] - 0.5*(sigma2[1]-sigma2[0])
    iface_v[1:-1] = 0.5*(sigma2[:-1]+sigma2[1:])
    iface_v[-1]   = sigma2[-1] + 0.5*(sigma2[-1]-sigma2[-2])
    si        = ds.createVariable('Interface', 'f8', ('Interface',))
    si[:]     = iface_v
    si.units  = 'kg m-3'

# Verificar IMEDIATAMENTE apos criar
import os
fsize = os.path.getsize(fname)
ds2   = nc.Dataset(fname)
n_lay = len(ds2['Layer'][:])
n_sig = len(ds2['sigma2'][:])
ds2.close()

print(f"Arquivo: {fname}")
print(f"  tamanho: {fsize} bytes")
print(f"  Layer:  {n_lay}  (deve ser {NK})")
print(f"  sigma2: {n_sig}  (deve ser {NK})")
print(f"  {'OK' if n_lay==NK and n_sig==NK else 'ERRO!'}")
EOF
#

echo    " =================================================="
echo    " Verificar antes de rodar ./hycom1_75_800m_mod.nc ./vgrid_75_2m.nc"
echo    " =================================================="

#
python3 -c "
import netCDF4 as nc
for f in ['./INPUT/hycom1_75_800m_mod.nc',
          './INPUT/hycom1_75_800m.nc',
          './INPUT/layer_coord.nc']:
    ds = nc.Dataset(f)
    print(f'{f}: Layer={len(ds[\"Layer\"][:])}  sigma2={len(ds[\"sigma2\"][:])}')
    ds.close()
"
echo    " =================================================="
echo    " Confirmar que nao e symlink e tem tamanho correto ./hycom1_75_800m_mod.nc ./vgrid_75_2m.nc"
echo    " =================================================="

ls -lh ./INPUT/hycom1_75_800m_mod.nc
python3 -c "
import netCDF4 as nc
ds = nc.Dataset('./INPUT/hycom1_75_800m_mod.nc')
print('Layer:', len(ds['Layer'][:]))
print('sigma2:', len(ds['sigma2'][:]))
ds.close()
"
