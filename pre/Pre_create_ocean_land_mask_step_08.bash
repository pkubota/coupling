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




echo "-------------------0-----------------------"
rm ./INPUT/ocean_mask.nc
rm ./INPUT/land_mask.nc

python3 << 'EOF'
import numpy as np
import netCDF4 as nc

# Ler topog
topog = nc.Dataset('./INPUT/topog.nc')
depth = np.array(topog['depth'][:])
topog.close()

ny, nx = depth.shape
print(f"topog shape: ny={ny}, nx={nx}")
print(f"celulas oceano: {(depth>0).sum()}")
print(f"celulas terra:  {(depth==0).sum()}")

if (depth == 0).sum() == 0:
    print("AVISO: topog nao tem terra ? mascaras podem estar erradas!")

# Mascaras
ocean_mask = (depth > 0).astype(np.int32)  # 1=oceano, 0=terra
land_mask  = (depth == 0).astype(np.int32) # 1=terra,  0=oceano

for fname, data, varname, descr in [
    ('./INPUT/ocean_mask.nc', ocean_mask, 'mask',
     'ocean wet mask: 1=ocean, 0=land'),
    ('./INPUT/land_mask.nc',  land_mask,  'mask',
     'land mask: 1=land, 0=ocean'),
]:
    with nc.Dataset(fname, 'w', format='NETCDF3_64BIT_OFFSET') as ds:
        ds.createDimension('ny', ny)
        ds.createDimension('nx', nx)
        v            = ds.createVariable(varname, 'i4', ('ny','nx'))
        v[:]         = data
        v.long_name  = descr
        v.valid_min  = np.int32(0)
        v.valid_max  = np.int32(1)

    ones = (data == 1).sum()
    zeros = (data == 0).sum()
    print(f"{fname}: 1s={ones}  0s={zeros}")

print("\n./INPUT/ocean_mask.nc e ./INPUT/land_mask.nc criados.")

# Verificar
for fname in ['./INPUT/ocean_mask.nc', './INPUT/land_mask.nc']:
    ds = nc.Dataset(fname)
    m  = ds['mask'][:]
    print(f"\n{fname}:")
    print(f"  shape: {m.shape}")
    print(f"  unicos: {np.unique(m)}")
    print(f"  sum: {m.sum()}")
    ds.close()
EOF
#```
#
#Se o `grep` do `MOM_input` mostrar um nome diferente de variável,
# ajuste o `varname` no script. O mais comum no MOM6 é:
#```
#INPUTDIR/ocean_mask.nc  ?  variavel 'mask'  com 1=oceano, 0=terra


echo "-------------------0-----------------------"


python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Carregar todos os arquivos
topog  = nc.Dataset('./INPUT/topog.nc')
depth  = np.array(topog['depth'][:])
topog.close()

om = nc.Dataset('./INPUT/ocean_mask.nc')
lm = nc.Dataset('./INPUT/land_mask.nc')

# Detectar nome da variavel automaticamente
om_var = [v for v in om.variables if v not in ('lat','lon','nx','ny')][0]
lm_var = [v for v in lm.variables if v not in ('lat','lon','nx','ny')][0]

ocean_mask = np.array(om[om_var][:])
land_mask  = np.array(lm[lm_var][:])
om.close()
lm.close()

print("=" * 55)
print("=== 1. ESTRUTURA ===")
print(f"  topog shape:      {depth.shape}")
print(f"  ocean_mask shape: {ocean_mask.shape}  var='{om_var}'")
print(f"  land_mask shape:  {land_mask.shape}   var='{lm_var}'")
print(f"  shapes OK: {ocean_mask.shape == depth.shape == land_mask.shape}")

print("\n=== 2. VALORES ===")
print(f"  ocean_mask unicos: {np.unique(ocean_mask)}  (deve ser [0 1])")
print(f"  land_mask  unicos: {np.unique(land_mask)}   (deve ser [0 1])")
print(f"  ocean_mask NaN:    {np.isnan(ocean_mask.astype(float)).sum()}")
print(f"  land_mask  NaN:    {np.isnan(land_mask.astype(float)).sum()}")

print("\n=== 3. CONSISTENCIA ENTRE SI ===")
soma      = ocean_mask + land_mask
sobrepos  = (soma == 2).sum()   # celula ocean E land ao mesmo tempo
vazio     = (soma == 0).sum()   # celula nem ocean nem land
correto   = (soma == 1).sum()   # cada celula eh so ocean ou so land
print(f"  ocean+land=2 (sobreposicao): {sobrepos}  (deve ser 0)")
print(f"  ocean+land=0 (vazio):        {vazio}      (deve ser 0)")
print(f"  ocean+land=1 (correto):      {correto}    (deve ser {depth.size})")
print(f"  complementares: {sobrepos == 0 and vazio == 0}")

print("\n=== 4. CONSISTENCIA COM TOPOG ===")
ocean_from_topog = (depth > 0).astype(np.int32)
land_from_topog  = (depth == 0).astype(np.int32)
diff_om = (ocean_mask != ocean_from_topog).sum()
diff_lm = (land_mask  != land_from_topog).sum()
print(f"  ocean_mask difere do topog: {diff_om} celulas  (deve ser 0)")
print(f"  land_mask  difere do topog: {diff_lm} celulas  (deve ser 0)")

print("\n=== 5. COBERTURA GLOBAL ===")
n_ocean = ocean_mask.sum()
n_land  = land_mask.sum()
n_total = ocean_mask.size
print(f"  celulas oceano: {n_ocean:6d} / {n_total} ({100*n_ocean/n_total:.1f}%)")
print(f"  celulas terra:  {n_land:6d} / {n_total} ({100*n_land/n_total:.1f}%)")
# Globo real: ~71% oceano, ~29% terra
print(f"  % oceano esperada: ~60-75% para grade 1 grau")

print("\n=== 6. REGIOES CONHECIDAS ===")
hgrid      = nc.Dataset('./INPUT/ocean_hgrid.nc')
lat2d      = np.array(hgrid['y'][1::2, 1::2])
lon2d      = np.array(hgrid['x'][1::2, 1::2])
hgrid.close()
lon2d_norm = ((lon2d + 180) % 360) - 180

checks_geo = [
    ("Amazonia (5S,60W)",      (-8,-2),  (-65,-55), 0, "terra"),
    ("Pacifico central (0,180W)", (-5,5),(-180,-170),1, "oceano"),
    ("Antartica (85S,0)",      (-90,-80),(-10,10),   0, "terra/gelo"),
    ("Atlantico Norte (50N,30W)",(45,55),(-35,-25),  1, "oceano"),
    ("Sahara (20N,10E)",       (15,25),  (5,20),     0, "terra"),
    ("Indico (10S,70E)",       (-15,-5), (65,75),    1, "oceano"),
]

for nome, lat_r, lon_r, esperado, tipo in checks_geo:
    mask = ((lat2d >= lat_r[0]) & (lat2d <= lat_r[1]) &
            (lon2d_norm >= lon_r[0]) & (lon2d_norm <= lon_r[1]))
    if mask.sum() > 0:
        vals     = ocean_mask[mask]
        dominante = 1 if vals.mean() > 0.5 else 0
        ok       = dominante == esperado
        print(f"  {'?' if ok else '?'} {nome}: "
              f"{vals.mean()*100:.0f}% oceano  (esperado: {tipo})")

print("\n=== 7. CHECKLIST MOM6+SIS2 ===")
checks = {
    "shapes iguais ao topog":       ocean_mask.shape == depth.shape,
    "apenas 0s e 1s (ocean)":       set(np.unique(ocean_mask)) <= {0,1},
    "apenas 0s e 1s (land)":        set(np.unique(land_mask))  <= {0,1},
    "complementares (soma=1)":      sobrepos == 0 and vazio == 0,
    "consistente com topog (ocean)":diff_om == 0,
    "consistente com topog (land)": diff_lm == 0,
    "sem NaN":                      not np.isnan(ocean_mask.astype(float)).any(),
    "cobertura oceano 50-80%":      0.50 < n_ocean/n_total < 0.80,
}
all_ok = True
for check, result in checks.items():
    print(f"  {'?' if result else '? PROBLEMA'}  {check}")
    if not result:
        all_ok = False

# Figura
fig, axes = plt.subplots(1, 3, figsize=(16, 4))
axes[0].imshow(depth,      origin='lower', aspect='auto', cmap='Blues')
axes[0].set_title('topog depth')
axes[1].imshow(ocean_mask, origin='lower', aspect='auto', cmap='Blues', vmin=0, vmax=1)
axes[1].set_title('ocean_mask (1=oceano)')
axes[2].imshow(land_mask,  origin='lower', aspect='auto', cmap='Greens', vmin=0, vmax=1)
axes[2].set_title('land_mask (1=terra)')
for ax in axes:
    ax.axis('off')
plt.tight_layout()
plt.savefig('check_masks.png', dpi=120)
print(f"\nFigura: check_masks.png")
print(f"\n{'PRONTO PARA MOM6+SIS2' if all_ok else 'CORRIGIR ITENS COM X'}")
EOF



#Os critérios críticos para o MOM6+SIS2 são:
#Criterio                                 Por que e importante
#Shapes iguais ao topog                   FMS rejeita se nao bater
#Apenas 0s e 1s                           modelo trata como inteiro booleano
#Complementares                           cada celula deve ser exatamente ocean OU land
#Consistente com topog                    mascara deve derivar da mesma batimetria
#Amazonia=terra, Pacífico=oceano          valida geografia basica
