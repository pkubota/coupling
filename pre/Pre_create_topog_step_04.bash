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
################################################################################
# PASSO 1 - Visão Geral e Dependencias
################################################################################
#
#Projeto MONAN - INPE/CPTEC
#Referencia: configuracao OM_1deg do GFDL
#Conteudo: topog.nc            · ocean_hgrid.nc          · ocean_mosaic.nc   · grid_spec.nc       · 
#          vgrid_75_2m.nc      · layer_coord.nc          · hycom1_75_800m.nc · tidal_amplitude.nc · 
#          KH_background_2d.nc · seawifs_1998-2006_smoothed_2X.nc ·
#MOM_channels · atmos/land mosaic exchange grids

#   1. Visão Geral e Dependencias

#Os arquivos de INPUT do OM_1deg sao gerados por tres ferramentas principais:
#Ferramenta        |                       Uso                          | Repositorio
#-----------------------------------------------------------------------------------------------------------
#FRE-NCtools       |  Grades, mosaicos, topografia                      | github.com/NOAA-GFDL/FRE-NCtools
#Python            |  (xarray/netCDF4) Grades verticais, KH, clorofila  | pip install xarray netCDF4 scipy
#NCO               |  (ncrcat, ncap2) Pos-processamento NetCDF          | nco.sourceforge.net
#ESMF              |  Exchange grids (opcional)                         | earthsystemmodeling.org

#1.1 Instalacao do FRE-NCtools
# Dependencias (modulos no ian05)
#module load intel/2022 openmpi/4.1 netcdf/4.9 hdf5/1.12
# Clonar e compilar
#git clone https://github.com/NOAA-GFDL/FRE-NCtools.git
#cd FRE-NCtools
#autoreconf -i
#./configure --prefix=$HOME/fre-nctools \
#CC=icc FC=ifort \
#CFLAGS="-O2" FCFLAGS="-O2"
#make -j8 && make install
# Adicionar ao PATH
#export PATH=$HOME/fre-nctools/bin:$PATH

## Ordem correta de geração de TODOS os arquivos
#```
#1. make_hgrid          -> ocean_hgrid.nc
#2. python vgrid.py     -> vgrid_75_2m.nc  (necessário para make_topog)
#3. make_solo_mosaic    -> ocean_mosaic.nc
#4. make_topog          -> topog.nc
#5. make_solo_mosaic    -> atmos_mosaic.nc
#6. make_coupler_mosaic -> grid_spec.nc + exchange grids
#7. python restante     -> KH, tidal, seawifs, masks, channels

################################################################################
# PASSO 2 - Grade Horizontal: ocean_hgrid.nc
################################################################################
#O arquivo ocean_hgrid.nc define a grade supergrid 2x (720x360 pontos de interface
# para uma grade 1-grau com 360x180
#celulas T). Gerado com make_hgrid.
#
#    2.1 Geracao com make_hgrid (tripolar - padrao OM4)
#
# Grade tripolar 1 grau (padrão GFDL OM4)
#                    Usage of make_hgrid

#  make_hgrid --grid_type grid_type --my_grid_file my_grid_file
#                 --nxbnds nxbnds --nybnds nybnds
#                 --xbnds x(1),...,x(nxbnds) --ybnds y(1),...,y(nybnds)
#                 --nlon nlon(1),...nlon(nxbnds-1)
#                 --nlat nlat(1),...nlat(nybnds-1)
#                 --dlon dlon(1),...dlon(nxbnds)
#                 --dlat dlat(1),...dlat(nybnds)
#                 --lat_join lat_join --num_lon num_lon --nratio nratio
#                 --simple_dx simple_dx --simple_dy simple_dy
#                 --grid_name gridname --center center --verbose --shift_fac #
#                 --do_schmidt --stretch_fac # --target_lon # --target_lat #
#                 --do_cube_transform
#                 --nest_grids nests
#                 --parent_tile parent_tile(1),...parent_tile(nests-1)
#                 --refine_ratio refine_ratio(1),...refine_ratio(nests-1)
#                 --halo #
#                 --istart_nest istart_nest(1),...istart_nest(nests-1)
#                 --iend_nest iend_nest(1),...iend_nest(nests-1)
#                 --jstart_nest jstart_nest(1),...jstart_nest(nests-1)
#                 --jend_nest jend_nest(1),...jend_nest(nests-1)
#                 --great_circle_algorithm --out_halo #
# -- 1. Grade horizontal tripolar 
#O erro  que --nlon e --nlat precisam ser pares quando center nao e none.
#Arquivos de entrada MOM6+SIS2 OM_1deg
################################################################################
# PASSO 4 - Topografia ocean
################################################################################
echo "--------------------------0---------------------------"
echo "     4.  Topografia ocean  topog.nc                   "
echo "--------------------------0---------------------------"
#3. Topografia: topog.nc e topo_edits_011818.nc
#   Define a profundidade do oceano em cada celula. 
#   Requer dado batimetrico bruto (GEBCO).
#3.1 Baixar GEBCO
#    GEBCO 2023 (15 arc-second, ~8GB)(passwrd -> pykDCS1975) -> requer registro gratuito em:
#    https://www.gebco.net/data_and_products/gridded_bathymetry_data/
#    wget  -c --timeout=60 -o wget.log -t 80 \
#    --dot-style=mega --load-cookies ~/.urs_cookies \
#    --save-cookies ~/.urs_cookies --auth-no-challenge=on \
#    --keep-session-cookies https://www.bodc.ac.uk/data/open_download/gebco/gebco_2023/zip/ \
#    -O GEBCO_2023.zip
#    https://dap.ceda.ac.uk/bodc/gebco/global/gebco_2025/sub_ice_topography_bathymetry/netcdf/gebco_2025_sub_ice_topo.zip?download=1
#    Versão de baixa resolução para teste (30 arc-second):
#    https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_30_second_grid/
rm  ./INPUT/topog.nc
python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import uniform_filter

# - 1. Ler pontos T do hgrid (supergrid: pegar 1 de cada 2) 

# Ler pontos T EXATOS do ocean_hgrid (supergrid -> grade T)
hgrid = nc.Dataset('./INPUT/ocean_hgrid.nc')
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
hgrid.close()

ny, nx = lon2d.shape
print(f"Grade T: ny={ny}, nx={nx}")
print(f"lon: {lon2d.min():.1f} a {lon2d.max():.1f}")
print(f"lat: {lat2d.min():.1f} a {lat2d.max():.1f}")

# -- 2. Normalizar longitudes para -180 a 180
lon2d_norm = ((lon2d + 180) % 360) - 180

# -- 3. Ler GEBCO 
print("Lendo GEBCO...")
gebco = nc.Dataset('GEBCO/GEBCO_2025_sub_ice_fixed.nc')
lat_g = np.array(gebco['lat'][:])
lon_g = np.array(gebco['lon'][:])
elev  = np.array(gebco['elevation'][:], dtype=np.float32)
gebco.close()

print("GEBCO: negativo=oceano, positivo=terra")
print(f"GEBCO elevation: min={elev.min():.0f} max={elev.max():.0f}")
print(f"GEBCO: lat {lat_g[0]:.1f}:{lat_g[-1]:.1f}, lon {lon_g[0]:.1f}:{lon_g[-1]:.1f}")

# Garantir lat crescente
if lat_g[0] > lat_g[-1]:
    lat_g = lat_g[::-1]
    elev  = elev[::-1, :]
# -- 4. Interpolar com RegularGridInterpolator
 
print("Criando interpolador...")
interp = RegularGridInterpolator(
    (lat_g, lon_g),
    elev,
    method='linear',
    bounds_error=False,
    fill_value=np.nan
)

print("Interpolando...")
pts        = np.column_stack([lat2d.flatten(), lon2d_norm.flatten()])
depth_flat = interp(pts)
depth_2d   = depth_flat.reshape(ny, nx)

# -- 5. Tratar NaN (pontos fora do dominio GEBCO)
 
nan_mask = np.isnan(depth_2d)
print(f"Pontos NaN: {nan_mask.sum()} de {ny*nx}")
depth_2d = np.where(nan_mask, 0.0, depth_2d)

# -- 6. Converter elevacao -> profundidade 
# oceano = elev < 0  -> depth = -elev (positivo)
# terra  = elev >= 0 -> depth = 0
depth_val  = np.where(depth_2d < 0, -depth_2d, 0.0)
mask_ocean = depth_val > 0

# -- 7. Suavizar 

depth_smooth = uniform_filter(depth_val, size=3)
depth_smooth = np.where(mask_ocean, np.maximum(depth_smooth, 10.0), 0.0)

print(f"Prof max:      {depth_smooth.max():.1f} m")
print(f"Celulas oceano:{mask_ocean.sum()} / {ny*nx}")

# -- 8. Salvar topog.nc 
with nc.Dataset('./INPUT/topog_bruto.nc', 'w', format='NETCDF3_64BIT_OFFSET') as ds:
    ds.createDimension('ny', ny)
    ds.createDimension('nx', nx)
    dep            = ds.createVariable('depth', 'f4', ('ny','nx'), fill_value=0.0)
    dep[:]         = depth_smooth.astype(np.float32)
    dep.units      = 'm'
    dep.long_name  = 'ocean depth'

print(f"./INPUT/topog.nc criado: ny={ny}, nx={nx}")

EOF

###############################################################
#  Identificar pontos problematicos
###############################################################

python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import shutil

# Ler grade e topografia
hgrid  = nc.Dataset('./INPUT/ocean_hgrid.nc')
lon_T  = np.array(hgrid['x'][1::2, 1::2])
lat_T  = np.array(hgrid['y'][1::2, 1::2])
hgrid.close()

shutil.copy('./INPUT/topog_bruto.nc', './INPUT/topog_bruto.nc.bak')
topog  = nc.Dataset('./INPUT/topog_bruto.nc.bak')
depth  = np.array(topog['depth'][:])
topog.close()

# Identificar pontos problematicos
# 6 pontos: lat < -80, depth <= 10m (MINIMUM_DEPTH original = 9.5m)
#prob = (lat_T < -80) & (depth > 0) & (depth <= 10.5)
prob =  (depth > 0) & (depth <= 10.5)
print(f"Pontos problematicos: {prob.sum()}")

j_idx, i_idx = np.where(prob)
for k in range(len(j_idx)):
    j, i = j_idx[k], i_idx[k]
    print(f"  j={j} i={i} lon={lon_T[j,i]:.2f} "
          f"lat={lat_T[j,i]:.3f} depth={depth[j,i]:.2f}m -> 0 (terra)")

# Setar como terra (depth=0)
depth[prob] = 0.0
print(f"\n{prob.sum()} pontos setados como terra")

# Salvar topog corrigido
with nc.Dataset('./INPUT/topog.nc', 'w',
                format='NETCDF3_64BIT_OFFSET') as ds:
    ds.createDimension('ny', depth.shape[0])
    ds.createDimension('nx', depth.shape[1])
    d          = ds.createVariable('depth', 'f4', ('ny','nx'),
                                   fill_value=0.0)
    d[:]       = depth
    d.units    = 'm'
    d.long_name = 'Ocean bottom depth'
    ds.history = f'{prob.sum()} pontos rasos antarticos removidos'

# Verificar
ds2   = nc.Dataset('./INPUT/topog.nc')
dep2  = np.array(ds2['depth'][:])
ds2.close()
print(f"\nVerificacao:")
print(f"  depth min={dep2[dep2>0].min():.2f}m (deve ser > 10m)")
print(f"  pontos ainda problematicos: "
      f"{((lat_T<-80) & (dep2>0) & (dep2<=10.5)).sum()}")
EOF



###############################################################
# plot TOPO (batimetria)
###############################################################
python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

kh_file  = nc.Dataset('./INPUT/ocean_hgrid.nc')
kh       = kh_file['area'][1::2, 1::2]
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
axes[0].set_title('hgrid_area')
plt.colorbar(axes[0].images[0], ax=axes[0], label='m2')
axes[1].imshow(depth, origin='lower', aspect='auto', cmap='Blues')
axes[1].set_title('topog depth')
plt.colorbar(axes[1].images[0], ax=axes[1], label='m')
plt.tight_layout()
plt.savefig('check_TOPO.png', dpi=100)
print("\n  Figura salva: check_TOPO.png")
EOF



echo ####################################
echo "end make_topo topog.nc"
echo ####################################
echo "--------------------------0---------------------------"
echo "     Check Topografia ocean  topog.nc                   "
echo "--------------------------0---------------------------"



python3 << 'EOF1'
import netCDF4 as nc
import numpy as np

hgrid = nc.Dataset('./INPUT/ocean_hgrid.nc')
print("Dimensoes do hgrid:")
for d in hgrid.dimensions:
    print(f"  {d} = {len(hgrid.dimensions[d])}")

print("\nVariaveis e shapes:")
for v in ['x','y','dx','dy']:
    print(f"  {v}: {hgrid[v].shape}")

x = hgrid['x'][:]
y = hgrid['y'][:]
print(f"\nx shape: {x.shape}")
print(f"x min/max: {x.min():.2f} / {x.max():.2f}")
print(f"y min/max: {y.min():.2f} / {y.max():.2f}")
hgrid.close()
EOF1

rm ./INPUT/topog_bruto.nc.bak ./INPUT/topog_bruto.nc
ncdump -h ./INPUT/topog.nc


python3 << 'EOF'
import numpy as np
import netCDF4 as nc

# Ler topog
ds_in  = nc.Dataset('./INPUT/topog.nc')
depth  = np.array(ds_in['depth'][:])
ds_in.close()

# h2 = rugosidade topografica subgrade (tipicamente ~(30m)^2 sobre fundo rugoso)
# h2 = profundidade ao quadrado normalizada (roughness amplitude)
# Para configuracao simples: h2 = (depth * 0.01)^2  (1% da profundidade)
# Pontos de terra: h2 = 0
h2 = np.where(depth > 0, (depth * 0.01)**2, 0.0)

print(f"depth: min={depth[depth>0].min():.1f}  max={depth.max():.1f}")
print(f"h2:    min={h2[h2>0].min():.4f}  max={h2.max():.2f}  (m^2)")

# Adicionar h2 ao topog.nc
import shutil
shutil.copy('./INPUT/topog.nc', './INPUT/topog.nc.bak_h2')

with nc.Dataset('./INPUT/topog.nc', 'a') as ds:
    if 'h2' not in ds.variables:
        vh2          = ds.createVariable('h2', 'f4', ('ny','nx'))
        vh2[:]       = h2.astype(np.float32)
        vh2.units    = 'm2'
        vh2.long_name = 'Squared sub-grid topographic roughness amplitude'
        print("Variavel h2 adicionada ao topog.nc")
    else:
        ds['h2'][:]  = h2.astype(np.float32)
        print("Variavel h2 atualizada em topog.nc")

# Verificar
ds2 = nc.Dataset('./INPUT/topog.nc')
print(f"\ntopog.nc vars: {list(ds2.variables.keys())}")
h2v = np.array(ds2['h2'][:])
print(f"h2: min={h2v[h2v>0].min():.4f}  max={h2v.max():.2f}")
ds2.close()
EOF
rm ./INPUT/topog.nc.bak_h2
ncdump -h ./INPUT/topog.nc

