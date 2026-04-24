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
#1. Grade horizontal
rm  INPUT/ocean_hgrid.nc
${path_fre}/make_hgrid \
     --grid_type tripolar_grid \
     --nxbnds 2 \
     --nybnds 7 \
     --xbnds -300,60 \
     --ybnds -90,-65,-30,0,30,60,90 \
     --nlon 360 \
     --nlat 14,72,56,56,56,62 \
     --grid_name ocean_hgrid \
     --center c_cell 
 mv ocean_hgrid.nc INPUT/ocean_hgrid.nc
              #Os --nlat somam 104+48+40+40+48+136 = 416 ? supergrid ny=417 
              #Os --nlon  360*2 soma 720
  echo "---------------check ocean_hgrid.nc-----------------" 
  
  
python3 << 'EOF1'
import netCDF4 as nc
import numpy as np

hgrid = nc.Dataset('INPUT/ocean_hgrid.nc')
x = hgrid['x'][:]
y = hgrid['y'][:]
hgrid.close()

# supergrid (317,361): pontos de interface sao os pares, centros sao impares
# Testar os dois casos:
print("Caso A - indices impares [1::2, 1::2]:")
xA = x[1::2, 1::2]
print(f"  shape: {xA.shape}")   # deve ser (158, 180)

print("Caso B - indices pares [0::2, 0::2]:")
xB = x[0::2, 0::2]
print(f"  shape: {xB.shape}")   # deve ser (159, 181)

# O correto para grade T de ny=316, nx=360:
# ny=316 celulas -> 317 interfaces -> pontos T nao batem com [1::2]
# Pontos T estao em todos os pontos da grade ny/nx (nao supergrid)
print("\nCaso C - grade completa ny x nx (sem supergrid):")
print(f"  x[0::1, 0::1] shape: {x[0::1,0::1].shape}")

# Para ny=316, nx=360: pegar linhas 0..315, colunas 0..359
print("Caso D - ny=316, nx=360 direto:")
xD = x[:316, :360]
print(f"  shape: {xD.shape}")
print(f"  lon min/max: {xD.min():.2f} / {xD.max():.2f}")
print(f"  lat min/max: {y[:316,:360].min():.2f} / {y[:316,:360].max():.2f}")
EOF1


echo "---------------plot ocean_hgrid.nc-----------------" 


python3 << 'EOF'
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ds     = nc.Dataset('INPUT/ocean_hgrid.nc')
area   = np.array(ds['area'][:])
x      = np.array(ds['x'][:])
y      = np.array(ds['y'][:])
ds.close()

# Pontos T da grade do modelo
area_T = area[1::2, 1::2]   # (158, 180)
lon_T  = x[1::2, 1::2]
lat_T  = y[1::2, 1::2]

# Normalizar lon para -180:180
lon_T  = ((lon_T + 180) % 360) - 180

area_km2 = area_T / 1e6   #  m*m-> km*km

print(f"area_T shape: {area_T.shape}")
print(f"lon:  {lon_T.min():.1f} a {lon_T.max():.1f}")
print(f"lat:  {lat_T.min():.1f} a {lat_T.max():.1f}")
print(f"area: {area_km2.min():.0f} a {area_km2.max():.0f} km2")

# -- Figura com 3 paineis -
fig = plt.figure(figsize=(18, 14))
fig.suptitle('INPUT/ocean_hgrid.nc - Campo de Area (pontos T)', fontsize=14, fontweight='bold')

# - Painel 1: Mapa global -
ax1 = fig.add_subplot(2, 2, (1,2),
      projection=ccrs.Robinson(central_longitude=0))
ax1.set_global()
ax1.add_feature(cfeature.LAND, color='lightgray', zorder=2)
ax1.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=3)
ax1.gridlines(linewidth=0.3, color='gray', alpha=0.5)

sc = ax1.pcolormesh(lon_T, lat_T, area_km2,
                    transform=ccrs.PlateCarree(),
                    cmap='plasma', vmin=0, vmax=area_km2.max())
cb = plt.colorbar(sc, ax=ax1, orientation='horizontal',
                  pad=0.04, shrink=0.8, label='Area (km2)')
ax1.set_title('Area das celulas T - Global (projecao Robinson)', fontsize=11)

# -- Painel 2: Perfil latitudinal -
ax2 = fig.add_subplot(2, 2, 3)
lat_1d   = lat_T[:, 0]
area_min = area_km2.min(axis=1)
area_max = area_km2.max(axis=1)
area_med = area_km2.mean(axis=1)

ax2.fill_betweenx(lat_1d, area_min, area_max,
                  alpha=0.3, color='steelblue', label='min-max')
ax2.plot(area_med, lat_1d, color='steelblue', lw=2, label='media zonal')
ax2.set_xlabel('Area (km2)', fontsize=10)
ax2.set_ylabel('Latitude (graus)', fontsize=10)
ax2.set_title('Variacao latitudinal da area', fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.axhline(0, color='k', lw=0.8, ls='--')

# Anotar valores tipicos
for lat_ref in [-60, -30, 0, 30, 60]:
    j   = np.argmin(np.abs(lat_1d - lat_ref))
    val = area_med[j]
    ax2.annotate(f'{val:.0f} km2', xy=(val, lat_ref),
                 xytext=(val + area_med.max()*0.05, lat_ref),
                 fontsize=7, va='center', color='navy')

# - Painel 3: Histograma -
ax3 = fig.add_subplot(2, 2, 4)
ax3.hist(area_km2.flatten(), bins=60, color='steelblue',
         edgecolor='white', linewidth=0.3)
ax3.axvline(area_km2.mean(), color='red', lw=1.5,
            label=f'm2dia = {area_km2.mean():.0f} km2')
ax3.axvline(np.median(area_km2), color='orange', lw=1.5,
            label=f'mediana = {np.median(area_km2):.0f} km2')
ax3.set_xlabel('Area (km2)', fontsize=10)
ax3.set_ylabel('Frequencia', fontsize=10)
ax3.set_title('Distribuicao das areas das celulas T', fontsize=11)
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)

# Estatisticas no painel (89)
stats = (f"min  = {area_km2.min():.0f} km2\n"
         f"max  = {area_km2.max():.0f} km2\n"
         f"mean = {area_km2.mean():.0f} km2\n"
         f"soma = {area_km2.sum():.2e} km2")
ax3.text(0.97, 0.97, stats, transform=ax3.transAxes,
         fontsize=8, va='top', ha='right',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('area_hgrid.png', dpi=150, bbox_inches='tight')
print("Figura salva: area_hgrid.png")
EOF
echo "---------------ocean_hgrid.nc-----------------" 


