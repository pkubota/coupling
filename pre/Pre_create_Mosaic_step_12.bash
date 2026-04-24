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
#cd INPUT
################################################################################
# PASSO 4 - mosaico ocean_mosaic
################################################################################
echo "--------------------------0---------------------------"
echo "     2. Mosaico do oceano                             "
echo "--------------------------0---------------------------"
# ocean_mosaic.nc + atmos_mosaic.nc
rm  INPUT/ocean_mosaic.nc
ls ./INPUT/ocean_hgrid.nc
ncdump -h ./INPUT/ocean_hgrid.nc 
# -- 2. Mosaico do oceano 
${path_fre}/make_solo_mosaic \
  --num_tiles 1 \
  --dir ./INPUT/ \
  --mosaic_name ocean_mosaic \
  --tile_file ocean_hgrid.nc \
  --periodx 360
echo "--------------------------0---------------------------"
echo "    End Mosaico do oceano                             "
echo "--------------------------0---------------------------"
  
mv ocean_mosaic.nc  INPUT/ocean_mosaic.nc
# deve mostrar: ny=180, nx=360, variavel depth
################################################################################
# PASSO 6 - Mosaico da atmosfera (forcamento JRA = mesma grade) -----
################################################################################

# -- 2.3. Mosaico da atmosfera (forcamento JRA = mesma grade) -----
#
cp ./INPUT/ocean_hgrid.nc ./INPUT/atmos_hgrid.nc

rm ./INPUT/atmos_mosaic.nc
${path_fre}/make_solo_mosaic \
  --num_tiles 1 \
  --dir ./INPUT/ \
  --mosaic_name atmos_mosaic \
  --tile_file atmos_hgrid.nc  \
  --periodx 360
mv atmos_mosaic.nc ./INPUT/atmos_mosaic.nc
echo "--------------------------0---------------------------"
echo "    End Mosaico do Atmos                              "
echo "--------------------------0---------------------------"

# -- 2.3. Mosaico da atmosfera (forcamento JRA = mesma grade) -----
#
cp ./INPUT/atmos_hgrid.nc ./INPUT/land_hgrid.nc
rm ./INPUT/land_mosaic.nc
${path_fre}/make_solo_mosaic \
  --num_tiles 1 \
  --dir ./INPUT/ \
  --mosaic_name land_mosaic \
  --tile_file land_hgrid.nc  \
  --periodx 360
echo "--------------------------0---------------------------"
echo "    End Mosaico do Land                               "
echo "--------------------------0---------------------------"
mv land_mosaic.nc ./INPUT/land_mosaic.nc
#atmos_mosaic_tile1Xocean_mosaic_tile1.nc -> .datasets/OM_1deg/INPUT/atmos_mosaic_tile1Xocean_mosaic_tile1.nc
#atmos_mosaic_tile1Xland_mosaic_tile1.nc ->  .datasets/OM_1deg/INPUT/atmos_mosaic_tile1Xland_mosaic_tile1.nc
################################################################################
# PASSO 6 - Coupler mosaic (gera grid_spec + exchange grids)
################################################################################
cp ./INPUT/atmos_hgrid.nc  .
cp ./INPUT/ocean_hgrid.nc  .
cp ./INPUT/land_hgrid.nc   .
cp ./INPUT/topog.nc   .
cp ./INPUT/atmos_mosaic.nc .
cp ./INPUT/land_mosaic.nc .
cp ./INPUT/ocean_mosaic.nc .
rm  ./INPUT/ocean_mask.nc
rm  ./INPUT/land_mask.nc
rm  ./INPUT/grid_spec.nc
rm  ./INPUT/atmos_mosaic_tile1Xocean_mosaic_tile1.nc
rm  ./INPUT/atmos_mosaic_tile1Xatmos_mosaic_tile1.nc

# -- 2.4. Coupler mosaic (gera grid_spec + exchange grids)
${path_fre}/make_coupler_mosaic \
  --atmos_mosaic atmos_mosaic.nc \
  --land_mosaic  atmos_mosaic.nc \
  --ocean_mosaic ocean_mosaic.nc \
  --ocean_topog  topog.nc \
  --mosaic_name  grid_spec \
  --check \
  --verbose
mv ocean_mask.nc   ./INPUT/
mv land_mask.nc    ./INPUT/
mv grid_spec.nc    ./INPUT/
mv atmos_mosaic_tile1Xocean_mosaic_tile1.nc  ./INPUT/
mv atmos_mosaic_tile1Xatmos_mosaic_tile1.nc  ./INPUT/
rm ./atmos_hgrid.nc
rm ./ocean_hgrid.nc
rm ./land_hgrid.nc 
rm ./topog.nc
rm ./atmos_mosaic.nc
rm ./land_mosaic.nc
rm ./ocean_mosaic.nc

echo "--------------------------0---------------------------"
echo "    End Mosaico do Coupler                               "
echo "--------------------------0---------------------------"
## O que isso gera (para o MOM6+SIS2)
#```
#grid_spec.nc                               <- lido pelo coupler FMS
#ocean_mosaic.nc                            <- lido pelo MOM6
#atmos_mosaic.nc                            <- lido pelo data_override (JRA)
#atmos_mosaic_tile1Xland_mosaic_tile1.nc    <- troca atmos<->terra
#atmos_mosaic_tile1Xocean_mosaic_tile1.nc   <- troca atmos<->oceano/gelo
#land_mosaic_tile1Xocean_mosaic_tile1.nc    <- fronteira terra<->oceano
# O problema eh que todos os *_dir estao com "./" 
# - isso significa que o FMS/MOM6 
#   vai procurar os arquivos referenciados no diretorio de trabalho (rundir),
#   nao em INPUT/.
#
# Opcao 1 - Recriar o grid_spec.nc com os diretorios corretos
# A forma mais limpa. Adicionar --dir INPUT/ (ou o path que preferir) 
# no make_coupler_mosaic nao eh suportado diretamente, mas você pode 
# patcher o arquivo NetCDF via Python:
################################################################################

###########################################################################################
#
#Você precisa corrigir 4 arquivos no total:
#
#Arquivo             Campo a corrigir                     Valor atual      Valor necessário
#grid_spec.nc        atm_mosaic_dir, lnd_mosaic_dir,         ./              INPUT/
#                    ocn_mosaic_dir, ocn_topog_dir        
#ocean_mosaic.nc     gridlocation                            ./              INPUT/
#atmos_mosaic.nc     gridlocation                            ./              INPUT/
#land_mosaic.nc      gridlocation                            ./              INPUT/
###########################################################################################
#
python3 << 'EOF'
import netCDF4 as nc
import numpy as np

NEW_DIR = "INPUT/"

def set_string_var(f, var_name, new_value):
    var = f.variables[var_name]
    str_len = var.shape[-1]
    padded = new_value.ljust(str_len)
    char_array = np.array(list(padded), dtype="S1")
    if var.ndim == 1:
        var[:] = char_array
    else:
        for i in range(var.shape[0]):
            var[i, :] = char_array

# grid_spec.nc
with nc.Dataset("INPUT/grid_spec.nc", "r+") as f:
    for v in ["atm_mosaic_dir", "lnd_mosaic_dir", "ocn_mosaic_dir", "ocn_topog_dir"]:
        set_string_var(f, v, NEW_DIR)
    print("grid_spec.nc: OK")
EOF
#Nota sobre o atmos_mosaic.nc
#No seu grid_spec.nc, tanto atm_mosaic_file quanto lnd_mosaic_file apontam para
# atmos_mosaic.nc - ou seja, não existe land_mosaic.nc separado,
# o que é normal para configurações ATM-only ou quando terra e atmosfera
#  compartilham a mesma grade. O script ja trata isso.
################################################################################
