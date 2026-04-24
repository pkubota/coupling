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
module load cdo

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
# PASSO 1 - Visao Geral e Dependencias
################################################################################
#
#Projeto MONAN - INPE/CPTEC
#Referencia: configuracao OM_1deg do GFDL
#Conteudo: topog.nc            - ocean_hgrid.nc          - ocean_mosaic.nc   - grid_spec.nc       -
#          vgrid_75_2m.nc      - layer_coord.nc          - hycom1_75_800m.nc - tidal_amplitude.nc -
#          KH_background_2d.nc - seawifs_1998-2006_smoothed_2X.nc -
#          MOM_channels         - atmos/land mosaic exchange grids

#   1. Visao Geral e Dependencias

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

## Ordem correta de geracao de TODOS os arquivos
#```
#1. make_hgrid          -> ocean_hgrid.nc
#2. python vgrid.py     -> vgrid_75_2m.nc  (necessario para make_topog)
#3. make_solo_mosaic    -> ocean_mosaic.nc
#4. make_topog          -> topog.nc
#5. make_solo_mosaic    -> atmos_mosaic.nc
#6. make_coupler_mosaic -> grid_spec.nc + exchange grids
#7. python restante     -> KH, tidal, seawifs, masks, channels


################################################################################
# PASSO 3 - Grades Verticais:  : vgrid_75_2m.nc, layer_coord.nc,hycom1_75_800m.nc
################################################################################
echo "---------------make_vgrid -----------------" 
echo "---------------vgrid.nc-----------------" 

# Comando para OM_1deg com 75 camadas, dz_min=2m
# Estrategia: 3 regioes
#   0-10m:    camadas finas (dz=2m)  ->  10 supergrid points
#   10-200m:  transicao              ->  20 supergrid points  
#   200-6500m: camadas grossas       -> 120 supergrid points
# Total supergrid = 150 = 2*75
rm INPUT/vgrid_75_2m.nc
${path_fre}/make_vgrid \
    --nbnds 4 \
    --bnds 0,10,200,6500 \
    --nz  10,20,120 \
    --center c_cell \
    --grid_name vgrid_75_2m \
    2>&1
mv vgrid_75_2m.nc  INPUT/vgrid_75_2m.nc
echo "---------------post-processing  vgrid.nc-----------------" 
echo "O make_vgrid do FRE-NCtools gera apenas zeta(nzv).     "
echo "O arquivo 'oficial' com zw, zt, dz e um formato        "
echo "estendido que precisa ser gerado por pos-processamento."

cat<<'EOF'>create_vgrid_post.py
import netCDF4 as nc
import numpy as np

# le o arquivo do make_vgrid
src = nc.Dataset("INPUT/vgrid_75_2m.nc")
zeta = src["zeta"][:]   # (151,) = supergrid vertical
src.close()

# calcula as variaveis derivadas
nzv = len(zeta)          # 151
nz  = (nzv - 1) // 2    # 75  (camadas)

# interfaces: indices pares do supergrid
zw = zeta[0::2]          # (76,) = nz_iface

# centros: indices impares
zt = zeta[1::2]          # (75,) = nz

# espessuras
dz = np.diff(zw)         # (75,)

print(f"nzv={nzv}  nz={nz}  nz_iface={len(zw)}")
print(f"zw: {zw[:3]} ... {zw[-3:]}")
print(f"zt: {zt[:3]} ... {zt[-3:]}")
print(f"dz: {dz[:3]} ... {dz[-3:]}")

# escreve o arquivo no formato esperado pelo MOM6
dst = nc.Dataset("INPUT/vgrid_75_2m_mom6.nc", "w", format="NETCDF4_CLASSIC")
dst.createDimension("nzv",      nzv)
dst.createDimension("nz",       nz)
dst.createDimension("nz_iface", nz + 1)

v = dst.createVariable("zeta", "f8", ("nzv",))
v.units    = "m"
v.positive = "down"
v.long_name = "Vertical grid supergrid positions"
v[:] = zeta

v = dst.createVariable("zw", "f8", ("nz_iface",))
v.units    = "m"
v.positive = "down"
v.long_name = "Depth of layer interfaces"
v[:] = zw

v = dst.createVariable("zt", "f8", ("nz",))
v.units    = "m"
v.positive = "down"
v.long_name = "Depth of layer centers"
v[:] = zt

v = dst.createVariable("dz", "f8", ("nz",))
v.units    = "m"
v.long_name = "Layer thickness"
v[:] = dz

dst.nk    = nz
dst.title = f"Vertical grid for MOM6 OM_1deg {nz} layers"
dst.close()
print("OK ? vgrid_75_2m_mom6.nc criado")
EOF

python3 create_vgrid_post.py
