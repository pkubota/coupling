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
module load cray-python/3.11.7
#   module load esmf/8.4.0         # ESMF (ajuste versao disponivel no seu Cray)
module load cdo2/2.5.4
module load ncl/6.2.2
module load nco/5.3.6
module list
###########################################################################################
export FC=ftn
export CC=cc
###########################################################################################
export PATH=/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/automake/automake-1.16.5/bin:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/bin:\
/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/bin:$PATH
export LD_LIBRARY_PATH=/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/libssl/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/automake/automake-1.16.5/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/lib:\
/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/lib:$LD_LIBRARY_PATH
path_fre=/p/projetos/monan_atm/paulo.kubota/coupler/lib_couplers/fms/FRE-NCtools/build/bin
###########################################################################################
#
# Otima pergunta -> esse e um ponto critico no acoplamento MPAS-A -> MOM6+SIS2 
#(via NUOPC/ESMF). Vou detalhar o fluxo completo para gerar o ocean_static.nc 
#e derivar as mascaras corretamente.
#
#Fluxo: topog.nc / ocean_mask.nc -> ocean_static.nc com mascaras para MOM6+SIS2
###########################################################################################
#1. Entenda o que o MOM6 espera
###########################################################################################
#O arquivo ocean_static.nc (tambem chamado de ocean_hgrid.nc + ocean_topog.nc em algumas configs) precisa conter:
#
#Variavel                       Significado
#wet                            1=oceano, 0=terra (em pontos T da grade B/C)
#mask                           idem, às vezes usado como alias
#tmask (ou mask_rho)            mascara nos pontos tracer (T) -> necessaria para SIS2 e diag
#depth_ocean                    batimetria positiva para baixo (metros)
#geolat, geolon                 coordenadas geograficas
#
###########################################################################################
#2. Ponto de partida: seus arquivos
###########################################################################################
#
#topog.nc         -> contem depth_ocean (ou "depth", "topo", "h")
#ocean_mask.nc    -> contem mascara oceanica (0/1)
#land_mask.nc     -> contem mascara terrestre (0/1, inverso)
#
###########################################################################################
#3. Script Python completo: gerar ocean_static.nc
###########################################################################################
rm -f ./INPUT/ocean_static.nc
python3 << 'EOF'
#######!/usr/bin/env python3
#"""
#Gera ocean_static.nc para acoplamento MPAS-MOM6+SIS2
#Deriva: wet -> mask -> tmask a partir de topog.nc / ocean_mask.nc / land_mask.nc
#"""

import numpy as np
import xarray as xr

# -- 1. Carrega arquivos de entrada ------------------------------------------
topo   = xr.open_dataset("./INPUT/topog.nc")
omask  = xr.open_dataset("./INPUT/ocean_mask.nc")
lmask  = xr.open_dataset("./INPUT/land_mask.nc")

# -- 2. Extrai batimetria -----------------------------------------------------
# Tenta nomes comuns; ajuste conforme seu arquivo
for vname in ["depth_ocean", "depth", "topo", "h", "deptho"]:
    if vname in topo:
        depth = topo[vname].values.astype(np.float64)
        break

# -- 3. Deriva 'wet' (mascara oceanica primaria) ------------------------------
# Estrategia 1: a partir da batimetria
wet_from_topo = np.where(depth > 0.0, 1.0, 0.0)

# Estrategia 2: a partir do ocean_mask.nc (mais confiavel se disponivel)
for vname in ["wet", "mask", "ocean_mask", "omask", "frocean"]:
    if vname in omask:
        wet_from_omask = omask[vname].values.astype(np.float64)
        break

# Estrategia 3: complemento do land_mask
for vname in ["land_mask", "lmask","mask", "frlnd", "landfrac"]:
    if vname in lmask:
        wet_from_lmask = 1.0 - lmask[vname].values.astype(np.float64)
        break

# -- 4. Mascara final: intersecao (oceano em TODOS os criterios) --------------
# Use logica AND para consistencia entre os tres:
wet = np.where(
    (wet_from_topo > 0) & (wet_from_omask > 0.5) & (wet_from_lmask > 0.5),
    1.0, 0.0
)

# Garante que celulas com depth <= 0 sejam terra
wet = np.where(depth > 0.0, wet, 0.0)
# Garante que celulas com depth > 0 mas marcadas como terra fiquem como terra
depth = np.where(wet > 0.5, depth, 0.0)

# -- 5. mask = wet (alias exigido por alguns componentes NUOPC) --------------
mask = wet.copy()

# -- 6. tmask: mascara nos pontos T (tracer) para MOM6/SIS2 ------------------
# tmask == wet para grade de pontos T (Arakawa B ou C)
# Para grade B: tmask e definida nos cantos -> aqui mantemos igual ao wet
# Para SIS2: tmask deve ser inteiro (0 ou 1) sem NaN
tmask = wet.astype(np.int32)

# -- 7. Coordenadas ----------------------------------------------------------
# Tenta pegar lat/lon do proprio topog.nc ou omask
lat = None
lon = None
for ds in [topo, omask, lmask]:
    for latname in ["lat", "geolat", "latitude", "y"]:
        if latname in ds:
            lat = ds[latname].values
            break
    for lonname in ["lon", "geolon", "longitude", "x"]:
        if lonname in ds:
            lon = ds[lonname].values
            break
    if lat is not None and lon is not None:
        break

# -- 8. Monta o Dataset ocean_static.nc --------------------------------------
coords = {}
if lat is not None and lon is not None:
    if lat.ndim == 1:
        coords = {"lat": (["yh"], lat), "lon": (["xh"], lon)}
        dims_2d = ["yh", "xh"]
    else:
        coords = {
            "geolat": (["yh", "xh"], lat),
            "geolon": (["yh", "xh"], lon),
        }
        dims_2d = ["yh", "xh"]
else:
    ny, nx = wet.shape
    dims_2d = ["yh", "xh"]

ds_out = xr.Dataset(
    {
        "wet": xr.DataArray(
            wet, dims=dims_2d,
            attrs={"long_name": "ocean fraction at T-cell centers",
                   "units": "none", "valid_range": [0., 1.]}
        ),
        "mask": xr.DataArray(
            mask, dims=dims_2d,
            attrs={"long_name": "ocean mask (1=ocean, 0=land)",
                   "units": "none"}
        ),
        "tmask": xr.DataArray(
            tmask, dims=dims_2d,
            attrs={"long_name": "T-point ocean mask for MOM6/SIS2",
                   "units": "none", "flag_values": "0 1",
                   "flag_meanings": "land ocean"}
        ),
        "depth_ocean": xr.DataArray(
            depth, dims=dims_2d,
            attrs={"long_name": "ocean depth at T-cell centers",
                   "units": "m", "positive": "down"}
        ),
    },
    coords=coords,
    attrs={
        "title": "MOM6+SIS2 ocean static file",
        "history": "Generated by make_ocean_static.py",
        "source_topog": "./INPUT/topog.nc",
        "source_omask": "./INPUT/ocean_mask.nc",
        "source_lmask": "./INPUT/land_mask.nc",
        "Conventions": "CF-1.7",
    }
)

# Adiciona geolat/geolon como variaveis de dados tambem (padrao MOM6)
if lat is not None and lon is not None:
    if lat.ndim == 2:
        ds_out["geolat"] = xr.DataArray(
            lat, dims=dims_2d,
            attrs={"long_name": "geographic latitude of T-cell centers",
                   "units": "degrees_north"}
        )
        ds_out["geolon"] = xr.DataArray(
            lon, dims=dims_2d,
            attrs={"long_name": "geographic longitude of T-cell centers",
                   "units": "degrees_east"}
        )

# -- 9. Salva -----------------------------------------------------------------
encoding = {v: {"zlib": True, "complevel": 4} for v in ds_out.data_vars}
ds_out.to_netcdf("./INPUT/ocean_static.nc", encoding=encoding)
print("- ./INPUT/ocean_static.nc gerado com sucesso!")
print(ds_out)
EOF
###########################################################################################
#4. Verificação com NCO/CDO
###########################################################################################
# Verifica estrutura
ncdump -h ./INPUT/ocean_static.nc

# Confere valores unicos de wet (deve ser so 0 e 1)
ncap2 -s 'print(wet.min(), wet.max())' ./INPUT/ocean_static.nc

# Checa consistencia: depth>0 onde wet=1
ncap2 -s 'check=depth_ocean*wet' ./INPUT/ocean_static.nc ./INPUT/check.nc
ncwa -a yh,xh -y min ./INPUT/check.nc  # minimo deve ser 0

# Plot rapido da tmask
python3 -c "
import xarray as xr, matplotlib.pyplot as plt
ds = xr.open_dataset('./INPUT/ocean_static.nc')
ds['tmask'].plot(cmap='Blues')
plt.title('tmask - MOM6/SIS2')
plt.savefig('tmask_check.png', dpi=150)
"
###########################################################################################
#5. Armadilhas comuns no acoplamento MPAS-MOM6+SIS2
###########################################################################################
#
#Problema                          | Causa                                     |  Solucao 
#tmask tem NaN                     | wet veio de float com fill_value          |  fillna(0) antes de salvar
#Máscara invertida                 | land_mask já é 0=terra                    |  use 1 - land_mask
#Desacordo MPAS<->MOM6             | grades diferentes (Voronoi vs retangular) |  interpole com ESMF_RegridWeightGen
#SIS2 nao le tmask                 | nome de variável diferente                |  confira input.nml -> mask_table
#wet /= 0/1 (valores fracionários) | dado de fração oceânica                   |  aplique threshold: wet = wet > 0.5
###########################################################################################
#6. Referência no MOM_input / SIS_input
###########################################################################################
#fortran! MOM_input
#INPUTDIR = "INPUT/"
#TOPO_FILE = "ocean_static.nc"  ! ou topog.nc separado
#TOPO_VARNAME = "depth_ocean"
#MASK_VARNAME = "wet"           ! MOM6 usa "wet"
#
#
#! SIS_input  
#THICKNESS_UNITS = "m"
#! SIS2 deriva tmask internamente a partir de wet
### Resumo do pipeline
#```
#topog.nc ----------|
#ocean_mask.nc -----|-- wet (AND lógico) --- mask (alias)
#land_mask.nc ------|                    ---- tmask (int32)
#                                        ---- depth_ocean (zerado em terra)
#                                        |
#                              ocean_static.nc  - MOM6 + SIS2
#
