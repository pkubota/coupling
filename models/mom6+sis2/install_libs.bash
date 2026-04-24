#!/bin/bash +x
################################################################################
# OPCAO B - Guia de CompilaCAo MOM6-examples com mkmf no Cray/GNU
#
# RepositOrio: https://github.com/NOAA-GFDL/MOM6-examples
#
# Este guia segue o workflow oficial GFDL (mkmf) adaptado para
# o ambiente Cray XC/EX com PrgEnv-gnu.
#
# Estrutura final apOs este processo:
#   MOM6-examples/
#   |__ build/
#       |__ gnu/
#           |__ shared/repro/           <- libfms.a
#           |__ ice_ocean_SIS2/repro/   <- executavel MOM6 standalone (teste)
#           |__ nuopc_cap/repro/        <- libmom6_nuopc.a (cap NUOPC para MPAS)
################################################################################
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
module load autoconf/2.72
module load libfabric/1.22.0
module load cray-pals/1.6.1
module list
export PATH=/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/bin:$PATH
export LD_LIBRARY_PATH=/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/libssl/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/lib:\
/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/lib:$LD_LIBRARY_PATH
#

# Verifique os modulos disponiveis com:
#   module avail esmf
#   module avail cray-netcdf
#
# IMPORTANTE no Cray: o compilador deve ser sempre "ftn" (wrapper), nao gfortran

################################################################################
# PASSO 1 - Clonar o repositOrio com submOdulos
################################################################################
#
# cd /p/projetos/monan_atm/paulo.kubota/coupler/compiler
# git clone --recursive https://github.com/NOAA-GFDL/MOM6-examples.git
 cd MOM6-examples
 mom6_example_dir=`pwd`
#
# Verificar submOdulos (FMS, MOM6, SIS2, mkmf devem aparecer):
# git submodule status
#
# Se algum submOdulo estiver vazio:
# git submodule update --init --recursive
cd ../
################################################################################
# PASSO 2 - Criar template mkmf para Cray/GNU
################################################################################
#
# O mkmf usa templates .mk para configurar compiladores por sistema.
# Copie o template GNU existente e adapte para o Cray:
 cd ${mom6_example_dir}
# cp src/mkmf/templates/ncrc5-gcc.mk \
#    src/mkmf/templates/cray-gnu.mk
 cd ..
#
# Edite src/mkmf/templates/cray-gnu.mk conforme abaixo:

# ---- CONTEUDO DO TEMPLATE cray-gnu.mk ----
# (salve este conteudo em src/mkmf/templates/cray-gnu.mk)
#
# # Template para Cray XC/EX com PrgEnv-gnu
# # Uso: mkmf -t cray-gnu.mk ...
#
export FC=ftn
export CC=cc
export LD=ftn
export NetCDF_ROOT=$(nc-config --prefix)
# FC = ftn
# CC = cc
# LD = ftn
#
# # Flags de compilaCAo reprodutivel (repro)
# FFLAGS_REPRO = -O2 -fbacktrace -ffree-line-length-none \
#                -fno-second-underscore
#
# # Flags de debug
# FFLAGS_DEBUG = -O0 -g -fbacktrace -fcheck=all -fcheck=bounds \
#                -ffpe-trap=invalid,zero,overflow -Wall
#
# # Flags OpenMP (se necessario)
# FFLAGS_OMP   = -fopenmp
#
# # O wrapper ftn do Cray injeta automaticamente:
# #   - MPI (cray-mpich)
# #   - NetCDF (cray-netcdf)
# #   - HDF5 (cray-hdf5)
# # Por isso LIBS pode ficar vazio ou com apenas extras
# LIBS =
#
# # Definicoes de pre-processamento
# CPPDEFS = -Duse_libMPI -Duse_netCDF -DSPMD
#
# # Modo padrao
# ifdef REPRO
#   FFLAGS = $(FFLAGS_REPRO) $(CPPDEFS)
# else ifdef DEBUG
#   FFLAGS = $(FFLAGS_DEBUG) $(CPPDEFS)
# else
#   FFLAGS = $(FFLAGS_REPRO) $(CPPDEFS)
# endif
# ---- FIM DO TEMPLATE ----
################################################################################
# PASSO 3 - Compilar FMS (biblioteca compartilhada base)
################################################################################
#
# O FMS e a infraestrutura base do MOM6. Deve ser compilado primeiro.
#
cd ${mom6_example_dir}
# Crie o diretOrio de build e compile:
mkdir  ${mom6_example_dir}/build

 TEMPLATE_MK="/p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/MOM6-examples/src/mkmf/templates/cray-gnu.mk"
 MKMF="$(pwd)/src/mkmf/bin/mkmf"
 LIST_PATHS="$(pwd)/src/mkmf/bin/list_paths"
#
 mkdir -p  ${mom6_example_dir}/build/gnu/shared/repro
 cd ${mom6_example_dir}/build/gnu/shared/repro
#
# # Gera path_names com todos os arquivos fonte do FMS
 rm -f path_names
 ${LIST_PATHS} -l ../../../../src/FMS
#
# # Gera o Makefile via mkmf
 ${MKMF} \
     -t ${TEMPLATE_MK} \
     -p libfms.a \
     -c "-Duse_libMPI -Duse_netCDF -DSPMD" \
     path_names
#
# # Compila
 make NETCDF=3 REPRO=1 libfms.a -j 8 2>&1 | tee make_fms.log

#
# # Verificar:
 ls -la libfms.a   # deve existir e ter tamanho > 0

 cd ${mom6_example_dir}

test_ice_ocean_SIS2=1

if [ $test_ice_ocean_SIS2 -eq 1 ]; then
 echo 'ice_ocean_SIS2'
 # Commands to execute if the condition is true
################################################################################
# PASSO 4 - Compilar MOM6 standalone (ice_ocean_SIS2) - para teste inicial
################################################################################
#
# Isso verifica se o MOM6 compila corretamente antes de tentar o cap NUOPC.
#
 mkdir -p ${mom6_example_dir}/build/gnu/ice_ocean_SIS2/repro
 cd ${mom6_example_dir}/build/gnu/ice_ocean_SIS2/repro
#
 rm -f path_names

 ${LIST_PATHS} -l              \
      ./ \
     ../../../../src/MOM6/config_src/{infra/FMS2,memory/dynamic_symmetric,drivers/FMS_cap,external} \
     ../../../../src/SIS2/config_src/dynamic_symmetric             \
     ../../../../src/MOM6/src/{*,*/*}/                              \
     ../../../../src/atmos_null                                    \
     ../../../../src/land_null                                     \
     ../../../../src/coupler                                       \
     ../../../../src/{ice_param,icebergs/src,SIS2,FMS/coupler,FMS/include}
#
 ${MKMF} \
     -t ${TEMPLATE_MK} \
     -o "-I../../shared/repro  -I../../ice_ocean_SIS2/repro " \
     -p MOM6 \
     -l "-L../../shared/repro -lfms" \
     -c "-Duse_libMPI -Duse_netCDF -DSPMD -Duse_AM3_physics -D_USE_LEGACY_LAND_ " \
     path_names
#
 cd ${mom6_example_dir}/build/gnu/ice_ocean_SIS2/repro
 make REPRO=1 MOM6 -j 8 2>&1 | tee make_mom6_ice_ocean_SIS2.log
echo "[2/3] MOM6 ice_ocean_SIS2: OK"
#

 ls -la MOM6   # executavel standalone

 cd ${mom6_example_dir}
else
  # Commands to execute if the condition is false

################################################################################
# PASSO 5 - Compilar MOM6 standalone (ocean_only) - para teste inicial
################################################################################
 echo 'ocean_only'
#
# Isso verifica se o MOM6 compila corretamente antes de tentar o cap NUOPC.
#
 mkdir -p ${mom6_example_dir}/build/gnu/ocean_only/repro
 cd ${mom6_example_dir}/build/gnu/ocean_only/repro
#
 rm -f path_names
 ${LIST_PATHS} -l \
     ../../../../src/MOM6/config_src/infra/FMS2 \
     ../../../../src/MOM6/config_src/memory/dynamic_symmetric \
     ../../../../src/MOM6/config_src/drivers/solo_driver \
     ../../../../src/MOM6/config_src/external \
     ../../../../src/MOM6/src/{*,*/*}
#
 ${MKMF} \
     -t ${TEMPLATE_MK} \
     -o "-I../../shared/repro" \
     -p MOM6 \
     -l "-L../../shared/repro -lfms" \
     -c "-Duse_libMPI -Duse_netCDF -DSPMD" \
     path_names
#
 cd ${mom6_example_dir}/build/gnu/ocean_only/repro
 make REPRO=1 MOM6 -j 8 2>&1 | tee make_mom6_ocean_only.log
echo "[2/3] MOM6 ocean_only: OK"
#
 ls -la MOM6   # executavel standalone

 cd ${mom6_example_dir}
fi



 cd ${mom6_example_dir}

################################################################################
# PASSO 6 - Compilar MOM6 com o cap NUOPC (para acoplamento com MPAS)
################################################################################
#
# Esta e a etapa principal para o acoplamento MPAS-MOM6.
# Usa config_src/drivers/nuopc_cap em vez de solo_driver.
# Inclui ESMF no build.
#
 ESMF_DIR="/p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf"
 ESMF_INC="${ESMF_DIR}/include"
 ESMF_MOD="${ESMF_DIR}/mod/modg/Linux"
 ESMF_LIB="${ESMF_DIR}/lib/libg/Linux"
 #ESMF_FLAGS="$(${ESMF_DIR}/bin/bing/Linux)"
 ESMF_FLAGS="$(${ESMF_DIR}/bin/esmf_libs 2>/dev/null || echo "-L${ESMF_DIR}/lib -lesmf")"

#
 mkdir -p ${mom6_example_dir}/build/gnu/nuopc_cap/repro
 cd ${mom6_example_dir}/build/gnu/nuopc_cap/repro
#
 rm -f path_names
#
# # Inclui: FMS (infraestrutura moderna), memOria dinamica,
#     ../../../../src/MOM6/config_src/drivers/nuopc_cap \
# #         driver nuopc_cap, cOdigo externo e todo o src MOM6
# 
# Inclui tambem o cap adaptado (MOM6_cap.F90 da OpCAo A)
# echo "../../../../../src/time_utils.F90" >> path_names
# echo "../../../../../src/MOM_ocean_model_nuopc.F90" >> path_names
# echo "../../../../../src/NUOPC_Driver.F90" >> path_names
# echo "../../../../../src/MOM6_cap.F90" >> path_names
#
 ${LIST_PATHS} -l \
     ../../../../src/MOM6/config_src/infra/FMS2 \
     ../../../../src/MOM6/config_src/memory/dynamic_symmetric \
     ../../../../src/MOM6/config_src/external \
     ../../../../../../../caps/ocean/*.F90               \
     ../../../../src/MOM6/src/{*,*/*}
#
 ${MKMF} \
     -t ${TEMPLATE_MK} \
     -o "-I../../shared/repro -I${ESMF_INC} -I${ESMF_MOD} -I${mom6_example_dir}/src/MOM6/src/framework" \
     -p libmom6_nuopc.a \
     -l "-L../../shared/repro -lfms ${ESMF_FLAGS}" \
     -c "-Duse_libMPI -Duse_netCDF -DSPMD -DUSE_ESMF_NUOPC -fcheck=all" \
     path_names
#
 make REPRO=1 libmom6_nuopc.a -j 8 2>&1 | tee make_mom6_nuopc.log
#
 ls -la libmom6_nuopc.a   # biblioteca do cap NUOPC

cd ${mom6_example_dir}
echo ""
echo "============================================"
echo " Build completo!"
echo ""
echo " Resultados:"
echo "   FMS    : build/gnu/shared/repro/libfms.a"
echo "   MOM6   : build/gnu/ocean_only/repro/MOM6"
echo "   MOM62  : build/gnu/ice_ocean_SIS2/repro/MOM6"
echo "   NUOPC  : build/gnu/nuopc_cap/repro/libmom6_nuopc.a"
echo ""
echo " PrOximo passo: linkar libmom6_nuopc.a ao driver MPAS-MOM6"
echo "============================================"
cp build/gnu/nuopc_cap/repro/*.a      ../lib/nuopc/ 
cp build/gnu/ice_ocean_SIS2/repro/*.o ../lib/mom6/ 
cp build/gnu/shared/repro/libfms.a    ../lib/fms/ 

cp build/gnu/nuopc_cap/repro/*.mod       ../mod/nuopc/ 
cp build/gnu/ice_ocean_SIS2/repro/*.mod  ../mod/mom6/ 
cp build/gnu/shared/repro/*.mod          ../mod/fms/ 
