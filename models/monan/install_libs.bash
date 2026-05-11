#!/bin/bash +x
path_monan="/p/projetos/monan_adm/daniel.massaru/Acopladores/scripts_CD-CT/sources/MONAN-Model_feature/monan_2.0.0/src"
rm -rf lib
mkdir -p lib/framework
cp -rf ${path_monan}/framework/*.a     lib/framework/
cp -rf ${path_monan}/framework/*.o     lib/framework/
 core_atmosphere/physics
mkdir -p lib/core_atmosphere
cp -rf ${path_monan}/core_atmosphere/*.a     lib/core_atmosphere/
cp -rf ${path_monan}/core_atmosphere/*.o     lib/core_atmosphere/

mkdir -p  lib/core_atmosphere/physics
cp -rf ${path_monan}/core_atmosphere/physics/*.a     lib/core_atmosphere/physics/
cp -rf ${path_monan}/core_atmosphere/physics/*.o     lib/core_atmosphere/physics/

mkdir -p lib/operators
cp -rf ${path_monan}/operators/*.a     lib/operators/
cp -rf ${path_monan}/operators/*.o     lib/operators/

mkdir -p lib/external/SMIOL
cp -rf ${path_monan}/external/SMIOL/*.a     lib/external/SMIOL/
cp -rf ${path_monan}/external/SMIOL/*.o     lib/external/SMIOL/

###########################################################################
rm -rf include
mkdir -p include/framework
cp -rf ${path_monan}/framework/*.h     include/framework/
cp -rf ${path_monan}/framework/*.inc   include/framework/
cp -rf ${path_monan}/framework/*.mod   include/framework/

mkdir -p include/core_atmosphere
cp -rf ${path_monan}/core_atmosphere/*.mod   include/core_atmosphere/
cp -rf ${path_monan}/core_atmosphere/*.h     include/core_atmosphere/
cp -rf ${path_monan}/core_atmosphere/*.inc   include/core_atmosphere/

mkdir -p  include/core_atmosphere/physics
cp -rf ${path_monan}/core_atmosphere/physics/*.mod   include/core_atmosphere/physics/
cp -rf ${path_monan}/core_atmosphere/physics/*.h     include/core_atmosphere/physics/
cp -rf ${path_monan}/core_atmosphere/physics/*.inc   include/core_atmosphere/physics/

mkdir -p include/operators
cp -rf ${path_monan}/operators/*.mod   include/operators/
cp -rf ${path_monan}/operators/*.h     include/operators/
cp -rf ${path_monan}/operators/*.inc   include/operators/

mkdir -p include/external/SMIOL
cp -rf ${path_monan}/external/SMIOL/*.mod   include/external/SMIOL/
cp -rf ${path_monan}/external/SMIOL/*.h     include/external/SMIOL/
cp -rf ${path_monan}/external/SMIOL/*.inc   include/external/SMIOL/













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
