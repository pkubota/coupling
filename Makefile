# ==============================================================================
# Makefile - Acoplamento MONAN-A 2.0 / NUOPC-ESMF 8.9.1
# Versao : 9.1 - Correcao de link: simbolos MPAS undefined (mpas_log, mpas_pool, etc.)
#
# Historico:
#   v8.0 - src/ em arvore: caps/atmos/, caps/ocean/, mediator/, driver/, main/
#   v8.1 - DATOCN_cap.F90 -> DOCN_cap.F90
#   v8.2 - B-51/B-52/B-53: escalabilidade validada 4/128/512 PETs
#   v9.0 - Merge: caminhos hardcoded do original + estrutura robusta do draft
#   v9.1 - FIX LINK: 3 bugs corrigidos:
#          Bug 1: MPAS_LIBDIR definido como flags -L (devia ser caminho puro)
#          Bug 2: MPAS_OBJS com wildcard duplo $(MPAS_LIBDIR)/lib/... (incorreto)
#          Bug 3: LDFLAGS_ALL sem -lframework -ldycore -lphys -lops -lsmiol
#                 -> todos os simbolos __mpas_*_MOD_* ficavam undefined
#
# Prerequisitos (carregar ANTES de make):
#   module purge
#   module load PrgEnv-gnu craype-x86-turin
#   module load cray-hdf5/1.14.3.3
#   module load cray-netcdf/4.9.0.15
#   module load cray-parallel-netcdf/1.12.3.15
#   module load cray-pals
#   source setenv-monan-gnu.bash   # define ESMFMKFILE, MPAS_DIR
#
# Uso tipico:
#   source setenv-monan-gnu.bash
#   make clean && make
#   make test NPES=4
#   make test NPES=128
#   make printenv
#   make check_stubs
# ==============================================================================
#==============================================================================#
# Makefile para o acoplamento MONAN + DATM + MED + MOM6+SIS2                    #
#/p/projetos/monan_adm/daniel.massaru/Acopladores/NUOPC-MPAS-Integrado-V4.2
#==============================================================================#
ESMFMKFILE := /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/lib/libg/Linux/esmf.mk
MPAS_DIR := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/monan
MPAS_DIR_LOCAL:= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3
# ------------------------------------------------------------------------------
# Validacao de variaveis obrigatorias de ambiente
# (devem vir de setenv-monan-gnu.bash)
# ------------------------------------------------------------------------------
ifndef ESMFMKFILE
  $(error ESMFMKFILE nao definido. Execute: source setenv-monan-gnu.bash)
endif

include $(ESMFMKFILE)

ifndef MPAS_DIR
  $(error MPAS_DIR nao definido. Execute: source setenv-monan-gnu.bash)
endif

# ------------------------------------------------------------------------------
# Caminhos das bibliotecas externas (Jaci/INPE - do Makefile original)
# Sobrescreva via linha de comando se necessario:
#   make MOM6_LIBDIR=/outro/caminho/lib/mom6
# ------------------------------------------------------------------------------

# -- MOM6 + SIS2 ---------------------------------------------------------------
MOM6_LIBDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/lib/mom6
MOM6_INCDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/include/mom6
MOM6_MODDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/mod/mom6

# -- FMS -----------------------------------------------------------------------
FMS_LIBDIR  ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/lib/fms
FMS_INCDIR  ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/include/fms
FMS_MODDIR  ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/mod/fms

# -- NUOPC cap MOM6 ------------------------------------------------------------
NUOPC_LIBDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/lib/nuopc
NUOPC_INCDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/include/nuopc
NUOPC_MODDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/mom6+sis2/mod/nuopc

# -- ESMF (caminhos explicitos alem do esmf.mk) --------------------------------
ESMF_DIR    ?= /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf
ESMF_INC    := -I$(ESMF_DIR)/mod/modg/Linux \
               -I$(ESMF_DIR)/include \
               -I$(ESMF_DIR)/esmf-8.9.0/mod/modO/Linux.gfortran.64.mpich2.default
ESMF_LIB    := -L$(ESMF_DIR)/lib/libg/Linux -lesmf

# -- MOAB (regridding) ---------------------------------------------------------
MOAB_DIR    ?= /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/libmoab
MOAB_INC    := -I$(MOAB_DIR)/include
MOAB_LIB    := -L$(MOAB_DIR)/lib -lMOAB

# -- PNetCDF -------------------------------------------------------------------
PNETCDF_DIR ?= /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/pnetcdf
PNETCDF_INC := -I$(PNETCDF_DIR)/include
PNETCDF_LIB := -L$(PNETCDF_DIR)/lib -lpnetcdf

# -- HDF5 privado (mesma versao com que libfms.a foi compilado) ----------------
# B-ABI-02: em runtime o Cray disponibiliza cray-hdf5-parallel que tem versao
# diferente. Usar RPATH para fixar a versao correta no executavel.
HDF5_DIR    ?= /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/hdf5

# -- NetCDF (C + Fortran) ------------------------------------------------------
NETCDF_DIR  ?= /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/netcdf
NETCDF_INC  := -I$(NETCDF_DIR)/include
NETCDF_LIB  := -L$(NETCDF_DIR)/lib -lnetcdf -lnetcdff -lnetcdf_c++4

# ------------------------------------------------------------------------------
# Validacao de caminhos criticos
# ------------------------------------------------------------------------------
ifndef MOM6_LIBDIR
  $(error MOM6_LIBDIR nao definido)
endif
ifndef FMS_LIBDIR
  $(error FMS_LIBDIR nao definido)
endif

# ------------------------------------------------------------------------------
# Arvore de diretorios do projeto (estrutura v8.0)
#/p/projetos/monan_adm/daniel.massaru/Acopladores/scripts_CD-CT/sources/MONAN-Model_feature/monan_2.0.0/src
# ------------------------------------------------------------------------------
#SRC_DIR      := $(MPAS_DIR_LOCAL)/src
SRC_DIR      := $(MPAS_DIR_LOCAL)
ATM_DIR      := $(SRC_DIR)/caps/atmos
OCN_DIR      := $(SRC_DIR)/caps/ocean
MEDIATOR_DIR := $(SRC_DIR)/mediator
DRIVER_DIR   := $(SRC_DIR)/driver
MAIN_DIR     := $(SRC_DIR)/main
BINDIR       := bin
DATAOUTDIR   := diag_export
LOGSDIR      := logs
BUILDDIR     := build
OBJDIR       := $(BUILDDIR)/obj
MODDIR       := $(BUILDDIR)/mod

# ------------------------------------------------------------------------------
# Executavel
# ------------------------------------------------------------------------------
TARGET := $(BINDIR)/esmApp

# ------------------------------------------------------------------------------
# Compilador (herdado do esmf.mk)
# ------------------------------------------------------------------------------
FC := $(ESMF_F90COMPILER)  -g -fcheck=all

# ------------------------------------------------------------------------------
# Includes do MONAN-A 2.0
# MPAS_DIR := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.3/models/monan
# ------------------------------------------------------------------------------
MPAS_MOD_DIR ?= $(MPAS_DIR)/include/framework

MPAS_INCLUDE := -I$(MPAS_MOD_DIR)                    \
                -I$(MPAS_DIR)/include/framework           \
                -I$(MPAS_DIR)/include/core_atmosphere     \
                -I$(MPAS_DIR)/include/core_atmosphere/physics     \
                -I$(MPAS_DIR)/include/operators           \
                -I$(MPAS_DIR)/include/external/SMIOL

MPAS_LIBDIR :=  $(MPAS_DIR)/lib

# ------------------------------------------------------------------------------
# F90FLAGS
# Combina: flags do esmf.mk + MPAS + MOM6/FMS/NUOPC mods + moddir proprio
#
# B-ABI-01: HÁ DOIS CONJUNTOS DE FLAGS NECESSÁRIOS POR INCOMPATIBILIDADE ABI:
#
#   1. F90FLAGS_MPAS: COM -fdefault-real-8 (e -fdefault-double-8)
#      Usado para arquivos do MPAS e que dependem de mod do MPAS.
#      MPAS é compilado com -fdefault-real-8 (default GFortran MPAS),
#      então qualquer arquivo que use mod do MPAS precisa da mesma flag.
#
#   2. F90FLAGS_FMS: SEM -fdefault-real-8 (e SEM -fdefault-double-8)
#      Usado para arquivos do cap MOM6/FMS.
#      FMS e MOM6 foram compilados via cray-gnu.mk SEM essas flags
#      (default GFDL). Compilar mom_cap.F90 COM -fdefault-real-8 causa
#      layout de struct (FmsNetcdfDomainFile_t) diferente entre caller
#      e biblioteca, manifestando-se como SIGSEGV em get_variable_num_dimensions
#      durante set_grid_metrics_from_mosaic.
#
#   3. Arquivos intermediários (esm.F90, MED_cap.F90, DATM_cap.F90,
#      ocn_comp_NUOPC.F90, esmApp.F90, mpas_cap_utils.F90, mpas_cap_config.F90):
#      compilar SEM -fdefault-real-8. Esses arquivos tocam ambos os mundos
#      mas declaram tipos com kind explícito (ESMF_KIND_R8, ESMF_KIND_R4) e
#      não fazem uso direto do "real" default. Compilar sem a flag mantém
#      compatibilidade com FMS e força declarações explícitas.
# ------------------------------------------------------------------------------

F90FLAGS_COMMON := $(ESMF_F90COMPILEOPTS)       \
                   $(ESMF_F90COMPILEPATHS)      \
                   $(ESMF_F90COMPILEFREECPP)    \
                   $(ESMF_F90COMPILECPPFLAGS)   \
                   $(MPAS_INCLUDE)              \
                   $(ESMF_INC)                  \
                   $(MOAB_INC)                  \
                   $(PNETCDF_INC)               \
                   $(NETCDF_INC)                \
                   -I$(MOM6_MODDIR)             \
                   -I$(FMS_MODDIR)              \
                   -I$(NUOPC_MODDIR)            \
                   -I$(MODDIR)                  \
                   -J$(MODDIR)                  \
                   -fconvert=big-endian         \
                   -ffree-form                  \
                   -ffree-line-length-none      \
                   -O2 -fPIC -m64

# Para arquivos que dependem de mods do MPAS (MPAS é compilado com -fdefault-real-8)
F90FLAGS_MPAS  := $(F90FLAGS_COMMON)            \
                  -fdefault-real-8              \
                  -fdefault-double-8

# Para arquivos do cap MOM6/FMS (FMS foi compilado sem -fdefault-real-8)
F90FLAGS_FMS   := $(F90FLAGS_COMMON)

# Alias antigo para retrocompatibilidade (era F90FLAGS sem distinção).
# Mantemos apontando para a variante MPAS por enquanto, mas o IDEAL é que
# nenhuma regra use $(F90FLAGS) ? sempre escolher F90FLAGS_MPAS ou F90FLAGS_FMS.
F90FLAGS := $(F90FLAGS_MPAS)
#            -fallow-argument-mismatch       \
#            -ffpe-summary=none              \
#            -O2 -g                          \
#            -fPIC                           \
#            -m64                            \
#            -Wall                           \
#            -Wno-unused-dummy-argument      \
#            -Wno-unused-variable

# ------------------------------------------------------------------------------
# Bibliotecas MOM6: objetos precompilados (exclui mains proprios do MOM6)
# ------------------------------------------------------------------------------
MOM6_OBJS := $(filter-out                  \
    $(MOM6_LIBDIR)/MOM_main.o              \
    $(MOM6_LIBDIR)/coupler_main.o          \
    $(MOM6_LIBDIR)/MOM_driver.o,           \
    $(wildcard $(MOM6_LIBDIR)/*.o))

# Objetos precompilados do MONAN-A (mesmo mecanismo do MOM6_OBJS).
# Se o MONAN for instalado como .a (static libs), use MPAS_LIBS abaixo.
# Se for como .o individuais em lib/framework/, lib/core_atmosphere/ etc.,
# o wildcard captura tudo automaticamente.
MPAS_OBJS := $(wildcard                                  \
    $(MPAS_LIBDIR)/framework/*.o                         \
    $(MPAS_LIBDIR)/core_atmosphere/*.o                   \
    $(MPAS_LIBDIR)/core_atmosphere/physics/*.o           \
    $(MPAS_LIBDIR)/operators/*.o                         \
    $(MPAS_LIBDIR)/external/SMIOL/*.o)

# Fallback: se o MONAN for instalado como bibliotecas estaticas (.a),
# MPAS_OBJS ficara vazio e MPAS_LIBS sera usado no lugar.
# Descomentar se necessario:
# MPAS_LIBS := -L$(MPAS_LIBDIR)/framework          \
#              -L$(MPAS_LIBDIR)/core_atmosphere     \
#              -L$(MPAS_LIBDIR)/operators           \
#              -L$(MPAS_LIBDIR)/external/SMIOL      \
#              -Wl,--start-group                    \
#              -lframework -ldycore -lphys        \
#              -lops -lsmiolf -lsmiol             \
#              -Wl,--end-group
# ------------------------------------------------------------------------------
# LDFLAGS_ALL - ordem de link critica:
#   1. objetos proprios (ALL_OBJS, passados no link)
#   2. objetos MOM6 precompilados
#   3. libs NUOPC/FMS/ESMF
#   4. libs NetCDF/PNetCDF/MOAB
#   5. libs sistema
# ------------------------------------------------------------------------------
SYS_LIBS  := -lz -ldl -lm
OMP_LIBS  := -lgomp

# Bibliotecas estaticas do MONAN-A (usadas quando MPAS_OBJS estiver vazio
# ou para garantir resolucao de simbolos nao capturados pelos .o).
# Ordem critica: framework antes de dycore/phys (dycore depende de framework).
MPAS_LINK := -L$(MPAS_LIBDIR)/framework                 \
              -L$(MPAS_LIBDIR)/core_atmosphere           \
              -L$(MPAS_LIBDIR)/core_atmosphere/physics   \
              -L$(MPAS_LIBDIR)/operators                 \
              -L$(MPAS_LIBDIR)/external/SMIOL            \
              -Wl,--start-group                          \
                -lframework -ldycore -lphys              \
                -lops -lsmiolf -lsmiol                   \
              -Wl,--end-group

LDFLAGS_ALL := $(MOM6_OBJS)                    \
               $(MPAS_OBJS)                     \
               $(MPAS_LINK)                     \
               -L$(NUOPC_LIBDIR) -lmom6_nuopc  \
               -L$(FMS_LIBDIR)   -lfms         \
               $(ESMF_LIB)                     \
               $(ESMF_F90LINKOPTS)             \
               $(ESMF_F90LINKPATHS)            \
               $(ESMF_F90LINKRPATHS)           \
               $(ESMF_F90ESMFLINKLIBS)         \
               $(MOAB_LIB)                     \
               $(PNETCDF_LIB)                  \
               $(NETCDF_LIB)                   \
               -lpioc                          \
               $(SYS_LIBS)                     \
               $(OMP_LIBS)                     \
               -Wl,-rpath,$(FMS_LIBDIR)        \
               -Wl,-rpath,$(NUOPC_LIBDIR)      \
               -Wl,-rpath,$(NETCDF_DIR)/lib    \
               -Wl,-rpath,$(HDF5_DIR)/lib

# ==============================================================================
# Fontes e objetos
#
# Grafo de dependencias - 6 camadas:
#
#  Camada 0 - sem dependencia interna de projeto:
#    mpas_cap_utils    mpas_cap_config    mpas_atm_types
#
#  Camada 1 - dependem de camada 0:
#    mpas_atm_model  -> mpas_atm_types + mpas_cap_config + mpas_cap_utils
#    mpas_cap_netcdf -> mpas_cap_config + mpas_cap_utils
#
#  Camada 2a - dependem de camada 1 (ATM):
#    mpas_atm_wrappers -> mpas_atm_model
#    mpas_cap_methods  -> mpas_atm_types + mpas_atm_model + mpas_cap_utils
#
#  Camada 2b - OCN interna (dependem so de ESMF/FMS/MOM6):
#    time_utils        -> FMS + ESMF
#    mom_cap_methods   -> ESMF
#    mom_cap_time      -> ESMF + mom_cap_methods
#
#  Camada 2c - caps independentes (so ESMF + mpas_cap_config):
#    DATM_cap          -> mpas_cap_config
#    ocn_comp_NUOPC    -> mpas_cap_config   (DOCN)
#
#  Camada 3 - caps principais:
#    mpas_cap -> mpas_cap_{config,methods,netcdf,utils}
#             + mpas_atm_{types,model,wrappers}
#    mom_cap  -> mom_cap_methods + mom_cap_time + time_utils
#             + mpas_cap_utils + MOM6/FMS/NUOPC libs
#    MED_cap  -> mpas_cap_utils
#
#  Camada 4 - driver NUOPC:
#    esm -> mpas_cap + DATM_cap + MED_cap + mom_cap + ocn_comp_NUOPC
#         + mpas_cap_config + mpas_cap_utils
#
#  Camada 5 - aplicativo principal:
#    esmApp -> esm + mpas_cap_config + mpas_cap_utils
# ==============================================================================

# -- ATM (src/caps/atmos/) -----------------------------------------------------
ATM_SRCS := $(ATM_DIR)/mpas_atm_types.F90      \
            $(ATM_DIR)/mpas_atm_model.F90       \
            $(ATM_DIR)/mpas_atm_wrappers.F90    \
            $(ATM_DIR)/mpas_cap_config.F90      \
            $(ATM_DIR)/mpas_cap_utils.F90       \
            $(ATM_DIR)/mpas_cap_methods.F90     \
            $(ATM_DIR)/mpas_cap_netcdf.F90      \
            $(ATM_DIR)/mpas_cap.F90             \
            $(ATM_DIR)/DATM_cap.F90

# -- OCN (src/caps/ocean/) -----------------------------------------------------
OCN_SRCS := $(OCN_DIR)/time_utils.F90          \
            $(OCN_DIR)/mom_cap_methods.F90      \
            $(OCN_DIR)/mom_cap_time.F90         \
            $(OCN_DIR)/mom_cap.F90              \
            $(OCN_DIR)/ocn_comp_NUOPC.F90

# -- Mediador, Driver, Main ----------------------------------------------------
MEDIATOR_SRCS := $(MEDIATOR_DIR)/MED_cap.F90
DRIVER_SRCS   := $(DRIVER_DIR)/esm.F90
MAIN_SRCS     := $(MAIN_DIR)/esmApp.F90

# -- Objetos derivados ---------------------------------------------------------
ATM_OBJS      := $(patsubst $(ATM_DIR)/%.F90,      $(OBJDIR)/%.o, $(ATM_SRCS))
OCN_OBJS      := $(patsubst $(OCN_DIR)/%.F90,      $(OBJDIR)/%.o, $(OCN_SRCS))
MEDIATOR_OBJS := $(patsubst $(MEDIATOR_DIR)/%.F90, $(OBJDIR)/%.o, $(MEDIATOR_SRCS))
DRIVER_OBJS   := $(patsubst $(DRIVER_DIR)/%.F90,   $(OBJDIR)/%.o, $(DRIVER_SRCS))
MAIN_OBJS     := $(patsubst $(MAIN_DIR)/%.F90,     $(OBJDIR)/%.o, $(MAIN_SRCS))
ALL_OBJS      := $(ATM_OBJS) $(OCN_OBJS) $(MEDIATOR_OBJS) $(DRIVER_OBJS) $(MAIN_OBJS)

# ==============================================================================
# Alvos principais
# ==============================================================================
.PHONY: all dirs test clean distclean printenv check_stubs help

all: $(TARGET)
	@echo ""
	@echo "  OK: executavel '$(TARGET)' gerado."
	@echo ""

dirs:
	@mkdir -p $(OBJDIR) $(MODDIR) $(BINDIR) $(DATAOUTDIR) $(LOGSDIR) \
	           diag_import diag_import/postproc $(DATAOUTDIR)/postproc \
	           $(ATM_DIR) $(OCN_DIR) $(MEDIATOR_DIR) $(DRIVER_DIR) $(MAIN_DIR)

$(TARGET): $(ALL_OBJS) | dirs
	$(FC) -o $@ $^ $(LDFLAGS_ALL)

# ==============================================================================
# Camada 0 - sem dependencia interna de projeto
# Arquivos MPAS (precisam de -fdefault-real-8 via F90FLAGS_MPAS)
# ==============================================================================

$(OBJDIR)/mpas_cap_utils.o: $(ATM_DIR)/mpas_cap_utils.F90 | dirs
	$(FC) $(F90FLAGS_MPAS) -c -o $@ $<

$(OBJDIR)/mpas_cap_config.o: $(ATM_DIR)/mpas_cap_config.F90 | dirs
	$(FC) $(F90FLAGS_MPAS) -c -o $@ $<

$(OBJDIR)/mpas_atm_types.o: $(ATM_DIR)/mpas_atm_types.F90 | dirs
	$(FC) $(F90FLAGS_MPAS) -c -o $@ $<

# ==============================================================================
# Camada 1 - dependem de camada 0
# ==============================================================================

$(OBJDIR)/mpas_atm_model.o: $(ATM_DIR)/mpas_atm_model.F90  \
                              $(OBJDIR)/mpas_atm_types.o    \
                              $(OBJDIR)/mpas_cap_config.o   \
                              $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS_MPAS) -c -o $@ $<

$(OBJDIR)/mpas_cap_netcdf.o: $(ATM_DIR)/mpas_cap_netcdf.F90 \
                              $(OBJDIR)/mpas_cap_config.o    \
                              $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS_MPAS) -c -o $@ $<

# ==============================================================================
# Camada 2a - ATM: dependem de camada 1
# ==============================================================================

$(OBJDIR)/mpas_atm_wrappers.o: $(ATM_DIR)/mpas_atm_wrappers.F90 \
                                $(OBJDIR)/mpas_atm_model.o
	$(FC) $(F90FLAGS_MPAS) -c -o $@ $<

$(OBJDIR)/mpas_cap_methods.o: $(ATM_DIR)/mpas_cap_methods.F90 \
                               $(OBJDIR)/mpas_atm_types.o     \
                               $(OBJDIR)/mpas_atm_model.o     \
                               $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS_MPAS) -c -o $@ $<

# ==============================================================================
# Camada 2b - OCN interna: time_utils, mom_cap_methods, mom_cap_time
# B-ABI-01: SEM -fdefault-real-8 (compatibilidade ABI com libfms.a)
# ==============================================================================

$(OBJDIR)/time_utils.o: $(OCN_DIR)/time_utils.F90 | dirs
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

$(OBJDIR)/mom_cap_methods.o: $(OCN_DIR)/mom_cap_methods.F90 | dirs
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

$(OBJDIR)/mom_cap_time.o: $(OCN_DIR)/mom_cap_time.F90 \
                           $(OBJDIR)/mom_cap_methods.o
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

# ==============================================================================
# Camada 2c - Caps independentes: DATM e DOCN (ocn_comp_NUOPC)
# B-ABI-01: SEM -fdefault-real-8.
# ==============================================================================

$(OBJDIR)/DATM_cap.o: $(ATM_DIR)/DATM_cap.F90      \
                      $(OBJDIR)/mpas_cap_config.o
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

$(OBJDIR)/ocn_comp_NUOPC.o: $(OCN_DIR)/ocn_comp_NUOPC.F90 \
                              $(OBJDIR)/mpas_cap_config.o
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

# ==============================================================================
# Camada 3 - Caps principais: mpas_cap, mom_cap, MED_cap
# ==============================================================================

# mpas_cap: depende de tudo do MPAS. F90FLAGS_MPAS obrigatorio.
$(OBJDIR)/mpas_cap.o: $(ATM_DIR)/mpas_cap.F90        \
                      $(OBJDIR)/mpas_cap_config.o     \
                      $(OBJDIR)/mpas_atm_types.o      \
                      $(OBJDIR)/mpas_atm_model.o      \
                      $(OBJDIR)/mpas_atm_wrappers.o   \
                      $(OBJDIR)/mpas_cap_methods.o    \
                      $(OBJDIR)/mpas_cap_netcdf.o     \
                      $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS_MPAS) -c -o $@ $<

# mom_cap: depende de MOM6/FMS. B-ABI-01: F90FLAGS_FMS obrigatorio.
# Era a CAUSA do crash em set_grid_metrics_from_mosaic.
$(OBJDIR)/mom_cap.o: $(OCN_DIR)/mom_cap.F90           \
                     $(OBJDIR)/mom_cap_methods.o       \
                     $(OBJDIR)/mom_cap_time.o          \
                     $(OBJDIR)/time_utils.o            \
                     $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

# MED_cap: troca campos via state, sem acesso direto a "real" do MPAS.
$(OBJDIR)/MED_cap.o: $(MEDIATOR_DIR)/MED_cap.F90      \
                     $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

# ==============================================================================
# Camada 4 - Driver NUOPC (src/driver/)
# ==============================================================================

$(OBJDIR)/esm.o: $(DRIVER_DIR)/esm.F90              \
                 $(OBJDIR)/mpas_cap.o                \
                 $(OBJDIR)/DATM_cap.o                \
                 $(OBJDIR)/MED_cap.o                 \
                 $(OBJDIR)/mom_cap.o                 \
                 $(OBJDIR)/ocn_comp_NUOPC.o          \
                 $(OBJDIR)/mpas_cap_config.o         \
                 $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

# ==============================================================================
# Camada 5 - Aplicativo principal (src/main/)
# ==============================================================================

$(OBJDIR)/esmApp.o: $(MAIN_DIR)/esmApp.F90           \
                    $(OBJDIR)/esm.o                  \
                    $(OBJDIR)/mpas_cap_config.o      \
                    $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS_FMS) -c -o $@ $<

# ==============================================================================
# Execucao local de teste (sem PBS)
# ==============================================================================
NPES ?= 4

test: $(TARGET)
	mpirun -np $(NPES) $(TARGET)

# ==============================================================================
# Limpeza
# ==============================================================================
clean:
	rm -rf $(BUILDDIR) $(BINDIR) $(DATAOUTDIR) diag_import $(LOGSDIR)
	rm -f *.stdout log.atmosphere.*.out

distclean: clean
	rm -f *.pbs esmApp-integrado.pbs

# ==============================================================================
# Diagnostico - variaveis de build
# ==============================================================================
printenv:
	@echo "======================================================================"
	@echo " NUOPC-MPAS-Integrado - variaveis de build (Makefile v9.1)"
	@echo "======================================================================"
	@echo "FC           = $(FC)"
	@echo "ESMFMKFILE   = $(ESMFMKFILE)"
	@echo "MPAS_DIR     = $(MPAS_DIR)"
	@echo "MPAS_MOD_DIR = $(MPAS_MOD_DIR)"
	@echo "MPAS_LIBDIR  = $(MPAS_LIBDIR)"
	@echo "TARGET       = $(TARGET)"
	@echo ""
	@echo "MOM6_LIBDIR  = $(MOM6_LIBDIR)"
	@echo "MOM6_MODDIR  = $(MOM6_MODDIR)"
	@echo "FMS_LIBDIR   = $(FMS_LIBDIR)"
	@echo "FMS_MODDIR   = $(FMS_MODDIR)"
	@echo "NUOPC_LIBDIR = $(NUOPC_LIBDIR)"
	@echo "NUOPC_MODDIR = $(NUOPC_MODDIR)"
	@echo "ESMF_DIR     = $(ESMF_DIR)"
	@echo "MOAB_DIR     = $(MOAB_DIR)"
	@echo "PNETCDF_DIR  = $(PNETCDF_DIR)"
	@echo "NETCDF_DIR   = $(NETCDF_DIR)"
	@echo ""
	@echo "ATM_DIR      = $(ATM_DIR)"
	@echo "OCN_DIR      = $(OCN_DIR)"
	@echo "MEDIATOR_DIR = $(MEDIATOR_DIR)"
	@echo "DRIVER_DIR   = $(DRIVER_DIR)"
	@echo "MAIN_DIR     = $(MAIN_DIR)"
	@echo "OBJDIR       = $(OBJDIR)"
	@echo "MODDIR       = $(MODDIR)"
	@echo ""
	@echo "ATM_SRCS:"
	@echo "  $(ATM_SRCS)" | tr ' ' '\n' | grep '\S' | sed 's/^/    /'
	@echo ""
	@echo "OCN_SRCS:"
	@echo "  $(OCN_SRCS)" | tr ' ' '\n' | grep '\S' | sed 's/^/    /'
	@echo ""
	@echo "ALL_OBJS:"
	@echo "  $(ALL_OBJS)" | tr ' ' '\n' | grep '\S' | sed 's/^/    /'
	@echo ""
	@echo "LDFLAGS_ALL:"
	@echo "  $(LDFLAGS_ALL)" | tr ' ' '\n' | grep '\S' | sed 's/^/    /'
	@echo "======================================================================"

# ==============================================================================
# Verifica ausencia de stubs esmf_time_f90
# (MONAN-A deve ser compilado sem src/external/esmf_time_f90/)
# ==============================================================================
check_stubs: $(TARGET)
	@echo "=== Verificando stubs esmf_time_f90 em $(TARGET) ==="
	@if nm $(TARGET) 2>/dev/null | grep -i 'esmf_timeset' | grep -qiv 'timesetdefault'; then \
	    echo "  FALHA   : Stubs esmf_time_f90 detectados no binario."; \
	    echo "  Execute : cd $$MPAS_DIR && make gfortran-xd2000 EXTRA_FFLAGS=\"-I$${ESMF_MOD}\""; \
	    exit 1; \
	elif nm $(TARGET) 2>/dev/null | grep -qi 'esmf_timesetdefault'; then \
	    echo "  OK      : ESMF 8.9.1 real confirmado."; \
	else \
	    echo "  AVISO   : Simbolo esmf_timeset* ausente no binario."; \
	fi

# ==============================================================================
# Ajuda
# ==============================================================================
help:
	@echo ""
	@echo "Uso: make [alvo] [VAR=valor]"
	@echo ""
	@echo "Estrutura do projeto (v8.0+):"
	@echo "  src/caps/atmos/  - cap ATM + interface MONAN-A (9 fontes)"
	@echo "  src/caps/ocean/  - cap OCN: mom_cap, ocn_comp_NUOPC + auxiliares"
	@echo "  src/mediator/    - mediador MED_cap.F90"
	@echo "  src/driver/      - driver ESM (esm.F90)"
	@echo "  src/main/        - aplicativo principal (esmApp.F90)"
	@echo "  bin/             - executavel gerado (esmApp)"
	@echo "  diag_export/     - campos ATM exportados (monan_export_*.nc)"
	@echo "  diag_import/     - campos OCN importados (docn_import_*.nc)"
	@echo "  logs/            - arquivos de log ESMF"
	@echo "  build/obj/       - objetos intermediarios"
	@echo "  build/mod/       - modulos Fortran gerados"
	@echo ""
	@echo "Alvos:"
	@echo "  all          compila bin/esmApp  (padrao)"
	@echo "  dirs         cria a arvore de diretorios"
	@echo "  test         executa: mpirun -np NPES bin/esmApp"
	@echo "  clean        remove build/, bin/, diag_export/, diag_import/, logs/"
	@echo "  distclean    idem clean + scripts PBS"
	@echo "  printenv     mostra variaveis de compilacao e link"
	@echo "  check_stubs  verifica ausencia de stubs ESMF"
	@echo "  help         esta mensagem"
	@echo ""
	@echo "Variaveis sobrescrevieis (linha de comando ou ambiente):"
	@echo "  NPES=N            processos MPI para 'make test'   (padrao: 4)"
	@echo "  MPAS_MOD_DIR=...  diretorio dos .mod do MONAN-A"
	@echo "  MOM6_LIBDIR=...   diretorio lib/ do MOM6+SIS2"
	@echo "  MOM6_MODDIR=...   diretorio mod/ do MOM6+SIS2"
	@echo "  FMS_LIBDIR=...    diretorio lib/ do FMS"
	@echo "  NETCDF_DIR=...    raiz do NetCDF"
	@echo "  PNETCDF_DIR=...   raiz do PNetCDF"
	@echo "  MOAB_DIR=...      raiz do MOAB"
	@echo ""
