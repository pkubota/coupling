# ==============================================================================
# Makefile - Acoplamento MONAN-A 2.0 / NUOPC-ESMF 8.9.1
# Versao : 9.0 - Merge Makefile_original + Makefile_massaru v8.2
#
# Historico:
#   v8.0 - src/ em arvore: caps/atmos/, caps/ocean/, mediator/, driver/, main/
#   v8.1 - DATOCN_cap.F90 -> DOCN_cap.F90
#   v8.2 - B-51/B-52/B-53: escalabilidade validada 4/128/512 PETs
#   v9.0 - Merge: caminhos hardcoded do original + estrutura robusta do draft
#          Adicionado: ocn_comp_NUOPC.F90, DATM_cap.F90
#          Adicionado: mom_cap_methods.F90, mom_cap_time.F90, time_utils.F90
#          Corrigido:  grafo de dependencias completo (5 camadas + OCN interna)
#          MOM6/FMS/MOAB/PNetCDF/NetCDF: caminhos do original (Jaci/INPE)
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
MPAS_DIR := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.2/models/monan
MPAS_DIR_LOCAL:= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.2
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
MOM6_LIBDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/lib/mom6
MOM6_INCDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/include/mom6
MOM6_MODDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/mod/mom6

# -- FMS -----------------------------------------------------------------------
FMS_LIBDIR  ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/lib/fms
FMS_INCDIR  ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/include/fms
FMS_MODDIR  ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/mod/fms

# -- NUOPC cap MOM6 ------------------------------------------------------------
NUOPC_LIBDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/lib/nuopc
NUOPC_INCDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/include/nuopc
NUOPC_MODDIR ?= /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/mod/nuopc

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
FC := $(ESMF_F90COMPILER)

# ------------------------------------------------------------------------------
# Includes do MONAN-A 2.0
# MPAS_DIR := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.2/models/monan
# ------------------------------------------------------------------------------
MPAS_MOD_DIR ?= $(MPAS_DIR)/include/framework

MPAS_INCLUDE := -I$(MPAS_MOD_DIR)                    \
                -I$(MPAS_DIR)/include/framework           \
                -I$(MPAS_DIR)/include/core_atmosphere     \
                -I$(MPAS_DIR)/include/operators           \
                -I$(MPAS_DIR)/include/external/SMIOL

MPAS_LIBDIR :=  -L$(MPAS_DIR)/lib/framework           \
                -L$(MPAS_DIR)/lib/core_atmosphere     \
                -L$(MPAS_DIR)/lib/operators           \
                -L$(MPAS_DIR)/lib/external/SMIOL

# ------------------------------------------------------------------------------
# F90FLAGS
# Combina: flags do esmf.mk + MPAS + MOM6/FMS/NUOPC mods + moddir proprio
# FCFLAGS := $(ESMF_F90COMPILEOPTS) \
#            $(ESMF_F90COMPILEPATHS) \
# 	   -O2 \
# 	   -fPIC \
# 	   -m64 \
# 	   -ffree-line-length-none \
#          $(ESMF_F90COMPILEFREECPP) \
# 	   $(ESMF_F90COMPILECPPFLAGS)
# ------------------------------------------------------------------------------
F90FLAGS := $(ESMF_F90COMPILEOPTS)          \
            $(ESMF_F90COMPILEPATHS)         \
            $(ESMF_F90COMPILEFREECPP)       \
            $(ESMF_F90COMPILECPPFLAGS)      \
            $(MPAS_INCLUDE)                 \
            $(ESMF_INC)                     \
            $(MOAB_INC)                     \
            $(PNETCDF_INC)                  \
            $(NETCDF_INC)                   \
            -I$(MOM6_MODDIR)                \
            -I$(FMS_MODDIR)                 \
            -I$(NUOPC_MODDIR)               \
            -I$(MODDIR)                     \
            -J$(MODDIR)                     \
            -ffree-form                     \
            -ffree-line-length-none         \
            -fopenmp                        \
            -fallow-argument-mismatch       \
            -ffpe-summary=none              \
            -O2 -g                          \
            -fPIC                           \
            -m64                            \
            -Wall                           \
            -Wno-unused-dummy-argument      \
            -Wno-unused-variable

# ------------------------------------------------------------------------------
# Bibliotecas MOM6: objetos precompilados (exclui mains proprios do MOM6)
# ------------------------------------------------------------------------------
MOM6_OBJS := $(filter-out                  \
    $(MOM6_LIBDIR)/MOM_main.o              \
    $(MOM6_LIBDIR)/coupler_main.o          \
    $(MOM6_LIBDIR)/MOM_driver.o,           \
    $(wildcard $(MOM6_LIBDIR)/*.o))

MPAS_OBJS := $(wildcard                    \
    $(MPAS_LIBDIR)/lib/framework/*.o       \
    $(MPAS_LIBDIR)/lib/core_atmosphere/*.o \
    $(MPAS_LIBDIR)/lib/operators/*.o       \
    $(MPAS_LIBDIR)/lib/external/SMIOL/*.o)
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

LDFLAGS_ALL := $(MOM6_OBJS)                    \
               $(MPAS_OBJS)                     \
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
               $(OMP_LIBS)

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
# ==============================================================================

$(OBJDIR)/mpas_cap_utils.o: $(ATM_DIR)/mpas_cap_utils.F90 | dirs
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/mpas_cap_config.o: $(ATM_DIR)/mpas_cap_config.F90 | dirs
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/mpas_atm_types.o: $(ATM_DIR)/mpas_atm_types.F90 | dirs
	$(FC) $(F90FLAGS) -c -o $@ $<

# ==============================================================================
# Camada 1 - dependem de camada 0
# ==============================================================================

$(OBJDIR)/mpas_atm_model.o: $(ATM_DIR)/mpas_atm_model.F90  \
                              $(OBJDIR)/mpas_atm_types.o    \
                              $(OBJDIR)/mpas_cap_config.o   \
                              $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/mpas_cap_netcdf.o: $(ATM_DIR)/mpas_cap_netcdf.F90 \
                              $(OBJDIR)/mpas_cap_config.o    \
                              $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS) -c -o $@ $<

# ==============================================================================
# Camada 2a - ATM: dependem de camada 1
# ==============================================================================

$(OBJDIR)/mpas_atm_wrappers.o: $(ATM_DIR)/mpas_atm_wrappers.F90 \
                                $(OBJDIR)/mpas_atm_model.o
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/mpas_cap_methods.o: $(ATM_DIR)/mpas_cap_methods.F90 \
                               $(OBJDIR)/mpas_atm_types.o     \
                               $(OBJDIR)/mpas_atm_model.o     \
                               $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS) -c -o $@ $<

# ==============================================================================
# Camada 2b - OCN interna: time_utils, mom_cap_methods, mom_cap_time
# Dependem apenas de ESMF e FMS (bibliotecas externas ? sem .o interno)
# ==============================================================================

$(OBJDIR)/time_utils.o: $(OCN_DIR)/time_utils.F90 | dirs
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/mom_cap_methods.o: $(OCN_DIR)/mom_cap_methods.F90 | dirs
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/mom_cap_time.o: $(OCN_DIR)/mom_cap_time.F90 \
                           $(OBJDIR)/mom_cap_methods.o
	$(FC) $(F90FLAGS) -c -o $@ $<

# ==============================================================================
# Camada 2c - Caps independentes: DATM e DOCN (ocn_comp_NUOPC)
# Dependem de mpas_cap_config (camada 0)
# ==============================================================================

$(OBJDIR)/DATM_cap.o: $(ATM_DIR)/DATM_cap.F90      \
                      $(OBJDIR)/mpas_cap_config.o
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/ocn_comp_NUOPC.o: $(OCN_DIR)/ocn_comp_NUOPC.F90 \
                              $(OBJDIR)/mpas_cap_config.o
	$(FC) $(F90FLAGS) -c -o $@ $<

# ==============================================================================
# Camada 3 - Caps principais: mpas_cap, mom_cap, MED_cap
# ==============================================================================

$(OBJDIR)/mpas_cap.o: $(ATM_DIR)/mpas_cap.F90        \
                      $(OBJDIR)/mpas_cap_config.o     \
                      $(OBJDIR)/mpas_atm_types.o      \
                      $(OBJDIR)/mpas_atm_model.o      \
                      $(OBJDIR)/mpas_atm_wrappers.o   \
                      $(OBJDIR)/mpas_cap_methods.o    \
                      $(OBJDIR)/mpas_cap_netcdf.o     \
                      $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/mom_cap.o: $(OCN_DIR)/mom_cap.F90           \
                     $(OBJDIR)/mom_cap_methods.o       \
                     $(OBJDIR)/mom_cap_time.o          \
                     $(OBJDIR)/time_utils.o            \
                     $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS) -c -o $@ $<

$(OBJDIR)/MED_cap.o: $(MEDIATOR_DIR)/MED_cap.F90      \
                     $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS) -c -o $@ $<

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
	$(FC) $(F90FLAGS) -c -o $@ $<

# ==============================================================================
# Camada 5 - Aplicativo principal (src/main/)
# ==============================================================================

$(OBJDIR)/esmApp.o: $(MAIN_DIR)/esmApp.F90           \
                    $(OBJDIR)/esm.o                  \
                    $(OBJDIR)/mpas_cap_config.o      \
                    $(OBJDIR)/mpas_cap_utils.o
	$(FC) $(F90FLAGS) -c -o $@ $<

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
	@echo " NUOPC-MPAS-Integrado - variaveis de build (Makefile v9.0)"
	@echo "======================================================================"
	@echo "FC           = $(FC)"
	@echo "ESMFMKFILE   = $(ESMFMKFILE)"
	@echo "MPAS_DIR     = $(MPAS_DIR)"
	@echo "MPAS_MOD_DIR = $(MPAS_MOD_DIR)"
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
