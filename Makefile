#==============================================================================#
# Makefile para o acoplamento MPAS + DATM + MED + MOM6+SIS2                    #
#==============================================================================#
ESMFMKFILE := /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf/lib/libg/Linux/esmf.mk
#/p/projetos/monan_atm/paulo.kubota/coupler/system_coupler_arquitetura_00/MOM6-examples/build/gnu/ice_ocean_SIS2
MOM6_LIBDIR := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/lib/mom6
MOM6_INCDIR := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/include/mom6
MOM6_MODDIR := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/mod/mom6
#/p/projetos/monan_atm/paulo.kubota/coupler/system_coupler_arquitetura_00/MOM6-examples/build/gnu/shared/repro
FMS_LIBDIR  := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/lib/fms
FMS_INCDIR  := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/include/fms
FMS_MODDIR  := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/mod/fms
#/p/projetos/monan_atm/paulo.kubota/coupler/system_coupler_arquitetura_00/MOM6-examples/build/gnu/nuopc_cap/repro
NUOPC_LIBDIR   := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/lib/nuopc
NUOPC_INCDIR   := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/include/nuopc
NUOPC_MODDIR   := /p/projetos/monan_atm/paulo.kubota/coupler/coupling_0.0.0/models/mom6+sis2/mod/nuopc
#NETCDF_DIR := /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/lib
include $(ESMFMKFILE)

# --- Validacao ---
ifndef ESMFMKFILE
  $(error "ESMFMKFILE nao definido")
endif
ifndef MOM6_LIBDIR
  $(error "MOM6_LIBDIR nao definido")
endif
ifndef FMS_LIBDIR
  $(error "FMS_LIBDIR nao definido")
endif

# Compilador e flags
FC = $(ESMF_F90COMPILER) -g -fcheck=all
FCFLAGS := $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) -O2 -fPIC -m64 -ffree-line-length-none \
           $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS)


LINKER  := $(ESMF_F90LINKER)
LFLAGS  := $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) \
           $(ESMF_F90LINKRPATHS) $(ESMF_EXTERNAL_LIBS)

# --- Paths de include ---
FCFLAGS += -I$(MOM6_MODDIR)  # -I$(MOM6_INCDIR)
FCFLAGS += -I$(NUOPC_MODDIR) # -I$(NUOPC_INCDIR)
FCFLAGS += -I$(FMS_MODDIR)   # -I$(FMS_MODDIR)
# --- Objetos MOM6 (exclui main e driver proprios do MOM6) ---
MOM6_OBJS := $(filter-out \
    $(MOM6_LIBDIR)/MOM_main.o \
    $(MOM6_LIBDIR)/coupler_main.o \
    $(MOM6_LIBDIR)/MOM_driver.o, \
    $(wildcard $(MOM6_LIBDIR)/*.o))

# --- Bibliotecas ---
LIBS := $(MOM6_OBJS)                       \
        -L$(NUOPC_LIBDIR) -lmom6_nuopc        \
        -L$(FMS_LIBDIR) -lfms              \
        $(ESMF_F90ESMFLINKLIBS)

ifdef NETCDF_DIR
  FCFLAGS += -I$(NETCDF_DIR)/include
  LIBS    += -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf
endif

# Diretorios ESMF
ESMF_DIR = /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/esmf
ESMF_INC = -I$(ESMF_DIR)/mod/modg/Linux -I$(ESMF_DIR)/include -I$(ESMF_DIR)/esmf-8.9.0/mod/modO/Linux.gfortran.64.mpich2.default
ESMF_LIB = -L$(ESMF_DIR)/lib/libg/Linux -lesmf

# Diretorios MOAB
MOAB_DIR = /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/libmoab
MOAB_INC = -I$(MOAB_DIR)/include
MOAB_LIB = -L$(MOAB_DIR)/lib -lMOAB

# Diretorios MOAB
PNETCDF_DIR = /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/pnetcdf
PNETCDF_INC = -I$(PNETCDF_DIR)/include
PNETCDF_LIB = -L$(PNETCDF_DIR)/lib -lpnetcdf
# Diretorios MOAB
NETCDF_DIR = /p/projetos/monan_adm/paulo.kubota/home/lib/lib_gnucray/netcdf
NETCDF_INC = -I$(NETCDF_DIR)/include
NETCDF_LIB = -L$(NETCDF_DIR)/lib  -lnetcdf -lnetcdff -lnetcdf_c++4

# Todos os includes
INCLUDES = $(ESMF_INC) $(MOAB_INC) $(PNETCDF_INC) $(NETCDF_INC)

# Todas as bibliotecas
LIBS := $(LIBS) $(ESMF_LIB) $(MOAB_LIB) $(PNETCDF_LIB) $(NETCDF_LIB) -lpioc
# Arquivos fonte
SOURCES = esmApp.F90 \
          esm.F90 \
          caps/MPAS_cap.F90 \
          caps/DATM_cap.F90 \
          mediator/MED_cap.F90

# Objetos
OBJECTS = esmApp.o esm.o MPAS_cap.o DATM_cap.o MED_cap.o

# Executavel
EXEC = exec/esmApp.x

# Regra principal
all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(FC) -o $@ $^ $(LIBS)

%.o: %.F90
	$(FC) -c $(FCFLAGS) $(INCLUDES) $< -o $@

esmApp.o: main/esmApp.F90 esm.o
	$(FC) -c $(FCFLAGS) $(INCLUDES) $< -o $@

esm.o: driver/esm.F90 MPAS_cap.o DATM_cap.o MED_cap.o
	$(FC) -c $(FCFLAGS) $(INCLUDES) $< -o $@

MPAS_cap.o: caps/atmos/MPAS_cap.F90
	$(FC) -c $(FCFLAGS) $(INCLUDES) $< -o $@

DATM_cap.o: caps/atmos/DATM_cap.F90
	$(FC) -c $(FCFLAGS) $(INCLUDES) $< -o $@

MED_cap.o: mediator/MED_cap.F90
	$(FC) -c $(FCFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -f *.o *.mod $(EXEC)

.PHONY: all clean
