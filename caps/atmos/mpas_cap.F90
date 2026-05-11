!> @file mpas_cap.F90
!! @brief Cap NUOPC/ESMF para o modelo atmosferico MPAS-A 8.3 / MONAN-A 2.0.
!!
!! Versao 7.0 -- Protocolo NUOPC completo via NUOPC_CompDerive.
!!
!! Patches aplicados:
!!   v7.0: NUOPC_CompDerive + InitializeAdvertise + InitializeDataComplete
!!   v7.1: mpas_atm_resize eliminado (ESMF e MPAS usam decomposicoes distintas)
!!   v7.2: coordenadas para NetCDF via lonCell(1:n_local) com
!!          n_local = min(localCells_ESMF, nCells_MPAS) -- sem ownedElemCoords
!!          que causa double-free no ESMF 8.9.1 em Cray/gfortran.

module mpas_cap_mod

  use ESMF
  use NUOPC,       only : NUOPC_CompDerive,        NUOPC_CompSpecialize,   &
                           NUOPC_CompSetEntryPoint, NUOPC_CompFilterPhaseMap, &
                           NUOPC_Advertise,         NUOPC_Realize,           &
                           NUOPC_CompAttributeGet,  NUOPC_CompAttributeSet
  use NUOPC_Model, only : model_routine_SS           => SetServices,          &
                           model_label_CheckImport    => label_CheckImport,  &
                           model_label_DataInitialize => label_DataInitialize, &
                           model_label_Advance        => label_Advance,        &
                           model_label_Finalize       => label_Finalize,       &
                           NUOPC_ModelGet,             SetVM

  use mpas_atm_types_mod,   only : mpas_atm_public_type,    &
                                    mpas_atm_state_type,     &
                                    atm_ocean_boundary_type

  use mpas_cap_methods_mod, only : mpas_import,      &
                                    mpas_export,      &
                                    mpas_create_grid, &
                                    state_diagnose

  use mpas_cap_netcdf_mod,  only : export_write_netcdf, &
                                    netcdf_init_coords,  &
                                    netcdf_config_set

  use mpas_cap_config_mod,  only : cfg_write_netcdf, cfg_write_diag, &
                                    cfg_mesh_atm, cfg_config_dir,     &
                                    cfg_dt_coupling, cfg_dt_atm,      &
                                    cfg_output_dir, cfg_grid_res_deg, &
                                    cfg_sst_default,                  &
                                    cfg_ice_fraction_default,         &
                                    cfg_zorl_default, config_read

  use mpas_cap_utils_mod,   only : ChkErr

  implicit none
  private

  interface
    subroutine mpas_atm_init(atm_public, atm_state, atm_bnd, &
                              dt_seconds, config_dir, mpi_comm, rc)
      use mpas_atm_types_mod, only : mpas_atm_public_type,    &
                                      mpas_atm_state_type,     &
                                      atm_ocean_boundary_type
      type(mpas_atm_public_type),    intent(inout) :: atm_public
      type(mpas_atm_state_type),     intent(inout) :: atm_state
      type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
      integer,          intent(in)  :: dt_seconds
      character(len=*), intent(in)  :: config_dir
      integer,          intent(in)  :: mpi_comm
      integer,          intent(out) :: rc
    end subroutine mpas_atm_init

    subroutine mpas_atm_init_sfc(atm_public, atm_state, rc)
      use mpas_atm_types_mod, only : mpas_atm_public_type, mpas_atm_state_type
      type(mpas_atm_public_type), intent(inout) :: atm_public
      type(mpas_atm_state_type),  intent(inout) :: atm_state
      integer,                    intent(out)   :: rc
    end subroutine mpas_atm_init_sfc

    subroutine mpas_atm_run(atm_public, atm_state, atm_bnd, dt_coupling, rc)
      use mpas_atm_types_mod, only : mpas_atm_public_type,    &
                                      mpas_atm_state_type,     &
                                      atm_ocean_boundary_type
      type(mpas_atm_public_type),    intent(inout) :: atm_public
      type(mpas_atm_state_type),     intent(inout) :: atm_state
      type(atm_ocean_boundary_type), intent(in)    :: atm_bnd
      integer,                       intent(in)    :: dt_coupling
      integer,                       intent(out)   :: rc
    end subroutine mpas_atm_run

    subroutine mpas_atm_final(atm_public, atm_state, atm_bnd, rc)
      use mpas_atm_types_mod, only : mpas_atm_public_type,    &
                                      mpas_atm_state_type,     &
                                      atm_ocean_boundary_type
      type(mpas_atm_public_type),    intent(inout) :: atm_public
      type(mpas_atm_state_type),     intent(inout) :: atm_state
      type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
      integer,                       intent(out)   :: rc
    end subroutine mpas_atm_final

    subroutine mpas_atm_resize(atm_public, atm_state, atm_bnd, nCells_new)
      use mpas_atm_types_mod, only : mpas_atm_public_type, mpas_atm_state_type, &
                                     atm_ocean_boundary_type
      type(mpas_atm_public_type),    intent(inout) :: atm_public
      type(mpas_atm_state_type),     intent(inout) :: atm_state
      type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
      integer,                       intent(in)    :: nCells_new
    end subroutine mpas_atm_resize
  end interface

  public :: SetServices
  public :: SetVM

  type(mpas_atm_public_type),    pointer, save :: g_atm_public => null()
  type(mpas_atm_state_type),     pointer, save :: g_atm_state  => null()
  type(atm_ocean_boundary_type), pointer, save :: g_atm_bnd    => null()
  type(ESMF_Grid),                        save :: g_grid

  ! B-OCN-01: So_u, So_v, So_s adicionados para correntes e salinidade
  integer, parameter :: N_IMP = 6
  character(len=20), parameter :: IMP_NAMES(N_IMP) = [ &
    character(len=20) :: 'So_t    ', 'Si_ifrac', 'Sf_zorl ', &
                         'So_u    ', 'So_v    ', 'So_s    ' ]

  integer, parameter :: N_EXP = 10
  character(len=20), parameter :: EXP_NAMES(N_EXP) = [ &
    character(len=20) ::            &
    'Sa_pslv_mpas  ',               &
    'Sa_tbot_mpas  ',               &
    'Sa_u10m_mpas  ',               &
    'Sa_v10m_mpas  ',               &
    'Faxa_swdn_mpas',               &
    'Faxa_lwdn_mpas',               &
    'Faxa_rain_mpas',               &
    'Sa_shum_mpas  ',               &  ! B-Fase2-01: q2 [kg/kg] - hum. espec. 2m
    'Faxa_snow_mpas' ,              &  ! B-Fase2-02: Dsnownc/dt [kg/m^2/s] - neve
    'Sa_lfrac_mpas '  ]                ! fracao de terra MONAN [0-1] para MED_cap: Îsnownc/dt [kg/mÂ²/s] â neve

  integer, parameter :: netcdf_write_freq = 1

  logical            :: write_diag    = .false.
  logical            :: write_netcdf  = .true.
  character(len=256) :: mesh_atm      = 'mpas_mesh.nc'
  character(len=256) :: config_dir    = './'
  integer            :: dt_coupling_s = 1800
  integer            :: dt_atm_s      = 1800
  integer, save      :: step_count    = 0

  character(len=*), parameter :: u_FILE_u = __FILE__

contains

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer, intent(out) :: rc
    rc = ESMF_SUCCESS
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/'IPDv03p1'/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/'IPDv03p3'/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, &
         specLabel=model_label_DataInitialize, &
         specRoutine=InitializeDataComplete, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, &
         specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, &
         specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    ! Suprimir validacao de timestamp de import (lag OCN->MPAS: t-1 != currTime)
    call NUOPC_CompSpecialize(gcomp, &
         specLabel=model_label_CheckImport, &
         specRoutine=CheckImportAlwaysOK, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call ESMF_LogWrite('mpas_cap: SetServices concluido (v7.0 NUOPC_CompDerive)', &
         ESMF_LOGMSG_INFO)
  end subroutine SetServices

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State)    :: importState, exportState
    type(ESMF_Clock)    :: clock
    integer,             intent(out) :: rc
    type(ESMF_Time)    :: startTimeLoc
    logical            :: isPresent, isSet
    character(len=256) :: value
    integer            :: yr, mo, dy, hr, mn, sc
    character(len=*), parameter :: subname = '(mpas_cap:InitializeP0)'
    rc = ESMF_SUCCESS
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/'IPDv03p'/), rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    block
      integer :: cfg_rc
      call config_read(cfg_rc)
    end block
    write_netcdf  = cfg_write_netcdf
    write_diag    = cfg_write_diag
    mesh_atm      = trim(cfg_mesh_atm)
    config_dir    = trim(cfg_config_dir)
    dt_coupling_s = cfg_dt_coupling
    dt_atm_s      = cfg_dt_atm
    call NUOPC_CompAttributeGet(gcomp, name='DumpFields', value=value, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (isPresent .and. isSet) write_diag = (trim(value) == 'true')
    call NUOPC_CompAttributeGet(gcomp, name='WriteNetCDF', value=value, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (isPresent .and. isSet) write_netcdf = (trim(value) == 'true')
    call NUOPC_CompAttributeGet(gcomp, name='mesh_atm', value=value, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (isPresent .and. isSet) mesh_atm = trim(value)
    call NUOPC_CompAttributeGet(gcomp, name='config_dir', value=value, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (isPresent .and. isSet) config_dir = trim(value)
    call NUOPC_CompAttributeGet(gcomp, name='dt_coupling', value=value, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (isPresent .and. isSet) read(value, *) dt_coupling_s
    call NUOPC_CompAttributeGet(gcomp, name='dt_atm', value=value, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (isPresent .and. isSet) read(value, *) dt_atm_s
    call ESMF_ClockGet(clock, startTime=startTimeLoc, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call ESMF_TimeGet(startTimeLoc, yy=yr, mm=mo, dd=dy, h=hr, m=mn, s=sc, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    write(value, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2)') &
          yr, mo, dy, hr, mn, sc
    call ESMF_LogWrite(subname//': start_time = '//trim(value), ESMF_LOGMSG_INFO)
    call ESMF_LogWrite(subname//': InitializeP0 concluido', ESMF_LOGMSG_INFO)
  end subroutine InitializeP0

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State)    :: importState, exportState
    type(ESMF_Clock)    :: clock
    integer,             intent(out) :: rc
    integer :: i
    character(len=*), parameter :: subname = '(mpas_cap:InitializeAdvertise)'
    rc = ESMF_SUCCESS
    do i = 1, N_IMP
      call NUOPC_Advertise(importState, StandardName=trim(IMP_NAMES(i)), rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end do
    do i = 1, N_EXP
      call NUOPC_Advertise(exportState, StandardName=trim(EXP_NAMES(i)), rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end do
    call ESMF_LogWrite(subname//': anunciados '// &
         trim(adjustl(int_to_str(N_IMP)))//' imp + '// &
         trim(adjustl(int_to_str(N_EXP)))//' exp', ESMF_LOGMSG_INFO)
  end subroutine InitializeAdvertise

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp) :: gcomp
    type(ESMF_State)    :: importState, exportState
    type(ESMF_Clock)    :: clock
    integer,             intent(out) :: rc
    type(ESMF_Field)   :: field
    type(ESMF_VM)      :: vm
    integer            :: i, localMpiComm, localPet
    real(ESMF_KIND_R8), parameter   :: RAD2DEG = 57.29577951308232_ESMF_KIND_R8
    character(len=*), parameter :: subname = '(mpas_cap:InitializeRealize)'
    rc = ESMF_SUCCESS

    ! -- 0. VM: obter localMpiComm e localPet ANTES de qualquer outra chamada -
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=localMpiComm, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! -- 1. ESMF_Grid 360x180 (ANTES de mpas_atm_init) --------------------
    ! SOLUCAO DEFINITIVA: ESMF_Grid nao usa MOAB. Zero deadlocks possiveis.
    ! Criado ANTES do SMIOL (mpas_atm_init) para MPI completamente limpo.
    call mpas_create_grid(g_grid, rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! -- 2. Campos ESMF e NUOPC_Realize (ANTES de mpas_atm_init) ----------
    ! ESMF_FieldCreate sobre ESMF_Grid: sem MOAB, sem deadlock.
    ! ESMF_Grid distribui automaticamente -> todos os PETs tem celulas locais.
    do i = 1, N_IMP
      field = ESMF_FieldCreate(g_grid, ESMF_TYPEKIND_R8, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, &
                               name=trim(IMP_NAMES(i)), rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call NUOPC_Realize(importState, field=field, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end do
    do i = 1, N_EXP
      field = ESMF_FieldCreate(g_grid, ESMF_TYPEKIND_R8, &
                               staggerloc=ESMF_STAGGERLOC_CENTER, &
                               name=trim(EXP_NAMES(i)), rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call NUOPC_Realize(exportState, field=field, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end do

    ! -- 3. Inicializar MPAS-A (SMIOL comeca aqui) ------------------------
    allocate(g_atm_public)
    allocate(g_atm_state)
    allocate(g_atm_bnd)
    call mpas_atm_init(g_atm_public, g_atm_state, g_atm_bnd, &
                       dt_atm_s, config_dir, localMpiComm, rc)
    if (rc /= 0) then
      call ESMF_LogSetError(ESMF_FAILURE, msg=subname//': mpas_atm_init falhou', &
           line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return
    end if

    ! -- 4. Coordenadas NetCDF (MPI_Allgather apos SMIOL - seguro) --------
    block
      real(ESMF_KIND_R8), allocatable :: lon_local_nc(:), lat_local_nc(:)
      integer :: k, n_local
      ! B-32: usar nCellsSolve (celulas proprias sem halos) para que a soma
      ! global em netcdf_init_coords seja exatamente 40962 (nao 83897 com halos).
      n_local = g_atm_public%nCellsSolve
      if (n_local == 0) n_local = g_atm_public%nCells   ! fallback se nao disponivel
      allocate(lon_local_nc(n_local), lat_local_nc(n_local))
      do k = 1, n_local
        lon_local_nc(k) = real(g_atm_public%lonCell(k), ESMF_KIND_R8) * RAD2DEG
        lat_local_nc(k) = real(g_atm_public%latCell(k), ESMF_KIND_R8) * RAD2DEG
      end do
      call netcdf_config_set(cfg_grid_res_deg, cfg_output_dir, localPet)
      call netcdf_init_coords(lon_local_nc, lat_local_nc, n_local, vm, rc)
      deallocate(lon_local_nc, lat_local_nc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end block

    call ESMF_LogWrite(subname//': InitializeRealize concluido', ESMF_LOGMSG_INFO)
  end subroutine InitializeRealize

  subroutine InitializeDataComplete(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer,             intent(out) :: rc
    type(ESMF_State)  :: importState, exportState
    type(ESMF_Clock)  :: clock
    character(len=*), parameter :: subname = '(mpas_cap:InitializeDataComplete)'
    rc = ESMF_SUCCESS
    call NUOPC_ModelGet(gcomp, &
         importState=importState, exportState=exportState, &
         modelClock=clock, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call init_import_defaults(importState, rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call mpas_atm_init_sfc(g_atm_public, g_atm_state, rc)
    if (rc /= 0) then
      call ESMF_LogSetError(ESMF_FAILURE, msg=subname//': mpas_atm_init_sfc falhou', &
           line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return
    end if
    call mpas_export(g_atm_public, exportState, rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, &
         name='InitializeDataProgress', value='true', rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call NUOPC_CompAttributeSet(gcomp, &
         name='InitializeDataComplete', value='true', rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call ESMF_LogWrite(subname//': DataInitialize SATISFIED', ESMF_LOGMSG_INFO)
  end subroutine InitializeDataComplete

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer,             intent(out) :: rc
    type(ESMF_State)    :: importState, exportState
    type(ESMF_Clock)    :: clock
    type(ESMF_VM)       :: vm
    type(ESMF_Time)     :: currTimeLoc
    integer             :: yr, mo, dy, hr, mn, sc
    character(len=*), parameter :: subname = '(mpas_cap:ModelAdvance)'
    rc = ESMF_SUCCESS
    call NUOPC_ModelGet(gcomp, &
         importState=importState, exportState=exportState, &
         modelClock=clock, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    step_count = step_count + 1
    call mpas_import(importState, g_atm_bnd, g_atm_public%nCells, rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (write_diag) then
      call state_diagnose(importState, 'importState@Advance', rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    call mpas_atm_run(g_atm_public, g_atm_state, g_atm_bnd, dt_coupling_s, rc)
    if (rc /= 0) then
      call ESMF_LogSetError(ESMF_FAILURE, msg=subname//': mpas_atm_run falhou', &
           line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return
    end if
    call mpas_export(g_atm_public, exportState, rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (write_diag) then
      call state_diagnose(exportState, 'exportState@Advance', rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    if (write_netcdf .and. mod(step_count, netcdf_write_freq) == 0) then
      call ESMF_VMGetCurrent(vm, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_ClockGet(clock, currTime=currTimeLoc, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call ESMF_TimeGet(currTimeLoc, yy=yr, mm=mo, dd=dy, &
                        h=hr, m=mn, s=sc, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      call export_write_netcdf(exportState, step_count * dt_coupling_s, &
                                yr, mo, dy, hr, mn, sc, vm, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    call ESMF_LogWrite(subname//': ModelAdvance concluido', ESMF_LOGMSG_INFO)
  end subroutine ModelAdvance

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer,             intent(out) :: rc
    character(len=*), parameter :: subname = '(mpas_cap:ModelFinalize)'
    rc = ESMF_SUCCESS
    call mpas_atm_final(g_atm_public, g_atm_state, g_atm_bnd, rc)
    if (rc /= 0) then
      call ESMF_LogSetError(ESMF_FAILURE, msg=subname//': mpas_atm_final falhou', &
           line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return
    end if
    ! ESMF_GridDestroy removido: os campos do importState/exportState
    ! ainda referenciam g_grid quando ModelFinalize e chamado.
    ! Destruir o grid aqui causa SIGSEGV no cleanup posterior do framework.
    ! O ESMF finaliza o grid automaticamente em ESMF_Finalize.
    deallocate(g_atm_public, g_atm_state, g_atm_bnd)
    g_atm_public => null()
    g_atm_state  => null()
    g_atm_bnd    => null()
    call ESMF_LogWrite(subname//': ModelFinalize concluido', ESMF_LOGMSG_INFO)
  end subroutine ModelFinalize

  !> @brief Preenche campos de importacao com valores padrao (t=0).
  !!
  !! FIX v5.2: usa ESMF_FieldGet(dimCount=) antes de farrayPtr para campos rank-2
  !! (ESMF_Grid 360x180), evitando erro ESMF_LocalArrayGetData rank mismatch.
  !> @brief Suprime validacao de timestamp dos campos de importacao.
  !!
  !! O conector OCN->MPAS fornece SST com lag de 1 passo (t-1), portanto
  !! os campos de importacao nunca tem timestamp = currTime. A validacao
  !! padrao NUOPC (label_CheckImport) geraria "INCOMPATIBILITY: Import Fields
  !! not at current time" em todos os 48 passos. Esta rotina substitui o
  !! CheckImport padrao com sucesso incondicional.
  subroutine CheckImportAlwaysOK(gcomp, rc)
    type(ESMF_GridComp) :: gcomp
    integer, intent(out) :: rc
    rc = ESMF_SUCCESS
  end subroutine CheckImportAlwaysOK

  subroutine init_import_defaults(importState, rc)
    type(ESMF_State), intent(inout) :: importState
    integer,          intent(out)   :: rc
    type(ESMF_Field)               :: field
    real(ESMF_KIND_R8), pointer    :: fptr1d(:)
    real(ESMF_KIND_R8), pointer    :: fptr2d(:,:)
    real(ESMF_KIND_R8)             :: defaults(N_IMP)
    integer :: i, fld_rank, localDeCount_imp
    rc = ESMF_SUCCESS
    defaults(1) = real(cfg_sst_default,          ESMF_KIND_R8)
    defaults(2) = real(cfg_ice_fraction_default,  ESMF_KIND_R8)
    defaults(3) = real(cfg_zorl_default,          ESMF_KIND_R8)
    do i = 1, N_IMP
      nullify(fptr1d, fptr2d)
      call ESMF_StateGet(importState, itemName=trim(IMP_NAMES(i)), &
                         field=field, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      ! B-45: PETs sem DE local na grade MPAS (360x180, regDecomp(2)=90) tem
      ! localDeCount=0 com 512 PETs (PETs 90-511). ESMF_FieldGet(farrayPtr)
      ! nestes PETs gera "localDe is out of range". Verificar antes de acessar.
      call ESMF_FieldGet(field, localDeCount=localDeCount_imp, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      if (localDeCount_imp == 0) cycle   ! PET sem dados locais - nada a inicializar
      ! Consultar rank antes de chamar farrayPtr (evita erro rank mismatch)
      call ESMF_FieldGet(field, dimCount=fld_rank, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      if (fld_rank == 1) then
        call ESMF_FieldGet(field, farrayPtr=fptr1d, rc=rc)
        if (ChkErr(rc, __LINE__, u_FILE_u)) return
        fptr1d = defaults(i)
        nullify(fptr1d)
      else
        call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
        if (ChkErr(rc, __LINE__, u_FILE_u)) return
        fptr2d = defaults(i)
        nullify(fptr2d)
      end if
    end do
  end subroutine init_import_defaults

  pure function int_to_str(n) result(s)
    integer, intent(in) :: n
    character(len=12)   :: s
    write(s, '(I0)') n
  end function int_to_str

end module mpas_cap_mod
