!==============================================================================!
! ocn_comp_NUOPC.F90 - Componente Oceânico de Dados NUOPC (DOCN)              !
!                                                                              !
! Baseado em: DOCN_cap.F90 (v1.0, GT Acoplamento MONAN / INPE/CGCT/DIMNT)    !
!                                                                              !
! Papel na arquitetura:                                                        !
!   Componente OCN alternativo ao MOM6 (mom_cap.F90). Usado em dois cenários: !
!     1. Desenvolvimento / testes: substitui MOM6 por campos prescritivos.    !
!     2. Acoplamento MPAS puro: fornece SST/gelo ao MPAS sem integrar o OCN.  !
!                                                                              !
! Campos exportados (para MED e para conector OCN?MPAS):                      !
!   So_t      temperatura da superfície do mar (SST)    [K]                   !
!   Si_ifrac  fração de gelo marinho                    [0?1]                 !
!   Sf_zorl   comprimento de rugosidade oceânica        [m]                   !
!   So_s      salinidade superficial (opcional)         [psu]                 !
!   So_u      corrente superficial zonal                [m/s]                 !
!   So_v      corrente superficial meridional           [m/s]                 !
!   So_omask  máscara oceânica (1=oceano, 0=terra)      [-]                   !
!                                                                              !
! Campos importados (do mediador MED?OCN ? recebidos mas não processados):    !
!   Foxx_taux, Foxx_tauy, Foxx_sen, Foxx_evap, Foxx_lwnet,                   !
!   Foxx_swnet_vdr, Foxx_swnet_vdf, Foxx_swnet_idr, Foxx_swnet_idf,          !
!   Faxa_rain, Faxa_snow, Sa_pslv, Si_ifrac, So_duu10n                       !
!                                                                              !
! Modos de operação (seleção via nuopc.input &nuopc_datocn):                  !
!   datocn_mode = 'netcdf'  ? lê SST de arquivo NetCDF com interp. temporal   !
!   datocn_mode = 'stub'    ? SST sintética constante (fallback/testes)       !
!                                                                              !
! Estratégia de leitura paralela:                                              !
!   PET0 lê o campo global inteiro do NetCDF e faz broadcast via              !
!   ESMF_VMBroadcast. Cada PET copia o seu subdomínio local.                  !
!   Adequado para grids até ~1440×1080 (OISST 0.25°): ~12 MB/campo/snapshot. !
!                                                                              !
! Arquivo NetCDF esperado (OISST v2.1 compatível):                            !
!   dims  : lon(1440), lat(720), time(N)                                      !
!   vars  : sst(lon,lat,time) [°C], aice(lon,lat,time) [0?1]                 !
!   Nota  : SST é convertida de °C ? K internamente (+273.15).                !
!           Se o arquivo já estiver em K, ajuste SST_CELSIUS_TO_K = 0.0.      !
!                                                                              !
! Nota sobre So_omask:                                                         !
!   O MED_cap.F90 usa So_omask para reconciliação de máscara costeira.        !
!   Em modo stub, omask = 1.0 (oceano puro ? sem informação de terra).        !
!   Em modo netcdf, omask é derivado da máscara de fill do arquivo SST        !
!   (pontos onde SST > 1e10 são tratados como terra ? omask = 0).            !
!                                                                              !
! Referência de design: DOCN_cap.F90, DATM_cap.F90 (JRA55).                  !
! Versão 1.0 ? GT Acoplamento MONAN / INPE/CGCT/DIMNT ? Maio 2026.           !
!==============================================================================!

module ocn_comp_NUOPC

  use ESMF
  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSetEntryPoint
  use ESMF, only: ESMF_GridCompGetInternalState, ESMF_GridCompSetInternalState
  use ESMF, only: ESMF_State, ESMF_StateGet
  use ESMF, only: ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet
  use ESMF, only: ESMF_Grid, ESMF_GridCreate1PeriDim, ESMF_GridAddCoord
  use ESMF, only: ESMF_GridGetCoord, ESMF_GridGet
  use ESMF, only: ESMF_Clock, ESMF_ClockGet
  use ESMF, only: ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF, only: ESMF_TimeInterval, ESMF_TimeIntervalSet, ESMF_TimeIntervalGet
  use ESMF, only: ESMF_METHOD_INITIALIZE, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_TYPEKIND_R8, ESMF_KIND_R8, ESMF_KIND_I8
  use ESMF, only: ESMF_INDEX_GLOBAL, ESMF_COORDSYS_SPH_DEG
  use ESMF, only: ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGERR_PASSTHRU
  use ESMF, only: ESMF_LogFoundError, ESMF_LogWrite, ESMF_LOGMSG_INFO
  use ESMF, only: ESMF_LOGMSG_WARNING, ESMF_LOGMSG_ERROR
  use ESMF, only: ESMF_VM, ESMF_VMGetGlobal, ESMF_VMGetCurrent
  use ESMF, only: ESMF_VMGet, ESMF_VMBroadcast
  use ESMF, only: ESMF_CALKIND_GREGORIAN

  use netcdf

  use NUOPC, only: NUOPC_CompDerive, NUOPC_CompSpecialize, NUOPC_CompSetEntryPoint
  use NUOPC, only: NUOPC_CompFilterPhaseMap, NUOPC_Advertise, NUOPC_Realize
  use NUOPC, only: NUOPC_SetTimestamp, NUOPC_CompAttributeSet
  use NUOPC_Model, &
    model_routine_SS           => SetServices,          &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_Advance        => label_Advance,        &
    model_label_Finalize       => label_Finalize,       &
    model_label_SetClock       => label_SetClock
  use NUOPC_Model, only: NUOPC_ModelGet, SetVM

  use mpas_cap_config_mod, only: cfg_sst_default,           &
                                  cfg_ice_fraction_default,  &
                                  cfg_zorl_default,          &
                                  cfg_datocn_mode,           &
                                  cfg_datocn_sst_file,       &
                                  cfg_datocn_ice_file,       &
                                  cfg_datocn_cur_file,       &
                                  cfg_datocn_nx,             &
                                  cfg_datocn_ny,             &
                                  cfg_datocn_dt_data,        &
                                  cfg_datocn_epoch_year,     &
                                  cfg_datocn_epoch_month,    &
                                  cfg_datocn_epoch_day,      &
                                  cfg_datocn_sst_varname,    &
                                  cfg_datocn_ice_varname,    &
                                  cfg_datocn_cur_u_varname,  &
                                  cfg_datocn_cur_v_varname,  &
                                  cfg_write_import_diag,     &
                                  cfg_import_diag_dir,       &
                                  cfg_datocn_ice_pct,        &
                                  cfg_grid_res_deg

  implicit none
  private

  public :: SetServices
  public :: SetVM

  !----------------------------------------------------------------------------
  ! Constantes físicas
  !----------------------------------------------------------------------------
  ! OISST v2.1 armazena SST em °C. Ajuste para 0.0 se o arquivo já for em K.
  real(ESMF_KIND_R8), parameter :: SST_CELSIUS_TO_K = 273.15_ESMF_KIND_R8

  ! Rugosidade oceânica padrão (Charnock, 1955 ? valor típico para oceano aberto)
  real(ESMF_KIND_R8), parameter :: ZORL_DEFAULT     = 0.001_ESMF_KIND_R8  ! [m]

  ! Máscara oceânica: 1.0 = oceano, 0.0 = terra
  real(ESMF_KIND_R8), parameter :: OMASK_OCEAN = 1.0_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: OMASK_LAND  = 0.0_ESMF_KIND_R8

  ! Limiar para fill values e outliers de correntes oceânicas
  real(ESMF_KIND_R8), parameter :: FILL_THRESHOLD  = 1.0e10_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: CURR_MAX_VALID   = 10.0_ESMF_KIND_R8   ! [m/s]

  !----------------------------------------------------------------------------
  ! Campos exportados (OCN ? MED e OCN ? MPAS via conector)
  !----------------------------------------------------------------------------
  integer, parameter :: N_EXP = 7
  character(len=32), parameter :: EXP_NAMES(N_EXP) = [ &
    "So_t    ", &  ! SST [K]
    "Si_ifrac", &  ! Fração de gelo [0-1]
    "Sf_zorl ", &  ! Rugosidade [m]
    "So_s    ", &  ! Salinidade superficial [psu]   (padrão 35 psu)
    "So_u    ", &  ! Corrente zonal [m/s]            (padrão 0.0)
    "So_v    ", &  ! Corrente meridional [m/s]       (padrão 0.0)
    "So_omask" ]   ! Máscara oceânica [-]            (1=oceano, 0=terra)

  !----------------------------------------------------------------------------
  ! Campos importados (MED ? OCN ? recebidos mas não integrados neste componente)
  ! Nota B-39: em modo 'stub' estes campos NÃO são anunciados para evitar que
  ! o conector MED?OCN tente criar RouteHandles bilineares que falham com DEs
  ! de largura 1 quando petCount > ny/2. Ver InitializeAdvertise.
  !----------------------------------------------------------------------------
  integer, parameter :: N_IMP = 14
  character(len=32), parameter :: IMP_NAMES(N_IMP) = [ &
    "Foxx_taux     ", "Foxx_tauy     ", "Foxx_sen      ", "Foxx_evap     ", &
    "Foxx_lwnet    ", "Foxx_swnet_vdr", "Foxx_swnet_vdf", "Foxx_swnet_idr", &
    "Foxx_swnet_idf", "Faxa_rain     ", "Faxa_snow     ", "Sa_pslv       ", &
    "Si_ifrac      ", "So_duu10n     " ]

  character(len=*), parameter :: u_FILE_u = __FILE__

  !----------------------------------------------------------------------------
  ! Estado interno do componente
  !----------------------------------------------------------------------------
  type :: OCN_DOCN_InternalState
    type(ESMF_Grid) :: grid
    ! Campos oceânicos interpolados (subdomínio local do PET)
    real(ESMF_KIND_R8), pointer :: sst(:,:)    => null()  ! SST [K]
    real(ESMF_KIND_R8), pointer :: aice(:,:)   => null()  ! fração de gelo [0-1]
    real(ESMF_KIND_R8), pointer :: sss(:,:)    => null()  ! salinidade [psu]
    real(ESMF_KIND_R8), pointer :: uocn(:,:)   => null()  ! corrente zonal [m/s]
    real(ESMF_KIND_R8), pointer :: vocn(:,:)   => null()  ! corrente meridional [m/s]
    real(ESMF_KIND_R8), pointer :: omask(:,:)  => null()  ! máscara oceânica [-]
    logical :: initialized = .false.
  end type OCN_DOCN_InternalState

  type :: OCN_DOCN_InternalStateWrapper
    type(OCN_DOCN_InternalState), pointer :: wrap => null()
  end type OCN_DOCN_InternalStateWrapper

contains

  !============================================================================
  ! SetServices ? registra fases IPDv03 e especializa ModelAdvance
  !============================================================================
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Fase 0: filtro de protocolo IPDv03
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Fase IPDv03p1: anúncio de campos
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Fase IPDv03p3: criação de grade e realização de campos
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Fase IPDv03p7: preenchimento inicial dos campos (DataInitialize)
    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_DataInitialize, &
      specRoutine=InitializeDataComplete, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Avanço: leitura/interpolação e preenchimento do exportState
    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('OCN_DOCN: SetServices concluido (modo='// &
      trim(cfg_datocn_mode)//')', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !============================================================================
  ! InitializeP0 ? filtra protocolo para IPDv03
  !============================================================================
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine InitializeP0

  !============================================================================
  ! InitializeAdvertise ? anuncia campos de SST/gelo/corrente para MED e MPAS
  !
  ! B-39: em modo 'stub' os campos de importação NÃO são anunciados para evitar
  ! que o conector NUOPC MED?OCN tente criar um RouteHandle bilinear que falha
  ! com: "not supported on Grids that contain a DE of width 1"
  ! quando petCount > ny/2.
  ! Em modo 'netcdf': campos importados são anunciados normalmente.
  !============================================================================
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    integer :: i, n_imp_active

    rc = ESMF_SUCCESS

    ! B-39: stub não usa fluxos do mediador ? suprimir importações
    if (trim(cfg_datocn_mode) == 'stub') then
      n_imp_active = 0
      call ESMF_LogWrite( &
        'OCN_DOCN stub: importacoes MED->OCN suprimidas (B-39)', &
        ESMF_LOGMSG_INFO)
    else
      n_imp_active = N_IMP
    end if

    do i = 1, n_imp_active
      call NUOPC_Advertise(importState, &
        StandardName=trim(IMP_NAMES(i)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do

    ! Campos exportados sempre anunciados (incluindo So_omask para MED)
    do i = 1, N_EXP
      call NUOPC_Advertise(exportState, &
        StandardName=trim(EXP_NAMES(i)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do

    call ESMF_LogWrite('OCN_DOCN: InitializeAdvertise concluido (' &
      //trim(int_str(N_EXP))//' exp, ' &
      //trim(int_str(n_imp_active))//' imp ativos de ' &
      //trim(int_str(N_IMP))//' disponiveis, modo=' &
      //trim(cfg_datocn_mode)//')', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !============================================================================
  ! InitializeRealize ? cria grade regular lat/lon e realiza campos
  !
  ! Grade configurável via nuopc.input (&nuopc_datocn):
  !   datocn_nx = 1440  (OISST 0.25°)  ou  360 (1.0°)
  !   datocn_ny =  720  (OISST 0.25°)  ou  180 (1.0°)
  !
  ! B-40: em modo 'stub' usa a MESMA grade do MPAS (cfg_grid_res_deg) para
  ! que o conector OCN?MPAS use ESMF_FieldBundleRedist (sem MOAB/bilinear).
  !
  ! B-57: decomposição 2D via sqrt(petCount) tiles por dimensão para evitar
  ! tiles com razão de aspecto extrema que causa deadlock no MOAB.
  !============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_Grid)   :: grid
    type(ESMF_VM)     :: vm
    integer           :: nx, ny, petCount
    real(ESMF_KIND_R8)              :: dx, dy
    real(ESMF_KIND_R8), pointer     :: coordX(:,:), coordY(:,:)
    integer                         :: localDeCount, lde
    integer                         :: i, j
    type(OCN_DOCN_InternalStateWrapper) :: iswrap
    type(OCN_DOCN_InternalState), pointer :: is

    rc = ESMF_SUCCESS

    !--------------------------------------------------------------------------
    ! B-40: seleção de resolução da grade conforme modo de operação
    !--------------------------------------------------------------------------
    if (trim(cfg_datocn_mode) == 'stub') then
      ! Grade idêntica à criada pelo mpas_create_grid (evita regridding MOAB)
      nx = nint(360.0_ESMF_KIND_R8 / real(cfg_grid_res_deg, ESMF_KIND_R8))
      ny = nint(180.0_ESMF_KIND_R8 / real(cfg_grid_res_deg, ESMF_KIND_R8))
    else
      ! Modo netcdf: resolução nativa do dado oceânico (OISST 0.25° = 1440×720)
      nx = cfg_datocn_nx
      ny = cfg_datocn_ny
    end if
    dx = 360.0_ESMF_KIND_R8 / real(nx, ESMF_KIND_R8)
    dy = 180.0_ESMF_KIND_R8 / real(ny, ESMF_KIND_R8)

    !--------------------------------------------------------------------------
    ! B-57: decomposição 2D com tiles quadradas (sqrt(petCount) por dimensão)
    !--------------------------------------------------------------------------
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    block
      integer :: nx_target, nx_max, ny_tiles, regDecomp_2d(2)
      nx_target      = max(1, nint(sqrt(real(petCount))))
      nx_max         = min(nx_target, nx / 2)
      ny_tiles       = (petCount + nx_max - 1) / nx_max
      regDecomp_2d(1) = min(nx_max, petCount)
      regDecomp_2d(2) = max(1, ny_tiles)

      grid = ESMF_GridCreate1PeriDim( &
        minIndex  = (/1, 1/),               &
        maxIndex  = (/nx, ny/),             &
        regDecomp = regDecomp_2d,           &
        indexflag = ESMF_INDEX_GLOBAL,      &
        coordSys  = ESMF_COORDSYS_SPH_DEG, &
        rc        = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end block

    call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! B-53: loop sobre DEs locais (com regDecomp 2D pode haver >1 DE/PET)
    call ESMF_GridGet(grid, localDeCount=localDeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    do lde = 0, localDeCount - 1
      call ESMF_GridGetCoord(grid, coordDim=1, localDE=lde, &
        staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordX, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      do j = lbound(coordX,2), ubound(coordX,2)
        do i = lbound(coordX,1), ubound(coordX,1)
          coordX(i,j) = (real(i,ESMF_KIND_R8) - 1.0_ESMF_KIND_R8)*dx &
                      + dx*0.5_ESMF_KIND_R8
        end do
      end do

      call ESMF_GridGetCoord(grid, coordDim=2, localDE=lde, &
        staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordY, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      do j = lbound(coordY,2), ubound(coordY,2)
        do i = lbound(coordY,1), ubound(coordY,1)
          coordY(i,j) = -90.0_ESMF_KIND_R8 &
            + (real(j,ESMF_KIND_R8) - 1.0_ESMF_KIND_R8)*dy + dy*0.5_ESMF_KIND_R8
        end do
      end do
    end do

    ! Realiza campos importados (somente em modo netcdf ? B-39)
    if (trim(cfg_datocn_mode) /= 'stub') then
      call RealizeFields(importState, grid, IMP_NAMES, N_IMP, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    ! Realiza campos exportados
    call RealizeFields(exportState, grid, EXP_NAMES, N_EXP, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Inicializa estado interno
    allocate(iswrap%wrap)
    is             => iswrap%wrap
    is%grid        = grid
    is%initialized = .false.

    call ESMF_GridCompSetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('OCN_DOCN: InitializeRealize concluido (grade ' &
      //trim(int_str(nx))//'x'//trim(int_str(ny))// &
      ', modo='//trim(cfg_datocn_mode)//')', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize

  !============================================================================
  ! InitializeDataComplete ? IPDv03p7: popula exportState e sinaliza conclusão
  !============================================================================
  subroutine InitializeDataComplete(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    type(ESMF_State)               :: exportState
    type(ESMF_Field)               :: field
    type(ESMF_Clock)               :: clock_idc
    type(ESMF_Time)                :: startTime_idc
    integer                        :: fieldCount, i
    character(len=64), allocatable :: fieldNameList(:)
    real(ESMF_KIND_R8), pointer    :: fptr(:,:)

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Preencher exportState com valores iniciais fisicamente consistentes
    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (fieldCount > 0) then
      allocate(fieldNameList(fieldCount))
      call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      do i = 1, fieldCount
        call ESMF_StateGet(exportState, itemName=trim(fieldNameList(i)), &
          field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return

        select case (trim(fieldNameList(i)))
          case ('So_t')
            ! SST padrão do namelist (&nuopc_atm_bnd) ? ex: 290 K ou 298 K
            fptr = real(cfg_sst_default, ESMF_KIND_R8)
          case ('Si_ifrac')
            fptr = real(cfg_ice_fraction_default, ESMF_KIND_R8)
          case ('Sf_zorl')
            fptr = ZORL_DEFAULT
          case ('So_s')
            fptr = 35.0_ESMF_KIND_R8    ! salinidade média global [psu]
          case ('So_u', 'So_v')
            fptr = 0.0_ESMF_KIND_R8     ! correntes em repouso
          case ('So_omask')
            fptr = OMASK_OCEAN          ! stub: oceano puro (sem máscara de terra)
          case default
            fptr = 0.0_ESMF_KIND_R8
        end select
        nullify(fptr)
      end do
      deallocate(fieldNameList)
    end if

    ! Atualizar timestamps: NUOPC_ModelBase verifica que campos do exportState
    ! estejam em currTime. Sem NUOPC_SetTimestamp o MPAS rejeita com
    ! "time mismatch" na fase InitializeIPDv03p8.
    call ESMF_GridCompGet(gcomp, clock=clock_idc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_ClockGet(clock_idc, startTime=startTime_idc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do i = 1, fieldCount
      call ESMF_StateGet(exportState, itemName=trim(fieldNameList(i)), &
        field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_SetTimestamp(field, startTime_idc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do
    deallocate(fieldNameList)

    call NUOPC_CompAttributeSet(gcomp, &
      name="InitializeDataProgress", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(gcomp, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('OCN_DOCN: InitializeDataComplete SATISFIED', &
      ESMF_LOGMSG_INFO)

  end subroutine InitializeDataComplete

  !============================================================================
  ! ModelAdvance ? lê campos oceânicos do NetCDF e popula exportState
  !
  ! Modo 'netcdf': lê SST e gelo do arquivo com interpolação temporal linear.
  !   Dados esperados (OISST v2.1 ou equivalente):
  !     sst_file: sst(lon,lat,time) em °C, dt=24h (diário)
  !     ice_file: aice(lon,lat,time) em [0-1], dt=24h (diário)
  !     cur_file: uo(lon,lat,time) e vo(lon,lat,time) em m/s (opcional)
  !
  ! Modo 'stub': mantém valores constantes do namelist (compatível com
  !   comportamento anterior do ocn_comp_NUOPC.F90 original).
  !============================================================================
  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    type(ESMF_State)         :: importState, exportState
    type(ESMF_Clock)         :: clock
    type(ESMF_Time)          :: currTime, nextTime
    type(ESMF_TimeInterval)  :: dt
    type(ESMF_Field)         :: field
    type(OCN_DOCN_InternalStateWrapper) :: iswrap
    type(OCN_DOCN_InternalState), pointer :: is
    real(ESMF_KIND_R8), pointer :: fptr(:,:)
    integer                  :: i1, i2, j1, j2
    integer                  :: year, month, day, hour, minu, sec
    integer                  :: fieldCount, k
    character(len=64), allocatable :: fieldNameList(:)
    character(len=256) :: msg

    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    is => iswrap%wrap

    call NUOPC_ModelGet(gcomp, modelClock=clock, &
      importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=dt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    nextTime = currTime + dt

    call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, &
      h=hour, m=minu, s=sec, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    write(msg,'(A,I4,5(A,I2.2))') 'OCN_DOCN: avancando para ', year, '-', &
      month, '-', day, ' ', hour, ':', minu, ':', sec
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! Obtém limites locais do subdomínio a partir do campo So_t
    call ESMF_StateGet(exportState, itemName="So_t", field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    i1 = lbound(fptr,1); i2 = ubound(fptr,1)
    j1 = lbound(fptr,2); j2 = ubound(fptr,2)
    nullify(fptr)

    ! Aloca buffers locais na primeira chamada
    if (.not. associated(is%sst)) then
      allocate(is%sst  (i1:i2, j1:j2))
      allocate(is%aice (i1:i2, j1:j2))
      allocate(is%sss  (i1:i2, j1:j2))
      allocate(is%uocn (i1:i2, j1:j2))
      allocate(is%vocn (i1:i2, j1:j2))
      allocate(is%omask(i1:i2, j1:j2))
    end if

    !--------------------------------------------------------------------------
    ! Despacho de modo
    !--------------------------------------------------------------------------
    if (trim(cfg_datocn_mode) == 'netcdf') then

      ! SST ? B-56: nome de variável configurável (padrão OISST: 'sst')
      call ReadOcnFieldInterp(gcomp, trim(cfg_datocn_sst_file), &
        trim(cfg_datocn_sst_varname), &
        currTime, cfg_datocn_nx, cfg_datocn_ny, is%sst, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="OCN_DOCN: falha ao ler SST de "//trim(cfg_datocn_sst_file), &
        line=__LINE__, file=__FILE__)) return
      is%sst = is%sst + SST_CELSIUS_TO_K

      ! Derivar máscara oceânica a partir do SST pré-conversão
      ! (fills > FILL_THRESHOLD indicam terra/gelo permanente)
      is%omask = merge(OMASK_LAND, OMASK_OCEAN, &
        abs(is%sst - SST_CELSIUS_TO_K) > FILL_THRESHOLD)

      ! Gelo marinho
      call ReadOcnFieldInterp(gcomp, trim(cfg_datocn_ice_file), &
        trim(cfg_datocn_ice_varname), &
        currTime, cfg_datocn_nx, cfg_datocn_ny, is%aice, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="OCN_DOCN: falha ao ler aice de "//trim(cfg_datocn_ice_file), &
        line=__LINE__, file=__FILE__)) return
      ! Conversão % ? fração (OISST armazena icec em percentagem 0?100)
      if (cfg_datocn_ice_pct) is%aice = is%aice / 100.0_ESMF_KIND_R8
      ! Clamping físico [0,1]
      is%aice = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, is%aice))

      ! Correntes superficiais (opcional ? arquivo pode não existir)
      if (len_trim(cfg_datocn_cur_file) > 0) then
        call ReadOcnFieldInterp(gcomp, trim(cfg_datocn_cur_file), &
          trim(cfg_datocn_cur_u_varname), &
          currTime, cfg_datocn_nx, cfg_datocn_ny, is%uocn, rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_LogWrite('OCN_DOCN: AVISO: falha uo ? corrente zonal = 0', &
            ESMF_LOGMSG_WARNING)
          is%uocn = 0.0_ESMF_KIND_R8; rc = ESMF_SUCCESS
        else
          ! Fill value OSCAR = -999 ? limiar |v| >= 10 m/s captura fills e
          ! valores fisicamente impossíveis (correntes oceânicas: 0.01?3 m/s)
          where (abs(is%uocn) >= CURR_MAX_VALID) is%uocn = 0.0_ESMF_KIND_R8
        end if

        call ReadOcnFieldInterp(gcomp, trim(cfg_datocn_cur_file), &
          trim(cfg_datocn_cur_v_varname), &
          currTime, cfg_datocn_nx, cfg_datocn_ny, is%vocn, rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_LogWrite( &
            'OCN_DOCN: AVISO: falha vo ? corrente meridional = 0', &
            ESMF_LOGMSG_WARNING)
          is%vocn = 0.0_ESMF_KIND_R8; rc = ESMF_SUCCESS
        else
          where (abs(is%vocn) >= CURR_MAX_VALID) is%vocn = 0.0_ESMF_KIND_R8
        end if
      else
        is%uocn = 0.0_ESMF_KIND_R8
        is%vocn = 0.0_ESMF_KIND_R8
      end if

      ! Salinidade: climatologia constante (sem arquivo de dado configurado)
      is%sss = 35.0_ESMF_KIND_R8

    else
      ! Modo stub: valores constantes do namelist
      is%sst   = real(cfg_sst_default,          ESMF_KIND_R8)
      is%aice  = real(cfg_ice_fraction_default,  ESMF_KIND_R8)
      is%sss   = 35.0_ESMF_KIND_R8
      is%uocn  = 0.0_ESMF_KIND_R8
      is%vocn  = 0.0_ESMF_KIND_R8
      is%omask = OMASK_OCEAN   ! stub: oceano puro
    end if

    !--------------------------------------------------------------------------
    ! Escreve campos no exportState
    !--------------------------------------------------------------------------
    call PutField(exportState, "So_t",     is%sst,   rc); if (rc /= ESMF_SUCCESS) return
    call PutField(exportState, "Si_ifrac", is%aice,  rc); if (rc /= ESMF_SUCCESS) return
    call PutField(exportState, "So_s",     is%sss,   rc); if (rc /= ESMF_SUCCESS) return
    call PutField(exportState, "So_u",     is%uocn,  rc); if (rc /= ESMF_SUCCESS) return
    call PutField(exportState, "So_v",     is%vocn,  rc); if (rc /= ESMF_SUCCESS) return
    call PutField(exportState, "So_omask", is%omask, rc); if (rc /= ESMF_SUCCESS) return
    call FillFieldConst(exportState, "Sf_zorl", ZORL_DEFAULT, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Diagnóstico opcional: escrita NetCDF dos campos preparados a cada passo.
    ! Ativado com write_import_diag=.true. em &nuopc_datocn no nuopc.input.
    if (cfg_write_import_diag) then
      call WriteDOCNDiag(gcomp, currTime, cfg_datocn_nx, cfg_datocn_ny, rc)
      if (rc /= ESMF_SUCCESS) then
        call ESMF_LogWrite( &
          'OCN_DOCN: AVISO: WriteDOCNDiag falhou ? continuando', &
          ESMF_LOGMSG_WARNING)
        rc = ESMF_SUCCESS
      end if
    end if

    !--------------------------------------------------------------------------
    ! Atualiza timestamps de todos os campos exportados para nextTime
    !--------------------------------------------------------------------------
    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do k = 1, fieldCount
      call ESMF_StateGet(exportState, itemName=trim(fieldNameList(k)), &
        field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_SetTimestamp(field, nextTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do
    deallocate(fieldNameList)

    call ESMF_LogWrite('OCN_DOCN: ModelAdvance concluido (modo='// &
      trim(cfg_datocn_mode)//')', ESMF_LOGMSG_INFO)

  end subroutine ModelAdvance

  !============================================================================
  ! ReadOcnFieldInterp ? interpolação temporal linear entre snapshots diários
  !
  ! Estratégia paralela: PET0 lê o campo global e faz broadcast via
  ! ESMF_VMBroadcast. Cada PET copia o seu subdomínio local.
  !
  ! B-54: lê ntime do arquivo para suportar arquivos anuais com ciclo correto.
  ! B-55: variáveis de escopo corrigidas para evitar "Expected INTEGER".
  ! B-59: verifica ordem de eixos para detectar (lon,lat,time) vs (lat,lon,time).
  !============================================================================
  subroutine ReadOcnFieldInterp(gcomp, filename, varname, currTime, &
                                 nx, ny, array, rc)
    type(ESMF_GridComp),  intent(in)    :: gcomp
    character(len=*),     intent(in)    :: filename
    character(len=*),     intent(in)    :: varname
    type(ESMF_Time),      intent(in)    :: currTime
    integer,              intent(in)    :: nx, ny
    real(ESMF_KIND_R8),   pointer       :: array(:,:)
    integer,              intent(out)   :: rc

    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: epochTime
    type(ESMF_TimeInterval) :: dt_since_epoch
    integer(ESMF_KIND_I8)   :: sec_since_epoch
    integer                 :: tidx0, tidx1
    integer                 :: ntime, ncid_nt, dimid_nt, nc_rc_nt
    real(ESMF_KIND_R8)      :: alpha
    integer(ESMF_KIND_I8)   :: dt_data_i8
    real(ESMF_KIND_R8)      :: f0_data(nx,ny), f1_data(nx,ny)
    real(ESMF_KIND_R8), allocatable :: buf_global(:)
    integer :: i1, i2, j1, j2, i, j, localPet
    character(len=256) :: msg

    rc         = ESMF_SUCCESS
    dt_data_i8 = int(cfg_datocn_dt_data, ESMF_KIND_I8)

    allocate(buf_global(nx*ny))

    call ESMF_VMGetGlobal(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    f0_data = 0.0_ESMF_KIND_R8
    f1_data = 0.0_ESMF_KIND_R8

    call ESMF_TimeSet(epochTime,               &
      yy=cfg_datocn_epoch_year,                &
      mm=cfg_datocn_epoch_month,               &
      dd=cfg_datocn_epoch_day,                 &
      h=0, m=0, s=0,                          &
      calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    dt_since_epoch = currTime - epochTime
    call ESMF_TimeIntervalGet(dt_since_epoch, s_i8=sec_since_epoch, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (sec_since_epoch < 0_ESMF_KIND_I8) then
      call ESMF_LogWrite( &
        'OCN_DOCN ReadOcnFieldInterp: currTime anterior ao epochTime!', &
        ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    ! B-54: ler ntime do arquivo para ciclo correto em arquivos anuais
    nc_rc_nt = nf90_open(filename, NF90_NOWRITE, ncid_nt)
    if (nc_rc_nt == NF90_NOERR) then
      nc_rc_nt = nf90_inq_dimid(ncid_nt, 'time',  dimid_nt)
      if (nc_rc_nt /= NF90_NOERR) &
        nc_rc_nt = nf90_inq_dimid(ncid_nt, 'Time',  dimid_nt)
      if (nc_rc_nt /= NF90_NOERR) &
        nc_rc_nt = nf90_inq_dimid(ncid_nt, 'TIME',  dimid_nt)
      if (nc_rc_nt == NF90_NOERR) then
        nc_rc_nt = nf90_inquire_dimension(ncid_nt, dimid_nt, len=ntime)
      else
        ntime = huge(ntime)
      end if
      nc_rc_nt = nf90_close(ncid_nt)
    else
      ntime = huge(ntime)
    end if

    tidx0 = mod(int(sec_since_epoch / real(dt_data_i8, ESMF_KIND_R8)), ntime) + 1
    tidx1 = mod(tidx0, ntime) + 1
    alpha = real(mod(sec_since_epoch, dt_data_i8), ESMF_KIND_R8) / &
            real(dt_data_i8, ESMF_KIND_R8)
    alpha = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, alpha))

    ! PET0 lê os dois snapshots e interpola
    if (localPet == 0) then
      call ReadGlobalField(filename, varname, tidx0, nx, ny, f0_data, rc)
      if (rc /= ESMF_SUCCESS) return
      call ReadGlobalField(filename, varname, tidx1, nx, ny, f1_data, rc)
      if (rc /= ESMF_SUCCESS) return
      f0_data    = f0_data + alpha * (f1_data - f0_data)
      buf_global = reshape(f0_data, [nx*ny])
    end if

    ! Broadcast do campo global interpolado para todos os PETs
    call ESMF_VMBroadcast(vm, bcstData=buf_global, count=nx*ny, &
      rootPet=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Copia subdomínio local
    i1 = lbound(array,1); i2 = ubound(array,1)
    j1 = lbound(array,2); j2 = ubound(array,2)
    do j = j1, j2
      do i = i1, i2
        array(i,j) = buf_global((j-1)*nx + i)
      end do
    end do
    deallocate(buf_global)

    ! B-55b: formato corrigido (5 strings antes do primeiro inteiro)
    write(msg,'(A,A,A,A,A,I5,A,I5,A,F6.4)') &
      'OCN_DOCN: interp ', trim(varname), ' [', trim(filename), &
      '] tidx0=', tidx0, ' tidx1=', tidx1, ' alpha=', alpha
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

  end subroutine ReadOcnFieldInterp

  !============================================================================
  ! ReadGlobalField ? lê um snapshot NetCDF global (chamado apenas em PET0)
  !
  ! B-59: verifica ordem de eixos do arquivo. Espera (lon,lat,time) em Fortran.
  ! Se dim1 /= nx, o arquivo tem ordem diferente (ex: OSCAR: lat,lon,time) ?
  ! aborta com mensagem clara indicando o comando ncpdq para corrigir.
  !============================================================================
  subroutine ReadGlobalField(filename, varname, tidx, nx, ny, array, rc)
    character(len=*),    intent(in)  :: filename
    character(len=*),    intent(in)  :: varname
    integer,             intent(in)  :: tidx
    integer,             intent(in)  :: nx, ny
    real(ESMF_KIND_R8),  intent(out) :: array(nx,ny)
    integer,             intent(out) :: rc

    integer :: ncid, varid, start(3), count_arr(3), nc_rc

    rc    = ESMF_SUCCESS
    nc_rc = nf90_open(filename, NF90_NOWRITE, ncid)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("OCN_DOCN ReadGlobalField: falha ao abrir " &
        //trim(filename)//": "//trim(nf90_strerror(nc_rc)), &
        ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    nc_rc = nf90_inq_varid(ncid, varname, varid)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("OCN_DOCN ReadGlobalField: variavel nao encontrada: " &
        //trim(varname)//" em "//trim(filename), ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
    end if

    ! B-59: verificar ordem dos eixos (lon deve ser dim1 com tamanho nx)
    block
      integer :: ndims_var, dimids(4), dim1_size, nc_rc_dim
      character(len=64) :: dim1_name
      nc_rc_dim = nf90_inquire_variable(ncid, varid, ndims=ndims_var, &
        dimids=dimids)
      if (nc_rc_dim == NF90_NOERR .and. ndims_var >= 2) then
        nc_rc_dim = nf90_inquire_dimension(ncid, dimids(1), &
          name=dim1_name, len=dim1_size)
        if (nc_rc_dim == NF90_NOERR .and. dim1_size /= nx) then
          call ESMF_LogWrite( &
            "OCN_DOCN ReadGlobalField: ERRO B-59 ? ordem de eixos incompativel!" &
            //" dim1='"//trim(dim1_name)//"' size="//trim(int_str(dim1_size)) &
            //" esperado nx="//trim(int_str(nx)) &
            //". Corrija com: ncpdq -a time,latitude,longitude " &
            //trim(filename)//" arquivo_corrigido.nc", &
            ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
        end if
      end if
    end block

    start     = [1, 1, tidx]
    count_arr = [nx, ny, 1]
    nc_rc = nf90_get_var(ncid, varid, array, start=start, count=count_arr)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("OCN_DOCN ReadGlobalField: falha ao ler " &
        //trim(varname)//": "//trim(nf90_strerror(nc_rc)), &
        ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
    end if

    nc_rc = nf90_close(ncid)

  end subroutine ReadGlobalField

  !============================================================================
  ! RealizeFields ? cria e realiza um array de campos numa ESMF_Grid
  !============================================================================
  subroutine RealizeFields(state, grid, names, n, rc)
    type(ESMF_State),  intent(inout) :: state
    type(ESMF_Grid),   intent(in)    :: grid
    character(len=32), intent(in)    :: names(:)
    integer,           intent(in)    :: n
    integer,           intent(out)   :: rc

    type(ESMF_Field) :: field
    integer          :: i

    rc = ESMF_SUCCESS
    do i = 1, n
      field = ESMF_FieldCreate(grid=grid, typekind=ESMF_TYPEKIND_R8, &
        staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(names(i)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Realize(state, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do

  end subroutine RealizeFields

  !============================================================================
  ! PutField ? copia array 2D local para campo do exportState
  !============================================================================
  subroutine PutField(state, name, array, rc)
    type(ESMF_State),    intent(inout) :: state
    character(len=*),    intent(in)    :: name
    real(ESMF_KIND_R8),  intent(in)    :: array(:,:)
    integer,             intent(out)   :: rc

    type(ESMF_Field)            :: field
    real(ESMF_KIND_R8), pointer :: fptr(:,:)

    rc = ESMF_SUCCESS
    call ESMF_StateGet(state, itemName=trim(name), field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="OCN_DOCN PutField: "//trim(name), &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    fptr = array
    nullify(fptr)

  end subroutine PutField

  !============================================================================
  ! FillFieldConst ? preenche campo do State com valor escalar constante
  !============================================================================
  subroutine FillFieldConst(state, name, value, rc)
    type(ESMF_State),    intent(inout) :: state
    character(len=*),    intent(in)    :: name
    real(ESMF_KIND_R8),  intent(in)    :: value
    integer,             intent(out)   :: rc

    type(ESMF_Field)            :: field
    real(ESMF_KIND_R8), pointer :: fptr(:,:)

    rc = ESMF_SUCCESS
    call ESMF_StateGet(state, itemName=trim(name), field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    fptr = value
    nullify(fptr)

  end subroutine FillFieldConst

  !============================================================================
  ! int_str ? converte inteiro para string (utilitário local)
  !============================================================================
  pure function int_str(n) result(s)
    integer, intent(in) :: n
    character(len=12) :: s
    write(s, '(I0)') n
  end function int_str

  !============================================================================
  ! WriteDOCNDiag ? escrita diagnóstica dos campos preparados a cada passo
  !
  ! B-58v2: PET0 relê e reinterporla diretamente do arquivo fonte (sem
  ! MPI_Gatherv, sem ambiguidade de layout espacial entre tiles de PETs).
  ! B-59b: cur_file usa tidx_cur independente do tidx do SST para suportar
  ! arquivos de corrente com ntime diferente do arquivo SST.
  !
  ! Gera: <cfg_import_diag_dir>/docn_import_YYYYMMDD_HHMMSS.nc
  ! Ativado via nuopc.input: write_import_diag = .true.
  !============================================================================
  subroutine WriteDOCNDiag(gcomp, currTime, nx, ny, rc)
    use mpi
    type(ESMF_GridComp), intent(in)  :: gcomp
    type(ESMF_Time),     intent(in)  :: currTime
    integer,             intent(in)  :: nx, ny
    integer,             intent(out) :: rc

    type(ESMF_VM)  :: vm
    type(ESMF_Time):: epochTime
    type(ESMF_TimeInterval) :: dt_since_epoch
    integer(ESMF_KIND_I8)   :: sec_since_epoch, dt_data_i8
    integer :: localPet, mpiComm, mpiErr
    integer :: yy, mm, dd, hh, mn, ss
    character(len=256) :: fname, dname
    character(len=19)  :: tstamp
    integer :: ncid_r, ncid_w, ncstat
    integer :: varid_src
    integer :: varid_sst, varid_ice, varid_u, varid_v
    integer :: varid_lat, varid_lon, dimid_lon, dimid_lat
    integer :: ntime,     dimid_nt,     tidx0,     tidx1
    integer :: ntime_i,   dimid_nt_i,   tidx0_i,   tidx1_i
    integer :: ntime_cur, dimid_nt_cur, tidx0_cur, tidx1_cur
    real(ESMF_KIND_R8) :: alpha, alpha_i, fill_val
    integer :: nlon_diag, nlat_diag, i, j
    real(ESMF_KIND_R8), allocatable :: lon_ax(:), lat_ax(:)
    real(ESMF_KIND_R8), allocatable :: f0(:,:), f1(:,:), fout(:,:)
    real(ESMF_KIND_R8), allocatable :: ice0(:,:), ice1(:,:), iceout(:,:)
    real(ESMF_KIND_R8), allocatable :: uout(:,:), vout(:,:)

    rc       = ESMF_SUCCESS
    fill_val = -9999.0_ESMF_KIND_R8

    call ESMF_VMGetCurrent(vm, rc=rc); if (rc /= ESMF_SUCCESS) return
    call ESMF_VMGet(vm, localPet=localPet, &
      mpiCommunicator=mpiComm, rc=rc)
    if (rc /= ESMF_SUCCESS) return

    ! Sincronizar todos os PETs antes da escrita do PET0
    call MPI_Barrier(mpiComm, mpiErr)
    if (localPet /= 0) return

    call ESMF_TimeGet(currTime, yy=yy, mm=mm, dd=dd, &
      h=hh, m=mn, s=ss, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    write(tstamp,'(I4.4,I2.2,I2.2,A,I2.2,I2.2,I2.2)') yy,mm,dd,'_',hh,mn,ss

    ! Recalcula tidx/alpha (mesmo algoritmo de ReadOcnFieldInterp)
    call ESMF_TimeSet(epochTime, yy=cfg_datocn_epoch_year, &
      mm=cfg_datocn_epoch_month, dd=cfg_datocn_epoch_day, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    dt_since_epoch = currTime - epochTime
    call ESMF_TimeIntervalGet(dt_since_epoch, s_i8=sec_since_epoch, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    dt_data_i8 = int(cfg_datocn_dt_data, ESMF_KIND_I8)

    ! ntime e tidx do SST
    ncstat = nf90_open(trim(cfg_datocn_sst_file), NF90_NOWRITE, ncid_r)
    if (ncstat /= NF90_NOERR) then
      call ESMF_LogWrite('OCN_DOCN WriteDOCNDiag: falha ao abrir '// &
        trim(cfg_datocn_sst_file), ESMF_LOGMSG_WARNING)
      return
    end if
    ncstat = nf90_inq_dimid(ncid_r, 'time', dimid_nt)
    if (ncstat /= NF90_NOERR) ncstat = nf90_inq_dimid(ncid_r, 'Time', dimid_nt)
    if (ncstat == NF90_NOERR) then
      ncstat = nf90_inquire_dimension(ncid_r, dimid_nt, len=ntime)
    else
      ntime = huge(ntime)
    end if
    tidx0 = mod(int(sec_since_epoch / real(dt_data_i8, ESMF_KIND_R8)), ntime) + 1
    tidx1 = mod(tidx0, ntime) + 1
    alpha = real(mod(sec_since_epoch, dt_data_i8), ESMF_KIND_R8) &
          / real(dt_data_i8, ESMF_KIND_R8)
    alpha = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, alpha))

    allocate(f0(nx,ny), f1(nx,ny), fout(nx,ny))
    allocate(uout(nx,ny), vout(nx,ny))
    uout = 0.0_ESMF_KIND_R8; vout = 0.0_ESMF_KIND_R8

    ncstat = nf90_inq_varid(ncid_r, trim(cfg_datocn_sst_varname), varid_src)
    if (ncstat == NF90_NOERR) then
      ncstat = nf90_get_var(ncid_r, varid_src, f0, &
        start=[1,1,tidx0], count=[nx,ny,1])
      ncstat = nf90_get_var(ncid_r, varid_src, f1, &
        start=[1,1,tidx1], count=[nx,ny,1])
      fout = (1.0_ESMF_KIND_R8 - alpha)*f0 + alpha*f1
      fout = fout + SST_CELSIUS_TO_K
      where (abs(f0) > FILL_THRESHOLD .or. abs(f1) > FILL_THRESHOLD) &
        fout = fill_val
    else
      fout = fill_val
    end if
    ncstat = nf90_close(ncid_r)

    ! Gelo marinho
    allocate(ice0(nx,ny), ice1(nx,ny), iceout(nx,ny))
    iceout = fill_val
    ncstat = nf90_open(trim(cfg_datocn_ice_file), NF90_NOWRITE, ncid_r)
    if (ncstat == NF90_NOERR) then
      ncstat = nf90_inq_dimid(ncid_r, 'time', dimid_nt_i)
      if (ncstat /= NF90_NOERR) &
        ncstat = nf90_inq_dimid(ncid_r, 'Time', dimid_nt_i)
      if (ncstat == NF90_NOERR) then
        ncstat = nf90_inquire_dimension(ncid_r, dimid_nt_i, len=ntime_i)
      else
        ntime_i = ntime
      end if
      tidx0_i = mod(int(sec_since_epoch / real(dt_data_i8, ESMF_KIND_R8)), &
        ntime_i) + 1
      tidx1_i = mod(tidx0_i, ntime_i) + 1
      alpha_i = alpha
      ncstat = nf90_inq_varid(ncid_r, trim(cfg_datocn_ice_varname), varid_src)
      if (ncstat == NF90_NOERR) then
        ncstat = nf90_get_var(ncid_r, varid_src, ice0, &
          start=[1,1,tidx0_i], count=[nx,ny,1])
        ncstat = nf90_get_var(ncid_r, varid_src, ice1, &
          start=[1,1,tidx1_i], count=[nx,ny,1])
        iceout = (1.0_ESMF_KIND_R8 - alpha_i)*ice0 + alpha_i*ice1
        if (cfg_datocn_ice_pct) iceout = iceout / 100.0_ESMF_KIND_R8
        iceout = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, iceout))
        where (abs(ice0) > FILL_THRESHOLD .or. abs(ice1) > FILL_THRESHOLD) &
          iceout = fill_val
      end if
      ncstat = nf90_close(ncid_r)
    end if

    ! B-59b: correntes com tidx independente (cur_file pode ter ntime=1)
    if (len_trim(cfg_datocn_cur_file) > 0) then
      ncstat = nf90_open(trim(cfg_datocn_cur_file), NF90_NOWRITE, ncid_r)
      if (ncstat == NF90_NOERR) then
        ncstat = nf90_inq_dimid(ncid_r, 'time', dimid_nt_cur)
        if (ncstat /= NF90_NOERR) &
          ncstat = nf90_inq_dimid(ncid_r, 'Time', dimid_nt_cur)
        if (ncstat == NF90_NOERR) then
          ncstat = nf90_inquire_dimension(ncid_r, dimid_nt_cur, len=ntime_cur)
        else
          ntime_cur = 1
        end if
        tidx0_cur = mod(int(sec_since_epoch / real(dt_data_i8, ESMF_KIND_R8)), &
          ntime_cur) + 1
        tidx1_cur = mod(tidx0_cur, ntime_cur) + 1

        ncstat = nf90_inq_varid(ncid_r, trim(cfg_datocn_cur_u_varname), varid_src)
        if (ncstat == NF90_NOERR) then
          ncstat = nf90_get_var(ncid_r, varid_src, f0, &
            start=[1,1,tidx0_cur], count=[nx,ny,1])
          ncstat = nf90_get_var(ncid_r, varid_src, f1, &
            start=[1,1,tidx1_cur], count=[nx,ny,1])
          uout = (1.0_ESMF_KIND_R8 - alpha)*f0 + alpha*f1
          where (abs(uout) >= CURR_MAX_VALID .or. &
                 abs(f0)   >= CURR_MAX_VALID .or. &
                 abs(f1)   >= CURR_MAX_VALID) uout = fill_val
        end if

        ncstat = nf90_inq_varid(ncid_r, trim(cfg_datocn_cur_v_varname), varid_src)
        if (ncstat == NF90_NOERR) then
          ncstat = nf90_get_var(ncid_r, varid_src, f0, &
            start=[1,1,tidx0_cur], count=[nx,ny,1])
          ncstat = nf90_get_var(ncid_r, varid_src, f1, &
            start=[1,1,tidx1_cur], count=[nx,ny,1])
          vout = (1.0_ESMF_KIND_R8 - alpha)*f0 + alpha*f1
          where (abs(vout) >= CURR_MAX_VALID .or. &
                 abs(f0)   >= CURR_MAX_VALID .or. &
                 abs(f1)   >= CURR_MAX_VALID) vout = fill_val
        end if
        ncstat = nf90_close(ncid_r)
      end if
    end if

    ! Escrita do NetCDF de diagnóstico
    nlon_diag = nx; nlat_diag = ny
    dname = trim(cfg_import_diag_dir)
    call execute_command_line('mkdir -p '//trim(dname), wait=.true.)
    fname = trim(dname)//'/docn_import_'//trim(tstamp)//'.nc'

    allocate(lon_ax(nlon_diag), lat_ax(nlat_diag))
    do i = 1, nlon_diag
      lon_ax(i) = real(i-1, ESMF_KIND_R8) * (360.0_ESMF_KIND_R8/nlon_diag)
    end do
    do j = 1, nlat_diag
      lat_ax(j) = -90.0_ESMF_KIND_R8 &
        + real(j-1, ESMF_KIND_R8) * (180.0_ESMF_KIND_R8/(nlat_diag-1))
    end do

    ncstat = nf90_create(trim(fname), NF90_CLOBBER, ncid_w)
    if (ncstat /= NF90_NOERR) then
      call ESMF_LogWrite('OCN_DOCN WriteDOCNDiag: nf90_create falhou: '// &
        trim(nf90_strerror(ncstat)), ESMF_LOGMSG_WARNING)
      goto 99
    end if

    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'Conventions',  'CF-1.8')
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'title', &
      'OCN_DOCN exportState ? SST/gelo interpolados por passo')
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'institution', &
      'INPE/CGCT/DIMNT ? GT Acoplamento MONAN')
    write(tstamp,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
      yy,'-',mm,'-',dd,'T',hh,':',mn,':',ss
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'valid_time',   trim(tstamp))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'datocn_mode',  trim(cfg_datocn_mode))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'sst_source',   trim(cfg_datocn_sst_file))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'ice_source',   trim(cfg_datocn_ice_file))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'tidx0',        tidx0)
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'alpha',        real(alpha,4))

    ncstat = nf90_def_dim(ncid_w, 'lon', nlon_diag, dimid_lon)
    ncstat = nf90_def_dim(ncid_w, 'lat', nlat_diag, dimid_lat)

    ncstat = nf90_def_var(ncid_w, 'lon', NF90_DOUBLE, [dimid_lon], varid_lon)
    ncstat = nf90_put_att(ncid_w, varid_lon, 'units', 'degrees_east')
    ncstat = nf90_def_var(ncid_w, 'lat', NF90_DOUBLE, [dimid_lat], varid_lat)
    ncstat = nf90_put_att(ncid_w, varid_lat, 'units', 'degrees_north')

    ncstat = nf90_def_var(ncid_w, 'So_t',     NF90_DOUBLE, &
      [dimid_lon,dimid_lat], varid_sst)
    ncstat = nf90_put_att(ncid_w, varid_sst, 'units', 'K')
    ncstat = nf90_put_att(ncid_w, varid_sst, '_FillValue', fill_val)

    ncstat = nf90_def_var(ncid_w, 'Si_ifrac', NF90_DOUBLE, &
      [dimid_lon,dimid_lat], varid_ice)
    ncstat = nf90_put_att(ncid_w, varid_ice, 'units', '1')
    ncstat = nf90_put_att(ncid_w, varid_ice, '_FillValue', fill_val)

    ncstat = nf90_def_var(ncid_w, 'So_u', NF90_DOUBLE, &
      [dimid_lon,dimid_lat], varid_u)
    ncstat = nf90_put_att(ncid_w, varid_u, 'units', 'm/s')
    ncstat = nf90_put_att(ncid_w, varid_u, '_FillValue', fill_val)

    ncstat = nf90_def_var(ncid_w, 'So_v', NF90_DOUBLE, &
      [dimid_lon,dimid_lat], varid_v)
    ncstat = nf90_put_att(ncid_w, varid_v, 'units', 'm/s')
    ncstat = nf90_put_att(ncid_w, varid_v, '_FillValue', fill_val)

    ncstat = nf90_enddef(ncid_w)
    ncstat = nf90_put_var(ncid_w, varid_lon, lon_ax)
    ncstat = nf90_put_var(ncid_w, varid_lat, lat_ax)
    ncstat = nf90_put_var(ncid_w, varid_sst, fout)
    ncstat = nf90_put_var(ncid_w, varid_ice, iceout)
    ncstat = nf90_put_var(ncid_w, varid_u,   uout)
    ncstat = nf90_put_var(ncid_w, varid_v,   vout)
    ncstat = nf90_close(ncid_w)

    call ESMF_LogWrite('OCN_DOCN WriteDOCNDiag: '//trim(fname)// &
      ' (tidx0='//trim(int_str(tidx0))//' B-58v2)', ESMF_LOGMSG_INFO)

    99 continue
    if (allocated(f0))    deallocate(f0, f1, fout)
    if (allocated(ice0))  deallocate(ice0, ice1, iceout)
    if (allocated(uout))  deallocate(uout, vout)
    if (allocated(lon_ax)) deallocate(lon_ax, lat_ax)

  end subroutine WriteDOCNDiag

end module ocn_comp_NUOPC
