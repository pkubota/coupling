!==============================================================================!
! DOCN_cap.F90 — Data Ocean NUOPC Component (MOM6 forçado por dados)        !
!                                                                              !
! Analogia exata com DATM_cap.F90 (Data Atmosphere / JRA55) porém para o      !
! componente oceânico: lê campos de SST e gelo marinho de arquivos NetCDF      !
! (ex: OISST v2.1 diário, OSTIA, MOM6 restart) e os exporta para o mediador  !
! MED_cap e para o cap atmosférico MPAS (condições de contorno de superfície). !
!                                                                              !
! Campos exportados (para MED e para OCN→MPAS):                               !
!   So_t      temperatura da superfície do mar (SST)    [K]                   !
!   Si_ifrac  fração de gelo marinho                    [0–1]                 !
!   Sf_zorl   comprimento de rugosidade oceânica        [m]                   !
!   So_s      salinidade superficial do mar (opcional)  [psu]                 !
!   So_u      corrente superficial zonal                [m/s]                 !
!   So_v      corrente superficial meridional           [m/s]                 !
!                                                                              !
! Campos importados (do mediador MED→OCN — recebidos mas não processados):    !
!   Foxx_taux, Foxx_tauy, Foxx_sen, Foxx_evap, Foxx_lwnet,                   !
!   Foxx_swnet_vdr, Foxx_swnet_vdf, Foxx_swnet_idr, Foxx_swnet_idf,          !
!   Faxa_rain, Faxa_snow, Sa_pslv, Si_ifrac, So_duu10n                       !
!                                                                              !
! Modos de operação (seleção via nuopc.input &nuopc_datocn):                  !
!   datocn_mode = 'netcdf'   — lê SST de arquivo NetCDF com interp. temporal  !
!   datocn_mode = 'stub'     — SST sintética (fallback, compatível v1)        !
!                                                                              !
! Estratégia de leitura paralela (igual ao DATM):                             !
!   PET0 lê o campo global inteiro do NetCDF e faz broadcast via              !
!   ESMF_VMBroadcast. Cada PET copia o seu subdomínio local.                  !
!   Adequado para grids até ~1440×1080 (OISST 0.25°): ~12 MB/campo/snapshot. !
!                                                                              !
! Arquivo NetCDF esperado (OISST v2.1 compatível):                            !
!   dims  : lon(1440), lat(720), time(N)                                      !
!   vars  : sst(lon,lat,time) [°C], aice(lon,lat,time) [0–1]                  !
!   Nota  : SST é convertida de °C → K internamente (+273.15).                !
!           Se o arquivo já estiver em K, ajuste SST_CELSIUS_TO_K = 0.0.      !
!                                                                              !
! Referência de design: DATM_cap.F90 (JRA55), AtmOcnMedPetListProto/ESMF.    !
! Versão 1.0 — GT Acoplamento MONAN / INPE/CGCT/DIMNT — Abril 2026.          !
!==============================================================================!

module DOCN_cap_mod

  use ESMF
  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSetEntryPoint
  use ESMF, only: ESMF_GridCompGetInternalState, ESMF_GridCompSetInternalState
  use ESMF, only: ESMF_State, ESMF_StateGet
  use ESMF, only: ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet
  use ESMF, only: ESMF_Grid, ESMF_GridCreate1PeriDim, ESMF_GridAddCoord, &
                  ESMF_GridGetCoord
  use ESMF, only: ESMF_Clock, ESMF_ClockGet
  use ESMF, only: ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF, only: ESMF_TimeInterval, ESMF_TimeIntervalSet, ESMF_TimeIntervalGet
  use ESMF, only: ESMF_METHOD_INITIALIZE, ESMF_STAGGERLOC_CENTER
  use ESMF, only: ESMF_TYPEKIND_R8, ESMF_KIND_R8, ESMF_KIND_I8
  use ESMF, only: ESMF_INDEX_GLOBAL, ESMF_COORDSYS_SPH_DEG
  use ESMF, only: ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGERR_PASSTHRU
  use ESMF, only: ESMF_LogFoundError, ESMF_LogWrite, ESMF_LOGMSG_INFO
  use ESMF, only: ESMF_VM, ESMF_VMGetGlobal, ESMF_VMGet, ESMF_VMBroadcast
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

  use mpas_cap_config_mod, only: cfg_sst_default,          &
                                  cfg_ice_fraction_default, &
                                  cfg_zorl_default,         &
                                  cfg_datocn_mode,          &
                                  cfg_datocn_sst_file,      &
                                  cfg_datocn_ice_file,      &
                                  cfg_datocn_cur_file,      &
                                  cfg_datocn_nx,            &
                                  cfg_datocn_ny,            &
                                  cfg_datocn_dt_data,       &
                                  cfg_datocn_epoch_year,    &
                                  cfg_datocn_epoch_month,   &
                                  cfg_datocn_epoch_day,     &
                                  cfg_datocn_sst_varname,   &
                                  cfg_datocn_ice_varname,   &
                                  cfg_datocn_cur_u_varname, &
                                  cfg_datocn_cur_v_varname, &
                                  cfg_write_import_diag,    &
                                  cfg_import_diag_dir,      &
                                  cfg_datocn_ice_pct,       &
                                  cfg_grid_res_deg

  implicit none
  private

  public :: SetServices
  public :: SetVM

  ! ── Conversão de unidades ──────────────────────────────────────────────────
  ! OISST v2.1 armazena SST em °C. Ajuste para 0.0 se o arquivo já for em K.
  real(ESMF_KIND_R8), parameter :: SST_CELSIUS_TO_K = 273.15_ESMF_KIND_R8

  ! ── Rugosidade oceânica padrão ─────────────────────────────────────────────
  real(ESMF_KIND_R8), parameter :: ZORL_DEFAULT = 0.001_ESMF_KIND_R8  ! [m]

  ! ── Campos exportados (OCN → MED e OCN → MPAS) ───────────────────────────
  integer, parameter :: N_EXP = 6
  character(len=32), parameter :: EXP_NAMES(N_EXP) = [ &
    "So_t    ", &  ! SST [K]
    "Si_ifrac", &  ! Fracao de gelo [0-1]
    "Sf_zorl ", &  ! Rugosidade [m]
    "So_s    ", &  ! Salinidade superficial [psu]  (opcional - padrao 35 psu)
    "So_u    ", &  ! Corrente zonal [m/s]          (opcional - padrao 0.0)
    "So_v    " ]   ! Corrente meridional [m/s]     (opcional - padrao 0.0)

  ! ── Campos importados (MED → OCN) ─────────────────────────────────────────
  integer, parameter :: N_IMP = 14
  character(len=32), parameter :: IMP_NAMES(N_IMP) = [ &
    "Foxx_taux     ", "Foxx_tauy     ", "Foxx_sen      ", "Foxx_evap     ", &
    "Foxx_lwnet    ", "Foxx_swnet_vdr", "Foxx_swnet_vdf", "Foxx_swnet_idr", &
    "Foxx_swnet_idf", "Faxa_rain     ", "Faxa_snow     ", "Sa_pslv       ", &
    "Si_ifrac      ", "So_duu10n     " ]

  character(len=*), parameter :: u_FILE_u = __FILE__

  !----------------------------------------------------------------------------
  ! Estado interno do DOCN
  !----------------------------------------------------------------------------
  type :: DOCN_InternalState
    type(ESMF_Grid) :: grid
    ! Campos oceânicos interpolados (subdomínio local do PET)
    real(ESMF_KIND_R8), pointer :: sst(:,:)   => null()  ! SST [K]
    real(ESMF_KIND_R8), pointer :: aice(:,:)  => null()  ! fração de gelo [0-1]
    real(ESMF_KIND_R8), pointer :: sss(:,:)   => null()  ! salinidade [psu]
    real(ESMF_KIND_R8), pointer :: uocn(:,:)  => null()  ! corrente zonal [m/s]
    real(ESMF_KIND_R8), pointer :: vocn(:,:)  => null()  ! corrente meridional [m/s]
    logical :: initialized = .false.
  end type DOCN_InternalState

  type :: DOCN_InternalStateWrapper
    type(DOCN_InternalState), pointer :: wrap => null()
  end type DOCN_InternalStateWrapper

contains

  !=============================================================================
  ! SetServices — registra fases IPDv03 e especializa ModelAdvance
  !=============================================================================
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer,              intent(out)   :: rc

    rc = ESMF_SUCCESS

    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_DataInitialize, &
      specRoutine=InitializeDataComplete, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('DOCN: SetServices concluido', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !=============================================================================
  ! InitializeP0 — filtra protocolo para IPDv03
  !=============================================================================
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer,              intent(out)   :: rc

    rc = ESMF_SUCCESS
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine InitializeP0

  !=============================================================================
  ! InitializeAdvertise — anuncia campos de SST/gelo/corrente para o MED e MPAS
  !
  ! B-39 (fix): em modo 'stub' o componente oceânico não anuncia campos de
  ! importação (fluxos do mediador). Sem anúncios no importState, o conector
  ! NUOPC MED→OCN não encontra campos a conectar e NÃO cria RouteHandles,
  ! eliminando assim a tentativa de bilinear regridding que falha com:
  !   "not supported on Grids that contain a DE of width 1"
  ! Esse erro é causado pelo grid interno do MED (fonte do regridding),
  ! cujas DEs têm largura 1 em alguma dimensão com 128 PETs. Como o stub
  ! ignora completamente todos os fluxos recebidos, suprimir os anúncios
  ! é semanticamente correto: o stub não precisa de dados do mediador.
  !
  ! Em modo 'netcdf': todos os campos importados são anunciados normalmente
  ! para que o mediador possa enviar os fluxos calculados.
  !=============================================================================
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer,              intent(out)   :: rc

    integer :: i
    integer :: n_imp_active

    rc = ESMF_SUCCESS

    ! B-39: stub nao usa fluxos do mediador — suprimir importacoes para evitar
    ! que o conector MED->OCN tente criar um RouteHandle bilinear que falha.
    ! Em modo netcdf, anunciar normalmente (DOCN precisa dos fluxos).
    if (trim(cfg_datocn_mode) == 'stub') then
      n_imp_active = 0
      call ESMF_LogWrite( &
        'DOCN stub: importacoes MED->OCN suprimidas (B-39: evita bilinear DE width 1)', &
        ESMF_LOGMSG_INFO)
    else
      n_imp_active = N_IMP
    end if

    do i = 1, n_imp_active
      call NUOPC_Advertise(importState, StandardName=trim(IMP_NAMES(i)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    end do

    ! Campos exportados para MED e para OCN->MPAS (sempre anunciados)
    do i = 1, N_EXP
      call NUOPC_Advertise(exportState, StandardName=trim(EXP_NAMES(i)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    end do

    call ESMF_LogWrite('DOCN: InitializeAdvertise concluido (' &
      //trim(adjustl(int_str(N_EXP)))//' exp, ' &
      //trim(adjustl(int_str(n_imp_active)))//' imp ativos de '&
      //trim(adjustl(int_str(N_IMP)))//' disponiveis)', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !=============================================================================
  ! InitializeRealize — cria grade regular lat/lon e realiza campos
  !
  ! Grade configurável via nuopc.input (&nuopc_datocn):
  !   datocn_nx = 1440  (OISST 0.25°)   ou  360 (1.0°)
  !   datocn_ny =  720  (OISST 0.25°)   ou  180 (1.0°)
  ! Coordenadas: lon=[0.125..359.875], lat=[-89.875..89.875] (centros de célula).
  !=============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer,              intent(out)   :: rc

    type(ESMF_Grid)   :: grid
    type(ESMF_VM)     :: vm
    integer           :: nx, ny, i, j, petCount
    real(ESMF_KIND_R8)              :: dx, dy
    real(ESMF_KIND_R8), pointer     :: coordX(:,:), coordY(:,:)
    type(DOCN_InternalStateWrapper) :: iswrap
    type(DOCN_InternalState), pointer :: is

    rc = ESMF_SUCCESS

    ! B-40 (fix): seleção da grade do DOCN conforme o modo de operação.
    !
    ! MOAB DEADLOCK com OCN→MPAS connector:
    !   A grade padrão do DOCN (1440×720) é DIFERENTE da grade de acoplamento
    !   do MPAS (360×180). O conector OCN→MPAS tentaria fazer bilinear regridding
    !   de 1440×720 → 360×180 usando MOAB (habilitado nesta instalação do ESMF 8.9.1).
    !   MOAB realiza operações MPI coletivas durante a criação dos pesos bilineares.
    !   Com 128 PETs todos entram na mesma barreira e o programa trava sem mensagem de erro.
    !   Confirmado: PET000.esmApp.log mostra last entry = "InitializeIPDvXp08 intro"
    !   após o qual ESMF_FieldBundleRegridStorePair/MOAB entraria em deadlock.
    !
    ! SOLUÇÃO (B-40):
    !   Em modo 'stub', usar a MESMA grade do MPAS (360×180 para grid_res_deg=1.0°).
    !   Quando fonte (OCN) e destino (MPAS) têm o mesmo DistGrid, o conector NUOPC
    !   usa ESMF_FieldBundleRedist (redistribuição) em vez de regridding bilinear.
    !   Redistribuição não usa MOAB → sem deadlock.
    !
    ! Para modo 'netcdf': mantém nx=1440, ny=720 (grade nativa do OISST).
    if (trim(cfg_datocn_mode) == 'stub') then
      ! Grade idêntica à criada por mpas_create_grid (cfg_grid_res_deg = 1.0° → 360×180)
      nx = nint(360.0_ESMF_KIND_R8 / real(cfg_grid_res_deg, ESMF_KIND_R8))
      ny = nint(180.0_ESMF_KIND_R8 / real(cfg_grid_res_deg, ESMF_KIND_R8))
    else
      ! Modo netcdf: resolução nativa do dado oceânico (OISST 0.25° → 1440×720)
      nx = cfg_datocn_nx
      ny = cfg_datocn_ny
    end if
    dx = 360.0_ESMF_KIND_R8 / real(nx, ESMF_KIND_R8)
    dy = 180.0_ESMF_KIND_R8 / real(ny, ESMF_KIND_R8)

    ! B-38 (fix): obter petCount para definir decomposicao explicitamente.
    ! ESMF_GridCreate1PeriDim sem regDecomp usa decomposicao default que,
    ! em ESMF 8.9.1, pode gerar DEs de largura 1 na dimensao latitudinal
    ! quando petCount > ny/2. Isso causa falha no regridding bilinear
    ! do conector MED->OCN com o erro:
    !   "not supported on Grids that contain a DE of width 1"
    ! Fix: forcar decomposicao 1D exclusivamente na dimensao longitudinal
    ! (dim 1 = periodica): cada PET recebe nx/petCount colunas e TODAS
    ! as ny linhas de latitude. Com nx=1440 e petCount=128:
    !   1440/128 = 11.25 -> min 11 colunas/PET >> 1 → OK para bilinear.
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_VMGet(vm, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! B-57 (fix B-46/B-52): regDecomp 2D com tiles quadradas — evita strips extremos.
    !
    ! PROBLEMA com B-46 (regDecomp=(/1, min(petCount,ny/2)/)):
    !   Decompõe em latitude → PETs vazios com petCount>ny/2.
    ! PROBLEMA com B-52 (nx_max=nx/2):
    !   netcdf 1440×720 a 512 PETs → regDecomp=(/512,1/) → 2-3 cols×720 rows
    !   → aspecto 256:1 → MOAB trava em ESMF_FieldBundleRegridStore.
    !
    ! SOLUCAO B-57: sqrt(petCount) tiles por dimensão.
    !   nx_tiles_target = nint(sqrt(N)) → aspecto ≈ 1.
    !   nx_max = min(target, nx/2) → garante col ≥ 2.
    !
    !   N=4:   sqrt=2  → nx_max=2   regDecomp=(/2,2/)=4    aspecto 0.5:1 ✓
    !   N=128: sqrt=11 → nx_max=11  regDecomp=(/11,12/)=132 aspecto 0.5:1 ✓
    !   N=512: stub(360×180)  → regDecomp=(/23,23/)=529  15col× 7row ✓
    !   N=512: netcdf(1440×720)→ regDecomp=(/23,23/)=529  62col×31row ✓

    block
      integer :: nx_tiles_target, nx_max, ny_tiles, regDecomp_2d(2)
      integer :: localDeCount_docn, lde_docn
      nx_tiles_target = max(1, nint(sqrt(real(petCount))))
      nx_max          = min(nx_tiles_target, nx / 2)
      ny_tiles        = (petCount + nx_max - 1) / nx_max
      regDecomp_2d(1) = min(nx_max, petCount)
      regDecomp_2d(2) = max(1, ny_tiles)

      grid = ESMF_GridCreate1PeriDim( &
        minIndex  = (/1, 1/),                &
        maxIndex  = (/nx, ny/),              &
        regDecomp = regDecomp_2d,            &
        indexflag = ESMF_INDEX_GLOBAL,       &
        coordSys  = ESMF_COORDSYS_SPH_DEG,  &
        rc        = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      ! B-53: loop sobre DEs locais — com regDecomp 2D e DEs>petCount,
      ! 17 PETs a 512 PETs têm localDeCount=2; GridGetCoord exige localDE=.
      call ESMF_GridGet(grid, localDeCount=localDeCount_docn, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      do lde_docn = 0, localDeCount_docn - 1
        call ESMF_GridGetCoord(grid, coordDim=1, localDE=lde_docn, &
          staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordX, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do j = lbound(coordX,2), ubound(coordX,2)
          do i = lbound(coordX,1), ubound(coordX,1)
            coordX(i,j) = (real(i,ESMF_KIND_R8) - 1.0_ESMF_KIND_R8)*dx + dx*0.5_ESMF_KIND_R8
          end do
        end do
        call ESMF_GridGetCoord(grid, coordDim=2, localDE=lde_docn, &
          staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordY, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        do j = lbound(coordY,2), ubound(coordY,2)
          do i = lbound(coordY,1), ubound(coordY,1)
            coordY(i,j) = -90.0_ESMF_KIND_R8 &
              + (real(j,ESMF_KIND_R8) - 1.0_ESMF_KIND_R8)*dy + dy*0.5_ESMF_KIND_R8
          end do
        end do
      end do  ! lde_docn
    end block

    ! Importados — somente em modo netcdf (B-39: stub nao anuncia importacoes)
    if (trim(cfg_datocn_mode) /= 'stub') then
      call RealizeFields(importState, grid, IMP_NAMES, N_IMP, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    ! Exportados
    call RealizeFields(exportState, grid, EXP_NAMES, N_EXP, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Estado interno
    allocate(iswrap%wrap)
    is             => iswrap%wrap
    is%grid        = grid
    is%initialized = .false.

    call ESMF_GridCompSetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('DOCN: InitializeRealize concluido (grade ' &
      //trim(adjustl(int_str(nx)))//'x' &
      //trim(adjustl(int_str(ny)))//')', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize

  !=============================================================================
  ! InitializeDataComplete — IPDv03p7: popula exportState e sinaliza conclusao
  !=============================================================================
  subroutine InitializeDataComplete(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer,              intent(out)   :: rc

    type(ESMF_State)               :: exportState
    type(ESMF_Field)               :: field
    type(ESMF_Clock)               :: clock_idc
    type(ESMF_Time)                :: startTime_idc
    integer                        :: fieldCount, i
    integer                        :: fieldCount_ts
    character(len=64), allocatable :: fieldNameList(:)
    character(len=64), allocatable :: fldNames_ts(:)
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
            ! SST default do namelist (&nuopc_atm_bnd) ou 290 K se stub
            fptr = real(cfg_sst_default, ESMF_KIND_R8)
          case ('Si_ifrac')
            fptr = real(cfg_ice_fraction_default, ESMF_KIND_R8)
          case ('Sf_zorl')
            fptr = ZORL_DEFAULT
          case ('So_s')
            fptr = 35.0_ESMF_KIND_R8   ! salinidade media global [psu]
          case ('So_u', 'So_v')
            fptr = 0.0_ESMF_KIND_R8    ! correntes em repouso
          case default
            fptr = 0.0_ESMF_KIND_R8
        end select
        nullify(fptr)
      end do
      deallocate(fieldNameList)
    end if

    ! Atualizar timestamps: NUOPC_ModelBase verifica que os campos no
    ! importState do componente seguinte estejam no currTime do clock.
    ! Sem NUOPC_SetTimestamp os campos ficam em t=0 e o MPAS rejeita.
    call ESMF_GridCompGet(gcomp, clock=clock_idc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    ! NUOPC_SetTimestamp recebe ESMF_Time, nao ESMF_Clock
    call ESMF_ClockGet(clock_idc, startTime=startTime_idc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_StateGet(exportState, itemCount=fieldCount_ts, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    allocate(fldNames_ts(fieldCount_ts))
    call ESMF_StateGet(exportState, itemNameList=fldNames_ts, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do i = 1, fieldCount_ts
      call ESMF_StateGet(exportState, itemName=trim(fldNames_ts(i)), &
        field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_SetTimestamp(field, startTime_idc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do
    deallocate(fldNames_ts)

    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataProgress", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete",  value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('DOCN: InitializeDataComplete SATISFIED', ESMF_LOGMSG_INFO)

  end subroutine InitializeDataComplete

  !=============================================================================
  ! ModelAdvance — le campos oceanicos do NetCDF e popula exportState
  !
  ! Modo 'netcdf': le SST e gelo do arquivo com interpolacao temporal linear.
  ! Modo 'stub':   mantem valores constantes (compatibilidade com mom_cap.F90).
  !
  ! Dados esperados (OISST v2.1 ou equivalente):
  !   sst_file: sst(lon,lat,time) em °C, dt=24h (diario)
  !   ice_file: aice(lon,lat,time) em [0-1], dt=24h (diario)
  !   cur_file: uo(lon,lat,time) e vo(lon,lat,time) em m/s (opcional)
  !=============================================================================
  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer,              intent(out)   :: rc

    type(ESMF_State)         :: importState, exportState
    type(ESMF_Clock)         :: clock
    type(ESMF_Time)          :: currTime, nextTime
    type(ESMF_TimeInterval)  :: dt
    type(ESMF_Field)         :: field
    type(DOCN_InternalStateWrapper) :: iswrap
    type(DOCN_InternalState), pointer :: is
    real(ESMF_KIND_R8), pointer :: fptr(:,:)
    integer                  :: i, j, i1, i2, j1, j2
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

    write(msg,'(A,I4,5(A,I2.2))') 'DOCN: avancando para ', year, '-', &
      month, '-', day, ' ', hour, ':', minu, ':', sec
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! Obtem limites locais do subdominio a partir do primeiro campo exportado
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
    end if

    ! ── Despacho de modo ──────────────────────────────────────────────────────
    if (trim(cfg_datocn_mode) == 'netcdf') then

      ! B-56: nomes de variavel configuráveis via nuopc.input (datocn_*_varname).
      ! OISST v2.1: sst_varname='sst'  ice_varname='icec'
      ! Sem B-56: 'aice' hardcoded → "variavel nao encontrada: aice" (OISST).
      call ReadOcnFieldInterp(gcomp, trim(cfg_datocn_sst_file), &
        trim(cfg_datocn_sst_varname), &
        currTime, cfg_datocn_nx, cfg_datocn_ny, is%sst,  rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg="DOCN: falha ao ler SST", &
        line=__LINE__, file=__FILE__)) return
      ! Conversao °C → K (OISST armazena em °C)
      is%sst = is%sst + SST_CELSIUS_TO_K

      call ReadOcnFieldInterp(gcomp, trim(cfg_datocn_ice_file), &
        trim(cfg_datocn_ice_varname), &
        currTime, cfg_datocn_nx, cfg_datocn_ny, is%aice, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg="DOCN: falha ao ler aice", &
        line=__LINE__, file=__FILE__)) return
      ! Conversão % → fração: OISST armazena icec em percentagem (0–100).
      ! cfg_datocn_ice_pct=.true. (padrão OISST); .false. para produtos em fração.
      if (cfg_datocn_ice_pct) is%aice = is%aice / 100.0_ESMF_KIND_R8
      ! Clamping fisico: fracao de gelo em [0,1]
      is%aice = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, is%aice))

      ! Correntes superficiais (opcional: arquivos podem nao existir)
      if (len_trim(cfg_datocn_cur_file) > 0) then
        call ReadOcnFieldInterp(gcomp, trim(cfg_datocn_cur_file), &
          trim(cfg_datocn_cur_u_varname), &
          currTime, cfg_datocn_nx, cfg_datocn_ny, is%uocn, rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_LogWrite('DOCN: AVISO: falha uo — corrente zonal = 0', &
            ESMF_LOGMSG_INFO)
          is%uocn = 0.0_ESMF_KIND_R8; rc = ESMF_SUCCESS
        else
          ! Fill value OSCAR = -999.0 — limiar |v|>10 m/s captura fills oceânicos
          ! e valores fisicamente impossíveis (correntes reais: 0.01–3 m/s).
          where (abs(is%uocn) >= 10.0_ESMF_KIND_R8) is%uocn = 0.0_ESMF_KIND_R8
        end if
        call ReadOcnFieldInterp(gcomp, trim(cfg_datocn_cur_file), &
          trim(cfg_datocn_cur_v_varname), &
          currTime, cfg_datocn_nx, cfg_datocn_ny, is%vocn, rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_LogWrite('DOCN: AVISO: falha vo — corrente meridional = 0', &
            ESMF_LOGMSG_INFO)
          is%vocn = 0.0_ESMF_KIND_R8; rc = ESMF_SUCCESS
        else
          where (abs(is%vocn) >= 10.0_ESMF_KIND_R8) is%vocn = 0.0_ESMF_KIND_R8
        end if
      else
        is%uocn = 0.0_ESMF_KIND_R8
        is%vocn = 0.0_ESMF_KIND_R8
      end if

      ! Salinidade: sem arquivo de dado, usar climatologia constante
      is%sss = 35.0_ESMF_KIND_R8

    else
      ! Modo stub (fallback): valores constantes do namelist
      is%sst  = real(cfg_sst_default,          ESMF_KIND_R8)
      is%aice = real(cfg_ice_fraction_default,  ESMF_KIND_R8)
      is%sss  = 35.0_ESMF_KIND_R8
      is%uocn = 0.0_ESMF_KIND_R8
      is%vocn = 0.0_ESMF_KIND_R8
    end if

    ! Escreve campos no exportState
    call PutField(exportState, "So_t",    is%sst,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Si_ifrac",is%aice, rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "So_s",    is%sss,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "So_u",    is%uocn, rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "So_v",    is%vocn, rc); if (rc/=ESMF_SUCCESS) return

    ! Diagnóstico: escrita NetCDF dos campos lidos/preparados a cada passo.
    ! Ativado com write_import_diag=.true. em &nuopc_datocn no nuopc.input.
    ! Gera: diag_import/docn_import_YYYYMMDD_HHMMSS.nc (grade DOCN, 1°×1°)
    ! Lido por: postproc_monan_import.py  (validação de SST/gelo vs fonte)
    if (cfg_write_import_diag) then
      call WriteDOCNDiag(gcomp, currTime, cfg_datocn_nx, cfg_datocn_ny, rc)
      if (rc /= ESMF_SUCCESS) then
        call ESMF_LogWrite('DOCN: AVISO: WriteDOCNDiag falhou — continuando', &
          ESMF_LOGMSG_WARNING)
        rc = ESMF_SUCCESS
      end if
    end if

    ! Sf_zorl: rugosidade constante (funcao de amplitude de onda nao modelada aqui)
    call FillFieldConst(exportState, "Sf_zorl", ZORL_DEFAULT, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Atualizar timestamps de todos os campos exportados
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

    call ESMF_LogWrite('DOCN: ModelAdvance concluido (modo='// &
      trim(cfg_datocn_mode)//')', ESMF_LOGMSG_INFO)

  end subroutine ModelAdvance

  !=============================================================================
  ! ReadOcnFieldInterp — interpolacao temporal linear entre snapshots diarios
  !
  ! Identica em estrutura a ReadJRAFieldInterp do DATM_cap.F90.
  ! Estrategia paralela: PET0 le campo global e broadcast via ESMF_VMBroadcast.
  !
  ! Parametros de epoch e dt_data configurados em &nuopc_datocn:
  !   datocn_epoch_year, datocn_epoch_month, datocn_epoch_day
  !   datocn_dt_data  (segundos entre snapshots; 86400 para diario)
  !=============================================================================
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
    integer                 :: ntime, ncid_nt, dimid_nt, nc_rc_nt  ! B-55 (fix B-54 scope)
    real(ESMF_KIND_R8)      :: alpha
    integer(ESMF_KIND_I8)   :: dt_data_i8
    real(ESMF_KIND_R8)      :: f0_data(nx,ny), f1_data(nx,ny)
    real(ESMF_KIND_R8), allocatable :: buf_global(:)
    integer :: i1, i2, j1, j2, i, j, localPet
    character(len=256) :: msg

    rc      = ESMF_SUCCESS
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

    ! Epoch: primeiro snapshot do arquivo
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
      call ESMF_LogWrite('DOCN ReadOcnFieldInterp: currTime anterior ao ' &
        //'epochTime do dado oceanico!', ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; return
    end if

    ! Indices de interpolacao temporal (base-1)
    ! B-54: ler ntime do arquivo para suportar arquivos anuais (ex: OISST 2026
    !   com 365 registros). Sem isso, tidx0=16281 >> 365 → "Index exceeds
    !   dimension bound". Com ntime: tidx0 = mod(dias_epoch, ntime) + 1.
    ! B-55a (fix B-54 scope): variaveis ntime/ncid_nt/dimid_nt/nc_rc_nt
    !   declaradas na subrotina (nao mais em block construct) para que
    !   tidx0/tidx1 calculados aqui sejam visiveis no write(msg,...) abaixo.
    nc_rc_nt = nf90_open(filename, NF90_NOWRITE, ncid_nt)
    if (nc_rc_nt == NF90_NOERR) then
      nc_rc_nt = nf90_inq_dimid(ncid_nt, 'time', dimid_nt)
      if (nc_rc_nt /= NF90_NOERR) &
        nc_rc_nt = nf90_inq_dimid(ncid_nt, 'Time', dimid_nt)
      if (nc_rc_nt /= NF90_NOERR) &
        nc_rc_nt = nf90_inq_dimid(ncid_nt, 'TIME', dimid_nt)
      if (nc_rc_nt == NF90_NOERR) then
        nc_rc_nt = nf90_inquire_dimension(ncid_nt, dimid_nt, len=ntime)
      else
        ntime = huge(ntime)   ! dim nao encontrada: sem clamping
      end if
      nc_rc_nt = nf90_close(ncid_nt)
    else
      ntime = huge(ntime)     ! arquivo nao abriu: ReadGlobalField reportara
    end if
    tidx0 = mod(int(sec_since_epoch / real(dt_data_i8, ESMF_KIND_R8)), ntime) + 1
    tidx1 = mod(tidx0, ntime) + 1   ! ciclo: ultimo registro volta ao 1
    alpha = real(mod(sec_since_epoch, dt_data_i8), ESMF_KIND_R8) / &
            real(dt_data_i8, ESMF_KIND_R8)
    alpha = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, alpha))

    ! PET0 le os dois snapshots e interpola
    if (localPet == 0) then
      call ReadGlobalField(filename, varname, tidx0, nx, ny, f0_data, rc)
      if (rc /= ESMF_SUCCESS) return
      call ReadGlobalField(filename, varname, tidx1, nx, ny, f1_data, rc)
      if (rc /= ESMF_SUCCESS) return

      ! Interpolacao temporal linear in-place
      f0_data = f0_data + alpha * (f1_data - f0_data)
      buf_global = reshape(f0_data, [nx*ny])
    end if

    ! Broadcast do campo global interpolado
    call ESMF_VMBroadcast(vm, bcstData=buf_global, count=nx*ny, rootPet=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Cada PET copia o seu subdominio local
    i1 = lbound(array,1); i2 = ubound(array,1)
    j1 = lbound(array,2); j2 = ubound(array,2)
    do j = j1, j2
      do i = i1, i2
        array(i,j) = buf_global((j-1)*nx + i)
      end do
    end do

    deallocate(buf_global)

    ! B-55b (fix): formato original (A,A,A,A,I5,A,I5,A,F6.4) tinha 4 A antes
    ! do primeiro I5, mas a chamada passa 5 strings antes de tidx0 →
    ! '] tidx0=' caía no I5 → "Expected INTEGER, got CHARACTER".
    ! Correto: (A,A,A,A,A,I5,A,I5,A,F6.4) — 5 A + I5 + A + I5 + A + F6.4.
    write(msg,'(A,A,A,A,A,I5,A,I5,A,F6.4)') &
      'DOCN: interp ', trim(varname), ' [', trim(filename), &
      '] tidx0=', tidx0, ' tidx1=', tidx1, ' alpha=', alpha
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

  end subroutine ReadOcnFieldInterp

  !=============================================================================
  ! ReadGlobalField — le um snapshot NetCDF global (chamado apenas em PET0)
  !=============================================================================
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
      call ESMF_LogWrite("ReadGlobalField DOCN: falha ao abrir " &
        //trim(filename)//": "//trim(nf90_strerror(nc_rc)), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; return
    end if

    nc_rc = nf90_inq_varid(ncid, varname, varid)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("ReadGlobalField DOCN: variavel nao encontrada: " &
        //trim(varname), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
    end if

    ! B-59: verificar ordem dos eixos do arquivo NetCDF.
    ! DOCN espera (lon, lat, time) em ordem Fortran = (time, lat, lon) em C/NetCDF.
    ! OSCAR NRT usa (time, lon, lat) em C = (lat, lon, time) em Fortran.
    ! → ReadGlobalField leria count=[nx,ny,1]=[1440,720,1] em (lat,lon,time):
    !   dim1=lat tem 719 pontos → count[1]=1440 > 719 → leitura fora dos limites.
    ! Verificação: se dim1_size ≠ nx → eixos incompatíveis → abortar com mensagem clara.
    block
      integer :: ndims_var, dimids(4), dim1_size, nc_rc_dim
      character(len=64) :: dim1_name
      nc_rc_dim = nf90_inquire_variable(ncid, varid, ndims=ndims_var, dimids=dimids)
      if (nc_rc_dim == NF90_NOERR .and. ndims_var >= 2) then
        nc_rc_dim = nf90_inquire_dimension(ncid, dimids(1), name=dim1_name, len=dim1_size)
        if (nc_rc_dim == NF90_NOERR .and. dim1_size /= nx) then
          call ESMF_LogWrite( &
            "ReadGlobalField DOCN: ERRO B-59 — ordem de eixos incompativel! "// &
            "Arquivo "//trim(filename)//" tem dim1='"//trim(dim1_name)// &
            "' com tamanho "//trim(adjustl(int2str_rg(dim1_size)))// &
            " mas DOCN espera nx="//trim(adjustl(int2str_rg(nx)))//". "// &
            "Execute prepare_cur_file.sh para transpor: "// &
            "ncpdq -a time,latitude,longitude arquivo.nc arquivo_corrigido.nc", &
            ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          nc_rc = nf90_close(ncid)
          return
        end if
      end if
    end block

    start     = [1, 1, tidx]
    count_arr = [nx, ny, 1]
    nc_rc     = nf90_get_var(ncid, varid, array, start=start, count=count_arr)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("ReadGlobalField DOCN: falha ao ler " &
        //trim(varname)//": "//trim(nf90_strerror(nc_rc)), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
    end if

    nc_rc = nf90_close(ncid)

  contains
    function int2str_rg(n) result(s)
      integer, intent(in) :: n
      character(len=12) :: s
      write(s,'(I12)') n
    end function int2str_rg

  end subroutine ReadGlobalField

  !=============================================================================
  ! RealizeFields — cria e realiza um array de campos numa ESMF_Grid
  !=============================================================================
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

  !=============================================================================
  ! PutField — copia array 2D local para campo do exportState
  !=============================================================================
  subroutine PutField(state, name, array, rc)
    type(ESMF_State),    intent(inout) :: state
    character(len=*),    intent(in)    :: name
    real(ESMF_KIND_R8),  intent(in)    :: array(:,:)
    integer,             intent(out)   :: rc

    type(ESMF_Field)            :: field
    real(ESMF_KIND_R8), pointer :: fptr(:,:)

    rc = ESMF_SUCCESS
    call ESMF_StateGet(state, itemName=trim(name), field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="PutField DOCN: "//trim(name), &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    fptr = array
    nullify(fptr)

  end subroutine PutField

  !=============================================================================
  ! FillFieldConst — preenche campo do State com valor escalar constante
  !=============================================================================
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

  !=============================================================================
  ! int_str — converte inteiro para string (utilitario local)
  !=============================================================================
  pure function int_str(n) result(s)
    integer, intent(in) :: n
    character(len=12) :: s
    write(s, '(I0)') n
  end function int_str

  !=============================================================================
  ! WriteDOCNDiag — escrita diagnóstica dos campos lidos pelo DOCN a cada passo
  !
  ! Gera: <cfg_import_diag_dir>/docn_import_YYYYMMDD_HHMMSS.nc
  ! Grade: grade DOCN regulada em 1°×1° (nx=cfg_datocn_nx, ny=cfg_datocn_ny)
  !        interpolada para 360×180 (cfg_grid_res_deg) por coleta de PET0.
  ! Ativado via nuopc.input: write_import_diag = .true.
  !
  ! Lido por postproc_monan_import.py para validação de:
  !   - Interpolação temporal (comparar vs arquivos fonte OISST)
  !   - Conversão de unidades (°C→K verificado por So_t_raw)
  !   - Clamping de gelo [0,1]
  !   - Correntes (So_u, So_v)
  !=============================================================================
  !=============================================================================
  ! WriteDOCNDiag v2 — PET0 rele e reinterporla diretamente do arquivo fonte.
  !
  ! PROBLEMA da versão anterior (MPI_Gatherv):
  !   MPI_Gatherv concatena os dados locais de cada PET em ordem de PE, não em
  !   ordem espacial (lon,lat). reshape(recvbuf, [nx,ny]) reconstrói um array
  !   geograficamente incorreto (dados embaralhados entre PE tiles).
  !   Além disso, sst_glob(1:360,1:181) recortava apenas os primeiros 360 de
  !   1440 colunas do OISST, cobrindo apenas ~25% do domínio global.
  !
  ! SOLUÇÃO (B-58 v2): PET0 rele e reinterporla o campo global diretamente do
  !   arquivo NetCDF fonte, usando o mesmo algoritmo B-54. Sem MPI, sem gather,
  !   sem ambiguidade de layout. O resultado é exatamente o campo global que o
  !   DOCN calculou — validação perfeita da interpolação temporal.
  !   Outros PETs executam MPI_Barrier e retornam.
  !=============================================================================
  subroutine WriteDOCNDiag(gcomp, currTime, nx, ny, rc)
    use netcdf
    use mpi
    type(ESMF_GridComp),  intent(in)  :: gcomp
    type(ESMF_Time),      intent(in)  :: currTime
    integer,              intent(in)  :: nx, ny
    integer,              intent(out) :: rc

    type(ESMF_VM)  :: vm
    type(ESMF_Time):: epochTime
    type(ESMF_TimeInterval) :: dt_since_epoch
    integer(ESMF_KIND_I8)   :: sec_since_epoch, dt_data_i8
    integer :: localPet, mpiComm, mpiErr
    integer :: yy, mm, dd, hh, mn, ss
    character(len=256) :: fname, dname
    character(len=19)  :: tstamp
    integer :: ncid_r, ncid_w
    integer :: varid_src, ncstat
    integer :: varid_sst, varid_ice, varid_u, varid_v
    integer :: varid_lat, varid_lon, dimid_lon, dimid_lat
    integer :: ntime, dimid_nt, tidx0, tidx1
    integer :: ntime_i, dimid_nt_i, tidx0_i, tidx1_i
    integer :: ntime_cur, dimid_nt_cur, tidx0_cur, tidx1_cur  ! B-59b: cur_file independente
    real(ESMF_KIND_R8) :: alpha, alpha_i, fill_val
    integer :: nlon_diag, nlat_diag, i, j
    real(ESMF_KIND_R8), allocatable :: lon_ax(:), lat_ax(:)
    real(ESMF_KIND_R8), allocatable :: f0(:,:), f1(:,:), fout(:,:)
    real(ESMF_KIND_R8), allocatable :: ice0(:,:), ice1(:,:), iceout(:,:)
    real(ESMF_KIND_R8), allocatable :: uout(:,:), vout(:,:)

    rc = ESMF_SUCCESS
    fill_val = -9999.0_ESMF_KIND_R8

    call ESMF_VMGetCurrent(vm, rc=rc); if (rc /= ESMF_SUCCESS) return
    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=mpiComm, rc=rc)
    if (rc /= ESMF_SUCCESS) return

    ! Sincronizar — todos os PETs chegam aqui antes da escrita do PET0
    call MPI_Barrier(mpiComm, mpiErr)
    if (localPet /= 0) return   ! apenas PET0 executa o restante

    call ESMF_TimeGet(currTime, yy=yy, mm=mm, dd=dd, h=hh, m=mn, s=ss, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    write(tstamp,'(I4.4,I2.2,I2.2,A,I2.2,I2.2,I2.2)') yy,mm,dd,'_',hh,mn,ss

    ! ── Calcular tidx0, tidx1, alpha (mesmo algoritmo de ReadOcnFieldInterp/B-54) ──
    call ESMF_TimeSet(epochTime, yy=cfg_datocn_epoch_year, &
      mm=cfg_datocn_epoch_month, dd=cfg_datocn_epoch_day, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    dt_since_epoch = currTime - epochTime
    call ESMF_TimeIntervalGet(dt_since_epoch, s_i8=sec_since_epoch, rc=rc)
    if (rc /= ESMF_SUCCESS) return
    dt_data_i8 = int(cfg_datocn_dt_data, ESMF_KIND_I8)

    ! ntime da SST
    ncstat = nf90_open(trim(cfg_datocn_sst_file), NF90_NOWRITE, ncid_r)
    if (ncstat /= NF90_NOERR) then
      call ESMF_LogWrite('WriteDOCNDiag: falha ao abrir '// &
        trim(cfg_datocn_sst_file)//': '//trim(nf90_strerror(ncstat)), &
        ESMF_LOGMSG_WARNING)
      return
    end if
    ncstat = nf90_inq_dimid(ncid_r, 'time', dimid_nt)
    if (ncstat /= NF90_NOERR) ncstat = nf90_inq_dimid(ncid_r, 'Time', dimid_nt)
    if (ncstat /= NF90_NOERR) ncstat = nf90_inq_dimid(ncid_r, 'TIME', dimid_nt)
    if (ncstat == NF90_NOERR) then
      ncstat = nf90_inquire_dimension(ncid_r, dimid_nt, len=ntime)
    else
      ntime = huge(ntime)
    end if
    tidx0 = mod(int(sec_since_epoch / real(dt_data_i8, ESMF_KIND_R8)), ntime) + 1
    tidx1 = mod(tidx0, ntime) + 1
    alpha = real(mod(sec_since_epoch, dt_data_i8), ESMF_KIND_R8) / real(dt_data_i8, ESMF_KIND_R8)
    alpha = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, alpha))

    ! Ler dois snapshots SST e interpolar
    allocate(f0(nx,ny), f1(nx,ny), fout(nx,ny))
    allocate(uout(nx,ny), vout(nx,ny))
    uout = 0.0_ESMF_KIND_R8; vout = 0.0_ESMF_KIND_R8
    ncstat = nf90_inq_varid(ncid_r, trim(cfg_datocn_sst_varname), varid_src)
    if (ncstat == NF90_NOERR) then
      ncstat = nf90_get_var(ncid_r, varid_src, f0, start=[1,1,tidx0], count=[nx,ny,1])
      ncstat = nf90_get_var(ncid_r, varid_src, f1, start=[1,1,tidx1], count=[nx,ny,1])
      fout = (1.0_ESMF_KIND_R8 - alpha)*f0 + alpha*f1
      fout = fout + 273.15_ESMF_KIND_R8   ! conversao °C → K
      ! Mascarar terra
      where (abs(f0) > 1.0e10_ESMF_KIND_R8 .or. abs(f1) > 1.0e10_ESMF_KIND_R8) &
        fout = fill_val
    else
      fout = fill_val
    end if
    ncstat = nf90_close(ncid_r)

    ! ntime do gelo
    ncstat = nf90_open(trim(cfg_datocn_ice_file), NF90_NOWRITE, ncid_r)
    allocate(ice0(nx,ny), ice1(nx,ny), iceout(nx,ny))
    iceout = fill_val
    if (ncstat == NF90_NOERR) then
      ncstat = nf90_inq_dimid(ncid_r, 'time', dimid_nt_i)
      if (ncstat /= NF90_NOERR) ncstat = nf90_inq_dimid(ncid_r, 'Time', dimid_nt_i)
      if (ncstat == NF90_NOERR) then
        ncstat = nf90_inquire_dimension(ncid_r, dimid_nt_i, len=ntime_i)
      else
        ntime_i = ntime
      end if
      tidx0_i = mod(int(sec_since_epoch / real(dt_data_i8, ESMF_KIND_R8)), ntime_i) + 1
      tidx1_i = mod(tidx0_i, ntime_i) + 1
      alpha_i  = alpha
      ncstat = nf90_inq_varid(ncid_r, trim(cfg_datocn_ice_varname), varid_src)
      if (ncstat == NF90_NOERR) then
        ncstat = nf90_get_var(ncid_r, varid_src, ice0, start=[1,1,tidx0_i], count=[nx,ny,1])
        ncstat = nf90_get_var(ncid_r, varid_src, ice1, start=[1,1,tidx1_i], count=[nx,ny,1])
        iceout = (1.0_ESMF_KIND_R8 - alpha_i)*ice0 + alpha_i*ice1
        if (cfg_datocn_ice_pct) iceout = iceout / 100.0_ESMF_KIND_R8
        iceout = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, iceout))
        where (abs(ice0) > 1.0e10_ESMF_KIND_R8 .or. abs(ice1) > 1.0e10_ESMF_KIND_R8) &
          iceout = fill_val
      end if
      ncstat = nf90_close(ncid_r)
    end if

    ! ── Correntes superficiais (opcional) ─────────────────────────────────────
    ! CORREÇÃO B-59b: WriteDOCNDiag usava tidx0/tidx1 calculados para o SST
    ! (ntime=365) para ler o OISST_cur.nc (ntime=1) → nf90_get_var falhava
    ! silenciosamente com start=[1,1,88] → f0/f1 com dados SST residuais
    ! (>10 K → fill_val) → So_u=So_v=fill no diag (mas ModelAdvance OK).
    ! Fix: ler ntime do próprio cur_file e calcular tidx0_cur/tidx1_cur.
    if (len_trim(cfg_datocn_cur_file) > 0) then
      ncstat = nf90_open(trim(cfg_datocn_cur_file), NF90_NOWRITE, ncid_r)
      if (ncstat == NF90_NOERR) then
        ! Calcular tidx específico para o cur_file (ntime pode ser 1 ou 365)
        ncstat = nf90_inq_dimid(ncid_r, 'time', dimid_nt_cur)
        if (ncstat /= NF90_NOERR) ncstat = nf90_inq_dimid(ncid_r, 'Time', dimid_nt_cur)
        if (ncstat == NF90_NOERR) then
          ncstat = nf90_inquire_dimension(ncid_r, dimid_nt_cur, len=ntime_cur)
        else
          ntime_cur = 1
        end if
        tidx0_cur = mod(int(sec_since_epoch / real(dt_data_i8, ESMF_KIND_R8)), ntime_cur) + 1
        tidx1_cur = mod(tidx0_cur, ntime_cur) + 1

        ! So_u
        ncstat = nf90_inq_varid(ncid_r, trim(cfg_datocn_cur_u_varname), varid_src)
        if (ncstat == NF90_NOERR) then
          ncstat = nf90_get_var(ncid_r, varid_src, f0, start=[1,1,tidx0_cur], count=[nx,ny,1])
          ncstat = nf90_get_var(ncid_r, varid_src, f1, start=[1,1,tidx1_cur], count=[nx,ny,1])
          uout = (1.0_ESMF_KIND_R8 - alpha)*f0 + alpha*f1
          where (abs(uout) >= 10.0_ESMF_KIND_R8 .or. &
                 abs(f0)   >= 10.0_ESMF_KIND_R8 .or. &
                 abs(f1)   >= 10.0_ESMF_KIND_R8) uout = fill_val
        else
          uout = 0.0_ESMF_KIND_R8
        end if
        ! So_v
        ncstat = nf90_inq_varid(ncid_r, trim(cfg_datocn_cur_v_varname), varid_src)
        if (ncstat == NF90_NOERR) then
          ncstat = nf90_get_var(ncid_r, varid_src, f0, start=[1,1,tidx0_cur], count=[nx,ny,1])
          ncstat = nf90_get_var(ncid_r, varid_src, f1, start=[1,1,tidx1_cur], count=[nx,ny,1])
          vout = (1.0_ESMF_KIND_R8 - alpha)*f0 + alpha*f1
          where (abs(vout) >= 10.0_ESMF_KIND_R8 .or. &
                 abs(f0)   >= 10.0_ESMF_KIND_R8 .or. &
                 abs(f1)   >= 10.0_ESMF_KIND_R8) vout = fill_val
        else
          vout = 0.0_ESMF_KIND_R8
        end if
        ncstat = nf90_close(ncid_r)
      else
        uout = 0.0_ESMF_KIND_R8; vout = 0.0_ESMF_KIND_R8
      end if
    else
      ! cur_file vazio: correntes zero (esperado — sem arquivo de correntes configurado)
      uout = 0.0_ESMF_KIND_R8; vout = 0.0_ESMF_KIND_R8
    end if

    ! ── Escrever NetCDF de diagnóstico (grade DOCN original, não subconjunto) ──
    nlon_diag = nx
    nlat_diag = ny
    dname = trim(cfg_import_diag_dir)
    call execute_command_line('mkdir -p '//trim(dname), wait=.true.)
    fname = trim(dname)//'/docn_import_'//trim(tstamp)//'.nc'

    allocate(lon_ax(nlon_diag), lat_ax(nlat_diag))
    ! lon nativo OISST: 0° → 359.75° (centro de célula do primeiro ponto = 0°).
    ! NÃO usar -180→179.75: o dado lido de nf90_get_var(start=[1,1,tidx]) começa
    ! em lon=0 do OISST — se lon_ax começar em -180 os continentes ficam
    ! deslocados 180° no mapa (bug visual). O postproc_monan_import.py faz o
    ! roll 0→360 para -180→180 no momento da plotagem (sem alterar o arquivo).
    do i = 1, nlon_diag
      lon_ax(i) = (real(i-1, ESMF_KIND_R8)) * (360.0_ESMF_KIND_R8/nlon_diag)
    end do
    do j = 1, nlat_diag
      lat_ax(j) = -90.0_ESMF_KIND_R8 + (j-1) * (180.0_ESMF_KIND_R8/(nlat_diag-1))
    end do

    ncstat = nf90_create(trim(fname), NF90_CLOBBER, ncid_w)
    if (ncstat /= NF90_NOERR) then
      call ESMF_LogWrite('WriteDOCNDiag: nf90_create falhou: '// &
        trim(nf90_strerror(ncstat)), ESMF_LOGMSG_WARNING)
      goto 99
    end if

    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'Conventions',  'CF-1.8')
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'title', &
      'DOCN importState — SST/gelo interpolados por passo (campo global)')
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'institution', &
      'INPE/CGCT/DIMNT — GT Acoplamento MONAN')
    write(tstamp,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
      yy,'-',mm,'-',dd,'T',hh,':',mn,':',ss
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'valid_time',   trim(tstamp))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'datocn_mode',  trim(cfg_datocn_mode))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'sst_source',   trim(cfg_datocn_sst_file))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'ice_source',   trim(cfg_datocn_ice_file))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'sst_varname',  trim(cfg_datocn_sst_varname))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'ice_varname',  trim(cfg_datocn_ice_varname))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'ice_pct',      merge('true ', 'false', cfg_datocn_ice_pct))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'tidx0',        tidx0)
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'tidx1',        tidx1)
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'alpha',        real(alpha,4))
    ncstat = nf90_put_att(ncid_w, NF90_GLOBAL, 'method', &
      'PET0 direct re-read (B-58v2) — grid='//trim(merge('1440x720','360x180 ', &
       trim(cfg_datocn_mode)=='netcdf')))

    ncstat = nf90_def_dim(ncid_w, 'lon', nlon_diag, dimid_lon)
    ncstat = nf90_def_dim(ncid_w, 'lat', nlat_diag, dimid_lat)

    ncstat = nf90_def_var(ncid_w, 'lon',  NF90_DOUBLE, [dimid_lon], varid_lon)
    ncstat = nf90_put_att(ncid_w, varid_lon, 'units', 'degrees_east')
    ncstat = nf90_def_var(ncid_w, 'lat',  NF90_DOUBLE, [dimid_lat], varid_lat)
    ncstat = nf90_put_att(ncid_w, varid_lat, 'units', 'degrees_north')

    ncstat = nf90_def_var(ncid_w, 'So_t',     NF90_DOUBLE, [dimid_lon,dimid_lat], varid_sst)
    ncstat = nf90_put_att(ncid_w, varid_sst, 'long_name', 'SST interpolada (OISST→NUOPC)')
    ncstat = nf90_put_att(ncid_w, varid_sst, 'units', 'K')
    ncstat = nf90_put_att(ncid_w, varid_sst, 'valid_min', 250.0_ESMF_KIND_R8)
    ncstat = nf90_put_att(ncid_w, varid_sst, 'valid_max', 315.0_ESMF_KIND_R8)
    ncstat = nf90_put_att(ncid_w, varid_sst, '_FillValue', fill_val)

    ncstat = nf90_def_var(ncid_w, 'Si_ifrac', NF90_DOUBLE, [dimid_lon,dimid_lat], varid_ice)
    ncstat = nf90_put_att(ncid_w, varid_ice, 'long_name', 'Fracao de gelo marinho')
    ncstat = nf90_put_att(ncid_w, varid_ice, 'units', '1')
    ncstat = nf90_put_att(ncid_w, varid_ice, 'valid_min', 0.0_ESMF_KIND_R8)
    ncstat = nf90_put_att(ncid_w, varid_ice, 'valid_max', 1.0_ESMF_KIND_R8)
    ncstat = nf90_put_att(ncid_w, varid_ice, '_FillValue', fill_val)

    ncstat = nf90_def_var(ncid_w, 'So_u', NF90_DOUBLE, [dimid_lon,dimid_lat], varid_u)
    ncstat = nf90_put_att(ncid_w, varid_u, 'long_name', 'Corrente zonal (zero se cur_file vazio)')
    ncstat = nf90_put_att(ncid_w, varid_u, 'units', 'm/s')
    ncstat = nf90_put_att(ncid_w, varid_u, '_FillValue', fill_val)

    ncstat = nf90_def_var(ncid_w, 'So_v', NF90_DOUBLE, [dimid_lon,dimid_lat], varid_v)
    ncstat = nf90_put_att(ncid_w, varid_v, 'long_name', 'Corrente meridional')
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
    call ESMF_LogWrite('WriteDOCNDiag: '//trim(fname)//' (tidx0='// &
      trim(adjustl(int2str(tidx0)))//' alpha='// &
      trim(adjustl(real2str(real(alpha,4))))//') [B-58v2]', ESMF_LOGMSG_INFO)

    99 continue
    if (allocated(f0))    deallocate(f0, f1, fout)
    if (allocated(ice0))  deallocate(ice0, ice1, iceout)
    if (allocated(uout))  deallocate(uout, vout)
    if (allocated(lon_ax)) deallocate(lon_ax, lat_ax)

  contains
    function int2str(n) result(s)
      integer, intent(in) :: n
      character(len=12) :: s
      write(s,'(I12)') n
    end function int2str
    function real2str(x) result(s)
      real, intent(in) :: x
      character(len=12) :: s
      write(s,'(F8.4)') x
    end function real2str

  end subroutine WriteDOCNDiag

end module DOCN_cap_mod
