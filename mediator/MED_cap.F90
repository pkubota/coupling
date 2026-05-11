!==============================================================================!
! MED_cap.F90 - NUOPC Mediator com fallback MPAS -> DATM                       !
!==============================================================================!
!                                                                              !
! O mediador recebe campos atmosfericos de duas fontes:                        !
!   - MPAS (primario): campos com sufixo _mpas                                 !
!   - DATM (fallback): campos sem sufixo                                       !
!                                                                              !
! Logica de fallback:                                                          !
!   1. Tenta obter campos do MPAS (Sa_u10m_mpas, etc.)                         !
!   2. Se disponivel, usa MPAS                                                 !
!   3. Se nao, usa DATM (Sa_u10m, etc.)                                        !
!                                                                              !
! Recebe SST do OCN e calcula bulk NCAR.                                       !
! Exporta fluxos Foxx_* para o OCN.                                            !
!                                                                              !
! CORRECOES APLICADAS:                                                         !
!   1. So_t (SST) realizado na grade OCN (era ATM - bug critico)               !
!   2. InitializeDataComplete usa NUOPC_MediatorGet em vez de GridCompGet      !
!   3. RegridOrCopy agora tem ramo else que copia direto se rh nao criado      !
!   4. Busca por Sa_u10m_mpas (MPAS) em vez de Sa_u10m (DATM) no IDC           !
!==============================================================================!
module MED_cap_mod
  use ESMF
  use ESMF, only: ESMF_State, ESMF_StateGet
  use NUOPC, only: NUOPC_CompDerive, NUOPC_CompSpecialize, NUOPC_CompSetEntryPoint
  use NUOPC, only: NUOPC_CompFilterPhaseMap, NUOPC_Advertise, NUOPC_Realize
  use NUOPC, only: NUOPC_SetTimestamp, NUOPC_CompAttributeSet
  use NUOPC, only: NUOPC_CompAttributeGet, NUOPC_CompAttributeAdd
  use NUOPC_Mediator, only: med_routine_SS          => SetServices
  use NUOPC_Mediator, only: med_label_DataInitialize => label_DataInitialize
  use NUOPC_Mediator, only: med_label_Advance        => label_Advance
  use NUOPC_Mediator, only: med_label_CheckImport    => label_CheckImport
  use NUOPC_Mediator, only: NUOPC_MediatorGet
  use netcdf
  implicit none
  private
  public :: SetServices

  !----------------------------------------------------------------------------
  ! Constantes fisicas (Large & Yeager 2009)
  !----------------------------------------------------------------------------
  real(ESMF_KIND_R8), parameter :: rho_air    = 1.225_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: Cd_neut    = 1.3e-3_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: Ch_neut    = 1.0e-3_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: Ce_neut    = 1.15e-3_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: Cp_air     = 1004.67_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: L_evap     = 2.501e6_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: T_freeze   = 273.15_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: eps_q      = 0.622_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: es_coef_a  = 611.2_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: es_coef_b  = 17.67_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: es_coef_c  = 243.5_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: sigma_sb   = 5.67e-8_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: albedo_ocn = 0.06_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: SST_default = 290.0_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: f_vis_dir  = 0.285_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: f_vis_dif  = 0.285_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: f_nir_dir  = 0.215_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: f_nir_dif  = 0.215_ESMF_KIND_R8

  ! Limiar de fracao de terra para decisao terra/oceano no MONAN.
  ! Ponto e oceano se lndfrac < LFRAC_OCEAN_THRESHOLD.
  ! Valor 0.5 = conservador (maioria terra -> tratar como terra).
  real(ESMF_KIND_R8), parameter :: LFRAC_OCEAN_THRESHOLD = 0.5_ESMF_KIND_R8

  ! Limiar de binarizacao da mascara apos regrid NEAREST_STOD.
  real(ESMF_KIND_R8), parameter :: MASK_BIN_THRESHOLD    = 0.5_ESMF_KIND_R8

  !----------------------------------------------------------------------------
  ! Estado interno do mediador
  !----------------------------------------------------------------------------
  type :: MED_InternalState
    ! Grade ATM (regular 640x320 para calculo do bulk)
    type(ESMF_Grid) :: atm_grid

    ! Grade OCN (para campos exportados ao oceano)
    type(ESMF_Grid) :: ocn_grid

    ! Campos internos na grade ATM
    type(ESMF_Field) :: f_taux_atm, f_tauy_atm, f_sen_atm, f_evap_atm
    type(ESMF_Field) :: f_lwnet_atm, f_swvdr_atm, f_swvdf_atm
    type(ESMF_Field) :: f_swidr_atm, f_swidf_atm
    type(ESMF_Field) :: f_rain_atm, f_snow_atm, f_pslv_atm
    type(ESMF_Field) :: f_ifrac_atm, f_duu10n_atm, f_sst_atm
    type(ESMF_Field) :: f_sst_atm_tmp 
    ! RouteHandles
    type(ESMF_RouteHandle) :: rh_atm2ocn   ! ATM -> OCN
    type(ESMF_RouteHandle) :: rh_ocn2atm_nearest
    type(ESMF_RouteHandle) :: rh_ocn2atm   ! OCN -> ATM

    ! -----------------------------------------------------------------------
    ! NOVO: Campos de mascara oceano-continente
    ! f_ocn_mask_ocn : mascara 0/1 lida do MOM6 (ocean_static.nc, var "wet")
    !                  na grade OCN nativa do MOM6.
    ! f_ocn_mask_atm : mascara regridada OCN->ATM via NEAREST_STOD e
    !                  binarizada; usada para zerar fluxos bulk em pontos terra.
    ! mask_loaded    : flag - true apos LoadOceanMask ter sido executada com
    !                  sucesso.
    ! -----------------------------------------------------------------------
    type(ESMF_Field) :: f_ocn_mask_ocn
    type(ESMF_Field) :: f_ocn_mask_atm
    logical          :: mask_loaded = .false.


    ! Mascara oceano/continente
    !PK real(ESMF_KIND_R8), allocatable :: ocn_mask_atm(:,:)
    
    logical :: rh_created    = .false.
    logical :: use_mpas_atm  = .false.  ! controlado por atributo NUOPC "use_mpas_atm"

    ! Diagnosticos NetCDF: SST na grade ATM e na malha Voronoi MPAS
    ! diag_step  : contador de passos (incrementado a cada MediatorAdvance)
    ! diag_dir   : diretorio de saida (atributo NUOPC "med_diag_dir", default "./diag")
    ! diag_freq  : frequencia de escrita em passos (atributo "med_diag_freq", default 1)
    ! diag_enabled: .true. se diagnostico habilitado (atributo "med_diag_sst", default "false")
    integer          :: diag_step    = 0
    integer          :: diag_freq    = 1
    logical          :: diag_enabled = .false.
    character(len=256) :: diag_dir   = './diag'

  end type MED_InternalState

  type :: MED_InternalStateWrapper
    type(MED_InternalState), pointer :: wrap => null()
  end type MED_InternalStateWrapper

  ! Campos de import do MPAS (primario) - com sufixo _mpas
  integer, parameter :: n_import_mpas = 10   ! +1 para Sa_lfrac_mpas (fracao terra MONAN)
  character(len=32), parameter :: import_mpas_names(n_import_mpas) = [ &
    "Sa_u10m_mpas  ", "Sa_v10m_mpas  ", "Sa_tbot_mpas  ", "Sa_shum_mpas  ", "Sa_pslv_mpas  ", &
    "Faxa_swdn_mpas", "Faxa_lwdn_mpas", "Faxa_rain_mpas", "Faxa_snow_mpas", &
    "Sa_lfrac_mpas " ]   ! NOVO: fracao de terra do MONAN para reconciliacao de mascara costeira

  ! Campos de import do DATM (fallback) - sem sufixo
  integer, parameter :: n_import_datm = 9
  character(len=32), parameter :: import_datm_names(n_import_datm) = [ &
    "Sa_u10m   ", "Sa_v10m   ", "Sa_tbot   ", "Sa_shum   ", "Sa_pslv   ", &
    "Faxa_swdn ", "Faxa_lwdn ", "Faxa_rain ", "Faxa_snow "]

  ! Campos de export para o OCN
  integer, parameter :: n_export = 14
  character(len=32), parameter :: export_names(n_export) = [ &
    "Foxx_taux     ", "Foxx_tauy     ", "Foxx_sen      ", "Foxx_evap     ", "Foxx_lwnet    ", &
    "Foxx_swnet_vdr", "Foxx_swnet_vdf", "Foxx_swnet_idr", "Foxx_swnet_idf", &
    "Faxa_rain     ", "Faxa_snow     ", "Sa_pslv       ", "Si_ifrac      ", "So_duu10n     " ]

contains

  !============================================================================
  ! SetServices
  !============================================================================
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call NUOPC_CompDerive(gcomp, med_routine_SS, rc=rc)
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

    call NUOPC_CompSpecialize(gcomp, specLabel=med_label_DataInitialize, &
      specRoutine=InitializeDataComplete, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=med_label_Advance, &
      specRoutine=MediatorAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=med_label_CheckImport, &
      specRoutine=CheckImportNoop, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

  !============================================================================
  ! CheckImportNoop
  !============================================================================
  subroutine CheckImportNoop(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    rc = ESMF_SUCCESS
    call ESMF_LogWrite('MED: CheckImport desabilitado (no-op)', ESMF_LOGMSG_INFO)
  end subroutine CheckImportNoop

  !============================================================================
  ! InitializeP0
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
  ! InitializeAdvertise
  !============================================================================
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    integer :: n
    logical                         :: isPresent, isSet
    character(len=8)                :: attr_val
    ! use_mpas_atm lido aqui apenas para log; o valor persistente fica no estado interno
    ! criado em InitializeRealize.
    logical          :: use_mpas_atm_local
    !PK ! --- ADICIONAR ao final de InitializeAdvertise, antes do LogWrite final ---
    !PK 
    !PK type(MED_InternalStateWrapper) :: iswrap
    !PK type(MED_InternalState), pointer :: is
    integer          :: localrc

    rc = ESMF_SUCCESS

    !PK allocate(iswrap%wrap)
    !PK is => iswrap%wrap

    ! CORRECAO 3: InternalState NAO e alocado aqui.
    ! A versao anterior alocava iswrap%wrap em InitializeAdvertise E depois
    ! em InitializeRealize, causando double-allocation e memory leak.
    ! O InternalState e criado uma unica vez, em InitializeRealize.

    use_mpas_atm_local = .false.
    !PK ! Inicializar todos os campos logicos do InternalState
    !PK is%use_mpas_atm = use_mpas_atm_advertise
    !PK is%rh_created   = .false.
    !PK 
    !PK call ESMF_GridCompSetInternalState(gcomp, iswrap, rc)
    !PK if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !PK   line=__LINE__, file=__FILE__)) return    !--- Le atributo use_mpas_atm definido pelo driver em esm.F90 ---
    ! Valores aceitos: "true" ou "false" (default: "false" = usa DATM)
    call NUOPC_CompAttributeGet(gcomp, name="use_mpas_atm", &
      value=attr_val, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) then
      !PK use_mpas_atm_advertise = (trim(attr_val) == "true")
      use_mpas_atm_local = (trim(attr_val) == "true")
    end if
    !PK if (use_mpas_atm_advertise) then
    if (use_mpas_atm_local) then
      call ESMF_LogWrite('MED: use_mpas_atm=true (MPAS/MONAN como fonte primaria)', &
        ESMF_LOGMSG_INFO)
    else
      call ESMF_LogWrite('MED: use_mpas_atm=false (DATM como fonte)', &
        ESMF_LOGMSG_INFO)
    end if

    !PK ! Anuncia campos de import conforme a fonte atmosferica configurada.
    !PK ! CRITICO: o NUOPC aborta em IPDv03p6 se um campo anunciado nao tiver
    !PK! conector ativo. Por isso MPAS e DATM sao anunciados exclusivamente.
    !PK if (use_mpas_atm_advertise) then
    ! CRITICO: anunciar MPAS ou DATM exclusivamente.
    ! NUOPC aborta em IPDv03p6 se um campo anunciado nao tiver conector ativo.
    ! -----------------------------------------------------------------------
    ! Anunciar campos de import:
    ! Modo MPAS: 9 campos atmosfericos + Sa_lfrac_mpas (fracao terra MONAN)
    ! Modo DATM: 9 campos atmosfericos sem sufixo
    !
    ! Sa_lfrac_mpas e opcional ? o mediador tolera ausencia dele em Advance,
    ! mas o anuncio aqui e necessario para o NUOPC criar o conector se o
    ! MONAN_cap exportar o campo.
    ! -----------------------------------------------------------------------
    if (use_mpas_atm_local) then
      ! Modo MPAS: anuncia campos _mpas (fornecidos pelo MPAS_cap)
      ! Sa_lfrac_mpas e opcional: nao abortar se o MONAN nao o oferecer
      do n = 1, n_import_mpas
        call NUOPC_Advertise(importState, StandardName=trim(import_mpas_names(n)), &
          TransferOfferGeomObject="cannot provide", &
          SharePolicyField="share", rc=localrc)
        ! Sa_lfrac_mpas e opcional: nao abortar se o MONAN nao o oferecer
        if (trim(import_mpas_names(n)) == "Sa_lfrac_mpas") then
          if (localrc /= ESMF_SUCCESS) then
            call ESMF_LogWrite( &
              'MED: AVISO - Sa_lfrac_mpas nao anunciado (MONAN nao exporta lndfrac)', &
              ESMF_LOGMSG_WARNING)
          end if
        else
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) then
            rc = localrc; return
          end if
        end if
      end do
    else
      ! Modo DATM: anuncia campos sem sufixo (fornecidos pelo DATM_cap)
      ! SharePolicyField="share" evita bondLevel ambiguo para Faxa_rain/snow
      ! que aparecem tanto no importState quanto no exportState do MED.
      do n = 1, n_import_datm
        call NUOPC_Advertise(importState, StandardName=trim(import_datm_names(n)), &
          TransferOfferGeomObject="cannot provide", &
          SharePolicyField="share", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      end do
    end if

    ! Advertise So_t (SST do OCN) - sempre presente (conector OCN->MED ativo nos dois modos)
    !           So_t (SST) e So_omask (mascara terra/oceano) do OCN - sempre presentes
    call NUOPC_Advertise(importState, StandardName="So_t", &
      TransferOfferGeomObject="cannot provide", &
      SharePolicyField="share", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! So_omask: mascara binaria MOM6 (1=oceano, 0=terra)
    ! Usada em InitializeRealize para LoadOceanMask via So_omask OCN->MED
    call NUOPC_Advertise(importState, StandardName="So_omask", &
      TransferOfferGeomObject="cannot provide", &
      SharePolicyField="share", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Advertise campos de export para o OCN
    do n = 1, n_export
      call NUOPC_Advertise(exportState, StandardName=trim(export_names(n)), &
        TransferOfferGeomObject="will provide", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do

    call ESMF_LogWrite('MED: InitializeAdvertise concluido', ESMF_LOGMSG_INFO)
  end subroutine InitializeAdvertise

  !============================================================================
  ! InitializeRealize
  ! CORRECAO 1: So_t (SST) realizado na grade OCN, nao na ATM.
  !   O campo So_t vem do componente OCN (grade ocn_grid). Realiza-lo na
  !   atm_grid fazia com que o routehandle OCN->ATM tivesse src e dst na
  !   mesma grade, tornando o regrid incorreto.
  !============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_Grid)  :: atm_grid, ocn_grid
    type(ESMF_Field) :: tmp_field
    type(MED_InternalStateWrapper) :: iswrap
    type(MED_InternalState), pointer :: is
    integer :: nx_atm, ny_atm, nx_ocn, ny_ocn, i, j, n
    integer :: nxp_ocn, nyp_ocn
    real(ESMF_KIND_R8), pointer :: coordX(:,:), coordY(:,:)
    integer :: ncid, varid, dimid
    ! BUG FIX: aloca (nyp, nxp) para espelhar o layout NetCDF do supergrid.
    ! NetCDF: x(nyp,nxp) dim0=lat(nyp) dim1=lon(nxp).
    ! nf90_get_var (col-major Fortran): ocn_lon(i,j) = x_nc[i-1, j-1]
    ! => ocn_lon(i,j) = longitude do ponto (lat_sg=i-1, lon_sg=j-1).
    ! Acesso ponto T do OCN (ig,jg): ocn_lon(2*jg-1, 2*ig-1)
    !   (jg = indice lat, ig = indice lon do OCN, ambos base-1 globais)
    real(ESMF_KIND_R8), allocatable :: ocn_lon(:,:), ocn_lat(:,:),wet_g(:,:)
    real(ESMF_KIND_R8) :: dlon
    integer, pointer :: mask_ptr(:,:) => null()
    logical             :: isPresent, isSet,file_ok
    character(len=256)  :: attr_val
    character(len=256)  :: ocn_hgrid_file, ocn_mask_file
    character(len=256)  :: line_tag
    integer :: clbnd(2), cubnd(2)   ! compute lower/upper bound do tile MPI
    integer             :: localrc
    integer :: ig, jg               ! indices globais OCN (base-1, lon e lat)
      
    rc = ESMF_SUCCESS

    ! CORRECAO 3 (cont.): InternalState e alocado UMA UNICA VEZ aqui.
    ! InitializeAdvertise nao aloca mais o IS; alocamos diretamente.
    ! Recuperar estado interno
    allocate(iswrap%wrap)
    is => iswrap%wrap
    is%rh_created   = .false.
    is%use_mpas_atm = .false.
    is%mask_loaded   = .false.

    ! Ler use_mpas_atm ANTES de realizar os campos.
    !PK call ESMF_GridCompGetInternalState(gcomp, iswrap, rc)
    !PK if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !PK line=__LINE__, file=__FILE__)) return
    !PK is => iswrap%wrap
    call NUOPC_CompAttributeGet(gcomp, name="use_mpas_atm", &
      value=attr_val, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) is%use_mpas_atm = (trim(attr_val) == "true")

    ! ---------------------------------------------------------------------------
    ! Dimensoes da grade ATM intermediaria do mediador (para calculo do bulk NCAR)
    ! Lidas do atributo NUOPC "med_nx_atm"/"med_ny_atm" definido pelo driver.
    ! Fallback: 640x320 (resolucao ~0.5 grau, adequada para qualquer config MPAS).
    ! ---------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------
    ! Dimensoes da grade ATM intermediaria (regular 640x320 para calculo do bulk)
    nx_atm = 640
    ny_atm = 320
    call NUOPC_CompAttributeGet(gcomp, name="med_nx_atm", &
      value=attr_val, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) read(attr_val, *) nx_atm

    call NUOPC_CompAttributeGet(gcomp, name="med_ny_atm", &
      value=attr_val, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) read(attr_val, *) ny_atm

    write(*,'(A,2I6)') 'MED: grade ATM intermediaria (bulk) nx_atm x ny_atm = ', nx_atm, ny_atm
    call ESMF_LogWrite('MED: grade ATM intermediaria configurada', ESMF_LOGMSG_INFO)

    ! Arquivo do supergrid OCN
    ! ---------------------------------------------------------------------------
    ! Dimensoes da grade OCN: NIGLOBAL x NJGLOBAL lidos de ocean_hgrid.nc
    ! O arquivo ocean_hgrid.nc (supergrid FRE-NCtools) tem dimensoes:
    !   nxp = NIGLOBAL + 1  (pontos de borda do supergrid)
    !   nyp = NJGLOBAL + 1 
    ! Portanto: nx_ocn = nxp - 1,  ny_ocn = nyp - 1
    ! O caminho do arquivo e lido do atributo NUOPC "ocn_hgrid_file" (obrigatorio).
    ! ---------------------------------------------------------------------------
    ocn_hgrid_file = 'INPUT/ocean_hgrid.nc'  ! default
    call NUOPC_CompAttributeGet(gcomp, name="ocn_hgrid_file", &
      value=ocn_hgrid_file, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (.not. (isPresent .and. isSet)) ocn_hgrid_file = 'INPUT/ocean_hgrid.nc'

    ! NOVO: Arquivo de mascara OCN (ocean_static.nc, variavel "wet")
    ocn_mask_file = 'INPUT/ocean_static.nc'
    call NUOPC_CompAttributeGet(gcomp, name="ocn_mask_file", &
      value=ocn_mask_file, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (.not. (isPresent .and. isSet)) ocn_mask_file = 'INPUT/ocean_static.nc'

    ! Configurar diagnosticos NetCDF de SST
    is%diag_enabled = .false.
    call NUOPC_CompAttributeGet(gcomp, name="med_diag_sst", &
      value=attr_val, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) is%diag_enabled = (trim(attr_val) == "true")

    is%diag_dir  = './diag'
    call NUOPC_CompAttributeGet(gcomp, name="med_diag_dir", &
      value=is%diag_dir, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (.not. (isPresent .and. isSet)) is%diag_dir = './diag'

    is%diag_freq = 1
    call NUOPC_CompAttributeGet(gcomp, name="med_diag_freq", &
      value=attr_val, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (isPresent .and. isSet) read(attr_val, *) is%diag_freq

    if (is%diag_enabled) then
      write(line_tag,'(A,A,A,I0,A)') 'MED: diagnostico SST habilitado -> ', &
        trim(is%diag_dir), ' (freq=', is%diag_freq, ' passos)'
      call ESMF_LogWrite(trim(line_tag), ESMF_LOGMSG_INFO)
      call execute_command_line('mkdir -p '//trim(is%diag_dir), wait=.true.)
    else
      call ESMF_LogWrite('MED: diagnostico SST DESABILITADO '// &
        '(definir med_diag_sst=true no driver para habilitar)', ESMF_LOGMSG_INFO)
    end if

    call ESMF_LogWrite('MED: lendo dimensoes OCN de '//trim(ocn_hgrid_file), ESMF_LOGMSG_INFO)

    ! Ler dimensoes do supergrid
    rc = nf90_open(trim(ocn_hgrid_file), NF90_NOWRITE, ncid)
    if (rc /= NF90_NOERR) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="MED: falha ao abrir "//trim(ocn_hgrid_file)//": "//trim(nf90_strerror(rc)), &
        line=__LINE__, file=__FILE__, rcToReturn=rc); return
    end if
    rc = nf90_inq_dimid(ncid, "nxp", dimid)
    if (rc /= NF90_NOERR) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="MED: dimensao 'nxp' nao encontrada em "//trim(ocn_hgrid_file), &
        line=__LINE__, file=__FILE__, rcToReturn=rc); return
    end if
    rc = nf90_inquire_dimension(ncid, dimid, len=nxp_ocn)
    rc = nf90_inq_dimid(ncid, "nyp", dimid)
    if (rc /= NF90_NOERR) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="MED: dimensao 'nyp' nao encontrada em "//trim(ocn_hgrid_file), &
        line=__LINE__, file=__FILE__, rcToReturn=rc); return
    end if
    rc = nf90_inquire_dimension(ncid, dimid, len=nyp_ocn)
    rc = nf90_close(ncid)
    rc = ESMF_SUCCESS

    !--------------------------------------------------------------------------
    ! Criar grade ATM intermediaria (para calculo do bulk NCAR)
    ! Dimensoes lidas do atributo med_nx_atm/med_ny_atm (default 640x320)
    !--------------------------------------------------------------------------
    atm_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), maxIndex=(/nx_atm, ny_atm/), &
      indexflag=ESMF_INDEX_GLOBAL, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha ao criar grade ATM", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_GridAddCoord(atm_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Preencher coordenadas da grade ATM
    call ESMF_GridGetCoord(atm_grid, coordDim=1, staggerloc=ESMF_STAGGERLOC_CENTER, &
      farrayPtr=coordX, &
      computationalLBound=clbnd, computationalUBound=cubnd, rc=rc)
    do jg = clbnd(2), cubnd(2)
      do ig = clbnd(1), cubnd(1)
         coordX(ig, jg) = (ig-1)*(360.0_ESMF_KIND_R8/nx_atm) + &
                      (360.0_ESMF_KIND_R8/nx_atm)*0.5_ESMF_KIND_R8
      end do
    end do

    call ESMF_GridGetCoord(atm_grid, coordDim=2, staggerloc=ESMF_STAGGERLOC_CENTER, &
      farrayPtr=coordY, &
      computationalLBound=clbnd, computationalUBound=cubnd, rc=rc)
    do jg = clbnd(2), cubnd(2)
      do ig = clbnd(1), cubnd(1)
        coordY(ig, jg) = -90.0_ESMF_KIND_R8 + (jg-1)*(180.0_ESMF_KIND_R8/ny_atm) + &
                         (180.0_ESMF_KIND_R8/ny_atm)*0.5_ESMF_KIND_R8
      end do
    end do
    !--------------------------------------------------------------------------
    ! Criar grade OCN (nx_ocn x ny_ocn lidos de ocean_hgrid.nc)
    ! Coordenadas reais lidas do supergrid (pontos T = indices pares do supergrid)
    !--------------------------------------------------------------------------
    ! example 
    ! Criar grade OCN tripolar (180x158) - simplificada
    ! Criar grade OCN (180x158)
    !--------------------------------------------------------------------------
    ! Pontos T (tracer) = supergrid / 2  (FRE-NCtools: nxp = NIGLOBAL+1)
    nx_ocn = (nxp_ocn - 1)/2
    ny_ocn = (nyp_ocn - 1)/2

    write(*,'(A,2I6)') 'MED: grade OCN lida de ocean_hgrid.nc: nx_ocn x ny_ocn = ', nx_ocn, ny_ocn
    call ESMF_LogWrite('MED: dimensoes OCN lidas com sucesso', ESMF_LOGMSG_INFO)
    !PK ocn_grid = ESMF_GridCreateNoPeriDim(minIndex=(/1,1/), maxIndex=(/nx_ocn, ny_ocn/), &
    !PK   indexflag=ESMF_INDEX_GLOBAL, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)

    ocn_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                      maxIndex=(/nx_ocn, ny_ocn/), &
                                      indexflag=ESMF_INDEX_GLOBAL, &
                                      coordSys=ESMF_COORDSYS_SPH_DEG, &
                                      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha ao criar grade OCN", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_GridAddCoord(ocn_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Preencher coordenadas da grade OCN (simplificada - lat/lon regulares)
    !--------------------------------------------------------------------------
    ! Ler coordenadas reais do supergrid ocean_hgrid.nc
    ! O supergrid tem dimensoes (nyp, nxp) = (2*ny+1, 2*nx+1).
    ! Pontos T (center) = indices pares do supergrid (base 0):
    !   lon_T(i,j) = x_sg(2*j, 2*i)   i=0..nx-1, j=0..ny-1
    !   lat_T(i,j) = y_sg(2*j, 2*i)
    ! Em Fortran (base 1): x_sg(2*i-1, 2*j-1)
    ! BUG FIX (transposicao nxp/nyp):
    !   NetCDF layout de x,y: (nyp, nxp)  -> dim0=lat, dim1=lon
    !   nf90_get_var preenche col-major Fortran:
    !     ocn_lon(i,j) = x_nc[i-1, j-1]  (i percorre nyp=lat, j percorre nxp=lon)
    !   Portanto alocamos (nyp_ocn, nxp_ocn) para que:
    !     ocn_lon(lat_sg, lon_sg) = longitude no supergrid
    !   Acesso ao ponto T(ig, jg) [ig=col_lon base-1, jg=row_lat base-1]:
    !     lon_T = ocn_lon(2*jg-1, 2*ig-1)   <- jg indexa dim1=lat, ig indexa dim2=lon
    !     lat_T = ocn_lat(2*jg-1, 2*ig-1)
    !
    ! BUG FIX (indices globais MPI):
    !   farrayPtr retorna a fatia LOCAL do processo MPI.
    !   Em runs paralelas, lbound(coordX,1) pode ser > 1 (indice global).
    !   Usamos ESMF_GridGetCoord com computationalLBound/UBound para obter
    !   os indices globais ig, jg e mapear corretamente no supergrid global.
    ! -----------------------------------------------------------------------
    allocate(ocn_lon(nxp_ocn, nyp_ocn))   ! (nxp, nyp) em ordem Fortran col-major
    allocate(ocn_lat(nxp_ocn, nyp_ocn))

    rc = nf90_open(trim(ocn_hgrid_file), NF90_NOWRITE, ncid)
    if (rc /= NF90_NOERR) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="MED: falha ao reabrir "//trim(ocn_hgrid_file)//" para coordenadas", &
        line=__LINE__, file=__FILE__, rcToReturn=rc); return
    end if
    ! Ler variavel x (longitude) - shape no NetCDF: (nyp, nxp) -> em Fortran: (nxp, nyp)
    rc = nf90_inq_varid(ncid, "x", varid)
    if (rc /= NF90_NOERR) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="MED: variavel 'x' nao encontrada em "//trim(ocn_hgrid_file), &
        line=__LINE__, file=__FILE__, rcToReturn=rc); return
    end if
    rc = nf90_get_var(ncid, varid, ocn_lon)
    ! Ler variavel y (latitude) - shape no NetCDF: (nyp, nxp) -> em Fortran: (nxp, nyp)
    rc = nf90_inq_varid(ncid, "y", varid)
    if (rc /= NF90_NOERR) then
      call ESMF_LogSetError(ESMF_FAILURE, &
        msg="MED: variavel 'y' nao encontrada em "//trim(ocn_hgrid_file), &
        line=__LINE__, file=__FILE__, rcToReturn=rc); return
    end if
    rc = nf90_get_var(ncid, varid, ocn_lat)
    rc = nf90_close(ncid)
    rc = ESMF_SUCCESS

    ! Preencher coordenadas da grade OCN (simplificada - lat/lon regulares)
    ! Preencher coordenadas da grade OCN a partir do supergrid (pontos T)
    ! Supergrid base-1 Fortran: ponto T(i,j) -> supergrid(2*i-1, 2*j-1)
    call ESMF_GridGetCoord(ocn_grid, coordDim=1, staggerloc=ESMF_STAGGERLOC_CENTER, &
          farrayPtr=coordX, &
          computationalLBound=clbnd, computationalUBound=cubnd, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) return
    do j = clbnd(2), cubnd(2)    ! j = jg global (lat)
       do i = clbnd(1), cubnd(1)  ! i = ig global (lon)
          !PK coordX(i,j) = (i-1) * (360.0_ESMF_KIND_R8/nx_ocn)
          !PK coordX(i,j) = ocn_lon(2*i-1, 2*j-1)
          ! offset local: coordX e 1-based localmente; global -> local: i-clbnd(1)+1
          coordX(i, j) = ocn_lon(2*i-1, 2*j-1)
       end do
    end do
            !PK 
            !PK       do j = 1, cubnd(2)-clbnd(2)
            !PK         jg = clbnd(2) + j - 1          ! indice global j no OCN grid (1-based)
            !PK         do i = 1, cubnd(1)-clbnd(1)
            !PK           ig = clbnd(1) + i - 1        ! indice global i no OCN grid (1-based)
            !PK           ! Ponto T do supergrid (base-1): supergrid(2*ig-1, 2*jg-1)
            !PK          coordX(i,j) = ocn_lon(2*ig-1, 2*jg-1)
            !PK         end do
            !PK       end do
      
    call ESMF_GridGetCoord(ocn_grid, coordDim=2, staggerloc=ESMF_STAGGERLOC_CENTER, &
       farrayPtr=coordY, &
       computationalLBound=clbnd, computationalUBound=cubnd, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, file=__FILE__)) return
    do j = clbnd(2), cubnd(2)
       do i = clbnd(1), cubnd(1)
          !coordY(i,j) = -90.0_ESMF_KIND_R8 + (j-1)*(180.0_ESMF_KIND_R8/ny_ocn) + &
          !              (180.0_ESMF_KIND_R8/ny_ocn)/2.0_ESMF_KIND_R8
          !coordY(i,j) = ocn_lat(2*i-1, 2*j-1)
          coordY(i, j) = ocn_lat(2*i-1, 2*j-1)
       end do
    end do
            !PK      do j = 1, cubnd(2)-clbnd(2)
            !PK        jg = clbnd(2) + j - 1
            !PK        do i = 1, cubnd(1)-clbnd(1)
            !PK          ig = clbnd(1) + i - 1
            !PK          coordY(i,j) = ocn_lat(2*ig-1, 2*jg-1)
            !PK        end do
            !PK      end do
    ! Log de verificacao dos cantos do dominio OCN

    write(line_tag,'(A,4F10.4)') &
      'MED: OCN T-pt cantos lon: (1,1)(nx,1)(1,ny)(nx,ny) = ', &
      ocn_lon(1,1), ocn_lon(2*nx_ocn-1,1), ocn_lon(1,2*ny_ocn-1), ocn_lon(2*nx_ocn-1,2*ny_ocn-1)
    call ESMF_LogWrite(trim(line_tag), ESMF_LOGMSG_INFO)
    write(line_tag,'(A,4F10.4)') &
      'MED: OCN T-pt cantos lat: (1,1)(nx,1)(1,ny)(nx,ny) = ', &
      ocn_lat(1,1), ocn_lat(2*nx_ocn-1,1), ocn_lat(1,2*ny_ocn-1), ocn_lat(2*nx_ocn-1,2*ny_ocn-1)
    call ESMF_LogWrite(trim(line_tag), ESMF_LOGMSG_INFO)

   ! Apos as linhas de ESMF_GridGetCoord/fill para ocn_grid (apos linha ~663)
   ! e ANTES de deallocate(ocn_lon, ocn_lat):
    call ESMF_GridAddItem(ocn_grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha GridAddItem MASK OCN", &
        line=__LINE__, file=__FILE__)) return
    call ESMF_GridGetItem(ocn_grid, itemflag=ESMF_GRIDITEM_MASK, &
        staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=mask_ptr, &
        computationalLBound=clbnd, computationalUBound=cubnd, rc=rc)
    ! Ler wet do ocean_static.nc (mesmo arquivo da LoadOceanMask)
    ! wet_g(ny, nx) ? layout NetCDF: dim0=lat, dim1=lon
    allocate(wet_g(nx_ocn, ny_ocn))
    wet_g = 1.0_ESMF_KIND_R8
    localrc = nf90_open(trim(ocn_mask_file), NF90_NOWRITE, ncid)
!PK    if (localrc == NF90_NOERR) then
!PK      localrc = nf90_inq_varid(ncid, "wet", varid)
!PK      if (localrc == NF90_NOERR) localrc = nf90_get_var(ncid, varid, wet_g)
!PK      localrc = nf90_close(ncid)
!PK    end if
!
    if (localrc == NF90_NOERR) then
      ! Tentar "wet" (MOM6/ocean_static.nc) depois "mask" (FRE ocean_mask.nc)
      localrc = nf90_inq_varid(ncid, "wet", varid)
!PK      if (localrc /= NF90_NOERR) &
!PK        localrc = nf90_inq_varid(ncid, "mask", varid)

      if (localrc == NF90_NOERR) then
        localrc = nf90_get_var(ncid, varid, wet_g)
        if (localrc == NF90_NOERR) then
          file_ok = .true.
        else
          call ESMF_LogWrite( &
            'MED LoadOceanMask: AVISO - falha lendo wet/mask; usando wet=1', &
            ESMF_LOGMSG_WARNING)
          wet_g = 1.0_ESMF_KIND_R8
        end if
      else
        call ESMF_LogWrite( &
          'MED LoadOceanMask: AVISO - var wet/mask nao encontrada; usando wet=1', &
          ESMF_LOGMSG_WARNING)
      end if
      do jg = clbnd(2), cubnd(2)
         do ig = clbnd(1), cubnd(1)
            ! wet=0 -> terra  -> mask=0 (ESMF ignora onde mask=srcMaskValues)
            ! wet=1 -> oceano -> mask=1
          if (wet_g(ig, jg) > 0.5_ESMF_KIND_R8) then
             mask_ptr(ig, jg) = 1
          else
             mask_ptr(ig, jg) = 0
          end if
        end do
      end do


      localrc = nf90_close(ncid)
    else
      call ESMF_LogWrite( &
        'MED LoadOceanMask: AVISO - '//trim(ocn_mask_file)//' nao aberto; usando wet=1', &
        ESMF_LOGMSG_WARNING)
    end if


    deallocate(wet_g)
    deallocate(ocn_lon, ocn_lat)

    !--------------------------------------------------------------------------
    ! Realizar campos de import conforme a fonte atmosferica configurada.
    ! Espelha exatamente o que foi anunciado em InitializeAdvertise.
    !--------------------------------------------------------------------------
    if (is%use_mpas_atm) then
      do n = 1, n_import_mpas
        tmp_field = ESMF_FieldCreate(grid=atm_grid, typekind=ESMF_TYPEKIND_R8, &
          staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(import_mpas_names(n)), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call NUOPC_Realize(importState, field=tmp_field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      end do
    else
      do n = 1, n_import_datm
        tmp_field = ESMF_FieldCreate(grid=atm_grid, typekind=ESMF_TYPEKIND_R8, &
          staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(import_datm_names(n)), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
        call NUOPC_Realize(importState, field=tmp_field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=__FILE__)) return
      end do
    end if

    ! So_t (SST) na grade OCN
    !--------------------------------------------------------------------------
    ! Realizar So_t (SST) na grade ATM (placeholder)
    ! CORRECAO 1: So_t (SST) realizado na grade OCN (era atm_grid - bug critico)
    ! O campo So_t vem do oceano, portanto sua grade nativa e ocn_grid.
    ! Realiza-lo na atm_grid causava conflito ao criar o routehandle OCN->ATM.
    !--------------------------------------------------------------------------
    tmp_field = ESMF_FieldCreate(grid=ocn_grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, name="So_t", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Realize(importState, field=tmp_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! So_omask: mascara terra/oceano do MOM6, recebida via conector OCN->MED
    ! Grade nativa e ocn_grid; realizada aqui para o conector poder conectar
    tmp_field = ESMF_FieldCreate(grid=ocn_grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, name="So_omask", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha criar campo So_omask", &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Realize(importState, field=tmp_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha realizar So_omask", &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! Realizar campos de export na grade OCN
    !--------------------------------------------------------------------------
    do n = 1, n_export
      tmp_field = ESMF_FieldCreate(grid=ocn_grid, typekind=ESMF_TYPEKIND_R8, &
        staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(export_names(n)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Realize(exportState, field=tmp_field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do

    !--------------------------------------------------------------------------
    ! Inicializar estado interno (IS ja alocado acima - apenas atribuir grades)
    !--------------------------------------------------------------------------
    !PK allocate(iswrap%wrap)
    !PK is => iswrap%wrap
    is%atm_grid      = atm_grid
    is%ocn_grid      = ocn_grid
    !PK is%rh_created    = .false.

    !PK !--- Le atributo use_mpas_atm e armazena no estado interno ---
    !PK is%use_mpas_atm = .false.
    !PK call NUOPC_CompAttributeGet(gcomp, name="use_mpas_atm", &
    !PK   value=attr_val, isPresent=isPresent, isSet=isSet, rc=rc)
    !PK if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !PK   line=__LINE__, file=__FILE__)) return
    !PK if (isPresent .and. isSet) is%use_mpas_atm = (trim(attr_val) == "true")

    ! Criar campos internos na grade ATM
    call CreateInternalField(is%f_taux_atm,   atm_grid, "med_taux",   rc)
    call CreateInternalField(is%f_tauy_atm,   atm_grid, "med_tauy",   rc)
    call CreateInternalField(is%f_sen_atm,    atm_grid, "med_sen",    rc)
    call CreateInternalField(is%f_evap_atm,   atm_grid, "med_evap",   rc)
    call CreateInternalField(is%f_lwnet_atm,  atm_grid, "med_lwnet",  rc)
    call CreateInternalField(is%f_swvdr_atm,  atm_grid, "med_swvdr",  rc)
    call CreateInternalField(is%f_swvdf_atm,  atm_grid, "med_swvdf",  rc)
    call CreateInternalField(is%f_swidr_atm,  atm_grid, "med_swidr",  rc)
    call CreateInternalField(is%f_swidf_atm,  atm_grid, "med_swidf",  rc)
    call CreateInternalField(is%f_rain_atm,   atm_grid, "med_rain",   rc)
    call CreateInternalField(is%f_snow_atm,   atm_grid, "med_snow",   rc)
    call CreateInternalField(is%f_pslv_atm,   atm_grid, "med_pslv",   rc)
    call CreateInternalField(is%f_ifrac_atm,  atm_grid, "med_ifrac",  rc)
    call CreateInternalField(is%f_duu10n_atm, atm_grid, "med_duu10n", rc)
    ! f_sst_atm: campo de SST interpolado para a grade ATM (destino do OCN->ATM)
    call CreateInternalField(is%f_sst_atm,    atm_grid, "med_sst",    rc)
    call CreateInternalField(is%f_sst_atm_tmp,    atm_grid, "med_sst_tmp",    rc)

    ! Zerar campos internos
    call ZeroInternalField(is%f_taux_atm,   rc)
    call ZeroInternalField(is%f_tauy_atm,   rc)
    call ZeroInternalField(is%f_sen_atm,    rc)
    call ZeroInternalField(is%f_evap_atm,   rc)
    call ZeroInternalField(is%f_lwnet_atm,  rc)
    call ZeroInternalField(is%f_swvdr_atm,  rc)
    call ZeroInternalField(is%f_swvdf_atm,  rc)
    call ZeroInternalField(is%f_swidr_atm,  rc)
    call ZeroInternalField(is%f_swidf_atm,  rc)
    call ZeroInternalField(is%f_rain_atm,   rc)
    call ZeroInternalField(is%f_snow_atm,   rc)
    call ZeroInternalField(is%f_pslv_atm,   rc)
    call ZeroInternalField(is%f_ifrac_atm,  rc)
    call ZeroInternalField(is%f_duu10n_atm, rc)
    ! Inicializa SST com valor padrao (nao zero, para evitar bulk erratico no t=0)
    !PK call ZeroInternalField(is%f_sst_atm,    rc)
    call FillInternalField(is%f_sst_atm, SST_default, rc)
    call FillInternalField(is%f_sst_atm_tmp, SST_default, rc)

    ! -----------------------------------------------------------------------
    ! NOVO: Carregar mascara oceano-continente do MOM6.
    ! LoadOceanMask le ocean_static.nc (variavel "wet"), realiza regrid
    ! NEAREST_STOD OCN->ATM, binariza e armazena em is%f_ocn_mask_atm. Nao e fatal se falhar.
    ! Deve ser chamado APOS is%ocn_grid e is%atm_grid estarem definidos.
    ! -----------------------------------------------------------------------
    call LoadOceanMask(gcomp, is, trim(ocn_mask_file), nx_ocn, ny_ocn, rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite( &
        'MED: AVISO - LoadOceanMask falhou; prosseguindo sem mascara explicita', &
        ESMF_LOGMSG_WARNING)
      rc = ESMF_SUCCESS   ! nao abortar - mascara e auxiliar, nao critica
    end if

    call ESMF_GridCompSetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('MED: InitializeRealize concluido', ESMF_LOGMSG_INFO)
  end subroutine InitializeRealize

  !============================================================================
  ! LoadOceanMask
  !
  ! Le a mascara terra/oceano do MOM6 (variavel "wet" em ocean_static.nc ou
  ! "mask" em ocean_mask.nc), cria campos ESMF nas grades OCN e ATM,
  ! interpola via NEAREST_STOD (binariza antes do regrid nao e necessario pois
  ! NEAREST_STOD preserva o valor do vizinho mais proximo, mas binarizamos
  ! depois por seguranca numerica) e armazena o resultado em
  ! is%f_ocn_mask_atm para uso nos loops bulk do MediatorAdvance.
  !
  ! Layout NetCDF de "wet": (ny_ocn, nx_ocn) ? dim0=lat, dim1=lon.
  ! nf90_get_var col-major Fortran: wet_g(i,j) = wet_nc[i-1, j-1]
  !   => wet_g(jg, ig) = mascara no ponto (lat_idx=jg-1, lon_idx=ig-1)
  !   => acesso: mptr_ocn(local_i, local_j) = wet_g(jg, ig)
  !
  ! Por que NEAREST_STOD e nao BILINEAR para a mascara?
  !   A mascara e binaria (0=terra, 1=oceano). BILINEAR produziria valores
  !   intermediarios (ex.: 0.4) que tornariam ambiguo o status do ponto.
  !   NEAREST_STOD propaga o valor 0 ou 1 do ponto mais proximo, mantendo
  !   a mascara binaria apos a interpolacao (antes da binarizacao explicita).
  !============================================================================
 subroutine LoadOceanMask(gcomp, is, mask_file, nx_ocn_arg, ny_ocn_arg, rc)
    type(ESMF_GridComp),     intent(in)    :: gcomp
    type(MED_InternalState), intent(inout) :: is
    character(len=*),        intent(in)    :: mask_file
    integer,                 intent(in)    :: nx_ocn_arg   ! nx da grade OCN (dim lon)
    integer,                 intent(in)    :: ny_ocn_arg   ! ny da grade OCN (dim lat)
    integer,                 intent(out)   :: rc

    integer :: ncid, varid
    integer :: ig, jg, iloc, jloc          ! indices globais e locais
    integer :: clbnd(2), cubnd(2)          ! bounds globais do tile MPI local
    real(ESMF_KIND_R8), pointer     :: mptr_ocn(:,:), mptr_atm(:,:)
    integer, pointer :: mask_ptr(:,:) => null()
    ! wet_g(ny, nx): mascara global lida do NetCDF.
    ! Alocada com (ny_ocn_arg, nx_ocn_arg) para espelhar dim0=lat, dim1=lon.
    real(ESMF_KIND_R8), allocatable :: wet_g(:,:)
    type(ESMF_RouteHandle) :: rh_mask
    integer                :: localrc
    logical                :: file_ok

    rc = ESMF_SUCCESS
    call ESMF_LogWrite('MED: LoadOceanMask: lendo '//trim(mask_file), ESMF_LOGMSG_INFO)

    ! ------------------------------------------------------------------
    ! 1. Criar campos de mascara nas grades OCN e ATM
    ! ------------------------------------------------------------------
    is%f_ocn_mask_ocn = ESMF_FieldCreate(grid=is%ocn_grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, name="med_ocn_mask_ocn", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="MED LoadOceanMask: falha criando f_ocn_mask_ocn", &
      line=__LINE__, file=__FILE__)) return

    is%f_ocn_mask_atm = ESMF_FieldCreate(grid=is%atm_grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, name="med_ocn_mask_atm", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="MED LoadOceanMask: falha criando f_ocn_mask_atm", &
      line=__LINE__, file=__FILE__)) return

    ! ------------------------------------------------------------------
    ! 2. Ler mascara do NetCDF para wet_g global (ny_ocn x nx_ocn)
    !    Fallback seguro: se qualquer passo falhar, wet_g fica com 1.0
    !    (tudo oceano) e o processamento continua sem abortar.
    ! ------------------------------------------------------------------
    allocate(wet_g(nx_ocn_arg, ny_ocn_arg))
    wet_g = 1.0_ESMF_KIND_R8   ! default: tudo oceano
    file_ok = .false.

    localrc = nf90_open(trim(mask_file), NF90_NOWRITE, ncid)
    if (localrc == NF90_NOERR) then
      ! Tentar "wet" (MOM6/ocean_static.nc) depois "mask" (FRE ocean_mask.nc)
      localrc = nf90_inq_varid(ncid, "wet", varid)
      if (localrc /= NF90_NOERR) &
        localrc = nf90_inq_varid(ncid, "mask", varid)

      if (localrc == NF90_NOERR) then
        localrc = nf90_get_var(ncid, varid, wet_g)
        if (localrc == NF90_NOERR) then
          file_ok = .true.
        else
          call ESMF_LogWrite( &
            'MED LoadOceanMask: AVISO - falha lendo wet/mask; usando wet=1', &
            ESMF_LOGMSG_WARNING)
          wet_g = 1.0_ESMF_KIND_R8
        end if
      else
        call ESMF_LogWrite( &
          'MED LoadOceanMask: AVISO - var wet/mask nao encontrada; usando wet=1', &
          ESMF_LOGMSG_WARNING)
      end if
      localrc = nf90_close(ncid)
    else
      call ESMF_LogWrite( &
        'MED LoadOceanMask: AVISO - '//trim(mask_file)//' nao aberto; usando wet=1', &
        ESMF_LOGMSG_WARNING)
    end if

    if (file_ok) then
      call ESMF_LogWrite('MED LoadOceanMask: mascara lida com sucesso', ESMF_LOGMSG_INFO)
    end if

    ! ------------------------------------------------------------------
    ! 3. Preencher campo OCN com wet_g (MPI-safe via indices globais).
    !
    ! farrayPtr e a fatia LOCAL do processo. computationalLBound/UBound
    ! retornam os indices GLOBAIS do tile local (base-1).
    ! Mapeamento:
    !   ig = indice global de lon (coluna)  [clbnd(1)..cubnd(1)]
    !   jg = indice global de lat (linha)   [clbnd(2)..cubnd(2)]
    !   wet_g(jg, ig) = mascara nesse ponto (layout NetCDF: dim0=lat, dim1=lon)
    !   iloc = ig - clbnd(1) + 1   (indice local no farrayPtr, base-1)
    !   jloc = jg - clbnd(2) + 1
    ! ------------------------------------------------------------------
    call ESMF_FieldGet(is%f_ocn_mask_ocn, farrayPtr=mptr_ocn, &
      computationalLBound=clbnd, computationalUBound=cubnd, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) then
      deallocate(wet_g); return
    end if

    do jg = clbnd(2), cubnd(2)
      jloc = jg - clbnd(2) + 1
      do ig = clbnd(1), cubnd(1)
        iloc = ig - clbnd(1) + 1
        mptr_ocn(ig, jg) = wet_g(ig, jg)   
        !mptr_ocn(iloc, jloc) = wet_g(ig, jg)     ! era mptr_ocn(ig,jg)=wet_g(ig,jg)
      end do
    end do
    deallocate(wet_g)   ! liberado aqui em TODOS os caminhos de sucesso

    ! ------------------------------------------------------------------
    ! 4. Regrid mascara OCN -> ATM via NEAREST_STOD
    !    Pontos ATM fora do dominio OCN ficam com valor inicial (0=terra).
    ! ------------------------------------------------------------------
    call ESMF_FieldGet(is%f_ocn_mask_atm, farrayPtr=mptr_atm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    mptr_atm = 0.0_ESMF_KIND_R8   ! inicializar: pontos nao mapeados = terra

    call ESMF_FieldRegridStore( &
      srcField       = is%f_ocn_mask_ocn,              &
      dstField       = is%f_ocn_mask_atm,              &
      routehandle    = rh_mask,                        &
      regridmethod   = ESMF_REGRIDMETHOD_NEAREST_STOD, &
      unmappedaction = ESMF_UNMAPPEDACTION_IGNORE,     &
      rc             = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="MED LoadOceanMask: falha FieldRegridStore NEAREST_STOD OCN->ATM", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_FieldRegrid(is%f_ocn_mask_ocn, is%f_ocn_mask_atm, rh_mask, &
      zeroregion=ESMF_REGION_TOTAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="MED LoadOceanMask: falha FieldRegrid mascara OCN->ATM", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_RouteHandleDestroy(rh_mask, nogarbage=.true., rc=rc)
    rc = ESMF_SUCCESS

    ! ------------------------------------------------------------------
    ! 5. Binarizar: > MASK_BIN_THRESHOLD -> 1.0 (oceano); resto -> 0.0
    !    Seguranca numerica: NEAREST_STOD deve preservar 0/1 mas pode
    !    haver imprecisao de ponto flutuante no limite.
    ! ------------------------------------------------------------------
    call ESMF_FieldGet(is%f_ocn_mask_atm, farrayPtr=mptr_atm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    where (mptr_atm > MASK_BIN_THRESHOLD)
      mptr_atm = 1.0_ESMF_KIND_R8
    elsewhere
      mptr_atm = 0.0_ESMF_KIND_R8
    end where

    is%mask_loaded = .true.
    call ESMF_LogWrite( &
      'MED: mascara OCN carregada e interpolada para grade ATM', &
      ESMF_LOGMSG_INFO)

  end subroutine LoadOceanMask

  !============================================================================
  ! ReconcileCoastalMask
  !
  ! Reconcilia a mascara MOM6 de terra/oceano do MED com a fracao de terra
  ! exportada pelo MONAN (Sa_lfrac_mpas).
  ! No-op se mask_loaded=.false. ou Sa_lfrac_mpas ausente.
  !
  ! Estrategia conservadora (intersecao):
  !   Um ponto e tratado como OCEANO somente se:
  !     (a) MOM6 diz oceano  (is%f_ocn_mask_atm = 1)  E
  !     (b) MONAN diz oceano (lndfrac < LFRAC_OCEAN_THRESHOLD)
  !
  ! Isso evita que o MED calcule fluxos bulk usando dados atmosfericos de
  ! pontos continentais do MONAN e os envie para pontos oceanicos do MOM6,
  ! o que causaria erros de balanco de energia na linha de costa.
  !
  ! O campo mascara e modificado IN PLACE a cada chamada. A mascara base
  ! (lida do MOM6 em LoadOceanMask) e preservada em f_ocn_mask_ocn e pode
  ! ser re-regridada se necessario ? mas como a linha de costa do MONAN nao
  ! muda durante a simulacao, a reconciliacao converge apos o primeiro Advance.
  !
  ! NOTA: Se Sa_lfrac_mpas nao estiver disponivel (modo DATM ou MONAN nao
  ! exportando lndfrac), a rotina retorna sem modificar a mascara e emite
  ! um aviso. Nao e fatal.
  !============================================================================
  subroutine ReconcileCoastalMask(is, importState, rc)
    type(MED_InternalState), intent(inout) :: is
    type(ESMF_State),        intent(inout) :: importState
    integer,                 intent(out)   :: rc

    type(ESMF_Field)            :: f_lndfrac
    real(ESMF_KIND_R8), pointer :: lndfrac(:,:), mask_atm(:,:)
    integer :: i, j
    integer :: localrc
    integer :: n_land_masked, n_total

    rc = ESMF_SUCCESS

    if (.not. is%mask_loaded) return

    ! Verificar existencia de Sa_lfrac_mpas via itemNameList (evita ERROR no log ESMF)
    block
      integer :: item_count, n_item
      character(len=64), allocatable :: item_names(:)
      logical :: found_it
      found_it = .false.
      call ESMF_StateGet(importState, itemCount=item_count, rc=localrc)
      if (localrc == ESMF_SUCCESS .and. item_count > 0) then
        allocate(item_names(item_count))
        call ESMF_StateGet(importState, itemNameList=item_names, rc=localrc)
        if (localrc == ESMF_SUCCESS) then
          do n_item = 1, item_count
            if (trim(item_names(n_item)) == "Sa_lfrac_mpas") then
              found_it = .true.; exit
            end if
          end do
        end if
        deallocate(item_names)
      end if
      if (.not. found_it) then
        call ESMF_LogWrite( &
          'MED ReconcileCoastalMask: Sa_lfrac_mpas ausente; mascara MOM6 mantida', &
          ESMF_LOGMSG_INFO)
        return
      end if
    end block

    ! Tentar Sa_lfrac_mpas (MONAN) ? opcional
    call ESMF_StateGet(importState, itemName="Sa_lfrac_mpas", &
      field=f_lndfrac, rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      call ESMF_LogWrite( &
        'MED ReconcileCoastalMask: Sa_lfrac_mpas nao disponivel; mascara MOM6 mantida sem alteracao', &
        ESMF_LOGMSG_WARNING)
      return
    end if

    call ESMF_FieldGet(f_lndfrac, farrayPtr=lndfrac, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_FieldGet(is%f_ocn_mask_atm, farrayPtr=mask_atm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Intersecao conservadora: zera mascara onde MONAN ve terra
    n_land_masked = 0
    n_total = (ubound(mask_atm,1)-lbound(mask_atm,1)+1) * &
              (ubound(mask_atm,2)-lbound(mask_atm,2)+1)

    do j = lbound(mask_atm,2), ubound(mask_atm,2)
      do i = lbound(mask_atm,1), ubound(mask_atm,1)
        ! Ponto e terra no MONAN mas oceano no MOM6: tratar como terra
        if (lndfrac(i,j) >= LFRAC_OCEAN_THRESHOLD .and. &
            mask_atm(i,j) > MASK_BIN_THRESHOLD) then
          mask_atm(i,j) = 0.0_ESMF_KIND_R8
          n_land_masked = n_land_masked + 1
        end if
      end do
    end do

    if (n_land_masked > 0) then
      write(*,'(A,I8,A,I8,A)') &
        'MED ReconcileCoastalMask: ', n_land_masked, &
        ' de ', n_total, ' pontos ATM reclassificados como terra pela mascara MONAN'
      call ESMF_LogWrite('MED: reconciliacao costeira MOM6/MONAN aplicada', ESMF_LOGMSG_INFO)
    end if

  end subroutine ReconcileCoastalMask

  !============================================================================
  ! CreateInternalField
  !============================================================================
  subroutine CreateInternalField(field, grid, name, rc)
    type(ESMF_Field), intent(out) :: field
    type(ESMF_Grid),  intent(in)  :: grid
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: rc

    field = ESMF_FieldCreate(grid=grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(name), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="MED CreateInternalField: "//trim(name), &
      line=__LINE__, file=__FILE__)) return
  end subroutine CreateInternalField

  !============================================================================
  ! ZeroInternalField
  !============================================================================
  subroutine ZeroInternalField(field, rc)
    type(ESMF_Field), intent(inout) :: field
    integer,          intent(out)   :: rc

    real(ESMF_KIND_R8), pointer :: fptr(:,:)
    rc = ESMF_SUCCESS
    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    fptr = 0.0_ESMF_KIND_R8

  end subroutine ZeroInternalField

  !============================================================================
  ! FillInternalField - preenche campo com valor constante
  !============================================================================
  subroutine FillInternalField(field, value, rc)
    type(ESMF_Field),   intent(inout) :: field
    real(ESMF_KIND_R8), intent(in)    :: value
    integer,            intent(out)   :: rc

    real(ESMF_KIND_R8), pointer :: fptr(:,:)
    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    fptr = value

  end subroutine FillInternalField

  !============================================================================
  ! InitializeDataComplete - cria routehandles
  ! CORRECAO 2: usa NUOPC_MediatorGet em vez de ESMF_GridCompGet para obter
  !   importState/exportState, que e a API correta para mediadores NUOPC.
  ! CORRECAO 4: busca Sa_u10m_mpas (MPAS, grade ATM) para obter a grade ATM,
  !   em vez de Sa_u10m (DATM), que pode nao estar presente se o DATM nao
  !   tiver sido conectado ainda. Usa fallback para Sa_u10m caso necessario.
  !============================================================================
  subroutine InitializeDataComplete(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    type(ESMF_State)         :: importState, exportState
    type(ESMF_Clock)         :: clock
    type(ESMF_Field)         :: atm_field, ocn_field, exp_field,ocn_mask_field
    type(ESMF_Grid)          :: atm_grid!PK
    type(MED_InternalStateWrapper) :: iswrap
    type(MED_InternalState), pointer :: is
    integer :: fieldCount, i, localrc
    character(len=64), allocatable :: fieldNameList(:)
    real(ESMF_KIND_R8), pointer :: fptr(:,:)
    integer, target  :: maskValues(1)

    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    is => iswrap%wrap

!PK    call ESMF_GridCompGet(gcomp, importState=importState, exportState=exportState, rc=rc)
!PK    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!PK      line=__LINE__, file=__FILE__)) return

    ! CORRECAO 2: NUOPC_MediatorGet e a API correta para mediadores
    call NUOPC_MediatorGet(gcomp, mediatorClock=clock, &
      importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Obtem campo de referencia para a grade ATM conforme o modo ativo.
    ! use_mpas_atm ja esta no estado interno (lido em InitializeRealize).
    if (is%use_mpas_atm) then
      call ESMF_StateGet(importState, itemName="Sa_u10m_mpas", &
        field=atm_field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="MED IDC: Sa_u10m_mpas nao encontrado", &
        line=__LINE__, file=__FILE__)) return
    else
      call ESMF_StateGet(importState, itemName="Sa_u10m", &
        field=atm_field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="MED IDC: Sa_u10m nao encontrado", &
        line=__LINE__, file=__FILE__)) return
    end if

    ! Obter campo de export para o OCN (Foxx_taux esta na grade OCN)
    call ESMF_StateGet(exportState, itemName="Foxx_taux", field=exp_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha Foxx_taux", &
      line=__LINE__, file=__FILE__)) return

    ! Routehandle ATM -> OCN (BILINEAR para fluxos continuos)
    ! Criar routehandle ATM -> OCN
    ! CORRECAO 5: BILINEAR e o metodo correto para interpolacao de fluxos
    ! continuos (momentum, calor sensivel, latente, SW, LW).
    ! NEAREST_STOD era usado anteriormente mas nao conserva energia e
    ! introduz discontinuidades na costa.
    call ESMF_FieldRegridStore( &
      srcField       = is%f_taux_atm,   &
      dstField       = exp_field,       &
      routehandle    = is%rh_atm2ocn,   &
      regridmethod   = ESMF_REGRIDMETHOD_BILINEAR, &
      unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
      rc             = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="MED: falha FieldRegridStore ATM->OCN", &
      line=__LINE__, file=__FILE__)) return

    ! Criar routehandle OCN -> ATM
    ! Routehandle OCN -> ATM (So_t ja esta na grade OCN)
    ! Criar routehandle OCN -> ATM
    ! So_t esta agora corretamente na grade OCN (ver InitializeRealize)

    call ESMF_StateGet(importState, itemName="So_t", field=ocn_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha So_t", &
      line=__LINE__, file=__FILE__)) return
    ! limpar campo destino 
    call ESMF_FieldGet(is%f_sst_atm_tmp, farrayPtr=fptr, rc=rc)
    fptr = 0.0_ESMF_KIND_R8
    ! limpar campo destino 
    call ESMF_FieldGet(is%f_sst_atm, farrayPtr=fptr, rc=rc)
    fptr = 0.0_ESMF_KIND_R8

    ! interpolar
    maskValues(1) = 0   ! valor que indica "ignorar este ponto"
    call ESMF_StateGet(importState, itemName="So_omask", field=ocn_mask_field, rc=rc)

    ! Routehandle principal: bilinear para domínio regular
    call ESMF_FieldRegridStore( &
      srcField       = ocn_field,       &
      dstField       = is%f_sst_atm,    &
      routehandle    = is%rh_ocn2atm,   &
      regridmethod   = ESMF_REGRIDMETHOD_BILINEAR, & !ESMF_REGRIDMETHOD_PATCH  ! mais suave que BILINEAR nas bordas
      srcMaskValues  = maskValues,         &  ! <- ignora pontos onde omask=0
      polemethod     = ESMF_POLEMETHOD_NONE,      &  ! <= media no polo
      unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
      rc             = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="MED: falha FieldRegridStore OCN->ATM", &
      line=__LINE__, file=__FILE__)) return
    ! Routehandle auxiliar: nearest para preencher pontos não mapeados (seam/polo)
    call ESMF_FieldRegridStore( &
      srcField       = ocn_field,                          &
      dstField       = is%f_sst_atm,                       &
      routehandle    = is%rh_ocn2atm_nearest,              &
      regridmethod   = ESMF_REGRIDMETHOD_NEAREST_STOD,     &
      srcMaskValues  = maskValues,                         &
      unmappedaction = ESMF_UNMAPPEDACTION_IGNORE,         &
      rc             = rc)
    is%rh_created = .true.

    ! Inicializar exportState com valores fisicamente razoaveis
    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    if (fieldCount > 0) then
      allocate(fieldNameList(fieldCount))
      call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
      do i = 1, fieldCount
        call ESMF_StateGet(exportState, itemName=trim(fieldNameList(i)), &
          field=exp_field, rc=rc)
        call ESMF_FieldGet(exp_field, farrayPtr=fptr, rc=rc)
        select case(trim(fieldNameList(i)))
          case('Sa_pslv')
            fptr = 101325.0_ESMF_KIND_R8
          case default
            fptr = 0.0_ESMF_KIND_R8
        end select
      end do
      deallocate(fieldNameList)
    end if

    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataProgress", value="true", rc=rc)
    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)

    call ESMF_LogWrite('MED: InitializeDataComplete SATISFIED', ESMF_LOGMSG_INFO)
  end subroutine InitializeDataComplete

  !============================================================================
  ! MediatorAdvance - com fallback MPAS -> DATM
  !
  ! Fluxo de execucao:
  !   1. Obter campos ATM (MPAS primario / DATM fallback)
  !   2. Regrid SST OCN -> ATM
  !   3. [NOVO] ReconcileCoastalMask: atualizar mascara com lndfrac MONAN
  !   4. Calcular bulk NCAR com guarda de mascara em cada loop
  !   5. Regrid fluxos ATM -> OCN e exportar
  !   6. Atualizar timestamps
  !============================================================================
  subroutine MediatorAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    type(ESMF_State)         :: importState, exportState
    type(ESMF_Clock)         :: clock
    type(ESMF_Time)          :: currTime, nextTime
    type(ESMF_TimeInterval)  :: dt
    type(ESMF_Field)         :: field
    type(MED_InternalStateWrapper) :: iswrap
    type(MED_InternalState), pointer :: is

    ! Campos da ATM do MPAS (primario)
    real(ESMF_KIND_R8), pointer :: uas_mpas(:,:)  => null()
    real(ESMF_KIND_R8), pointer :: vas_mpas(:,:)  => null()
    real(ESMF_KIND_R8), pointer :: tas_mpas(:,:)  => null()
    real(ESMF_KIND_R8), pointer :: shum_mpas(:,:) => null()
    real(ESMF_KIND_R8), pointer :: psl_mpas(:,:)  => null()
    real(ESMF_KIND_R8), pointer :: swdn_mpas(:,:) => null()
    real(ESMF_KIND_R8), pointer :: lwdn_mpas(:,:) => null()
    real(ESMF_KIND_R8), pointer :: rain_mpas(:,:) => null()
    real(ESMF_KIND_R8), pointer :: snow_mpas(:,:) => null()
    logical :: mpas_available

    ! Campos da ATM do DATM (fallback)
    real(ESMF_KIND_R8), pointer :: uas_datm(:,:)  => null()
    real(ESMF_KIND_R8), pointer :: vas_datm(:,:)  => null()
    real(ESMF_KIND_R8), pointer :: tas_datm(:,:)  => null()
    real(ESMF_KIND_R8), pointer :: shum_datm(:,:) => null()
    real(ESMF_KIND_R8), pointer :: psl_datm(:,:)  => null()
    real(ESMF_KIND_R8), pointer :: swdn_datm(:,:) => null()
    real(ESMF_KIND_R8), pointer :: lwdn_datm(:,:) => null()
    real(ESMF_KIND_R8), pointer :: rain_datm(:,:) => null()
    real(ESMF_KIND_R8), pointer :: snow_datm(:,:) => null()

    ! Campos finais (alias para MPAS ou DATM)
    ! Alias para o conjunto ativo
    real(ESMF_KIND_R8), pointer :: uas(:,:), vas(:,:), tas(:,:), shum(:,:)
    real(ESMF_KIND_R8), pointer :: psl(:,:), swdn(:,:), lwdn(:,:)
    real(ESMF_KIND_R8), pointer :: rain(:,:), snow(:,:)

    ! Campos calculados
    real(ESMF_KIND_R8), pointer :: sst(:,:), fptr(:,:),sst_bil(:,:),sst_nst(:,:)

    ! NOVO: Ponteiro para mascara oceano-continente interpolada para grade ATM
    real(ESMF_KIND_R8), pointer :: mask_atm(:,:) => null()

    real(ESMF_KIND_R8) :: wspd, qsat, sst_eff
    integer :: i, j, i1, i2, j1, j2
    integer :: fieldCount, k
    character(len=64), allocatable :: fieldNameList(:)
    character(len=256) :: msg

    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    is => iswrap%wrap

    call NUOPC_MediatorGet(gcomp, mediatorClock=clock, &
      importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=dt, rc=rc)
    nextTime = currTime + dt

    !==========================================================================
    ! 1. TENTAR OBTER CAMPOS DA ATM DO MPAS (PRIMARIO)
    !==========================================================================
    ! use_mpas_atm vem do atributo NUOPC definido em esm.F90.
    ! Se false, pula a tentativa e vai direto ao DATM.
    mpas_available = is%use_mpas_atm

    if (mpas_available) then
      call GetFieldPtrOptional(importState, "Sa_u10m_mpas", uas_mpas, rc)
      if (rc /= ESMF_SUCCESS) mpas_available = .false.
    end if

    if (mpas_available) then
      call GetFieldPtrOptional(importState, "Sa_v10m_mpas",   vas_mpas,  rc)
      call GetFieldPtrOptional(importState, "Sa_tbot_mpas",   tas_mpas,  rc)
      call GetFieldPtrOptional(importState, "Sa_shum_mpas",   shum_mpas, rc)
      call GetFieldPtrOptional(importState, "Sa_pslv_mpas",   psl_mpas,  rc)
      call GetFieldPtrOptional(importState, "Faxa_swdn_mpas", swdn_mpas, rc)
      call GetFieldPtrOptional(importState, "Faxa_lwdn_mpas", lwdn_mpas, rc)
      call GetFieldPtrOptional(importState, "Faxa_rain_mpas", rain_mpas, rc)
      call GetFieldPtrOptional(importState, "Faxa_snow_mpas", snow_mpas, rc)

      if (.not. (associated(uas_mpas)  .and. associated(vas_mpas)  .and. &
                 associated(tas_mpas)  .and. associated(shum_mpas) .and. &
                 associated(psl_mpas)  .and. associated(swdn_mpas) .and. &
                 associated(lwdn_mpas) .and. associated(rain_mpas) .and. &
                 associated(snow_mpas))) then
        mpas_available = .false.
      end if
    end if
    !==========================================================================
    ! 2. SE MPAS NAO DISPONIVEL, USAR DATM (FALLBACK)
    !==========================================================================
    if (.not. mpas_available) then
      call GetFieldPtr(importState, "Sa_u10m",   uas_datm,  rc); if (rc/=ESMF_SUCCESS) return
      call GetFieldPtr(importState, "Sa_v10m",   vas_datm,  rc); if (rc/=ESMF_SUCCESS) return
      call GetFieldPtr(importState, "Sa_tbot",   tas_datm,  rc); if (rc/=ESMF_SUCCESS) return
      call GetFieldPtr(importState, "Sa_shum",   shum_datm, rc); if (rc/=ESMF_SUCCESS) return
      call GetFieldPtr(importState, "Sa_pslv",   psl_datm,  rc); if (rc/=ESMF_SUCCESS) return
      call GetFieldPtr(importState, "Faxa_swdn", swdn_datm, rc); if (rc/=ESMF_SUCCESS) return
      call GetFieldPtr(importState, "Faxa_lwdn", lwdn_datm, rc); if (rc/=ESMF_SUCCESS) return
      call GetFieldPtr(importState, "Faxa_rain", rain_datm, rc); if (rc/=ESMF_SUCCESS) return
      call GetFieldPtr(importState, "Faxa_snow", snow_datm, rc); if (rc/=ESMF_SUCCESS) return

      uas  => uas_datm;  vas  => vas_datm;  tas  => tas_datm
      shum => shum_datm; psl  => psl_datm;  swdn => swdn_datm
      lwdn => lwdn_datm; rain => rain_datm; snow => snow_datm

      call ESMF_LogWrite('MED: Usando DATM (JRA55) como fonte atmosferica (fallback)', &
        ESMF_LOGMSG_INFO)
    else
      uas  => uas_mpas;  vas  => vas_mpas;  tas  => tas_mpas
      shum => shum_mpas; psl  => psl_mpas;  swdn => swdn_mpas
      lwdn => lwdn_mpas; rain => rain_mpas; snow => snow_mpas

      call ESMF_LogWrite('MED: Usando MPAS/MONAN como fonte atmosferica primaria', &
        ESMF_LOGMSG_INFO)
    end if

    i1 = lbound(uas,1); i2 = ubound(uas,1)
    j1 = lbound(uas,2); j2 = ubound(uas,2)

    !==========================================================================
    ! 3. SST: regrid OCN -> ATM (So_t esta agora na grade OCN)
    !==========================================================================
    if (is%rh_created) then
      call ESMF_StateGet(importState, itemName="So_t", field=field, rc=rc)
      ! limpar campo destino
      call ESMF_FieldGet(is%f_sst_atm_tmp, farrayPtr=fptr, rc=rc)
      fptr = 0.0_ESMF_KIND_R8
      ! limpar campo destino
      call ESMF_FieldGet(is%f_sst_atm, farrayPtr=fptr, rc=rc)
      fptr = 0.0_ESMF_KIND_R8
      ! interpolar
      ! 1. Bilinear para o dominio principal
      call ESMF_FieldRegrid(field      , &
                           is%f_sst_atm, &
                           is%rh_ocn2atm, &
                           zeroregion=ESMF_REGION_EMPTY, rc=rc)
      ! 2. Nearest apenas para pontos que ficaram zero (seam line e polo)
      ! Usar campo temporário para nao sobrescrever o bilinear onde ele funcionou
      call ESMF_FieldRegrid(field, &
                            is%f_sst_atm_tmp, &
                            is%rh_ocn2atm_nearest, &
                            zeroregion=ESMF_REGION_EMPTY, rc=rc)
      ! 3. Merge: onde bilinear deu zero, usar nearest
      call ESMF_FieldGet(is%f_sst_atm    , farrayPtr=sst_bil, rc=rc)
      call ESMF_FieldGet(is%f_sst_atm_tmp, farrayPtr=sst_nst, rc=rc)
      where (sst_bil == 0.0_ESMF_KIND_R8) sst_bil = sst_nst
      
      sst=>sst_bil
    else
      ! Routehandles nao criados: usa SST padrao (ja preenchido em InitializeRealize)
      call ESMF_FieldGet(is%f_sst_atm, farrayPtr=sst, rc=rc)
      !PK call GetFieldPtr(importState, "So_t", sst, rc)
    end if

    !==========================================================================
    ! 3b. NOVO: Reconciliar mascara costeira MONAN x MOM6
    !    Chamado a cada Advance pois lndfrac pode variar com o modelo ATM.
    !    A rotina e no-op se is%mask_loaded=.false. ou Sa_lfrac_mpas ausente.
    !==========================================================================
    call ReconcileCoastalMask(is, importState, rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite( &
        'MED: AVISO - ReconcileCoastalMask falhou; continuando sem reconciliacao costeira', &
        ESMF_LOGMSG_WARNING)
      rc = ESMF_SUCCESS
    end if

    ! Obter ponteiro para a mascara (null se nao disponivel)
    nullify(mask_atm)
    if (is%mask_loaded) then
      call ESMF_FieldGet(is%f_ocn_mask_atm, farrayPtr=mask_atm, rc=rc)
      if (rc /= ESMF_SUCCESS) nullify(mask_atm)
    end if

    !PK !==========================================================================
    !PK ! 3c. Diagnostico NetCDF: SST na grade ATM e na malha Voronoi MPAS
    !PK !     Escrito a cada diag_freq passos se med_diag_sst="true"
    !PK !==========================================================================
    !PK is%diag_step = is%diag_step + 1
    !PK if (is%diag_enabled .and. mod(is%diag_step, is%diag_freq) == 0) then
    !PK   call WriteDiagSST(gcomp, is, importState, currTime, is%diag_step, rc)
    !PK   if (rc /= ESMF_SUCCESS) then
    !PK     call ESMF_LogWrite('MED: AVISO - WriteDiagSST falhou no passo '// &
    !PK       trim(adjustl(transfer(is%diag_step, ' '))), ESMF_LOGMSG_WARNING)
    !PK     rc = ESMF_SUCCESS
    !PK   end if
    !PK end if    
    !==========================================================================
    ! 4. CALCULAR BULK NCAR
    !
    ! NOVO: Cada loop inclui guarda de mascara:
    !   if (associated(mask_atm)) then
    !     if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
    !       fptr(i,j) = 0.0; cycle
    !     end if
    !   end if
    !
    ! Isso garante que fluxos calculados sobre pontos continentais sejam
    ! explicitamente zerados antes do regrid ATM->OCN, evitando contaminacao
    ! de pontos costeiros do MOM6 com dados de superficie terrestre do MONAN.
    !==========================================================================

    ! --- Taux ---
    call ESMF_FieldGet(is%f_taux_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2
      do i = i1, i2
        if (associated(mask_atm)) then
          if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
            fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
          end if
        end if
        wspd = sqrt(uas(i,j)**2 + vas(i,j)**2) + 1.0e-10_ESMF_KIND_R8
        fptr(i,j) = rho_air * Cd_neut * wspd * uas(i,j)
      end do
    end do

    ! --- Tauy ---
    call ESMF_FieldGet(is%f_tauy_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2
      do i = i1, i2
        if (associated(mask_atm)) then
          if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
            fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
          end if
        end if
        wspd = sqrt(uas(i,j)**2 + vas(i,j)**2) + 1.0e-10_ESMF_KIND_R8
        fptr(i,j) = rho_air * Cd_neut * wspd * vas(i,j)
      end do
    end do

    ! --- Calor sensivel ---
    call ESMF_FieldGet(is%f_sen_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2
      do i = i1, i2
        if (associated(mask_atm)) then
          if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
            fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
          end if
        end if
        wspd = sqrt(uas(i,j)**2 + vas(i,j)**2) + 1.0e-10_ESMF_KIND_R8
        sst_eff = merge(sst(i,j), SST_default, &
          sst(i,j) > 271.0_ESMF_KIND_R8 .and. sst(i,j) < 308.0_ESMF_KIND_R8)
        fptr(i,j) = rho_air * Cp_air * Ch_neut * wspd * (tas(i,j) - sst_eff)
      end do
    end do

    ! --- Evaporacao ---
    call ESMF_FieldGet(is%f_evap_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2
      do i = i1, i2
        if (associated(mask_atm)) then
          if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
            fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
          end if
        end if
        wspd = sqrt(uas(i,j)**2 + vas(i,j)**2) + 1.0e-10_ESMF_KIND_R8
        sst_eff = merge(sst(i,j), SST_default, &
          sst(i,j) > 271.0_ESMF_KIND_R8 .and. sst(i,j) < 308.0_ESMF_KIND_R8)
        qsat = eps_q * es_coef_a * &
          exp(es_coef_b*(sst_eff-T_freeze)/(sst_eff-T_freeze+es_coef_c)) / &
          max(psl(i,j), 1.0_ESMF_KIND_R8)
        fptr(i,j) = rho_air * Ce_neut * wspd * (shum(i,j) - qsat)
      end do
    end do

    ! --- Radiacao de onda longa liquida ---
    call ESMF_FieldGet(is%f_lwnet_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2
      do i = i1, i2
        if (associated(mask_atm)) then
          if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
            fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
          end if
        end if
        sst_eff = merge(sst(i,j), SST_default, &
          sst(i,j) > 271.0_ESMF_KIND_R8 .and. sst(i,j) < 308.0_ESMF_KIND_R8)
        fptr(i,j) = lwdn(i,j) - sigma_sb * sst_eff**4
      end do
    end do

    ! --- Bandas de onda curta (SW) ---
    ! SW: nao precisa de SST, mas aplica mascara para evitar fluxos em terra
    call ESMF_FieldGet(is%f_swvdr_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2; do i = i1, i2
      if (associated(mask_atm)) then
        if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
          fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
        end if
      end if
      fptr(i,j) = swdn(i,j) * (1.0_ESMF_KIND_R8 - albedo_ocn) * f_vis_dir
    end do; end do

    call ESMF_FieldGet(is%f_swvdf_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2; do i = i1, i2
      if (associated(mask_atm)) then
        if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
          fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
        end if
      end if
      fptr(i,j) = swdn(i,j) * (1.0_ESMF_KIND_R8 - albedo_ocn) * f_vis_dif
    end do; end do

    call ESMF_FieldGet(is%f_swidr_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2; do i = i1, i2
      if (associated(mask_atm)) then
        if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
          fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
        end if
      end if
      fptr(i,j) = swdn(i,j) * (1.0_ESMF_KIND_R8 - albedo_ocn) * f_nir_dir
    end do; end do

    call ESMF_FieldGet(is%f_swidf_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2; do i = i1, i2
      if (associated(mask_atm)) then
        if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
          fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
        end if
      end if
      fptr(i,j) = swdn(i,j) * (1.0_ESMF_KIND_R8 - albedo_ocn) * f_nir_dif
    end do; end do

    ! --- Precipitacao e pressao: propagar mascara (nao zerados em terra
    !     pois MOM6 espera estes campos mesmo em pontos costeiros; porem
    !     zerados para consistencia com a mascara se o ponto for interior) ---
    ! Rain, snow, psl
    call ESMF_FieldGet(is%f_rain_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2; do i = i1, i2
      if (associated(mask_atm)) then
        if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
          fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
        end if
      end if
      fptr(i,j) = rain(i,j)
    end do; end do

    call ESMF_FieldGet(is%f_snow_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2; do i = i1, i2
      if (associated(mask_atm)) then
        if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
          fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
        end if
      end if
      fptr(i,j) = snow(i,j)
    end do; end do

    call ESMF_FieldGet(is%f_pslv_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2; do i = i1, i2
      if (associated(mask_atm)) then
        if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
          fptr(i,j) = 101325.0_ESMF_KIND_R8; cycle  ! padrao atmosferico em terra
        end if
      end if
      fptr(i,j) = psl(i,j)
    end do; end do

    ! --- duu10n ---
    call ESMF_FieldGet(is%f_duu10n_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2
      do i = i1, i2
        if (associated(mask_atm)) then
          if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
            fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
          end if
        end if
        fptr(i,j) = uas(i,j)**2 + vas(i,j)**2
      end do
    end do

    ! --- Fracao de gelo (threshold em SST) ---
    call ESMF_FieldGet(is%f_ifrac_atm, farrayPtr=fptr, rc=rc)
    do j = j1, j2
      do i = i1, i2
        if (associated(mask_atm)) then
          if (mask_atm(i,j) < MASK_BIN_THRESHOLD) then
            fptr(i,j) = 0.0_ESMF_KIND_R8; cycle
          end if
        end if
        sst_eff = merge(sst(i,j), SST_default, &
          sst(i,j) > 271.0_ESMF_KIND_R8 .and. sst(i,j) < 308.0_ESMF_KIND_R8)
        fptr(i,j) = merge(1.0_ESMF_KIND_R8, 0.0_ESMF_KIND_R8, &
          sst_eff <= 271.35_ESMF_KIND_R8)
      end do
    end do

    !==========================================================================
    ! 5. DIAGNOSTICO NetCDF de SST
    !    Escrito APOS o bulk, para que f_sst_atm ja contenha a SST real
    !    interpolada do OCN (regrid OCN->ATM feito na secao 3).
    !    Incrementa o contador e escreve a cada diag_freq passos.
    !==========================================================================
    is%diag_step = is%diag_step + 1
    if (is%diag_enabled .and. mod(is%diag_step, is%diag_freq) == 0) then
      call WriteDiagSST(gcomp, is, importState, currTime, is%diag_step, rc)
      if (rc /= ESMF_SUCCESS) then
        write(msg,'(A,I0)') 'MED: AVISO - WriteDiagSST falhou no passo ', is%diag_step
        call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_WARNING)
        rc = ESMF_SUCCESS
      end if
    end if    
    !==========================================================================
    ! 6. REGRID E EXPORTA PARA O OCEANO
    ! CORRECAO 3: RegridOrCopy agora tem ramo else explicito: se routehandles
    !   nao estiverem criados, copia direto da grade ATM interna para a grade
    !   OCN do exportState via ESMF_FieldSMM (ou copia simples). Isso evita
    !   que os campos exportados permanecam zerados silenciosamente.
    !==========================================================================
    call RegridOrCopy(is%f_taux_atm,   exportState, "Foxx_taux",      is, rc)
    call RegridOrCopy(is%f_tauy_atm,   exportState, "Foxx_tauy",      is, rc)
    call RegridOrCopy(is%f_sen_atm,    exportState, "Foxx_sen",       is, rc)
    call RegridOrCopy(is%f_evap_atm,   exportState, "Foxx_evap",      is, rc)
    call RegridOrCopy(is%f_lwnet_atm,  exportState, "Foxx_lwnet",     is, rc)
    call RegridOrCopy(is%f_swvdr_atm,  exportState, "Foxx_swnet_vdr", is, rc)
    call RegridOrCopy(is%f_swvdf_atm,  exportState, "Foxx_swnet_vdf", is, rc)
    call RegridOrCopy(is%f_swidr_atm,  exportState, "Foxx_swnet_idr", is, rc)
    call RegridOrCopy(is%f_swidf_atm,  exportState, "Foxx_swnet_idf", is, rc)
    call RegridOrCopy(is%f_rain_atm,   exportState, "Faxa_rain",      is, rc)
    call RegridOrCopy(is%f_snow_atm,   exportState, "Faxa_snow",      is, rc)
    call RegridOrCopy(is%f_pslv_atm,   exportState, "Sa_pslv",        is, rc)
    call RegridOrCopy(is%f_ifrac_atm,  exportState, "Si_ifrac",       is, rc)
    call RegridOrCopy(is%f_duu10n_atm, exportState, "So_duu10n",      is, rc)

    ! Atualizar timestamps do exportState
    call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(exportState, itemNameList=fieldNameList, rc=rc)
    do k = 1, fieldCount
      call ESMF_StateGet(exportState, itemName=trim(fieldNameList(k)), &
        field=field, rc=rc)
      call NUOPC_SetTimestamp(field, nextTime, rc=rc)
    end do
    deallocate(fieldNameList)

    call ESMF_LogWrite('MED: MediatorAdvance concluido', ESMF_LOGMSG_INFO)
  end subroutine MediatorAdvance

  !============================================================================
  ! GetFieldPtr - obtem ponteiro para campo (falha se nao existir)
  !============================================================================
  subroutine GetFieldPtr(state, name, ptr, rc)
    type(ESMF_State),            intent(in)    :: state
    character(len=*),            intent(in)    :: name
    real(ESMF_KIND_R8), pointer, intent(inout) :: ptr(:,:)
    integer,                     intent(out)   :: rc

    type(ESMF_Field) :: field
    !PK type(ESMF_StateItem_Flag) :: itemFlag
    integer :: localrc

    rc = ESMF_SUCCESS
    nullify(ptr)

    ! Verificar se o campo existe no state
    call ESMF_StateGet(state, trim(name), field, rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      rc = ESMF_FAILURE; return
    end if

    call ESMF_FieldGet(field, farrayPtr=ptr, rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      rc = ESMF_FAILURE; return
    end if

  end subroutine GetFieldPtr

  !============================================================================
  ! GetFieldPtrOptional - obtem ponteiro, nao falha se nao existir
  !============================================================================
  subroutine GetFieldPtrOptional(state, name, ptr, rc)
    type(ESMF_State),            intent(in)    :: state
    character(len=*),            intent(in)    :: name
    real(ESMF_KIND_R8), pointer, intent(inout) :: ptr(:,:)
    integer,                     intent(out)   :: rc

    type(ESMF_Field) :: field
    integer :: localrc

    rc = ESMF_SUCCESS
    nullify(ptr)

    ! Tentar obter o campo - se falhar, apenas retorna rc=ESMF_FAILURE
    call ESMF_StateGet(state, trim(name), field, rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      rc = ESMF_FAILURE; return
    end if

    call ESMF_FieldGet(field, farrayPtr=ptr, rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      rc = ESMF_FAILURE; return
    end if

    rc = ESMF_SUCCESS

  end subroutine GetFieldPtrOptional

  !============================================================================
  ! RegridOrCopy
  ! CORRECAO 3: ramo else adicionado para o caso rh_created = .false.
  !   Sem o else, os campos exportados ao OCN ficavam zerados/inalterados
  !   silenciosamente quando os routehandles nao tinham sido criados, o que
  !   causava fluxos incorretos no primeiro passo ou em caso de erro na IDC.
  !   Com o else, faz regrid on-the-fly via ESMF_FieldRegridStore temporario.
  !============================================================================
  subroutine RegridOrCopy(src_field, dst_state, dst_name, is, rc)
    type(ESMF_Field),        intent(inout) :: src_field
    type(ESMF_State),        intent(inout) :: dst_state
    character(len=*),        intent(in)    :: dst_name
    type(MED_InternalState), intent(inout) :: is
    integer,                 intent(out)   :: rc

    type(ESMF_Field) :: dst_field
    type(ESMF_RouteHandle) :: rh_tmp
    real(ESMF_KIND_R8), pointer :: dst_ptr(:,:)

    rc = ESMF_SUCCESS

    call ESMF_StateGet(dst_state, itemName=trim(dst_name), field=dst_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="RegridOrCopy: "//trim(dst_name), &
      line=__LINE__, file=__FILE__)) return

    if (is%rh_created) then
      call ESMF_FieldRegrid(src_field, dst_field, is%rh_atm2ocn, &
        zeroregion=ESMF_REGION_TOTAL, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="RegridOrCopy: falha no regrid de "//trim(dst_name), &
        line=__LINE__, file=__FILE__)) return

      ! Sanitizar NaNs
      call ESMF_FieldGet(dst_field, farrayPtr=dst_ptr, rc=rc)
      where (dst_ptr /= dst_ptr) dst_ptr = 0.0_ESMF_KIND_R8

    else
      ! Routehandle ainda nao disponivel: cria um temporario para este campo
      ! CORRECAO 5 (cont.): BILINEAR tambem no fallback
      call ESMF_FieldRegridStore( &
        srcField       = src_field,    &
        dstField       = dst_field,    &
        routehandle    = rh_tmp,       &
        regridmethod   = ESMF_REGRIDMETHOD_BILINEAR, &
        unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
        rc             = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="RegridOrCopy fallback: falha store "//trim(dst_name), &
        line=__LINE__, file=__FILE__)) return

      call ESMF_FieldRegrid(src_field, dst_field, rh_tmp, &
        zeroregion=ESMF_REGION_TOTAL, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="RegridOrCopy fallback: falha regrid "//trim(dst_name), &
        line=__LINE__, file=__FILE__)) return

      call ESMF_RouteHandleDestroy(rh_tmp, nogarbage=.true., rc=rc)

      call ESMF_FieldGet(dst_field, farrayPtr=dst_ptr, rc=rc)
      where (dst_ptr /= dst_ptr) dst_ptr = 0.0_ESMF_KIND_R8
    end if

  end subroutine RegridOrCopy
  !============================================================================
  ! WriteDiagSST
  !
  ! Escreve dois arquivos NetCDF de diagnostico por chamada:
  !
  !   med_sst_atm_NNNNN.nc  : SST na grade ATM intermediaria (lat-lon regular)
  !                            Variaveis: lon(nx,ny), lat(nx,ny), sst(nx,ny)
  !
  !   med_sst_voronoi_NNNNN.nc : SST exportada para a malha Voronoi do MPAS
  !                              Obtida do exportState do MPAS (campo "So_t" ou
  !                              campo de SST no importState do MPAS se disponivel)
  !                              Variaveis: lon(ncells), lat(ncells), sst(ncells)
  !
  ! Ambos os arquivos incluem tempo como variavel escalar (segundos desde epoch).
  ! A escrita e feita apenas no PET 0 (root) via coleta com ESMF_FieldGather.
  !
  ! NOTA sobre a malha Voronoi:
  !   O ESMF nao exporta diretamente as coordenadas do Mesh. Para obte-las,
  !   tentamos ler latCell/lonCell do campo "Sa_u10m_mpas" (que esta na grade
  !   MPAS). Se nao disponivel, escrevemos apenas os valores sem coordenadas.
  !
  ! @param gcomp     componente mediador
  ! @param is        estado interno (contem f_sst_atm e grade ATM)
  ! @param importState estado de import (para obter SST na grade MPAS)
  ! @param currTime  tempo atual (para timestamp no arquivo)
  ! @param step      numero do passo (usado no nome do arquivo)
  ! @param rc        codigo de retorno
  !============================================================================
  ! WriteDiagSST -> VERSAO CORRIGIDA
  !
  ! Escreve diagnostico de SST em NetCDF a cada diag_freq passos.
  !   Arquivo 1: med_sst_atm_NNNNN.nc  ? SST na grade ATM 2D (lat-lon regular)
  !   Arquivo 2: med_sst_voronoi_NNNNN.nc ? SST regridada para a malha MPAS
  !
  ! Correcoes em relacao a versao anterior:
  !   - ESMF_FieldGet sem argumentos keyword invalidos (computationalLBound
  !     nao e valido nessa forma; usar ESMF_GridGet + lbound/ubound do ptr)
  !   - ESMF_GridGet sem minIndex/maxIndex (usar totalLBound/totalUBound)
  !   - ESMF_VMAllReduce com array Fortran explicitamente alocado
  !   - Sem gridToFieldMap (argumento invalido em ESMF_FieldGet)
  !   - Sem block construct (nao suportado em gfortran antigo)
  !============================================================================
  ! WriteDiagSST ? VERSAO CORRIGIDA
  !
  ! Escreve diagnostico de SST em NetCDF a cada diag_freq passos.
  !   Arquivo 1: med_sst_atm_NNNNN.nc  ? SST na grade ATM 2D (lat-lon regular)
  !   Arquivo 2: med_sst_voronoi_NNNNN.nc ? SST regridada para a malha MPAS
  !
  ! Correcoes em relacao a versao anterior:
  !   - ESMF_FieldGet sem argumentos keyword invalidos (computationalLBound
  !     nao e valido nessa forma; usar ESMF_GridGet + lbound/ubound do ptr)
  !   - ESMF_GridGet sem minIndex/maxIndex (usar totalLBound/totalUBound)
  !   - ESMF_VMAllReduce com array Fortran explicitamente alocado
  !   - Sem gridToFieldMap (argumento invalido em ESMF_FieldGet)
  !   - Sem block construct (nao suportado em gfortran antigo)
  !============================================================================
  !============================================================================
  ! WriteDiagSST
  !
  ! Escreve diagnostico de SST em dois arquivos NetCDF por chamada:
  !
  !   med_sst_atm_NNNNN.nc     : SST na grade ATM intermediaria (lat-lon regular)
  !   med_sst_voronoi_NNNNN.nc : SST interpolada para a malha Voronoi MPAS
  !
  ! ESTRATEGIA DE I/O (compativel com ESMF 8.9.0 + ESMF_MOAB=enabled):
  !
  !   Arquivo ATM: usa ESMF_FieldWrite (paralelo, MPI-aware, sem gather manual).
  !     Evita completamente ESMF_FieldGather para a grade ATM 2D.
  !     Gera um arquivo NetCDF com variaveis lon/lat/sst com dimensoes globais.
  !
  !   Arquivo Voronoi: usa ESMF_FieldGather com buffer alocado em TODOS os PETs
  !     com tamanho global obtido via ESMF_VMBroadcast do PET 0.
  !     A escrita NetCDF e feita somente no PET 0.
  !============================================================================
  subroutine WriteDiagSST(gcomp, is, importState, currTime, step, rc)
    type(ESMF_GridComp),     intent(in)    :: gcomp
    type(MED_InternalState), intent(inout) :: is
    type(ESMF_State),        intent(in)    :: importState
    type(ESMF_Time),         intent(in)    :: currTime
    integer,                 intent(in)    :: step
    integer,                 intent(out)   :: rc

    integer :: localPet, petCount, localrc
    integer :: nx_g, ny_g
    integer :: ncid, dimid_x, dimid_y, dimid_n
    integer :: varid_lon, varid_lat, varid_sst, varid_time
    integer :: yy, mm, dd, hh, mn, ss
    real(ESMF_KIND_R8) :: time_sec
    character(len=512) :: filename, filename_nc
    character(len=8)   :: step_str

    ! Buffers para arquivo Voronoi
    ! Buffers globais (alocados apenas no PET 0)
    real(ESMF_KIND_R8), allocatable :: sst2d(:,:), lon2d(:,:), lat2d(:,:)
    real(ESMF_KIND_R8), allocatable :: sst1d(:), lon1d(:), lat1d(:)
    real(ESMF_KIND_R8), allocatable :: elem_coords(:)

    ! Campos e grids ESMF
    type(ESMF_Field)       :: f_sst_mpas, f_sst_voronoi, f_lon_atm, f_lat_atm
    type(ESMF_Mesh)        :: mpas_mesh
    type(ESMF_RouteHandle) :: rh_atm2mpas
    type(ESMF_VM)          :: vm
    real(ESMF_KIND_R8), pointer :: cx(:,:), cy(:,:)
    integer :: n_local, n_global_arr(1), n_out_arr(1), n_global
    integer :: k

    rc = ESMF_SUCCESS

    !-------------------------------------------------------------------
    ! Obter VM e identificar PET 0
    !-------------------------------------------------------------------
    call ESMF_VMGetGlobal(vm, rc=localrc)
    if (localrc /= ESMF_SUCCESS) return
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=localrc)
    if (localrc /= ESMF_SUCCESS) return

    !-------------------------------------------------------------------
    ! Tempo atual
    !-------------------------------------------------------------------
    call ESMF_TimeGet(currTime, yy=yy, mm=mm, dd=dd, h=hh, m=mn, s=ss, rc=localrc)
    if (localrc /= ESMF_SUCCESS) return
    time_sec = real(ss + 60*mn + 3600*hh, ESMF_KIND_R8)
    write(step_str, '(I5.5)') step

    !=======================================================================
    ! ARQUIVO 1: SST na grade ATM via ESMF_FieldWrite (sem gather manual)
    ! ESMF_FieldWrite e MPI-aware: todos os PETs participam e o resultado
    ! e um NetCDF paralelo com os dados globais. Nao requer alocacao de
    ! buffer global nem ESMF_FieldGather.
    !=======================================================================

    ! Garantir que o diretorio existe (chamada coletiva no PET 0)
    if (localPet == 0) &
      call execute_command_line('mkdir -p '//trim(is%diag_dir), wait=.true.)
    call ESMF_VMBarrier(vm, rc=localrc)

    write(filename_nc, '(A,"/med_sst_atm_",A,".nc")') trim(is%diag_dir), trim(step_str)

    ! ESMF_FieldWrite: escreve campo distribuido para arquivo NetCDF
    ! O arquivo resultante tem dimensoes globais e e legivel por xarray/ncdump
    call ESMF_FieldWrite(is%f_sst_atm, &
      fileName       = trim(filename_nc), &
      variableName   = "sst", &
      status         = ESMF_FILESTATUS_REPLACE, &
      timeslice      = 1, &
      rc             = localrc)
    if (localrc /= ESMF_SUCCESS) then
      call ESMF_LogWrite('MED WriteDiagSST: ESMF_FieldWrite SST ATM falhou; skip', &
        ESMF_LOGMSG_WARNING)
      goto 200
    end if

    ! Escrever tambem os campos de coordenadas no mesmo arquivo
    ! Criar campos temporarios lon/lat na grade ATM
    f_lon_atm = ESMF_FieldCreate(is%atm_grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, name="lon_atm_diag", rc=localrc)
    if (localrc /= ESMF_SUCCESS) goto 200
    f_lat_atm = ESMF_FieldCreate(is%atm_grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, name="lat_atm_diag", rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      call ESMF_FieldDestroy(f_lon_atm, rc=localrc); goto 200
    end if

    call ESMF_FieldGet(f_lon_atm, farrayPtr=cx, rc=localrc)
    if (localrc == ESMF_SUCCESS) &
      call ESMF_GridGetCoord(is%atm_grid, coordDim=1, &
        staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=cx, rc=localrc)

    call ESMF_FieldGet(f_lat_atm, farrayPtr=cy, rc=localrc)
    if (localrc == ESMF_SUCCESS) &
      call ESMF_GridGetCoord(is%atm_grid, coordDim=2, &
        staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=cy, rc=localrc)

    call ESMF_FieldWrite(f_lon_atm, &
      fileName     = trim(filename_nc), &
      variableName = "lon", &
      status       = ESMF_FILESTATUS_OLD, &
      rc           = localrc)
    call ESMF_FieldWrite(f_lat_atm, &
      fileName     = trim(filename_nc), &
      variableName = "lat", &
      status       = ESMF_FILESTATUS_OLD, &
      rc           = localrc)

    call ESMF_FieldDestroy(f_lon_atm, rc=localrc)
    call ESMF_FieldDestroy(f_lat_atm, rc=localrc)

    call ESMF_LogWrite('MED WriteDiagSST: '//trim(filename_nc)//' escrito (ESMF_FieldWrite)', &
      ESMF_LOGMSG_INFO)

    !=======================================================================
    ! ARQUIVO 2: SST na malha Voronoi do MPAS
    ! Regrid bilinear ATM Grid -> MPAS Mesh + gather com buffer global
    ! As dimensoes globais sao obtidas via ESMF_VMBroadcast do PET 0
    !=======================================================================
    200 continue

    ! Verificar existencia de Sa_u10m_mpas via itemNameList
    block
      integer :: item_count, n_item
      character(len=64), allocatable :: item_names(:)
      logical :: found_it
      found_it = .false.
      call ESMF_StateGet(importState, itemCount=item_count, rc=localrc)
      if (localrc == ESMF_SUCCESS .and. item_count > 0) then
        allocate(item_names(item_count))
        call ESMF_StateGet(importState, itemNameList=item_names, rc=localrc)
        if (localrc == ESMF_SUCCESS) then
          do n_item = 1, item_count
            if (trim(item_names(n_item)) == "Sa_u10m_mpas") then
              found_it = .true.; exit
            end if
          end do
        end if
        deallocate(item_names)
      end if
      if (.not. found_it) then
        call ESMF_LogWrite('MED WriteDiagSST: Sa_u10m_mpas ausente; skip Voronoi', &
          ESMF_LOGMSG_INFO)
        goto 300
      end if
    end block

    call ESMF_StateGet(importState, itemName="Sa_u10m_mpas", &
      field=f_sst_mpas, rc=localrc)
    if (localrc /= ESMF_SUCCESS) goto 300

    call ESMF_FieldGet(f_sst_mpas, mesh=mpas_mesh, rc=localrc)
    if (localrc /= ESMF_SUCCESS) goto 300

    f_sst_voronoi = ESMF_FieldCreate(mesh=mpas_mesh, &
      typekind=ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
      name="sst_voronoi_diag", rc=localrc)
    if (localrc /= ESMF_SUCCESS) goto 300

    ! interpolar
    call ESMF_FieldRegridStore( &
      srcField       = is%f_sst_atm,              &
      dstField       = f_sst_voronoi,             &
      routehandle    = rh_atm2mpas,               &
      regridmethod   = ESMF_REGRIDMETHOD_BILINEAR, &
      unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
      rc             = localrc)
    if (localrc /= ESMF_SUCCESS) then
      call ESMF_FieldDestroy(f_sst_voronoi, rc=localrc); goto 300
    end if

    call ESMF_FieldRegrid(is%f_sst_atm, f_sst_voronoi, rh_atm2mpas, &
      zeroregion=ESMF_REGION_TOTAL, rc=localrc)
    call ESMF_RouteHandleDestroy(rh_atm2mpas, nogarbage=.true., rc=localrc)

    ! Obter numero total de elementos via reducao MPI
    call ESMF_MeshGet(mpas_mesh, numOwnedElements=n_local, rc=localrc)
    if (localrc /= ESMF_SUCCESS) n_local = 0
    n_global_arr(1) = n_local
    call ESMF_VMAllReduce(vm, n_global_arr, n_out_arr, 1, &
      ESMF_REDUCE_SUM, rc=localrc)
    n_global = merge(n_out_arr(1), 0, localrc == ESMF_SUCCESS)

    if (n_global <= 0) then
      call ESMF_FieldDestroy(f_sst_voronoi, rc=localrc); goto 300
    end if

    ! Obter tamanho global via PET 0 e broadcast para todos os PETs
    ! ESMF_FieldGather exige que sst1d tenha tamanho global em TODOS os PETs
    ! PET 0: ja tem n_global correto; outros PETs: recebem via broadcast
    allocate(sst1d(n_global))    ! todos os PETs
    if (localPet == 0) then
      allocate(lon1d(n_global), lat1d(n_global))
    end if

    call ESMF_FieldGather(f_sst_voronoi, sst1d, rootPet=0, rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      deallocate(sst1d)
      if (localPet == 0) deallocate(lon1d, lat1d)
      call ESMF_FieldDestroy(f_sst_voronoi, rc=localrc); goto 300
    end if

    ! Coordenadas e escrita NetCDF apenas no PET 0
    if (localPet == 0) then
      allocate(elem_coords(2 * n_global))
      call ESMF_MeshGet(mpas_mesh, ownedElemCoords=elem_coords, rc=localrc)
      if (localrc == ESMF_SUCCESS) then
        do k = 1, n_global
          lon1d(k) = elem_coords(2*k - 1)
          lat1d(k) = elem_coords(2*k)
        end do
      else
        lon1d = 0.0_ESMF_KIND_R8; lat1d = 0.0_ESMF_KIND_R8
      end if
      deallocate(elem_coords)

      write(filename_nc, '(A,"/med_sst_voronoi_",A,".nc")') &
        trim(is%diag_dir), trim(step_str)

      localrc = nf90_create(trim(filename_nc), NF90_CLOBBER, ncid)
      if (localrc == NF90_NOERR) then
        localrc = nf90_def_dim(ncid, "ncells", n_global, dimid_n)
        localrc = nf90_put_att(ncid, NF90_GLOBAL, "title", &
          "MED SST na malha Voronoi MPAS/MONAN")
        localrc = nf90_put_att(ncid, NF90_GLOBAL, "step", step)
        write(filename, '(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":",I2.2,":",I2.2)') &
          yy, mm, dd, hh, mn, ss
        localrc = nf90_put_att(ncid, NF90_GLOBAL, "valid_time", trim(filename))

        localrc = nf90_def_var(ncid, "time", NF90_DOUBLE, varid=varid_time)
        localrc = nf90_put_att(ncid, varid_time, "units", "seconds within day")
        localrc = nf90_def_var(ncid, "lon", NF90_DOUBLE, dimids=(/dimid_n/), varid=varid_lon)
        localrc = nf90_put_att(ncid, varid_lon, "units",     "degrees_east")
        localrc = nf90_put_att(ncid, varid_lon, "long_name", "longitude celula Voronoi")
        localrc = nf90_def_var(ncid, "lat", NF90_DOUBLE, dimids=(/dimid_n/), varid=varid_lat)
        localrc = nf90_put_att(ncid, varid_lat, "units",     "degrees_north")
        localrc = nf90_put_att(ncid, varid_lat, "long_name", "latitude celula Voronoi")
        localrc = nf90_def_var(ncid, "sst", NF90_DOUBLE, dimids=(/dimid_n/), varid=varid_sst)
        localrc = nf90_put_att(ncid, varid_sst, "units",         "K")
        localrc = nf90_put_att(ncid, varid_sst, "long_name",     "sea surface temperature")
        localrc = nf90_put_att(ncid, varid_sst, "standard_name", "sea_surface_temperature")
        localrc = nf90_put_att(ncid, varid_sst, "_FillValue",    0.0_ESMF_KIND_R8)
        localrc = nf90_put_att(ncid, varid_sst, "coordinates",   "lon lat")
        localrc = nf90_enddef(ncid)

        localrc = nf90_put_var(ncid, varid_time, time_sec)
        localrc = nf90_put_var(ncid, varid_lon,  lon1d)
        localrc = nf90_put_var(ncid, varid_lat,  lat1d)
        localrc = nf90_put_var(ncid, varid_sst,  sst1d)
        localrc = nf90_close(ncid)

        write(filename_nc, '(A,"/med_sst_voronoi_",A,".nc")') &
          trim(is%diag_dir), trim(step_str)
        call ESMF_LogWrite('MED WriteDiagSST: '//trim(filename_nc)//' escrito', &
          ESMF_LOGMSG_INFO)
      end if

      deallocate(lon1d, lat1d)
    end if

    deallocate(sst1d)
    call ESMF_FieldDestroy(f_sst_voronoi, rc=localrc)

    300 continue
    rc = ESMF_SUCCESS

  end subroutine WriteDiagSST
  
end module MED_cap_mod
