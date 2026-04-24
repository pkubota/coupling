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
!   4. Busca por Sa_u10m_mpas (MPAS) em vez de Sa_u10m (DATM) no IDC          !
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

  ! Limiar de binarizacao da mascara apas regrid NEAREST_STOD.
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

    ! RouteHandles
    type(ESMF_RouteHandle) :: rh_atm2ocn   ! ATM -> OCN
    type(ESMF_RouteHandle) :: rh_ocn2atm   ! OCN -> ATM

    ! Mascara oceano/continente
    !PK real(ESMF_KIND_R8), allocatable :: ocn_mask_atm(:,:)
    
    logical :: rh_created    = .false.
    logical :: use_mpas_atm  = .false.  ! controlado por atributo NUOPC "use_mpas_atm"
  end type MED_InternalState

  type :: MED_InternalStateWrapper
    type(MED_InternalState), pointer :: wrap => null()
  end type MED_InternalStateWrapper

  ! Campos de import do MPAS (primario) - com sufixo _mpas
  integer, parameter :: n_import_mpas = 10   ! +1 para Sa_lfrac_mpas
  character(len=32), parameter :: import_mpas_names(n_import_mpas) = [ &
    "Sa_u10m_mpas  ", "Sa_v10m_mpas  ", "Sa_tbot_mpas  ", "Sa_shum_mpas  ", "Sa_pslv_mpas  ", &
    "Faxa_swdn_mpas", "Faxa_lwdn_mpas", "Faxa_rain_mpas", "Faxa_snow_mpas", &
    "Sa_lfrac_mpas " ]   ! NOVO: fracao de terra do MONAN para reconciliacao de mascara

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
    !PK ! Inicializar todos os campos lógicos do InternalState
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
      do n = 1, n_import_mpas
        call NUOPC_Advertise(importState, StandardName=trim(import_mpas_names(n)), &
          TransferOfferGeomObject="cannot provide", &
          SharePolicyField="share", rc=localrc)
        ! Sa_lfrac_mpas e opcional: nao abortar se o MONAN nao o oferecer
        if (localrc /= ESMF_SUCCESS .and. &
            trim(import_mpas_names(n)) == "Sa_lfrac_mpas") then
          call ESMF_LogWrite( &
            'MED: AVISO - Sa_lfrac_mpas nao anunciado (MONAN nao exporta lndfrac)', &
            ESMF_LOGMSG_WARNING)
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
    call NUOPC_Advertise(importState, StandardName="So_t", &
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
    real(ESMF_KIND_R8), allocatable :: ocn_lon(:,:), ocn_lat(:,:)
    logical             :: isPresent, isSet
    character(len=256)    :: attr_val
    character(len=256)  :: ocn_hgrid_file, ocn_mask_file
    character(len=256)  :: line_tag
    integer :: clbnd(2), cubnd(2)   ! compute lower/upper bound (global index)
    integer :: ig, jg               ! indices globais no grid OCN
      
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
      farrayPtr=coordX, rc=rc)
    do j = lbound(coordX,2), ubound(coordX,2)
      do i = lbound(coordX,1), ubound(coordX,1)
        coordX(i,j) = (i-1)*(360.0_ESMF_KIND_R8/nx_atm) + &
                      (360.0_ESMF_KIND_R8/nx_atm)*0.5_ESMF_KIND_R8
      end do
    end do

    call ESMF_GridGetCoord(atm_grid, coordDim=2, staggerloc=ESMF_STAGGERLOC_CENTER, &
      farrayPtr=coordY, rc=rc)
    do j = lbound(coordY,2), ubound(coordY,2)
      do i = lbound(coordY,1), ubound(coordY,1)
        coordY(i,j) = -90.0_ESMF_KIND_R8 + (j-1)*(180.0_ESMF_KIND_R8/ny_atm) + &
                      (180.0_ESMF_KIND_R8/ny_atm)/2.0_ESMF_KIND_R8
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

    ocn_grid = ESMF_GridCreateNoPeriDim(minIndex=(/1,1/), maxIndex=(/nx_ocn, ny_ocn/), &
      indexflag=ESMF_INDEX_GLOBAL, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha ao criar grade OCN", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_GridAddCoord(ocn_grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! Ler coordenadas reais do supergrid ocean_hgrid.nc
    ! O supergrid tem dimensoes (nyp, nxp) = (2*ny+1, 2*nx+1).
    ! Pontos T (center) = indices pares do supergrid (base 0):
    !   lon_T(i,j) = x_sg(2*j, 2*i)   i=0..nx-1, j=0..ny-1
    !   lat_T(i,j) = y_sg(2*j, 2*i)
    ! Em Fortran (base 1): x_sg(2*i-1, 2*j-1)
    !--------------------------------------------------------------------------
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

    ! -----------------------------------------------------------------------
    ! NOVO: Carregar mascara oceano-continente do MOM6.
    ! LoadOceanMask le ocean_static.nc (variavel "wet"), realiza regrid
    ! NEAREST_STOD OCN->ATM, binariza e armazena em is%f_ocn_mask_atm.
    ! Deve ser chamado APOS is%ocn_grid e is%atm_grid estarem definidos.
    ! -----------------------------------------------------------------------
    call LoadOceanMask(gcomp, is, trim(ocn_mask_file), nx_ocn, ny_ocn, rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite( &
        'MED: AVISO - LoadOceanMask falhou; prosseguindo sem mascara explicita', &
        ESMF_LOGMSG_WARNING)
      rc = ESMF_SUCCESS   ! nao abortar ? mascara e auxiliar, nao critica
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
      'MED: mascara oceano-continente carregada e interpolada para grade ATM', &
      ESMF_LOGMSG_INFO)

  end subroutine LoadOceanMask

  !============================================================================
  ! ReconcileCoastalMask
  !
  ! Reconcilia a mascara de terra/oceano do MED com a fracao de terra
  ! exportada pelo MONAN (Sa_lfrac_mpas).
  !
  ! Estrategia conservadora (interseção):
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

    ! Interseção conservadora: zera mascara onde MONAN vê terra
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
        'MED ReconcileCoastalMask: ', n_land_masked, ' de ', n_total, &
        ' pontos ATM reclassificados como terra pela mascara MONAN'
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
    type(ESMF_Field)         :: atm_field, ocn_field, exp_field
    type(ESMF_Grid)          :: atm_grid!PK
    type(MED_InternalStateWrapper) :: iswrap
    type(MED_InternalState), pointer :: is
    integer :: fieldCount, i, localrc
    character(len=64), allocatable :: fieldNameList(:)
    real(ESMF_KIND_R8), pointer :: fptr(:,:)

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

    ! Routehandle OCN -> ATM (So_t ja esta na grade OCN)
    ! Criar routehandle OCN -> ATM
    ! So_t esta agora corretamente na grade OCN (ver InitializeRealize)
    call ESMF_StateGet(importState, itemName="So_t", field=ocn_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MED: falha So_t", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_FieldRegridStore( &
      srcField       = ocn_field,       &
      dstField       = is%f_sst_atm,    &
      routehandle    = is%rh_ocn2atm,   &
      regridmethod   = ESMF_REGRIDMETHOD_BILINEAR, &
      unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
      rc             = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="MED: falha FieldRegridStore OCN->ATM", &
      line=__LINE__, file=__FILE__)) return

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

    ! Campos da ATM  do MPAS (primario)
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
    real(ESMF_KIND_R8), pointer :: sst(:,:), fptr(:,:)

    ! NOVO: ponteiro para mascara oceano-continente interpolada para ATM
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
    ! 1. OBTER CAMPOS DA ATM DO MPAS (PRIMARIO)
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

      call ESMF_LogWrite('MED: Usando MPAS/MONAN como fonte atmosferica', &
        ESMF_LOGMSG_INFO)
    end if

    i1 = lbound(uas,1); i2 = ubound(uas,1)
    j1 = lbound(uas,2); j2 = ubound(uas,2)

    !==========================================================================
    ! 3. SST: regrid OCN -> ATM (So_t esta agora na grade OCN)
    !==========================================================================
    if (is%rh_created) then
      call ESMF_StateGet(importState, itemName="So_t", field=field, rc=rc)
      call ESMF_FieldRegrid(field, is%f_sst_atm, is%rh_ocn2atm, &
        zeroregion=ESMF_REGION_TOTAL, rc=rc)
      call ESMF_FieldGet(is%f_sst_atm, farrayPtr=sst, rc=rc)
    else
      ! Routehandles nao criados: usa SST padrao (ja preenchido em InitializeRealize)
      call ESMF_FieldGet(is%f_sst_atm, farrayPtr=sst, rc=rc)
      !PK call GetFieldPtr(importState, "So_t", sst, rc)
    end if

    !==========================================================================
    ! 3. NOVO: Reconciliar mascara costeira MONAN x MOM6
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
    if (is%mask_loaded) then
      call ESMF_FieldGet(is%f_ocn_mask_atm, farrayPtr=mask_atm, rc=rc)
      if (rc /= ESMF_SUCCESS) nullify(mask_atm)
    end if
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

    ! --- Bandas de onda curta ---
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
    ! 5. REGRID E EXPORTA PARA O OCEANO
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

end module MED_cap_mod
