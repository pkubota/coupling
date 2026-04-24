!==============================================================================!
! esm.F90 - ESM Driver Component (MPAS + DATM + MOM6+SIS2)                     !
!==============================================================================!
! Arquitetura:                                                                 !
!   MPAS (primario) + DATM (fallback) -> MED -> OCN                           !
!                                                                              !
! Conectores (5 no total):                                                     !
!   MPAS -> MED   : campos atmosfericos do MPAS (primario)                     !
!   DATM -> MED   : campos atmosfericos do JRA55 (fallback)                    !
!   OCN  -> MED   : SST/So_s/So_u/So_v para bulk formula (lag 1 step)          !
!   MED  -> OCN   : fluxos calculados (Foxx_*, Faxa_*, etc.)                   !
!   OCN  -> MPAS  : SST/correntes diretamente para pool sfc_input do MPAS      !
!                                                                              !
! Nota sobre OCN->MPAS:                                                        !
!   O MPAS precisa de SST antes de cada timestep de fisica de superficie.      !
!   O conector OCN->MPAS usa regrid bilinear MOM6_grid -> MPAS_mesh.           !
!   Os campos conectados sao: So_t, So_s, So_u, So_v.                          !
!   O MPAS_cap anunciou esses campos com "cannot provide" + "share",           !
!   portanto o conector reusa a geometria já realizada pelo MOM_cap.            !
!==============================================================================!
module ESM
  !use ESMF
  use ESMF,  only: ESMF_GridComp,ESMF_Clock,ESMF_Time,ESMF_GridComp
  use ESMF,  only: ESMF_LogFoundError,ESMF_GridCompGet,ESMF_ClockGet,ESMF_TimeGet
  use ESMF,  only: ESMF_LOGERR_PASSTHRU,ESMF_LOGMSG_INFO,ESMF_SUCCESS,ESMF_LogWrite
  use NUOPC, only: NUOPC_FreeFormatCreate, NUOPC_FreeFormat, NUOPC_FreeFormatDestroy
  use NUOPC, only: NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd
  use NUOPC, only: NUOPC_FieldDictionarySetAutoAdd
  use NUOPC, only: NUOPC_CompDerive, NUOPC_CompSpecialize
  use NUOPC_Driver, &
    driver_routine_SS             => SetServices,            &
    driver_label_SetModelServices => label_SetModelServices, &
    driver_label_SetRunSequence   => label_SetRunSequence
  ! Connector NUOPC padrao
  use NUOPC_Connector, only: CPL_SetServices => SetServices
  ! Caps dos componentes
  use MPAS_cap_mod,  only: MPAS_SetServices => SetServices
  use DATM_cap_mod,  only: DATM_SetServices => SetServices
  use MED_cap_mod,   only: MED_SetServices  => SetServices
  use MOM_cap_mod,   only: OCN_SetServices  => SetServices
  implicit none
  private
  public :: SetServices

  ! Rotulos dos componentes
  character(len=*), parameter :: MPAS_LABEL = "MPAS"
  character(len=*), parameter :: DATM_LABEL = "DATM"
  character(len=*), parameter :: MED_LABEL  = "MED"
  character(len=*), parameter :: OCN_LABEL  = "OCN"

  !----------------------------------------------------------------------------
  ! CHAVE DE CONTROLE: escolhe a fonte atmosferica em tempo de compilacao.
  !
  !   use_mpas_atm = .true.  ? MPAS-ATM como fonte primária (producao)
  !   use_mpas_atm = .false. ? DATM/JRA55 como fonte (desenvolvimento/fallback)
  !
  ! Esta flag e passada ao MED_cap via atributo NUOPC "use_mpas_atm" e
  ! tambem controla quais componentes e conectores sao ativados na RunSequence.
  !----------------------------------------------------------------------------
  logical, parameter :: use_mpas_atm = .false.

contains

  !============================================================================
  ! SetServices
  !============================================================================
  subroutine SetServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    rc = ESMF_SUCCESS

    call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(driver, &
      specLabel=driver_label_SetModelServices, &
      specRoutine=SetModelServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(driver, &
      specLabel=driver_label_SetRunSequence, &
      specRoutine=SetRunSequence, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

  !============================================================================
  ! SetModelServices - registra componentes e conectores
  !============================================================================
  subroutine SetModelServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    type(ESMF_GridComp)  :: mpasComp, datmComp, medComp, ocnComp
    integer              :: petCount, i
    integer, allocatable :: petList(:)

    type(ESMF_Clock)  :: driverClock
    type(ESMF_Time)   :: stopTime_local
    integer           :: syy, smm, sdd, sh, sm, ss
    integer           :: stop_ymd_int, stop_tod_int
    character(len=16) :: stop_ymd_str, stop_tod_str

    rc = ESMF_SUCCESS

    ! AutoAdd necessario para nomes customizados (Sa_u10m_mpas, So_t, etc.)
    call NUOPC_FieldDictionarySetAutoAdd(.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompGet(driver, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    allocate(petList(petCount))
    petList = [(i-1, i=1,petCount)]

    !--------------------------------------------------------------------------
    ! MPAS (primario) - registrado apenas quando use_mpas_atm=.true.
    !--------------------------------------------------------------------------
    if (use_mpas_atm) then
      call NUOPC_DriverAddComp(driver,                          &
        compLabel              = MPAS_LABEL,                    &
        compSetServicesRoutine = MPAS_SetServices,              &
        petList                = petList,                       &
        comp                   = mpasComp,                      &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call NUOPC_CompAttributeSet(mpasComp, name="Verbosity", value="high",  rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_CompAttributeSet(mpasComp, name="DumpFields", value="false", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    !--------------------------------------------------------------------------
    ! DATM (fallback) - registrado apenas quando use_mpas_atm=.false.
    !--------------------------------------------------------------------------
    if (.not. use_mpas_atm) then
      call NUOPC_DriverAddComp(driver,                          &
        compLabel              = DATM_LABEL,                    &
        compSetServicesRoutine = DATM_SetServices,              &
        petList                = petList,                       &
        comp                   = datmComp,                      &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call NUOPC_CompAttributeSet(datmComp, name="Verbosity", value="high",  rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_CompAttributeSet(datmComp, name="DumpFields", value="false", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    !--------------------------------------------------------------------------
    ! MED: Mediador
    !--------------------------------------------------------------------------
    call NUOPC_DriverAddComp(driver,                          &
      compLabel              = MED_LABEL,                     &
      compSetServicesRoutine = MED_SetServices,               &
      petList                = petList,                       &
      comp                   = medComp,                       &
      rc                     = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompAttributeSet(medComp, name="Verbosity", value="high", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--- Passa flag use_mpas_atm para o MED como atributo NUOPC ---
    ! O MED_cap le este atributo em InitializeAdvertise e InitializeRealize
    ! para decidir qual fonte atmosferica usar em MediatorAdvance.
    call NUOPC_CompAttributeAdd(medComp, attrList=(/"use_mpas_atm"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (use_mpas_atm) then
      call NUOPC_CompAttributeSet(medComp, name="use_mpas_atm", value="true",  rc=rc)
    else
      call NUOPC_CompAttributeSet(medComp, name="use_mpas_atm", value="false", rc=rc)
    end if
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! OCN: MOM6+SIS2
    !--------------------------------------------------------------------------
    call NUOPC_DriverAddComp(driver,                          &
      compLabel              = OCN_LABEL,                     &
      compSetServicesRoutine = OCN_SetServices,               &
      petList                = petList,                       &
      comp                   = ocnComp,                       &
      rc                     = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompAttributeSet(ocnComp, name="Verbosity",     value="high",  rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(ocnComp, name="DumpFields",    value="false", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(ocnComp, name="ProfileMemory", value="false", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(ocnComp, name="timeStampValidation", &
      value="false", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(ocnComp, name="restart_n", value="0", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--- Propaga stopTime para o OCN ---
    call ESMF_GridCompGet(driver, clock=driverClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_ClockGet(driverClock, stopTime=stopTime_local, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_TimeGet(stopTime_local, yy=syy, mm=smm, dd=sdd, &
      h=sh, m=sm, s=ss, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    stop_ymd_int = syy*10000 + smm*100 + sdd
    stop_tod_int = sh*3600   + sm*60   + ss
    write(stop_ymd_str,'(i0)') stop_ymd_int
    write(stop_tod_str,'(i0)') stop_tod_int

    call NUOPC_CompAttributeSet(ocnComp, name="stop_ymd", &
      value=trim(stop_ymd_str), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(ocnComp, name="stop_tod", &
      value=trim(stop_tod_str), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! Conectores condicionais: MPAS ou DATM -> MED
    !--------------------------------------------------------------------------
    if (use_mpas_atm) then
      ! MPAS -> MED : campos atmosfericos da malha Voronoi
      call NUOPC_DriverAddComp(driver,                        &
        srcCompLabel           = MPAS_LABEL,                  &
        dstCompLabel           = MED_LABEL,                   &
        compSetServicesRoutine = CPL_SetServices,             &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      ! DATM -> MED : campos JRA55 (grade regular 640x320)
      call NUOPC_DriverAddComp(driver,                        &
        srcCompLabel           = DATM_LABEL,                  &
        dstCompLabel           = MED_LABEL,                   &
        compSetServicesRoutine = CPL_SetServices,             &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    !--------------------------------------------------------------------------
    ! Conector MED -> OCN : fluxos bulk -> MOM6 (sempre ativo)
    !--------------------------------------------------------------------------
    call NUOPC_DriverAddComp(driver,                          &
      srcCompLabel           = MED_LABEL,                     &
      dstCompLabel           = OCN_LABEL,                     &
      compSetServicesRoutine = CPL_SetServices,               &
      rc                     = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! Conector OCN -> MED : SST/correntes -> bulk do mediador (sempre ativo)
    !--------------------------------------------------------------------------
    call NUOPC_DriverAddComp(driver,                          &
      srcCompLabel           = OCN_LABEL,                     &
      dstCompLabel           = MED_LABEL,                     &
      compSetServicesRoutine = CPL_SetServices,               &
      rc                     = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! Conector OCN -> MPAS : SST lag-1 -> pool sfc_input do MPAS
    !   Ativo apenas quando use_mpas_atm=.true.
    !   Campos: So_t, So_s, So_u, So_v
    !   Regrid: MOM6_grid (estruturado) -> MPAS_mesh (Voronoi), bilinear
    !--------------------------------------------------------------------------
    if (use_mpas_atm) then
      call NUOPC_DriverAddComp(driver,                        &
        srcCompLabel           = OCN_LABEL,                   &
        dstCompLabel           = MPAS_LABEL,                  &
        compSetServicesRoutine = CPL_SetServices,             &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    deallocate(petList)

    if (use_mpas_atm) then
      call ESMF_LogWrite( &
        'ESM: modo MPAS-ATM ? 4 conectores (MPAS->MED, MED->OCN, OCN->MED, OCN->MPAS)', &
        ESMF_LOGMSG_INFO)
    else
      call ESMF_LogWrite( &
        'ESM: modo DATM ? 3 conectores (DATM->MED, MED->OCN, OCN->MED)', &
        ESMF_LOGMSG_INFO)
    end if
  end subroutine SetModelServices

  !============================================================================
  ! SetRunSequence
  !
  !  Sequencia por passo de acoplamento (3h = 10800 s):
  !
  !  1. OCN -> MPAS : SST do passo anterior -> sfc_input (lag 1 passo)
  !  2. MPAS run    : fisica + dinâmica com SST atualizado
  !  3. DATM run    : le JRA55 (fallback)
  !  4. MPAS -> MED : campos atmosfericos do MPAS
  !  5. DATM -> MED : campos atmosfericos JRA55 (fallback)
  !  6. OCN  -> MED : SST para bulk formula no mediador
  !  7. MED  run    : calcula fluxos bulk NCAR
  !  8. MED  -> OCN : fluxos calculados -> MOM6
  !  9. OCN  run    : integra MOM6 com fluxos recebidos
  !
  !  Nota sobre a ordem OCN->MPAS antes de MPAS run:
  !    O SST transferido e do passo t-1 (lag de 1 intervalo de acoplamento),
  !    que e o comportamento padrao em sistemas acoplados atmosfera-oceano.
  !    Isso e equivalente ao "ocean lag" do CESM e de outros sistemas NUOPC.
  !    O primeiro passo usa o SST da condicao inicial do MPAS (namelist).
  !============================================================================
  subroutine SetRunSequence(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    type(NUOPC_FreeFormat) :: runSeqFF

    rc = ESMF_SUCCESS

    if (use_mpas_atm) then
      !------------------------------------------------------------------------
      ! Modo MPAS-ATM (use_mpas_atm = .true.)
      !
      !  1. OCN -> MPAS : SST lag t-1 -> pool sfc_input (antes do run MPAS)
      !  2. MPAS        : dinâmica + fisica com SST atualizado
      !  3. MPAS -> MED : campos atmosfericos da malha Voronoi
      !  4. OCN  -> MED : SST -> mediador para bulk
      !  5. MED         : calcula fluxos bulk NCAR (usa campos MPAS)
      !  6. MED  -> OCN : fluxos -> MOM6
      !  7. OCN         : integra MOM6
      !------------------------------------------------------------------------
      runSeqFF = NUOPC_FreeFormatCreate(stringList=(/ &
        "@10800            ",  &
        "  OCN -> MPAS     ",  &
        "  MPAS            ",  &
        "  MPAS -> MED     ",  &
        "  OCN  -> MED     ",  &
        "  MED             ",  &
        "  MED  -> OCN     ",  &
        "  OCN             ",  &
        "@                 " /), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_LogWrite('ESM: RunSequence modo MPAS-ATM', ESMF_LOGMSG_INFO)

    else
      !------------------------------------------------------------------------
      ! Modo DATM/JRA55 (use_mpas_atm = .false.)
      !
      !  1. DATM        : le JRA55 e exporta campos brutos
      !  2. DATM -> MED : campos JRA55 -> mediador
      !  3. OCN  -> MED : SST -> mediador para bulk
      !  4. MED         : calcula fluxos bulk NCAR (usa campos DATM)
      !  5. MED  -> OCN : fluxos -> MOM6
      !  6. OCN         : integra MOM6
      !------------------------------------------------------------------------
      runSeqFF = NUOPC_FreeFormatCreate(stringList=(/ &
        "@10800            ",  &
        "  DATM            ",  &
        "  DATM -> MED     ",  &
        "  OCN  -> MED     ",  &
        "  MED             ",  &
        "  MED  -> OCN     ",  &
        "  OCN             ",  &
        "@                 " /), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_LogWrite('ESM: RunSequence modo DATM/JRA55', ESMF_LOGMSG_INFO)

    end if

    call NUOPC_DriverIngestRunSequence(driver, runSeqFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_FreeFormatDestroy(runSeqFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine SetRunSequence

end module ESM
