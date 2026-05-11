!==============================================================================!
! esm.F90 - ESM Driver Component (MPAS + DATM + MOM6+SIS2)                     !
!==============================================================================!
! Arquitetura:                                                                 !
!   ATM (MPAS ou DATM) -> MED -> OCN (MOM6 real ou DOCN)                      !
!                                                                              !
! Conectores (até 5 ativos por configuração):                                  !
!   MPAS -> MED   : campos atmosfericos do MPAS (primário)                    !
!   DATM -> MED   : campos atmosfericos do JRA55 (fallback)                   !
!   OCN  -> MED   : SST/So_s/So_u/So_v para bulk formula (lag 1 step)         !
!   MED  -> OCN   : fluxos calculados (Foxx_*, Faxa_*, etc.)                  !
!   OCN  -> MPAS  : SST/correntes diretamente para pool sfc_input do MPAS     !
!                                                                              !
! Chaves de controle (runtime via nuopc.input &nuopc_esm):                    !
!                                                                              !
!   use_mpas_atm  : .true.  -> MPAS-ATM como fonte atmosferica primaria       !
!                   .false. -> DATM/JRA55 como fonte (desenvolvimento/fallback)!
!                                                                              !
!   use_mom6_ocn  : .true.  -> MOM6+SIS2 integrado (mom_cap.F90)              !
!                   .false. -> DOCN dados prescritos (ocn_comp_NUOPC.F90)      !
!                                                                              !
! Combinações válidas e casos de uso:                                          !
!                                                                              !
!   use_mpas_atm=F  use_mom6_ocn=F  ? DATM + DOCN                             !
!     Desenvolvimento puro: sem modelos activos, tudo prescritivo.             !
!     Valida o mediador e os conectores de forma isolada.                      !
!                                                                              !
!   use_mpas_atm=F  use_mom6_ocn=T  ? DATM + MOM6  (padrão atual)             !
!     MOM6 integrado com forçamento atmosférico JRA55. Produção sem MPAS.      !
!                                                                              !
!   use_mpas_atm=T  use_mom6_ocn=F  ? MPAS + DOCN                             !
!     MPAS com SST prescrita. Útil para ajuste de parametrização atmosférica   !
!     sem o custo do MOM6 e sem precisar de restart oceânico.                  !
!                                                                              !
!   use_mpas_atm=T  use_mom6_ocn=T  ? MPAS + MOM6  (acoplamento completo)     !
!     Configuração alvo do GT Acoplamento MONAN. Ambos os modelos ativos,      !
!     troca bidirecional de fluxos e SST via MED.                              !
!                                                                              !
! Nota sobre OCN->MPAS:                                                        !
!   O MPAS precisa de SST antes de cada timestep de fisica de superficie.      !
!   O conector OCN->MPAS usa regrid bilinear MOM6_grid -> MPAS_mesh.           !
!   Os campos conectados sao: So_t, So_s, So_u, So_v.                          !
!   O MPAS_cap anunciou esses campos com "cannot provide" + "share",           !
!   portanto o conector reusa a geometria já realizada pelo MOM_cap.            !
!==============================================================================!
module ESM

  use ESMF,  only: ESMF_GridComp,ESMF_Clock,ESMF_Time,ESMF_GridComp
  use ESMF,  only: ESMF_LogFoundError,ESMF_FAILURE,ESMF_GridCompGet
  use ESMF,  only: ESMF_ClockGet,ESMF_TimeGet
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

  ! Caps dos componentes atmosfericos
  use MPAS_cap_mod,   only: MPAS_SetServices => SetServices
  use DATM_cap_mod,   only: DATM_SetServices => SetServices

  ! Cap do mediador
  use MED_cap_mod,    only: MED_SetServices  => SetServices

  ! Caps oceanicos
  ! NOTA B-60: MOM_cap_mod NAO pode aparecer em "use" no topo do modulo ESM
  ! quando use_mom6_ocn=F. O FMS2 chama MOM_infra_init() dentro do
  ! MOM_cap SetServices e grid_init() tenta abrir INPUT/ocean_mosaic.nc
  ! mesmo sem NUOPC_DriverAddComp ter sido chamado para o MOM_cap.
  ! Solucao: use encapsulado em subrotina separada chamada apenas quando
  ! cfg_use_mom6_ocn=.true. (RegisterMOM6 abaixo).
  ! DOCN: sempre seguro de importar (nao tem FMS/mosaic dependency).
  !PK use MOM_cap_mod,    only: MOM6_SetServices => SetServices   ! MOM6+SIS2 integrado
  use ocn_comp_NUOPC, only: DOCN_SetServices => SetServices   ! DOCN dados prescritos

  ! Configuracao runtime: le use_mpas_atm e use_mom6_ocn do nuopc.input
  use mpas_cap_config_mod, only:  cfg_dt_coupling,   &
                                  config_read,       &
                                  cfg_use_mpas_atm,  &
                                  cfg_use_mom6_ocn

  implicit none
  private
  public :: SetServices

  !----------------------------------------------------------------------------
  ! Rótulos dos componentes NUOPC
  !----------------------------------------------------------------------------
  character(len=*), parameter :: MPAS_LABEL = "MPAS"
  character(len=*), parameter :: DATM_LABEL = "DATM"
  character(len=*), parameter :: MED_LABEL  = "MED"
  character(len=*), parameter :: OCN_LABEL  = "OCN"

  !============================================================================
  ! CHAVE 1 ? Fonte atmosférica       (runtime ? lida de nuopc.input)
  ! CHAVE 2 ? Componente oceânico     (runtime ? lida de nuopc.input)
  !
  ! Ambas as flags são lidas em SetModelServices via config_read() antes de
  ! qualquer registro de componente. Os defaults abaixo nunca entram em vigor
  ! durante execucao normal: config_read() sobrescreve com o valor do arquivo.
  ! Sao mantidos apenas como documentacao dos defaults do modulo config.
  !
  !   use_mpas_atm = .false.  (default: DATM + MOM6, comportamento v1.0)
  !   use_mom6_ocn = .true.
  !
  ! Para mudar entre configuracoes: edite nuopc.input, nao este arquivo.
  !============================================================================
  ! (sem parametros aqui ? variaveis lidas de mpas_cap_config_mod)

  !----------------------------------------------------------------------------
  ! Diagnóstico de SST do mediador
  !   med_diag_sst  : .true. habilita escrita NetCDF de SST
  !   med_diag_dir  : diretório de saída (criado automaticamente pelo MED_cap)
  !   med_diag_freq : frequência em passos de acoplamento
  !                   (1 = todo passo 3h; 4 = a cada 12h; 8 = a cada 24h)
  !----------------------------------------------------------------------------
  logical,          parameter :: med_diag_sst  = .true.
  character(len=*), parameter :: med_diag_dir  = "./diag"
  integer,          parameter :: med_diag_freq = 1

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
    character(len=4)  :: diag_freq_str   ! para converter med_diag_freq -> string
    integer           :: cfg_rc   ! rc de config_read (nao fatal)

    rc = ESMF_SUCCESS

    !--------------------------------------------------------------------------
    ! Leitura do nuopc.input ? deve ocorrer ANTES de qualquer registro de
    ! componente pois use_mpas_atm e use_mom6_ocn controlam quais caps
    ! e conectores sao registrados nesta subrotina.
    !
    ! config_read() e seguro antes de MPI_Init (B-29): usa apenas I/O Fortran.
    ! Se o arquivo nao for encontrado, os defaults do modulo sao mantidos
    ! (use_mpas_atm=F, use_mom6_ocn=T) e a execucao continua normalmente.
    !--------------------------------------------------------------------------
    call config_read(cfg_rc)
    if (cfg_rc == 2) then
      ! Erro fatal de leitura (arquivo corrompido ou permissao negada)
      call ESMF_LogWrite( &
        'ESM: config_read falhou (rc=2) -- abortando SetModelServices', &
        ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; return
    end if

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
    ! ATM: MPAS (primario) - registrado apenas quando cfg_use_mpas_atm=.true.
    !--------------------------------------------------------------------------
    if (cfg_use_mpas_atm) then
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
    ! ATM: DATM (fallback) - registrado apenas quando cfg_use_mpas_atm=.false.
    !--------------------------------------------------------------------------
    if (.not. cfg_use_mpas_atm) then
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
    ! MED: Mediador (sempre ativo)
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
    ! O MED_cap le este atributo em InitializeAdvertise e InitializeRealize.
    ! para decidir qual fonte atmosferica usar em MediatorAdvance.
    call NUOPC_CompAttributeAdd(medComp, attrList=(/"use_mpas_atm"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    if (cfg_use_mpas_atm) then
      call NUOPC_CompAttributeSet(medComp, name="use_mpas_atm", value="true",  rc=rc)
    else
      call NUOPC_CompAttributeSet(medComp, name="use_mpas_atm", value="false", rc=rc)
    end if
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! -----------------------------------------------------------------------
    ! NOVO: Atributos de diagnostico de SST do mediador
    !
    ! med_diag_sst  : "true" habilita escrita de NetCDF de SST a cada passo
    ! med_diag_dir  : diretorio de saida (criado automaticamente pelo MED_cap)
    ! med_diag_freq : frequencia em passos de acoplamento
    !                 Exemplos para passo de acoplamento de 3h:
    !                   "1"  = todo passo (3h)
    !                   "4"  = a cada 12h
    !                   "8"  = a cada 24h
    !
    ! Arquivos gerados:
    !   <med_diag_dir>/med_sst_atm_NNNNN.nc     ? SST na grade ATM 640x320
    !   <med_diag_dir>/med_sst_voronoi_NNNNN.nc ? SST na malha Voronoi MPAS
    !                                              (apenas se use_mpas_atm=.true.)
    !
    ! Para DESABILITAR: mudar med_diag_sst = .false. no topo deste arquivo.
    ! -----------------------------------------------------------------------
    call NUOPC_CompAttributeAdd(medComp, &
      attrList=(/"med_diag_sst ", "med_diag_dir ", "med_diag_freq"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (med_diag_sst) then
      call NUOPC_CompAttributeSet(medComp, name="med_diag_sst", value="true", rc=rc)
    else
      call NUOPC_CompAttributeSet(medComp, name="med_diag_sst", value="false", rc=rc)
    end if
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompAttributeSet(medComp, name="med_diag_dir", &
      value=trim(med_diag_dir), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    write(diag_freq_str, '(I0)') med_diag_freq
    call NUOPC_CompAttributeSet(medComp, name="med_diag_freq", &
      value=trim(diag_freq_str), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (med_diag_sst) then
      call ESMF_LogWrite('ESM: diagnostico SST MED habilitado -> '// &
        trim(med_diag_dir), ESMF_LOGMSG_INFO)
    else
      call ESMF_LogWrite('ESM: diagnostico SST MED desabilitado', ESMF_LOGMSG_INFO)
    end if

    !--------------------------------------------------------------------------
    ! OCN: MOM6+SIS2 real ou DOCN dados prescritos ? chave cfg_use_mom6_ocn
    !--------------------------------------------------------------------------
    if (cfg_use_mom6_ocn) then

      !------------------------------------------------------------------------
      ! MOM6+SIS2 integrado (mom_cap.F90)
      ! B-60: registrado via subrotina wrapper RegisterMOM6 para evitar que
      ! o "use MOM_cap_mod" no topo do modulo dispare MOM_infra_init/grid_init
      ! quando use_mom6_ocn=F e os arquivos INPUT/ocean_mosaic.nc nao existem.
      !------------------------------------------------------------------------
      call NUOPC_DriverAddComp(driver,                          &
        compLabel              = OCN_LABEL,                     &
        compSetServicesRoutine = MOM6_SetServices_wrapper,      &
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

      ! Propaga stopTime para o MOM6 (necessario para o restart automatico)
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

      call ESMF_LogWrite('ESM: OCN = MOM6+SIS2 integrado (mom_cap.F90)', &
        ESMF_LOGMSG_INFO)

    else

      !------------------------------------------------------------------------
      ! DOCN dados prescritos (ocn_comp_NUOPC.F90)
      !------------------------------------------------------------------------
      call NUOPC_DriverAddComp(driver,                          &
        compLabel              = OCN_LABEL,                     &
        compSetServicesRoutine = DOCN_SetServices,              &
        petList                = petList,                       &
        comp                   = ocnComp,                       &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call NUOPC_CompAttributeSet(ocnComp, name="Verbosity",  value="high",  rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_CompAttributeSet(ocnComp, name="DumpFields", value="false", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_LogWrite('ESM: OCN = DOCN dados prescritos (ocn_comp_NUOPC.F90)', &
        ESMF_LOGMSG_INFO)

    end if

    !--------------------------------------------------------------------------
    ! Conectores condicionais: ATM -> MED ou DATM -> MED
    !--------------------------------------------------------------------------
    if (cfg_use_mpas_atm) then
      ! MPAS -> MED: campos atmosfericos da malha Voronoi
      call NUOPC_DriverAddComp(driver,                          &
        srcCompLabel           = MPAS_LABEL,                    &
        dstCompLabel           = MED_LABEL,                     &
        compSetServicesRoutine = CPL_SetServices,               &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    else
      ! DATM -> MED: campos JRA55 (grade regular 640x320)
      call NUOPC_DriverAddComp(driver,                          &
        srcCompLabel           = DATM_LABEL,                    &
        dstCompLabel           = MED_LABEL,                     &
        compSetServicesRoutine = CPL_SetServices,               &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    !--------------------------------------------------------------------------
    ! Conector MED -> OCN : fluxos bulk -> MOM6 (sempre ativo)
    !--------------------------------------------------------------------------
    call NUOPC_DriverAddComp(driver,                            &
      srcCompLabel           = MED_LABEL,                       &
      dstCompLabel           = OCN_LABEL,                       &
      compSetServicesRoutine = CPL_SetServices,                 &
      rc                     = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! Conector OCN -> MED: SST/correntes -> bulk do mediador (sempre ativo) (sempre registrado)
    !--------------------------------------------------------------------------
    call NUOPC_DriverAddComp(driver,                            &
      srcCompLabel           = OCN_LABEL,                       &
      dstCompLabel           = MED_LABEL,                       &
      compSetServicesRoutine = CPL_SetServices,                 &
      rc                     = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! Conector OCN -> MPAS : SST lag-1 -> pool sfc_input do MPAS
    ! Ativo apenas quando cfg_use_mpas_atm=.true.
    !   Campos: So_t, So_s, So_u, So_v
    !   Regrid: MOM6_grid (estruturado) -> MPAS_mesh (Voronoi), bilinear
    !--------------------------------------------------------------------------
    if (cfg_use_mpas_atm) then
      call NUOPC_DriverAddComp(driver,                          &
        srcCompLabel           = OCN_LABEL,                     &
        dstCompLabel           = MPAS_LABEL,                    &
        compSetServicesRoutine = CPL_SetServices,               &
        rc                     = rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    deallocate(petList)

    if (cfg_use_mpas_atm) then
      call ESMF_LogWrite( &
        'ESM: modo MPAS-ATM ? 4 conectores (MPAS->MED, MED->OCN, OCN->MED, OCN->MPAS)', &
        ESMF_LOGMSG_INFO)
    else
      call ESMF_LogWrite( &
        'ESM: modo DATM ? 3 conectores (DATM->MED, MED->OCN, OCN->MED)', &
        ESMF_LOGMSG_INFO)
    end if
    !--------------------------------------------------------------------------
    ! Log do modo de operacao ativo
    !--------------------------------------------------------------------------
    call LogModoAtivo()

  end subroutine SetModelServices

  !============================================================================
  ! SetRunSequence
  !
  !  Sequencia por passo de acoplamento (3h = 10800 s):
  ! A RunSequence é idêntica para MOM6 e DOCN ? o rótulo "OCN" é o mesmo em
  ! ambos os casos. A diferença está apenas em quais campos são de fato
  ! transferidos pelos conectores (determinado pelos anúncios dos caps).
  !
  ! Sequência por passo de acoplamento (dt_coupling = 10800 s = 3 h):
  !
  !  Modo MPAS + OCN (qualquer):
  !    OCN ? MPAS  : SST lag t-1 ? sfc_input (antes do run MPAS)
  !    MPAS        : dinâmica + física com SST atualizado
  !    MPAS ? MED  : campos atmosféricos da malha Voronoi
  !    OCN  ? MED  : SST/correntes ? mediador para bulk
  !    MED         : calcula fluxos bulk NCAR
  !    MED  ? OCN  : fluxos ? OCN (integrado pelo MOM6 ou ignorados pelo DOCN)
  !    OCN         : integra MOM6 / avança DOCN (lê próximo snapshot)
  !
  !  Modo DATM + OCN (qualquer):
  !    DATM        : lê JRA55 e exporta campos brutos
  !    DATM ? MED  : campos JRA55 ? mediador
  !    OCN  ? MED  : SST/correntes ? mediador para bulk
  !    MED         : calcula fluxos bulk NCAR
  !    MED  ? OCN  : fluxos ? OCN
  !    OCN         : integra MOM6 / avança DOCN
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
    character(len=18)       :: dt_str   ! "@NNNN" do RunSequence
    character(len=6)        :: dt_val   ! valor numerico

    rc = ESMF_SUCCESS

    ! B-60: usar dt_coupling do nuopc.input em vez de hardcode
    write(dt_val, '(I0)') cfg_dt_coupling
    dt_str = '@' // trim(adjustl(dt_val))
    ! Pad para 18 chars (NUOPC_FreeFormat exige strings do mesmo comprimento)
    dt_str = dt_str // repeat(' ', 18 - len_trim(dt_str))

    if (cfg_use_mpas_atm) then
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
        dt_str               ,  &
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

      call ESMF_LogWrite('ESM: RunSequence modo MPAS-ATM + OCN ('// &
        merge('MOM6','DOCN', cfg_use_mom6_ocn)//')', ESMF_LOGMSG_INFO)

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
        dt_str               ,  &
        "  DATM            ",  &
        "  DATM -> MED     ",  &
        "  OCN  -> MED     ",  &
        "  MED             ",  &
        "  MED  -> OCN     ",  &
        "  OCN             ",  &
        "@                 " /), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return

      call ESMF_LogWrite('ESM: RunSequence modo DATM + OCN ('// &
        merge('MOM6','DOCN', cfg_use_mom6_ocn)//')', ESMF_LOGMSG_INFO)

    end if

    call NUOPC_DriverIngestRunSequence(driver, runSeqFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_FreeFormatDestroy(runSeqFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine SetRunSequence

  !============================================================================
  ! MOM6_SetServices_wrapper ? B-60: encapsula o use MOM_cap_mod
  !
  ! O FMS2 inicializa o framework (MOM_infra_init, grid_init) dentro do
  ! SetServices do MOM_cap, abrindo INPUT/ocean_mosaic.nc imediatamente.
  ! Se esse use estivesse no topo do modulo ESM, o linker resolveria os
  ! simbolos e o FMS seria ativado mesmo com use_mom6_ocn=F.
  !
  ! Ao colocar o "use MOM_cap_mod" aqui ? dentro de uma subrotina chamada
  ! APENAS quando cfg_use_mom6_ocn=.true. ? o FMS so e inicializado quando
  ! o MOM6 esta realmente sendo usado.
  !
  ! Assinatura NUOPC obrigatoria: sem intent nos argumentos ESMF.
  !============================================================================
  subroutine MOM6_SetServices_wrapper(gcomp, rc)
    use MOM_cap_mod, only: MOM6_SS => SetServices
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    call MOM6_SS(gcomp, rc)
  end subroutine MOM6_SetServices_wrapper

  !============================================================================
  ! LogModoAtivo ? escreve no log ESMF a configuração completa ativa
  !============================================================================
  subroutine LogModoAtivo()

    character(len=128) :: msg

    ! Linha de resumo compacta para grep rápido nos logs
    write(msg, '(A,L1,A,L1,A)') &
      'ESM: use_mpas_atm=', cfg_use_mpas_atm, &
      '  use_mom6_ocn=',    cfg_use_mom6_ocn, &
      ' ??????????????????????????'
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! Descrição expandida
    if (cfg_use_mpas_atm .and. cfg_use_mom6_ocn) then
      call ESMF_LogWrite( &
        'ESM: CONFIGURACAO = MPAS-ATM + MOM6+SIS2 (acoplamento completo)', &
        ESMF_LOGMSG_INFO)
      call ESMF_LogWrite( &
        'ESM: conectores ativos: MPAS->MED, OCN->MED, MED->OCN, OCN->MPAS', &
        ESMF_LOGMSG_INFO)

    else if (cfg_use_mpas_atm .and. (.not. cfg_use_mom6_ocn)) then
      call ESMF_LogWrite( &
        'ESM: CONFIGURACAO = MPAS-ATM + DOCN (SST prescrita)', &
        ESMF_LOGMSG_INFO)
      call ESMF_LogWrite( &
        'ESM: conectores ativos: MPAS->MED, OCN->MED, MED->OCN, OCN->MPAS', &
        ESMF_LOGMSG_INFO)
      call ESMF_LogWrite( &
        'ESM: NOTA: MED->OCN registrado mas DOCN stub suprime importacoes (B-39)', &
        ESMF_LOGMSG_INFO)

    else if ((.not. cfg_use_mpas_atm) .and. cfg_use_mom6_ocn) then
      call ESMF_LogWrite( &
        'ESM: CONFIGURACAO = DATM/JRA55 + MOM6+SIS2', &
        ESMF_LOGMSG_INFO)
      call ESMF_LogWrite( &
        'ESM: conectores ativos: DATM->MED, OCN->MED, MED->OCN', &
        ESMF_LOGMSG_INFO)

    else
      call ESMF_LogWrite( &
        'ESM: CONFIGURACAO = DATM/JRA55 + DOCN (tudo prescritivo)', &
        ESMF_LOGMSG_INFO)
      call ESMF_LogWrite( &
        'ESM: conectores ativos: DATM->MED, OCN->MED, MED->OCN', &
        ESMF_LOGMSG_INFO)
      call ESMF_LogWrite( &
        'ESM: NOTA: MED->OCN registrado mas DOCN stub suprime importacoes (B-39)', &
        ESMF_LOGMSG_INFO)
    end if

  end subroutine LogModoAtivo

end module ESM
