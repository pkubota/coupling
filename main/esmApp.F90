!==============================================================================!
! esmApp.F90 - Programa principal do acoplamento MPAS + DATM + MOM6+SIS2       !
!==============================================================================!
! Arquitetura:                                                                 !
!   MPAS (primario) + DATM (fallback) -> MED (bulk NCAR) -> MOM6+SIS2          !
!                                                                              !
! O MED implementa fallback: usa MPAS quando disponivel, senao usa DATM.       !
!                                                                              !
! Parametros de tempo (start_date, stop_date, dt_coupling) e chaves de        !
! configuracao (use_mpas_atm, use_mom6_ocn) sao lidos de nuopc.input via      !
! mpas_cap_config_mod antes de ESMF_Initialize.                                !
!                                                                              !
! Estrutura do arquivo nuopc.input:                                            !
!   &nuopc_esm                                                                 !
!     use_mpas_atm = .true.                                                    !
!     use_mom6_ocn = .true.                                                    !
!   /                                                                          !
!   &nuopc_driver                                                              !
!     start_date  = '2016-01-01'    ! formato YYYY-MM-DD                       !
!     stop_date   = '2016-01-05'    ! formato YYYY-MM-DD                       !
!     dt_coupling = 10800           ! passo do acoplador [s]                   !
!     dt_atm      = 60              ! passo interno do MPAS [s]                !
!     log_dir     = 'logs'                                                     !
!   /                                                                          !
!                                                                              !
! Sequencia de inicializacao:                                                  !
!   1. config_read()         -- le nuopc.input (puro Fortran, sem MPI)         !
!   2. ESMF_Initialize()     -- inicia ESMF + MPI                              !
!   3. config_print()        -- imprime config ativa (apenas PET 0)            !
!   4. ESMF_TimeSet/Clock    -- usa cfg_start_date / cfg_stop_date             !
!   5. ESM_SetServices       -- registra caps com cfg_use_mpas_atm/ocn         !
!==============================================================================!
program esmApp

  use ESMF, only: ESMF_GridComp, ESMF_State, ESMF_Clock, ESMF_TimeInterval
  use ESMF, only: ESMF_Time, ESMF_VM
  use ESMF, only: ESMF_Initialize, ESMF_VMGetGlobal, ESMF_VMGet
  use ESMF, only: ESMF_ClockCreate, ESMF_TimeIntervalSet, ESMF_TimeSet
  use ESMF, only: ESMF_ClockGet, ESMF_ClockAdvance, ESMF_ClockIsStopTime
  use ESMF, only: ESMF_ClockDestroy, ESMF_ClockPrint
  use ESMF, only: ESMF_StateCreate, ESMF_StateDestroy
  use ESMF, only: ESMF_GridCompCreate, ESMF_GridCompSetServices
  use ESMF, only: ESMF_GridCompInitialize, ESMF_GridCompRun, ESMF_GridCompFinalize
  use ESMF, only: ESMF_GridCompDestroy
  use ESMF, only: ESMF_LogFoundError, ESMF_Finalize
  use ESMF, only: ESMF_STATEINTENT_EXPORT, ESMF_STATEINTENT_IMPORT
  use ESMF, only: ESMF_END_ABORT, ESMF_LOGKIND_MULTI
  use ESMF, only: ESMF_LOGERR_PASSTHRU, ESMF_CALKIND_GREGORIAN

  use NUOPC, only: NUOPC_FieldDictionarySetAutoAdd

  use ESM, only: ESM_SetServices => SetServices

  ! Configuracao runtime: todas as datas e flags vem do nuopc.input
  use mpas_cap_config_mod, only: config_read,        &
                                  config_print,        &
                                  config_parse_date,   &
                                  cfg_start_date,      &
                                  cfg_stop_date,       &
                                  cfg_dt_coupling,     &
                                  cfg_log_dir,         &
                                  cfg_use_mpas_atm,    &
                                  cfg_use_mom6_ocn

  implicit none

  type(ESMF_GridComp)     :: esmComp
  type(ESMF_State)        :: importState, exportState
  type(ESMF_Clock)        :: clock
  type(ESMF_TimeInterval) :: timeStep
  type(ESMF_Time)         :: startTime, stopTime
  type(ESMF_VM)           :: vm

  integer :: localPet, petCount
  integer :: rc, urc, cfg_rc
  logical :: isStopTime

  ! Componentes das datas (extraidos de cfg_start_date / cfg_stop_date)
  integer :: yy_start, mm_start, dd_start
  integer :: yy_stop,  mm_stop,  dd_stop

  !----------------------------------------------------------------------------
  ! PASSO 1 ? Leitura de nuopc.input
  !
  ! DEVE ocorrer ANTES de ESMF_Initialize (que faz MPI_Init).
  ! config_read usa apenas I/O Fortran puro (sem MPI, sem ESMF).
  ! Se o arquivo nao for encontrado, os defaults do modulo sao usados e
  ! a execucao continua normalmente (rc=1 = nao fatal).
  !----------------------------------------------------------------------------
  call config_read(cfg_rc)
  if (cfg_rc == 2) then
    write(*,'(A)') '[ERRO] esmApp: config_read falhou (rc=2). Abortando.'
    stop 1
  end if

  !----------------------------------------------------------------------------
  ! PASSO 2 ? Inicializa ESMF
  !----------------------------------------------------------------------------
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI,          &
                       defaultCalKind=ESMF_CALKIND_GREGORIAN,   &
                       defaultLogFileName='esm_acoplamento.log', &
                       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call NUOPC_FieldDictionarySetAutoAdd(.false.)

  call ESMF_VMGetGlobal(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !----------------------------------------------------------------------------
  ! PASSO 3 ? Imprime configuracao ativa (apenas PET 0)
  !
  ! Feito apos ESMF_VMGet para ter localPet disponivel (B-29).
  !----------------------------------------------------------------------------
  if (localPet == 0) then
    call config_print()
    write(*,'(A)') '======================================================='
    write(*,'(A)') ' Acoplamento MPAS (primario) + DATM (fallback) + MOM6+SIS2'
    write(*,'(A)') '======================================================='
    write(*,'(A,I4)')   '  PETs          = ', petCount
    write(*,'(A,A)')    '  start_date    = ', trim(cfg_start_date)
    write(*,'(A,A)')    '  stop_date     = ', trim(cfg_stop_date)
    write(*,'(A,I7,A)') '  dt_coupling   = ', cfg_dt_coupling, ' s'
    write(*,'(A,L1)')   '  use_mpas_atm  = ', cfg_use_mpas_atm
    write(*,'(A,L1)')   '  use_mom6_ocn  = ', cfg_use_mom6_ocn
    write(*,'(A,A)')    '  log_dir       = ', trim(cfg_log_dir)
    if (cfg_use_mpas_atm .and. cfg_use_mom6_ocn) then
      write(*,'(A)') '  Modo: MPAS + MOM6 (acoplamento completo)'
    else if (cfg_use_mpas_atm) then
      write(*,'(A)') '  Modo: MPAS + DOCN (SST prescrita)'
    else if (cfg_use_mom6_ocn) then
      write(*,'(A)') '  Modo: DATM + MOM6'
    else
      write(*,'(A)') '  Modo: DATM + DOCN (tudo prescritivo)'
    end if
    write(*,'(A)') '======================================================='
  end if

  !----------------------------------------------------------------------------
  ! PASSO 4 ? Configura relogio ESMF com datas do nuopc.input
  !
  ! config_parse_date converte 'YYYY-MM-DD' -> (yy, mm, dd).
  ! As horas e minutos do experimento original (hh=1, mn=30) sao lidos
  ! do formato estendido 'YYYY-MM-DD HH:MM' se fornecido; caso contrario,
  ! assume meia-noite (hh=0, mn=0, s=0).
  !
  ! NOTA: o formato atual de cfg_start_date e 'YYYY-MM-DD' (10 chars).
  ! Para suporte a hora/minuto, expanda para 'YYYY-MM-DD HH:MM' (16 chars)
  ! no mpas_cap_config.F90 e ajuste config_parse_date adequadamente.
  ! Enquanto isso, hh e mn sao zero (meia-noite) automaticamente.
  !----------------------------------------------------------------------------
  call config_parse_date(cfg_start_date, yy_start, mm_start, dd_start, rc)
  if (rc /= 0) then
    if (localPet == 0) write(*,'(A,A)') &
      '[ERRO] esmApp: formato invalido em start_date: ', trim(cfg_start_date)
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end if

  call config_parse_date(cfg_stop_date, yy_stop, mm_stop, dd_stop, rc)
  if (rc /= 0) then
    if (localPet == 0) write(*,'(A,A)') &
      '[ERRO] esmApp: formato invalido em stop_date: ', trim(cfg_stop_date)
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end if

  call ESMF_TimeIntervalSet(timeStep, s=cfg_dt_coupling, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_TimeSet(startTime,                         &
    yy=yy_start, mm=mm_start, dd=dd_start,            &
    h=0, m=0, s=0,                                    &
    calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_TimeSet(stopTime,                          &
    yy=yy_stop, mm=mm_stop, dd=dd_stop,               &
    h=0, m=0, s=0,                                    &
    calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  clock = ESMF_ClockCreate(timeStep=timeStep,      &
                           startTime=startTime,    &
                           stopTime=stopTime,      &
                           name="ESM_ApplicationClock", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (localPet == 0) write(*,'(A)') '[OK] Relogio criado'

  !----------------------------------------------------------------------------
  ! PASSO 5 ? Cria driver ESM
  !----------------------------------------------------------------------------
  esmComp = ESMF_GridCompCreate(name="ESM_Driver", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_GridCompSetServices(esmComp, ESM_SetServices, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !----------------------------------------------------------------------------
  ! PASSO 6 ? Estados import/export do driver raiz (vazios)
  !----------------------------------------------------------------------------
  importState = ESMF_StateCreate(name="ESM_Import", &
    stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  exportState = ESMF_StateCreate(name="ESM_Export", &
    stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !----------------------------------------------------------------------------
  ! PASSO 7 ? Inicializacao (MPAS, DATM, MED, OCN + conectores)
  !----------------------------------------------------------------------------
  if (localPet == 0) write(*,'(A)') '[  ] Inicializando componentes...'

  call ESMF_GridCompInitialize(esmComp,               &
    importState=importState, exportState=exportState, &
    clock=clock, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (localPet == 0) write(*,'(A)') '[OK] Driver inicializado'

  !----------------------------------------------------------------------------
  ! PASSO 8 ? Loop temporal
  !----------------------------------------------------------------------------
  if (localPet == 0) write(*,'(A)') '[  ] Iniciando integracao...'

  do
    isStopTime = ESMF_ClockIsStopTime(clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (isStopTime) exit

    call ESMF_GridCompRun(esmComp,                  &
                          importState=importState,  &
                          exportState=exportState,  &
                          clock=clock,              &
                          userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_ClockAdvance(clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end do

  if (localPet == 0) write(*,'(A)') '[OK] Integracao concluida'

  !----------------------------------------------------------------------------
  ! PASSO 9 ? Finalizacao
  !----------------------------------------------------------------------------
  if (localPet == 0) write(*,'(A)') '[  ] Finalizando...'

  call ESMF_GridCompFinalize(esmComp,                &
    importState=importState, exportState=exportState,&
    clock=clock, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (localPet == 0) write(*,'(A)') '[OK] Finalizado com sucesso'

  !----------------------------------------------------------------------------
  ! PASSO 10 ? Limpeza
  !----------------------------------------------------------------------------
  call ESMF_ClockDestroy(clock, rc=rc)
  call ESMF_StateDestroy(importState, rc=rc)
  call ESMF_StateDestroy(exportState, rc=rc)
  call ESMF_GridCompDestroy(esmComp, rc=rc)
  call ESMF_Finalize(rc=rc)

end program esmApp
