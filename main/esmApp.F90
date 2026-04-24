!==============================================================================!
! esmApp.F90 - Programa principal do acoplamento MPAS + DATM + MOM6+SIS2       !
!==============================================================================!
! Arquitetura:                                                                 !
!   MPAS (primario) + DATM (fallback) -> MED (bulk NCAR) -> MOM6+SIS2          !
!                                                                              !
! O MED implementa fallback: usa MPAS quando disponivel, senao usa DATM        !
!==============================================================================!
program esmApp
  !use ESMF
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
  implicit none

  type(ESMF_GridComp)     :: esmComp
  type(ESMF_State)        :: importState, exportState
  type(ESMF_Clock)        :: clock
  type(ESMF_TimeInterval) :: timeStep
  type(ESMF_Time)         :: startTime, stopTime
  type(ESMF_VM)           :: vm
  integer                 :: localPet, petCount
  integer                 :: rc, urc
  logical                 :: isStopTime
  integer                 :: yy_start, mm_start, dd_start, hh_start, mn_start
  integer                 :: yy_stop,  mm_stop,  dd_stop,  hh_stop,  mn_stop
  integer                 :: dt_coupling_s

  !----------------------------------------------------------------------------
  ! Parametros do experimento
  !----------------------------------------------------------------------------
  yy_start      = 2016 ;  mm_start = 1  ;  dd_start = 1 ; hh_start = 1 ; mn_start = 30
  yy_stop       = 2016 ;  mm_stop  = 1  ;  dd_stop  = 5 ; hh_stop  = 1 ; mn_stop  = 30
  dt_coupling_s = 10800    ! passo do driver = 3h

  !----------------------------------------------------------------------------
  ! 1. Inicializa ESMF
  !----------------------------------------------------------------------------
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, &
                       defaultCalKind=ESMF_CALKIND_GREGORIAN, &
                       defaultLogFileName='esm_acoplamento.log', &
                       rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !PK REMOVER ou comentar esta linha:
  !PK call NUOPC_FieldDictionarySetAutoAdd(.true.)

  !PK Substituir por:
  call NUOPC_FieldDictionarySetAutoAdd(.false.)

  call ESMF_VMGetGlobal(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (localPet == 0) then
    write(*,'(A)') '========================================================='
    write(*,'(A)') ' Acoplamento MPAS (primario) + DATM (fallback) + MOM6+SIS2'
    write(*,'(A)') '========================================================='
    write(*,'(A,I4)')                      '  PETs        = ', petCount
    write(*,'(A,I4.4,"-",I2.2,"-",I2.2)') '  Inicio: ', yy_start, mm_start, dd_start
    write(*,'(A,I4.4,"-",I2.2,"-",I2.2)') '  Fim:    ', yy_stop,  mm_stop,  dd_stop
    write(*,'(A,I7," s")')                 '  dt_cpl  = ', dt_coupling_s
    write(*,'(A)') '  Arquitetura: MPAS + DATM -> MED -> OCN'
    write(*,'(A)') '  Fallback: MED usa MPAS se disponivel, senao DATM'
    write(*,'(A)') '========================================================='
  end if

  !----------------------------------------------------------------------------
  ! 2. Configura relogio
  !----------------------------------------------------------------------------
  call ESMF_TimeIntervalSet(timeStep, s=dt_coupling_s, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_TimeSet(startTime, yy=yy_start, mm=mm_start, dd=dd_start, &
    h=hh_start, m=mn_start, s=0, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  call ESMF_TimeSet(stopTime, yy=yy_stop, mm=mm_stop, dd=dd_stop, &
    h=hh_stop, m=mn_stop, s=0, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  clock = ESMF_ClockCreate(timeStep=timeStep, startTime=startTime, &
                           stopTime=stopTime, &
                           name="ESM_ApplicationClock", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (localPet == 0) write(*,'(A)') '[OK] Relogio criado'

  !----------------------------------------------------------------------------
  ! 3. Cria driver ESM
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
  ! 4. Estados import/export do driver raiz (vazios)
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
  ! 5. Inicializacao (MPAS, DATM, MED, OCN + conectores)
  !----------------------------------------------------------------------------
  if (localPet == 0) write(*,'(A)') '[  ] Inicializando componentes...'

  call ESMF_GridCompInitialize(esmComp, importState=importState, &
    exportState=exportState, clock=clock, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (localPet == 0) write(*,'(A)') '[OK] Driver inicializado'

  !----------------------------------------------------------------------------
  ! 6. Loop temporal
  ! CORRECAO: rc de ClockIsStopTime agora e verificado explicitamente.
  !----------------------------------------------------------------------------
  if (localPet == 0) write(*,'(A)') '[  ] Iniciando integracao...'

  do
    isStopTime = ESMF_ClockIsStopTime(clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (isStopTime) exit

    call ESMF_GridCompRun(esmComp,                &
                          importState=importState, &
                          exportState=exportState, &
                          clock=clock,             &
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
  ! 7. Finalizacao
  !----------------------------------------------------------------------------
  if (localPet == 0) write(*,'(A)') '[  ] Finalizando...'

  call ESMF_GridCompFinalize(esmComp, importState=importState, &
    exportState=exportState, clock=clock, userRc=urc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  if (localPet == 0) write(*,'(A)') '[OK] Finalizado com sucesso'

  !----------------------------------------------------------------------------
  ! 8. Limpeza
  !----------------------------------------------------------------------------
  call ESMF_ClockDestroy(clock, rc=rc)
  call ESMF_StateDestroy(importState, rc=rc)
  call ESMF_StateDestroy(exportState, rc=rc)
  call ESMF_GridCompDestroy(esmComp, rc=rc)
  call ESMF_Finalize(rc=rc)

end program esmApp
