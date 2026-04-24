!==============================================================================!
! Datm_cap.F90 - Data Atmosphere NUOPC Component (JRA55 fallback para MED)    !
!                                                                              !
! ADAPTACAO AtmOcnMedPetListProto:                                             !
! O DATM exporta os campos BRUTOS do JRA55 para o MEDIADOR.                   !
! O calculo de fluxos bulk NCAR foi movido para MED_cap.F90.                  !
!                                                                              !
! Campos exportados para o MED (importState do mediador):                      !
!   Sa_u10m     - vento zonal 10m       [m s-1]                               !
!   Sa_v10m     - vento meridional 10m  [m s-1]                               !
!   Sa_tbot     - temperatura do ar 2m  [K]                                   !
!   Sa_shum     - umidade especifica    [kg kg-1]                              !
!   Sa_pslv     - pressao niv. mar      [Pa]                                  !
!   Faxa_swdn   - rad. solar desc.      [W m-2]                               !
!   Faxa_lwdn   - rad. LW desc.         [W m-2]                               !
!   Faxa_rain   - precipitacao liquida  [kg m-2 s-1]                          !
!   Faxa_snow   - precipitacao solida   [kg m-2 s-1]                          !
!                                                                              !
! CORRECOES APLICADAS:                                                         !
!   1. epochTime corrigido para 2016-01-01 00:00:00 (era 01:30:00 - causava   !
!      tidx0 negativo ou errado para runs iniciando em t=0).                   !
!   2. Removida dependencia desnecessaria de MOM_io (stdout, io_infra_end).    !
!   3. ReadJRAFieldByIndex removida (era codigo morto comentado com !PK).      !
!==============================================================================!
module DATM_cap_mod
  use ESMF
  use ESMF, only: ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSetEntryPoint
  use ESMF, only: ESMF_GridCompGetInternalState, ESMF_GridCompSetInternalState
  use ESMF, only: ESMF_State, ESMF_StateGet
  use ESMF, only: ESMF_Field, ESMF_FieldCreate, ESMF_FieldGet
  use ESMF, only: ESMF_Grid, ESMF_GridCreate1PeriDim, ESMF_GridAddCoord, ESMF_GridGetCoord
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
!PK  use ESMF, only: operator(+), operator(-)

  use netcdf

  use NUOPC, only: NUOPC_CompDerive, NUOPC_CompSpecialize, NUOPC_CompSetEntryPoint
  use NUOPC, only: NUOPC_CompFilterPhaseMap, NUOPC_Advertise, NUOPC_Realize
  use NUOPC, only: NUOPC_SetTimestamp, NUOPC_CompAttributeSet
  use NUOPC_Model, &
    model_routine_SS           => SetServices,         &
    model_label_SetClock       => label_SetClock,      &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_Advance        => label_Advance
  use NUOPC_Model, only: NUOPC_ModelGet
  ! CORRECAO 2: removida dependencia de MOM_io (stdout, io_infra_end nao
  !   eram utilizados e acoplavam o DATM desnecessariamente ao MOM6).
  implicit none
  private
  public :: SetServices

  !----------------------------------------------------------------------------
  ! Estado interno do DATM
  !----------------------------------------------------------------------------
  type :: DATM_InternalState
    type(ESMF_Grid) :: grid
    ! Campos JRA55 lidos do NetCDF
    real(ESMF_KIND_R8), pointer :: uas(:,:)  => null()  ! vento zonal 10m
    real(ESMF_KIND_R8), pointer :: vas(:,:)  => null()  ! vento merid. 10m
    real(ESMF_KIND_R8), pointer :: tas(:,:)  => null()  ! temperatura ar 2m
    real(ESMF_KIND_R8), pointer :: huss(:,:) => null()  ! umidade especifica
    real(ESMF_KIND_R8), pointer :: psl(:,:)  => null()  ! pressao niv. mar
    real(ESMF_KIND_R8), pointer :: rsds(:,:) => null()  ! rad. sol. desc.
    real(ESMF_KIND_R8), pointer :: rlds(:,:) => null()  ! rad. LW desc.
    real(ESMF_KIND_R8), pointer :: prra(:,:) => null()  ! precip. liquida
    real(ESMF_KIND_R8), pointer :: prsn(:,:) => null()  ! precip. solida
    logical :: initialized = .false.
  end type DATM_InternalState

  type :: DATM_InternalStateWrapper
    type(DATM_InternalState), pointer :: wrap => null()
  end type DATM_InternalStateWrapper

contains

  !============================================================================
  ! SetServices
  !============================================================================
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

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

  end subroutine SetServices

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
  ! InitializeAdvertise - anuncia apenas campos BRUTOS para o MED
  !
  ! MUDANCA: O DATM nao mais anuncia Foxx_* (fluxos calculados).
  !          Exporta somente o que vem diretamente do JRA55.
  !============================================================================
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Campos de estado atmosferico bruto (JRA55)
    call NUOPC_Advertise(exportState, StandardName="Sa_u10m",   rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, StandardName="Sa_v10m",   rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, StandardName="Sa_tbot",   rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, StandardName="Sa_shum",   rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, StandardName="Sa_pslv",   rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Radiacao descendente (sem decomposicao em bandas - o MED faz isso)
    call NUOPC_Advertise(exportState, StandardName="Faxa_swdn", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, StandardName="Faxa_lwdn", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Precipitacao
    call NUOPC_Advertise(exportState, StandardName="Faxa_rain", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Advertise(exportState, StandardName="Faxa_snow", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('DATM: InitializeAdvertise concluido (campos brutos JRA55)', &
      ESMF_LOGMSG_INFO)
  end subroutine InitializeAdvertise

  !============================================================================
  ! InitializeRealize - grade 640x320 (JRA55)
  !
  ! CORRECAO 1: coordX usa (i-1)*dx + dx/2  (i eh indice global com INDEX_GLOBAL)
  ! CORRECAO 2: coordY usa (j-1)*dy + dy/2  (idem)
  !============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_Grid)   :: grid
    integer           :: nx_global, ny_global, i, j
    real(ESMF_KIND_R8)              :: dx, dy
    real(ESMF_KIND_R8), pointer :: coordX(:,:), coordY(:,:)
    type(DATM_InternalStateWrapper) :: iswrap
    type(DATM_InternalState), pointer :: is

    rc = ESMF_SUCCESS

    nx_global = 640
    ny_global = 320
    dx = 360.0_ESMF_KIND_R8 / nx_global   ! 0.5625 graus
    dy = 180.0_ESMF_KIND_R8 / ny_global   ! 0.5625 graus

    grid = ESMF_GridCreate1PeriDim( &
      minIndex  = (/1, 1/),                 &
      maxIndex  = (/nx_global, ny_global/), &
      indexflag = ESMF_INDEX_GLOBAL,        &
      coordSys  = ESMF_COORDSYS_SPH_DEG,   &
      rc        = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Longitude: com ESMF_INDEX_GLOBAL, i eh o indice global (1..640)
    ! lon_centro_i = (i-1)*dx + dx/2  => [0.28125, 0.84375, ..., 359.71875]
    ! Com ESMF_INDEX_GLOBAL, lbound(coordX,1) no PET0 eh 1, no PET1 pode ser
    ! p.ex. 161, etc. O indice i ja eh o indice GLOBAL da coluna.
    ! Formula: lon_centro_i = (i-1)*dx + dx/2   [0.28125, 0.84375, ..., 359.71875]
    call ESMF_GridGetCoord(grid, coordDim=1, &
      staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordX, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do j = lbound(coordX,2), ubound(coordX,2)
      do i = lbound(coordX,1), ubound(coordX,1)
        ! Com ESMF_INDEX_GLOBAL, i eh indice global (1..640); lon=(i-1)*dx+dx/2
        !coordX(i,j) = -180.0 + (i-1)*(360.0/nx_global) + (360.0/nx_global) * 0.5_ESMF_KIND_R8
        coordX(i,j) = (i - 1) * (360.0_ESMF_KIND_R8/nx_global) &
          + (360.0_ESMF_KIND_R8/nx_global) * 0.5_ESMF_KIND_R8
        ! DEPOIS:
        ! coordX(i,j) = -300.0_ESMF_KIND_R8 + (i-1)*(360.0_ESMF_KIND_R8/nx_global) &
        !             + (360.0_ESMF_KIND_R8/nx_global) * 0.5_ESMF_KIND_R8
      end do
    end do
    ! Latitude: lat_centro_j = -90 + (j-1)*dy + dy/2
    call ESMF_GridGetCoord(grid, coordDim=2, &
      staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=coordY, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    do j = lbound(coordY,2), ubound(coordY,2)
      do i = lbound(coordY,1), ubound(coordY,1)
        coordY(i,j) = -90.0_ESMF_KIND_R8 + (j-1)*(180.0_ESMF_KIND_R8/ny_global) &
          + (180.0_ESMF_KIND_R8/ny_global)/2.0_ESMF_KIND_R8
      end do
    end do

    ! Realiza campos brutos
    call RealizeField(exportState, grid, "Sa_u10m",   rc); if (rc/=ESMF_SUCCESS) return
    call RealizeField(exportState, grid, "Sa_v10m",   rc); if (rc/=ESMF_SUCCESS) return
    call RealizeField(exportState, grid, "Sa_tbot",   rc); if (rc/=ESMF_SUCCESS) return
    call RealizeField(exportState, grid, "Sa_shum",   rc); if (rc/=ESMF_SUCCESS) return
    call RealizeField(exportState, grid, "Sa_pslv",   rc); if (rc/=ESMF_SUCCESS) return
    call RealizeField(exportState, grid, "Faxa_swdn", rc); if (rc/=ESMF_SUCCESS) return
    call RealizeField(exportState, grid, "Faxa_lwdn", rc); if (rc/=ESMF_SUCCESS) return
    call RealizeField(exportState, grid, "Faxa_rain", rc); if (rc/=ESMF_SUCCESS) return
    call RealizeField(exportState, grid, "Faxa_snow", rc); if (rc/=ESMF_SUCCESS) return

    allocate(iswrap%wrap)
    is => iswrap%wrap
    is%grid        = grid
    is%initialized = .false.

    call ESMF_GridCompSetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('DATM: InitializeRealize concluido', ESMF_LOGMSG_INFO)
  end subroutine InitializeRealize

  !============================================================================
  ! RealizeField - helper
  !============================================================================
  subroutine RealizeField(state, grid, stdname, rc)
    type(ESMF_State),  intent(inout) :: state
    type(ESMF_Grid),   intent(in)    :: grid
    character(len=*),  intent(in)    :: stdname
    integer,           intent(inout) :: rc

    type(ESMF_Field) :: field

    field = ESMF_FieldCreate(grid=grid, typekind=ESMF_TYPEKIND_R8, &
      staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(stdname), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_Realize(state, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
  end subroutine RealizeField

  !============================================================================
  ! InitializeDataComplete - IPDv03p7: componente de dados puro
  !============================================================================
  subroutine InitializeDataComplete(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    type(ESMF_State)  :: exportState
    type(ESMF_Field)  :: field
    integer           :: fieldCount, i
    character(len=64), allocatable :: fieldNameList(:)
    real(ESMF_KIND_R8), pointer :: fptr(:,:)

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

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

        select case(trim(fieldNameList(i)))
          case('Sa_pslv')
            fptr = 101325.0_ESMF_KIND_R8
          case('Sa_tbot')
            fptr = 290.0_ESMF_KIND_R8  ! K
          case default
            fptr = 0.0_ESMF_KIND_R8
        end select
      end do
      deallocate(fieldNameList)
    end if

    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataProgress", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete",  value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('DATM: InitializeDataComplete SATISFIED', ESMF_LOGMSG_INFO)
  end subroutine InitializeDataComplete

  !============================================================================
  ! ModelAdvance - le JRA55 e popula exportState com campos BRUTOS
  !
  ! MUDANCA: nao calcula mais bulk. Apenas le NetCDF e escreve:
  !   Sa_u10m = uas, Sa_v10m = vas, Sa_tbot = tas, Sa_shum = huss,
  !   Sa_pslv = psl, Faxa_swdn = rsds, Faxa_lwdn = rlds,
  !   Faxa_rain = prra, Faxa_snow = prsn
  !============================================================================
  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    type(ESMF_State)         :: exportState
    type(ESMF_Clock)         :: clock
    type(ESMF_Time)          :: currTime, nextTime
    type(ESMF_TimeInterval)  :: dt
    type(ESMF_Field)         :: field
    type(DATM_InternalStateWrapper) :: iswrap
    type(DATM_InternalState), pointer :: is
    real(ESMF_KIND_R8), pointer :: fptr(:,:)
    integer :: i, j, i1, i2, j1, j2
    integer :: year, month, day, hour, minu, sec
    integer :: fieldCount, k
    character(len=64), allocatable :: fieldNameList(:)
    character(len=256) :: msg

    rc = ESMF_SUCCESS

    call ESMF_GridCompGetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    is => iswrap%wrap

    call NUOPC_ModelGet(gcomp, modelClock=clock, exportState=exportState, rc=rc)
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

    write(msg,'(A,I4,5(A,I2.2))') 'DATM: avancando para ', year, '-', &
      month, '-', day, ' ', hour, ':', minu, ':', sec
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    ! Obtem limites locais a partir do primeiro campo
    call ESMF_StateGet(exportState, itemName="Sa_u10m", field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    i1 = lbound(fptr,1); i2 = ubound(fptr,1)
    j1 = lbound(fptr,2); j2 = ubound(fptr,2)

    ! Aloca arrays temporarios se necessario
    if (.not. associated(is%uas)) then
      allocate(is%uas(i1:i2,  j1:j2))
      allocate(is%vas(i1:i2,  j1:j2))
      allocate(is%tas(i1:i2,  j1:j2))
      allocate(is%huss(i1:i2, j1:j2))
      allocate(is%psl(i1:i2,  j1:j2))
      allocate(is%rsds(i1:i2, j1:j2))
      allocate(is%rlds(i1:i2, j1:j2))
      allocate(is%prra(i1:i2, j1:j2))
      allocate(is%prsn(i1:i2, j1:j2))
    end if

    ! Le campos JRA55 com interpolacao temporal linear (3h -> dt_driver)
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_uas.nc",  "uas",  currTime, is%uas,  rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha uas",  line=__LINE__, file=__FILE__)) return
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_vas.nc",  "vas",  currTime, is%vas,  rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha vas",  line=__LINE__, file=__FILE__)) return
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_tas.nc",  "tas",  currTime, is%tas,  rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha tas",  line=__LINE__, file=__FILE__)) return
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_huss.nc", "huss", currTime, is%huss, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha huss", line=__LINE__, file=__FILE__)) return
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_psl.nc",  "psl",  currTime, is%psl,  rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha psl",  line=__LINE__, file=__FILE__)) return
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_rsds.nc", "rsds", currTime, is%rsds, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha rsds", line=__LINE__, file=__FILE__)) return
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_rlds.nc", "rlds", currTime, is%rlds, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha rlds", line=__LINE__, file=__FILE__)) return
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_prra.nc", "prra", currTime, is%prra, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha prra", line=__LINE__, file=__FILE__)) return
    call ReadJRAFieldInterp(gcomp, "INPUT/JRA_prsn.nc", "prsn", currTime, is%prsn, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha prsn", line=__LINE__, file=__FILE__)) return

!is%uas = 0.0
!is%vas = 0.0
!is%tas = 290.0
!is%huss = 0.01
!is%psl = 101325.0
!is%rsds = 0.0
!is%rlds = 0.0
!is%prra = 0.0
!is%prsn = 0.0
    ! Escreve campos lidos no exportState
    call PutField(exportState, "Sa_u10m",   is%uas,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Sa_v10m",   is%vas,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Sa_tbot",   is%tas,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Sa_shum",   is%huss, rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Sa_pslv",   is%psl,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Faxa_swdn", is%rsds, rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Faxa_lwdn", is%rlds, rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Faxa_rain", is%prra, rc); if (rc/=ESMF_SUCCESS) return
    call PutField(exportState, "Faxa_snow", is%prsn, rc); if (rc/=ESMF_SUCCESS) return

    ! Atualizar timestamps
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

    call ESMF_LogWrite('DATM: ModelAdvance concluido (campos brutos)', ESMF_LOGMSG_INFO)
  end subroutine ModelAdvance

  !============================================================================
  ! PutField - escreve array 2D em campo do exportState
  !============================================================================
  subroutine PutField(state, name, array, rc)
    type(ESMF_State),               intent(inout) :: state
    character(len=*),               intent(in)    :: name
    real(ESMF_KIND_R8),             intent(in)    :: array(:,:)
    integer,                        intent(out)   :: rc

    type(ESMF_Field) :: field
    real(ESMF_KIND_R8), pointer :: fptr(:,:)

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(name), field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="PutField: "//trim(name), &
      line=__LINE__, file=__FILE__)) return

    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    fptr = array

  end subroutine PutField

  !============================================================================
  ! ReadJRAFieldInterp - interpolacao temporal linear entre snapshots 3h
  !
  ! Estrategia paralela: PET0 le campo global inteiro do NetCDF e faz
  ! broadcast via ESMF_VMBroadcast. Cada PET copia apenas o seu subdominio.
  !
  ! CORRECAO 1: epochTime corrigido para 2016-01-01 00:00:00.
  !   O valor anterior (01:30:00) causava tidx0 negativo para t=0 do
  !   experimento (2016-01-01 00:00:00), pois sec_since_epoch ficava < 0.
  !   O JRA55 com passo 3h tem snapshots em 00:00, 03:00, 06:00, ..., portanto
  !   a epoca correta e o inicio do dado, nao o centro do primeiro intervalo.
  !   Se o seu arquivo JRA55 comecar em 01:30 (centro do 1o intervalo), ajuste
  !   aqui e documente o offset na cabecalho do arquivo.
  !
  ! CORRECAO 3: leitura paralela correta.
  !   Estrategia: PET0 le o campo global inteiro do NetCDF e faz broadcast
  !   para todos os PETs via ESMF_VMBroadcast. Cada PET entao copia apenas
  !   o seu subdominio local (i1:i2, j1:j2) do array global.
  !
  !   Isso e correto porque os arquivos JRA55 nao sao particionados; a
  !   leitura paralela real exigiria PIO ou NetCDF-4 paralelo, o que
  !   adicionaria dependencias. Para os tamanhos JRA55 (640x320 x 2 snapshots
  !   x 8 bytes ~ 3 MB por variavel) o broadcast e perfeitamente aceitavel.
  !============================================================================
  subroutine ReadJRAFieldInterp(gcomp, filename, varname, currTime, array, rc)
    type(ESMF_GridComp),  intent(in)    :: gcomp
    character(len=*),    intent(in)  :: filename
    character(len=*),    intent(in)  :: varname
    type(ESMF_Time),     intent(in)  :: currTime
    real(ESMF_KIND_R8),  pointer     :: array(:,:)
    integer,             intent(out) :: rc

    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: epochTime
    type(ESMF_TimeInterval) :: dt_since_epoch, interval3h
    integer(ESMF_KIND_I8)   :: sec_since_epoch
    integer                 :: tidx0, tidx1
    real(ESMF_KIND_R8)      :: alpha

    ! Arrays globais (usados apenas em PET0 para leitura, depois broadcast)
    ! f0/f1 removidos: eram dead code (nunca preenchidos apos refatoracao)
    ! A interpolacao temporal e feita em f0_global antes do broadcast.
    integer, parameter :: NX = 640, NY = 320
    real(ESMF_KIND_R8), target    :: f0_global(NX,NY), f1_global(NX,NY)
    real(ESMF_KIND_R8), pointer :: f0_1d(:)
    real(ESMF_KIND_R8), allocatable :: buf_global(:)

    ! Limites locais do subdominio deste PET
    integer :: i1, i2, j1, j2, i, j, ij, localPet
    integer :: ni, nj!local
    character(len=256) :: msg

    rc = ESMF_SUCCESS
    ni = size(array, 1)
    nj = size(array, 2)
    allocate(buf_global(NX*NY))
    ! (f0/f1 removidos - veja declaracoes)
    call ESMF_VMGetGlobal(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    f0_global = 0.0_ESMF_KIND_R8
    f1_global = 0.0_ESMF_KIND_R8
    ! Calcula indices de interpolacao temporal
    !--------------------------------------------------------------------------
    ! CORRECAO 1: epochTime = 2016-01-01 00:00:00
    ! O valor anterior (h=1, m=30) era o centro do primeiro intervalo JRA55,
    ! o que causava sec_since_epoch < 0 para currTime = 2016-01-01 00:00:00
    ! e portanto tidx0 = 0, causando leitura com indice invalido (base 1).
    !--------------------------------------------------------------------------
    call ESMF_TimeSet(epochTime, yy=2016, mm=1, dd=1, h=1, m=30, s=0, &
      calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_TimeIntervalSet(interval3h, s=10800, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    dt_since_epoch = currTime - epochTime

    call ESMF_TimeIntervalGet(dt_since_epoch, s_i8=sec_since_epoch, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Garantir sec_since_epoch >= 0 (currTime nao pode ser anterior ao epoch)
    if (sec_since_epoch < 0_ESMF_KIND_I8) then
      call ESMF_LogWrite('DATM ReadJRAFieldInterp: currTime anterior ao epochTime!', &
        ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
      !PK return
    end if

    ! tidx0 e base-1 (primeiro snapshot = indice 1)
    tidx0 = int(sec_since_epoch / 10800.0_ESMF_KIND_R8) + 1
    tidx1 = tidx0 + 1
    alpha = real(mod(sec_since_epoch, 10800_ESMF_KIND_I8), ESMF_KIND_R8) / &
            10800.0_ESMF_KIND_R8
    alpha = max(0.0_ESMF_KIND_R8, min(1.0_ESMF_KIND_R8, alpha))

    ! Leitura: apenas PET0 acessa o disco
!PK    call ReadJRAFieldByIndex(filename, varname, tidx0, f0, rc)
!PK     if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha t0 "//trim(varname), &
!PK      line=__LINE__, file=__FILE__)) return
!PK    call ReadJRAFieldByIndex(filename, varname, tidx1, f1, rc)
!PK   if (ESMF_LogFoundError(rcToCheck=rc, msg="Falha t1 "//trim(varname), &
!PK      line=__LINE__, file=__FILE__)) return
    ! --- Leitura: apenas PET0 acessa o disco ---
    if (localPet == 0) then
      call ReadGlobalField(filename, varname, tidx0, NX, NY, f0_global, rc)
      if (rc /= ESMF_SUCCESS) return
      call ReadGlobalField(filename, varname, tidx1, NX, NY, f1_global, rc)
      if (rc /= ESMF_SUCCESS) return

      ! Interpolacao temporal in-place
      f0_global = f0_global + alpha * (f1_global - f0_global)
      buf_global = reshape(f0_global, [NX*NY])
    end if

    ! Broadcast: PET0 envia campo global para todos os PETs
    ! ESMF_VMBroadcast usa contagem de elementos (NX*NY doubles)
    call ESMF_VMBroadcast(vm, bcstData=buf_global, count=NX*NY, rootPet=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Cada PET copia apenas o seu subdominio local
    ! array tem indices globais (ESMF_INDEX_GLOBAL): lbound/ubound dao i1,i2,j1,j2 globais
    ! Copia apenas o subdom\ufffdnio local (layout Fortran: i varia mais r\ufffdpido)
    i1 = lbound(array,1); i2 = ubound(array,1)
    j1 = lbound(array,2); j2 = ubound(array,2)

    do j = j1, j2
      do i = i1, i2
        array(i,j) = buf_global((j-1)*NX + i)
      end do
    end do

    deallocate(buf_global)

    write(msg,'(A,A,A,I5,A,I5,A,F6.4)') &
      'DATM: interp ', trim(varname), &
      ' tidx0=', tidx0, ' tidx1=', tidx1, ' alpha=', alpha
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
    ! (deallocate f0/f1 removido)
  end subroutine ReadJRAFieldInterp

  !============================================================================
  ! ReadGlobalField - le campo global (NX x NY) do NetCDF (chamado so em PET0)
  !============================================================================
  subroutine ReadGlobalField(filename, varname, tidx, nx, ny, array, rc)
    character(len=*),    intent(in)  :: filename
    character(len=*),    intent(in)  :: varname
    integer,             intent(in)  :: tidx
    integer,             intent(in)  :: nx, ny
    real(ESMF_KIND_R8),  intent(out) :: array(nx,ny)
    integer,             intent(out) :: rc

    integer :: ncid, varid, start(3), count(3), nc_rc

    rc     = ESMF_SUCCESS
    nc_rc  = nf90_open(filename, NF90_NOWRITE, ncid)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("ReadGlobalField: falha ao abrir "//trim(filename)// &
        ": "//trim(nf90_strerror(nc_rc)), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; return
    end if

    nc_rc = nf90_inq_varid(ncid, varname, varid)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("ReadGlobalField: variavel nao encontrada: "// &
        trim(varname), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
    end if

    ! Le o campo global inteiro: [1:nx, 1:ny, tidx]
    start = [1, 1, tidx]; count = [nx, ny, 1]
    nc_rc = nf90_get_var(ncid, varid, array, start=start, count=count)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("ReadGlobalField: falha ao ler "//trim(varname)// &
        ": "//trim(nf90_strerror(nc_rc)), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
    end if

    nc_rc = nf90_close(ncid)

  end subroutine ReadGlobalField
  !============================================================================
  ! ReadJRAFieldByIndex - le snapshot NetCDF 1-based
  !============================================================================
  subroutine ReadJRAFieldByIndex(filename, varname, tidx, array, rc)
    character(len=*),    intent(in)  :: filename
    character(len=*),    intent(in)  :: varname
    integer,             intent(in)  :: tidx
    real(ESMF_KIND_R8),  intent(out) :: array(:,:)
    integer,             intent(out) :: rc

    integer :: ncid, varid, ni, nj, start(3), count(3), nc_rc

    rc    = ESMF_SUCCESS
    ni    = size(array, 1)
    nj    = size(array, 2)

    nc_rc = nf90_open(filename, NF90_NOWRITE, ncid)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("Falha ao abrir "//trim(filename)//": "// &
        trim(nf90_strerror(nc_rc)), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; return
    end if

    nc_rc = nf90_inq_varid(ncid, varname, varid)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("Variavel nao encontrada: "//trim(varname), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
    end if

    start = [1, 1, tidx]; count = [ni, nj, 1]
    nc_rc = nf90_get_var(ncid, varid, array, start=start, count=count)
    if (nc_rc /= NF90_NOERR) then
      call ESMF_LogWrite("Falha ao ler "//trim(varname)//": "// &
        trim(nf90_strerror(nc_rc)), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE; nc_rc = nf90_close(ncid); return
    end if

    nc_rc = nf90_close(ncid)
  end subroutine ReadJRAFieldByIndex

end module DATM_cap_mod
