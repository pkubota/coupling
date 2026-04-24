!==============================================================================!
! MPAS_cap.F90 - MPAS Atmosphere NUOPC Component (Primario)                    !
!==============================================================================!
! Exporta campos atmosfericos da malha Voronoi MPAS para o mediador.           !
!                                                                              !
! Campos exportados:                                                           !
!   Sa_u10m_mpas   - vento zonal 10m       [m s-1]                            !
!   Sa_v10m_mpas   - vento meridional 10m  [m s-1]                            !
!   Sa_tbot_mpas   - temperatura do ar 2m  [K]                                !
!   Sa_shum_mpas   - umidade especifica    [kg kg-1]                           !
!   Sa_pslv_mpas   - pressao niv. mar      [Pa]                               !
!   Faxa_swdn_mpas - rad. solar descendente [W m-2]                           !
!   Faxa_lwdn_mpas - rad. longa desc.      [W m-2]                            !
!   Faxa_rain_mpas - precipitacao liquida  [kg m-2 s-1]                       !
!   Faxa_snow_mpas - precipitacao solida   [kg m-2 s-1]                       !
!                                                                              !
! Malha ESMF construida a partir de mpas_scrip.nc:                            !
!   - nos:           grid_corner_lon / grid_corner_lat  (nCells x nCorners)   !
!   - conectividade: um elemento por celula, vertices locais                  !
!   - elementCoords: grid_center_lon / grid_center_lat  (nCells)              !
!     OBRIGATORIO para ESMF_MeshCreateDual (regrid MESHLOC_ELEMENT)           !
!==============================================================================!

module MPAS_cap_mod

  use ESMF
  use NUOPC, only: NUOPC_CompDerive, NUOPC_CompSpecialize, NUOPC_CompSetEntryPoint
  use NUOPC, only: NUOPC_CompFilterPhaseMap, NUOPC_Advertise, NUOPC_Realize
  use NUOPC, only: NUOPC_SetTimestamp, NUOPC_CompAttributeSet
  use NUOPC_Model, only: model_routine_SS => SetServices
  use NUOPC_Model, only: model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model, only: model_label_Advance => label_Advance
  use NUOPC_Model, only: NUOPC_ModelGet

  use netcdf

  implicit none
  private
  public :: SetServices

  real(ESMF_KIND_R8), parameter :: PI      = 3.14159265358979323846_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: RAD2DEG = 180.0_ESMF_KIND_R8 / PI

  type :: MPAS_InternalState
    type(ESMF_Mesh) :: mesh
    integer :: nCells
    real(ESMF_KIND_R8), allocatable :: u10(:), v10(:), t2m(:), q2(:)
    real(ESMF_KIND_R8), allocatable :: psl(:), swdn(:), lwdn(:)
    real(ESMF_KIND_R8), allocatable :: prra(:), prsn(:)
    logical :: initialized = .false.
  end type MPAS_InternalState

  type :: MPAS_InternalStateWrapper
    type(MPAS_InternalState), pointer :: wrap => null()
  end type MPAS_InternalStateWrapper

  integer, parameter :: n_export = 9
  character(len=32), parameter :: export_names(n_export) = [ &
    "Sa_u10m_mpas  ", "Sa_v10m_mpas  ", "Sa_tbot_mpas  ", "Sa_shum_mpas  ", "Sa_pslv_mpas  ", &
    "Faxa_swdn_mpas", "Faxa_lwdn_mpas", "Faxa_rain_mpas", "Faxa_snow_mpas" ]

contains

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

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
      specRoutine=InitializeDataComplete, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
  end subroutine SetServices

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    rc = ESMF_SUCCESS
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
  end subroutine InitializeP0

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    integer :: n
    rc = ESMF_SUCCESS
    do n = 1, n_export
      call NUOPC_Advertise(exportState, StandardName=trim(export_names(n)), &
        TransferOfferGeomObject="will provide", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do
    call ESMF_LogWrite('MPAS: InitializeAdvertise concluido', ESMF_LOGMSG_INFO)
  end subroutine InitializeAdvertise

  !============================================================================
  ! InitializeRealize
  !
  ! Le mpas_scrip.nc e constroi ESMF_Mesh com:
  !   - nos:           grid_corner_lon/lat  -> nodeCoords
  !   - elementos:     uma celula por elemento, nCorners vertices locais
  !   - elementCoords: grid_center_lon/lat  -> centros dos elementos
  !                    OBRIGATORIO: sem este argumento, ESMF_MeshCreateDual
  !                    falha com "Creation of a dual mesh requires element
  !                    coordinates" ao fazer regrid de MESHLOC_ELEMENT.
  !============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_Mesh)  :: mesh
    type(ESMF_Field) :: field
    type(MPAS_InternalStateWrapper) :: iswrap
    type(MPAS_InternalState), pointer :: is

    ! Variaveis NetCDF
    integer :: ncid, ncrc, dimid
    integer :: varid_clon, varid_clat   ! grid_corner_lon / grid_corner_lat
    integer :: varid_elon, varid_elat   ! grid_center_lon / grid_center_lat

    ! Dimensoes
    integer :: nCells, nCorners, nNodes

    ! Arrays lidos do SCRIP
    ! CORRECAO 1: SCRIP armazena (nCorners, nCells) - primeira dimensao e nCorners.
    ! A declaracao anterior (nCells, nCorners) invertia os indices e corromperia
    ! a malha: celula i receberia os vertices de outra celula apos o reshape.
    real(ESMF_KIND_R8), allocatable :: corner_lon(:,:)   ! (nCorners, nCells)
    real(ESMF_KIND_R8), allocatable :: corner_lat(:,:)   ! (nCorners, nCells)
    real(ESMF_KIND_R8), allocatable :: center_lon(:)     ! (nCells)
    real(ESMF_KIND_R8), allocatable :: center_lat(:)     ! (nCells)

    ! Arrays para ESMF_MeshCreate
    integer,            allocatable :: nodeIds(:)
    real(ESMF_KIND_R8), allocatable :: nodeCoords(:)     ! (2*nNodes)
    integer,            allocatable :: nodeOwners(:)
    integer,            allocatable :: elemIds(:)
    integer,            allocatable :: elemTypes(:)
    integer,            allocatable :: elemConn(:)       ! (nCells*nCorners)
    real(ESMF_KIND_R8), allocatable :: elemCoords(:)     ! (2*nCells)

    integer :: i, j, n
    character(len=256) :: msg
    logical :: file_exists

    rc = ESMF_SUCCESS

    !--------------------------------------------------------------------------
    ! 1. Verificar arquivo
    !--------------------------------------------------------------------------
    inquire(file="mpas_scrip.nc", exist=file_exists)
    if (.not. file_exists) then
      call ESMF_LogWrite('MPAS: arquivo mpas_scrip.nc nao encontrado!', &
        ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    call ESMF_LogWrite('MPAS: lendo malha de mpas_scrip.nc', ESMF_LOGMSG_INFO)

    !--------------------------------------------------------------------------
    ! 2. Abrir e ler dimensoes
    !--------------------------------------------------------------------------
    ncrc = nf90_open("mpas_scrip.nc", NF90_NOWRITE, ncid)
    if (ncrc /= NF90_NOERR) then
      write(msg,*) 'MPAS: falha ao abrir mpas_scrip.nc, erro=', ncrc
      call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    ncrc = nf90_inq_dimid(ncid, "grid_size", dimid)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao ler dimensao grid_size', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if
    ncrc = nf90_inquire_dimension(ncid, dimid, len=nCells)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao obter tamanho grid_size', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    ncrc = nf90_inq_dimid(ncid, "grid_corners", dimid)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao ler dimensao grid_corners', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if
    ncrc = nf90_inquire_dimension(ncid, dimid, len=nCorners)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao obter tamanho grid_corners', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    write(msg,'(A,I8,A,I2,A)') 'MPAS: lendo ', nCells, &
      ' celulas com ', nCorners, ' vertices'
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    !--------------------------------------------------------------------------
    ! 3. Obter varids
    !--------------------------------------------------------------------------
    ncrc = nf90_inq_varid(ncid, "grid_corner_lon", varid_clon)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao encontrar grid_corner_lon', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if
    ncrc = nf90_inq_varid(ncid, "grid_corner_lat", varid_clat)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao encontrar grid_corner_lat', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if
    ncrc = nf90_inq_varid(ncid, "grid_center_lon", varid_elon)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao encontrar grid_center_lon', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if
    ncrc = nf90_inq_varid(ncid, "grid_center_lat", varid_elat)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao encontrar grid_center_lat', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    !--------------------------------------------------------------------------
    ! 4. Alocar e ler arrays
    !--------------------------------------------------------------------------
    ! CORRECAO 1 (cont.): alocacao com ordem correta (nCorners, nCells)
    allocate(corner_lon(nCorners, nCells), corner_lat(nCorners, nCells))
    allocate(center_lon(nCells),           center_lat(nCells))

    ! CORRECAO 2: leitura em bloco unico em vez de loop celula-a-celula.
    ! O loop anterior gerava O(nCells) chamadas nf90_get_var (~40 mil para
    ! x1.40962), o que e extremamente lento para I/O paralelo.
    call ESMF_LogWrite('MPAS: lendo coordenadas dos vertices (leitura bulk)...', &
      ESMF_LOGMSG_INFO)
    ncrc = nf90_get_var(ncid, varid_clon, corner_lon)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao ler grid_corner_lon', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if
    ncrc = nf90_get_var(ncid, varid_clat, corner_lat)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao ler grid_corner_lat', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    ! Centros: array 1D, lido de uma vez
    ncrc = nf90_get_var(ncid, varid_elon, center_lon)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao ler grid_center_lon', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if
    ncrc = nf90_get_var(ncid, varid_elat, center_lat)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: erro ao ler grid_center_lat', ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE; return
    end if

    ncrc = nf90_close(ncid)
    if (ncrc /= NF90_NOERR) then
      call ESMF_LogWrite('MPAS: aviso - erro ao fechar arquivo', ESMF_LOGMSG_WARNING)
    end if

    call ESMF_LogWrite('MPAS: dados lidos com sucesso', ESMF_LOGMSG_INFO)

    !--------------------------------------------------------------------------
    ! 5. Converter radianos -> graus se necessario
    !    Criterio: max|corner_lon| <= 2*pi + 0.1 => radianos
    !--------------------------------------------------------------------------
    if (maxval(abs(corner_lon)) <= 2.0_ESMF_KIND_R8*PI + 0.1_ESMF_KIND_R8) then
      call ESMF_LogWrite('MPAS: convertendo coordenadas de radianos para graus', &
        ESMF_LOGMSG_INFO)
      corner_lon = corner_lon * RAD2DEG
      corner_lat = corner_lat * RAD2DEG
      center_lon = center_lon * RAD2DEG
      center_lat = center_lat * RAD2DEG
    else
      call ESMF_LogWrite('MPAS: coordenadas ja estao em graus', ESMF_LOGMSG_INFO)
    end if

    !--------------------------------------------------------------------------
    ! 6. Montar arrays de nos (cantos de cada celula - sem compartilhamento)
    !--------------------------------------------------------------------------
    nNodes = nCells * nCorners
    allocate(nodeIds(nNodes))
    allocate(nodeCoords(2*nNodes))
    allocate(nodeOwners(nNodes))

    do i = 1, nCells
      do j = 1, nCorners
        n = (i-1)*nCorners + j
        nodeIds(n)        = n
        ! CORRECAO 1 (cont.): corner_lon/lat agora e (nCorners,nCells) -> indice (j,i)
        nodeCoords(2*n-1) = corner_lon(j, i)
        nodeCoords(2*n  ) = corner_lat(j, i)
        nodeOwners(n)     = 0
      end do
    end do

    !--------------------------------------------------------------------------
    ! 7. Montar conectividade de elementos (cada celula tem nCorners vertices)
    !--------------------------------------------------------------------------
    allocate(elemIds(nCells))
    allocate(elemTypes(nCells))
    allocate(elemConn(nCells * nCorners))

    do i = 1, nCells
      elemIds(i)   = i
      elemTypes(i) = nCorners
      do j = 1, nCorners
        elemConn((i-1)*nCorners + j) = (i-1)*nCorners + j
      end do
    end do

    !--------------------------------------------------------------------------
    ! 8. Montar coordenadas de elemento (centros das celulas)
    !    OBRIGATORIO: sem este argumento, ESMF_MeshCreateDual falha com:
    !      "Creation of a dual mesh requires element coordinates"
    !    O dual mesh e criado internamente pelo ESMF ao fazer regrid de
    !    campos em MESHLOC_ELEMENT (caso do conector MPAS->MED).
    !--------------------------------------------------------------------------
    allocate(elemCoords(2*nCells))
    do i = 1, nCells
      elemCoords(2*i-1) = center_lon(i)
      elemCoords(2*i  ) = center_lat(i)
    end do

    !--------------------------------------------------------------------------
    ! 9. Criar malha ESMF
    !--------------------------------------------------------------------------
    call ESMF_LogWrite('MPAS: criando malha ESMF...', ESMF_LOGMSG_INFO)

    mesh = ESMF_MeshCreate(                   &
      parametricDim  = 2,                     &
      spatialDim     = 2,                     &
      coordSys       = ESMF_COORDSYS_SPH_DEG, &
      nodeIds        = nodeIds,               &
      nodeCoords     = nodeCoords,            &
      nodeOwners     = nodeOwners,            &
      elementIds     = elemIds,               &
      elementTypes   = elemTypes,             &
      elementConn    = elemConn,              &
      elementCoords  = elemCoords,            &
      rc             = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="MPAS: ESMF_MeshCreate falhou", &
      line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite('MPAS: malha criada com sucesso', ESMF_LOGMSG_INFO)

    !--------------------------------------------------------------------------
    ! 10. Realizar campos de exportacao na malha MPAS
    !--------------------------------------------------------------------------
    do n = 1, n_export
      field = ESMF_FieldCreate(mesh=mesh, typekind=ESMF_TYPEKIND_R8, &
        meshloc=ESMF_MESHLOC_ELEMENT, name=trim(export_names(n)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, &
        msg="MPAS: falha ao criar campo "//trim(export_names(n)), &
        line=__LINE__, file=__FILE__)) return
      call NUOPC_Realize(exportState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end do

    !--------------------------------------------------------------------------
    ! 11. Inicializar estado interno
    !--------------------------------------------------------------------------
    allocate(iswrap%wrap)
    is => iswrap%wrap
    is%mesh        = mesh
    is%nCells      = nCells
    is%initialized = .false.

    allocate(is%u10(nCells), is%v10(nCells), is%t2m(nCells), is%q2(nCells))
    allocate(is%psl(nCells), is%swdn(nCells), is%lwdn(nCells))
    allocate(is%prra(nCells), is%prsn(nCells))

    call ESMF_GridCompSetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    !--------------------------------------------------------------------------
    ! 12. Limpeza
    !--------------------------------------------------------------------------
    deallocate(corner_lon, corner_lat, center_lon, center_lat)
    deallocate(nodeIds, nodeCoords, nodeOwners)
    deallocate(elemIds, elemTypes, elemConn, elemCoords)

    call ESMF_LogWrite('MPAS: InitializeRealize concluido', ESMF_LOGMSG_INFO)
  end subroutine InitializeRealize

  subroutine InitializeDataComplete(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    type(ESMF_State)  :: exportState
    type(ESMF_Field)  :: field
    real(ESMF_KIND_R8), pointer :: fptr(:)
    integer :: fieldCount, i
    character(len=64), allocatable :: fieldNameList(:)

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
        fptr = 0.0_ESMF_KIND_R8
      end do
      deallocate(fieldNameList)
    end if
    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataProgress", &
      value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", &
      value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_LogWrite('MPAS: InitializeDataComplete SATISFIED', ESMF_LOGMSG_INFO)
  end subroutine InitializeDataComplete

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    type(ESMF_State)         :: exportState
    type(ESMF_Clock)         :: clock
    type(ESMF_Time)          :: currTime, nextTime
    type(ESMF_TimeInterval)  :: dt
    type(MPAS_InternalStateWrapper) :: iswrap
    type(MPAS_InternalState), pointer :: is
    type(ESMF_Field) :: field
    character(len=256) :: msg
    integer :: year, month, day, hour, minute, second
    integer :: fieldCount, k
    character(len=64), allocatable :: fieldNameList(:)

    rc = ESMF_SUCCESS
    call ESMF_GridCompGetInternalState(gcomp, iswrap, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    is => iswrap%wrap

    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=dt, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    nextTime = currTime + dt

    call ESMF_TimeGet(currTime, yy=year, mm=month, dd=day, &
      h=hour, m=minute, s=second, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    write(msg,'(A,I4,5(A,I2.2))') 'MPAS: avancando para ', year, '-', &
      month, '-', day, ' ', hour, ':', minute, ':', second
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    !=========================================================================
    ! AQUI VOCE INTEGRA COM O MODELO MPAS
    ! Por enquanto, valores constantes para teste:
    !=========================================================================
    is%u10  = 5.0_ESMF_KIND_R8
    is%v10  = 5.0_ESMF_KIND_R8
    is%t2m  = 290.0_ESMF_KIND_R8
    is%q2   = 0.01_ESMF_KIND_R8
    is%psl  = 101325.0_ESMF_KIND_R8
    is%swdn = 0.0_ESMF_KIND_R8
    is%lwdn = 300.0_ESMF_KIND_R8
    is%prra = 0.0_ESMF_KIND_R8
    is%prsn = 0.0_ESMF_KIND_R8

    call PutField1D(exportState, "Sa_u10m_mpas",   is%u10,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField1D(exportState, "Sa_v10m_mpas",   is%v10,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField1D(exportState, "Sa_tbot_mpas",   is%t2m,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField1D(exportState, "Sa_shum_mpas",   is%q2,   rc); if (rc/=ESMF_SUCCESS) return
    call PutField1D(exportState, "Sa_pslv_mpas",   is%psl,  rc); if (rc/=ESMF_SUCCESS) return
    call PutField1D(exportState, "Faxa_swdn_mpas", is%swdn, rc); if (rc/=ESMF_SUCCESS) return
    call PutField1D(exportState, "Faxa_lwdn_mpas", is%lwdn, rc); if (rc/=ESMF_SUCCESS) return
    call PutField1D(exportState, "Faxa_rain_mpas", is%prra, rc); if (rc/=ESMF_SUCCESS) return
    call PutField1D(exportState, "Faxa_snow_mpas", is%prsn, rc); if (rc/=ESMF_SUCCESS) return

    ! Atualizar timestamps dos campos exportados
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

    call ESMF_LogWrite('MPAS: ModelAdvance concluido', ESMF_LOGMSG_INFO)
  end subroutine ModelAdvance

  subroutine PutField1D(state, name, array, rc)
    type(ESMF_State),   intent(inout) :: state
    character(len=*),   intent(in)    :: name
    real(ESMF_KIND_R8), intent(in)    :: array(:)
    integer,            intent(out)   :: rc
    type(ESMF_Field) :: field
    real(ESMF_KIND_R8), pointer :: fptr(:)
    rc = ESMF_SUCCESS
    call ESMF_StateGet(state, itemName=trim(name), field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="PutField1D: "//trim(name), &
      line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return
    fptr = array
  end subroutine PutField1D

end module MPAS_cap_mod
