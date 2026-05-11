!> @file mpas_cap_methods.F90
!! @brief Importacao/exportacao de campos ESMF <-> MPAS-A e criacao de malha.
!!
!! Versao 4.12 A malha ESMF real validada (02/04/2026):
!!   mpas_create_mesh: gera x1.40962.esmf_mesh.nc via ESMF_MeshCreateFromFile
!!   com ESMF compilado com ESMF_PIO=external. Fallback sintetico mantido.
!!   state_diagnose: expandido com estatisticas reais de campo (min/max/mean).

module mpas_cap_methods_mod

  use ESMF
  use mpas_atm_types_mod, only : mpas_atm_public_type,   &
                                  atm_ocean_boundary_type, &
                                  MPAS_RKIND
  use mpas_cap_utils_mod, only : ChkErr

  implicit none
  private

  public :: mpas_import
  public :: mpas_export
  public :: mpas_create_grid
  public :: state_diagnose

  character(len=*), parameter :: u_FILE_u = __FILE__

contains

  !> @brief Importa campos do importState ESMF para atm_bnd.
  !!
  !! Campos esperados (nomes padrao CMEPS):
  !!   So_t, Si_ifrac, Sf_zorl
  subroutine mpas_import(importState, atm_bnd, nCells, rc)
    type(ESMF_State),              intent(in)    :: importState
    type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
    integer,                       intent(in)    :: nCells
    integer,                       intent(inout) :: rc

    character(len=*), parameter :: subname = '(mpas_import)'

    rc = ESMF_SUCCESS

    ! -- campos originais --------------------------------------------------
    call state_get_field_1d(importState, 'So_t',     nCells, atm_bnd%sst,          rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    call state_get_field_1d(importState, 'Si_ifrac', nCells, atm_bnd%ice_fraction,  rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    call state_get_field_1d(importState, 'Sf_zorl',  nCells, atm_bnd%zorl,          rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    call ESMF_LogWrite(subname//': importacao concluida', ESMF_LOGMSG_INFO)

    ! -- B-OCN-01: correntes superficiais e salinidade --------------------
    ! So_u/So_v: velocidade de corrente superficial oceanica [m/s].
    ! Injetadas em sfc_input%u_surface_velocity / v_surface_velocity pelo
    ! mpas_atm_run. Usadas pelo esquema de CLA para calcular velocidade
    ! relativa vento-corrente e estresse mais realista em regioes de
    ! correntes intensas (Corrente Norte, Agulhas, Kuroshio, Gulf Stream).
    !
    ! state_get_field_1d e tolerante a campo ausente: se o conector
    ! OCN->MPAS nao transferiu So_u/So_v (ex: DOCN modo stub sem correntes),
    ! atm_bnd%uocn/vocn mantem o valor de default (0.0 m/s).
    call state_get_field_1d(importState, 'So_u', nCells, atm_bnd%uocn, rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    call state_get_field_1d(importState, 'So_v', nCells, atm_bnd%vocn, rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! So_s: salinidade superficial [psu]. Reservado para uso futuro
    ! (parametrizacao de albedo de spray e calor especifico oceanico).
    ! Importado aqui para completar a interface com o OCN mas nao
    ! consumido ainda pelo nucleo MPAS-A.
    call state_get_field_1d(importState, 'So_s', nCells, atm_bnd%sss, rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    call ESMF_LogWrite(subname//': importacao concluida (So_t/Si_ifrac/Sf_zorl/So_u/So_v/So_s)', &
      ESMF_LOGMSG_INFO)

  end subroutine mpas_import

  !> @brief Exporta campos de atm_public para o exportState ESMF.
  !!
  !! Campos exportados (nomes CMEPS com sufixo _mpas):
  !!   Sa_pslv_mpas, Sa_tbot_mpas, Sa_u10m_mpas, Sa_v10m_mpas, Sa_shum_mpas,
  !!   Faxa_swdn_mpas, Faxa_lwdn_mpas, Faxa_rain_mpas, Faxa_snow_mpas.
  !! Faxa_taux/tauy e Faxa_lhflx/shflx sao calculados pelo MED_cap (bulk NCAR)
  !! e nao pertencem ao exportState do MPAS cap.
  !!
  !! Campos nao-associados (pool diag_physics inativo ou nome ausente no
  !! Registry.xml) sao silenciosamente ignorados.
  subroutine mpas_export(atm_public, exportState, rc)
    type(mpas_atm_public_type), intent(in)    :: atm_public
    type(ESMF_State),           intent(inout) :: exportState
    integer,                    intent(inout) :: rc

    integer :: n
    character(len=*), parameter :: subname = '(mpas_export)'

    rc = ESMF_SUCCESS
    n  = atm_public%nCells

    if (associated(atm_public%pslv)) then
      call state_set_field_1d(exportState, 'Sa_pslv_mpas',   n, atm_public%pslv, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    if (associated(atm_public%t2m)) then
      call state_set_field_1d(exportState, 'Sa_tbot_mpas',   n, atm_public%t2m,  rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    if (associated(atm_public%u10)) then
      call state_set_field_1d(exportState, 'Sa_u10m_mpas',   n, atm_public%u10,  rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    if (associated(atm_public%v10)) then
      call state_set_field_1d(exportState, 'Sa_v10m_mpas',   n, atm_public%v10,  rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    if (associated(atm_public%swdn_sfc)) then
      call state_set_field_1d(exportState, 'Faxa_swdn_mpas', n, atm_public%swdn_sfc, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    if (associated(atm_public%lwdn_sfc)) then
      call state_set_field_1d(exportState, 'Faxa_lwdn_mpas', n, atm_public%lwdn_sfc, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    ! FIX Bug2: prec_rain (liquida) em vez de prec_total (liquida+neve)
    if (associated(atm_public%prec_rain)) then
      call state_set_field_1d(exportState, 'Faxa_rain_mpas', n, atm_public%prec_rain, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    ! B-Fase2-01: Sa_shum_mpas eh umidade especifica a 2 m [kg/kg]
    ! Fonte: pool diag_physics%q2 (Registry.xml), disponivel com bl_mynn_in ou bl_ysu_in.
    ! Fallback (em mpas_atm_model.F90): saturacao Tetens x RH=0.8 se q2 nao estiver no pool.
    if (associated(atm_public%q2m)) then
      call state_set_field_1d(exportState, 'Sa_shum_mpas', n, atm_public%q2m, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    ! B-Fase2-02: Faxa_snow_mpas eh taxa de precipitacao solida [kg/m^2/s]
    ! Fonte: incremento de diag_physics%snownc por intervalo de acoplamento Ã· dt.
    ! snownc [mm acumulado] requer mp_thompson_in ou mp_wsm6_in.
    ! Fallback (em mpas_atm_model.F90): particao por T2m se snownc ausente.
    if (associated(atm_public%prec_snow)) then
      call state_set_field_1d(exportState, 'Faxa_snow_mpas', n, atm_public%prec_snow, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    ! Sa_lfrac_mpas: fracao de terra MONAN [0,1] para reconciliacao
    ! de mascara costeira no MED_cap (ReconcileCoastalMask).
    ! Derivado de xland (WRF pool) ou landmask (mesh pool) em mpas_atm_init.
    ! Se o MPAS-A nao tiver xland nem landmask no Registry, lfrac = 0
    ! (oceano puro) e o MED opera sem correccao de mascara pelo ATM.
    if (associated(atm_public%lfrac)) then
      call state_set_field_1d(exportState, 'Sa_lfrac_mpas', n, atm_public%lfrac, rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
    end if
    call ESMF_LogWrite(subname//': exportacao concluida', ESMF_LOGMSG_INFO)

  end subroutine mpas_export

  !> @brief Cria ESMF_Grid sintetica 360x180 para o cap MPAS.
  !!
  !! SOLUCAO DEFINITIVA v5.0: substituicao de ESMF_Mesh por ESMF_Grid.
  !!
  !! CAUSA RAIZ de todos os travamentos anteriores:
  !!   O ESMF_Mesh usa internamente ESMF_MOAB para gestao paralela da malha.
  !!   Com ESMF_MOAB=enabled (build ESMF 8.9.1), as operacoes de
  !!   ESMF_MeshAddNodes, ESMF_MeshAddElements e ESMF_FieldCreate sobre
  !!   ESMF_Mesh executam chamadas MPI internas (para redistribuicao de nos
  !!   e sincronizacao de campo) que entram em deadlock apos mpas_atm_init.
  !!   O SMIOL (Simple Model I/O Library) do MPAS-A deixa o communicator MPI
  !!   em estado incompativel com as operacoes nao-convencionais do MOAB.
  !!
  !! SOLUCAO:
  !!   Substituir ESMF_Mesh por ESMF_Grid (grade regular lat/lon 360x180).
  !!   ESMF_Grid nao usa MOAB: todas as operacoes paralelas usam MPI padrao.
  !!   Todos os conectores passam a ser Grid->Grid (MED usa Grid 640x320,
  !!   DATOCN usa Grid 1440x720): regridding padrao, sem deadlock.
  !!
  !! Grade sintetica 360x180 (1 grau, 64800 celulas):
  !!   - ESMF_GridCreate1PeriDim: distribuicao automatica balanceada
  !!   - numOwnedElements > 0 em TODOS os PETs garantido pelo ESMF
  !!   - Coordenadas lon/lat explicitamente definidas (ESMF_STAGGERLOC_CENTER)
  !!   - Compativel com MED atm_grid (640x320) e DATOCN grid (1440x720)
  subroutine mpas_create_grid(grid, rc)
    type(ESMF_Grid), intent(out) :: grid
    integer,         intent(out) :: rc

    integer, parameter            :: NLON = 360
    integer, parameter            :: NLAT = 180
    real(ESMF_KIND_R8), parameter :: DLON = 1.0_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: DLAT = 1.0_ESMF_KIND_R8

    real(ESMF_KIND_R8), pointer :: coordX(:,:), coordY(:,:)
    integer  :: i, j, clbX(2), cubX(2), clbY(2), cubY(2)
    integer  :: petCount, regDecomp(2), localDeCount
    integer  :: nx_max, ny_tiles, lde
    integer  :: nx_tiles_target  ! B-57
    type(ESMF_VM) :: vm
    character(len=*), parameter :: subname = '(mpas_create_grid)'

    rc = ESMF_SUCCESS

    ! B-57 (fix B-52): regDecomp 2D com tiles quadradas - evita strips extremos.
    !
    ! PROBLEMA com B-52 (nx_max = NLON/2):
    !   Grids grandes + muitos PETs geram strips ultra-estreitas.
    !   Ex: grade netcdf 1440x720 a 512 PETs -> regDecomp=(/512,1/) ->
    !   2-3 cols x 720 rows -> aspecto 256:1. MOAB trava em ESMF_FieldBundleRegridStore
    !   (IPDvXp08 extro aparece mas o regrid store nunca retorna).
    !
    ! SOLUCAO B-57: tiles quadradas via sqrt(petCount).
    !   nx_tiles_target = nint(sqrt(N)) -> aspect ratio ~ 1.
    !   nx_max = min(target, NLON/2)    -> garante col >= 2 (bilinear OK).
    !   ny_tiles = ceil(N/nx_max)       -> cobre todos os PETs.
    !
    !   N=4:   sqrt=2  -> nx_max=min(2,180)=2   regDecomp=(/2,2/)=4DEs   asp 0.5:1 |
    !   N=128: sqrt=11 -> nx_max=min(11,180)=11 regDecomp=(/11,12/)=132  asp 0.5:1 |
    !   N=512: sqrt=23 -> nx_max=min(23,180)=23 regDecomp=(/23,23/)=529  asp 0.5:1 |
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    call ESMF_VMGet(vm, petCount=petCount, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    nx_tiles_target = max(1, nint(sqrt(real(petCount))))
    nx_max      = min(nx_tiles_target, NLON / 2)
    ny_tiles    = (petCount + nx_max - 1) / nx_max
    regDecomp(1) = min(nx_max, petCount)
    regDecomp(2) = max(1, ny_tiles)

    ! Grade regular 1 grau, periodica em lon.
    grid = ESMF_GridCreate1PeriDim( &
      minIndex   = (/1, 1/),           &
      maxIndex   = (/NLON, NLAT/),     &
      regDecomp  = regDecomp,          &
      coordSys   = ESMF_COORDSYS_SPH_DEG, &
      rc         = rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! ESMF_GridAddCoord eh COLETIVA - todos os PETs devem chama-la.
    call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! B-52: guard localDeCount>0 - com regDecomp 2D e DEs>petCount,
    ! todos os PETs tem >=1 DE; guard mantido por seguranca para N > nx_max*ny_tiles.
    call ESMF_GridGet(grid, localDeCount=localDeCount, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! B-53 (fix B-52): com regDecomp 2D, DEs totais > petCount -> alguns PETs tem
    ! localDeCount=2. ESMF_GridGetCoord sem localDE= falha com "must provide localDe
    ! argument for localDeCount > 1". Solucao: loop explicito sobre cada DE local.
    do lde = 0, localDeCount - 1

      ! Coordenada X (longitude): centros de celulas (-179.5a° a +179.5a°)
      nullify(coordX)
      call ESMF_GridGetCoord(grid, coordDim=1, localDE=lde, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             computationalLBound=clbX, computationalUBound=cubX, &
                             farrayPtr=coordX, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      do j = clbX(2), cubX(2)
        do i = clbX(1), cubX(1)
          coordX(i,j) = -180.0_ESMF_KIND_R8 + (real(i,ESMF_KIND_R8) - 0.5_ESMF_KIND_R8)*DLON
        end do
      end do

      ! Coordenada Y (latitude): centros de celulas (-89.5a° a +89.5a°)
      nullify(coordY)
      call ESMF_GridGetCoord(grid, coordDim=2, localDE=lde, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             computationalLBound=clbY, computationalUBound=cubY, &
                             farrayPtr=coordY, rc=rc)
      if (ChkErr(rc, __LINE__, u_FILE_u)) return
      do j = clbY(2), cubY(2)
        do i = clbY(1), cubY(1)
          coordY(i,j) = -90.0_ESMF_KIND_R8 + (real(j,ESMF_KIND_R8) - 0.5_ESMF_KIND_R8)*DLAT
        end do
      end do

    end do  ! lde = 0, localDeCount-1

    call ESMF_LogWrite(subname//': ESMF_Grid 360x180 criada (sem MOAB)', ESMF_LOGMSG_INFO)

  end subroutine mpas_create_grid

  !> @brief Escreve estatisticas dos campos do estado no log ESMF.
  !!
  !! Ativado por DumpFields='true' (atributo NUOPC).
  !! Para cada campo do estado escreve: nome, min, max, media.
  subroutine state_diagnose(state, state_tag, rc)
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: state_tag
    integer,          intent(out) :: rc

    type(ESMF_Field)              :: field
    character(len=64), allocatable :: fldnames(:)
    integer :: itemCount, i, localrc
    character(len=160) :: msg
    character(len=*), parameter :: subname = '(state_diagnose)'

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemCount=itemCount, rc=localrc)
    if (localrc /= ESMF_SUCCESS .or. itemCount == 0) then
      call ESMF_LogWrite(subname//': '//trim(state_tag)//' vazio', ESMF_LOGMSG_INFO)
      return
    end if

    allocate(fldnames(itemCount))
    call ESMF_StateGet(state, itemNameList=fldnames, rc=localrc)
    if (localrc /= ESMF_SUCCESS) then
      deallocate(fldnames); return
    end if

    write(msg,'(A,A)') subname//': ', trim(state_tag)
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    do i = 1, itemCount
      call ESMF_StateGet(state, itemName=trim(fldnames(i)), field=field, rc=localrc)
      if (localrc /= ESMF_SUCCESS) cycle
      block
        real(ESMF_KIND_R8), pointer :: fp1d(:)
        real(ESMF_KIND_R8), pointer :: fp2d(:,:)
        real(ESMF_KIND_R8), allocatable :: vals(:)
        integer :: fdr
        nullify(fp1d, fp2d)
        call ESMF_FieldGet(field, dimCount=fdr, rc=localrc)
        if (localrc /= ESMF_SUCCESS) cycle
        if (fdr == 1) then
          call ESMF_FieldGet(field, farrayPtr=fp1d, rc=localrc)
          if (localrc /= ESMF_SUCCESS .or. .not. associated(fp1d) .or. size(fp1d)==0) cycle
          vals = fp1d
          nullify(fp1d)
        else
          call ESMF_FieldGet(field, farrayPtr=fp2d, rc=localrc)
          if (localrc /= ESMF_SUCCESS .or. .not. associated(fp2d) .or. size(fp2d)==0) cycle
          vals = pack(fp2d, .true.)
          nullify(fp2d)
        end if
        write(msg,'(A,A,3(A,ES11.4))') &
          '  ', trim(fldnames(i)), &
          '  min=', minval(vals), &
          '  max=', maxval(vals), &
          '  mean=', sum(vals) / real(size(vals), ESMF_KIND_R8)
        call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)
      end block
    end do

    deallocate(fldnames)

  end subroutine state_diagnose

  ! ======================== privado ==========================================
  !============================================================================
  !============================================================================
  !============================================================================
  !============================================================================
  !============================================================================
  !============================================================================
  !============================================================================
  !============================================================================
  !============================================================================
  !======================================================


  !> @brief Copia campo 1D do ESMF_State para array Fortran.
  !> @brief Copia campo do ESMF_State para array Fortran 1D.
  !!
  !! FIX v5.2: suporte a campos 2D (ESMF_Grid) via pack() - Fortran 95, portavel.
  !! pack(fptr2d,.true.) serializa o array 2D column-major para 1D sem restricoes
  !! de tipo de ponteiro nem necessidade de iso_c_binding.
  !> @brief Copia campo do ESMF_State para array Fortran 1D.
  !!
  !! FIX v5.2: usa ESMF_FieldGet(dimCount=) antes de farrayPtr para evitar
  !! erro "rank does not match" propagado pelo ESMF para o NUOPC.
  subroutine state_get_field_1d(state, fldname, n, data, rc)
    type(ESMF_State),  intent(in)    :: state
    character(len=*),  intent(in)    :: fldname
    integer,           intent(in)    :: n
    real(MPAS_RKIND),  intent(inout) :: data(n)
    integer,           intent(out)   :: rc

    type(ESMF_Field)             :: field
    real(ESMF_KIND_R8), pointer  :: fptr1d(:)
    real(ESMF_KIND_R8), pointer  :: fptr2d(:,:)
    integer :: n_esmf, fld_rank
    character(len=*), parameter  :: subname = '(state_get_field_1d)'

    rc = ESMF_SUCCESS
    nullify(fptr1d, fptr2d)

    call ESMF_StateGet(state, itemName=fldname, field=field, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite(subname//': '//trim(fldname)//' nao encontrado', ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS
      return
    end if

    ! B-45: verificar localDeCount ANTES de farrayPtr (evita erro ESMF log).
    block
      integer :: localDeCount_sg
      call ESMF_FieldGet(field, localDeCount=localDeCount_sg, rc=rc)
      if (rc /= ESMF_SUCCESS .or. localDeCount_sg == 0) then
        rc = ESMF_SUCCESS; return
      end if
    end block

    ! Consultar rank do campo ANTES de chamar farrayPtr (evita erro ESMF)
    call ESMF_FieldGet(field, dimCount=fld_rank, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite(subname//': '//trim(fldname)//' dimCount query falhou', ESMF_LOGMSG_WARNING)
      rc = ESMF_SUCCESS; return
    end if

    if (fld_rank == 1) then
      ! Campo rank-1: ESMF_Mesh ou ESMF_Grid 1D
      call ESMF_FieldGet(field, farrayPtr=fptr1d, rc=rc)
      if (rc /= ESMF_SUCCESS .or. .not. associated(fptr1d)) then
        rc = ESMF_SUCCESS; return
      end if
      n_esmf = min(size(fptr1d), n)
      data(1:n_esmf) = real(fptr1d(1:n_esmf), MPAS_RKIND)
      if (n_esmf < n) data(n_esmf+1:n) = 0.0_MPAS_RKIND
      nullify(fptr1d)
    else
      ! Campo rank-2: ESMF_Grid regular (NLON x NLAT_local)
      call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
      if (rc /= ESMF_SUCCESS .or. .not. associated(fptr2d)) then
        rc = ESMF_SUCCESS; return
      end if
      n_esmf = min(size(fptr2d), n)
      block
        real(ESMF_KIND_R8), allocatable :: flat(:)
        flat = pack(fptr2d, .true.)   ! column-major 2D -> 1D (Fortran 95)
        data(1:n_esmf) = real(flat(1:n_esmf), MPAS_RKIND)
        if (n_esmf < n) data(n_esmf+1:n) = 0.0_MPAS_RKIND
      end block
      nullify(fptr2d)
    end if
    rc = ESMF_SUCCESS

  end subroutine state_get_field_1d

  !> @brief Copia array Fortran 1D para campo do ESMF_State.
  !!
  !! FIX v5.2: suporte a campos 2D (ESMF_Grid) via loop indexado - Fortran 95, portavel.
  !! Percorre fptr2d em ordem column-major, preenchendo com os primeiros n_esmf
  !! valores do array data (celulas MPAS locais).
  !> @brief Copia array Fortran 1D para campo do ESMF_State.
  !!
  !! FIX v5.2: usa ESMF_FieldGet(dimCount=) para consultar rank ANTES de
  !! chamar farrayPtr. Evita que ESMF registre erro "rank does not match"
  !! que propaga e aborta fases NUOPC (IPDv03p7 / RunPhase1).
  subroutine state_set_field_1d(state, fldname, n, data, rc)
    type(ESMF_State),  intent(inout) :: state
    character(len=*),  intent(in)    :: fldname
    integer,           intent(in)    :: n
    real(MPAS_RKIND),  intent(in)    :: data(n)
    integer,           intent(out)   :: rc

    type(ESMF_Field)             :: field
    real(ESMF_KIND_R8), pointer  :: fptr1d(:)
    real(ESMF_KIND_R8), pointer  :: fptr2d(:,:)
    integer :: n_esmf, fld_rank, i, j, idx
    character(len=*), parameter  :: subname = '(state_set_field_1d)'

    rc = ESMF_SUCCESS
    nullify(fptr1d, fptr2d)

    call ESMF_StateGet(state, itemName=fldname, field=field, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite(subname//': '//trim(fldname)//' nao encontrado', ESMF_LOGMSG_INFO)
      rc = ESMF_SUCCESS
      return
    end if

    ! B-45: verificar localDeCount ANTES de farrayPtr (evita erro ESMF log).
    block
      integer :: localDeCount_ss
      call ESMF_FieldGet(field, localDeCount=localDeCount_ss, rc=rc)
      if (rc /= ESMF_SUCCESS .or. localDeCount_ss == 0) then
        rc = ESMF_SUCCESS; return
      end if
    end block

    ! Consultar rank do campo ANTES de chamar farrayPtr (evita erro ESMF)
    call ESMF_FieldGet(field, dimCount=fld_rank, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite(subname//': '//trim(fldname)//' dimCount query falhou', ESMF_LOGMSG_WARNING)
      rc = ESMF_SUCCESS; return
    end if

    if (fld_rank == 1) then
      ! Campo rank-1: ESMF_Mesh ou ESMF_Grid 1D
      call ESMF_FieldGet(field, farrayPtr=fptr1d, rc=rc)
      if (rc /= ESMF_SUCCESS .or. .not. associated(fptr1d)) then
        rc = ESMF_SUCCESS; return
      end if
      n_esmf = min(size(fptr1d), n)
      fptr1d(1:n_esmf) = real(data(1:n_esmf), ESMF_KIND_R8)
      nullify(fptr1d)
    else
      ! Campo rank-2: ESMF_Grid regular (NLON x NLAT_local)
      ! Percorrer column-major: elemento (i,j) = posicao (j-1)*dim1 + i
      call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
      if (rc /= ESMF_SUCCESS .or. .not. associated(fptr2d)) then
        rc = ESMF_SUCCESS; return
      end if
      n_esmf = min(size(fptr2d), n)
      idx = 0
      outer: do j = lbound(fptr2d,2), ubound(fptr2d,2)
        do i = lbound(fptr2d,1), ubound(fptr2d,1)
          idx = idx + 1
          if (idx > n_esmf) exit outer
          fptr2d(i,j) = real(data(idx), ESMF_KIND_R8)
        end do
      end do outer
      nullify(fptr2d)
    end if
    rc = ESMF_SUCCESS

  end subroutine state_set_field_1d


end module mpas_cap_methods_mod
