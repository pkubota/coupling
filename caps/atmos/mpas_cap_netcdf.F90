!> @file mpas_cap_netcdf.F90
!! @brief Escrita dos campos exportados pelo MONAN-A em NetCDF com grade regular lat/lon.
!!
!! Versão 2.9 (13/04/2026) — GT Acoplamento MONAN / INPE/CGCT/DIMNT
!!
!! MUDANÇAS EM RELAÇÃO À v2.5:
!!
!!   Bug corrigido: timestamp duplo em export_write_netcdf (v2.6).
!!   Os parâmetros step e dt_s foram removidos da assinatura.
!!   O currTime passado por ModelRun (yr,mo,dy,hr,mn,sc via ESMF_ClockGet)
!!   é usado diretamente como timestamp do arquivo NetCDF.
!!   Resultado: arquivo monan_export_YYYYMMDD_HHMMSS.nc agora recebe o
!!   timestamp correto em vez de startTime + 2×step×dt_coupling.
!!
!! MUDANÇAS EM RELAÇÃO À v2.4:
!!
!!   Removida a lógica de conversão de campos acumulados (field_is_accumulated,
!!   acum_factor, ÷elapsed_s). Esta conversão foi movida para mpas_atm_model.F90
!!   (mpas_atm_run), que agora fornece campos já em unidades instantâneas:
!!
!!     Faxa_swdn = (acswdnb_N − acswdnb_{N−1}) / dt_coupling  [W/m²]
!!     Faxa_lwdn = (aclwdnb_N − aclwdnb_{N−1}) / dt_coupling  [W/m²]
!!     Faxa_prec = (rainnc_N + rainc_N − prev_N) / dt / 1000  [kg/m²/s]
!!     Faxa_taux = ρ_a · ust² · u10 / |V10|                   [N/m²]
!!     Faxa_tauy = ρ_a · ust² · v10 / |V10|                   [N/m²]
!!
!!   A divisão por elapsed_s (tempo total desde t=0) era incorreta para passos
!!   posteriores ao primeiro: produzia a média from t=0 em vez da média do
!!   intervalo de acoplamento corrente.
!!
!!   O limiar de outlier Faxa_taux/tauy foi reduzido de 1e4 para 10 N/m²
!!   (valor físico máximo realista de stress superficial). Com o cálculo
!!   correto via ust, não há mais lixo de memória nestes campos.
!!
!! ESTRUTURA DO ARQUIVO NetCDF GERADO (grade regular 1°×1°):
!!   dimensions  : lat(181), lon(360)
!!   variables   :
!!     double lat(lat)      [degrees_north, -90 a +90, passo 1°]
!!     double lon(lon)      [degrees_east, -180 a +179, passo 1°]
!!     double time          [escalar CF: seconds since start_time]
!!     double Sa_pslv(lat,lon)    [Pa,        instantâneo]
!!     double Sa_tbot(lat,lon)    [K,         instantâneo]
!!     double Sa_ubot(lat,lon)    [m s-1,     instantâneo]
!!     double Sa_vbot(lat,lon)    [m s-1,     instantâneo]
!!     double Faxa_swdn(lat,lon)  [W m-2,     média do intervalo de acoplamento]
!!     double Faxa_lwdn(lat,lon)  [W m-2,     média do intervalo de acoplamento]
!!     double Faxa_prec(lat,lon)  [kg m-2 s-1, média do intervalo de acoplamento]
!!     double Faxa_taux(lat,lon)  [N m-2,     ρ·ust²·u10/|V10|, outliers |v|>10 N/m² descartados]
!!     double Faxa_tauy(lat,lon)  [N m-2,     ρ·ust²·v10/|V10|, outliers |v|>10 N/m² descartados]
!!     double Faxa_lhflx(lat,lon) [W m-2,     instantâneo]
!!     double Faxa_shflx(lat,lon) [W m-2,     instantâneo]
!!
!! CONVENÇÃO DE DIMENSÕES (compatível com Python/netCDF4):
!!   Fortran: nf90_def_var([dimid_lon, dimid_lat]) → lon varia mais rápido
!!   Python:  nc.variables['Sa_tbot'][:] → shape (181, 360) = (nlat, nlon) ✓
!!
!! FLUXO MPI:
!!   Todos os PETs → MPI_Allgather (tamanhos locais)
!!               → MPI_Gatherv    (dados de campo → PET0)
!!   PET0        → voronoi_to_latlon (binning + conversão)
!!               → nf90_create / nf90_put_var / nf90_close

module mpas_cap_netcdf_mod

  use ESMF
  use mpi
  use netcdf
  use mpas_cap_utils_mod, only : ChkErr

  implicit none
  private

  ! ── Interface pública ──────────────────────────────────────────────────────
  public :: netcdf_init_coords    ! coleta coordenadas locais de todos os PETs
  public :: export_write_netcdf   ! interpola e escreve NetCDF com grade lat/lon
  public :: netcdf_config_set     ! configura grade e diretório a partir do namelist

  ! ── Grade regular de saída ────────────────────────────────────────────────
  ! Inicializada com 1°×1° (padrão do namelist &nuopc_netcdf).
  ! Atualizada via netcdf_config_set() antes de export_write_netcdf.
  integer,            save :: NLON     = 360   !  -180° a +179°
  integer,            save :: NLAT     = 181   !   -90° a  +90°
  real,               save :: GRID_RES = 1.0   !  resolução [°]
  real(ESMF_KIND_R8), save :: DLON = 1.0_ESMF_KIND_R8  ! passo em lon
  real(ESMF_KIND_R8), save :: DLAT = 1.0_ESMF_KIND_R8  ! passo em lat

  ! Fill value para pontos da grade sem nenhuma célula Voronoi
  real(ESMF_KIND_R8), parameter :: FILL_VALUE = -9.99e+33_ESMF_KIND_R8

  ! ── Coordenadas globais (módulo save — preenchidas por netcdf_init_coords)
  real(ESMF_KIND_R8), allocatable, save :: g_lon_global(:)  ! (nGlobal) graus
  real(ESMF_KIND_R8), allocatable, save :: g_lat_global(:)  ! (nGlobal) graus
  logical,                          save :: g_coords_ready = .false.

  ! ── Decomposição MPI salva (BUG FIX v2.8 — Bug 3) ───────────────────────
  ! Garante que allCounts/displs em export_write_netcdf sejam idênticos
  ! aos usados em netcdf_init_coords, evitando mapeamento geográfico errado.
  integer,                          save :: g_nLocal_saved  = 0
  integer,                          save :: g_nGlobal_saved = 0
  integer, allocatable,             save :: g_allCounts_saved(:)
  integer, allocatable,             save :: g_displs_saved(:)

  character(len=256), save :: OUTPUT_DIR = 'diag_export'
  character(len=*), parameter :: u_FILE_u   = __FILE__

contains

  !> @brief Configura os parâmetros de grade e diretório de saída a partir do namelist.
  !!
  !! Deve ser chamada em InitializeRealize antes de netcdf_init_coords.
  !! @param[in] res_deg   Resolução da grade em graus (ex: 1.0, 0.5, 0.25)
  !! @param[in] out_dir   Diretório de saída para os arquivos NetCDF
  !! @param[in] localPet  PET local do ESMF — suprime impressão em PETs > 0 (B-30)
  subroutine netcdf_config_set(res_deg, out_dir, localPet)
    real,             intent(in) :: res_deg
    character(len=*), intent(in) :: out_dir
    integer,          intent(in) :: localPet   ! B-30: guarda de rank

    GRID_RES   = res_deg
    DLON       = real(res_deg, ESMF_KIND_R8)
    DLAT       = real(res_deg, ESMF_KIND_R8)
    NLON       = nint(360.0 / res_deg)
    NLAT       = nint(180.0 / res_deg) + 1
    OUTPUT_DIR = trim(out_dir)

    ! B-30: sem guarda, N PETs × N chamadas = N² mensagens em stdout.
    ! Só PET 0 imprime; demais passam silenciosamente.
    if (localPet == 0) &
      write(*,'(A,F5.2,A,I0,A,I0,A,A)') &
        '[NetCDF] grade configurada: ', res_deg, '° -> NLON=', NLON, &
        ' NLAT=', NLAT, ' output_dir=', trim(OUTPUT_DIR)
  end subroutine netcdf_config_set

  ! ============================================================================
  !> Coleta coordenadas locais de todos os PETs via MPI_Gatherv e armazena no PET0.
  !!
  !! Deve ser chamada UMA VEZ em InitializeRealize do cap, após a malha ESMF
  !! estar disponível e ANTES da primeira chamada a export_write_netcdf.
  !!
  !! A decomposição MPI_Gatherv usada aqui é idêntica à de export_write_netcdf
  !! para os campos, garantindo que lon_global(i)/lat_global(i) corresponde
  !! exatamente ao dado(i) em recvBuf — eliminando o padrão de xadrez.
  !!
  !! Chamada idempotente (retorna imediatamente se já executada).
  !!
  !! Uso em mpas_cap.F90 (InitializeRealize), após ESMF_MeshGet:
  !!   call netcdf_init_coords(lon_local, lat_local, nLocalElem, vm, rc)
  subroutine netcdf_init_coords(lon_local, lat_local, nLocal, vm, rc)
    real(ESMF_KIND_R8), intent(in)    :: lon_local(:)   ! longitudes do PET (graus)
    real(ESMF_KIND_R8), intent(in)    :: lat_local(:)   ! latitudes  do PET (graus)
    integer,            intent(in)    :: nLocal          ! número de células locais
    type(ESMF_VM),      intent(in)    :: vm
    integer,            intent(inout) :: rc

    integer :: localPet, petCount, mpiComm, mpi_ierr, i, nGlobal
    integer, allocatable :: allCounts(:), displs(:)
    character(len=*), parameter :: subname = '(netcdf_init_coords)'

    rc = ESMF_SUCCESS
    if (g_coords_ready) return   ! idempotente

    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, &
                    mpiCommunicator=mpiComm, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! ── Reunir tamanhos locais de cada PET ────────────────────────────────
    allocate(allCounts(petCount))
    call MPI_Allgather(nLocal, 1, MPI_INTEGER, &
                       allCounts, 1, MPI_INTEGER, mpiComm, mpi_ierr)
    if (mpi_ierr /= MPI_SUCCESS) then
      call ESMF_LogSetError(ESMF_FAILURE, msg=subname//': MPI_Allgather falhou', &
           line=__LINE__, file=u_FILE_u, rcToReturn=rc)
      return
    end if
    nGlobal = sum(allCounts)

    allocate(displs(petCount))
    displs(1) = 0
    do i = 2, petCount
      displs(i) = displs(i-1) + allCounts(i-1)
    end do

    ! ── Alocar buffers (PETs >0 recebem array mínimo — argumento inativo) ─
    if (localPet == 0) then
      allocate(g_lon_global(nGlobal))
      allocate(g_lat_global(nGlobal))
    else
      allocate(g_lon_global(1))
      allocate(g_lat_global(1))
    end if

    ! ── Gather de lon e lat ────────────────────────────────────────────────
    call MPI_Gatherv(lon_local, nLocal, MPI_DOUBLE_PRECISION, &
                     g_lon_global, allCounts, displs, MPI_DOUBLE_PRECISION, &
                     0, mpiComm, mpi_ierr)
    if (mpi_ierr /= MPI_SUCCESS .and. localPet == 0) &
      write(*,'(A)') '[NetCDF] AVISO: MPI_Gatherv de lon_local falhou'

    call MPI_Gatherv(lat_local, nLocal, MPI_DOUBLE_PRECISION, &
                     g_lat_global, allCounts, displs, MPI_DOUBLE_PRECISION, &
                     0, mpiComm, mpi_ierr)
    if (mpi_ierr /= MPI_SUCCESS .and. localPet == 0) &
      write(*,'(A)') '[NetCDF] AVISO: MPI_Gatherv de lat_local falhou'

    ! BUG FIX v2.8: salvar decomposicao MPI para reuso em export_write_netcdf.
    ! Garante que recvBuf(i) corresponde a g_lon/lat_global(i) — sem desfase.
    g_nLocal_saved  = nLocal
    g_nGlobal_saved = nGlobal
    allocate(g_allCounts_saved(petCount))
    allocate(g_displs_saved(petCount))
    g_allCounts_saved = allCounts
    g_displs_saved    = displs

    deallocate(allCounts, displs)
    g_coords_ready = .true.

    if (localPet == 0) then
      write(*,'(A,I0,A)') &
        '[NetCDF] Coordenadas prontas: ', nGlobal, ' células (grade 1° pronta)'
      call ESMF_LogWrite(subname//': '//trim(int_to_str(nGlobal))// &
                         ' células — interpolação lat/lon ativa', ESMF_LOGMSG_INFO)
    end if

  end subroutine netcdf_init_coords

  ! ============================================================================
  !> Reúne campos do exportState, interpola para grade lat/lon e escreve NetCDF.
  !!
  !! Fluxo por campo (todos os PETs participam nas chamadas MPI coletivas):
  !!   1. MPI_Gatherv → recvBuf(nGlobal) no PET0
  !!   2. voronoi_to_latlon (PET0): nearest-neighbor binning + filtro de outliers
  !!   3. nf90_put_var: escreve grid(NLON,NLAT) → Python lê como (nlat,nlon)
  !!
  !! Todos os campos chegam já em unidades instantâneas de mpas_atm_model.F90
  !! (sem conversão de acumulados aqui).
  subroutine export_write_netcdf(exportState,               &
                                  elapsed_s,                &
                                  s_yr, s_mo, s_dy,         &
                                  s_hr, s_mn, s_sc,         &
                                  vm, rc)
    type(ESMF_State), intent(in)    :: exportState
    integer,          intent(in)    :: elapsed_s
    integer,          intent(in)    :: s_yr, s_mo, s_dy, s_hr, s_mn, s_sc
    type(ESMF_VM),    intent(in)    :: vm
    integer,          intent(inout) :: rc

    ! Locais
    integer :: localPet, petCount, mpiComm, mpi_ierr
    integer :: i, itemCount, nLocal, nGlobal
    integer :: ncid, varid, ncstat
    integer :: dimid_lat, dimid_lon
    integer :: varid_lat, varid_lon, varid_t
    integer :: c_yr, c_mo, c_dy, c_hr, c_mn, c_sc
    integer :: st_yr, st_mo, st_dy, st_hr, st_mn, st_sc

    character(len=64),  allocatable :: fldnames(:)
    integer,            allocatable :: allCounts(:), displs(:)
    real(ESMF_KIND_R8), allocatable :: sendBuf(:), recvBuf(:)

    real(ESMF_KIND_R8) :: grid_2d(NLON, NLAT)
    real(ESMF_KIND_R8) :: lat_axis(NLAT), lon_axis(NLON)
    real(ESMF_KIND_R8) :: time_val, othr

    type(ESMF_Field) :: field
    character(len=36) :: fname_base
    character(len=64) :: fname
    character(len=19) :: valid_time_iso, time_units_str
    integer :: cmd_stat
    character(len=*), parameter :: subname = '(export_write_netcdf)'

    rc = ESMF_SUCCESS

    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, &
                    mpiCommunicator=mpiComm, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! Coordenadas requeridas — preenchidas por netcdf_init_coords em InitializeRealize
    if (.not. g_coords_ready) then
      if (localPet == 0) write(*,'(A)') &
        '[NetCDF] ERRO: netcdf_init_coords nao foi chamado em InitializeRealize.'
      rc = ESMF_FAILURE
      return
    end if

    ! ── 0. Data/hora do passo ─────────────────────────────────────────────
    ! s_yr..s_sc = currTime (ESMF_ClockGet em ModelRun) → nome do arquivo.
    ! elapsed_s  = step_count * dt_coupling_s (calculado pelo chamador)
    !            → variável CF time: "elapsed_s seconds since startTime".
    c_yr = s_yr; c_mo = s_mo; c_dy = s_dy
    c_hr = s_hr; c_mn = s_mn; c_sc = s_sc

    ! ── 1. Inventário do exportState ──────────────────────────────────────
    call ESMF_StateGet(exportState, itemCount=itemCount, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return
    if (itemCount == 0) return

    allocate(fldnames(itemCount))
    call ESMF_StateGet(exportState, itemNameList=fldnames, rc=rc)
    if (ChkErr(rc, __LINE__, u_FILE_u)) return

    ! ── 2-3. Decomposição MPI reutilizada de netcdf_init_coords (BUG FIX v2.8)
    ! Antes: nLocal = size(fptr) [localCells_ESMF] podia diferir do nLocal
    ! usado em netcdf_init_coords [min(localCells_ESMF, nCells_MPAS)],
    ! gerando allCounts distintos e mapeamento geográfico errado no NetCDF.
    ! Agora: reutilizar exatamente os mesmos allCounts/displs do init_coords
    ! garante que recvBuf(i) corresponde a g_lon/lat_global(i).
    nLocal  = g_nLocal_saved
    nGlobal = g_nGlobal_saved
    allocate(allCounts(petCount), displs(petCount))
    allCounts = g_allCounts_saved
    displs    = g_displs_saved
    rc = ESMF_SUCCESS

    ! ── 4. Buffers ────────────────────────────────────────────────────────
    allocate(sendBuf(max(nLocal, 1)))
    if (localPet == 0) then
      allocate(recvBuf(max(nGlobal, 1)))
    else
      allocate(recvBuf(1))
    end if

    ! ── 5. PET0: criar e definir estrutura do arquivo NetCDF ──────────────
    if (localPet == 0) then
      fname_base     = datetime_to_fname(c_yr,c_mo,c_dy,c_hr,c_mn,c_sc)
      fname          = trim(OUTPUT_DIR)//'/'//trim(fname_base)
      valid_time_iso = datetime_to_iso(c_yr,c_mo,c_dy,c_hr,c_mn,c_sc)
      ! startTime = currTime - elapsed_s (para CF time_units "seconds since startTime")
      call datetime_add_seconds(s_yr,s_mo,s_dy,s_hr,s_mn,s_sc, -elapsed_s, &
                                 st_yr,st_mo,st_dy,st_hr,st_mn,st_sc)
      time_units_str = datetime_to_cf_base(st_yr,st_mo,st_dy,st_hr,st_mn,st_sc)
      time_val       = real(elapsed_s, ESMF_KIND_R8)

      call execute_command_line('mkdir -p '//trim(OUTPUT_DIR), exitstat=cmd_stat)

      ncstat = nf90_create(trim(fname), NF90_CLOBBER, ncid)
      if (ncstat /= NF90_NOERR) then
        write(*,'(3A)') '[NetCDF] ERRO ao criar ', trim(fname), &
          ': '//trim(nf90_strerror(ncstat))
        call ESMF_LogSetError(ESMF_FAILURE, msg=subname//': nf90_create falhou', &
             line=__LINE__, file=u_FILE_u, rcToReturn=rc)
        return
      end if

      ! ── Atributos globais CF-1.8 ────────────────────────────────────────
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions',    'CF-1.8')
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'title',          &
               'MPAS-A exportState — grade regular 1° lat/lon — NUOPC/CMEPS')
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'source',         &
               'NUOPC-MPAS-Integrado v5.2 (mpas_cap_netcdf_mod v2.7)')
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'institution',    &
               'INPE / CGCT / DIMNT')
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'valid_time',     &
               trim(valid_time_iso))
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'start_time',     &
               datetime_to_iso(s_yr,s_mo,s_dy,s_hr,s_mn,s_sc))
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'elapsed_time_s', elapsed_s)
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'ncells_global',  nGlobal)
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'petCount',       petCount)
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'grid_resolution','1.0 degree')
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'interp_method',  &
               'Nearest-neighbor binning, Voronoi x1.40962 (~120 km) -> 1 deg')
      ncstat = nf90_put_att(ncid, NF90_GLOBAL, 'processing_note', &
               'Faxa_swdn/lwdn/prec: media do intervalo de acoplamento (incremento/dt). ' // &
               'Faxa_taux/tauy: rho*ust^2*(u,v)/|V10|, outliers |v|>10 N/m2 descartados. ' // &
               'time: seconds since startTime (CF-1.8).')

      ! ── Dimensões ────────────────────────────────────────────────────────
      ! lat e lon — sem dimensão time (1 arquivo por passo)
      ncstat = nf90_def_dim(ncid, 'lat', NLAT, dimid_lat)
      if (ncstat /= NF90_NOERR) then
        write(*,'(A)') '[NetCDF] ERRO nf90_def_dim(lat): '//trim(nf90_strerror(ncstat))
        ncstat = nf90_close(ncid)
        call ESMF_LogSetError(ESMF_FAILURE, msg=subname//': nf90_def_dim falhou', &
             line=__LINE__, file=u_FILE_u, rcToReturn=rc); return
      end if
      ncstat = nf90_def_dim(ncid, 'lon', NLON, dimid_lon)

      ! ── Variáveis de coordenada ──────────────────────────────────────────
      ncstat = nf90_def_var(ncid, 'lat', NF90_DOUBLE, [dimid_lat], varid_lat)
      ncstat = nf90_put_att(ncid, varid_lat, 'long_name',     'latitude')
      ncstat = nf90_put_att(ncid, varid_lat, 'units',         'degrees_north')
      ncstat = nf90_put_att(ncid, varid_lat, 'standard_name', 'latitude')
      ncstat = nf90_put_att(ncid, varid_lat, 'axis',          'Y')

      ncstat = nf90_def_var(ncid, 'lon', NF90_DOUBLE, [dimid_lon], varid_lon)
      ncstat = nf90_put_att(ncid, varid_lon, 'long_name',     'longitude')
      ncstat = nf90_put_att(ncid, varid_lon, 'units',         'degrees_east')
      ncstat = nf90_put_att(ncid, varid_lon, 'standard_name', 'longitude')
      ncstat = nf90_put_att(ncid, varid_lon, 'axis',          'X')

      ! Variável time escalar CF
      ncstat = nf90_def_var(ncid, 'time', NF90_DOUBLE, varid_t)
      ncstat = nf90_put_att(ncid, varid_t, 'long_name', &
               'tempo da simulacao ao final do passo de acoplamento')
      ncstat = nf90_put_att(ncid, varid_t, 'units',     &
               'seconds since '//trim(time_units_str))
      ncstat = nf90_put_att(ncid, varid_t, 'calendar',  'gregorian')
      ncstat = nf90_put_att(ncid, varid_t, 'valid_time',trim(valid_time_iso))

      ! ── Variáveis dos campos (lon, lat) em Fortran column-major ──────────
      ! Python: nc['campo'][:] → shape (NLAT, NLON) = (181, 360)  ✓
      do i = 1, itemCount
        ncstat = nf90_def_var(ncid, trim(fldnames(i)), NF90_DOUBLE, &
                              [dimid_lon, dimid_lat], varid)
        if (ncstat /= NF90_NOERR) then
          write(*,'(3A)') '[NetCDF] AVISO nf90_def_var: ', &
            trim(fldnames(i)), ' — '//trim(nf90_strerror(ncstat))
          cycle
        end if
        ncstat = nf90_put_att(ncid, varid, 'long_name',     &
                 field_long_name(fldnames(i)))
        ncstat = nf90_put_att(ncid, varid, 'units',         &
                 field_units(fldnames(i)))
        ncstat = nf90_put_att(ncid, varid, 'standard_name', &
                 field_stdname(fldnames(i)))
        ncstat = nf90_put_att(ncid, varid, 'CMEPS_name',    trim(fldnames(i)))
        ncstat = nf90_put_att(ncid, varid, '_FillValue',    FILL_VALUE)
        ncstat = nf90_put_att(ncid, varid, 'missing_value', FILL_VALUE)
      end do

      ncstat = nf90_enddef(ncid)
      if (ncstat /= NF90_NOERR) then
        write(*,'(A)') '[NetCDF] ERRO nf90_enddef: '//trim(nf90_strerror(ncstat))
        ncstat = nf90_close(ncid)
        call ESMF_LogSetError(ESMF_FAILURE, msg=subname//': nf90_enddef falhou', &
             line=__LINE__, file=u_FILE_u, rcToReturn=rc); return
      end if

      ! ── Escrever eixos e time ────────────────────────────────────────────
      do i = 1, NLAT
        lat_axis(i) = -90.0_ESMF_KIND_R8 + real(i-1, ESMF_KIND_R8) * DLAT
      end do
      do i = 1, NLON
        lon_axis(i) = -180.0_ESMF_KIND_R8 + real(i-1, ESMF_KIND_R8) * DLON
      end do
      ncstat = nf90_put_var(ncid, varid_lat, lat_axis)
      ncstat = nf90_put_var(ncid, varid_lon, lon_axis)
      ncstat = nf90_put_var(ncid, varid_t,   time_val)

    end if   ! localPet == 0

    ! ── 6. Loop por campo: gather → interpolação → escrita ─────────────────
    ! MPI_Gatherv é coletivo — TODOS os PETs participam em cada iteração.
    do i = 1, itemCount

      ! FIX v2.9: dimCount guard antes de farrayPtr (campos rank-2 ESMF_Grid)
      ! ESMF_FieldGet(farrayPtr=rank1_ptr) em campo rank-2 gera erro ESMF.
      sendBuf(1:max(nLocal,1)) = 0.0_ESMF_KIND_R8   ! default: zeros
      call ESMF_StateGet(exportState, itemName=trim(fldnames(i)), &
                         field=field, rc=rc)
      if (rc == ESMF_SUCCESS) then
        block
          real(ESMF_KIND_R8), pointer     :: fptr1d(:)
          real(ESMF_KIND_R8), pointer     :: fptr2d(:,:)
          integer :: fld_rank
          nullify(fptr1d, fptr2d)
          call ESMF_FieldGet(field, dimCount=fld_rank, rc=rc)
          if (rc == ESMF_SUCCESS) then
            if (fld_rank == 1) then
              call ESMF_FieldGet(field, farrayPtr=fptr1d, rc=rc)
              if (rc == ESMF_SUCCESS .and. associated(fptr1d) .and. &
                  size(fptr1d) >= nLocal) &
                sendBuf(1:nLocal) = fptr1d(1:nLocal)
              if (associated(fptr1d)) nullify(fptr1d)
            else
              call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
              if (rc == ESMF_SUCCESS .and. associated(fptr2d) .and. &
                  size(fptr2d) >= nLocal) then
                block
                  real(ESMF_KIND_R8), allocatable :: flat(:)
                  flat = pack(fptr2d, .true.)
                  sendBuf(1:nLocal) = flat(1:nLocal)
                end block
              end if
              if (associated(fptr2d)) nullify(fptr2d)
            end if
          end if
        end block
      end if
      rc = ESMF_SUCCESS

      ! Gather coletivo → PET0
      call MPI_Gatherv(sendBuf, nLocal, MPI_DOUBLE_PRECISION, &
                       recvBuf, allCounts, displs,            &
                       MPI_DOUBLE_PRECISION, 0, mpiComm, mpi_ierr)
      if (mpi_ierr /= MPI_SUCCESS .and. localPet == 0) &
        write(*,'(3A)') '[NetCDF] AVISO MPI_Gatherv: ', trim(fldnames(i)), &
          ' — dados podem estar incompletos'

      ! PET0: interpola e escreve
      if (localPet == 0) then

        ! Limiar de outlier por campo
        othr = field_outlier_threshold(fldnames(i))

        ! Nearest-neighbor binning com filtro de outliers
        call voronoi_to_latlon(recvBuf(1:nGlobal),      &
                                g_lon_global(1:nGlobal), &
                                g_lat_global(1:nGlobal), &
                                nGlobal, grid_2d, othr)

        ncstat = nf90_inq_varid(ncid, trim(fldnames(i)), varid)
        if (ncstat == NF90_NOERR) then
          ncstat = nf90_put_var(ncid, varid, grid_2d)
          if (ncstat /= NF90_NOERR) &
            write(*,'(3A)') '[NetCDF] AVISO nf90_put_var: ', &
              trim(fldnames(i)), ' — '//trim(nf90_strerror(ncstat))
        end if

      end if
    end do   ! campos

    ! ── 7. PET0: fechar arquivo ───────────────────────────────────────────
    if (localPet == 0) then
      ncstat = nf90_close(ncid)
      if (ncstat == NF90_NOERR) then
        write(*,'(A,4A)') '[NetCDF] Escrito (', trim(valid_time_iso), ') → ', trim(fname), ''
        call ESMF_LogWrite(subname//': '//trim(fname)//' escrito', &
                           ESMF_LOGMSG_INFO)
      else
        write(*,'(A)') '[NetCDF] AVISO nf90_close: '//trim(nf90_strerror(ncstat))
      end if
    end if

    deallocate(fldnames, allCounts, displs, sendBuf, recvBuf)

  end subroutine export_write_netcdf

  ! ============================================================================
  !> Nearest-neighbor binning: mapeia células Voronoi para grade regular NLON×NLAT.
  !!
  !! Para cada célula k:
  !!   1. Descartar se |data_in(k)| > outlier_thr (fill value ou lixo de memória)
  !!   2. Normalizar longitude para [-180, 180)
  !!   3. Calcular índice Fortran (ilon, ilat) do ponto de grade mais próximo:
  !!        ilon = nint((lon + 180) / DLON) + 1   [1..NLON]
  !!        ilat = nint((lat +  90) / DLAT) + 1   [1..NLAT]
  !!   4. Acumular soma e contagem
  !!
  !! Após o loop:
  !!   grid_out(ilon,ilat) = soma / contagem   onde cnt > 0
  !!   grid_out(ilon,ilat) = FILL_VALUE        onde cnt = 0
  !!
  !! Spray adaptativo em longitude (ver detalhes nos comentários inline).
  !! Todos os campos chegam já em unidades corretas — sem conversão aqui.
  subroutine voronoi_to_latlon(data_in, lon_v, lat_v, n, grid_out, outlier_thr)
    real(ESMF_KIND_R8), intent(in)  :: data_in(n), lon_v(n), lat_v(n)
    integer,            intent(in)  :: n
    real(ESMF_KIND_R8), intent(out) :: grid_out(NLON, NLAT)
    real(ESMF_KIND_R8), intent(in)  :: outlier_thr

    ! Spray fixo em latitude (espaçamento uniforme)
    integer, parameter :: NSPAN_LAT = 1

    ! Constante para spray adaptativo em longitude
    ! nspan_lon = ceil(0.60 / (cos(lat) * DLON))   [DLON=1°]
    ! Fator 0.60 > 0.54 = margem de segurança para células ~120 km
    real(ESMF_KIND_R8), parameter :: CELL_HALF_DEG = 0.60_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: PI            = 3.14159265358979323846_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: DEG2RAD       = PI / 180.0_ESMF_KIND_R8

    real(ESMF_KIND_R8) :: acc(NLON, NLAT), val, lon_norm, cos_lat
    integer            :: cnt(NLON, NLAT), k, ilon_c, ilat_c, ilon2, ilat2
    integer            :: di, dj, nspan_lon

    acc      = 0.0_ESMF_KIND_R8
    cnt      = 0
    grid_out = FILL_VALUE

    do k = 1, n
      val = data_in(k)

      ! Filtrar fill values (-9.99e33) e outliers campo-específicos
      if (abs(val) > outlier_thr) cycle

      ! Normalizar longitude: qualquer convenção → [-180, 180)
      lon_norm = lon_v(k)
      do while (lon_norm >= 180.0_ESMF_KIND_R8)
        lon_norm = lon_norm - 360.0_ESMF_KIND_R8
      end do
      do while (lon_norm < -180.0_ESMF_KIND_R8)
        lon_norm = lon_norm + 360.0_ESMF_KIND_R8
      end do

      ! Ponto de grade central mais próximo (1-based)
      ilon_c = nint((lon_norm + 180.0_ESMF_KIND_R8) / DLON) + 1
      ilat_c = nint((lat_v(k) +  90.0_ESMF_KIND_R8) / DLAT) + 1
      ilon_c = min(max(ilon_c, 1), NLON)
      ilat_c = min(max(ilat_c, 1), NLAT)

      ! Spray longitude adaptativo: compensa convergência de meridianos nos polos.
      ! cos_lat mínimo = sin(0.5°) ≈ 0.009 para evitar nspan_lon → ∞ exatamente no polo.
      cos_lat   = max(cos(lat_v(k) * DEG2RAD), 0.009_ESMF_KIND_R8)
      nspan_lon = min(int(CELL_HALF_DEG / (cos_lat * DLON)) + 1, NLON/2)
      nspan_lon = max(nspan_lon, NSPAN_LAT)   ! pelo menos ±1

      ! Spray: acumular em janela (2*NSPAN_LAT+1) × (2*nspan_lon+1)
      do dj = -NSPAN_LAT, NSPAN_LAT
        ilat2 = min(max(ilat_c + dj, 1), NLAT)
        do di = -nspan_lon, nspan_lon
          ! Longitude: tratamento periódico (wrap ±180°)
          ilon2 = ilon_c + di
          if (ilon2 < 1)    ilon2 = ilon2 + NLON
          if (ilon2 > NLON) ilon2 = ilon2 - NLON
          acc(ilon2, ilat2) = acc(ilon2, ilat2) + val
          cnt(ilon2, ilat2) = cnt(ilon2, ilat2) + 1
        end do
      end do
    end do

    ! Média das contribuições; pontos sem nenhuma célula ficam com FILL_VALUE
    where (cnt > 0) grid_out = acc / real(cnt, ESMF_KIND_R8)

  end subroutine voronoi_to_latlon

  ! ============================================================================
  ! Funções auxiliares de metadados de campo
  ! ============================================================================

  !> Limiar de outlier por campo.
  !! Faxa_taux/tauy: máximo físico realista de stress superficial ≈ 3–5 N/m²
  !!   (furacão Cat.5: ~3 N/m²); limiar = 10 N/m² com margem.
  !!   Com cálculo correto via ust, não há mais lixo de memória — apenas
  !!   células com wind-shear extremo podem ultrapassar 5 N/m².
  !! Demais campos: 1e30 (captura apenas fill value -9.99e33).
  !> Limiar de outlier por campo (filtra lixo de memória e fill values).
  !!
  !! Faxa_taux/tauy: stress superficial máximo físico ≈ 3–5 N/m² (furacão Cat.5);
  !!   limiar = 10 N/m² com margem de segurança.
  !!
  !! Sa_u10m_mpas / Sa_v10m_mpas: vento 10 m. Valores > 10 m/s são NORMAIS
  !!   (jatos de baixos níveis, alísios fortes, ciclones extratropicais);
  !!   usar 10 m/s cortava 2.6% dos bins — justamente os de vento forte —
  !!   reduzindo σ_cap em 7% vs σ_standalone (Bug B-28).
  !!   Limiar correto: 150 m/s (fisicamente impossível → só filtra garbage).
  !!
  !! Demais campos: 1e30 (captura apenas fill value -9.99e33).
  pure real(ESMF_KIND_R8) function field_outlier_threshold(fname)
    character(len=*), intent(in) :: fname
    if (trim(fname) == 'Faxa_taux' .or. trim(fname) == 'Faxa_tauy') then
      field_outlier_threshold = 10.0_ESMF_KIND_R8    ! N/m²: stress max Cat.5
    else if (trim(fname) == 'Sa_u10m_mpas' .or. &
             trim(fname) == 'Sa_v10m_mpas') then
      field_outlier_threshold = 150.0_ESMF_KIND_R8   ! m/s: impossível fisicamente
    else
      field_outlier_threshold = 1.0e30_ESMF_KIND_R8  ! filtra só fill value
    end if
  end function field_outlier_threshold

  !> Unidades dos campos após conversão (corrige metadados do MPAS).
  function field_units(fname) result(units)
    character(len=*), intent(in) :: fname
    character(len=32) :: units
    select case (trim(fname))
      ! Nomes _mpas (cap NUOPC v3+, sufixo identifica fonte MPAS vs DATM)
      case ('Sa_pslv_mpas')                          ; units = 'Pa'
      case ('Sa_tbot_mpas')                          ; units = 'K'
      case ('Sa_u10m_mpas', 'Sa_v10m_mpas')          ; units = 'm s-1'
      case ('Sa_shum_mpas')                          ; units = 'kg kg-1'
      case ('Faxa_swdn_mpas', 'Faxa_lwdn_mpas')      ; units = 'W m-2'
      case ('Faxa_rain_mpas', 'Faxa_snow_mpas')      ; units = 'kg m-2 s-1'
      ! Nomes legado (sem _mpas) para compatibilidade retroativa
      case ('Sa_pslv')                               ; units = 'Pa'
      case ('Sa_tbot')                               ; units = 'K'
      case ('Sa_ubot', 'Sa_vbot')                    ; units = 'm s-1'
      case ('Faxa_swdn', 'Faxa_lwdn')                ; units = 'W m-2'
      case ('Faxa_prec')                             ; units = 'kg m-2 s-1'
      case ('Faxa_taux', 'Faxa_tauy')                ; units = 'N m-2'
      case ('Faxa_lhflx', 'Faxa_shflx')             ; units = 'W m-2'
      case default                                   ; units = '1'
    end select
  end function field_units

  !> Long name descritivo para cada campo CMEPS.
  function field_long_name(fname) result(lname)
    character(len=*), intent(in) :: fname
    character(len=96) :: lname
    select case (trim(fname))
      ! Nomes _mpas
      case ('Sa_pslv_mpas')    ; lname = 'Pressao ao nivel do mar'
      case ('Sa_tbot_mpas')    ; lname = 'Temperatura do ar a 2 m'
      case ('Sa_u10m_mpas')    ; lname = 'Vento zonal a 10 m'
      case ('Sa_v10m_mpas')    ; lname = 'Vento meridional a 10 m'
      case ('Sa_shum_mpas')    ; lname = 'Umidade especifica a 2 m'
      case ('Faxa_swdn_mpas')  ; lname = 'Radiacao SW descendente media no intervalo'
      case ('Faxa_lwdn_mpas')  ; lname = 'Radiacao LW descendente media no intervalo'
      case ('Faxa_rain_mpas')  ; lname = 'Precipitacao liquida media no intervalo'
      case ('Faxa_snow_mpas')  ; lname = 'Precipitacao solida (neve) media no intervalo'
      ! Nomes legado
      case ('Sa_pslv')    ; lname = 'Pressao ao nivel do mar'
      case ('Sa_tbot')    ; lname = 'Temperatura do ar a 2 m'
      case ('Sa_ubot')    ; lname = 'Vento zonal a 10 m'
      case ('Sa_vbot')    ; lname = 'Vento meridional a 10 m'
      case ('Faxa_swdn')  ; lname = 'Radiacao SW descendente media no intervalo de acoplamento'
      case ('Faxa_lwdn')  ; lname = 'Radiacao LW descendente media no intervalo de acoplamento'
      case ('Faxa_prec')  ; lname = 'Taxa de precipitacao total media no intervalo de acoplamento'
      case ('Faxa_taux')  ; lname = 'Tensao de cisalhamento zonal na superficie'
      case ('Faxa_tauy')  ; lname = 'Tensao de cisalhamento meridional na superficie'
      case ('Faxa_lhflx') ; lname = 'Fluxo de calor latente na superficie'
      case ('Faxa_shflx') ; lname = 'Fluxo de calor sensivel na superficie'
      case default        ; lname = trim(fname)
    end select
  end function field_long_name

  !> CF standard_name para campos CMEPS.
  function field_stdname(fname) result(sname)
    character(len=*), intent(in) :: fname
    character(len=80) :: sname
    select case (trim(fname))
      ! Nomes _mpas
      case ('Sa_pslv_mpas')    ; sname = 'air_pressure_at_mean_sea_level'
      case ('Sa_tbot_mpas')    ; sname = 'air_temperature'
      case ('Sa_u10m_mpas')    ; sname = 'eastward_wind'
      case ('Sa_v10m_mpas')    ; sname = 'northward_wind'
      case ('Sa_shum_mpas')    ; sname = 'specific_humidity'
      case ('Faxa_swdn_mpas')  ; sname = 'surface_downwelling_shortwave_flux_in_air'
      case ('Faxa_lwdn_mpas')  ; sname = 'surface_downwelling_longwave_flux_in_air'
      case ('Faxa_rain_mpas')  ; sname = 'rainfall_flux'
      case ('Faxa_snow_mpas')  ; sname = 'snowfall_flux'
      ! Nomes legado
      case ('Sa_pslv')    ; sname = 'air_pressure_at_mean_sea_level'
      case ('Sa_tbot')    ; sname = 'air_temperature'
      case ('Sa_ubot')    ; sname = 'eastward_wind'
      case ('Sa_vbot')    ; sname = 'northward_wind'
      case ('Faxa_swdn')  ; sname = 'surface_downwelling_shortwave_flux_in_air'
      case ('Faxa_lwdn')  ; sname = 'surface_downwelling_longwave_flux_in_air'
      case ('Faxa_prec')  ; sname = 'precipitation_flux'
      case ('Faxa_taux')  ; sname = 'surface_downward_eastward_stress'
      case ('Faxa_tauy')  ; sname = 'surface_downward_northward_stress'
      case ('Faxa_lhflx') ; sname = 'surface_upward_latent_heat_flux'
      case ('Faxa_shflx') ; sname = 'surface_upward_sensible_heat_flux'
      case default        ; sname = 'unknown'
    end select
  end function field_stdname

  ! ── Utilitários de conversão ──────────────────────────────────────────────

  function int_to_str(n) result(s)
    integer, intent(in) :: n; character(len=16) :: s
    write(s,'(I0)') n
  end function int_to_str

  ! ── Aritmética de calendário gregoriano proléptico ────────────────────────

  pure logical function is_leap_year(yr)
    integer, intent(in) :: yr
    is_leap_year = (mod(yr,4)==0 .and. mod(yr,100)/=0) .or. mod(yr,400)==0
  end function is_leap_year

  subroutine datetime_add_seconds(yr,mo,dy,hr,mn,sc, nadd, &
                                   yr_o,mo_o,dy_o,hr_o,mn_o,sc_o)
    integer, intent(in)  :: yr,mo,dy,hr,mn,sc, nadd
    integer, intent(out) :: yr_o,mo_o,dy_o,hr_o,mn_o,sc_o
    integer :: dim(12), tot, extra, dsec
    tot   = hr*3600 + mn*60 + sc + nadd
    extra = tot / 86400
    dsec  = mod(tot, 86400)
    if (dsec < 0) then; extra = extra - 1; dsec = dsec + 86400; end if
    sc_o = mod(dsec, 60); mn_o = mod(dsec/60, 60); hr_o = dsec/3600
    yr_o = yr; mo_o = mo; dy_o = dy + extra
    dim  = [31,28,31,30,31,30,31,31,30,31,30,31]
    if (is_leap_year(yr_o)) dim(2) = 29
    do while (dy_o > dim(mo_o))
      dy_o = dy_o - dim(mo_o); mo_o = mo_o + 1
      if (mo_o > 12) then
        mo_o = 1; yr_o = yr_o + 1
        dim = [31,28,31,30,31,30,31,31,30,31,30,31]
        if (is_leap_year(yr_o)) dim(2) = 29
      end if
    end do
  end subroutine datetime_add_seconds

  ! ── Formatadores de data/hora ─────────────────────────────────────────────

  function datetime_to_fname(yr,mo,dy,hr,mn,sc) result(s)
    integer, intent(in) :: yr,mo,dy,hr,mn,sc; character(len=36) :: s
    write(s,'(A,I4.4,2I2.2,A,3I2.2,A)') 'monan_export_',yr,mo,dy,'_',hr,mn,sc,'.nc'
  end function datetime_to_fname

  function datetime_to_iso(yr,mo,dy,hr,mn,sc) result(s)
    integer, intent(in) :: yr,mo,dy,hr,mn,sc; character(len=19) :: s
    write(s,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
      yr,'-',mo,'-',dy,'T',hr,':',mn,':',sc
  end function datetime_to_iso

  function datetime_to_cf_base(yr,mo,dy,hr,mn,sc) result(s)
    integer, intent(in) :: yr,mo,dy,hr,mn,sc; character(len=19) :: s
    write(s,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
      yr,'-',mo,'-',dy,' ',hr,':',mn,':',sc
  end function datetime_to_cf_base

end module mpas_cap_netcdf_mod
