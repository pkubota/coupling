!> @file mpas_cap_config.F90
!! @brief Leitura de namelists de configuracao para o sistema de acoplamento
!!        NUOPC-MPAS-Integrado (MONAN-A 2.0 / ESMF 8.9.1).
!!
!! Modulo central de configuracao: le o arquivo nuopc.input e expoe
!! as variaveis de configuracao via acessores publicos.
!!
!! Arquivo de entrada padrao: nuopc.input (no diretorio de execucao).
!! Caminho alternativo: variavel de ambiente NUOPC_INPUT.
!!
!! Estrutura do arquivo nuopc.input:
!!
!!   &nuopc_esm
!!     use_mpas_atm = .false.     ! .true.  = MPAS-ATM como fonte atmosferica primaria
!!                                ! .false. = DATM/JRA55 (fallback/desenvolvimento)
!!     use_mom6_ocn = .true.      ! .true.  = MOM6+SIS2 integrado (mom_cap.F90)
!!                                ! .false. = DOCN dados prescritos (ocn_comp_NUOPC.F90)
!!   /
!!
!!   &nuopc_driver
!!     start_date   = '2026-03-29'
!!     stop_date    = '2026-03-30'
!!     dt_coupling  = 10800
!!     dt_atm       = 60
!!     log_dir      = 'logs'
!!   /
!!
!!   &nuopc_atm
!!     mesh_atm     = 'mpas_mesh.nc'
!!     config_dir   = './'
!!     write_diag   = .false.
!!   /
!!
!!   &nuopc_netcdf
!!     write_netcdf = .true.
!!     output_dir   = 'diag_export'
!!     grid_res_deg = 1.0
!!   /
!!
!!   &nuopc_atm_bnd
!!     sst_default          = 298.0
!!     ice_fraction_default = 0.0
!!     zorl_default         = 0.01
!!   /
!!
!! Versao 2.0 -- GT Acoplamento MONAN / INPE/CGCT/DIMNT -- Maio 2026.
!!   v1.0: versao inicial (Abril 2026).
!!   v2.0: adicionado &nuopc_esm com use_mpas_atm e use_mom6_ocn (runtime).
!!         Estas flags eram parametros compile-time em esm.F90; agora sao
!!         lidas do nuopc.input antes de ESMF_Initialize e passadas ao
!!         driver via cfg_use_mpas_atm / cfg_use_mom6_ocn.

module mpas_cap_config_mod

  implicit none
  private

  ! -- Constantes ------------------------------------------------------------
  character(len=*), parameter, public :: CONFIG_FILE_DEFAULT = 'nuopc.input'
  integer, parameter :: UNITN = 42   ! unidade Fortran para leitura de namelist

  ! ==========================================================================
  ! Grupo &nuopc_esm -- chaves de configuracao do driver ESM (runtime)
  !
  ! Defaults conservadores: DATM + MOM6 (comportamento da v1.0).
  ! Altere APENAS via nuopc.input -- nao altere os defaults aqui para
  ! evitar recompilar ao mudar de configuracao entre experimentos.
  ! ==========================================================================
  logical, public, protected :: cfg_use_mpas_atm = .false.
  logical, public, protected :: cfg_use_mom6_ocn = .true.

  ! ==========================================================================
  ! Grupo &nuopc_driver -- parametros do driver NUOPC
  ! ==========================================================================
  character(len=10), public, protected :: cfg_start_date   = '2026-03-29'
  character(len=10), public, protected :: cfg_stop_date    = '2026-03-30'
  integer,           public, protected :: cfg_dt_coupling  = 10800  ! [s]
  integer,           public, protected :: cfg_dt_atm       = 60     ! [s]
  character(len=256),public, protected :: cfg_log_dir      = 'logs'

  ! ==========================================================================
  ! Grupo &nuopc_atm -- parametros do cap atmosferico MPAS
  ! ==========================================================================
  character(len=256),public, protected :: cfg_mesh_atm     = 'mpas_mesh.nc'
  character(len=256),public, protected :: cfg_config_dir   = './'
  logical,           public, protected :: cfg_write_diag   = .false.

  ! ==========================================================================
  ! Grupo &nuopc_netcdf -- parametros de saida NetCDF
  ! ==========================================================================
  logical,           public, protected :: cfg_write_netcdf = .true.
  character(len=256),public, protected :: cfg_output_dir   = 'diag_export'
  real,              public, protected :: cfg_grid_res_deg = 1.0    ! [deg]

  ! ==========================================================================
  ! Grupo &nuopc_atm_bnd -- defaults de condicao de contorno
  ! ==========================================================================
  real,              public, protected :: cfg_sst_default          = 298.0  ! [K]
  real,              public, protected :: cfg_ice_fraction_default = 0.0    ! [0-1]
  real,              public, protected :: cfg_zorl_default         = 0.01   ! [m]
  ! B-OCN-01: defaults para correntes oceanicas e salinidade
  real,              public, protected :: cfg_uocn_default         = 0.0    ! [m/s]
  real,              public, protected :: cfg_vocn_default         = 0.0    ! [m/s]
  real,              public, protected :: cfg_sss_default          = 35.0   ! [psu]

  ! ==========================================================================
  ! Grupo &nuopc_datm -- parametros do componente DATM (JRA55)
  !
  ! datm_epoch_date: data de inicio do arquivo JRA55 no formato 'YYYY-MM-DD'.
  !   Define o t=0 da interpolacao temporal linear entre snapshots.
  !   DEVE coincidir com o primeiro snapshot do arquivo NetCDF JRA55.
  !   Erro comum: usar a data de inicio da simulacao (cfg_start_date) em vez
  !   da data de inicio do arquivo -- isso causa tidx0 errado ou negativo.
  !
  ! datm_dt_data: passo temporal dos snapshots JRA55 [s].
  !   JRA55 horario = 3600; JRA55 3-horario (padrao) = 10800.
  ! ==========================================================================
  integer,           public, protected :: cfg_datm_epoch_year  = 2016
  integer,           public, protected :: cfg_datm_epoch_month = 1
  integer,           public, protected :: cfg_datm_epoch_day   = 1
  integer,           public, protected :: cfg_datm_dt_data     = 10800  ! [s] JRA55 3h
  character(len=16), public, protected :: cfg_datocn_mode        = "stub"
  integer,           public, protected :: cfg_datocn_nx          = 1440
  integer,           public, protected :: cfg_datocn_ny          = 720
  integer,           public, protected :: cfg_datocn_dt_data     = 86400
  integer,           public, protected :: cfg_datocn_epoch_year  = 1981
  integer,           public, protected :: cfg_datocn_epoch_month = 9
  integer,           public, protected :: cfg_datocn_epoch_day   = 1
  character(len=256),public, protected :: cfg_datocn_sst_file    = "INPUT/OISST_sst.nc"
  character(len=256),public, protected :: cfg_datocn_ice_file    = "INPUT/OISST_ice.nc"
  character(len=256),public, protected :: cfg_datocn_cur_file    = ""
  ! B-56: nomes das variaveis nos arquivos NetCDF (dependem do produto).
  !   OISST v2.1: sst_varname='sst'  ice_varname='icec'  cur='uo'/'vo'
  !   ERSSTv5   : sst_varname='sst'  ice_varname='sic'
  character(len=64), public, protected :: cfg_datocn_sst_varname  = "sst"
  character(len=64), public, protected :: cfg_datocn_ice_varname  = "icec"
  character(len=64), public, protected :: cfg_datocn_cur_u_varname= "uo"
  character(len=64), public, protected :: cfg_datocn_cur_v_varname= "vo"
  ! Unidade do campo de gelo: .false. = fracao [0,1]; .true. = percentagem [0,100].
  !   Verificar: ncdump -h ice_file.nc | grep "units\|scale_factor"
  logical,           public, protected :: cfg_datocn_ice_pct      = .false.
  ! Diagnostico de importacao: escreve NetCDF por passo
  logical,           public, protected :: cfg_write_import_diag   = .false.
  character(len=256),public, protected :: cfg_import_diag_dir     = "diag_import"

  ! -- API publica -----------------------------------------------------------
  public :: config_read
  public :: config_print
  public :: config_parse_date

contains

  !> @brief Le o arquivo nuopc.input e popula as variaveis de configuracao.
  !!
  !! Sequencia de busca do arquivo:
  !!   1. Argumento opcional file_path
  !!   2. Variavel de ambiente NUOPC_INPUT
  !!   3. nuopc.input no diretorio de execucao (CONFIG_FILE_DEFAULT)
  !!
  !! NOTA B-29: esta subrotina e chamada ANTES de ESMF_Initialize (que faz
  !! MPI_Init). NAO usar MPI aqui. O print de confirmacao foi relocado para
  !! esmApp.F90 apos ESMF_VMGet, onde localPet ja esta disponivel.
  !!
  !! @param[out] rc   0 = sucesso; 1 = arquivo nao encontrado; 2 = erro.
  !! @param[in]  file_path  Caminho alternativo (opcional).
  subroutine config_read(rc, file_path)
    integer,          intent(out)          :: rc
    character(len=*), intent(in), optional :: file_path

    ! -- Variaveis locais espelhando as de modulo --------------------------
    ! &nuopc_esm
    logical            :: use_mpas_atm, use_mom6_ocn
    ! &nuopc_driver
    character(len=10)  :: start_date, stop_date
    integer            :: dt_coupling, dt_atm
    character(len=256) :: log_dir
    ! &nuopc_atm
    character(len=256) :: mesh_atm, config_dir
    logical            :: write_diag
    ! &nuopc_netcdf
    logical            :: write_netcdf
    character(len=256) :: output_dir
    real               :: grid_res_deg
    ! &nuopc_atm_bnd
    real               :: sst_default, ice_fraction_default, zorl_default
    real               :: uocn_default, vocn_default, sss_default
    ! &nuopc_datm
    integer            :: datm_epoch_year, datm_epoch_month, datm_epoch_day
    integer            :: datm_dt_data
    ! &nuopc_datocn
    character(len=16)  :: datocn_mode
    integer            :: datocn_nx, datocn_ny, datocn_dt_data
    integer            :: datocn_epoch_year, datocn_epoch_month, datocn_epoch_day
    character(len=256) :: datocn_sst_file, datocn_ice_file, datocn_cur_file
    character(len=64)  :: datocn_sst_varname, datocn_ice_varname
    character(len=64)  :: datocn_cur_u_varname, datocn_cur_v_varname
    logical            :: datocn_ice_pct
    logical            :: write_import_diag
    character(len=256) :: import_diag_dir

    ! -- Declaracoes de namelists ------------------------------------------
    namelist /nuopc_esm/ use_mpas_atm, use_mom6_ocn

    namelist /nuopc_driver/ start_date, stop_date, dt_coupling, dt_atm, log_dir

    namelist /nuopc_atm/    mesh_atm, config_dir, write_diag

    namelist /nuopc_netcdf/ write_netcdf, output_dir, grid_res_deg

    namelist /nuopc_atm_bnd/ sst_default, ice_fraction_default, zorl_default, &
                              uocn_default, vocn_default, sss_default

    namelist /nuopc_datm/ datm_epoch_year, datm_epoch_month, datm_epoch_day, &
                           datm_dt_data

    namelist /nuopc_datocn/ datocn_mode, datocn_nx, datocn_ny, datocn_dt_data, &
                             datocn_epoch_year, datocn_epoch_month, datocn_epoch_day, &
                             datocn_sst_file, datocn_ice_file, datocn_cur_file, &
                             datocn_sst_varname, datocn_ice_varname, &
                             datocn_cur_u_varname, datocn_cur_v_varname, &
                             datocn_ice_pct, write_import_diag, import_diag_dir

    character(len=512) :: fpath
    logical :: file_exists
    integer :: ios

    rc = 0

    ! -- 1. Inicializar locais com os defaults do modulo ------------------
    use_mpas_atm         = cfg_use_mpas_atm
    use_mom6_ocn         = cfg_use_mom6_ocn
    start_date           = cfg_start_date
    stop_date            = cfg_stop_date
    dt_coupling          = cfg_dt_coupling
    dt_atm               = cfg_dt_atm
    log_dir              = cfg_log_dir
    mesh_atm             = cfg_mesh_atm
    config_dir           = cfg_config_dir
    write_netcdf         = cfg_write_netcdf
    write_diag           = cfg_write_diag
    output_dir           = cfg_output_dir
    grid_res_deg         = cfg_grid_res_deg
    sst_default          = cfg_sst_default
    ice_fraction_default = cfg_ice_fraction_default
    zorl_default         = cfg_zorl_default
    uocn_default         = cfg_uocn_default
    vocn_default         = cfg_vocn_default
    sss_default          = cfg_sss_default
    datm_epoch_year      = cfg_datm_epoch_year
    datm_epoch_month     = cfg_datm_epoch_month
    datm_epoch_day       = cfg_datm_epoch_day
    datm_dt_data         = cfg_datm_dt_data
    datocn_mode          = cfg_datocn_mode
    datocn_nx            = cfg_datocn_nx
    datocn_ny            = cfg_datocn_ny
    datocn_dt_data       = cfg_datocn_dt_data
    datocn_epoch_year    = cfg_datocn_epoch_year
    datocn_epoch_month   = cfg_datocn_epoch_month
    datocn_epoch_day     = cfg_datocn_epoch_day
    datocn_sst_file      = cfg_datocn_sst_file
    datocn_ice_file      = cfg_datocn_ice_file
    datocn_cur_file      = cfg_datocn_cur_file
    datocn_sst_varname   = cfg_datocn_sst_varname
    datocn_ice_varname   = cfg_datocn_ice_varname
    datocn_cur_u_varname = cfg_datocn_cur_u_varname
    datocn_cur_v_varname = cfg_datocn_cur_v_varname
    datocn_ice_pct       = cfg_datocn_ice_pct
    write_import_diag    = cfg_write_import_diag
    import_diag_dir      = cfg_import_diag_dir

    ! -- 2. Determinar caminho do arquivo ---------------------------------
    if (present(file_path) .and. len_trim(file_path) > 0) then
      fpath = trim(file_path)
    else
      call get_environment_variable('NUOPC_INPUT', fpath, status=ios)
      if (ios /= 0 .or. len_trim(fpath) == 0) fpath = CONFIG_FILE_DEFAULT
    end if

    inquire(file=trim(fpath), exist=file_exists)
    if (.not. file_exists) then
      write(*,'(A,A,A)') '[mpas_cap_config] AVISO: arquivo "', &
        trim(fpath), '" nao encontrado -- usando defaults.'
      rc = 1
      return
    end if

    ! -- 3. Abrir e ler cada grupo ----------------------------------------
    open(unit=UNITN, file=trim(fpath), status='old', action='read', &
         form='formatted', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') '[mpas_cap_config] ERRO: nao foi possivel abrir ', trim(fpath)
      rc = 2; return
    end if

    ! &nuopc_esm -- lido PRIMEIRO pois controla quais componentes registrar
    rewind(UNITN)
    read(UNITN, nml=nuopc_esm, iostat=ios)
    if (ios /= 0) write(*,'(A)') &
      '[mpas_cap_config] AVISO: &nuopc_esm ausente -- ' // &
      'use_mpas_atm=F use_mom6_ocn=T (defaults v1.0).'

    ! &nuopc_driver
    rewind(UNITN)
    read(UNITN, nml=nuopc_driver, iostat=ios)
    if (ios /= 0) write(*,'(A)') &
      '[mpas_cap_config] AVISO: &nuopc_driver ausente -- usando defaults.'

    ! &nuopc_atm
    rewind(UNITN)
    read(UNITN, nml=nuopc_atm, iostat=ios)
    if (ios /= 0) write(*,'(A)') &
      '[mpas_cap_config] AVISO: &nuopc_atm ausente -- usando defaults.'

    ! &nuopc_netcdf
    rewind(UNITN)
    read(UNITN, nml=nuopc_netcdf, iostat=ios)
    if (ios /= 0) write(*,'(A)') &
      '[mpas_cap_config] AVISO: &nuopc_netcdf ausente -- usando defaults.'

    ! &nuopc_atm_bnd
    rewind(UNITN)
    read(UNITN, nml=nuopc_atm_bnd, iostat=ios)
    if (ios /= 0) write(*,'(A)') &
      '[mpas_cap_config] AVISO: &nuopc_atm_bnd ausente -- usando defaults.'

    ! &nuopc_datm -- epoch e passo temporal do JRA55
    rewind(UNITN)
    read(UNITN, nml=nuopc_datm, iostat=ios)
    if (ios /= 0) write(*,'(A)') &
      '[mpas_cap_config] AVISO: &nuopc_datm ausente -- ' // &
      'epoch=2016-01-01, dt_data=10800s (defaults JRA55 3h).'

    ! &nuopc_datocn -- ANTES do close para evitar fort.42
    rewind(UNITN)
    read(UNITN, nml=nuopc_datocn, iostat=ios)
    if (ios /= 0) write(*,'(A)') &
      '[mpas_cap_config] AVISO: &nuopc_datocn ausente -- modo stub ativo.'

    close(UNITN)

    ! -- 4. Validacao basica ----------------------------------------------
    if (dt_coupling <= 0) then
      write(*,'(A,I0)') '[mpas_cap_config] ERRO: dt_coupling invalido: ', dt_coupling
      rc = 2; return
    end if
    if (dt_atm <= 0 .or. dt_atm > dt_coupling) then
      write(*,'(A,I0,A,I0)') '[mpas_cap_config] AVISO: dt_atm=', dt_atm, &
        ' deve ser <= dt_coupling=', dt_coupling
    end if
    if (mod(dt_coupling, dt_atm) /= 0) then
      write(*,'(A)') '[mpas_cap_config] AVISO: dt_coupling nao e multiplo de dt_atm.'
    end if
    if (grid_res_deg <= 0.0 .or. grid_res_deg > 10.0) then
      write(*,'(A)') '[mpas_cap_config] AVISO: grid_res_deg fora do intervalo (0,10]. Usando 1.0.'
      grid_res_deg = 1.0
    end if
    if (sst_default < 150.0 .or. sst_default > 350.0) then
      write(*,'(A,F7.2)') '[mpas_cap_config] AVISO: sst_default fora do intervalo fisico: ', &
        sst_default
    end if

    ! -- 5. Copiar para variaveis de modulo -------------------------------
    cfg_use_mpas_atm         = use_mpas_atm
    cfg_use_mom6_ocn         = use_mom6_ocn
    cfg_start_date           = start_date
    cfg_stop_date            = stop_date
    cfg_dt_coupling          = dt_coupling
    cfg_dt_atm               = dt_atm
    cfg_log_dir              = trim(log_dir)
    cfg_mesh_atm             = trim(mesh_atm)
    cfg_config_dir           = trim(config_dir)
    cfg_write_diag           = write_diag
    cfg_write_netcdf         = write_netcdf
    cfg_output_dir           = trim(output_dir)
    cfg_grid_res_deg         = grid_res_deg
    cfg_sst_default          = sst_default
    cfg_ice_fraction_default = ice_fraction_default
    cfg_zorl_default         = zorl_default
    cfg_uocn_default         = uocn_default
    cfg_vocn_default         = vocn_default
    cfg_sss_default          = sss_default
    cfg_datm_epoch_year      = datm_epoch_year
    cfg_datm_epoch_month     = datm_epoch_month
    cfg_datm_epoch_day       = datm_epoch_day
    cfg_datm_dt_data         = datm_dt_data
    cfg_datocn_mode          = trim(datocn_mode)
    cfg_datocn_nx            = datocn_nx
    cfg_datocn_ny            = datocn_ny
    cfg_datocn_dt_data       = datocn_dt_data
    cfg_datocn_epoch_year    = datocn_epoch_year
    cfg_datocn_epoch_month   = datocn_epoch_month
    cfg_datocn_epoch_day     = datocn_epoch_day
    cfg_datocn_sst_file      = trim(datocn_sst_file)
    cfg_datocn_ice_file      = trim(datocn_ice_file)
    cfg_datocn_cur_file      = trim(datocn_cur_file)
    cfg_datocn_sst_varname   = trim(datocn_sst_varname)
    cfg_datocn_ice_varname   = trim(datocn_ice_varname)
    cfg_datocn_cur_u_varname = trim(datocn_cur_u_varname)
    cfg_datocn_cur_v_varname = trim(datocn_cur_v_varname)
    cfg_datocn_ice_pct       = datocn_ice_pct
    cfg_write_import_diag    = write_import_diag
    cfg_import_diag_dir      = trim(import_diag_dir)

  end subroutine config_read

  !> @brief Imprime o resumo da configuracao ativa no stdout.
  subroutine config_print()
    write(*,'(/,A)') '=============================================='
    write(*,'(A)')   ' Configuracao NUOPC-MPAS-Integrado (nuopc.input)'
    write(*,'(A)')   '=============================================='
    write(*,'(A)')   '  [nuopc_esm]'
    write(*,'(2X,A,L1)') '  use_mpas_atm      = ', cfg_use_mpas_atm
    write(*,'(2X,A,L1)') '  use_mom6_ocn      = ', cfg_use_mom6_ocn
    write(*,'(A)') ''
    write(*,'(A)')   '  [nuopc_driver]'
    write(*,'(2X,A,A)')  '  start_date        = ', trim(cfg_start_date)
    write(*,'(2X,A,A)')  '  stop_date         = ', trim(cfg_stop_date)
    write(*,'(2X,A,I0)') '  dt_coupling       = ', cfg_dt_coupling
    write(*,'(2X,A,I0)') '  dt_atm            = ', cfg_dt_atm
    write(*,'(2X,A,A)')  '  log_dir           = ', trim(cfg_log_dir)
    write(*,'(A)') ''
    write(*,'(A)')   '  [nuopc_atm]'
    write(*,'(2X,A,A)')  '  mesh_atm          = ', trim(cfg_mesh_atm)
    write(*,'(2X,A,A)')  '  config_dir        = ', trim(cfg_config_dir)
    write(*,'(2X,A,L1)') '  write_diag        = ', cfg_write_diag
    write(*,'(A)') ''
    write(*,'(A)')   '  [nuopc_netcdf]'
    write(*,'(2X,A,L1)') '  write_netcdf      = ', cfg_write_netcdf
    write(*,'(2X,A,A)')    '  output_dir        = ', trim(cfg_output_dir)
    write(*,'(2X,A,F5.2)') '  grid_res_deg      = ', cfg_grid_res_deg
    write(*,'(A)') ''
    write(*,'(A)')   '  [nuopc_atm_bnd]'
    write(*,'(2X,A,F7.2)') '  sst_default       = ', cfg_sst_default
    write(*,'(2X,A,F5.3)') '  ice_frac_default  = ', cfg_ice_fraction_default
    write(*,'(2X,A,F6.4)') '  zorl_default      = ', cfg_zorl_default
    write(*,'(2X,A,F6.3)') '  uocn_default      = ', cfg_uocn_default
    write(*,'(2X,A,F6.3)') '  vocn_default      = ', cfg_vocn_default
    write(*,'(2X,A,F5.1)') '  sss_default       = ', cfg_sss_default
    write(*,'(A)') ''
    write(*,'(A)')   '  [nuopc_datm]'
    write(*,'(2X,A,I4.4,A,I2.2,A,I2.2)') '  datm_epoch_date   = ', &
      cfg_datm_epoch_year, '-', cfg_datm_epoch_month, '-', cfg_datm_epoch_day
    write(*,'(2X,A,I0,A)') '  datm_dt_data      = ', cfg_datm_dt_data, ' s'
    write(*,'(A)') ''
    write(*,'(A)')   '  [nuopc_datocn]'
    write(*,'(2X,A,A)')       '  datocn_mode       = ', trim(cfg_datocn_mode)
    write(*,'(2X,A,I5,A,I5)') '  grid nx x ny      = ', cfg_datocn_nx, ' x ', cfg_datocn_ny
    write(*,'(2X,A,A)')  '  sst_file          = ', trim(cfg_datocn_sst_file)
    write(*,'(2X,A,A)')  '  sst_varname       = ', trim(cfg_datocn_sst_varname)
    write(*,'(2X,A,A)')  '  ice_file          = ', trim(cfg_datocn_ice_file)
    write(*,'(2X,A,A)')  '  ice_varname       = ', trim(cfg_datocn_ice_varname)
    write(*,'(2X,A,L1)') '  ice_pct           = ', cfg_datocn_ice_pct
    write(*,'(2X,A,A)')  '  cur_file          = ', trim(cfg_datocn_cur_file)
    write(*,'(2X,A,L1)') '  write_import_diag = ', cfg_write_import_diag
    write(*,'(2X,A,A)')  '  import_diag_dir   = ', trim(cfg_import_diag_dir)
    write(*,'(A,/)') '=============================================='
  end subroutine config_print

  !> @brief Converte string 'YYYY-MM-DD' para componentes inteiros.
  !!
  !! @param[in]  date_str  String no formato 'YYYY-MM-DD'
  !! @param[out] yy, mm, dd  Componentes extraidos
  !! @param[out] rc  0 = sucesso; 1 = formato invalido
  subroutine config_parse_date(date_str, yy, mm, dd, rc)
    character(len=*), intent(in)  :: date_str
    integer,          intent(out) :: yy, mm, dd, rc
    integer :: ios

    rc = 0
    if (len_trim(date_str) < 10) then
      rc = 1; return
    end if

    read(date_str(1:4),  '(I4)', iostat=ios) yy
    if (ios /= 0) then; rc = 1; return; end if
    read(date_str(6:7),  '(I2)', iostat=ios) mm
    if (ios /= 0) then; rc = 1; return; end if
    read(date_str(9:10), '(I2)', iostat=ios) dd
    if (ios /= 0) then; rc = 1; return; end if

    if (mm < 1 .or. mm > 12 .or. dd < 1 .or. dd > 31) rc = 1
  end subroutine config_parse_date

end module mpas_cap_config_mod
