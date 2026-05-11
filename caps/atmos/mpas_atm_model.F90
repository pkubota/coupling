!> @file mpas_atm_model.F90
!! @brief Interface com o modelo atmosferico MPAS-A 8.3 / MONAN-A 2.0.
!!
!! Versao 5.1 - Interface MONAN-A 2.0: init/run/final/resize + buffers de acoplamento.
!!   Com -DMPAS_EXTERNAL_ESMF_LIB, mpas_timekeeping.F usa 'use ESMF' (externo).
!!   mpas_advance_stop_time controla o relogio INTERNO do MONAN-A (g_domain%%clock),
!!   independente do relogio ESMF do driver. Ambos sao necessarios.
!!
!! Sequencia de inicializacao do MONAN-A (confirmada via probe no Jaci):
!!   phase1(external_comm) -> atm_setup_core -> atm_setup_domain -> setup_log ->
!!   setup_namelist -> phase2 -> streamInfo -> define_packages -> setup_packages ->
!!   setup_decompositions -> setup_clock -> bootstrap_phase1 -> stream_mgr_init ->
!!   setup_immutable_streams -> xml_stream_parser -> bootstrap_phase2 -> core_init ->
!!   extracao de ponteiros zero-copy.
!!
!! Assinaturas confirmadas (probe_block_type.bash, Jaci):
!!   core_init     : function(domain, startTimeStamp) result(ierr)  [integer]
!!   core_run      : function(domain) result(ierr)                   [integer]
!!   core_finalize : function(domain) result(ierr)                   [integer]

module mpas_atm_model_mod

  ! Tipos publicos em modulo isolado (sem dependencia ESMF externa)
  use mpas_atm_types_mod, only : MPAS_RKIND,             &
                                  mpas_atm_public_type,   &
                                  mpas_atm_state_type,    &
                                  atm_ocean_boundary_type

  use mpas_kind_types,    only : RKIND
  use mpas_derived_types, only : domain_type, mpas_pool_type, MPAS_LOG_CRIT

  ! Duas fases confirmadas no probe (mpas_framework.F linhas 47/78/165)
  use mpas_timekeeping,   only : mpas_timekeeping_init,  &
                                  mpas_advance_stop_time
  use mpas_framework,     only : mpas_framework_init_phase1,  &
                                  mpas_framework_init_phase2,  &
                                  mpas_framework_finalize

  use mpas_domain_routines, only : mpas_allocate_domain

  use mpas_pool_routines, only : mpas_pool_get_array,    &
                                  mpas_pool_get_dimension, &
                                  mpas_pool_get_subpool,   &
                                  mpas_pool_get_config

  use mpas_bootstrapping, only : mpas_bootstrap_framework_phase1, &
                                  mpas_bootstrap_framework_phase2

  use mpas_stream_inquiry, only : MPAS_stream_inquiry_new_streaminfo

  use mpas_stream_manager, only : MPAS_stream_mgr_init,            &
                                   MPAS_stream_mgr_validate_streams
  use iso_c_binding,        only : c_loc, c_ptr, c_int, c_char

  use mpas_io,             only : MPAS_IO_PNETCDF

  ! mpas_log_write confirmado (probe mpas_log.F linha 480)
  use mpas_log,           only : mpas_log_write, mpas_log_info

  ! atm_setup_core: registra core_init/run/finalize em domain%core (secao 8 probe)
  use atm_core_interface, only : atm_setup_core, atm_setup_domain


  use mpas_cap_config_mod,  only : cfg_sst_default, &
                                    cfg_ice_fraction_default, &
                                    cfg_zorl_default, &
                                    cfg_sss_default, &
                                    cfg_uocn_default, &
                                    cfg_vocn_default

  implicit none
  private

  ! Ponteiros PRIVADOS de modulo
  type(domain_type), pointer, private, save :: g_domain   => null()
  integer,                    private, save :: g_mpi_comm = -1

  ! -- Ponteiros para arrays de pool (leitura em mpas_atm_run) ----------------
  ! Necessarios para computar incrementos de acumulados e stress superficial.
  ! Sao ponteiros para memoria do MPAS pool - NAO devem ser desalocados aqui.
  real(MPAS_RKIND), pointer, private, save :: g_pool_acswdnb(:) => null() ! J/m^2 acumulado
  real(MPAS_RKIND), pointer, private, save :: g_pool_aclwdnb(:) => null() ! J/m^2 acumulado
  real(MPAS_RKIND), pointer, private, save :: g_pool_rainnc(:)  => null() ! mm acumulado (estratiforme)
  real(MPAS_RKIND), pointer, private, save :: g_pool_rainc(:)   => null() ! mm acumulado (convectiva)
  real(MPAS_RKIND), pointer, private, save :: g_pool_ust(:)     => null() ! vel. atrito [m/s]
  ! Sa_lfrac_mpas: buffer e ponteiros para fracao de terra
  real(MPAS_RKIND), dimension(:), allocatable,target, save :: g_lfrac_buf
  real(MPAS_RKIND), pointer, private, save :: g_pool_xland(:)    => null()
  integer,          pointer, private, save :: g_pool_landmask(:) => null()
  real(MPAS_RKIND), pointer, private, save :: g_pool_snownc(:)  => null() ! mm acum. neve estrat.
  real(MPAS_RKIND), pointer, private, save :: g_pool_q2(:)      => null() ! hum. espec. 2m [kg/kg]

  ! -- Valores do passo anterior (para calculo de incrementos) ----------------
  real(MPAS_RKIND), allocatable, private, save :: g_prev_acswdnb(:) ! J/m^2
  real(MPAS_RKIND), allocatable, private, save :: g_prev_aclwdnb(:) ! J/m^2
  real(MPAS_RKIND), allocatable, private, save :: g_prev_precip(:)  ! mm (rainnc+rainc)

  ! -- Buffers de saida em unidades instantaneas (apontados por atm_public) ---
  ! atm_public%swdn_sfc, lwdn_sfc, prec_total, taux_sfc, tauy_sfc
  ! apontam para estes arrays apos mpas_atm_init.
  ! OBRIGAToRIO: atributo TARGET para que ptr => array seja valido em Fortran.
  real(MPAS_RKIND), allocatable, target, private, save :: g_swdn_inst(:)  ! W/m^2
  real(MPAS_RKIND), allocatable, target, private, save :: g_lwdn_inst(:)  ! W/m^2
  real(MPAS_RKIND), allocatable, target, private, save :: g_prec_inst(:)  ! kg/m^2/s
  real(MPAS_RKIND), allocatable, target, private, save :: g_taux_buf(:)       ! N/m^2
  real(MPAS_RKIND), allocatable, target, private, save :: g_tauy_buf(:)       ! N/m^2
  real(MPAS_RKIND), allocatable, target, private, save :: g_q2m_buf(:)        ! kg/kg
  real(MPAS_RKIND), allocatable, target, private, save :: g_prec_rain_buf(:)  ! kg/m^2/s
  real(MPAS_RKIND), allocatable, target, private, save :: g_prec_snow_buf(:)  ! kg/m^2/s
  real(MPAS_RKIND), allocatable,         private, save :: g_prev_snow(:)      ! mm acum.

  ! Densidade do ar Ã  superficie: constante de referencia para calculo de stress.
  ! Fonte: NIST, condicoes padrao (1013 hPa, 15Â°C). Erro < 5% na pratica.
  real(MPAS_RKIND), parameter, private :: RHO_AIR_SFC = 1.2_MPAS_RKIND  ! kg/m^3

  ! Velocidade minima para evitar divisao por zero no calculo de stress
  real(MPAS_RKIND), parameter, private :: VMIN = 0.1_MPAS_RKIND   ! m/s

  public :: mpas_atm_init
  public :: mpas_atm_run
  public :: mpas_atm_final
  public :: mpas_atm_init_sfc
  public :: mpas_atm_resize

contains

  !> @brief Obtem o nome do arquivo de malha do namelist.
  function mesh_filename_for_bootstrap(domain) result(fname)
    use mpas_derived_types, only : domain_type
    use mpas_pool_routines, only : mpas_pool_get_config
    use mpas_kind_types,    only : StrKIND
    type(domain_type), pointer, intent(in) :: domain
    character(len=256) :: fname
    character(len=StrKIND), pointer :: config_input_name => null()
    call mpas_pool_get_config(domain%configs, 'config_input_name', config_input_name)
    if (associated(config_input_name)) then
      fname = trim(config_input_name)
    else
      fname = 'x1.40962.init.nc'  ! fallback
    end if
  end function mesh_filename_for_bootstrap

  ! ============================================================================
  !> @brief Inicializa o MONAN-A 2.0.
  !!
  !! Sequencia baseada no probe e em mpas_subdriver.F linhas 202/257:
  !!
  !!   1. mpas_framework_init_phase1(dminfo, external_comm=mpi_comm)
  !!      Inicializa dmpar (MPI wrapper) com o comunicador da VM ESMF.
  !!
  !!   2. atm_setup_core(domain%core)
  !!      Registra os procedure pointers core_init/core_run/core_finalize.
  !!      (Deve ser chamado entre phase1 e phase2 - confirmado pelo probe secao 8)
  !!
  !!   3. mpas_framework_init_phase2(domain, calendar=...)
  !!      Le namelist.atmosphere, decompoe malha, aloca blocklist e pools,
  !!      inicializa clock MPAS-A com config_start_time do namelist.
  !!
  !!   4. Obtem startTimeStamp do clock via mpas_get_clock_time
  !!
  !!   5. ierr = core_init(domain, startTimeStamp)
  !!      Le init.nc, configura fisica e integrador. Retorna inteiro.
  !!
  !!   6. Extrai nCells do subpool 'mesh' via blocklist%structs
  !!
  !!   7. Liga ponteiros zero-copy via subpools structs%mesh, structs%diag
  !!
  !! @param mpi_comm  Comunicador MPI inteiro (extraido pelo cap da VM ESMF).
  ! ============================================================================
  subroutine mpas_atm_init(atm_public, atm_state, atm_bnd, &
                            dt_seconds, config_dir, mpi_comm, rc)

    type(mpas_atm_public_type),    intent(inout) :: atm_public
    type(mpas_atm_state_type),     intent(inout) :: atm_state
    type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
    integer,          intent(in)  :: dt_seconds
    character(len=*), intent(in)  :: config_dir
    integer,          intent(in)  :: mpi_comm
    integer,          intent(out) :: rc

    type(mpas_pool_type), pointer :: meshPool     => null()
    type(mpas_pool_type), pointer :: diagPool     => null()
    type(mpas_pool_type), pointer :: diagPhysPool => null()

    integer, pointer :: nCells_ptr      => null()
    integer, pointer :: nCellsSolve_ptr => null()   ! B-32: celulas proprias (sem halos)
    integer, pointer :: nVertLev_ptr    => null()
    integer          :: n, nSolve, ierr
    character(len=64) :: startTimeStamp
    character(len=256) :: msg

    rc = 0
    g_mpi_comm           = mpi_comm
    atm_state%mpi_comm   = mpi_comm
    atm_state%dt_seconds = dt_seconds
    atm_state%config_dir = trim(config_dir)

    ! ------------------------------------------------------------------
    ! Sequencia replicada de mpas_subdriver.F (confirmada pelo probe):
    !
    ! Aqui: g_domain â‰¡ domain_ptr; g_domain%core â‰¡ corelist.
    ! ------------------------------------------------------------------
    allocate(g_domain)
    nullify(g_domain%next)       ! linked-list: sem proximo domain

    allocate(g_domain%core)
    nullify(g_domain%core%next)  ! linked-list: sem proximo core

    ! Back-link: core%domainlist aponta para o domain
    g_domain%core%domainlist => g_domain

    ! ------------------------------------------------------------------
    ! 1. mpas_allocate_domain: aloca configs, packages, clock,
    !    streamManager, ioContext e faz nullify(blocklist).
    !    Tambem faz allocate(dom%dminfo) - mas phase1 vai re-alocar.
    ! ------------------------------------------------------------------
    call mpas_allocate_domain(g_domain)

    ! ------------------------------------------------------------------
    ! 2. Inicializa timekeeping do MONAN-A com calendario gregoriano.
    !    mpas_timekeeping_init cria os calendarios ESMF via ESMF_CalendarCreate.
    !    Com -DMPAS_EXTERNAL_ESMF_LIB, usa ESMF real (esmf_base%%this).
    ! ------------------------------------------------------------------
    call mpas_timekeeping_init('gregorian')

    ! ------------------------------------------------------------------
    ! 3. Phase 1: inicializa MPI com o comunicador da VM ESMF.
    ! ------------------------------------------------------------------
    nullify(g_domain%dminfo)
    call mpas_framework_init_phase1(g_domain%dminfo, external_comm=g_mpi_comm)

    ! ------------------------------------------------------------------
    ! 3. Registra procedure pointers do nucleo (APOS phase1, conforme
    !    mpas_subdriver.F). atm_setup_core recebe g_domain%core que e
    !    do tipo core_type - ja alocado acima.
    ! ------------------------------------------------------------------
    call atm_setup_core(g_domain%core)

    ! ------------------------------------------------------------------
    ! 4. atm_setup_domain: registra campos adicionais no domain_type
    !    (nomes de variaveis, streams, etc.) - chamado por mpas_subdriver
    !    apos atm_setup_core e antes de phase2.
    ! ------------------------------------------------------------------
    call atm_setup_domain(g_domain)

    ! ------------------------------------------------------------------
    ! 5. setup_log: inicializa o gerenciador de log do MPAS-A.
    !    DEVE ser chamado apos atm_setup_core (que registra o procedure
    !    pointer setup_log em g_domain%core) e apos phase1 (dminfo pronto).
    !    Qualquer mpas_log_write ANTES deste ponto -> g_domain%logInfo
    !    nao inicializado -> SIGSEGV.
    !
    !    B-35 (fix): removida chamada prematura a mpas_log_write que existia
    !    logo apos atm_setup_domain - era a causa raiz do SIGSEGV observado
    !    em todos os 128 ranks (backtrace: mpas_atm_model.F90:251).
    !    Mensagens de progresso anteriores a este ponto devem usar write(*,...).
    !
    !    Sequencia de mpas_subdriver.F:
    !      ierr = domain_ptr%core%setup_log(domain_ptr%logInfo, domain_ptr)
    ! ------------------------------------------------------------------
    ierr = g_domain%core%setup_log(g_domain%logInfo, g_domain)
    if (ierr /= 0) then
      write(*,'(A)') 'ERRO mpas_atm_init: setup_log falhou'
      rc = ierr; return
    end if

    ! ------------------------------------------------------------------
    ! 6. setup_namelist: le namelist.atmosphere para domain%configs.
    !    CRITICO: phase2 le config_pio_num_iotasks e config_pio_stride
    !    de domain%configs. Sem setup_namelist, configs esta vazio ->
    !    mpas_pool_get_config retorna null pointer -> SIGSEGV em phase2.
    !    mpas_subdriver.F:
    !      ierr = domain_ptr%core%setup_namelist(domain_ptr%configs,
    !               domain_ptr%namelist_filename, domain_ptr%dminfo)
    ! ------------------------------------------------------------------
    g_domain%namelist_filename = trim(atm_state%config_dir) // 'namelist.atmosphere'
    g_domain%streams_filename  = trim(atm_state%config_dir) // 'streams.atmosphere'

    ierr = g_domain%core%setup_namelist(g_domain%configs,         &
                                         g_domain%namelist_filename, &
                                         g_domain%dminfo)
    if (ierr /= 0) then
      call mpas_log_write('ERRO mpas_atm_init: setup_namelist falhou', &
                          messageType=MPAS_LOG_CRIT)
      rc = ierr; return
    end if
    call mpas_log_write('mpas_atm_init: setup_log + setup_namelist concluidos')

    ! ------------------------------------------------------------------
    ! 7. Phase 2: decompoe malha, aloca blocklist e pools,
    !    configura clock com config_start_time do namelist.
    !    Chamada sem calendar= pois setup_namelist ja leu config_calendar_type
    !    para domain%configs; phase2 o le de la quando calendar nao e passado.
    ! ------------------------------------------------------------------
    call mpas_framework_init_phase2(g_domain)
    call mpas_log_write('mpas_atm_init: phase2 concluida (I/O init, timekeeping)')

    ! ------------------------------------------------------------------
    ! 8. streamInfo: informacoes sobre streams (lido do XML).
    ! ------------------------------------------------------------------
    g_domain%streamInfo => MPAS_stream_inquiry_new_streaminfo()
    if (.not. associated(g_domain%streamInfo)) then
      call mpas_log_write('ERRO: streamInfo falhou', messageType=MPAS_LOG_CRIT)
      rc = 1; return
    end if
    if (g_domain%streamInfo%init(g_domain%dminfo%comm, g_domain%streams_filename) /= 0) then
      call mpas_log_write('ERRO: streamInfo%init falhou', messageType=MPAS_LOG_CRIT)
      rc = 1; return
    end if

    ! ------------------------------------------------------------------
    ! 9. define_packages / setup_packages / setup_decompositions / setup_clock
    ! ------------------------------------------------------------------
    ierr = g_domain%core%define_packages(g_domain%packages)
    if (ierr /= 0) then
      call mpas_log_write('ERRO: define_packages falhou', messageType=MPAS_LOG_CRIT)
      rc = ierr; return
    end if

    ierr = g_domain%core%setup_packages(g_domain%configs, g_domain%streamInfo, &
                                         g_domain%packages, g_domain%ioContext)
    if (ierr /= 0) then
      call mpas_log_write('ERRO: setup_packages falhou', messageType=MPAS_LOG_CRIT)
      rc = ierr; return
    end if

    ierr = g_domain%core%setup_decompositions(g_domain%decompositions)
    if (ierr /= 0) then
      call mpas_log_write('ERRO: setup_decompositions falhou', messageType=MPAS_LOG_CRIT)
      rc = ierr; return
    end if

    ierr = g_domain%core%setup_clock(g_domain%clock, g_domain%configs)
    if (ierr /= 0) then
      call mpas_log_write('ERRO: setup_clock falhou', messageType=MPAS_LOG_CRIT)
      rc = ierr; return
    end if
    call mpas_log_write('mpas_atm_init: packages + decomp + clock configurados')

    ! ------------------------------------------------------------------
    ! 10. mpas_bootstrap_framework_phase1: le malha, cria blocos,
    !     distribui dominio. Apos esta chamada, blocklist esta alocado.
    !     O filename do mesh e lido de config_input_name no namelist.
    ! ------------------------------------------------------------------
    call mpas_bootstrap_framework_phase1(g_domain, &
         trim(mesh_filename_for_bootstrap(g_domain)), MPAS_IO_PNETCDF)

    if (.not. associated(g_domain%blocklist)) then
      call mpas_log_write('ERRO: blocklist nulo apos bootstrap_phase1', &
                          messageType=MPAS_LOG_CRIT)
      rc = 1; return
    end if
    call mpas_log_write('mpas_atm_init: bootstrap_phase1 concluido (blocklist alocado)')

    ! ------------------------------------------------------------------
    ! 11. Configura stream manager e streams imutaveis.
    ! ------------------------------------------------------------------
    call MPAS_stream_mgr_init(g_domain%streamManager, g_domain%ioContext, &
                              g_domain%clock, g_domain%blocklist%allFields, &
                              g_domain%packages, g_domain%blocklist%allStructs)

    ierr = g_domain%core%setup_immutable_streams(g_domain%streamManager)
    if (ierr /= 0) then
      call mpas_log_write('ERRO: setup_immutable_streams falhou', messageType=MPAS_LOG_CRIT)
      rc = ierr; return
    end if

    ! ------------------------------------------------------------------
    ! 11b. xml_stream_parser: parseia streams.atmosphere e registra todas
    !      as streams dinamicas no stream manager.
    !      CRiTICO: sem esta chamada, as streams do namelist nao sao
    !      registradas -> reads retornam garbage -> crash na fisica.
    !      Interface C definida localmente (igual ao mpas_subdriver.F).
    ! ------------------------------------------------------------------
    block
      use iso_c_binding, only : c_loc, c_ptr, c_int, c_char
      interface
        subroutine xml_stream_parser(xmlname, mgr_p, comm, ierr) bind(c)
          use iso_c_binding, only : c_char, c_ptr, c_int
          character(kind=c_char), dimension(*), intent(in)    :: xmlname
          type(c_ptr),                          intent(inout) :: mgr_p
          integer(kind=c_int),                  intent(inout) :: comm
          integer(kind=c_int),                  intent(out)   :: ierr
        end subroutine xml_stream_parser
      end interface
      type(c_ptr)                            :: mgr_p
      integer(kind=c_int)                    :: c_comm, c_ierr
      character(kind=c_char,len=1), dimension(512) :: c_filename
      integer :: k, slen

      ! Converter streams_filename para C string
      slen = len_trim(g_domain%streams_filename)
      do k = 1, slen
        c_filename(k) = g_domain%streams_filename(k:k)
      end do
      c_filename(slen+1) = achar(0)  ! null terminator (C string)

#ifdef MPAS_USE_MPI_F08
      c_comm = g_domain%dminfo%comm%mpi_val
#else
      c_comm = g_domain%dminfo%comm
#endif
      mgr_p = c_loc(g_domain%streamManager)
      call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
      if (c_ierr /= 0) then
        call mpas_log_write('ERRO: xml_stream_parser falhou para streams.atmosphere', &
                            messageType=MPAS_LOG_CRIT)
        rc = 1; return
      end if
    end block

    call mpas_log_write('mpas_atm_init: xml_stream_parser concluido')

    ! Valida streams apos configuracao
    call MPAS_stream_mgr_validate_streams(g_domain%streamManager, ierr=ierr)
    if (ierr /= 0) then
      call mpas_log_write('ERRO: stream manager validation falhou', messageType=MPAS_LOG_CRIT)
      rc = 1; return
    end if
    call mpas_log_write('mpas_atm_init: streams validadas')

    ! ------------------------------------------------------------------
    ! 12. mpas_bootstrap_framework_phase2: finaliza alocacao de campos e halos.
    ! ------------------------------------------------------------------
    call mpas_bootstrap_framework_phase2(g_domain)
    call mpas_log_write('mpas_atm_init: bootstrap_phase2 concluido')

    ! ------------------------------------------------------------------
    ! 13. Inicializa o nucleo atmosferico (core_init).
    ! ------------------------------------------------------------------
    startTimeStamp = ''
    ierr = g_domain%core%core_init(g_domain, startTimeStamp)
    if (ierr /= 0) then
      write(msg,'(A,I0)') 'ERRO mpas_atm_init: core_init retornou ierr=', ierr
      call mpas_log_write(trim(msg), messageType=MPAS_LOG_CRIT)
      rc = ierr; return
    end if
    call mpas_log_write('mpas_atm_init: core_init concluido')

    ! ------------------------------------------------------------------
    ! 6. Extrai nCells do subpool 'mesh'
    !    Probe secao 7, mpas_atm_core.F linha 167:
    !    O campo do bloco e 'structs' (confirmado pelo probe).
    ! ------------------------------------------------------------------
    call mpas_pool_get_subpool(g_domain%blocklist%structs, 'mesh', meshPool)

    if (.not. associated(meshPool)) then
      write(*,'(A)') 'ERRO mpas_atm_init: subpool mesh nao encontrado em blocklist%structs'
      rc = 1; return
    end if

    call mpas_pool_get_dimension(meshPool, 'nCells',      nCells_ptr)
    call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve_ptr)  ! B-32
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLev_ptr)

    if (.not. associated(nCells_ptr)) then
      write(*,'(A)') 'ERRO mpas_atm_init: nCells nao encontrado no subpool mesh'
      rc = 1; return
    end if

    n = nCells_ptr

    ! B-33: merge() avalia AMBOS os argumentos (tsource e fsource) antes de
    ! aplicar a mascara - comportamento mandatorio do padrao Fortran (7.1.5.2).
    ! Se nCellsSolve_ptr for null(), a referencia implicita ao ponteiro em tsource
    ! gera SIGSEGV independentemente do valor de mask=associated(...).
    ! Correcao: usar if/else para evitar qualquer dereference quando null.
    if (associated(nCellsSolve_ptr)) then
      nSolve = nCellsSolve_ptr
    else
      nSolve = n
      write(*,'(A)') 'AVISO mpas_atm_init: nCellsSolve ausente no pool mesh - usando nCells'
    end if

    ! B-34: nVertLev_ptr pode ser null() se 'nVertLevels' nao existir no pool
    ! (e.g., nome divergente no Registry.xml de alguma versao). Desreferenciar
    ! um ponteiro null gera SIGSEGV. Guardar com associated() antes de usar.
    if (associated(nVertLev_ptr)) then
      atm_state%nVertLevels  = nVertLev_ptr
      atm_public%nVertLevels = nVertLev_ptr
    else
      write(*,'(A)') 'AVISO mpas_atm_init: nVertLevels ausente no pool mesh - usando default 55'
      atm_state%nVertLevels  = 55   ! default da fisica MONAN-A 2.0
      atm_public%nVertLevels = 55
    end if

    atm_state%nCells       = n
    atm_public%nCells      = n
    atm_public%nCellsSolve = nSolve   ! B-32: expoe para netcdf_init_coords

    ! ------------------------------------------------------------------
    ! 7a. Ponteiros zero-copy: geometria (subpool 'mesh')
    !     Confirmado: mpas_atm_core.F linha 437 usa 'areaCell' de mesh.
    ! ------------------------------------------------------------------
    call mpas_pool_get_array(meshPool, 'latCell',  atm_public%latCell)
    call mpas_pool_get_array(meshPool, 'lonCell',  atm_public%lonCell)
    call mpas_pool_get_array(meshPool, 'areaCell', atm_public%areaCell)

    if (.not. associated(atm_public%latCell)) then
      write(*,'(A)') 'ERRO mpas_atm_init: latCell nao encontrado no subpool mesh'
      rc = 1; return
    end if

    ! ------------------------------------------------------------------
    ! 7b. Ponteiros zero-copy: diagnosticos
    !
    !  No MONAN-A 2.0 os campos estao distribuidos em dois subpools:
    !
    !  subpool 'diag'         - variaveis termodinamicas e de radiacao:
    !    mslp, acswdnb, aclwdnb, rainnc, u10, v10
    !
    !  subpool 'diag_physics' - saidas de pacotes de CLP/superficie
    !    (ativo com bl_mynn_in=T ou bl_ysu_in=T):
    !    t2m, lh, hfx
    !
    !  Estrategia: busca em 'diag' primeiro; para qualquer campo ainda
    !  nulo, tenta 'diag_physics'. Cobre reorganizacoes do Registry.xml
    !  entre versoes sem exigir probe externo.
    !
    !  Confirmado nos logs: mslp encontrado em 'diag'; t2m, acswdnb,
    !  rainnc e lh retornavam nulo em 'diag' com mesoscale_reference_monan.
    ! ------------------------------------------------------------------

    ! Passa 1: subpool 'diag'
    call mpas_pool_get_subpool(g_domain%blocklist%structs, 'diag', diagPool)

    if (associated(diagPool)) then
      call mpas_pool_get_array(diagPool, 'mslp',    atm_public%pslv)     ! PSLV [Pa]
      call mpas_pool_get_array(diagPool, 'u10',     atm_public%u10)      ! U 10m [m/s]
      call mpas_pool_get_array(diagPool, 'v10',     atm_public%v10)      ! V 10m [m/s]
      ! Ponteiros privados para pools acumulados - nao expostos diretamente
      call mpas_pool_get_array(diagPool, 'acswdnb', g_pool_acswdnb)      ! J/m^2 acum.
      call mpas_pool_get_array(diagPool, 'aclwdnb', g_pool_aclwdnb)      ! J/m^2 acum.
      call mpas_pool_get_array(diagPool, 'rainnc',  g_pool_rainnc)       ! mm acum. (estrat.)
      call mpas_pool_get_array(diagPool, 'rainc',   g_pool_rainc)        ! mm acum. (conv.)
      call mpas_pool_get_array(diagPool, 'snownc',  g_pool_snownc)       ! mm acum. neve estrat.
      ! t2m, lh, hfx: tentativa em 'diag'
      call mpas_pool_get_array(diagPool, 't2m',     atm_public%t2m)
      call mpas_pool_get_array(diagPool, 'lh',      atm_public%lhflx)
      call mpas_pool_get_array(diagPool, 'hfx',     atm_public%shflx)
    else
      write(*,'(A)') 'AVISO mpas_atm_init: subpool diag nao encontrado em structs'
    end if

    ! Passa 2: subpool 'diag_physics' - fallback para campos de CLP/superficie
    ! No MONAN-A 2.0 com suite mesoscale_reference_monan, t2m/lh/hfx/ust estao aqui.
    call mpas_pool_get_subpool(g_domain%blocklist%structs, 'diag_physics', diagPhysPool)

    if (associated(diagPhysPool)) then
      ! Sobrescreve apenas ponteiros ainda nulos apos a busca em 'diag'
      if (.not. associated(atm_public%t2m))    &
        call mpas_pool_get_array(diagPhysPool, 't2m',     atm_public%t2m)
      if (.not. associated(atm_public%u10))    &
        call mpas_pool_get_array(diagPhysPool, 'u10',     atm_public%u10)
      if (.not. associated(atm_public%v10))    &
        call mpas_pool_get_array(diagPhysPool, 'v10',     atm_public%v10)
      if (.not. associated(g_pool_acswdnb))    &
        call mpas_pool_get_array(diagPhysPool, 'acswdnb', g_pool_acswdnb)
      if (.not. associated(g_pool_aclwdnb))    &
        call mpas_pool_get_array(diagPhysPool, 'aclwdnb', g_pool_aclwdnb)
      if (.not. associated(g_pool_rainnc))     &
        call mpas_pool_get_array(diagPhysPool, 'rainnc',  g_pool_rainnc)
      if (.not. associated(g_pool_rainc))      &
        call mpas_pool_get_array(diagPhysPool, 'rainc',   g_pool_rainc)
      if (.not. associated(g_pool_snownc))     &
        call mpas_pool_get_array(diagPhysPool, 'snownc',  g_pool_snownc)
      if (.not. associated(g_pool_q2))         &
        call mpas_pool_get_array(diagPhysPool, 'q2',      g_pool_q2)
      if (.not. associated(atm_public%lhflx))  &
        call mpas_pool_get_array(diagPhysPool, 'lh',      atm_public%lhflx)
      if (.not. associated(atm_public%shflx))  &
        call mpas_pool_get_array(diagPhysPool, 'hfx',     atm_public%shflx)
      ! Velocidade de atrito - necessaria para calcular stress superficial
      call mpas_pool_get_array(diagPhysPool, 'ust', g_pool_ust)
    end if

    call warn_if_null(atm_public%t2m,      't2m')
    call warn_if_null(atm_public%pslv,     'mslp')
    call warn_if_null(g_pool_acswdnb,      'acswdnb')
    call warn_if_null(g_pool_rainnc,       'rainnc')
    call warn_if_null(atm_public%lhflx,    'lh')
    if (.not. associated(g_pool_ust)) &
      write(*,'(A)') 'AVISO mpas_atm_init: ust nulo - taux/tauy serao zero'

    ! ------------------------------------------------------------------
    ! 7c. Alocar buffers de saida em unidades instantaneas e apontar
    !     atm_public para eles.
    !
    !  Os campos acumulados do MPAS (acswdnb, aclwdnb, rainnc, rainc)
    !  NaO podem ser expostos diretamente como Faxa_swdn/lwdn/prec porque:
    !    1. Sao acumulados desde t=0 - nao representam o intervalo de acoplamento.
    !    2. Dividir pelo tempo total (Ã· elapsed_s) da a media desde t=0, nao
    !       a media do ultimo intervalo - divergencia crescente ao longo do dia.
    !
    !  Solucao: em cada mpas_atm_run, computar:
    !    swdn_inst = (acswdnb_N - acswdnb_{N-1}) / dt_coupling  [W/m^2]
    !    lwdn_inst = (aclwdnb_N - aclwdnb_{N-1}) / dt_coupling  [W/m^2]
    !    prec_inst = (rainnc_N + rainc_N - prev_N) / dt / 1000  [kg/m^2/s]
    !
    !  Analogamente, taux/tauy nunca foram populados em nenhuma passada pelos
    !  pools - permanecem nulos -> Faxa_taux/tauy nunca sao exportados. Fix:
    !    taux = rho· ust^2 x u10 / max(|V10|, VMIN)  [N/m^2]
    !    tauy = rho· ust^2 x v10 / max(|V10|, VMIN)  [N/m^2]
    ! ------------------------------------------------------------------
    allocate(g_prev_acswdnb(n), g_prev_aclwdnb(n), g_prev_precip(n))
    allocate(g_swdn_inst(n), g_lwdn_inst(n), g_prec_inst(n))
    allocate(g_taux_buf(n), g_tauy_buf(n))
    allocate(g_q2m_buf(n), g_prec_rain_buf(n), g_prec_snow_buf(n))
    allocate(g_prev_snow(n))

    ! Sa_lfrac_mpas: derivar fracao de terra para exportar ao MED_cap
    ! Estrategia: tentar xland (WRF real, 1=terra 2=mar) -> fallback landmask
    allocate(g_lfrac_buf(n))
    g_lfrac_buf = 0.0_MPAS_RKIND   ! default: oceano puro
    call mpas_pool_get_array(diagPool, 'xland', g_pool_xland)
    if (associated(g_pool_xland)) then
      ! xland WRF: 1.0=terra, 2.0=agua -> lfrac: 1.0=terra, 0.0=oceano
      g_lfrac_buf(1:n) = max(0.0_MPAS_RKIND, &
                            2.0_MPAS_RKIND - g_pool_xland(1:n))
      call mpas_log_write('mpas_atm_init: lfrac <- xland (WRF)')
    else
      ! fallback: landmask do pool mesh (1=terra, 0=mar)
      call mpas_pool_get_array(meshPool, 'landmask', g_pool_landmask)
      if (associated(g_pool_landmask)) then
        g_lfrac_buf(1:n) = real(g_pool_landmask(1:n), MPAS_RKIND)
        call mpas_log_write('mpas_atm_init: lfrac <- landmask (mesh)')
      else
        call mpas_log_write('mpas_atm_init: AVISO - xland e landmask ' // &
          'ausentes; Sa_lfrac_mpas = 0 (oceano puro)')
      end if
    end if
    atm_public%lfrac => g_lfrac_buf

    ! Inicializar valores do passo anterior com estado t=0 (apos core_init)
    if (associated(g_pool_acswdnb)) then
      g_prev_acswdnb = g_pool_acswdnb(1:n)
    else
      g_prev_acswdnb = 0.0_MPAS_RKIND
    end if
    if (associated(g_pool_aclwdnb)) then
      g_prev_aclwdnb = g_pool_aclwdnb(1:n)
    else
      g_prev_aclwdnb = 0.0_MPAS_RKIND
    end if
    ! Precip total t=0: rainnc + rainc (podem ser nao-zero apos hot-start)
    if (associated(g_pool_rainnc) .and. associated(g_pool_rainc)) then
      g_prev_precip = g_pool_rainnc(1:n) + g_pool_rainc(1:n)
    else if (associated(g_pool_rainnc)) then
      g_prev_precip = g_pool_rainnc(1:n)
    else
      g_prev_precip = 0.0_MPAS_RKIND
    end if
    ! Neve acumulada t=0
    if (associated(g_pool_snownc)) then
      g_prev_snow = g_pool_snownc(1:n)
    else
      g_prev_snow = 0.0_MPAS_RKIND
    end if

    ! Buffers inicializados a zero (serao preenchidos no primeiro core_run)
    g_swdn_inst     = 0.0_MPAS_RKIND
    g_lwdn_inst     = 0.0_MPAS_RKIND
    g_prec_inst     = 0.0_MPAS_RKIND
    g_taux_buf      = 0.0_MPAS_RKIND
    g_tauy_buf      = 0.0_MPAS_RKIND
    g_q2m_buf       = 0.0_MPAS_RKIND
    g_prec_rain_buf = 0.0_MPAS_RKIND
    g_prec_snow_buf = 0.0_MPAS_RKIND

    ! Redirecionar atm_public para buffers computados (em vez de pool diretamente)
    atm_public%swdn_sfc   => g_swdn_inst
    atm_public%lwdn_sfc   => g_lwdn_inst
    atm_public%prec_total => g_prec_inst
    atm_public%taux_sfc   => g_taux_buf
    atm_public%tauy_sfc   => g_tauy_buf
    atm_public%q2m        => g_q2m_buf
    atm_public%prec_rain  => g_prec_rain_buf
    atm_public%prec_snow  => g_prec_snow_buf

    ! ------------------------------------------------------------------
    ! 8. Aloca arrays de propriedade deste modulo
    ! ------------------------------------------------------------------
    allocate(atm_bnd%sst(n), atm_bnd%ice_fraction(n), atm_bnd%zorl(n))
    allocate(atm_bnd%uocn(n), atm_bnd%vocn(n), atm_bnd%sss(n))
    atm_bnd%sst          = real(cfg_sst_default,          MPAS_RKIND)
    atm_bnd%ice_fraction = real(cfg_ice_fraction_default, MPAS_RKIND)
    atm_bnd%zorl         = real(cfg_zorl_default,         MPAS_RKIND)
    atm_bnd%uocn         = real(cfg_uocn_default,         MPAS_RKIND)
    atm_bnd%vocn         = real(cfg_vocn_default,         MPAS_RKIND)
    atm_bnd%sss          = real(cfg_sss_default,          MPAS_RKIND)

    atm_state%initialized = .true.
    ! B-32: nSolve = celulas proprias (sem halos); n = nCells total (com halos).
    ! netcdf_init_coords deve usar nSolve -> soma global = 40962.
    write(msg,'(A,I0,A,I0,A)') &
      'mpas_atm_init: OK (', nSolve, ' celulas proprias / ', n, ' com halos - SMIOL ativo)'
    call mpas_log_write(trim(msg))

  end subroutine mpas_atm_init

  ! ============================================================================
  subroutine mpas_atm_init_sfc(atm_public, atm_state, rc)
    type(mpas_atm_public_type), intent(inout) :: atm_public
    type(mpas_atm_state_type),  intent(inout) :: atm_state
    integer,                    intent(out)   :: rc
    rc = 0
    if (.not. atm_state%initialized) then
      write(*,'(A)') 'ERRO mpas_atm_init_sfc: modelo nao inicializado'
      rc = 1; return
    end if
    ! core_init ja preencheu o subpool diag com dados do init.nc via SMIOL.
    ! Os ponteiros zero-copy em atm_public ja contem dados validos.
    call mpas_log_write('mpas_atm_init_sfc: campos t=0 prontos (zero-copy)')
  end subroutine mpas_atm_init_sfc

  ! ============================================================================
  !> @brief Avanca o MONAN-A por um intervalo de acoplamento.
  !!
  !! Probe secao 4 / mpas_atm_core.F linha 605:
  !!   function atm_core_run(domain) result(ierr)
  !! core_run e INTEGER FUNCTION - retorna codigo de erro MPAS.
  !!
  !! I/O (history/restart) via SMIOL/smiolf ocorre automaticamente
  !! conforme alarmes definidos em streams.atmosphere.
  ! ============================================================================
  subroutine mpas_atm_run(atm_public, atm_state, atm_bnd, dt_coupling, rc)

    type(mpas_atm_public_type),    intent(inout) :: atm_public
    type(mpas_atm_state_type),     intent(inout) :: atm_state
    type(atm_ocean_boundary_type), intent(in)    :: atm_bnd
    integer,                       intent(in)    :: dt_coupling
    integer,                       intent(out)   :: rc

    type(mpas_pool_type), pointer   :: sfcInputPool => null()
    real(MPAS_RKIND), dimension(:), pointer :: sst_field  => null()
    real(MPAS_RKIND), dimension(:), pointer :: ice_field  => null()
    real(MPAS_RKIND), dimension(:), pointer :: zorl_field => null()
    real(MPAS_RKIND), dimension(:), pointer :: uocn_field => null()  !< corrente zonal
    real(MPAS_RKIND), dimension(:), pointer :: vocn_field => null()  !< corrente meridional
    integer :: n, ierr
    character(len=256) :: msg

    rc = 0
    n  = atm_state%nCells

    if (.not. atm_state%initialized .or. .not. associated(g_domain)) then
      write(*,'(A)') 'ERRO mpas_atm_run: modelo nao inicializado'
      rc = 1; return
    end if

    ! ------------------------------------------------------------------
    ! Injeta condicoes de fronteira no subpool 'sfc_input'
    ! Probe secao 7 / mpas_atm_core.F linha 553:
    ! Nomes Registry.xml: sst, iceAreaCell, znt
    ! ------------------------------------------------------------------
    call mpas_pool_get_subpool(g_domain%blocklist%structs, 'sfc_input', sfcInputPool)

    if (associated(sfcInputPool)) then
      call mpas_pool_get_array(sfcInputPool, 'sst',              sst_field)
      call mpas_pool_get_array(sfcInputPool, 'iceAreaCell',      ice_field)
      call mpas_pool_get_array(sfcInputPool, 'znt',              zorl_field)
      ! B-OCN-01: correntes superficiais e salinidade
      ! Nomes Registry.xml MPAS-A 8.x: u_surface_velocity, v_surface_velocity
      ! Disponivel no pool sfc_input quando config_frac_seaice=.true.
      ! e physics_suite inclui bl_mynn_surface ou sfc_sfclay.
      call mpas_pool_get_array(sfcInputPool, 'u_surface_velocity', uocn_field)
      call mpas_pool_get_array(sfcInputPool, 'v_surface_velocity', vocn_field)

      if (associated(sst_field)  .and. allocated(atm_bnd%sst)) &
           sst_field(1:n)  = atm_bnd%sst(1:n)
      if (associated(ice_field)  .and. allocated(atm_bnd%ice_fraction)) &
           ice_field(1:n)  = atm_bnd%ice_fraction(1:n)
      if (associated(zorl_field) .and. allocated(atm_bnd%zorl)) &
           zorl_field(1:n) = atm_bnd%zorl(1:n)
      ! Corrente relativa: injeta apenas se o pool tem o campo
      if (associated(uocn_field) .and. allocated(atm_bnd%uocn)) &
           uocn_field(1:n) = atm_bnd%uocn(1:n)
      if (associated(vocn_field) .and. allocated(atm_bnd%vocn)) &
           vocn_field(1:n) = atm_bnd%vocn(1:n)
    else
      write(*,'(A)') 'AVISO mpas_atm_run: subpool sfc_input nao encontrado em structs'
    end if

    call mpas_log_write('mpas_atm_run: sfc_input injetado')

    ! mpas_advance_stop_time: avanca o stop time do relogio MPAS interno
    ! por exatamente dt_coupling antes de core_run.
    ! Avanca o stop time do relogio interno do MONAN-A (g_domain%clock),
    ! independente do relogio ESMF do driver. Controla quantos passos
    ! internos (dt_atm) core_run integra por chamada a mpas_atm_run.
    call mpas_advance_stop_time(g_domain%clock, dt_coupling)

    ! ------------------------------------------------------------------
    ! Ativa mpas_log_info -> domain%logInfo antes de core_run.
    ! mpas_subdriver.F linha 414:
    ! Sem isso, mpas_log_write dentro de core_run derreferencia null -> SIGSEGV.
    ! ------------------------------------------------------------------
    if (associated(g_domain%logInfo)) mpas_log_info => g_domain%logInfo

    ! ------------------------------------------------------------------
    ! Avanca o nucleo: integra passos internos de dt_atm, escreve I/O
    ! via SMIOL conforme streams.atmosphere.
    ! core_run e INTEGER FUNCTION.
    ! ------------------------------------------------------------------
    ierr = g_domain%core%core_run(g_domain)
    if (ierr /= 0) then
      write(msg,'(A,I0)') 'ERRO mpas_atm_run: core_run retornou ierr=', ierr
      write(*,'(A)') trim(msg)
      call mpas_log_write(trim(msg))
      rc = ierr; return
    end if

    call mpas_log_write('mpas_atm_run: core_run concluido')

    ! ------------------------------------------------------------------
    ! Pos-processamento dos campos acumulados e stress superficial.
    !
    ! Os arrays do pool (g_pool_*) foram atualizados por core_run.
    ! Agora computamos os valores instantaneos para o intervalo de
    ! acoplamento e armazenamos nos buffers g_*_inst / g_taux_buf / g_tauy_buf
    ! que sao apontados por atm_public%swdn_sfc, lwdn_sfc, prec_total,
    ! taux_sfc, tauy_sfc (configurado em mpas_atm_init).
    !
    ! IMPORTANTE: usar real(dt_coupling, MPAS_RKIND) para evitar perda de
    ! precisao quando MPAS_RKIND = kind(1.0) (single precision).
    ! ------------------------------------------------------------------
    block
      real(MPAS_RKIND) :: dt_r, precip_now, spd
      integer          :: k
      dt_r = real(dt_coupling, MPAS_RKIND)

      ! -- SW e LW descendentes: incremento Ã· dt -> W/m^2 -------------
      if (associated(g_pool_acswdnb)) then
        do k = 1, n
          g_swdn_inst(k) = max((g_pool_acswdnb(k) - g_prev_acswdnb(k)) / dt_r, &
                               0.0_MPAS_RKIND)
        end do
        g_prev_acswdnb(1:n) = g_pool_acswdnb(1:n)
      end if

      if (associated(g_pool_aclwdnb)) then
        do k = 1, n
          g_lwdn_inst(k) = max((g_pool_aclwdnb(k) - g_prev_aclwdnb(k)) / dt_r, &
                               0.0_MPAS_RKIND)
        end do
        g_prev_aclwdnb(1:n) = g_pool_aclwdnb(1:n)
      end if

      ! -- Precipitacao total: (rainnc + rainc) incremento Ã· dt ------
      ! rainnc [mm] = precipitacao estratiforme acumulada
      ! rainc  [mm] = precipitacao convectiva acumulada (esquema GF/KF)
      ! 1 mm = 1 kg/m^2 -> taxa = D[mm] / dt [kg/m^2/s]
      do k = 1, n
        precip_now = 0.0_MPAS_RKIND
        if (associated(g_pool_rainnc)) precip_now = precip_now + g_pool_rainnc(k)
        if (associated(g_pool_rainc))  precip_now = precip_now + g_pool_rainc(k)
        g_prec_inst(k) = max((precip_now - g_prev_precip(k)) / dt_r, &
                              0.0_MPAS_RKIND)
      end do
      ! Atualizar acumulado anterior
      do k = 1, n
        g_prev_precip(k) = 0.0_MPAS_RKIND
        if (associated(g_pool_rainnc)) g_prev_precip(k) = g_prev_precip(k) + g_pool_rainnc(k)
        if (associated(g_pool_rainc))  g_prev_precip(k) = g_prev_precip(k) + g_pool_rainc(k)
      end do

      ! -- Precipitacao solida (neve): snownc incremento Ã· dt --------
      ! snownc [mm] = neve estratiforme acumulada (subconjunto de rainnc)
      ! Se snownc nao estiver disponivel, usa particao por temperatura:
      !   T < T_FREEZE -> tudo neve; caso contrario -> tudo chuva
      block
        real(MPAS_RKIND), parameter :: T_FREEZE = 273.15_MPAS_RKIND
        real(MPAS_RKIND) :: snow_now, delta_snow, delta_total
        do k = 1, n
          delta_total = g_prec_inst(k)
          if (associated(g_pool_snownc)) then
            snow_now = g_pool_snownc(k)
            delta_snow = max((snow_now - g_prev_snow(k)) / dt_r, 0.0_MPAS_RKIND)
            g_prec_snow_buf(k) = min(delta_snow, delta_total)
            g_prec_rain_buf(k) = max(delta_total - g_prec_snow_buf(k), 0.0_MPAS_RKIND)
          else if (associated(atm_public%t2m)) then
            ! Fallback: particao por temperatura
            if (atm_public%t2m(k) < T_FREEZE) then
              g_prec_snow_buf(k) = delta_total
              g_prec_rain_buf(k) = 0.0_MPAS_RKIND
            else
              g_prec_snow_buf(k) = 0.0_MPAS_RKIND
              g_prec_rain_buf(k) = delta_total
            end if
          else
            g_prec_rain_buf(k) = delta_total
            g_prec_snow_buf(k) = 0.0_MPAS_RKIND
          end if
        end do
        ! Atualizar acumulado anterior de neve
        if (associated(g_pool_snownc)) then
          g_prev_snow(1:n) = g_pool_snownc(1:n)
        end if
      end block

      ! -- Umidade especifica a 2m: q2 [kg/kg] -----------------------
      ! g_pool_q2 e ponteiro direto para o pool - sem buffer de incremento.
      ! Valor instantaneo -> valido para o instante corrente.
      if (associated(g_pool_q2)) then
        g_q2m_buf(1:n) = g_pool_q2(1:n)
      else if (associated(atm_public%t2m)) then
        ! Fallback: umidade de saturacao em T2m (Tetens) Ã— RH=0.8
        block
          real(MPAS_RKIND) :: es, qs
          real(MPAS_RKIND), parameter :: es0 = 611.2_MPAS_RKIND
          real(MPAS_RKIND), parameter :: a   = 17.67_MPAS_RKIND
          real(MPAS_RKIND), parameter :: b   = 243.5_MPAS_RKIND
          real(MPAS_RKIND), parameter :: eps = 0.622_MPAS_RKIND
          real(MPAS_RKIND), parameter :: p0  = 101325.0_MPAS_RKIND
          do k = 1, n
            es = es0 * exp(a*(atm_public%t2m(k)-273.15_MPAS_RKIND) / &
                           (b + atm_public%t2m(k)-273.15_MPAS_RKIND))
            qs = eps * es / (p0 - es)
            g_q2m_buf(k) = 0.8_MPAS_RKIND * qs   ! RH=80% como fallback
          end do
        end block
      end if

      ! -- Stress superficial: tau = rho· ust^2 x (u,v) / |V10| ---------
      ! Formula de Monin-Obukhov: CD = (ust/|V10|)^2
      ! -> tau_x = rho· CD x |V10| x u10 = rho· ust^2 x u10 / |V10|
      ! Direcao positiva: eastward (taux>0 quando vento vai para leste)
      if (associated(g_pool_ust) .and. &
          associated(atm_public%u10) .and. associated(atm_public%v10)) then
        do k = 1, n
          spd = sqrt(atm_public%u10(k)**2 + atm_public%v10(k)**2)
          spd = max(spd, VMIN)
          g_taux_buf(k) = RHO_AIR_SFC * g_pool_ust(k)**2 * atm_public%u10(k) / spd
          g_tauy_buf(k) = RHO_AIR_SFC * g_pool_ust(k)**2 * atm_public%v10(k) / spd
        end do
      end if

    end block

    atm_state%running = .true.

    nullify(sfcInputPool, sst_field, ice_field, zorl_field, uocn_field, vocn_field)

  end subroutine mpas_atm_run

  ! ============================================================================
  !> @brief Finaliza o MONAN-A.
  !!
  !! Probe secao 5 / mpas_atm_core.F linha 1027:
  !!   function atm_core_finalize(domain) result(ierr)
  !! Probe mpas_framework.F linha 165:
  !!   subroutine mpas_framework_finalize(dminfo, domain, io_system)
  !!   io_system e OPCIONAL (mpas_subdriver linha 474 omite).
  !!
  !! Sequencia obrigatoria com SMIOL:
  !!   nullify(ponteiros zero-copy) -> core_finalize -> mpas_framework_finalize
  !!   -> deallocate(domain)
  ! ============================================================================
  subroutine mpas_atm_final(atm_public, atm_state, atm_bnd, rc)

    type(mpas_atm_public_type),    intent(inout) :: atm_public
    type(mpas_atm_state_type),     intent(inout) :: atm_state
    type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
    integer,                       intent(out)   :: rc

    integer :: ierr

    rc = 0

    if (.not. atm_state%initialized) then
      call mpas_log_write('mpas_atm_final: nada a finalizar')
      return
    end if

    if (associated(g_domain)) then

      ! IMPORTANTE: core_finalize e mpas_framework_finalize sao OMITIDAS.
      !
      ! core_finalize do MPAS-A (compilado com -DMPAS_EXTERNAL_ESMF_LIB) destroi
      ! internamente objetos ESMF_Time e ESMF_Calendar que o framework NUOPC
      ! ainda precisa para cleanup dos conectores (RouteHandles) apos ModelFinalize.
      ! Chamar core_finalize dentro de ESMF_GridCompFinalize -> SIGSEGV.
      !
      ! Os streams SMIOL ja foram fechados automaticamente no ultimo core_run
      ! (streams.atmosphere define alarm de output/restart). O restart final
      ! pode ser obtido configurando output_alarm no streams.atmosphere.
      !
      ! mpas_framework_finalize tambem omitida pelos mesmos motivos.
      ! A memoria e liberada pelo SO no termino do processo MPI.
      !
      ! Apenas nulifica ponteiros para evitar dangling references:
      nullify(atm_public%latCell,    atm_public%lonCell,  atm_public%areaCell)
      nullify(atm_public%t2m,        atm_public%u10,      atm_public%v10)
      nullify(atm_public%pslv)
      nullify(atm_public%lhflx,      atm_public%shflx)
      nullify(atm_public%swdn_sfc,   atm_public%lwdn_sfc, atm_public%prec_total)
      nullify(atm_public%taux_sfc,   atm_public%tauy_sfc)
      nullify(atm_public%q2m,        atm_public%prec_rain, atm_public%prec_snow)
      nullify(g_pool_acswdnb, g_pool_aclwdnb, g_pool_rainnc, g_pool_rainc)
      nullify(g_pool_snownc,  g_pool_q2,       g_pool_ust)
      g_domain => null()

      write(*,'(A)') 'mpas_atm_final: ponteiros nulificados (ESMF preservado)'
    end if

    ! 4. Desaloca apenas arrays de propriedade deste modulo
    if (allocated(g_lfrac_buf))            deallocate(g_lfrac_buf)
    if (allocated(atm_bnd%sst))           deallocate(atm_bnd%sst)
    if (allocated(atm_bnd%ice_fraction))  deallocate(atm_bnd%ice_fraction)
    if (allocated(atm_bnd%zorl))          deallocate(atm_bnd%zorl)
    if (allocated(atm_bnd%uocn))          deallocate(atm_bnd%uocn)
    if (allocated(atm_bnd%vocn))          deallocate(atm_bnd%vocn)
    if (allocated(atm_bnd%sss))           deallocate(atm_bnd%sss)
    ! Buffers de saida computados (propriedade deste modulo)
    if (allocated(g_prev_acswdnb)) deallocate(g_prev_acswdnb)
    if (allocated(g_prev_aclwdnb)) deallocate(g_prev_aclwdnb)
    if (allocated(g_prev_precip))  deallocate(g_prev_precip)
    if (allocated(g_prev_snow))    deallocate(g_prev_snow)
    if (allocated(g_swdn_inst))    deallocate(g_swdn_inst)
    if (allocated(g_lwdn_inst))    deallocate(g_lwdn_inst)
    if (allocated(g_prec_inst))    deallocate(g_prec_inst)
    if (allocated(g_taux_buf))     deallocate(g_taux_buf)
    if (allocated(g_tauy_buf))     deallocate(g_tauy_buf)
    if (allocated(g_q2m_buf))      deallocate(g_q2m_buf)
    if (allocated(g_prec_rain_buf))deallocate(g_prec_rain_buf)
    if (allocated(g_prec_snow_buf))deallocate(g_prec_snow_buf)

    atm_state%initialized = .false.
    atm_state%running     = .false.

  end subroutine mpas_atm_final

  ! ============================================================================

  ! ============================================================================
  subroutine mpas_atm_resize(atm_public, atm_state, atm_bnd, nCells_new)
    type(mpas_atm_public_type),    intent(inout) :: atm_public
    type(mpas_atm_state_type),     intent(inout) :: atm_state
    type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
    integer,                       intent(in)    :: nCells_new
    character(len=256) :: msg

    if (nCells_new == atm_public%nCells) then
      write(msg,'(A,I0,A)') 'mpas_atm_resize: nCells=', nCells_new, ' consistente'
      call mpas_log_write(trim(msg))
      return
    end if

    write(msg,'(A,I0,A,I0)') 'mpas_atm_resize: AVISO ESMF nCells=', nCells_new, &
         ' difere de MPAS nCells=', atm_public%nCells
    call mpas_log_write(trim(msg))
    write(*,'(A)') trim(msg)

    ! -- Atualiza contadores de celulas em atm_public e atm_state ---------
    atm_public%nCells = nCells_new
    atm_state%nCells  = nCells_new

    ! -- Redimensiona atm_bnd ----------------------------------------------
    if (allocated(atm_bnd%sst))          deallocate(atm_bnd%sst)
    if (allocated(atm_bnd%ice_fraction)) deallocate(atm_bnd%ice_fraction)
    if (allocated(atm_bnd%zorl))         deallocate(atm_bnd%zorl)
    if (allocated(atm_bnd%uocn))         deallocate(atm_bnd%uocn)
    if (allocated(atm_bnd%vocn))         deallocate(atm_bnd%vocn)
    if (allocated(atm_bnd%sss))          deallocate(atm_bnd%sss)

    allocate(atm_bnd%sst(nCells_new))
    allocate(atm_bnd%ice_fraction(nCells_new))
    allocate(atm_bnd%zorl(nCells_new))
    allocate(atm_bnd%uocn(nCells_new))
    allocate(atm_bnd%vocn(nCells_new))
    allocate(atm_bnd%sss(nCells_new))
    atm_bnd%sst          = real(cfg_sst_default,          MPAS_RKIND)
    atm_bnd%ice_fraction = real(cfg_ice_fraction_default, MPAS_RKIND)
    atm_bnd%zorl         = real(cfg_zorl_default,         MPAS_RKIND)
    atm_bnd%uocn         = real(cfg_uocn_default,         MPAS_RKIND)
    atm_bnd%vocn         = real(cfg_vocn_default,         MPAS_RKIND)
    atm_bnd%sss          = real(cfg_sss_default,          MPAS_RKIND)

    ! -- Redimensiona buffers de modulo (acumulados e stress) -------------
    ! Estes arrays sao alocados em mpas_atm_init com tamanho = MPAS nCells.
    ! Se o numero de celulas mudou (particionamento ESMF diferente do MPAS),
    ! os buffers devem ser realocados para evitar acesso fora dos limites em
    ! mpas_atm_run (loop 1..n onde n = atm_state%nCells).
    if (allocated(g_prev_acswdnb)) then
      deallocate(g_prev_acswdnb, g_prev_aclwdnb, g_prev_precip)
      deallocate(g_prev_snow)
      deallocate(g_swdn_inst, g_lwdn_inst, g_prec_inst)
      deallocate(g_taux_buf, g_tauy_buf)
      deallocate(g_q2m_buf, g_prec_rain_buf, g_prec_snow_buf)

      allocate(g_prev_acswdnb(nCells_new), g_prev_aclwdnb(nCells_new), &
               g_prev_precip(nCells_new))
      allocate(g_prev_snow(nCells_new))
      allocate(g_swdn_inst(nCells_new), g_lwdn_inst(nCells_new), &
               g_prec_inst(nCells_new))
      allocate(g_taux_buf(nCells_new), g_tauy_buf(nCells_new))
      allocate(g_q2m_buf(nCells_new), g_prec_rain_buf(nCells_new), &
               g_prec_snow_buf(nCells_new))

      ! Reinicializar com estado atual dos pools (se disponiveis)
      if (associated(g_pool_acswdnb)) then
        g_prev_acswdnb(1:nCells_new) = g_pool_acswdnb(1:nCells_new)
      else
        g_prev_acswdnb = 0.0_MPAS_RKIND
      end if
      if (associated(g_pool_aclwdnb)) then
        g_prev_aclwdnb(1:nCells_new) = g_pool_aclwdnb(1:nCells_new)
      else
        g_prev_aclwdnb = 0.0_MPAS_RKIND
      end if
      if (associated(g_pool_rainnc) .and. associated(g_pool_rainc)) then
        g_prev_precip(1:nCells_new) = g_pool_rainnc(1:nCells_new) &
                                     + g_pool_rainc(1:nCells_new)
      else if (associated(g_pool_rainnc)) then
        g_prev_precip(1:nCells_new) = g_pool_rainnc(1:nCells_new)
      else
        g_prev_precip = 0.0_MPAS_RKIND
      end if
      g_swdn_inst     = 0.0_MPAS_RKIND
      g_lwdn_inst     = 0.0_MPAS_RKIND
      g_prec_inst     = 0.0_MPAS_RKIND
      g_taux_buf      = 0.0_MPAS_RKIND
      g_tauy_buf      = 0.0_MPAS_RKIND
      g_q2m_buf       = 0.0_MPAS_RKIND
      g_prec_rain_buf = 0.0_MPAS_RKIND
      g_prec_snow_buf = 0.0_MPAS_RKIND

      ! Redirecionar ponteiros de atm_public para os novos buffers
      atm_public%swdn_sfc   => g_swdn_inst
      atm_public%lwdn_sfc   => g_lwdn_inst
      atm_public%prec_total => g_prec_inst
      atm_public%taux_sfc   => g_taux_buf
      atm_public%tauy_sfc   => g_tauy_buf
      atm_public%q2m        => g_q2m_buf
      atm_public%prec_rain  => g_prec_rain_buf
      atm_public%prec_snow  => g_prec_snow_buf
    end if

    write(msg,'(A,I0,A)') 'mpas_atm_resize: buffers realocados para ', &
         nCells_new, ' celulas'
    call mpas_log_write(trim(msg))

  end subroutine mpas_atm_resize

  ! --- auxiliares privados ---------------------------------------------------

  subroutine warn_if_null(ptr, name)
    real(MPAS_RKIND), pointer, intent(in) :: ptr(:)
    character(len=*),          intent(in) :: name
    character(len=256) :: msg
    if (.not. associated(ptr)) then
      write(msg,'(A,A,A)') 'AVISO: ponteiro nulo "', trim(name), &
           '" - verificar Registry.xml e namelist'
      call mpas_log_write(trim(msg))
      write(*,'(A)') trim(msg)
    end if
  end subroutine warn_if_null

end module mpas_atm_model_mod
