!> @file mpas_atm_types.F90
!! @brief Tipos publicos do cap MONAN-A 2.0 - sem dependencia direta do ESMF.
!!
!! Versao 6.0 - Campos ajustados ao mediador MED_cap (ESMF/NUOPC 8.9.1):
!!   mpas_atm_public_type agora inclui:
!!     q2m       : umidade especifica a 2 m [kg/kg]  -> Sa_shum_mpas
!!     prec_rain : precipitacao liquida  [kg/m^2/s]   -> Faxa_rain_mpas
!!     prec_snow : precipitacao solida   [kg/m^2/s]   -> Faxa_snow_mpas
!!
!!   O campo prec_total eh mantido para compatibilidade, mas NaO eh exportado
!!   ao mediador. O mediador espera rain e snow separados.
!!
!! Depende apenas de mpas_kind_types (sem ESMF).
!! Usado por: mpas_atm_model_mod, mpas_cap_methods_mod, mpas_cap_mod.

module mpas_atm_types_mod

  use mpas_kind_types, only : RKIND

  implicit none
  private

  ! -- Parametro de kind ------------------------------------------------------
  integer, parameter, public :: MPAS_RKIND = RKIND

  ! -- Campos diagnosticos exportados pelo MPAS-A para o mediador ------------
  !
  ! Mapeamento cap -> mediador (nomes NUOPC com sufixo _mpas):
  !   u10       -> Sa_u10m_mpas    vento zonal      10 m [m/s]
  !   v10       -> Sa_v10m_mpas    vento meridional 10 m [m/s]
  !   t2m       -> Sa_tbot_mpas    temperatura       2 m [K]
  !   q2m       -> Sa_shum_mpas    umidade especifica 2 m [kg/kg]
  !   pslv      -> Sa_pslv_mpas    pressao ao nivel do mar [Pa]
  !   swdn_sfc  -> Faxa_swdn_mpas  radiacao SW descendente [W/m^2]
  !   lwdn_sfc  -> Faxa_lwdn_mpas  radiacao LW descendente [W/m^2]
  !   prec_rain -> Faxa_rain_mpas  precipitacao liquida    [kg/m^2/s]
  !   prec_snow -> Faxa_snow_mpas  precipitacao solida     [kg/m^2/s]
  !
  type, public :: mpas_atm_public_type
    integer :: nCells      = 0  !< celulas locais incluindo halos (para zero-copy)
    integer :: nCellsSolve = 0  !< celulas proprias sem halos (B-32 - para NetCDF/export)
    integer :: nVertLevels = 0

    ! -- Geometria (ponteiros zero-copy -> pool 'mesh') ---------------------
    real(MPAS_RKIND), pointer :: latCell(:)    => null()  !< lat [rad]
    real(MPAS_RKIND), pointer :: lonCell(:)    => null()  !< lon [rad]
    real(MPAS_RKIND), pointer :: areaCell(:)   => null()  !< Area [m^2]

    ! -- Vento e temperatura em baixa atmosfera ----------------------------
    real(MPAS_RKIND), pointer :: t2m(:)        => null()  !< T a 2 m [K]
    real(MPAS_RKIND), pointer :: q2m(:)        => null()  !< Hum. especifica 2 m [kg/kg]
    real(MPAS_RKIND), pointer :: u10(:)        => null()  !< U a 10 m [m/s]
    real(MPAS_RKIND), pointer :: v10(:)        => null()  !< V a 10 m [m/s]

    ! -- Pressao -----------------------------------------------------------
    real(MPAS_RKIND), pointer :: pslv(:)       => null()  !< PSLV [Pa]

    ! -- Radiacao (medias do intervalo de acoplamento) ---------------------
    real(MPAS_RKIND), pointer :: swdn_sfc(:)   => null()  !< SWdn [W/m^2]
    real(MPAS_RKIND), pointer :: lwdn_sfc(:)   => null()  !< LWdn [W/m^2]

    ! -- Precipitacao (separada em liquida e solida) -----------------------
    real(MPAS_RKIND), pointer :: prec_rain(:)  => null()  !< Prec. liquida [kg/m^2/s]
    real(MPAS_RKIND), pointer :: prec_snow(:)  => null()  !< Prec. solida  [kg/m^2/s]
    !> Campo legado: prec_rain + prec_snow. Mantido para compatibilidade interna.
    real(MPAS_RKIND), pointer :: prec_total(:) => null()  !< Prec. total [kg/m^2/s]

    ! -- Fluxos turbulentos de superficie ---------------------------------
    real(MPAS_RKIND), pointer :: taux_sfc(:)   => null()  !< tau_x [N/m^2]
    real(MPAS_RKIND), pointer :: tauy_sfc(:)   => null()  !< tau_y [N/m^2]
    real(MPAS_RKIND), pointer :: lhflx(:)      => null()  !< LH [W/m^2]
    real(MPAS_RKIND), pointer :: shflx(:)      => null()  !< SH [W/m^2]

    ! Sa_lfrac_mpas: fracao de terra [0,1] para reconciliacao de mascara
    ! costeira no MED_cap. Derivado de xland (WRF: 1=terra,2=mar) ou
    ! landmask (1=terra,0=mar) do pool MPAS.
    real(MPAS_RKIND), pointer :: lfrac(:)      => null()  !< fracao terra [0-1] (Sa_lfrac_mpas)
  end type mpas_atm_public_type

  ! -- Estado interno do cap (sem tipos ESMF) ---------------------------------
  type, public :: mpas_atm_state_type
    logical            :: initialized    = .false.
    logical            :: running        = .false.
    character(len=256) :: config_dir     = './'
    character(len=64)  :: calendar_type  = 'gregorian'
    integer            :: dt_seconds     = 1800
    integer            :: nCells         = 0
    integer            :: nVertLevels    = 55
    integer            :: mpi_comm       = -1
  end type mpas_atm_state_type

  ! -- Condicoes de contorno vindas do oceano (via mediador) -----------------
  !
  ! Mapeamento mediador -> campo do cap:
  !   (conector OCN->MPAS)
  !   So_t  -> sst          SST [K]
  !   Si_ifrac -> ice_fraction  fracao de gelo [0-1]
  !   Sf_zorl  -> zorl       comprimento de rugosidade [m]
  !
  !   So_t     -> sst          SST [K]
  !   Si_ifrac -> ice_fraction fracao de gelo [0-1]
  !   Sf_zorl  -> zorl         comprimento de rugosidade [m]
  !   So_u     -> uocn         corrente zonal superficial [m/s]
  !   So_v     -> vocn         corrente meridional superficial [m/s]
  !   So_s     -> sss          salinidade superficial [psu]
  !
  ! Impacto das correntes na CLA do MPAS:
  !   uocn/vocn entram na velocidade relativa vento-corrente usada pelo
  !   esquema de fluxo de superficie (bl_mynn ou sfc_sfclay). Sem elas,
  !   o estresse tau e calculado com vento absoluto, subestimando a
  !   transferencia de momento em regioes de correntes intensas (ex: CN, AG).
  !   sss e reservado para uso futuro (parametrizacao de albedo de spray).

  type, public :: atm_ocean_boundary_type
    real(MPAS_RKIND), allocatable :: sst(:)          !< SST [K]
    real(MPAS_RKIND), allocatable :: ice_fraction(:) !< fracao de gelo [0-1]
    real(MPAS_RKIND), allocatable :: zorl(:)         !< rugosidade [m]
    real(MPAS_RKIND), allocatable :: uocn(:)         !< corrente zonal [m/s]   (So_u)
    real(MPAS_RKIND), allocatable :: vocn(:)         !< corrente meridional [m/s] (So_v)
    real(MPAS_RKIND), allocatable :: sss(:)          !< salinidade superficial [psu] (So_s)
  end type atm_ocean_boundary_type

end module mpas_atm_types_mod
