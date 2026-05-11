!> @file mpas_atm_wrappers.F90
!! @brief Wrappers não-módulo para as rotinas públicas de mpas_atm_model_mod.
!!
!! mpas_cap.F90 declara estas rotinas via bloco interface externo e as chama
!! diretamente por símbolo (sem use-association do módulo), evitando conflitos
!! de dependência entre os módulos do cap e do modelo.
!!
!! Cada wrapper apenas repassa os argumentos para a rotina de módulo correspondente.

subroutine mpas_atm_init(atm_public, atm_state, atm_bnd, &
                          dt_seconds, config_dir, mpi_comm, rc)
  use mpas_atm_types_mod, only : mpas_atm_public_type,    &
                                  mpas_atm_state_type,     &
                                  atm_ocean_boundary_type
  use mpas_atm_model_mod, only : impl => mpas_atm_init
  implicit none
  type(mpas_atm_public_type),    intent(inout) :: atm_public
  type(mpas_atm_state_type),     intent(inout) :: atm_state
  type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
  integer,          intent(in)  :: dt_seconds
  character(len=*), intent(in)  :: config_dir
  integer,          intent(in)  :: mpi_comm
  integer,          intent(out) :: rc
  call impl(atm_public, atm_state, atm_bnd, dt_seconds, config_dir, mpi_comm, rc)
end subroutine mpas_atm_init

subroutine mpas_atm_init_sfc(atm_public, atm_state, rc)
  use mpas_atm_types_mod, only : mpas_atm_public_type, mpas_atm_state_type
  use mpas_atm_model_mod, only : impl => mpas_atm_init_sfc
  implicit none
  type(mpas_atm_public_type), intent(inout) :: atm_public
  type(mpas_atm_state_type),  intent(inout) :: atm_state
  integer,                    intent(out)   :: rc
  call impl(atm_public, atm_state, rc)
end subroutine mpas_atm_init_sfc

subroutine mpas_atm_run(atm_public, atm_state, atm_bnd, dt_coupling, rc)
  use mpas_atm_types_mod, only : mpas_atm_public_type,    &
                                  mpas_atm_state_type,     &
                                  atm_ocean_boundary_type
  use mpas_atm_model_mod, only : impl => mpas_atm_run
  implicit none
  type(mpas_atm_public_type),    intent(inout) :: atm_public
  type(mpas_atm_state_type),     intent(inout) :: atm_state
  type(atm_ocean_boundary_type), intent(in)    :: atm_bnd
  integer,                       intent(in)    :: dt_coupling
  integer,                       intent(out)   :: rc
  call impl(atm_public, atm_state, atm_bnd, dt_coupling, rc)
end subroutine mpas_atm_run

subroutine mpas_atm_final(atm_public, atm_state, atm_bnd, rc)
  use mpas_atm_types_mod, only : mpas_atm_public_type,    &
                                  mpas_atm_state_type,     &
                                  atm_ocean_boundary_type
  use mpas_atm_model_mod, only : impl => mpas_atm_final
  implicit none
  type(mpas_atm_public_type),    intent(inout) :: atm_public
  type(mpas_atm_state_type),     intent(inout) :: atm_state
  type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
  integer,                       intent(out)   :: rc
  call impl(atm_public, atm_state, atm_bnd, rc)
end subroutine mpas_atm_final

subroutine mpas_atm_resize(atm_public, atm_state, atm_bnd, nCells_new)
  use mpas_atm_types_mod, only : mpas_atm_public_type, mpas_atm_state_type, &
                                  atm_ocean_boundary_type
  use mpas_atm_model_mod, only : impl => mpas_atm_resize
  implicit none
  type(mpas_atm_public_type),    intent(inout) :: atm_public
  type(mpas_atm_state_type),     intent(inout) :: atm_state
  type(atm_ocean_boundary_type), intent(inout) :: atm_bnd
  integer,                       intent(in)    :: nCells_new
  call impl(atm_public, atm_state, atm_bnd, nCells_new)
end subroutine mpas_atm_resize
