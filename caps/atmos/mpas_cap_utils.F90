!> @file mpas_cap_utils.F90
!! @brief Utilitário de verificação de erros para o cap NUOPC do MPAS-A.

module mpas_cap_utils_mod

  use ESMF

  implicit none
  private

  public :: ChkErr

contains

  function ChkErr(rc, line, file) result(found_error)
    integer,          intent(in) :: rc
    integer,          intent(in) :: line
    character(len=*), intent(in) :: file
    logical :: found_error

    found_error = ESMF_LogFoundError(rcToCheck=rc,             &
                                     msg=ESMF_LOGERR_PASSTHRU, &
                                     line=line, file=file)
  end function ChkErr

end module mpas_cap_utils_mod
