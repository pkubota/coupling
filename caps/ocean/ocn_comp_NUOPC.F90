!==============================================================================!
! ocn_comp_NUOPC.F90 - Modulo ponte para o cap NUOPC do oceano (MOM6+SIS2)    !
!                                                                              !
! FIX: O modulo original apenas fazia "use MOM_cap_mod" sem re-exportar       !
! nenhum simbolo, tornando-o inutil como ponto de entrada alternativo.        !
!                                                                              !
! Este modulo re-exporta explicitamente SetServices do MOM_cap_mod para que   !
! outros drivers possam usar:                                                  !
!   use ocn_comp_NUOPC, only: OCN_SetServices => SetServices                  !
! em lugar de usar MOM_cap_mod diretamente.                                   !
!                                                                              !
! Uso tipico em um driver alternativo:                                         !
!   use ocn_comp_NUOPC, only: OCN_SetServices => SetServices                  !
!   call NUOPC_DriverAddComp(driver, compLabel="OCN",                         !
!        compSetServicesRoutine=OCN_SetServices, ...)                          !
!==============================================================================!

module ocn_comp_NUOPC

  use MOM_cap_mod, only: SetServices

  implicit none
  private

  !> Ponto de entrada publico do componente oceanico MOM6+SIS2
  !! Registra todas as fases NUOPC (InitializeP0, AdvertiseFields,
  !! RealizeFields, DataInitialize, ModelAdvance, Finalize).
  public :: SetServices

end module ocn_comp_NUOPC
