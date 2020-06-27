module SystemOfEquationsLblAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: sysofeqns_lbl_auxvar_type
     PetscReal :: gbh
     PetscReal :: gbv
     PetscReal :: gbc
   contains
     procedure, public :: Init => LblSOEAuxVarInit
  end type sysofeqns_lbl_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine LblSOEAuxVarInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize an auxiliary variable
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_lbl_auxvar_type) :: this

    this%gbh = 0.d0
    this%gbv = 0.d0
    this%gbc = 0.d0

  end subroutine LblSOEAuxVarInit

#endif

end module SystemOfEquationsLblAuxType
