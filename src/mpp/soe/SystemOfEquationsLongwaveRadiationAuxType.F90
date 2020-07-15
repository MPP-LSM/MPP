module SystemOfEquationsLongwaveRadAuxType

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

  type, public :: sysofeqns_longwave_auxvar_type
     PetscReal :: Iup
     PetscReal :: Idn
     PetscReal :: Iabs
   contains
     procedure, public :: Init => LongwaveRadSOEAuxVarInit
  end type sysofeqns_longwave_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine LongwaveRadSOEAuxVarInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize an auxiliary variable
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_longwave_auxvar_type) :: this

    this%Iup  = 0.d0
    this%Idn  = 0.d0
    this%Iabs = 0.d0

  end subroutine LongwaveRadSOEAuxVarInit

#endif

end module SystemOfEquationsLongwaveRadAuxType
