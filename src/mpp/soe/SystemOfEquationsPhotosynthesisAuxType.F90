module SystemOfEquationsPhotosynthesisAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use SystemOfEquationsLblAuxType, only : sysofeqns_lbl_auxvar_type
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  private

  type, public :: sysofeqns_photosynthesis_auxvar_type
     PetscReal :: ci
   contains
     procedure, public :: Init => LblSOEAuxVarInit
  end type sysofeqns_photosynthesis_auxvar_type

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
    class(sysofeqns_photosynthesis_auxvar_type) :: this

    call this%Init()
    
    this%ci  = 0.d0

  end subroutine LblSOEAuxVarInit

#endif

end module SystemOfEquationsPhotosynthesisAuxType
