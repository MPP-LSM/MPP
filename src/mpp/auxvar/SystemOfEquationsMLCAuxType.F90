
module SystemOfEquationsMlcAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                 , only : iulog
  use mpp_abortutils             , only : endrun
  use mpp_shr_log_mod            , only : errMsg => shr_log_errMsg
  use SystemOfEquationsBaseType  , only : sysofeqns_base_type
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscdm
  use petscdmda
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: sysofeqns_mlc_auxvar_type

     PetscReal :: tair
     PetscReal :: water_vapor
     PetscReal :: tleaf_sun
     PetscReal :: tleaf_shd

     PetscInt :: goveqn_id
     PetscInt :: condition_id

     PetscBool :: is_in                  ! [True/False] T = auxvar is for an internal control volume
     PetscBool :: is_bc                  ! [True/False] T = auxvar is for a boundary condition
     PetscBool :: is_ss                  ! [True/False] T = auxvar is for a source/sink condition
     PetscBool :: is_active              ! [True/False] T = auxvar is for a grid cell that is active

   contains
     procedure, public :: Init => MlcSOEAuxVarInit
  end type sysofeqns_mlc_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine MlcSOEAuxVarInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize an auxiliary variable
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_mlc_auxvar_type) :: this

    this%tair        = 0.d0
    this%water_vapor = 0.d0
    this%tleaf_sun   = 0.d0
    this%tleaf_shd   = 0.d0

    this%goveqn_id   = 0
    this%condition_id= 0

    this%is_in       = PETSC_FALSE
    this%is_bc       = PETSC_FALSE
    this%is_ss       = PETSC_FALSE
    this%is_active   = PETSC_TRUE

  end subroutine MlcSOEAuxVarInit

#endif

end module SystemOfEquationsMlcAuxType

