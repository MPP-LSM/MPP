module swv
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbShortwave , only : mpp_shortwave_type
  use ml_model_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: init_swv

contains

  !------------------------------------------------------------------------
  subroutine initialize(swv_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_shortwave_KSP
    !
    implicit none

    type(mpp_shortwave_type) :: swv_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call swv_mpp%Init       ()
    call swv_mpp%SetName    ('Shortgwave radiation model')
    call swv_mpp%SetID      (MPP_shortwave_KSP)
    call swv_mpp%SetMPIRank (iam)

  end subroutine initialize

  !------------------------------------------------------------------------
  subroutine init_swv(swv_mpp)
    !
    implicit none
    !
    type (mpp_shortwave_type) :: swv_mpp

    call initialize(swv_mpp)
    !call setup_meshes(swv_mpp)
    !call add_goveqn(swv_mpp)
    !call add_conditions_to_goveqns(swv_mpp)

    !call swv_mpp%AllocateAuxVars()
    !call swv_mpp%SetupProblem()
    !call set_parameters(swv_mpp)

  end subroutine init_swv

end module swv
