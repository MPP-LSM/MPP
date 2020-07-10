module photosynthesis_problem

  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbPhotosynthesis , only : mpp_photosynthesis_type
  use photosynthesis_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none


#include <petsc/finclude/petsc.h>

  type(mpp_photosynthesis_type):: phtsyn_mpp

  public :: run_photosynthesis_problem
  !public :: output_regression_photosynthesis_problem

contains

  !------------------------------------------------------------------------
  subroutine run_photosynthesis_problem(namelist_filename)
    !
    implicit none
    !
    character(len=256), optional :: namelist_filename
    !
    PetscReal          :: dt
    PetscBool          :: converged
    PetscInt           :: istep
    PetscInt           :: converged_reason
    PetscBool          :: flg
    PetscErrorCode     :: ierr

    ncair = 1;
    ntree = 1;

    call Init(phtsyn_mpp)

    dt = 0.0
    istep = 1
    call phtsyn_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine run_photosynthesis_problem

  !------------------------------------------------------------------------
  subroutine Init(phtsyn_mpp)
    !
    use photosynthesis_meshes     , only : setup_meshes
    use photosynthesis_parameters , only : set_parameters
    !
    implicit none
    !
    type(mpp_photosynthesis_type) :: phtsyn_mpp

    call initialize_mpp(phtsyn_mpp)

    call setup_meshes(phtsyn_mpp)

    call add_goveqn(phtsyn_mpp)

    call phtsyn_mpp%AllocateAuxVars()

    call phtsyn_mpp%SetupProblem()

    call set_parameters(phtsyn_mpp)

    call set_initial_condition(phtsyn_mpp)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(phtsyn_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_PHOTOSYNTHESIS_SNES
    !
    implicit none

    type(mpp_photosynthesis_type) :: phtsyn_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call phtsyn_mpp%Init       ()
    call phtsyn_mpp%SetName    ('Photosynthesis model')
    call phtsyn_mpp%SetID      (MPP_PHOTOSYNTHESIS_SNES)
    call phtsyn_mpp%SetMPIRank (iam)

  end subroutine Initialize_Mpp

  !------------------------------------------------------------------------
  subroutine add_goveqn(phtsyn_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_PHOTOSYNTHESIS
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_photosynthesis_type) :: phtsyn_mpp

    PHTSYN_GE = 1

    call phtsyn_mpp%AddGovEqnWithMeshRank(GE_PHOTOSYNTHESIS, 'Photosynthesis model', PHTSYN_MESH)
    
    call phtsyn_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqn

  !------------------------------------------------------------------------
  subroutine set_initial_condition(phtsyn_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_PHOTOSYNTHESIS_SNES
    !
    implicit none

    type(mpp_photosynthesis_type) :: phtsyn_mpp
    PetscScalar, pointer          :: ci_p(:)
    PetscErrorCode                :: ierr

    call VecGetArrayF90(phtsyn_mpp%soe%solver%soln, ci_p, ierr); CHKERRQ(ierr)
    ci_p(:) = 152.d0
    call VecRestoreArrayF90(phtsyn_mpp%soe%solver%soln, ci_p, ierr); CHKERRQ(ierr)

  end subroutine set_initial_condition

end module photosynthesis_problem
