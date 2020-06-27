module leafbndlyr_problem

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbLBL , only : mpp_lbl_type
  use lbl_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none


#include <petsc/finclude/petsc.h>

  type(mpp_lbl_type):: lbl_mpp

  public :: run_leafbndlyr_problem
  public :: output_regression_leafbndlyr_problem

contains

  !------------------------------------------------------------------------
  subroutine run_leafbndlyr_problem(namelist_filename)
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

    call Init(lbl_mpp)

    dt = 0.0
    istep = 1
    call lbl_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine run_leafbndlyr_problem

  !------------------------------------------------------------------------
  subroutine Init(lbl_mpp)
    !
    use lbl_meshes, only : setup_meshes
    use lbl_parameters, only : set_parameters
    !
    implicit none
    !
    type(mpp_lbl_type) :: lbl_mpp

    call initialize_mpp(lbl_mpp)

    call setup_meshes(lbl_mpp)

    call add_goveqn(lbl_mpp)

    call lbl_mpp%AllocateAuxVars()

    call lbl_mpp%SetupProblem()

    call set_parameters(lbl_mpp)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(lbl_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_LBL_KSP
    !
    implicit none

    type(mpp_lbl_type) :: lbl_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call lbl_mpp%Init       ()
    call lbl_mpp%SetName    ('Leaf boundary layer model')
    call lbl_mpp%SetID      (MPP_LBL_KSP)
    call lbl_mpp%SetMPIRank (iam)

  end subroutine Initialize_Mpp

  !------------------------------------------------------------------------
  subroutine add_goveqn(lbl_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_LEAF_BND_LAYER
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_lbl_type) :: lbl_mpp

    LBL_GE = 1

    call lbl_mpp%AddGovEqnWithMeshRank(GE_LEAF_BND_LAYER, 'leaf boundary layer', LBL_MESH)
    
    call lbl_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqn

  !------------------------------------------------------------------------
  subroutine output_regression_leafbndlyr_problem(filename_base, num_cells)
    !
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_HEAT
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_H2O
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_CO2
    use regression_mod                  , only : regression_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnLeafBoundaryLayer         , only : goveqn_leaf_bnd_lyr_type
    !
    implicit none
    !
    character(len=256)    :: filename_base
    PetscInt, intent(in)  :: num_cells
    !
    PetscInt              :: output, ieqn, ncells
    character(len=512)    :: filename
    character(len=64)     :: category
    character(len=64)     :: name
    PetscReal, pointer    :: data(:)
    type(regression_type) :: regression
    class(goveqn_base_type) , pointer :: goveq

    ncells = (nz_cair+1)*ncair

    allocate(data(ncells))

    call regression%Init(filename_base, num_cells)
    call regression%OpenOutput()
    
    do ieqn = 1,1

       call lbl_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

       select type(goveq)
       class is (goveqn_leaf_bnd_lyr_type)
          name = 'gbh'; category = 'general'
          call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_HEAT, ncells, data)
          call regression%WriteData(name, category, data)

          name = 'gbv'; category = 'general'
          call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_H2O, ncells, data)
          call regression%WriteData(name, category, data)

          name = 'gbc'; category = 'general'
          call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_CO2, ncells, data)

          call regression%WriteData(name, category, data)
       end select

    enddo

    call regression%CloseOutput()
    
    deallocate(data)

  end subroutine output_regression_leafbndlyr_problem

end module leafbndlyr_problem
