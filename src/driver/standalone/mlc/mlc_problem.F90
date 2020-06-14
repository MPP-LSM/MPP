module mlc_problem

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use mlc_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none

#include <petsc/finclude/petsc.h>

  public :: run_mlc_problem

contains

  !------------------------------------------------------------------------
  subroutine run_mlc_problem()

    use MultiPhysicsProbMLC, only : mpp_mlc_type

    type(mpp_mlc_type) :: mlc_mpp
    PetscReal          :: dt
    PetscBool          :: converged
    PetscInt           :: istep
    PetscInt           :: converged_reason
    PetscBool          :: flg
    PetscErrorCode     :: ierr

    ncair = 1;

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ncair',ncair,flg,ierr)

    call Init(mlc_mpp)
    !call mlc_mpp%soe%PrintInfo()

    dt = 5.0 * 60.d0 ! sec
    istep = 1
    call mlc_mpp%soe%SetDtime(dt)
    call mlc_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine run_mlc_problem

  !------------------------------------------------------------------------
  subroutine Init(mlc_mpp)

    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use mlc_conditions            , only : add_conditions_to_goveqns
    use mlc_meshes                , only : setup_meshes
    use mlc_parameters            , only : set_parameters

    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp(mlc_mpp)

    ! 2. Add all meshes needed for the MPP
    call setup_meshes(mlc_mpp)
    !call add_meshes(mlc_mpp)
    !call add_connection_sets_within_meshes(mlc_mpp)

    ! 3. Add govering equations
    call add_multiple_goveqns(mlc_mpp)

    ! 5. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns(mlc_mpp)
    !call add_conditions_to_goveqns(mlc_mpp)

    ! 6. Add internal coupling variables
    !call add_internal_coupling_vars(mlc_mpp)

    ! 7. Allocate auxvars
    call mlc_mpp%AllocateAuxVars(ncair)

    ! 8. Setup Problem
    call mlc_mpp%SetupProblem()

    ! 9.
    call set_parameters(mlc_mpp)

    ! 10. 
    call set_initial_conditions(mlc_mpp)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(mlc_mpp)
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_MLC_KSP
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call mlc_mpp%Init       ()
    call mlc_mpp%SetName    ('Multi Layer Canopy Model')
    call mlc_mpp%SetID      (MPP_MLC_KSP)
    call mlc_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp


  !------------------------------------------------------------------------
  subroutine add_multiple_goveqns(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_TEMP
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_VAPOR
    use MultiPhysicsProbConstants, only : GE_CANOPY_LEAF_TEMP
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp

    CAIR_TEMP_GE = 1
    CAIR_VAPR_GE = 2 
    CLEF_TEMP_SUN_GE = 3
    CLEF_TEMP_SHD_GE = 4

    call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_AIR_TEMP, 'Canopy air temperature', CAIR_MESH)
    
    call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_AIR_VAPOR, 'Canopy air vapor', CAIR_MESH)

    call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_LEAF_TEMP, 'Sunlit canopy', CLEF_MESH)

    call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_LEAF_TEMP, 'Shaded canopy', CLEF_MESH)

    call mlc_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_multiple_goveqns

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsMLCType  , only : sysofeqns_mlc_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscReal                            :: relhum, eref, esat, desatdt
    PetscInt                             :: p
    PetscInt                             :: nDM
    DM                         , pointer :: dms(:)
    Vec                        , pointer :: soln_subvecs(:)
    PetscReal                  , pointer :: v_p(:)
    PetscInt                             :: ii
    PetscViewer                          :: viewer
    PetscInt                             :: soe_auxvar_id
    PetscErrorCode                       :: ierr

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    p = 1
    soe%cturb%vcan(p) = soe%cturb%vref(p)
    soe%cturb%tcan(p) = soe%cturb%tref(p)

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(mlc_mpp%soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(mlc_mpp%soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(mlc_mpp%soe%solver%dm, &
         mlc_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)

       if (ii == CAIR_TEMP_GE .or. ii == CLEF_TEMP_SUN_GE .or. ii == CLEF_TEMP_SHD_GE) then
          v_p(:) = soe%cturb%tref(1)

       else if (ii == CAIR_VAPR_GE) then
          v_p(:) = soe%cturb%vref(1)          

       endif

       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(mlc_mpp%soe%solver%dm, &
         mlc_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)
    
    call VecCopy(mlc_mpp%soe%solver%soln, mlc_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(mlc_mpp%soe%solver%soln, mlc_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine set_initial_conditions

end module mlc_problem

