module lbl

  use mpp_varctl               , only : iulog
  use mpp_abortutils           , only : endrun
  use mpp_shr_log_mod          , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbLBL , only : mpp_lbl_type
  use ml_model_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: init_lbl
  public :: solve_lbl
  public :: extract_data_from_lbl

contains

  !------------------------------------------------------------------------
  subroutine initialize(lbl_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_LBL_KSP
    !
    implicit none

    class(mpp_lbl_type) :: lbl_mpp

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

  end subroutine initialize

  !------------------------------------------------------------------------
  subroutine setup_meshes(lbl_mpp)
    !
    use MeshType             , only : mesh_type
    use ml_model_meshes      , only : create_canopy_mesh_for_leaf
    use ml_model_global_vars , only : LBL_MESH
    !
    implicit none

    class(mpp_lbl_type) :: lbl_mpp

    class(mesh_type), pointer :: mesh

    LBL_MESH = 1
    call create_canopy_mesh_for_leaf(mesh)
    call lbl_mpp%SetNumMeshes(1)
    call lbl_mpp%AddMesh(LBL_MESH, mesh)

    deallocate(mesh)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine add_goveqns(lbl_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_LEAF_BND_LAYER
    use ml_model_global_vars      , only : LBL_MESH, LBL_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_lbl_type) :: lbl_mpp

    LBL_GE = 1
    call lbl_mpp%AddGovEqnWithMeshRank(GE_LEAF_BND_LAYER, 'LBL equation', LBL_MESH)
    call lbl_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine set_parameters(lbl_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsLBLType  , only : sysofeqns_lbl_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnLeafBoundaryLayer   , only : goveqn_leaf_bnd_lyr_type
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use ml_model_global_vars      , only : LBL_MESH, LBL_GE
    use ml_model_meshes           , only : nleaf
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_lbl_type)                   :: lbl_mpp
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    class(connection_set_type) , pointer :: cur_conn_set
    class(condition_type)      , pointer :: cur_cond
    PetscInt                             :: k, icol, icell, ileaf, iconn, sum_conn, nz, ncol
    PetscReal                  , parameter :: dleaf = 0.04d0 ! [m]

    call lbl_mpp%soe%SetPointerToIthGovEqn(LBL_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_leaf_bnd_lyr_type)

       ncol = ncair * ntree
       nz   = (ntop-nbot+1)

       icell = 0
       do ileaf = 1, nleaf
          do icol = 1, ncol
             do k = 1, nz
                icell = icell + 1
                cur_goveq%aux_vars_in(icell)%dleaf = dleaf    ! [m]
             end do
          end do
       end do

    end select

  end subroutine set_parameters

  !------------------------------------------------------------------------
  subroutine lbl_set_boundary_conditions(lbl_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsLBLType  , only : sysofeqns_lbl_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnLeafBoundaryLayer   , only : goveqn_leaf_bnd_lyr_type
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use ml_model_global_vars      , only : Tleaf_sun, Tleaf_shd, Tair
    use ml_model_utils            , only : get_value_from_condition
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_lbl_type)                   :: lbl_mpp
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    class(connection_set_type) , pointer :: cur_conn_set
    class(condition_type)      , pointer :: cur_cond
    PetscInt                             :: k, icol, icell, ileaf, icair, idx_data, sum_conn, nz, ncol

    call lbl_mpp%soe%SetPointerToIthGovEqn(LBL_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_leaf_bnd_lyr_type)

       ncol = ntree
       nz   = (ntop-nbot+1)

       icell = 0
       ileaf = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1
             ileaf = ileaf + 1
             cur_goveq%aux_vars_in(icell           )%tleaf = get_value_from_condition(Tleaf_sun, ileaf) ! [K]
             cur_goveq%aux_vars_in(icell + ncol*nz )%tleaf = get_value_from_condition(Tleaf_shd, ileaf) ! [K]
          end do
       end do

       icell = 0
       idx_data = 0
       do icair = 1, ncair
          do k = 1, nz_cair
              idx_data = idx_data + 1
              if (k >= nbot .and. k<=ntop) then
                 icell = icell + 1

                 cur_goveq%aux_vars_in(icell           )%tair  = get_value_from_condition(Tair, idx_data)  ! [K]
                 cur_goveq%aux_vars_in(icell + ncol*nz )%tair  = get_value_from_condition(Tair, idx_data)  ! [K]

                 cur_goveq%aux_vars_in(icell           )%wind  = get_value_from_condition(Uref, icair)  ! [m/s]
                 cur_goveq%aux_vars_in(icell + ncol*nz )%wind  = get_value_from_condition(Uref, icair)  ! [m/s]

                 cur_goveq%aux_vars_in(icell           )%patm  = get_value_from_condition(Pref, icair)  ! [Pa]
                 cur_goveq%aux_vars_in(icell + ncol*nz )%patm  = get_value_from_condition(Pref, icair)  ! [Pa]
              endif
            enddo
       enddo

    end select

  end subroutine lbl_set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine init_lbl(lbl_mpp)
    !
    implicit none
    !
    class(mpp_lbl_type) :: lbl_mpp

    call initialize(lbl_mpp)

    call setup_meshes(lbl_mpp)

    call add_goveqns(lbl_mpp)

    call lbl_mpp%AllocateAuxVars()

    call lbl_mpp%SetupProblem()

    call set_parameters(lbl_mpp)

  end subroutine init_lbl

  !------------------------------------------------------------------------
  subroutine extract_data_from_lbl(lbl_mpp)
    !
    ! !DESCRIPTION:
    !   Extracts following variables from the LBL model:
    !     - gbh
    !     - gbv
    !     - gbc
    !
    ! !USES:
    use ml_model_global_vars            , only : nbot, ntop, ncair, ntree, nz_cair
    use ml_model_global_vars            , only : gbh, gbv, gbc
    use ml_model_meshes                 , only : nleaf
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnLeafBoundaryLayer         , only : goveqn_leaf_bnd_lyr_type
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_HEAT
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_H2O
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_CO2
    use ml_model_utils                  , only : set_value_in_condition
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_lbl_type)               :: lbl_mpp
    !
    PetscScalar             , pointer :: soln_p(:)
    PetscInt                          :: idx_leaf, idx_data, idx_soil, idx_air
    PetscInt                          :: ileaf, icair, itree, k, ieqn, icell
    PetscInt                          :: ncells
    PetscReal               , pointer :: gbh_data(:), gbv_data(:), gbc_data(:)
    class(goveqn_base_type) , pointer :: goveq
    PetscErrorCode                    :: ierr

    ncells = ncair * ntree * (ntop - nbot + 1) * nleaf

    allocate(gbh_data(ncells))
    allocate(gbv_data(ncells))
    allocate(gbc_data(ncells))

    call lbl_mpp%soe%SetPointerToIthGovEqn(1, goveq)

    select type(goveq)
    class is (goveqn_leaf_bnd_lyr_type)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_HEAT, ncells, gbh_data)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_H2O , ncells, gbv_data)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_CO2 , ncells, gbc_data)
    end select

    do icell = 1, ncells
       call set_value_in_condition(gbh, icell, gbh_data(icell))
       call set_value_in_condition(gbv, icell, gbv_data(icell))
       call set_value_in_condition(gbc, icell, gbc_data(icell))
    end do

    deallocate(gbh_data)
    deallocate(gbv_data)
    deallocate(gbc_data)

  end subroutine extract_data_from_lbl
 
 !------------------------------------------------------------------------
  subroutine solve_lbl(lbl_mpp, istep, dt)
    !
    implicit none
    !
    class(mpp_lbl_type) :: lbl_mpp
    PetscInt                 :: istep
    PetscReal                :: dt
    !
    PetscBool                :: converged
    PetscInt                 :: converged_reason
    PetscBool                :: flg
    PetscErrorCode           :: ierr

    call lbl_set_boundary_conditions(lbl_mpp)

    call lbl_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine solve_lbl

end module lbl
