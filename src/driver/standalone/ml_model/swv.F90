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
  public :: swv_set_boundary_conditions
  public :: solve_swv

contains

  !------------------------------------------------------------------------
  subroutine initialize(swv_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_shortwave_KSP
    !
    implicit none

    class(mpp_shortwave_type) :: swv_mpp

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
  subroutine setup_meshes(swv_mpp)
    !
    use MeshType             , only : mesh_type
    use ml_model_meshes      , only : create_canopy_and_soil_mesh
    use ml_model_global_vars , only : SHORTWAVE_MESH
    !
    implicit none

    class(mpp_shortwave_type) :: swv_mpp

    class(mesh_type), pointer :: mesh

    SHORTWAVE_MESH = 1
    call create_canopy_and_soil_mesh(mesh)
    call swv_mpp%SetNumMeshes(1)
    call swv_mpp%AddMesh(SHORTWAVE_MESH, mesh)

    deallocate(mesh)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine add_goveqns(swv_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_SHORTWAVE
    use ml_model_global_vars      , only : SHORTWAVE_MESH, SHORTWAVE_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_shortwave_type) :: swv_mpp

    SHORTWAVE_GE = 1
    call swv_mpp%AddGovEqnWithMeshRank(GE_SHORTWAVE, 'Shortwave radiation model', SHORTWAVE_MESH)
    call swv_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(swv_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_SHORTWAVE
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use ml_model_global_vars      , only : SHORTWAVE_MESH, SHORTWAVE_GE
    use ConnectionSetType         , only : connection_set_type
    use ml_model_meshes           , only : create_connection_set_to_canopy_and_soil_mesh
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_shortwave_type)            :: swv_mpp
    class(connection_set_type) , pointer :: conn_set    

    call create_connection_set_to_canopy_and_soil_mesh(swv_mpp%meshes(SHORTWAVE_MESH), conn_set)

    call swv_mpp%soe%AddConditionInGovEqn( &
         SHORTWAVE_GE                 , &
         ss_or_bc_type = COND_BC      , &
         name = 'Atmospheric forcing' , &
         unit = 'K'                   , &
         cond_type = COND_DIRICHLET   , &
         conn_set = conn_set)
    
  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine set_parameters(swv_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsShortwaveType , only : sysofeqns_shortwave_type
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnShortwaveType            , only : goveqn_shortwave_type
    use ConditionType                  , only : condition_type
    use ConnectionSetType              , only : connection_set_type
    use ml_model_global_vars           , only : SHORTWAVE_MESH, SHORTWAVE_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type)               :: swv_mpp
    !
    class(goveqn_base_type)    , pointer   :: cur_goveq
    class(connection_set_type) , pointer   :: cur_conn_set
    class(condition_type)      , pointer   :: cur_cond
    PetscInt                               :: k, icol, icell, ileaf, iconn, sum_conn, nz, ncol

    PetscReal                  , parameter :: clumpfac  = 1.d0
    PetscReal                  , parameter :: Kb        = 0.577350269189626d0
    PetscReal                  , parameter :: td        = 0.913235689378651d0

    call swv_mpp%soe%SetPointerToIthGovEqn(SHORTWAVE_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_shortwave_type)

       ncol = ncair * ntree
       nz   = (ntop-nbot+1) + 1

       icell = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1

             if (k == 1) then
                cur_goveq%aux_vars_in(icell)%is_soil = PETSC_TRUE

                cur_goveq%aux_vars_in(icell)%soil_albedo_b(1) = 0.1d0 ! vis + direct
                cur_goveq%aux_vars_in(icell)%soil_albedo_b(2) = 0.2d0 ! nir + direct

                cur_goveq%aux_vars_in(icell)%soil_albedo_d(1) = 0.1d0 ! vis + diffuse
                cur_goveq%aux_vars_in(icell)%soil_albedo_d(2) = 0.2d0 ! nir + diffuse
             else
                cur_goveq%aux_vars_in(icell)%leaf_rho(1)   = 0.10d0
                cur_goveq%aux_vars_in(icell)%leaf_rho(2)   = 0.45d0

                cur_goveq%aux_vars_in(icell)%leaf_tau(1)   = 0.05d0
                cur_goveq%aux_vars_in(icell)%leaf_tau(2)   = 0.25d0

                cur_goveq%aux_vars_in(icell)%leaf_omega(1) = 0.55d0
                cur_goveq%aux_vars_in(icell)%leaf_omega(2) = 0.30d0

                cur_goveq%aux_vars_in(icell)%leaf_dlai        = dpai(k)
                cur_goveq%aux_vars_in(icell)%leaf_fraction(1) = fssh(k)
                cur_goveq%aux_vars_in(icell)%leaf_fraction(2) = 1.d0 - fssh(k)

                cur_goveq%aux_vars_in(icell)%leaf_tb    = exp(-Kb * dpai(k)   * clumpfac)
                cur_goveq%aux_vars_in(icell)%leaf_tbcum = exp(-Kb * cumlai(k) * clumpfac)
                cur_goveq%aux_vars_in(icell)%leaf_td    = td
             end if
          end do
       end do

    end select

  end subroutine set_parameters

  !------------------------------------------------------------------------
  subroutine swv_set_boundary_conditions(swv_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsShortwaveType , only : sysofeqns_shortwave_type
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnShortwaveType            , only : goveqn_shortwave_type
    use ConditionType                  , only : condition_type
    use ConnectionSetType              , only : connection_set_type
    use ml_model_utils                 , only : get_value_from_condition
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type)             :: swv_mpp
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    class(connection_set_type) , pointer :: cur_conn_set
    class(condition_type)      , pointer :: cur_cond
    PetscInt                             :: k, icol, icell, iconn, sum_conn, nz, ncol

    call swv_mpp%soe%SetPointerToIthGovEqn(SHORTWAVE_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_shortwave_type)

       ncol = ncair * ntree
       nz   = (ntop-nbot+1) + 1

       icell = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1
             if (k > 1) then
                cur_goveq%aux_vars_in(icell)%Iskyb(1) = get_value_from_condition(Iskyb_vis, icol)
                cur_goveq%aux_vars_in(icell)%Iskyb(2) = get_value_from_condition(Iskyb_nir, icol)
                cur_goveq%aux_vars_in(icell)%Iskyd(1) = get_value_from_condition(Iskyd_vis, icol)
                cur_goveq%aux_vars_in(icell)%Iskyd(2) = get_value_from_condition(Iskyd_nir, icol)

             end if
          end do
       end do

       sum_conn = 0
       cur_cond => cur_goveq%boundary_conditions%first
       do
          if (.not.associated(cur_cond)) exit

          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             icell = cur_conn_set%conn(iconn)%GetIDDn()

             cur_goveq%aux_vars_bc(sum_conn)%Iskyb(1) = get_value_from_condition(Iskyb_vis, sum_conn)
             cur_goveq%aux_vars_bc(sum_conn)%Iskyb(2) = get_value_from_condition(Iskyb_nir, sum_conn)
             cur_goveq%aux_vars_bc(sum_conn)%Iskyd(1) = get_value_from_condition(Iskyd_vis, sum_conn)
             cur_goveq%aux_vars_bc(sum_conn)%Iskyd(2) = get_value_from_condition(Iskyd_nir, sum_conn)

          enddo
          cur_cond => cur_cond%next
       enddo

    end select

  end subroutine swv_set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_number_of_leaves(swv_mpp, nleaves)
    !
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnShortwaveType       , only : goveqn_shortwave_type
    !
    implicit none
    !
    class(mpp_shortwave_type)        :: swv_mpp
    PetscInt                         :: nleaves
    !
    class(goveqn_base_type), pointer :: cur_goveq

    cur_goveq => swv_mpp%soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_shortwave_type)
          cur_goveq%nleaf = nleaves
       end select

       cur_goveq => cur_goveq%next
    end do

  end subroutine set_number_of_leaves

  !------------------------------------------------------------------------
  subroutine init_swv(swv_mpp)
    !
    implicit none
    !
    class(mpp_shortwave_type) :: swv_mpp

    call initialize(swv_mpp)

    call setup_meshes(swv_mpp)

    call add_goveqns(swv_mpp)

    call add_conditions_to_goveqns(swv_mpp)

    call set_number_of_leaves(swv_mpp, 2)

    call swv_mpp%AllocateAuxVars()

    call swv_mpp%SetupProblem()

    call set_parameters(swv_mpp)

  end subroutine init_swv

 !------------------------------------------------------------------------
  subroutine solve_swv(swv_mpp, istep, dt)
    !
    implicit none
    !
    class(mpp_shortwave_type) :: swv_mpp
    PetscInt                  :: istep
    PetscReal                 :: dt
    !
    PetscBool                 :: converged
    PetscInt                  :: converged_reason
    PetscBool                 :: flg
    PetscErrorCode            :: ierr

    call swv_set_boundary_conditions(swv_mpp)

    call swv_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine solve_swv

end module swv
