module lwv

  use mpp_varctl               , only : iulog
  use mpp_abortutils           , only : endrun
  use mpp_shr_log_mod          , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbLongwave , only : mpp_longwave_type
  use ml_model_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: init_lwv
  public :: lwv_set_boundary_conditions
  public :: solve_lwv
  public :: extract_data_from_lwv

contains

  !------------------------------------------------------------------------
  subroutine initialize(lwv_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_LONGWAVE_KSP
    !
    implicit none

    class(mpp_longwave_type) :: lwv_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call lwv_mpp%Init       ()
    call lwv_mpp%SetName    ('Longwave radiation model')
    call lwv_mpp%SetID      (MPP_LONGWAVE_KSP)
    call lwv_mpp%SetMPIRank (iam)

  end subroutine initialize

  !------------------------------------------------------------------------
  subroutine setup_meshes(lwv_mpp)
    !
    use MeshType             , only : mesh_type
    use ml_model_meshes      , only : create_canopy_and_soil_mesh
    use ml_model_global_vars , only : LONGWAVE_MESH
    !
    implicit none

    class(mpp_longwave_type) :: lwv_mpp

    class(mesh_type), pointer :: mesh

    LONGWAVE_MESH = 1
    call create_canopy_and_soil_mesh(mesh)
    call lwv_mpp%SetNumMeshes(1)
    call lwv_mpp%AddMesh(LONGWAVE_MESH, mesh)

    deallocate(mesh)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine add_goveqns(lwv_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_LONGWAVE
    use ml_model_global_vars      , only : LONGWAVE_MESH, LONGWAVE_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_longwave_type) :: lwv_mpp

    LONGWAVE_GE = 1
    call lwv_mpp%AddGovEqnWithMeshRank(GE_LONGWAVE, 'Longwave radiation model', LONGWAVE_MESH)
    call lwv_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(lwv_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_LONGWAVE
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use ml_model_global_vars      , only : LONGWAVE_MESH, LONGWAVE_GE
    use ConnectionSetType         , only : connection_set_type
    use ml_model_meshes           , only : create_connection_set_to_canopy_and_soil_mesh
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_longwave_type)            :: lwv_mpp
    class(connection_set_type) , pointer :: conn_set    

    call create_connection_set_to_canopy_and_soil_mesh(lwv_mpp%meshes(LONGWAVE_MESH), conn_set)

    call lwv_mpp%soe%AddConditionInGovEqn( &
         LONGWAVE_GE                 , &
         ss_or_bc_type = COND_BC      , &
         name = 'Atmospheric forcing' , &
         unit = 'K'                   , &
         cond_type = COND_DIRICHLET   , &
         conn_set = conn_set)
    
  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine set_parameters(lwv_mpp, dpai)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsLongwaveType , only : sysofeqns_longwave_type
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnLongwaveType            , only : goveqn_longwave_type
    use ConditionType                  , only : condition_type
    use ConnectionSetType              , only : connection_set_type
    use ml_model_global_vars           , only : LONGWAVE_MESH, LONGWAVE_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_longwave_type)               :: lwv_mpp
    PetscReal, pointer, intent(in)         :: dpai(:)
    !
    class(goveqn_base_type)    , pointer   :: cur_goveq
    class(connection_set_type) , pointer   :: cur_conn_set
    class(condition_type)      , pointer   :: cur_cond
    PetscInt                               :: k, icol, icell, ileaf, iconn, sum_conn, nz, ncol

    PetscReal               , parameter :: emleaf = 0.98d0
    PetscReal               , parameter :: emgrnd = 0.96d0
    PetscReal               , parameter :: Irsky  = 400.d0
    PetscReal               , parameter :: td = 0.915d0

    call lwv_mpp%soe%SetPointerToIthGovEqn(LONGWAVE_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_longwave_type)

       ncol = ncair * ntree
       nz   = (ntop-nbot+1) + 1

       icell = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1

             cur_goveq%aux_vars_in(icell)%trans      = td
             cur_goveq%aux_vars_in(icell)%leaf_rho   = 1.d0 - emleaf
             cur_goveq%aux_vars_in(icell)%leaf_tau   = 0.d0
             cur_goveq%aux_vars_in(icell)%leaf_emiss = emleaf

             if (k == 1) then
                cur_goveq%aux_vars_in(icell)%is_soil     = PETSC_TRUE
                cur_goveq%aux_vars_in(icell)%soil_emiss = emgrnd
             else
                if (cur_goveq%aux_vars_in(icell)%nleaf /= 2) then
                   write(iulog,*)'Longwave model: number of leaves is not 2'
                   call exit(0)
                end if
                ileaf = 1; cur_goveq%aux_vars_in(icell)%leaf_fssh(ileaf) = fssh(nbot + k - 2)
                ileaf = 2; cur_goveq%aux_vars_in(icell)%leaf_fssh(ileaf) = 1.d0 -fssh(nbot + k - 2)
                ileaf = 1; cur_goveq%aux_vars_in(icell)%leaf_dlai(ileaf) = dpai(k + nbot - 2)
                ileaf = 2; cur_goveq%aux_vars_in(icell)%leaf_dlai(ileaf) = dpai(k + nbot - 2)
             end if
          end do
       end do

    end select

  end subroutine set_parameters

  !------------------------------------------------------------------------
  subroutine lwv_set_boundary_conditions(lwv_mpp, istep, isubstep)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType     , only : sysofeqns_base_type
    use SystemOfEquationsLongwaveType , only : sysofeqns_longwave_type
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnLongwaveType            , only : goveqn_longwave_type
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    use ml_model_utils                , only : get_value_from_condition
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_longwave_type)              :: lwv_mpp
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    class(connection_set_type) , pointer :: cur_conn_set
    class(condition_type)      , pointer :: cur_cond
    PetscInt                             :: k, icol, icell, iconn, sum_conn, nz, ncol, leaf_count, ileaf
    PetscInt :: istep, isubstep

    call lwv_mpp%soe%SetPointerToIthGovEqn(LONGWAVE_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_longwave_type)

       ncol = ncair * ntree
       nz   = (ntop-nbot+1) + 1

       icell = 0
       leaf_count = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1

             if (k == 1) then
                cur_goveq%aux_vars_in(icell)%soil_temperature = get_value_from_condition(tg, icol)
             else
                leaf_count = leaf_count + 1
                ileaf = 1; cur_goveq%aux_vars_in(icell)%leaf_temperature(ileaf) = get_value_from_condition(Tleaf_sun, leaf_count)
                ileaf = 2; cur_goveq%aux_vars_in(icell)%leaf_temperature(ileaf) = get_value_from_condition(Tleaf_shd, leaf_count)
                cur_goveq%aux_vars_in(icell)%trans = leaf_td(nbot + k - 2)
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

             cur_goveq%aux_vars_bc(sum_conn)%Idn = get_value_from_condition(Irsky, sum_conn)

          enddo
          cur_cond => cur_cond%next
       enddo

    end select

  end subroutine lwv_set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_number_of_leaves(lwv_mpp, nleaves)
    !
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnLongwaveType        , only : goveqn_longwave_type
    !
    implicit none
    !
    class(mpp_longwave_type)         :: lwv_mpp
    PetscInt                         :: nleaves
    !
    class(goveqn_base_type), pointer :: cur_goveq

    cur_goveq => lwv_mpp%soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_longwave_type)
          cur_goveq%nleaf = nleaves
       end select

       cur_goveq => cur_goveq%next
    end do

  end subroutine set_number_of_leaves

  !------------------------------------------------------------------------
  subroutine init_lwv(lwv_mpp)
    !
    use ml_model_global_vars, only : dpai
    !
    implicit none
    !
    class(mpp_longwave_type) :: lwv_mpp

    call initialize(lwv_mpp)

    call setup_meshes(lwv_mpp)

    call add_goveqns(lwv_mpp)

    call add_conditions_to_goveqns(lwv_mpp)

    call set_number_of_leaves(lwv_mpp, 2)

    call lwv_mpp%AllocateAuxVars()

    call lwv_mpp%SetupProblem()

    call set_parameters(lwv_mpp, dpai)

  end subroutine init_lwv

  !------------------------------------------------------------------------
  subroutine extract_data_from_lwv(lwv_mpp, istep, isubstep)
    !
    ! !DESCRIPTION:
    !   Extracts following variables from the longwave radiation model:
    !     - Iabs_leaf:
    !     - Iabs_soil
    !
    ! !USES:
    use ml_model_global_vars      , only : nbot, ntop, ncair, ntree, nz_cair
    use ml_model_global_vars      , only : Labs_leaf_sun, Labs_leaf_shd, Labs_soil
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnLongwaveType        , only : goveqn_longwave_type
    use MultiPhysicsProbLongwave  , only : mpp_longwave_type
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_LEAF_ABSORBED_LONGWAVE_RAD_PER_LAI
    use MultiPhysicsProbConstants , only : VAR_SOIL_ABSORBED_LONGWAVE_RAD_PER_GROUND
    use ml_model_utils            , only : set_value_in_condition
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_longwave_type)           :: lwv_mpp
    PetscInt :: istep, isubstep
    !
    PetscInt                            :: idx_leaf, idx_data, idx_soil, idx_air
    PetscInt                            :: icair, itree, k, leaf_icell, soil_icell
    PetscInt                            :: ncells, count
    PetscReal               , pointer   :: Labs_leaf_data(:), Labs_soil_data(:)
    class(goveqn_base_type) , pointer   :: goveq
    PetscErrorCode                      :: ierr

    ncells = ncair * ntree * (ntop - nbot + 1 + 1) ! num of canopy airspace x num. of tree x num. of levels
    allocate(Labs_leaf_data(ncells))
    allocate(Labs_soil_data(ncells))

    call lwv_mpp%soe%SetPointerToIthGovEqn(1, goveq)

    select type(goveq)
    class is (goveqn_longwave_type)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_ABSORBED_LONGWAVE_RAD_PER_LAI   , ncells, Labs_leaf_data)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_SOIL_ABSORBED_LONGWAVE_RAD_PER_GROUND, ncells, Labs_soil_data)
    end select

    count = 0
    leaf_icell = 0
    soil_icell = 0
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, ntop - nbot + 1 + 1
             count = count + 1;
             if (k == 1) then
                soil_icell = soil_icell + 1
                call set_value_in_condition(Labs_soil, soil_icell, Labs_soil_data(count))
             else
                leaf_icell = leaf_icell + 1
                call set_value_in_condition(Labs_leaf_sun, leaf_icell, Labs_leaf_data(count))
                call set_value_in_condition(Labs_leaf_shd, leaf_icell, Labs_leaf_data(count))
            end if
          end do
       end do
    end do

    deallocate(Labs_leaf_data)
    deallocate(Labs_soil_data)

  end subroutine extract_data_from_lwv

 !------------------------------------------------------------------------
  subroutine solve_lwv(lwv_mpp, istep, isubstep, dt)
    !
    implicit none
    !
    class(mpp_longwave_type) :: lwv_mpp
    PetscInt                 :: istep, isubstep
    PetscReal                :: dt
    !
    PetscBool                :: converged
    PetscInt                 :: converged_reason
    PetscBool                :: flg
    PetscErrorCode           :: ierr

    call lwv_set_boundary_conditions(lwv_mpp, istep, isubstep)

    call lwv_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine solve_lwv

end module lwv
