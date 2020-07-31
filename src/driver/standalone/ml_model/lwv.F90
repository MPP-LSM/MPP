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
  subroutine set_parameters(lwv_mpp)

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
    !
    class(goveqn_base_type)    , pointer   :: cur_goveq
    class(connection_set_type) , pointer   :: cur_conn_set
    class(condition_type)      , pointer   :: cur_cond
    PetscInt                               :: k, icol, icell, ileaf, iconn, sum_conn, nz, ncol

    PetscReal               , parameter :: emleaf = 0.98d0
    PetscReal               , parameter :: emgrnd = 1.00d0
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
                ileaf = 1; cur_goveq%aux_vars_in(icell)%leaf_fraction(ileaf) = fssh(k)
                ileaf = 2; cur_goveq%aux_vars_in(icell)%leaf_fraction(ileaf) = 1.d0 -fssh(k)
             end if
          end do
       end do

    end select

  end subroutine set_parameters

  !------------------------------------------------------------------------
  subroutine lwv_set_boundary_conditions(lwv_mpp, Tsoil, Tleaf_sun, Tleaf_shd, Irsky)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType     , only : sysofeqns_base_type
    use SystemOfEquationsLongwaveType , only : sysofeqns_longwave_type
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnLongwaveType            , only : goveqn_longwave_type
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_longwave_type)              :: lwv_mpp
    PetscReal                  , pointer :: Tsoil(:)
    PetscReal                  , pointer :: Tleaf_sun(:)
    PetscReal                  , pointer :: Tleaf_shd(:)
    PetscReal                            :: Irsky
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    class(connection_set_type) , pointer :: cur_conn_set
    class(condition_type)      , pointer :: cur_cond
    PetscInt                             :: k, icol, icell, iconn, sum_conn, nz, ncol, leaf_count, ileaf

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
                cur_goveq%aux_vars_in(icell)%soil_temperature = Tsoil(icol) + 20.d0             
             else
                leaf_count = leaf_count + 1
                ileaf = 1; cur_goveq%aux_vars_in(icell)%leaf_temperature(ileaf) = Tleaf_sun(leaf_count) + 25.d0
                ileaf = 1; cur_goveq%aux_vars_in(icell)%leaf_temperature(ileaf) = Tleaf_shd(leaf_count) + 25.d0
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

             cur_goveq%aux_vars_bc(sum_conn)%Idn = Irsky

          enddo
          cur_cond => cur_cond%next
       enddo

    end select

  end subroutine lwv_set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine init_lwv(lwv_mpp)
    !
    implicit none
    !
    class(mpp_longwave_type) :: lwv_mpp

    call initialize(lwv_mpp)

    call setup_meshes(lwv_mpp)

    call add_goveqns(lwv_mpp)

    call add_conditions_to_goveqns(lwv_mpp)

    call lwv_mpp%AllocateAuxVars()

    call lwv_mpp%SetupProblem()

    call set_parameters(lwv_mpp)

  end subroutine init_lwv

end module lwv
