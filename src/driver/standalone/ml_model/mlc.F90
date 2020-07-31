module mlc

  use mpp_varctl               , only : iulog
  use mpp_abortutils           , only : endrun
  use mpp_shr_log_mod          , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use ml_model_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: init_mlc

contains

  !------------------------------------------------------------------------
  subroutine initialize(mlc_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_MLC_KSP
    !
    implicit none

    class(mpp_mlc_type) :: mlc_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call mlc_mpp%Init       ()
    call mlc_mpp%SetName    ('MLC radiation model')
    call mlc_mpp%SetID      (MPP_MLC_KSP)
    call mlc_mpp%SetMPIRank (iam)

  end subroutine initialize

  !------------------------------------------------------------------------
  subroutine setup_meshes(mlc_mpp)
    !
    use MeshType             , only : mesh_type
    use ml_model_meshes      , only : create_canopy_airspace_mesh
    use ml_model_meshes      , only : create_canopy_mesh
    use ml_model_global_vars , only : CAIR_MESH, CLEF_MESH
    !
    implicit none

    class(mpp_mlc_type) :: mlc_mpp

    class(mesh_type), pointer :: mesh

    call mlc_mpp%SetNumMeshes(2)

    CAIR_MESH = 1; call create_canopy_airspace_mesh(mesh); call mlc_mpp%AddMesh(CAIR_MESH, mesh)
    CLEF_MESH = 2; call create_canopy_mesh(mesh)         ; call mlc_mpp%AddMesh(CLEF_MESH, mesh)

    deallocate(mesh)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine setup_leaf2cair_map(mlc_mpp)
    !
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType , only : goveqn_cair_temp_type
    use GoveqnCanopyAirVaporType       , only : goveqn_cair_vapor_type
    use GoveqnCanopyLeafTemperatureType, only : goveqn_cleaf_temp_type
    !
    implicit none
    !
    class(mpp_mlc_type)          :: mlc_mpp
    !
    class(goveqn_base_type) , pointer :: goveq
    PetscInt                , pointer :: leaf2cair(:)
    PetscInt                          :: nleaf2cair, icair, itree, k, count, ieqn

    if (ntree == 1) then

       do ieqn = 1, 4
          call mlc_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

          select type(goveq)
          class is (goveqn_cair_temp_type)
             call goveq%SetDefaultLeaf2CAirMap()
          class is (goveqn_cair_vapor_type)
             call goveq%SetDefaultLeaf2CAirMap()
          class is (goveqn_cleaf_temp_type)
             call goveq%SetDefaultLeaf2CAirMap()
          end select
       end do
    else

       nleaf2cair = (nz_cair+1)*ncair*ntree

       allocate(leaf2cair(nleaf2cair))

       count = 0
       do icair = 1, ncair
          do itree = 1, ntree
             do k = 1, nz_cair + 1
                count = count + 1
                leaf2cair(count) = (nz_cair + 1)*(icair-1) + k
             end do
          end do
       end do

       do ieqn = 1, 4
          call mlc_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

          select type(goveq)
          class is (goveqn_cair_temp_type)
             call goveq%SetLeaf2CAirMap(leaf2cair, nleaf2cair)
          class is (goveqn_cair_vapor_type)
             call goveq%SetLeaf2CAirMap(leaf2cair, nleaf2cair)
          class is (goveqn_cleaf_temp_type)
             call goveq%SetLeaf2CAirMap(leaf2cair, nleaf2cair)
          end select
       end do

    end if

  end subroutine setup_leaf2cair_map

  !------------------------------------------------------------------------
  subroutine add_goveqns(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_CANOPY_AIR_TEMP
    use MultiPhysicsProbConstants , only : GE_CANOPY_AIR_VAPOR
    use MultiPhysicsProbConstants , only : GE_CANOPY_LEAF_TEMP
    use ml_model_global_vars      , only : CAIR_TEMP_GE, CAIR_VAPR_GE
    use ml_model_global_vars      , only : CLEF_TEMP_SUN_GE, CLEF_TEMP_SUN_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp

    CAIR_TEMP_GE     = 1; call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_AIR_TEMP , 'Canopy air temperature' , CAIR_MESH)
    CAIR_VAPR_GE     = 2; call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_AIR_VAPOR, 'Canopy air vapor'       , CAIR_MESH)
    CLEF_TEMP_SUN_GE = 3; call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_LEAF_TEMP, 'Sunlit canopy'          , CLEF_MESH)
    CLEF_TEMP_SHD_GE = 4; call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_LEAF_TEMP, 'Shaded canopy'          , CLEF_MESH)

    call mlc_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine init_mlc(mlc_mpp)
    !
    use mlc_conditions, only : mlc_add_conditions_to_goveqns
    use mlc_parameters, only : mlc_set_parameters
    !
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp

    call initialize(mlc_mpp)

    call setup_meshes(mlc_mpp)

    call add_goveqns(mlc_mpp)

    call setup_leaf2cair_map(mlc_mpp)

    call mlc_add_conditions_to_goveqns(mlc_mpp)

    call mlc_mpp%AllocateAuxVars(ncair)

    call mlc_mpp%SetupProblem()

    call mlc_set_parameters(mlc_mpp)

  end subroutine init_mlc

end module mlc
