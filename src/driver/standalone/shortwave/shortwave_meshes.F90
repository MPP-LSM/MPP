module shortwave_meshes

#include <petsc/finclude/petsc.h>

  use MultiPhysicsProbShortwave , only : mpp_shortwave_type
  use shortwave_global_vars
  use petscsys

  implicit none

  public :: setup_meshes

contains

  !------------------------------------------------------------------------
  subroutine setup_meshes(shortwave_mpp)
    !
    implicit none
    !
    type(mpp_shortwave_type) :: shortwave_mpp

    call add_meshes(shortwave_mpp)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine add_meshes(shortwave_mpp)
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type) :: shortwave_mpp

    call shortwave_mpp%SetNumMeshes(1)
    call add_canopy_mesh(shortwave_mpp)

  end subroutine add_meshes

  !------------------------------------------------------------------------
  subroutine add_canopy_mesh(shortwave_mpp)
    !
    use MultiPhysicsProbConstants , only : CONN_IN_Z_DIR
    use MeshType                  , only : mesh_type, MeshCreate
    use MultiPhysicsProbConstants , only : MESH_AGAINST_GRAVITY
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type) :: shortwave_mpp
    !
    class(mesh_type) , pointer :: mesh
    PetscReal        , pointer :: dx(:), dy(:), dz(:)
    PetscReal        , pointer :: xc(:), yc(:), zc(:)
    PetscReal        , pointer :: area(:)
    PetscInt                   :: vert_nconn
    PetscInt         , pointer :: vert_conn_id_up(:)   !
    PetscInt         , pointer :: vert_conn_id_dn(:)   !
    PetscReal        , pointer :: vert_conn_dist_up(:) !
    PetscReal        , pointer :: vert_conn_dist_dn(:) !
    PetscReal        , pointer :: vert_conn_area(:)    !
    PetscInt         , pointer :: vert_conn_type(:)    !
    PetscInt                   :: k, ncells, icair, icell, iconn, itree, offset
    PetscReal                  :: dz_cair

    ncells =(nz_cair+1)*ncair*ntree
    allocate(xc(ncells))
    allocate(yc(ncells))
    allocate(zc(ncells))
    allocate(dx(ncells))
    allocate(dy(ncells))
    allocate(dz(ncells))
    allocate(area(ncells))

    dz_cair = z_cair/nz_cair;

    xc(:) = 0.d0
    yc(:) = 0.d0
    dx(:) = 1.d0
    dy(:) = 1.d0
    dz(:) = dz_cair
    area(:) = 1.d0

    icell = 0    
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, nz_cair+1
             icell = icell + 1
             if (k == 1) then
                zc(icell) = 0.d0
             else if (k == 2) then
                zc(icell) = dz_cair/2.d0
             else
                zc(icell) = zc(icell-1) + dz_cair
             endif
          end do
       end do
    end do

    vert_nconn = nz_cair*ncair*ntree
    allocate(vert_conn_id_up(vert_nconn))
    allocate(vert_conn_id_dn(vert_nconn))
    allocate(vert_conn_dist_up(vert_nconn))
    allocate(vert_conn_dist_dn(vert_nconn))
    allocate(vert_conn_area(vert_nconn))
    allocate(vert_conn_type(vert_nconn))

    vert_conn_area(:) = 1.d0
    vert_conn_type(:) = CONN_VERTICAL
    vert_conn_dist_up(:) = dz_cair/2.d0
    vert_conn_dist_dn(:) = dz_cair/2.d0

    iconn = 0
    offset = 0
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, nz_cair
             iconn = iconn + 1
             vert_conn_id_up(iconn) = offset + k
             vert_conn_id_dn(iconn) = vert_conn_id_up(iconn) + 1
          end do
          offset = offset + (nz_cair + 1)
       end do
    enddo

    ! Canopy
    shortwave_MESH = 1

    allocate(mesh)
    call mesh%SetName('Canopy Air')
    call mesh%SetOrientation(MESH_AGAINST_GRAVITY)
    call mesh%SetDimensions(ncells, 0, nz_cair+1)

    call mesh%SetGeometricAttributes(VAR_XC   , xc   )
    call mesh%SetGeometricAttributes(VAR_YC   , yc   )
    call mesh%SetGeometricAttributes(VAR_ZC   , zc   )
    call mesh%SetGeometricAttributes(VAR_DX   , dx   )
    call mesh%SetGeometricAttributes(VAR_DY   , dy   )
    call mesh%SetGeometricAttributes(VAR_DZ   , dz   )
    call mesh%SetGeometricAttributes(VAR_AREA , area )

    call mesh%ComputeVolume()

    call mesh%intrn_conn_set_list%Init()
    call mesh%lateral_conn_set_list%Init()
    call mesh%conditions_conn_set_list%Init()

    call mesh%CreateAndAddConnectionSet(CONN_SET_INTERNAL,      &
         vert_nconn,  vert_conn_id_up, vert_conn_id_dn,         &
         vert_conn_dist_up, vert_conn_dist_dn,  vert_conn_area, &
         vert_conn_type)

    call shortwave_mpp%AddMesh(shortwave_MESH, mesh)
    call mesh%Clean()
    deallocate(mesh)

    deallocate(xc)
    deallocate(yc)
    deallocate(zc)
    deallocate(dx)
    deallocate(dy)
    deallocate(dz)
    deallocate(area)
    deallocate(vert_conn_id_up)
    deallocate(vert_conn_id_dn)
    deallocate(vert_conn_dist_up)
    deallocate(vert_conn_dist_dn)
    deallocate(vert_conn_area)
    deallocate(vert_conn_type)

  end subroutine add_canopy_mesh

end module shortwave_meshes
