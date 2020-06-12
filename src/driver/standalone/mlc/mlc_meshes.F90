module mlc_meshes

#include <petsc/finclude/petsc.h>

  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use mlc_global_vars
  use petscsys

  implicit none

  public :: setup_meshes

contains

  !------------------------------------------------------------------------
  subroutine setup_meshes(mlc_mpp)
    !
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp

    call add_meshes(mlc_mpp)
    call add_internal_connections(mlc_mpp)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine add_meshes(mlc_mpp)
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
    type(mpp_mlc_type) :: mlc_mpp
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
    PetscInt                   :: k, ncells, icair, icell, iconn
    PetscReal                  :: dz_cair

    call mlc_mpp%SetNumMeshes(2)

    ncells =(nz_cair+1)*ncair
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
       do k = 1, nz_cair
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

    vert_nconn = nz_cair*ncair
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

    !vert_conn_dist_up(1) = dz_cair/4.d0
    !vert_conn_dist_dn(1) = dz_cair/4.d0

    iconn = 0
    do icair = 1, ncair
       do k = 1, nz_cair
          iconn = iconn + 1
          vert_conn_id_up(iconn) = (nz_cair+1)*(icair-1) + k
          vert_conn_id_dn(iconn) = vert_conn_id_up(iconn) + 1
      end do
    enddo

    ! Canopy air space
    CAIR_MESH = 1

    allocate(mesh)
    call mesh%SetName('Canopy Air')
    call mesh%SetOrientation(MESH_AGAINST_GRAVITY)
    call mesh%SetDimensions((nz_cair+1)*ncair, 0, nz_cair+1)

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

    call mlc_mpp%AddMesh(CAIR_MESH, mesh)
    call mesh%Clean()
    deallocate(mesh)

    ! Canopy
    CLEF_MESH = 2

    allocate(mesh)
    call mesh%SetName('Canopy Air')
    call mesh%SetOrientation(MESH_AGAINST_GRAVITY)
    call mesh%SetDimensions((nz_cair+1)*ncair, 0, nz_cair+1)

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

    call mlc_mpp%AddMesh(CLEF_MESH, mesh)
    call mesh%Clean()
    deallocate(mesh)

  end subroutine add_meshes

  !------------------------------------------------------------------------
  subroutine add_internal_connections(mlc_mpp)
    !
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use MultiPhysicsProbConstants , only : CONN_SET_CONDITIONS
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    PetscInt            :: kk , nconn
    PetscInt  , pointer :: id_up(:)
    PetscInt  , pointer :: id_dn(:)
    PetscReal , pointer :: dist_up(:)
    PetscReal , pointer :: dist_dn(:)
    PetscReal , pointer :: area(:)
    PetscInt  , pointer :: itype(:)
    PetscReal , pointer :: unit_vec(:,:)

    nconn  = nz_cleaf

    allocate(id_up    (nconn   ))
    allocate(id_dn    (nconn   ))
    allocate(dist_up  (nconn   ))
    allocate(dist_dn  (nconn   ))
    allocate(area     (nconn   ))
    allocate(itype    (nconn   ))
    allocate(unit_vec (nconn,3 ))

    do kk = 1, nconn
       id_up(kk)      = 0
       id_dn(kk)      = kk
       dist_up(kk)    = z_cleaf/nz_cleaf/2.d0
       dist_dn(kk)    = z_cleaf/nz_cleaf/2.d0
       area(kk)       = 1.d0
       unit_vec(kk,1) = -1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_VERTICAL
    enddo

    CLEF_REGION_IN_CAIR_MESH = 1
    call mlc_mpp%CreateAndAddConnectionSet(CAIR_MESH, CONN_SET_CONDITIONS, &
         nconn,  id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    CAIR_REGION_IN_CLEF_MESH = 2
    unit_vec(:,1) = 1.d0
    call mlc_mpp%CreateAndAddConnectionSet(CLEF_MESH, CONN_SET_CONDITIONS, &
         nconn,  id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

  end subroutine add_internal_connections

end module mlc_meshes
