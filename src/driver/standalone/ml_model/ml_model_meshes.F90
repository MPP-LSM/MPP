module ml_model_meshes
  !
  use mpp_varctl           , only : iulog
  use mpp_abortutils       , only : endrun
  use mpp_shr_log_mod      , only : errMsg => shr_log_errMsg
  !
  implicit none
  !
#include <petsc/finclude/petsc.h>
  PetscInt, parameter :: CANOPY_AIR_MESH      = 1
  PetscInt, parameter :: CANOPY_MESH          = 2
  PetscInt, parameter :: CANOPY_AND_SOIL_MESH = 3

  public :: create_canopy_airspace_mesh
  public :: create_canopy_mesh
  public :: create_canopy_and_soil_mesh
  public :: create_connection_set_to_canopy_airspace_mesh
  public :: create_connection_set_to_canopy_and_soil_mesh

contains

  !------------------------------------------------------------------------
  subroutine create_canopy_airspace_mesh(mesh)
    !
    use MeshType, only : mesh_type
    !
    implicit none
    !
    class(mesh_type) , pointer :: mesh

    call create_mesh(mesh, CANOPY_AIR_MESH)

    call mesh%SetName('Canopy Air')

  end subroutine create_canopy_airspace_mesh

  !------------------------------------------------------------------------
  subroutine create_canopy_mesh(mesh)
    !
    use MeshType , only : mesh_type
    !
    implicit none
    !
    class(mesh_type) , pointer :: mesh

    call create_mesh(mesh, CANOPY_MESH)

    call mesh%SetName('Canopy')

  end subroutine create_canopy_mesh

  !------------------------------------------------------------------------
  subroutine create_canopy_and_soil_mesh(mesh)
    !
    use MeshType, only : mesh_type
    !
    implicit none
    !
    class(mesh_type) , pointer :: mesh

    call create_mesh(mesh, CANOPY_AND_SOIL_MESH)

    call mesh%SetName('Canopy with soil layer')

  end subroutine create_canopy_and_soil_mesh

  !------------------------------------------------------------------------
  subroutine create_mesh(mesh, mesh_id)
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
    use ml_model_global_vars      , only : ncair, nz_cair, dz_cair, ntop, nbot, ntree
    ! !ARGUMENTS
    implicit none
    !
    class(mesh_type) , pointer :: mesh
    PetscInt                   :: mesh_id
    !
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
    PetscInt                   :: k, ncells
    PetscInt                   :: icol, itree, icell
    PetscInt                   :: nz, ncol

    select case(mesh_id)
    case (CANOPY_AIR_MESH)
       nz   = nz_cair + 1
       ncol = ncair

    case (CANOPY_MESH)
       nz   = (ntop-nbot+1)
       ncol = ncair * ntree * 2

    case (CANOPY_AND_SOIL_MESH)
       nz   = (ntop-nbot+1) + 1
       ncol = ncair * ntree

    end select

    ncells = ncol * nz

    allocate(xc(ncells))
    allocate(yc(ncells))
    allocate(zc(ncells))
    allocate(dx(ncells))
    allocate(dy(ncells))
    allocate(dz(ncells))
    allocate(area(ncells))

    xc(:)   = 0.d0
    yc(:)   = 0.d0
    dx(:)   = 1.d0
    dy(:)   = 1.d0
    dz(:)   = dz_cair
    area(:) = 1.d0

    select case(mesh_id)
    case (CANOPY_AIR_MESH)
       icell = 0
       do icol = 1, ncol
          do k = 1, nz
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

    case (CANOPY_MESH)
       nz = (ntop-nbot+1)

       !write(*,*)'CANOPY_MESH:'
       icell = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1
             zc(icell) = (nbot + k - 1)* dz_cair + dz_cair/2.d0
             !write(*,*) icell,zc(icell)
          end do
       end do

    case (CANOPY_AND_SOIL_MESH)
       nz = (ntop-nbot+1) + 1

       !write(*,*)'CANOPY_AND_SOIL_MESH: '
       icell = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1
             if (k == 1) then
                zc(icell) = 0.d0
             else
                zc(icell) = (nbot + k - 2)* dz_cair + dz_cair/2.d0
             end if
             !write(*,*) icell,zc(icell)
                
          end do
       end do

    end select

    call SetupVerticalConnection( &
         ncol              , &
         nz-1              , &
         dz_cair           , &
         vert_nconn        , &
         vert_conn_id_up   , &
         vert_conn_id_dn   , &
         vert_conn_dist_up , &
         vert_conn_dist_dn , &
         vert_conn_area    , &
         vert_conn_type      &
         )

    allocate(mesh)
    call mesh%SetOrientation(MESH_AGAINST_GRAVITY)
    call mesh%SetDimensions(ncells, 0, nz)

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

  end subroutine create_mesh

  !------------------------------------------------------------------------
  subroutine SetupVerticalConnection( &
       ncol              , &
       nz                , &
       dz                , &
       vert_nconn        , &
       vert_conn_id_up   , &
       vert_conn_id_dn   , &
       vert_conn_dist_up , &
       vert_conn_dist_dn , &
       vert_conn_area    , &
       vert_conn_type)
    !
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    !
    implicit none
    !
    PetscInt                          :: ncol, nz
    PetscReal                         :: dz
    PetscInt            , intent(inout) :: vert_nconn
    PetscInt  , pointer , intent(inout) :: vert_conn_id_up(:)
    PetscInt  , pointer , intent(inout) :: vert_conn_id_dn(:)
    PetscReal , pointer , intent(inout) :: vert_conn_dist_up(:)
    PetscReal , pointer , intent(inout) :: vert_conn_dist_dn(:)
    PetscReal , pointer , intent(inout) :: vert_conn_area(:)
    PetscInt  , pointer , intent(inout) :: vert_conn_type(:)
    !
    PetscInt            :: k, ncells, icol, iconn

    vert_nconn = ncol*nz

    allocate(vert_conn_id_up(vert_nconn))
    allocate(vert_conn_id_dn(vert_nconn))
    allocate(vert_conn_dist_up(vert_nconn))
    allocate(vert_conn_dist_dn(vert_nconn))
    allocate(vert_conn_area(vert_nconn))
    allocate(vert_conn_type(vert_nconn))

    vert_conn_area(:)    = 1.d0
    vert_conn_type(:)    = CONN_VERTICAL
    vert_conn_dist_up(:) = dz/2.d0
    vert_conn_dist_dn(:) = dz/2.d0

    iconn = 0
    do icol = 1, ncol
       do k = 1, nz
          iconn = iconn + 1
          vert_conn_id_up(iconn) = (nz+1)*(icol-1) + k
          vert_conn_id_dn(iconn) = vert_conn_id_up(iconn) + 1
       end do
    enddo
  end subroutine SetupVerticalConnection

  !------------------------------------------------------------------------
  subroutine create_connection_set_to_canopy_airspace_mesh(mesh, conn_set)
    !
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : mesh_type
    use ml_model_global_vars      , only : ncair, nz_cair
    !
    implicit none
    !
    class(mesh_type)                    :: mesh
    class(connection_set_type), pointer :: conn_set

    call create_connection_set(mesh, conn_set, ncair, nz_cair + 1)

  end subroutine create_connection_set_to_canopy_airspace_mesh

  !------------------------------------------------------------------------
  subroutine create_connection_set_to_canopy_and_soil_mesh(mesh, conn_set)
    !
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : mesh_type
    use ml_model_global_vars      , only : ncair, ntree, ntop, nbot
    !
    implicit none
    !
    class(mesh_type)                    :: mesh
    class(connection_set_type), pointer :: conn_set

    call create_connection_set(mesh, conn_set, ncair*ntree, (ntop-nbot+1) + 1)

  end subroutine create_connection_set_to_canopy_and_soil_mesh

  !------------------------------------------------------------------------
  subroutine create_connection_set(mesh, conn_set, ncol, nz)
    !
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet, mesh_type
    use ml_model_global_vars      , only : dz_cair
    !
    implicit none
    !
    class(mesh_type)                     :: mesh
    class(connection_set_type) , pointer :: conn_set
    PetscInt                             :: ncol, nz
    !
    PetscInt                             :: iconn, nconn
    PetscReal                            :: dz
    PetscInt                   , pointer :: id_up(:)
    PetscInt                   , pointer :: id_dn(:)
    PetscReal                  , pointer :: dist_up(:)
    PetscReal                  , pointer :: dist_dn(:)
    PetscReal                  , pointer :: area(:)
    PetscInt                   , pointer :: itype(:)
    PetscReal                  , pointer :: unit_vec(:,:)

    nconn  = ncol

    allocate(id_up    (nconn   ))
    allocate(id_dn    (nconn   ))
    allocate(dist_up  (nconn   ))
    allocate(dist_dn  (nconn   ))
    allocate(area     (nconn   ))
    allocate(itype    (nconn   ))
    allocate(unit_vec (nconn,3 ))

    dz = dz_cair

    do iconn = 1, nconn
       id_up(iconn)      = 0
       id_dn(iconn)      = (nz)*iconn
       dist_up(iconn)    =  0.d0
       dist_dn(iconn)    = dz
       area(iconn)       = 1.d0
       unit_vec(iconn,1) = -1.d0
       unit_vec(iconn,2) = 0.d0
       unit_vec(iconn,3) = 0.d0
       itype(iconn)      = CONN_VERTICAL
    end do

    allocate(conn_set)
    call MeshCreateConnectionSet(mesh, &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, conn_set)

  end subroutine create_connection_set

  end module ml_model_meshes
