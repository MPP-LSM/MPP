module mpp_mesh_utils

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  use petscsys
  use mpp_varctl              , only : iulog
  use mpp_abortutils          , only : endrun
  use mpp_shr_log_mod         , only : errMsg => shr_log_errMsg

  implicit none

  public :: ComputeXYZCentroids
  public :: ComputeCentroids
  public :: ComputeInternalConnections
  public :: ComputeBoundaryDomainConnection
  public :: ComputeLeftBoundaryDomainConnection
  public :: ComputeRightBoundaryDomainConnection
  public :: ComputeTopBoundaryDomainConnection
  public :: ComputeBottomBoundaryDomainConnection
  public :: ComputeSouthBoundaryDomainConnection
  public :: ComputeNorthBoundaryDomainConnection

  interface ComputeXYZCentroids
     procedure ComputeXYZCentroids1
     procedure ComputeXYZCentroids2
     procedure ComputeXYZCentroids3
     procedure ComputeXYZCentroids4
  end interface ComputeXYZCentroids

  interface ComputeCentroids
     procedure ComputeCentroids1
     procedure ComputeCentroids2
  end interface ComputeCentroids

contains

  !------------------------------------------------------------------------
  subroutine ComputeXYZCentroids1( &
       nx, ny, nz,                  &
       dx, dy, dz,                  &
       x_min, y_min, z_min,         &
       xc, yc, zc)
    !
    ! !DESCRIPTION:
    ! Computes centroids along X, Y, and Z direction for a structured
    ! grid that has spatially homogeneous dx, dy, and dz. Returns values
    ! in a 1D array.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_XC, VAR_YC, VAR_ZC
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(in)           :: dx, dy, dz
    PetscReal, intent(in)           :: x_min, y_min, z_min
    PetscReal, intent(out), pointer :: xc(:), yc(:), zc(:)
    !
    PetscInt :: ii, jj, kk
    PetscInt :: count

    allocate(xc(nx*ny*nz))
    allocate(yc(nx*ny*nz))
    allocate(zc(nx*ny*nz))

    call ComputeCentroids(nx, ny, nz, dx, x_min, VAR_XC, xc)
    call ComputeCentroids(nx, ny, nz, dy, y_min, VAR_YC, yc)
    call ComputeCentroids(nx, ny, nz, dz, z_min, VAR_ZC, zc)

  end subroutine ComputeXYZCentroids1

  !------------------------------------------------------------------------
  subroutine ComputeXYZCentroids2( &
       nx, ny, nz,                  &
       dx, dy, dz,                  &
       xc, yc, zc)
    !
    ! !DESCRIPTION:
    ! Computes centroids along X, Y, and Z direction for a structured
    ! grid that has spatially homogeneous dx, dy, and dz. Returns values
    ! in a 1D array.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_XC, VAR_YC, VAR_ZC
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(in)           :: dx, dy, dz
    PetscReal, intent(out), pointer :: xc(:), yc(:), zc(:)

    call ComputeXYZCentroids( &
         nx = nx, ny = ny, nz = nz, &
         dx = dx, dy = dy, dz = dz, &
         x_min = 0.d0, y_min = 0.d0, z_min = 0.d0, &
         xc = xc, yc = yc, zc = zc)
    
  end subroutine ComputeXYZCentroids2

  !------------------------------------------------------------------------
  subroutine ComputeXYZCentroids3( &
       nx, ny, nz,                 &
       dx, dy, dz,                 &
       x_min, y_min, z_min,        &
       xc, yc, zc)
    !
    ! !DESCRIPTION:
    ! Same as ComputeXYZCentroids1D except returns values in a 1D array.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_XC, VAR_YC, VAR_ZC
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(in)           :: dx, dy, dz
    PetscReal, intent(in)           :: x_min, y_min, z_min
    PetscReal, intent(out), pointer :: xc(:,:,:), yc(:,:,:), zc(:,:,:)
    !
    PetscInt                        :: ii, jj, kk
    PetscReal, pointer              :: xc_1D(:), yc_1D(:), zc_1D(:)

    allocate(xc(nx,ny,nz))
    allocate(yc(nx,ny,nz))
    allocate(zc(nx,ny,nz))
    
    call ComputeXYZCentroids( &
         nx = nx, ny = ny, nz = nz, &
         dx = dx, dy = dy, dz = dz, &
         x_min = x_min, y_min = y_min, z_min = z_min, &
         xc = xc_1D, yc = yc_1D, zc = zc_1D)
    
    call Remap1Dto3D(nx, ny, nz, xc_1D, xc)
    call Remap1Dto3D(nx, ny, nz, yc_1D, yc)
    call Remap1Dto3D(nx, ny, nz, zc_1D, zc)

    deallocate(xc_1D)
    deallocate(yc_1D)
    deallocate(zc_1D)

  end subroutine ComputeXYZCentroids3

  !------------------------------------------------------------------------
  subroutine ComputeXYZCentroids4( &
       nx, ny, nz,                 &
       dx, dy, dz,                 &
       xc, yc, zc)
    !
    ! !DESCRIPTION:
    ! Same as ComputeXYZCentroids1D except returns values in a 1D array.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_XC, VAR_YC, VAR_ZC
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(in)           :: dx, dy, dz
    PetscReal, intent(out), pointer :: xc(:,:,:), yc(:,:,:), zc(:,:,:)
    !
    PetscInt                        :: ii, jj, kk

    call ComputeXYZCentroids( &
         nx = nx, ny = ny, nz = nz, &
         dx = dx, dy = dy, dz = dz, &
         x_min = 0.d0, y_min = 0.d0, z_min = 0.d0, &
         xc = xc, yc = yc, zc = zc)
    
  end subroutine ComputeXYZCentroids4

  !------------------------------------------------------------------------
  subroutine ComputeCentroids1(nx, ny, nz, &
       dd, dmin, idir, &
       dc)
    !
    ! !DESCRIPTION:
    ! Computes centroids along 'idir' direction for a structured
    ! grid that has spatially homogeneous dx, dy, and dz. Returns
    ! values in a 1D array.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_XC, VAR_YC, VAR_ZC
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(in)           :: dd, dmin
    PetscInt , intent(in)           :: idir
    PetscReal, intent(out), pointer :: dc(:)
    !
    PetscInt :: ii, jj, kk
    PetscInt :: index
    PetscInt :: count

    count = 0
    do kk = 1,nz
       do jj = 1,ny
          do ii = 1, nx
             count = count + 1

             select case (idir)
             case (VAR_XC)
                index = ii
             case (VAR_YC)
                index = jj
             case (VAR_ZC)
                index = kk
             end select

             dc(count) = dd/2.d0 + dd * (index - 1) + dmin
          end do
       end do
    end do

  end subroutine ComputeCentroids1

  !------------------------------------------------------------------------
  subroutine ComputeCentroids2(nx, ny, nz, &
       dd, dmin, idir, &
       dc)
    !
    ! !DESCRIPTION:
    ! Similar to ComputeCentroids1D except returns values in a 3D array.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_XC, VAR_YC, VAR_ZC
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(in)           :: dd, dmin
    PetscInt , intent(in)           :: idir
    PetscReal, intent(out), pointer :: dc(:,:,:)
    !
    PetscReal, pointer              :: dc_1D(:)

    allocate(dc_1D(nx*ny*nz))

    call ComputeCentroids1(nx, ny, nz, &
       dd, dmin, idir, dc_1D)

    call Remap1Dto3D(nx, ny, nz, dc_1D, dc)

    deallocate(dc_1D)

  end subroutine ComputeCentroids2

  !------------------------------------------------------------------------
  subroutine ComputeInternalConnections(nx, ny, nz, &
       dx, dy, dz, idir, &
       nconn, conn_id_up, conn_id_dn, conn_dist_up, conn_dist_dn, &
       conn_area, conn_type)
    !
    ! !DESCRIPTION:
    ! Sets up connection set along 'idir' direction
    ! for a structured grid that has spatially homogeneous
    ! dx, dy, and dz
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : CONN_IN_X_DIR
    use MultiPhysicsProbConstants, only : CONN_IN_Y_DIR
    use MultiPhysicsProbConstants, only : CONN_IN_Z_DIR
    use MultiPhysicsProbConstants, only : CONN_IN_XYZ_DIR
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(in)           :: dx, dy, dz
    PetscInt , intent(in)           :: idir
    PetscInt , intent(out)          :: nconn
    PetscInt , intent(out), pointer :: conn_id_up(:), conn_id_dn(:)
    PetscReal, intent(out), pointer :: conn_dist_up(:), conn_dist_dn(:)
    PetscReal, intent(out), pointer :: conn_area(:)
    PetscInt , intent(out), pointer :: conn_type(:)
    !
    PetscInt                        :: iconn

    select case(idir)
    case (CONN_IN_X_DIR)
       nconn = (nx-1)*ny*nz
       
    case (CONN_IN_Y_DIR)
       nconn = nx*(ny-1)*nz
       
    case (CONN_IN_Z_DIR)
       nconn = nx*ny*(nz-1)
       
    case (CONN_IN_XYZ_DIR)
       nconn = &
            (nx-1)* ny   *nz    + & ! connections in x-dir
            nx   *(ny-1)*nz     + & ! connections in y-dir
            nx   * ny   *(nz-1)     ! connections in z-dir
    case default
       write(iulog,*)'Unsupported idir'
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select
    
    allocate(conn_id_up   (nconn))
    allocate(conn_id_dn   (nconn))
    allocate(conn_dist_up (nconn))
    allocate(conn_dist_dn (nconn))
    allocate(conn_area    (nconn))
    allocate(conn_type    (nconn))

    iconn = 0
    select case(idir)
    case (CONN_IN_X_DIR)
       call ComputeIntConnAlongADirection( &
            nx, ny, nz, dx, dy, dz,        &
            iconn, CONN_IN_X_DIR,          &
            conn_id_up, conn_id_dn,        &
            conn_dist_up, conn_dist_dn,    &
            conn_area, conn_type)

    case (CONN_IN_Y_DIR)
       call ComputeIntConnAlongADirection( &
            nx, ny, nz, dx, dy, dz,        &
            iconn, CONN_IN_Y_DIR,          &
            conn_id_up, conn_id_dn,        &
            conn_dist_up, conn_dist_dn,    &
            conn_area, conn_type)
       
    case (CONN_IN_Z_DIR)
       call ComputeIntConnAlongADirection( &
            nx, ny, nz, dx, dy, dz,        &
            iconn, CONN_IN_Z_DIR,          &
            conn_id_up, conn_id_dn,        &
            conn_dist_up, conn_dist_dn,    &
            conn_area, conn_type)
       
    case (CONN_IN_XYZ_DIR)
       call ComputeIntConnAlongADirection( &
            nx, ny, nz, dx, dy, dz,        &
            iconn, CONN_IN_X_DIR,          &
            conn_id_up, conn_id_dn,        &
            conn_dist_up, conn_dist_dn,    &
            conn_area, conn_type)

       call ComputeIntConnAlongADirection( &
            nx, ny, nz, dx, dy, dz,        &
            iconn, CONN_IN_Y_DIR,          &
            conn_id_up, conn_id_dn,        &
            conn_dist_up, conn_dist_dn,    &
            conn_area, conn_type)

       call ComputeIntConnAlongADirection( &
            nx, ny, nz, dx, dy, dz,        &
            iconn, CONN_IN_Z_DIR,          &
            conn_id_up, conn_id_dn,        &
            conn_dist_up, conn_dist_dn,    &
            conn_area, conn_type)

    case default
       write(iulog,*)'Unsupported idir'
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine ComputeInternalConnections

  !------------------------------------------------------------------------
  subroutine ComputeIntConnAlongADirection(nx, ny, nz, &
       dx, dy, dz, &
       iconn, idir, &
       conn_id_up, conn_id_dn, conn_dist_up, conn_dist_dn, &
       conn_area, conn_type)
    !
    ! !DESCRIPTION:
    ! Sets up connection along either x, y, or z-direction
    ! for a structured grid that has spatially homogeneous
    ! dx, dy, and dz
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : CONN_IN_X_DIR
    use MultiPhysicsProbConstants, only : CONN_IN_Y_DIR
    use MultiPhysicsProbConstants, only : CONN_IN_Z_DIR
    use MultiPhysicsProbConstants, only : CONN_HORIZONTAL
    use MultiPhysicsProbConstants, only : CONN_VERTICAL
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(in)           :: dx, dy, dz
    PetscInt , intent(inout)        :: iconn
    PetscInt , intent(in)           :: idir
    PetscInt , intent(out), pointer :: conn_id_up(:), conn_id_dn(:)
    PetscReal, intent(out), pointer :: conn_dist_up(:), conn_dist_dn(:)
    PetscReal, intent(out), pointer :: conn_area(:)
    PetscInt , intent(out), pointer :: conn_type(:)
    !
    PetscInt, pointer               :: id(:,:,:)
    PetscInt                        :: count
    PetscInt                        :: ii, jj, kk
    PetscInt                        :: ii_offset, jj_offset, kk_offset
    PetscInt                        :: ii_max, jj_max, kk_max
    PetscInt                        :: ii_up, jj_up, kk_up
    PetscInt                        :: ii_dn, jj_dn, kk_dn
    PetscReal                       :: dist, area
    PetscInt                        :: itype

    ! Compute cell ids
    call ComputeCellID(nx, ny, nz, id)
    
    ! Default value for offset & max value for all direction
    ii_offset = 0; ii_max = nx
    jj_offset = 0; jj_max = ny
    kk_offset = 0; kk_max = nz

    ! Update value for offset & max value depending on the direction
    ! of connection set
    ! For a connection set in x-dir, ii index will loop over 1:nx-1 and
    ! cell ids of connections will be between ii and ii+1
    !
    select case (idir)
    case (CONN_IN_X_DIR)
       ii_offset = 1; ii_max = nx-1

    case (CONN_IN_Y_DIR)
       jj_offset = 1; jj_max = ny-1

    case (CONN_IN_Z_DIR)
       kk_offset = 1; kk_max = nz-1

    case default
       write(iulog,*)'Unsupported idir'
       call endrun(msg=errMsg(__FILE__,__LINE__))
   end select

    do ii = 1, ii_max
       do jj = 1, jj_max
          do kk = 1, kk_max

             iconn = iconn + 1
             ii_up = ii; ii_dn = ii + ii_offset
             jj_up = jj; jj_dn = jj + jj_offset
             kk_up = kk; kk_dn = kk + kk_offset

             conn_id_up(iconn)   = id(ii_up, jj_up, kk_up)
             conn_id_dn(iconn)   = id(ii_dn, jj_dn, kk_dn)

             select case (idir)
             case (CONN_IN_X_DIR)
                dist      = dx
                area      = dy*dz
                itype = CONN_HORIZONTAL

             case (CONN_IN_Y_DIR)
                dist      = dy
                area      = dx*dz
                itype = CONN_HORIZONTAL

             case (CONN_IN_Z_DIR)
                dist      = dz
                area      = dx*dy
                itype = CONN_VERTICAL
             end select

             conn_dist_up(iconn) = 0.5d0*dist
             conn_dist_dn(iconn) = 0.5d0*dist
             conn_area(iconn)    = area
             conn_type(iconn)    = CONN_HORIZONTAL !itype

          end do
       end do
    end do

    deallocate(id)

  end subroutine ComputeIntConnAlongADirection

  !------------------------------------------------------------------------
  subroutine AllocateMemoryForBoundaryConnection( &
       nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    !
    implicit none
    !

    PetscInt   :: nconn
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscInt  , pointer     :: tmp_3d   (:,:,:)
    PetscReal , pointer     :: unit_vec (:,:)

    allocate(id_up    (nconn ))
    allocate(id_dn    (nconn ))
    allocate(dist_up  (nconn ))
    allocate(dist_dn  (nconn ))
    allocate(area     (nconn ))
    allocate(itype    (nconn ))
    allocate(unit_vec (nconn,3 ))

    unit_vec(:,:) = 0.d0

  end subroutine AllocateMemoryForBoundaryConnection

  !------------------------------------------------------------------------
  subroutine ComputeBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    !
    implicit none
    !
    PetscInt  , intent(in)  :: nx, ny, nz
    PetscReal , intent(in)  :: dx, dy, dz
    PetscInt  , intent(out) :: nconn
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscInt  , pointer     :: tmp_3d   (:,:,:)
    PetscReal , pointer     :: unit_vec (:,:)
    !
    PetscInt                :: iconn
    PetscInt                :: count
  
    nconn = 0

    if (nx > 1) nconn = nconn + ny*nz*2
    if (ny > 1) nconn = nconn + nx*nz*2
    if (nz > 1) nconn = nconn + nx*ny*2

    call AllocateMemoryForBoundaryConnection( &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
        
    count = 0
    call ComputeXBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    
    call ComputeYBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    
    call ComputeZBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    
  end subroutine ComputeBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeLeftBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    !
    implicit none
    !
    PetscInt  , intent(in)  :: nx, ny, nz
    PetscReal , intent(in)  :: dx, dy, dz
    PetscInt  , intent(out) :: nconn
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscInt  , pointer     :: tmp_3d   (:,:,:)
    PetscReal , pointer     :: unit_vec (:,:)
    !
    PetscInt                :: iconn
    PetscInt                :: count
  
    nconn = ny*nz
        
    call AllocateMemoryForBoundaryConnection( &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    count = 0
    call ComputeXBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
         count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
         include_left_boundary =.true., &
         include_right_boundary=.false.)
    
  end subroutine ComputeLeftBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeRightBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    !
    implicit none
    !
    PetscInt  , intent(in)  :: nx, ny, nz
    PetscReal , intent(in)  :: dx, dy, dz
    PetscInt  , intent(out) :: nconn
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscInt  , pointer     :: tmp_3d   (:,:,:)
    PetscReal , pointer     :: unit_vec (:,:)
    !
    PetscInt                :: iconn
    PetscInt                :: count

    nconn = ny*nz
        
    call AllocateMemoryForBoundaryConnection( &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    count = 0
    call ComputeXBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
         count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
         include_left_boundary =.false., &
         include_right_boundary=.true.)
    
  end subroutine ComputeRightBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeSouthBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    !
    implicit none
    !
    PetscInt  , intent(in)  :: nx, ny, nz
    PetscReal , intent(in)  :: dx, dy, dz
    PetscInt  , intent(out) :: nconn
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscInt  , pointer     :: tmp_3d   (:,:,:)
    PetscReal , pointer     :: unit_vec (:,:)
    !
    PetscInt                :: iconn
    PetscInt                :: count
  
    nconn = nx*nz
        
    call AllocateMemoryForBoundaryConnection( &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    count = 0
    call ComputeYBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
         count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
         include_south_boundary =.true., &
         include_north_boundary=.false.)
    
  end subroutine ComputeSouthBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeNorthBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    !
    implicit none
    !
    PetscInt  , intent(in)  :: nx, ny, nz
    PetscReal , intent(in)  :: dx, dy, dz
    PetscInt  , intent(out) :: nconn
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscInt  , pointer     :: tmp_3d   (:,:,:)
    PetscReal , pointer     :: unit_vec (:,:)
    !
    PetscInt                :: iconn
    PetscInt                :: count

    nconn = nx*nz
        
    call AllocateMemoryForBoundaryConnection( &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    count = 0
    call ComputeYBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
         count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
         include_south_boundary =.false., &
         include_north_boundary=.true.)
    
  end subroutine ComputeNorthBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeBottomBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    !
    implicit none
    !
    PetscInt  , intent(in)  :: nx, ny, nz
    PetscReal , intent(in)  :: dx, dy, dz
    PetscInt  , intent(out) :: nconn
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscInt  , pointer     :: tmp_3d   (:,:,:)
    PetscReal , pointer     :: unit_vec (:,:)
    !
    PetscInt                :: iconn
    PetscInt                :: count
  
    nconn = ny*nz
        
    call AllocateMemoryForBoundaryConnection( &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    count = 0
    call ComputeZBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
         count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
         include_bottom_boundary =.true., &
         include_top_boundary    =.false.)
    
  end subroutine ComputeBottomBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeTopBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)
    !
    implicit none
    !
    PetscInt  , intent(in)  :: nx, ny, nz
    PetscReal , intent(in)  :: dx, dy, dz
    PetscInt  , intent(out) :: nconn
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscInt  , pointer     :: tmp_3d   (:,:,:)
    PetscReal , pointer     :: unit_vec (:,:)
    !
    PetscInt                :: iconn
    PetscInt                :: count

    nconn = ny*nz
        
    call AllocateMemoryForBoundaryConnection( &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    count = 0
    call ComputeZBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
         count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
         include_bottom_boundary =.true., &
         include_top_boundary    =.false.)
    
  end subroutine ComputeTopBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeXBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
       include_left_boundary, include_right_boundary)
    !
    use MultiPhysicsProbConstants, only : CONN_HORIZONTAL
    use MultiPhysicsProbConstants, only : CONN_VERTICAL
    !
    implicit none
    !
    PetscInt  , intent(in)  :: nx, ny, nz
    PetscReal , intent(in)  :: dx, dy, dz
    PetscInt  , intent(inout) :: count
    PetscInt  , pointer     :: id_up    (:)
    PetscInt  , pointer     :: id_dn    (:)
    PetscReal , pointer     :: dist_up  (:)
    PetscReal , pointer     :: dist_dn  (:)
    PetscReal , pointer     :: area     (:)
    PetscInt  , pointer     :: itype    (:)
    PetscReal , pointer     :: unit_vec (:,:)
    PetscBool , optional    :: include_left_boundary
    PetscBool , optional    :: include_right_boundary
    !
    PetscInt                :: ii, jj, kk
    PetscInt, pointer       :: id(:,:,:)
    PetscBool               :: add_left_boundary
    PetscBool               :: add_right_boundary

    add_left_boundary  = PETSC_TRUE
    add_right_boundary = PETSC_TRUE

    if (present(include_left_boundary )) add_left_boundary  = include_left_boundary
    if (present(include_right_boundary)) add_right_boundary = include_right_boundary

    ! Compute cell ids
    call ComputeCellID(nx, ny, nz, id)

    ! X-direction
    if (nx > 1) then
       id_up(:) = 0
       
       do kk = 1,nz
          do jj = 1,ny

             ! At the beginning
             if (add_left_boundary) then
                ii = 1;
                count             = count + 1
                id_dn(count)      = id(ii,jj,kk)
                dist_up(count)    = dx/2.d0*0.d0
                dist_dn(count)    = dx/2.d0
                area(count)       = dy*dz
                unit_vec(count,1) = 1.d0
                itype(count)      = CONN_HORIZONTAL
             endif

             ! At the end
             if (add_right_boundary) then
                ii = nx
                count             = count + 1
                id_dn(count)      = id(ii,jj,kk)
                dist_up(count)    = dx/2.d0*0.d0
                dist_dn(count)    = dx/2.d0
                area(count)       = dy*dz
                unit_vec(count,1) = -1.d0
                itype(count)      = CONN_HORIZONTAL
             endif
          end do
       end do
    end if

  end subroutine ComputeXBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeYBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
       include_south_boundary, include_north_boundary)
    !
    use MultiPhysicsProbConstants, only : CONN_HORIZONTAL
    use MultiPhysicsProbConstants, only : CONN_VERTICAL
    !
    implicit none
    !
    PetscInt  , intent(in)    :: nx, ny, nz
    PetscReal , intent(in)    :: dx, dy, dz
    PetscInt  , intent(inout) :: count
    PetscInt  , pointer       :: id_up    (:)
    PetscInt  , pointer       :: id_dn    (:)
    PetscReal , pointer       :: dist_up  (:)
    PetscReal , pointer       :: dist_dn  (:)
    PetscReal , pointer       :: area     (:)
    PetscInt  , pointer       :: itype    (:)
    PetscReal , pointer       :: unit_vec (:,:)
    PetscBool , optional      :: include_south_boundary
    PetscBool , optional      :: include_north_boundary
    !
    PetscInt                  :: ii, jj, kk
    PetscInt, pointer         :: id(:,:,:)
    PetscBool                 :: add_south_boundary
    PetscBool                 :: add_north_boundary

    add_south_boundary = PETSC_TRUE
    add_north_boundary = PETSC_TRUE

    if (present(include_south_boundary)) add_south_boundary = include_south_boundary
    if (present(include_north_boundary)) add_north_boundary = include_north_boundary

    ! Compute cell ids
    call ComputeCellID(nx, ny, nz, id)

    ! Y-direction
    if (ny > 1) then
       id_up(:) = 0
       
       do kk = 1,nz
          do ii = 1,nx

             ! At the beginning
             if (add_south_boundary) then
                jj = 1;
                count             = count + 1
                id_dn(count)      = id(ii,jj,kk)
                dist_up(count)    = 0.d0
                dist_dn(count)    = dy/2.d0
                area(count)       = dx*dz
                unit_vec(count,1) = 1.d0
                itype(count)      = CONN_HORIZONTAL
             endif

             ! At the end
             if (add_north_boundary) then
                jj = ny
                count             = count + 1
                id_dn(count)      = id(ii,jj,kk)
                dist_up(count)    = 0.d0
                dist_dn(count)    = dy/2.d0
                area(count)       = dx*dz
                unit_vec(count,1) = -1.d0
                itype(count)      = CONN_HORIZONTAL
             endif
          end do
       end do
    end if

  end subroutine ComputeYBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeZBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
       count, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, &
       include_bottom_boundary, include_top_boundary)
    !
    use MultiPhysicsProbConstants, only : CONN_HORIZONTAL
    use MultiPhysicsProbConstants, only : CONN_VERTICAL
    !
    implicit none
    !
    PetscInt       , intent(in)    :: nx, ny, nz
    PetscReal      , intent(in)    :: dx, dy, dz
    PetscInt       , intent(inout) :: count
    PetscInt       , pointer       :: id_up    (:)
    PetscInt       , pointer       :: id_dn    (:)
    PetscReal      , pointer       :: dist_up  (:)
    PetscReal      , pointer       :: dist_dn  (:)
    PetscReal      , pointer       :: area     (:)
    PetscInt       , pointer       :: itype    (:)
    PetscReal      , pointer       :: unit_vec (:,:)
    PetscBool      , optional      :: include_bottom_boundary
    PetscBool      , optional      :: include_top_boundary
    !
    PetscInt                       :: ii , jj, kk
    PetscInt       , pointer       :: id(:,:,:)
    PetscBool                      :: add_bottom_boundary
    PetscBool                      :: add_top_boundary

    add_bottom_boundary = PETSC_TRUE
    add_top_boundary    = PETSC_TRUE

    if (present(include_bottom_boundary)) add_bottom_boundary = include_bottom_boundary
    if (present(include_top_boundary   )) add_top_boundary    = include_top_boundary

    ! Compute cell ids
    call ComputeCellID(nx, ny, nz, id)

    ! Z-direction
    if (nz > 1) then
       id_up(:) = 0
       
       do jj = 1,ny
          do ii = 1,nx

             ! At the beginning
             if (add_bottom_boundary) then
                kk = 1;
                count             = count + 1
                id_dn(count)      = id(ii,jj,kk)
                dist_up(count)    = 0.d0
                dist_dn(count)    = dz/2.d0
                area(count)       = dx*dy
                unit_vec(count,1) = 1.d0
                itype(count)      = CONN_VERTICAL
             endif

             ! At the end
             if (add_top_boundary) then
                kk = nz
                count             = count + 1
                id_dn(count)      = id(ii,jj,kk)
                dist_up(count)    = 0.d0
                dist_dn(count)    = dz/2.d0
                area(count)       = dx*dy
                unit_vec(count,1) = -1.d0
                itype(count)      = CONN_VERTICAL
             endif
          end do
       end do
    end if

  end subroutine ComputeZBoundaryDomainConnection

  !------------------------------------------------------------------------
  subroutine ComputeCellID(nx, ny, nz, id)
    !
    implicit none
    !
    PetscInt, intent(in)           :: nx, ny, nz
    PetscInt, intent(out), pointer :: id(:,:,:)
    !
    PetscInt :: ii, jj, kk
    PetscInt :: count
    
    allocate(id(nx,ny,nz))
      
    ! Determine ID of control volumes
    count = 0
    do kk = 1, nz
       do jj = 1, ny
          do ii = 1, nx
             count = count + 1
             id(ii,jj,kk) = count
          end do
       end do
    end do

  end subroutine ComputeCellID

  !------------------------------------------------------------------------
  subroutine Remap1Dto3D(nx, ny, nz, data_1D, data_3D)
    !
    implicit none
    !
    PetscInt , intent(in)           :: nx, ny, nz
    PetscReal, intent(out), pointer :: data_1D(:)
    PetscReal, intent(out), pointer :: data_3D(:,:,:)
    !
    PetscInt :: ii, jj, kk
    PetscInt :: count
    
    count = 0
    do kk = 1, nz
       do jj = 1, ny
          do ii = 1, nx
             count = count + 1
             data_3D(ii,jj,kk) = data_1D(count)
          end do
       end do
    end do

  end subroutine Remap1Dto3D

#endif

end module mpp_mesh_utils
