
module MeshType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Mesh data type allocation
  ! --------------------------------------------------------

#include <petsc/finclude/petsc.h>
  use petscsys
  use mpp_mesh_utils

  ! !USES:
  use mpp_varctl         , only : iulog
  use mpp_abortutils     , only : endrun
  use mpp_shr_log_mod    , only : errMsg => shr_log_errMsg
  use ConnectionSetType  , only : connection_set_list_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: mesh_type
     character (len=256) :: name         ! indentifer for the mesh

     PetscInt            :: ncells_local ! Number of cells within the mesh
     PetscInt            :: ncells_ghost ! Number of ghost cells with the mesh
     PetscInt            :: ncells_all   ! Number of local+ghost cells with the mesh
     PetscInt            :: nlev         ! Number of cells in z
     PetscInt            :: itype        ! identifier

     ! centroid
     PetscReal, pointer  :: x(:)         ! [m]
     PetscReal, pointer  :: y(:)         ! [m]
     PetscReal, pointer  :: z(:)         ! [m]

     PetscReal, pointer  :: z_m(:)       ! [m]
     PetscReal, pointer  :: z_p(:)       ! [m]

     PetscReal, pointer  :: dx(:)        ! [m]
     PetscReal, pointer  :: dy(:)        ! [m]
     PetscReal, pointer  :: dz(:)        ! [m]
     PetscReal, pointer  :: area_xy(:)   ! [m^2]
     PetscReal, pointer  :: vol(:)       ! [m^3]

     PetscBool, pointer  :: is_active(:) ! [true/false]

     PetscInt            :: orientation

     type(connection_set_list_type) :: intrn_conn_set_list
     type(connection_set_list_type) :: lateral_conn_set_list
     type(connection_set_list_type) :: conditions_conn_set_list

   contains
     procedure, public :: Init
     procedure, public :: Clean
     procedure, public :: CreateFromCLMCols
     procedure, public :: ComputeVolume             => MeshComputeVolume
     procedure, public :: CreateAndAddConnectionSet => MeshCreateAndAddConnectionSet
     procedure, public :: SetDimensions             => MeshSetDimensions
     procedure, public :: SetGeometricAttributes    => MeshSetGeometricAttributes
     procedure, public :: SetGridCellFilter         => MeshSetGridCellFilter
     procedure, public :: SetID                     => MeshSetID
     procedure, public :: SetName                   => MeshSetName
     procedure, public :: SetOrientation            => MeshSetOrientation
     procedure, public :: Copy                      => MeshCopy
     procedure         :: AllocateMemory            => MeshAllocateMemory
  end type mesh_type

  public :: MeshCreate
  public :: MeshCreateConnectionSet

  interface MeshCreate
     module procedure MeshCreate1
     module procedure MeshCreate2
  end interface MeshCreate

  interface MeshCreateConnectionSet
    module procedure MeshCreateConnectionSet1
    module procedure MeshCreateConnectionSet2
  end interface
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this)
    !
    ! !DESCRIPTION:
    ! Initializes a mesh object
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type) :: this

    this%name         = ""
    this%ncells_local = 0
    this%ncells_ghost = 0
    this%ncells_all   = 0
    this%nlev         = 0
    this%itype        = -1
    this%orientation  = -1

    nullify(this%x         )
    nullify(this%y         )
    nullify(this%z         )

    nullify(this%z_m       )
    nullify(this%z_p       )

    nullify(this%dx        )
    nullify(this%dy        )
    nullify(this%dz        )
    nullify(this%area_xy   )
    nullify(this%vol       )

    nullify(this%is_active )

    call this%intrn_conn_set_list%Init()
    call this%lateral_conn_set_list%Init()
    call this%conditions_conn_set_list%Init()

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine MeshCopy(this, inp_mesh)
    !
    ! !DESCRIPTION:
    ! Copies a mesh
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type) :: this
    class(mesh_type) :: inp_mesh

    this%name         = inp_mesh%name
    this%ncells_local = inp_mesh%ncells_local
    this%ncells_ghost = inp_mesh%ncells_ghost
    this%ncells_all   = inp_mesh%ncells_all
    this%nlev         = inp_mesh%nlev
    this%itype        = inp_mesh%itype
    this%orientation  = inp_mesh%orientation

    call this%AllocateMemory()

    this%x(:)         = inp_mesh%x(:)
    this%y(:)         = inp_mesh%y(:)
    this%z(:)         = inp_mesh%z(:)

    this%z_m(:)       = inp_mesh%z_m(:)
    this%z_p(:)       = inp_mesh%z_P(:)

    this%dx(:)        = inp_mesh%dx(:)
    this%dy(:)        = inp_mesh%dy(:)
    this%dz(:)        = inp_mesh%dz(:)

    this%area_xy(:)   = inp_mesh%area_xy(:)
    this%vol(:)       = inp_mesh%vol(:)

    this%is_active(:) = inp_mesh%is_active(:)

    call this%intrn_conn_set_list%Copy(inp_mesh%intrn_conn_set_list)
    call this%lateral_conn_set_list%Copy(inp_mesh%lateral_conn_set_list)
    call this%conditions_conn_set_list%Copy(inp_mesh%conditions_conn_set_list)

  end subroutine MeshCopy

  !------------------------------------------------------------------------
  subroutine MeshCreate1(mesh, name, x_min, y_min, z_min, &
       x_max, y_max, z_max, nx, ny, nz, internal_conn_dir)
    !
    use MultiPhysicsProbConstants , only : MESH_AGAINST_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)  , intent(inout) :: mesh
    character(len =*) , intent(in)    :: name
    PetscReal         , intent(in)    :: x_min, y_min, z_min
    PetscReal         , intent(in)    :: x_max, y_max, z_max
    PetscInt          , intent(in)    :: nx, ny, nz
    PetscInt          , intent(in)    :: internal_conn_dir
    !
    PetscReal          :: dx, dy, dz
    PetscReal, pointer :: dx_1d(:), dy_1d(:), dz_1d(:)
    PetscReal, pointer :: xc_1d(:), yc_1d(:), zc_1d(:)
    PetscReal, pointer :: area_1d(:)
    PetscInt            :: vert_nconn
    PetscInt  , pointer :: vert_conn_id_up(:)   !
    PetscInt  , pointer :: vert_conn_id_dn(:)   !
    PetscReal , pointer :: vert_conn_dist_up(:) !
    PetscReal , pointer :: vert_conn_dist_dn(:) !
    PetscReal , pointer :: vert_conn_area(:)    !
    PetscInt  , pointer :: vert_conn_type(:)    !

    dx = (x_max - x_min)/nx
    dy = (y_max - y_min)/ny
    dz = (z_max - z_min)/nz
    
    call mesh%SetName(name)
    call mesh%SetOrientation(MESH_AGAINST_GRAVITY)
    call mesh%SetID(MESH_CLM_SOIL_COL)

    call mesh%SetDimensions(nx*ny*nz, 0, nz)
    
    call ComputeXYZCentroids( &
         nx, ny, nz, &
         dx, dy, dz, &
         x_min, y_min, z_min, &
         xc_1d, yc_1d, zc_1d)

    allocate(dx_1d  (nx*ny*nz))
    allocate(dy_1d  (nx*ny*nz))
    allocate(dz_1d  (nx*ny*nz))
    allocate(area_1d(nx*ny*nz))

    dx_1d   (:) = dx
    dy_1d   (:) = dy
    dz_1d   (:) = dz
    area_1d (:) = dx*dy

    call mesh%SetGeometricAttributes(VAR_XC   , xc_1d   )
    call mesh%SetGeometricAttributes(VAR_YC   , yc_1d   )
    call mesh%SetGeometricAttributes(VAR_ZC   , zc_1d   )    
    call mesh%SetGeometricAttributes(VAR_DX   , dx_1d   )
    call mesh%SetGeometricAttributes(VAR_DY   , dy_1d   )
    call mesh%SetGeometricAttributes(VAR_DZ   , dz_1d   )
    call mesh%SetGeometricAttributes(VAR_AREA , area_1d )

    call mesh%ComputeVolume()

    vert_nconn = 0
    call ComputeInternalConnections(nx, ny, nz,      &
       dx, dy, dz, internal_conn_dir,               &
       vert_nconn, vert_conn_id_up, vert_conn_id_dn, &
       vert_conn_dist_up, vert_conn_dist_dn,         &
       vert_conn_area, vert_conn_type)

    call mesh%intrn_conn_set_list%Init()
    call mesh%lateral_conn_set_list%Init()
    call mesh%conditions_conn_set_list%Init()

    call MeshCreateAndAddConnectionSet(mesh, CONN_SET_INTERNAL,           &
         vert_nconn,  vert_conn_id_up, vert_conn_id_dn,         &
         vert_conn_dist_up, vert_conn_dist_dn,  vert_conn_area, &
         vert_conn_type)

    deallocate(xc_1d   )
    deallocate(yc_1d   )
    deallocate(zc_1d   )
    deallocate(dx_1d   )
    deallocate(dy_1d   )
    deallocate(dz_1d   )
    deallocate(area_1d )

  end subroutine MeshCreate1

  !------------------------------------------------------------------------
  subroutine MeshCreate2(mesh, name, x_max, y_max, z_max, &
       nx, ny, nz, internal_conn_dir)
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)  , intent(inout) :: mesh
    character(len =*) , intent(in)    :: name
    PetscReal         , intent(in)    :: x_max, y_max, z_max
    PetscInt          , intent(in)    :: nx, ny, nz
    PetscInt          , intent(in)    :: internal_conn_dir

    call MeshCreate(mesh, name, &
         x_min =  0.d0, y_min =  0.d0, z_min =  0.d0, &
         x_max = x_max, y_max = y_max, z_max = z_max, &
         nx    =    nx, ny    =    ny, nz    =    nz, &
         internal_conn_dir = internal_conn_dir)

  end subroutine MeshCreate2

  !------------------------------------------------------------------------
  subroutine CreateFromCLMCols(this, begg, begc, endc, &
       ugrid, grc_landunit_indices, lun_coli, lun_colf, &
       discretization_type, ncols_ghost, &
       xc_col, yc_col, zc_col, zi, dz,  &
       area_col, grid_owner, col_itype, filter)
    !
    ! !DESCRIPTION:
    ! Creates a mesh from CLM column level data structure
    !
    ! !USES:
    use mpp_varpar                  , only : nlevgrnd
    use ConnectionSetType           , only : ConnectionSetNew
    use MultiPhysicsProbConstants   , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants   , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants   , only : DISCRETIZATION_THREE_DIM
    use MultiPhysicsProbConstants   , only : DISCRETIZATION_VERTICAL_WITH_SS
    use MultiPhysicsProbConstants   , only : CONN_VERTICAL
    use MultiPhysicsProbConstants   , only : CONN_HORIZONTAL
    use ConnectionSetType           , only : connection_set_type
    use mpp_varcon                  , only : istsoil
    use UnstructuredGridType        , only : ugrid_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)                :: this
    integer            , intent(in) :: begg
    integer            , intent(in) :: begc,endc
    type(ugrid_type)   , pointer    :: ugrid
    integer            , intent(in) :: grc_landunit_indices(:,:)
    integer            , intent(in) :: lun_coli(:)
    integer            , intent(in) :: lun_colf(:)
    integer            , intent(in) :: discretization_type
    integer            , intent(in) :: ncols_ghost
    PetscReal, pointer , intent(in) :: xc_col(:)
    PetscReal, pointer , intent(in) :: yc_col(:)
    PetscReal, pointer , intent(in) :: zc_col(:)
    PetscReal, pointer , intent(in) :: area_col(:)
    PetscReal, pointer , intent(in) :: zi(:,:)
    PetscReal, pointer , intent(in) :: dz(:,:)
    PetscInt, pointer  , intent(in) :: grid_owner(:)
    PetscInt, pointer  , intent(in) :: col_itype(:)
    PetscInt, pointer  , intent(in) :: filter(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                           :: c,j,l                              !indices
    PetscInt                           :: icell
    PetscInt                           :: iconn
    PetscInt                           :: nconn
    PetscInt                           :: nconn_vert
    PetscInt                           :: nconn_horz
    PetscInt                           :: id_up, id_dn
    PetscInt                           :: first_active_hydro_col_id
    PetscInt                           :: col_id
    class(connection_set_type),pointer :: conn_set
    PetscInt                           :: g_up, g_dn, iedge
    PetscInt                           :: c_idx_up, c_idx_dn
    PetscInt                           :: l_idx_up, l_idx_dn
    PetscInt                           :: ltype, ctype
    PetscInt                           :: tmp
    PetscReal                          :: dist_x, dist_y, dist_z, dist
    PetscReal                          :: dc, dv
    PetscErrorCode                     :: ierr

    call this%Init()

    this%name           = "Soil Mesh"
    this%itype          = MESH_CLM_SOIL_COL
    this%ncells_local   = (endc - begc + 1)*nlevgrnd
    this%ncells_ghost   = ncols_ghost*nlevgrnd
    this%ncells_all     = this%ncells_local + this%ncells_ghost
    this%nlev           = nlevgrnd
    this%orientation    = MESH_ALONG_GRAVITY

    allocate(this%x(this%ncells_all         ))
    allocate(this%y(this%ncells_all         ))
    allocate(this%z(this%ncells_all         ))

    allocate(this%z_m(this%ncells_all       ))
    allocate(this%z_p(this%ncells_all       ))

    allocate(this%dx(this%ncells_all        ))
    allocate(this%dy(this%ncells_all        ))
    allocate(this%dz(this%ncells_all        ))
    allocate(this%area_xy(this%ncells_all   ))
    allocate(this%vol(this%ncells_all       ))
    allocate(this%is_active(this%ncells_all ))


    !
    ! Populate location of cell centroids
    !

    first_active_hydro_col_id = -1
    do c = begc, endc
       if (filter(c) == 1) then
          if (first_active_hydro_col_id == -1) then
             first_active_hydro_col_id = c
             exit
          endif
       endif
    enddo

    if (first_active_hydro_col_id == -1) then
       write(iulog,*)'No active soil hydrology column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    icell = 0
    do c = begc, endc + ncols_ghost
       do j = 1, this%nlev

          icell = icell + 1
          this%x (icell) = xc_col(c)
          this%y (icell) = yc_col(c)

          this%dx(icell) = 1.0d0
          this%dy(icell) = 1.0d0

          if (filter(c) == 1) then
             col_id = c
             this%is_active(icell) = PETSC_TRUE
          else
             col_id = first_active_hydro_col_id
             this%is_active(icell) = PETSC_FALSE
          endif

          this%z(icell) = -0.5d0*(zi(col_id,j-1) + zi(col_id,j)) + zc_col(c)
          this%z_m(icell) = -zi(col_id,j-1) + zc_col(c)
          this%z_p(icell) = -zi(col_id,j  ) + zc_col(c)

          this%dz      (icell) = dz(col_id,j)
          this%area_xy (icell) = area_col(c)
          this%vol     (icell) = area_col(c)*this%dz(icell)
       enddo
    enddo

    nconn_horz = 0

    ! Determine number of horizontal connections
    if (discretization_type == DISCRETIZATION_THREE_DIM        .or. &
        discretization_type == DISCRETIZATION_VERTICAL_WITH_SS ) then
       !
       ! Sets up lateral connection between columns of type 'istsoil'
       !
       ! Assumptions:
       ! - There is only ONE 'istsoil' column per landunit per grid cell.
       ! - Grid cells that are laterally connected via cellsOnCell
       !   field defined in domain netcdf file have at least ONE
       !   column of 'istsoil' type

       nconn_horz = 0

       ltype = istsoil
       ctype = istsoil

       ! Determine number of lateral connections

       do icell = 1, ugrid%ngrid_local
          do iedge = 1, ugrid%maxEdges

             if (ugrid%gridsOnGrid_local(iedge,icell) > icell) then
                g_up = icell + begg - 1
                g_dn = ugrid%gridsOnGrid_local(iedge,icell) + begg - 1

                l_idx_up = grc_landunit_indices(ltype, g_up)
                l_idx_dn = grc_landunit_indices(ltype, g_dn)

                c_idx_up = -1
                c_idx_dn = -1

                do c = lun_coli(l_idx_up), lun_colf(l_idx_up)
                   if (col_itype(c) == ctype) then
                      if (c_idx_up /= -1) then
                         write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                              'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_up
                         call endrun(msg=errMsg(__FILE__, __LINE__))
                      endif
                      c_idx_up = c
                   endif
                enddo

                do c = lun_coli(l_idx_dn), lun_colf(l_idx_dn)
                   if (col_itype(c) == ctype) then
                      if (c_idx_dn /= -1) then
                         write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                              'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_dn
                         call endrun(msg=errMsg(__FILE__, __LINE__))
                      endif
                      c_idx_dn = c
                   endif
                enddo

                if (c_idx_up > -1 .and. c_idx_dn > -1) then
                   nconn_horz = nconn_horz + 1
                else
                   write(iulog,*)'CreateFromCLMCols: No column of ctype = ', ctype, &
                        ' found between following grid cells: ',g_up,g_dn
                   call endrun(msg=errMsg(__FILE__, __LINE__))
                endif
             endif
          enddo
       enddo

       nconn_horz = nconn_horz * this%nlev

    endif

    ! Number of vertical connections
    nconn_vert = (this%nlev - 1)*(endc - begc + 1)

    ! Total number of connections
    nconn = nconn_vert + nconn_horz

    conn_set => ConnectionSetNew(nconn)

    ! Set vertical connections
    iconn = 0
    do c = begc, endc
       do j = 1, this%nlev-1
          iconn = iconn + 1

          id_up = (c-begc)*this%nlev + j
          id_dn = id_up + 1

          call conn_set%conn(iconn)%SetType(CONN_VERTICAL)
          call conn_set%conn(iconn)%SetIDUp(id_up)
          call conn_set%conn(iconn)%SetIDDn(id_dn)

          call conn_set%conn(iconn)%SetArea(area_col(c))

          call conn_set%conn(iconn)%SetDistUp(0.5d0*this%dz(id_up))
          call conn_set%conn(iconn)%SetDistDn(0.5d0*this%dz(id_dn))

          call conn_set%conn(iconn)%SetDistUnitVec(0.d0, 0.d0, -1.d0)

       enddo
    enddo

    ! If it is not a 3D problem, then internal connections are done
    if (discretization_type /= DISCRETIZATION_THREE_DIM) then
       call this%intrn_conn_set_list%AddSet(conn_set)
       nullify(conn_set)
    endif

    if (discretization_type == DISCRETIZATION_THREE_DIM        .or. &
        discretization_type == DISCRETIZATION_VERTICAL_WITH_SS ) then

       if (discretization_type == DISCRETIZATION_VERTICAL_WITH_SS) then
          ! Set up a new connection set to hold horizontal connections
          conn_set => ConnectionSetNew(nconn_horz)
          iconn = 0
       else
          ! Contine adding horizontal connections in the internal
          ! connection list
       endif

       do icell = 1, ugrid%ngrid_local

          do iedge = 1, ugrid%maxEdges
             if (ugrid%gridsOnGrid_local(iedge,icell) > icell) then
                g_up = icell + begg - 1
                g_dn = ugrid%gridsOnGrid_local(iedge,icell) + begg - 1

                l_idx_up = grc_landunit_indices(ltype, g_up)
                l_idx_dn = grc_landunit_indices(ltype, g_dn)

                c_idx_up = -1
                c_idx_dn = -1

                do c = lun_coli(l_idx_up), lun_colf(l_idx_up)
                   if (col_itype(c) == ctype) then
                      if (c_idx_up /= -1) then
                         write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                              'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_up
                         call endrun(msg=errMsg(__FILE__, __LINE__))
                      endif
                      c_idx_up = c
                   endif
                enddo

                do c = lun_coli(l_idx_dn), lun_colf(l_idx_dn)
                   if (col_itype(c) == ctype) then
                      if (c_idx_dn /= -1) then
                         write(iulog,*)'CreateFromCLMCols: More than one column found for ' // &
                              'ctype = ', ctype, ' for ltype = ', ltype, ' in grid cell ', g_dn
                         call endrun(msg=errMsg(__FILE__, __LINE__))
                      endif
                      c_idx_dn = c
                   endif
                enddo

                if (c_idx_up > -1 .and. c_idx_dn > -1) then

                   if (grid_owner(g_up) > grid_owner(g_dn)) then
                      tmp      = g_up;
                      g_up     = g_dn
                      g_dn     = tmp

                      tmp      = l_idx_up;
                      l_idx_up = l_idx_dn
                      l_idx_dn = tmp

                      tmp      = c_idx_up;
                      c_idx_up = c_idx_dn
                      c_idx_dn = tmp
                   endif

                   dc = ugrid%dcOnGrid_local(iedge, icell)
                   dv = ugrid%dvOnGrid_local(iedge, icell)

                   do j = 1, this%nlev
                      iconn = iconn + 1

                      id_up = (c_idx_up - begc)*this%nlev + j
                      id_dn = (c_idx_dn - begc)*this%nlev + j

                      call conn_set%conn(iconn)%SetType(CONN_HORIZONTAL)
                      call conn_set%conn(iconn)%SetIDUp(id_up)
                      call conn_set%conn(iconn)%SetIDDn(id_dn)

                      this%is_active(id_up) = PETSC_TRUE
                      this%is_active(id_dn) = PETSC_TRUE

                      call conn_set%conn(iconn)%SetArea(this%dz(id_up)*dv)

                      dist_x = dc
                      dist_y = 0.d0
                      dist_z = this%z(id_dn) - this%z(id_up)
                      dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0

                      call conn_set%conn(iconn)%SetDistUp(0.5d0*dist)
                      call conn_set%conn(iconn)%SetDistDn(0.5d0*dist)

                      call conn_set%conn(iconn)%SetDistUnitVec(dist_x/dist, dist_y/dist, dist_z/dist)
                   enddo
                endif ! if (c_idx_up > -1 .and. c_idx_dn > -1)
             endif ! if (ugrid%gridsOnGrid_local(iedge,icell) > icell)
          enddo
       enddo

       if (discretization_type == DISCRETIZATION_VERTICAL_WITH_SS) then
          call this%lateral_conn_set_list%AddSet(conn_set)
          nullify(conn_set)

       else if (discretization_type == DISCRETIZATION_THREE_DIM) then
          call this%intrn_conn_set_list%AddSet(conn_set)
          nullify(conn_set)
       endif

    endif

  end subroutine CreateFromCLMCols

!------------------------------------------------------------------------
  subroutine MeshCreateConnectionSet1(mesh, region_itype, conn_set, ncells_local, &
       soil_top_cell_offset, use_clm_dist_to_interface, begc, endc, z, zi)
    !
    ! !DESCRIPTION:
    ! Creates a connection set for a mesh that stores:
    ! - ID of control volumes (CVs) upwind and downwind of a connection,
    ! - Area of the face shared between upwind and downwind CVs,
    ! - Unit normal vector from up to down control volume,
    ! - Distance from cell centroid to face centroid for upwind and downwind CV.
    !
    ! !USES:
    use ConnectionSetType           , only : ConnectionSetNew
    use MultiPhysicsProbConstants   , only : SOIL_TOP_CELLS, SOIL_BOTTOM_CELLS, SOIL_CELLS
    use MultiPhysicsProbConstants   , only : SNOW_TOP_CELLS, SNOW_BOTTOM_CELLS, SSW_TOP_CELLS, ALL_CELLS
    use MultiPhysicsProbConstants   , only : MESH_ALONG_GRAVITY, MESH_AGAINST_GRAVITY
    use ConnectionSetType           , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    type(mesh_type)                   :: mesh
    PetscInt                          :: region_itype
    type(connection_set_type),pointer :: conn_set
    PetscInt, intent(out)             :: ncells_local
    PetscInt,optional                 :: soil_top_cell_offset
    PetscBool, optional               :: use_clm_dist_to_interface
    PetscInt, optional                :: begc,endc
    PetscReal, optional, pointer      :: z(:,:)
    PetscReal, optional, pointer      :: zi(:,:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                          :: c,j
    PetscInt                          :: iconn
    PetscInt                          :: nconn
    PetscInt                          :: id_up, id_dn
    PetscInt                          :: ncols
    PetscInt                          :: offset
    PetscBool                         :: use_centroid_in_dist_computation
    PetscInt                          :: num_active_conn
    PetscInt , pointer                :: active_conn_ids(:)

    ncols = mesh%ncells_local/mesh%nlev

    offset = 0
    if (present(soil_top_cell_offset)) then
       offset = soil_top_cell_offset
    end if

    use_centroid_in_dist_computation = PETSC_FALSE
    if (present(use_clm_dist_to_interface)) then
       use_centroid_in_dist_computation = use_clm_dist_to_interface

       if (use_centroid_in_dist_computation) then
          if (.not.present(begc)) then
             write(iulog,*)'use_centroid_in_dist_computation=.true. but optional argument is absent: begc'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif

          if (.not.present(endc)) then
             write(iulog,*)'use_centroid_in_dist_computation=.true. but optional argument is absent: endc'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif

          if (.not.present(z)) then
             write(iulog,*)'use_centroid_in_dist_computation=.true. but optional argument is absent: z'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif

          if (.not.present(zi)) then
             write(iulog,*)'use_centroid_in_dist_computation=.true. but optional argument is absent: zi'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif
       endif
    endif

    select case (region_itype)
    case (SOIL_BOTTOM_CELLS, SOIL_TOP_CELLS, SNOW_TOP_CELLS, SNOW_BOTTOM_CELLS, SSW_TOP_CELLS)

       nconn = ncols
       conn_set => ConnectionSetNew(nconn)
       allocate(active_conn_ids(nconn))
       num_active_conn = 0

       iconn = 0
       do c = 1, ncols
          iconn = iconn + 1

          select case (mesh%orientation)
          case (MESH_AGAINST_GRAVITY)

             select case(region_itype)
             case (SOIL_TOP_CELLS, SNOW_TOP_CELLS, SSW_TOP_CELLS)
                id_up = -1
                id_dn = mesh%nlev*c + offset

                if (mesh%is_active(id_dn)) then
                   num_active_conn = num_active_conn + 1
                   active_conn_ids(num_active_conn) = iconn
                end if

                call conn_set%conn(iconn)%SetDistUnitVec(0.d0, 0.d0, -1.d0)

             case (SOIL_BOTTOM_CELLS, SNOW_BOTTOM_CELLS)

                id_up = -1
                id_dn = mesh%nlev*(c-1) + 1 + offset

                if (mesh%is_active(id_dn)) then
                   num_active_conn = num_active_conn + 1
                   active_conn_ids(num_active_conn) = iconn
                end if

                call conn_set%conn(iconn)%SetDistUnitVec(0.d0, 0.d0, 1.d0)
             end select

          case (MESH_ALONG_GRAVITY)
             select case(region_itype)
             case (SOIL_BOTTOM_CELLS, SNOW_BOTTOM_CELLS)
                id_up = -1
                id_dn = mesh%nlev*c + offset

                if (mesh%is_active(id_dn)) then
                   num_active_conn = num_active_conn + 1
                   active_conn_ids(num_active_conn) = iconn
                end if
                call conn_set%conn(iconn)%SetDistUnitVec(0.d0, 0.d0, 1.d0)

             case (SOIL_TOP_CELLS, SNOW_TOP_CELLS, SSW_TOP_CELLS)
                id_up = -1
                id_dn = mesh%nlev*(c-1) + 1 + offset

                if (mesh%is_active(id_dn)) then
                   num_active_conn = num_active_conn + 1
                   active_conn_ids(num_active_conn) = iconn
                end if

                call conn_set%conn(iconn)%SetDistUnitVec(0.d0, 0.d0, -1.d0)
             end select

          case default
             write(iulog,*)'MeshCreateConnectionSet: Unknown mesh%orientation'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select

          call conn_set%conn(iconn)%SetIdUp(id_up)
          call conn_set%conn(iconn)%SetIDDn(id_dn)
          call conn_set%conn(iconn)%SetArea(mesh%area_xy(id_dn))
          call conn_set%conn(iconn)%SetDistUp(0.0d0)
          if (use_centroid_in_dist_computation) then
             j = 0
             call conn_set%conn(iconn)%SetDistDn(z(begc+c-1,j+1) - zi(begc+c-1,j))
          else
             call conn_set%conn(iconn)%SetDistDn(0.5d0*mesh%dz(id_dn))
          endif
       enddo

       call conn_set%SetActiveConnections(num_active_conn, active_conn_ids)
       
       deallocate(active_conn_ids)

    case (SOIL_CELLS, ALL_CELLS)

       nconn = mesh%ncells_local
       conn_set => ConnectionSetNew(nconn)
       allocate(active_conn_ids(nconn))
       num_active_conn = 0

       iconn = 0
       do c = 1,ncols
          do j = 1,mesh%nlev
             iconn = iconn + 1

             id_up = -1
             id_dn = iconn

             call conn_set%conn(iconn)%SetIDUp(id_up)
             call conn_set%conn(iconn)%SetIDDn(id_dn)

             call conn_set%conn(iconn)%SetDistUnitVec(0.d0, 0.d0, 0.d0)

             call conn_set%conn(iconn)%SetArea(mesh%area_xy(id_dn))

             call conn_set%conn(iconn)%SetDistUp(0.0d0)
             call conn_set%conn(iconn)%SetDistDn(0.0d0)

             if (mesh%is_active(id_dn)) then
                num_active_conn = num_active_conn + 1
                active_conn_ids(num_active_conn) = iconn
             end if

          enddo
       enddo

       call conn_set%SetActiveConnections(num_active_conn, active_conn_ids)
       
       deallocate(active_conn_ids)

    case default

    end select

    ncells_local = nconn

  end subroutine MeshCreateConnectionSet1

  !------------------------------------------------------------------------
  subroutine MeshCreateConnectionSet2(this, nconn, id_up, id_dn, &
       dist_up, dist_dn, area, itype, unit_vec, conn_set, skip_id_check)
    !
    ! !DESCRIPTION:
    ! Creates a connection set based on information passed
    !
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetNew
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)                  :: this
    PetscInt                          :: nconn
    PetscInt, pointer                 :: id_up(:)
    PetscInt, pointer                 :: id_dn(:)
    PetscReal, pointer                :: dist_up(:)
    PetscReal, pointer                :: dist_dn(:)
    PetscReal, pointer                :: area(:)
    PetscInt, pointer                 :: itype(:)
    PetscReal, pointer, optional      :: unit_vec(:,:)
    PetscBool, optional               :: skip_id_check
    class(connection_set_type),pointer:: conn_set
    !
    ! !LOCAL VARIABLES:
    PetscInt                          :: iconn
    PetscReal                         :: dist_x
    PetscReal                         :: dist_y
    PetscReal                         :: dist_z
    PetscReal                         :: dist
    PetscBool                         :: check
    PetscInt                          :: num_active_conn
    PetscInt , pointer                :: active_conn_ids(:)

    check = PETSC_TRUE
    if (present(skip_id_check)) check = .not.skip_id_check

    conn_set => ConnectionSetNew(nconn)
    allocate(active_conn_ids(nconn))
    num_active_conn = 0

    do iconn = 1, nconn

       if (check .and. (id_up(iconn) > this%ncells_all)) then
          write(iulog,*)'Cell id up is greater than total number of cells '
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       if (id_dn(iconn) > this%ncells_all) then
          write(iulog,*)'Cell id down is greater than total number of cells '
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       call conn_set%conn(iconn)%SetIDUp(id_up(iconn))
       call conn_set%conn(iconn)%SetIDDn(id_dn(iconn))
       call conn_set%conn(iconn)%SetArea(area(iconn))

       call conn_set%conn(iconn)%SetDistUp(dist_up(iconn))
       call conn_set%conn(iconn)%SetDistDn(dist_dn(iconn))

       if (id_up(iconn) > 0 .and. id_dn(iconn) > 0) then
          if (this%is_active(id_up(iconn)) .and. this%is_active(id_dn(iconn))) then
             num_active_conn = num_active_conn + 1
             active_conn_ids(num_active_conn) = iconn
          end if
       elseif (id_up(iconn) > 0) then
          if (this%is_active(id_up(iconn))) then
             num_active_conn = num_active_conn + 1
             active_conn_ids(num_active_conn) = iconn
          end if
       else
          if (this%is_active(id_dn(iconn))) then
             num_active_conn = num_active_conn + 1
             active_conn_ids(num_active_conn) = iconn
          end if
       end if
       if (.not.present(unit_vec)) then
          dist_x = this%x(id_dn(iconn)) - this%x(id_up(iconn))
          dist_y = this%y(id_dn(iconn)) - this%y(id_up(iconn))
          dist_z = this%z(id_dn(iconn)) - this%z(id_up(iconn))
          dist   = (dist_x**2.d0 + dist_y**2.d0 + dist_z**2.d0)**0.5d0

          call conn_set%conn(iconn)%SetDistUnitVec(dist_x/dist, dist_y/dist, dist_z/dist)
       else
          call conn_set%conn(iconn)%SetDistUnitVec(unit_vec(iconn,1), unit_vec(iconn,2), unit_vec(iconn,3))
       endif

       call conn_set%conn(iconn)%SetType(itype(iconn))
    end do

    call conn_set%SetActiveConnections(num_active_conn, active_conn_ids)
       
    deallocate(active_conn_ids)

  end subroutine MeshCreateConnectionSet2

  !------------------------------------------------------------------------
  subroutine MeshSetName(this, name)
    !
    ! !DESCRIPTION:
    ! Set name of Mesh
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)  :: this
    character(len =*) :: name

    this%name = trim(name)

  end subroutine MeshSetName

  !------------------------------------------------------------------------
  subroutine MeshSetID(this, id)
    !
    ! !DESCRIPTION:
    ! Set ID of mesh
    !
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants, only : MESH_CLM_THERMAL_SOIL_COL
    use MultiPhysicsProbConstants, only : MESH_CLM_SNOW_COL
    use MultiPhysicsProbConstants, only : MESH_CLM_SSW_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_XYLEM_COL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)  :: this
    PetscInt          :: id

    select case(id)
    case (MESH_CLM_SOIL_COL, MESH_CLM_THERMAL_SOIL_COL, &
         MESH_CLM_SNOW_COL, MESH_CLM_SSW_COL, MESH_SPAC_ROOT_COL, &
         MESH_SPAC_XYLEM_COL)
       this%itype = id

    case default
       write(iulog,*) 'Attempting to set unsupported mesh ID. id = ', id
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine MeshSetID

  !------------------------------------------------------------------------
  subroutine MeshSetOrientation(this, orientation)
    !
    ! !DESCRIPTION:
    ! Set an indentifier that mesh along/against gravity.
    !
    use MultiPhysicsProbConstants, only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants, only : MESH_AGAINST_GRAVITY
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)  :: this
    PetscInt          :: orientation

    select case(orientation)
    case (MESH_ALONG_GRAVITY, MESH_AGAINST_GRAVITY)
       this%orientation = orientation

    case default
       write(iulog,*) 'Attempting to set invalid mesh orientation. ' // &
            ' orientation = ', orientation
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine MeshSetOrientation

  !------------------------------------------------------------------------
  subroutine MeshSetDimensions(this, ncells_local, ncells_ghost, nlev)
    !
    ! !DESCRIPTION:
    ! Set mesh dimensions
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)  :: this
    PetscInt          :: ncells_local
    PetscInt          :: ncells_ghost
    PetscInt          :: nlev

    if (ncells_local <= 0) then
       write(iulog,*) 'Invalid number of local cells. ncells_local = ',ncells_local
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (ncells_ghost < 0) then
       write(iulog,*) 'Invalid number of ghost cells. ncells_ghost = ',ncells_ghost
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

      if (nlev <= 0) then
       write(iulog,*) 'Invalid number of vertical layers. nlev = ',nlev
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    this%ncells_local = ncells_local
    this%ncells_ghost = ncells_ghost
    this%ncells_all   = this%ncells_local + this%ncells_ghost
    this%nlev         = nlev

    call this%AllocateMemory()

    ! Initialize values
    this%x(:)         = 0.d0
    this%y(:)         = 0.d0
    this%z(:)         = 0.d0

    this%z_m(:)       = 0.d0
    this%z_p(:)       = 0.d0

    this%dx(:)        = 0.d0
    this%dy(:)        = 0.d0
    this%dz(:)        = 0.d0

    this%area_xy(:)   = 0.d0
    this%vol(:)       = 0.d0
    this%is_active(:) = PETSC_TRUE

  end subroutine MeshSetDimensions

  !------------------------------------------------------------------------
  subroutine MeshAllocateMemory (this)
    !
    implicit none
    !
    class(mesh_type) :: this

    ! Allocate memory
    allocate (this%x         (this%ncells_all ))
    allocate (this%y         (this%ncells_all ))
    allocate (this%z         (this%ncells_all ))

    allocate (this%z_m       (this%ncells_all ))
    allocate (this%z_p       (this%ncells_all ))

    allocate (this%dx        (this%ncells_all ))
    allocate (this%dy        (this%ncells_all ))
    allocate (this%dz        (this%ncells_all ))

    allocate (this%area_xy   (this%ncells_all ))
    allocate (this%vol       (this%ncells_all ))
    allocate (this%is_active (this%ncells_all ))

  end subroutine MeshAllocateMemory

  !------------------------------------------------------------------------
  subroutine MeshDeallocateMemory (this)
    !
    implicit none
    !
    class(mesh_type) :: this

    ! Deallocate memory
    if (associated(this%x         )) deallocate (this%x)
    if (associated(this%y         )) deallocate (this%y)
    if (associated(this%z         )) deallocate (this%z)

    if (associated(this%z_m       )) deallocate (this%z_m)
    if (associated(this%z_p       )) deallocate (this%z_p)

    if (associated(this%dx        )) deallocate (this%dx)
    if (associated(this%dy        )) deallocate (this%dy)
    if (associated(this%dx        )) deallocate (this%dz)

    if (associated(this%area_xy   )) deallocate (this%area_xy  )
    if (associated(this%vol       )) deallocate (this%vol      )
    if (associated(this%is_active )) deallocate (this%is_active)

  end subroutine MeshDeallocateMemory

  !------------------------------------------------------------------------
  subroutine MeshSetGeometricAttributes(this, var_id, values)
    !
    ! !DESCRIPTION:
    ! Set geometric attributes of mesh
    !
    use MultiPhysicsProbConstants, only : VAR_XC
    use MultiPhysicsProbConstants, only : VAR_YC
    use MultiPhysicsProbConstants, only : VAR_ZC
    use MultiPhysicsProbConstants, only : VAR_DX
    use MultiPhysicsProbConstants, only : VAR_DY
    use MultiPhysicsProbConstants, only : VAR_DZ
    use MultiPhysicsProbConstants, only : VAR_AREA
    use MultiPhysicsProbConstants, only : VAR_VOLUME
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)   :: this
    PetscInt           :: var_id
    PetscReal, pointer :: values(:)
    !
    PetscInt :: icell

    if (size(values) /= this%ncells_all) then
       write(iulog,*) 'The size of geometric attributes does not match mesh dim. ' // &
            ' this%ncells_all = ', this%ncells_all, ' size(values) = ', size(values)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_id)
    case (VAR_XC)
       do icell = 1, this%ncells_all
          this%x(icell) = values(icell)
       enddo

    case (VAR_YC)
       do icell = 1, this%ncells_all
          this%y(icell) = values(icell)
       enddo

    case (VAR_ZC)
       do icell = 1, this%ncells_all
          this%z(icell) = values(icell)
       enddo

    case (VAR_DX)
       do icell = 1, this%ncells_all
          this%dx(icell) = values(icell)
       enddo

    case (VAR_DY)
       do icell = 1, this%ncells_all
          this%dy(icell) = values(icell)
       enddo

    case (VAR_DZ)
       do icell = 1, this%ncells_all
          this%dz(icell) = values(icell)
       enddo

    case (VAR_AREA)
       do icell = 1, this%ncells_all
          this%area_xy(icell) = values(icell)
       enddo

    case (VAR_VOLUME)
       do icell = 1, this%ncells_all
          this%vol(icell) = values(icell)
       enddo

    case default
       write(iulog,*) 'Attempting to set unknown geometric attributes. ' // &
            ' var_id = ', var_id
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine MeshSetGeometricAttributes

  !------------------------------------------------------------------------

  subroutine MeshSetGridCellFilter(this, filter)
    !
    ! !DESCRIPTION:
    ! Set grid cell active/inactive based on filter
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)  :: this
    PetscInt, pointer :: filter(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt          :: icell

    if (size(filter) /= this%ncells_all) then
       write(iulog,*) 'The size of filter does not match mesh dim. ' // &
            ' this%ncells_all = ', this%ncells_all, ' size(filter) = ', size(filter)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do icell = 1, this%ncells_all
       if (filter(icell) == 1) then
          this%is_active(icell) = PETSC_TRUE
       else
          this%is_active(icell) = PETSC_FALSE
       endif
    enddo

  end subroutine MeshSetGridCellFilter

  !------------------------------------------------------------------------

  subroutine MeshComputeVolume(this)
    !
    ! !DESCRIPTION:
    ! Compute volume of grid cells based on xy area and layer thickness
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)  :: this
    PetscInt, pointer :: filter(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt          :: icell

    do icell = 1, this%ncells_all
       if (this%is_active(icell)) then
          this%vol(icell) = this%area_xy(icell) * this%dz(icell)
       endif
    enddo

  end subroutine MeshComputeVolume

  !------------------------------------------------------------------------

  subroutine MeshCreateAndAddConnectionSet(this, conn_type, nconn, id_up, id_dn, &
       dist_up, dist_dn, area, itype, unit_vec)
    !
    ! !DESCRIPTION:
    ! Creates a connection set based on information passed
    !
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_SET_CONDITIONS
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetNew
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type)                   :: this
    PetscInt                           :: conn_type
    PetscInt                           :: nconn
    PetscInt, pointer                  :: id_up(:)
    PetscInt, pointer                  :: id_dn(:)
    PetscReal, pointer                 :: dist_up(:)
    PetscReal, pointer                 :: dist_dn(:)
    PetscReal, pointer                 :: area(:)
    PetscInt, pointer                  :: itype(:)
    PetscReal, pointer, optional       :: unit_vec(:,:)
    !
    ! !LOCAL VARIABLES:
    class(connection_set_type),pointer :: conn_set
    PetscInt                           :: iconn
    PetscReal                          :: dist_x
    PetscReal                          :: dist_y
    PetscReal                          :: dist_z
    PetscReal                          :: dist

    call MeshCreateConnectionSet(this, nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec=unit_vec, conn_set=conn_set)

    select case(conn_type)
    case (CONN_SET_INTERNAL)
       call this%intrn_conn_set_list%AddSet(conn_set)
    case (CONN_SET_LATERAL)
       call this%lateral_conn_set_list%AddSet(conn_set)
    case (CONN_SET_CONDITIONS)
       call this%conditions_conn_set_list%AddSet(conn_set)
    case default
       write(iulog,*)'Unknown connection set type = ',conn_type
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine MeshCreateAndAddConnectionSet

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !
    ! !DESCRIPTION:
    ! Release memory allocated to the mesh
    !
    ! !USES:
    implicit none
    !
    ! !ARGUMENTS
    class(mesh_type) :: this

    deallocate(this%x  )
    deallocate(this%y  )
    deallocate(this%z  )

    deallocate(this%z_m)
    deallocate(this%z_p)

    deallocate(this%dx )
    deallocate(this%dy )
    deallocate(this%dz )

    call this%intrn_conn_set_list%Destroy()
    call this%lateral_conn_set_list%Destroy()
    call this%conditions_conn_set_list%Destroy()

  end subroutine Clean

#endif

end module MeshType
