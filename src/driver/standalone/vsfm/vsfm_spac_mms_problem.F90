module vsfm_spac_mms_problem

!#define USE_VG
#include <petsc/finclude/petsc.h>
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  use MultiPhysicsProbVSFM , only : mpp_vsfm_type
  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : PRESSURE_REF
  use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
  use mpp_mesh_utils

  implicit none
  
  PetscInt              :: num_xylm   , num_root   , num_soil
  PetscReal             :: x_xylm_min , y_xylm_min , z_xylm_min
  PetscReal             :: x_xylm_max , y_xylm_max , z_xylm_max
  PetscReal             :: x_root_min , y_root_min , z_root_min
  PetscReal             :: x_root_max , y_root_max , z_root_max
  PetscReal             :: x_soil_min , y_soil_min , z_soil_min
  PetscReal             :: x_soil_max , y_soil_max , z_soil_max
  PetscReal             :: dx_xylm    , dy_xylm    , dz_xylm
  PetscReal             :: dx_root    , dy_root    , dz_root
  PetscReal             :: dx_soil    , dy_soil    , dz_soil
  PetscReal , pointer   :: xc_xylm(:) , yc_xylm(:) , zc_xylm(:)
  PetscReal , pointer   :: xc_root(:) , yc_root(:) , zc_root(:)
  PetscReal , pointer   :: xc_soil(:) , yc_soil(:) , zc_soil(:)

  PetscReal , parameter :: PI                  = 4.d0 * atan (1.d0)
  PetscReal , parameter :: root_area           = 1.d0!4.d0 * PI * 1.d-7 ! [m^2]
  PetscReal , parameter :: xylm_area           = 1.d0!0.013d0           ! [m^2]
  PetscInt              :: grid_factor

  PetscReal , parameter :: vis                 = 8.904156d-4       ! Pa s
  PetscReal , parameter :: rho                 = 1000.d0           ! kg m^{-3}

  PetscReal , parameter :: root_kmax           = 1.6d-6            ! Pa
  PetscReal , parameter :: root_phi50          = -2.5d6            ! Pa
  PetscReal , parameter :: root_phi88          = -0.5d6            ! Pa
  PetscReal , parameter :: root_c1             = 1.7d6             ! Pa
  PetscReal , parameter :: root_c2             = 3.0d0             ! -
  
  PetscReal , parameter :: xylm_kmax           = 1.6d-6            ! Pa
  PetscReal , parameter :: xylm_phi50          = -2.5d6            ! Pa
  PetscReal , parameter :: xylm_phi88          = -0.5d6            ! Pa
  !PetscReal , parameter :: xylm_phis50         = -20500.d0         ! Pa
  PetscReal , parameter :: xylm_phis50         = -0.91d6           ! Pa
  PetscReal , parameter :: xylm_c1             = 1.7d6             ! Pa
  PetscReal , parameter :: xylm_c2             = 3.0d0             ! -
  PetscReal , parameter :: xylm_c3             = 12.3d0            ! -
  PetscReal , parameter :: max_pet             = 2.d-4
  
  PetscInt  , parameter :: SOIL_XC               = 001
  PetscInt  , parameter :: SOIL_YC               = 002
  PetscInt  , parameter :: SOIL_ZC               = 003
  PetscInt  , parameter :: SOIL_PERMEABILITY     = 004
  PetscInt  , parameter :: SOIL_PRESSURE_BC      = 005
  PetscInt  , parameter :: SOIL_MASS_SOURCE      = 006
  PetscInt  , parameter :: SOIL_PRESSURE         = 007
  PetscInt  , parameter :: SOIL_POROSITY         = 008
  PetscInt  , parameter :: SOIL_SATFUNC_TYPE     = 009
  PetscInt  , parameter :: SOIL_SATFUNC_ALPHA    = 010
  PetscInt  , parameter :: SOIL_SATFUNC_LAMBDA   = 011
  PetscInt  , parameter :: SOIL_RES_SAT          = 012
  PetscInt  , parameter :: SOIL_CONDUCTANCE_VAL  = 013
  PetscInt  , parameter :: SOIL_INITIAL_PRESSURE = 014
  PetscInt  , parameter :: SOIL_BC_PRESSURE      = 015
  PetscInt  , parameter :: SOIL_LIQ_SAT          = 016
  PetscInt  , parameter :: SOIL_REL_PERM         = 017

  PetscInt  , parameter :: ROOT_XC               = 101
  PetscInt  , parameter :: ROOT_YC               = 102
  PetscInt  , parameter :: ROOT_ZC               = 103
  PetscInt  , parameter :: ROOT_PERMEABILITY     = 104
  PetscInt  , parameter :: ROOT_PRESSURE_BC      = 105
  PetscInt  , parameter :: ROOT_MASS_SOURCE      = 106
  PetscInt  , parameter :: ROOT_PRESSURE         = 107
  PetscInt  , parameter :: ROOT_POROSITY         = 108
  PetscInt  , parameter :: ROOT_SATFUNC_TYPE     = 109
  PetscInt  , parameter :: ROOT_SATFUNC_ALPHA    = 110
  PetscInt  , parameter :: ROOT_SATFUNC_LAMBDA   = 111
  PetscInt  , parameter :: ROOT_RES_SAT          = 112
  PetscInt  , parameter :: ROOT_WEIBULL_C        = 113
  PetscInt  , parameter :: ROOT_WEIBULL_D        = 114
  PetscInt  , parameter :: ROOT_CONDUCTANCE_VAL  = 115
  PetscInt  , parameter :: ROOT_INITIAL_PRESSURE = 116
  PetscInt  , parameter :: ROOT_BC_PRESSURE      = 117
  PetscInt  , parameter :: ROOT_LIQ_SAT          = 118
  PetscInt  , parameter :: ROOT_REL_PERM         = 119

  PetscInt  , parameter :: XYLM_XC               = 201
  PetscInt  , parameter :: XYLM_YC               = 202
  PetscInt  , parameter :: XYLM_ZC               = 203
  PetscInt  , parameter :: XYLM_PERMEABILITY     = 204
  PetscInt  , parameter :: XYLM_PRESSURE_BC      = 205
  PetscInt  , parameter :: XYLM_MASS_SOURCE      = 206
  PetscInt  , parameter :: XYLM_PRESSURE         = 207
  PetscInt  , parameter :: XYLM_POROSITY         = 208
  PetscInt  , parameter :: XYLM_SATFUNC_TYPE     = 209
  PetscInt  , parameter :: XYLM_SATFUNC_ALPHA    = 210
  PetscInt  , parameter :: XYLM_SATFUNC_LAMBDA   = 211
  PetscInt  , parameter :: XYLM_RES_SAT          = 212
  PetscInt  , parameter :: XYLM_WEIBULL_C        = 213
  PetscInt  , parameter :: XYLM_WEIBULL_D        = 214
  PetscInt  , parameter :: XYLM_SINK_EXPONENT    = 215
  PetscInt  , parameter :: XYLM_SINK_PRESSURE    = 216
  PetscInt  , parameter :: XYLM_CONDUCTANCE_VAL  = 217
  PetscInt  , parameter :: XYLM_INITIAL_PRESSURE = 218
  PetscInt  , parameter :: XYLM_BC_PRESSURE      = 219
  PetscInt  , parameter :: XYLM_PET              = 220
  PetscInt  , parameter :: XYLM_LIQ_SAT          = 221
  PetscInt  , parameter :: XYLM_REL_PERM         = 222

  PetscReal, pointer :: soil_root_flux(:)

  Public :: run_vsfm_spac_mms_problem

  
contains

  !------------------------------------------------------------------------
  subroutine run_vsfm_spac_mms_problem()
    !
    implicit none
    !
    PetscBool           :: converged
    PetscInt            :: converged_reason
    PetscBool           :: flg
    PetscErrorCode      :: ierr
    character(len=256)  :: true_soln_fname
    character(len=256)  :: source_fname
    character(len=256)  :: liq_sat_fname
    character(len=256)  :: rel_perm_fname

    grid_factor = 2
    true_soln_fname = ''
    source_fname    = ''
    liq_sat_fname   = ''
    rel_perm_fname  = ''
    
    call PetscOptionsGetInt( PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-grid_factor', grid_factor, flg, ierr)
    CHKERRQ(ierr)

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_true_solution', true_soln_fname, flg, ierr)

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_source', source_fname, flg, ierr)

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_liq_saturation', liq_sat_fname, flg, ierr)

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_rel_permeability', rel_perm_fname, flg, ierr)

    call Init()

    ! Preform Pre-StepDT operations
    call vsfm_mpp%soe%SetDtime(1.d0)
    call vsfm_mpp%soe%PreStepDT()
    !call vsfm_mpp%soe%PreSolve()
    !call vsfm_mpp%soe%PostSolve()

    call set_boundary_conditions()
    call set_source_sink_conditions()

    call vsfm_mpp%soe%StepDT(1.d0, 1, &
         converged, converged_reason, ierr); CHKERRQ(ierr)

    if (len(trim(adjustl(true_soln_fname)))>0) then
       call save_problem_variable(vsfm_mpp, true_soln_fname, &
       SOIL_PRESSURE, ROOT_PRESSURE, XYLM_PRESSURE)
    end if

    if (len(trim(adjustl(source_fname)))>0) then
       call save_problem_variable(vsfm_mpp, source_fname, &
       SOIL_MASS_SOURCE, ROOT_MASS_SOURCE, XYLM_MASS_SOURCE)
    end if

    if (len(trim(adjustl(liq_sat_fname)))>0) then
       call save_problem_variable(vsfm_mpp, liq_sat_fname, &
       SOIL_LIQ_SAT, ROOT_LIQ_SAT, XYLM_LIQ_SAT)
    end if

    if (len(trim(adjustl(rel_perm_fname)))>0) then
       call save_problem_variable(vsfm_mpp, rel_perm_fname, &
       SOIL_REL_PERM, ROOT_REL_PERM, XYLM_REL_PERM)
    end if

  end subroutine run_vsfm_spac_mms_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use SystemOfEquationsVSFMType , only : VSFMSOEUpdateConnections
    !
    implicit none
    !

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp()

    ! 2. Add all meshes needed for the MPP
    call add_multiple_meshes()

    ! 3. Add all governing equations
    call add_multiple_goveqns()

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns()

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars()

    ! 6. Setup the MPP
    call vsfm_mpp%SetupProblem()

    !call vsfm_mpp%soe%PrintInfo()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties() 

    ! 8. Set initial conditions
    call set_initial_conditions()

    call VSFMSOEUpdateConnections(vsfm_mpp%soe, MPP_VSFM_SNES_CLM)

  end subroutine Init


  !------------------------------------------------------------------------
  subroutine initialize_mpp()
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call vsfm_mpp%Init       ()
    call vsfm_mpp%SetName    ('Variably-Saturated-Flow-Model')
    call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_multiple_meshes()
    !
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_XYLEM_COL
    !
    implicit none

    x_xylm_min =  0.d0; x_xylm_max = 10.d0;
    x_root_min = -5.d0; x_root_max =  0.d0;
    x_soil_min = -5.d0; x_soil_max =  0.d0;

    num_xylm = 20*grid_factor;
    num_root = num_xylm/2;
    num_soil = num_root;
    
    !x_xylm_min =  0.d0; x_xylm_max = 5.d0;
    !x_xylm_min = -5.d0; x_xylm_max = 0.d0;
    !num_xylm = num_root;

    dx_xylm = (x_xylm_max - x_xylm_min)/num_xylm
    dx_root = (x_root_max - x_root_min)/num_root
    dx_soil = (x_soil_max - x_soil_min)/num_soil

    dy_xylm = 1.d0
    dy_root = 1.d0
    dy_soil = 1.d0

    dz_xylm = xylm_area
    dz_root = root_area
    dz_soil = 1.d0
    
    call vsfm_mpp%SetNumMeshes(3)

    call AddMesh(imesh = 1, name = 'Soil mesh', mesh_id = MESH_CLM_SOIL_COL    , &
         nx    = num_soil   , ny    = 1       , nz    = 1                      , &
         dx    = dx_soil    , dy    = dy_soil , dz    = dz_soil                , &
         x_min = x_soil_min , y_min = 0.d0    , z_min = 0.d0                   , &
         xc_1d = xc_soil    , yc_1d = yc_soil , zc_1d = zc_soil)
    
    call AddMesh(imesh = 2, name = 'Root mesh', mesh_id = MESH_SPAC_ROOT_COL   , &
         nx    = num_root   , ny    = 1       , nz    = 1                      , &
         dx    = dx_root    , dy    = dy_root , dz    = dz_root                , &
         x_min = x_root_min , y_min = 0.d0    , z_min = 0.d0                   , &
         xc_1d = xc_root    , yc_1d = yc_root , zc_1d = zc_root)
    
    call AddMesh(imesh = 3, name = 'Xylem Mesh', mesh_id = MESH_SPAC_XYLEM_COL , &
         nx    = num_xylm   , ny    = 1       , nz    = 1                      , &
         dx    = dx_xylm    , dy    = dy_xylm , dz    = dz_xylm                , &
         x_min = x_xylm_min , y_min = 0.d0    , z_min = 0.d0                   , &
         xc_1d = xc_xylm    , yc_1d = yc_xylm , zc_1d = zc_xylm)
    
  end subroutine add_multiple_meshes

  !------------------------------------------------------------------------
  subroutine AddMesh(imesh, name, mesh_id, nx, ny, nz, dx, dy, dz, &
       x_min, y_min, z_min, xc_1d, yc_1d, zc_1d)
    !
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_AGAINST_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : VAR_VOLUME
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_IN_XYZ_DIR
    use MeshType                  , only : mesh_type
    !
    implicit none
    !
    PetscInt                   :: id
    character(len =*)          :: name
    PetscInt                   :: mesh_id
    PetscInt                   :: nx     , ny, nz
    PetscReal                  :: dx, dy, dz
    PetscReal                  :: x_min, y_min, z_min
    PetscReal        , pointer :: xc_1d(:), yc_1d(:), zc_1d(:)
    !
    class(mesh_type) , pointer :: mesh
    PetscInt                   :: imesh, ncells
    PetscInt         , pointer :: filter(:)
    PetscReal        , pointer :: dx_1d(:), dy_1d(:), dz_1d(:)
    PetscReal        , pointer :: area(:)
    PetscReal        , pointer :: volume(:)
    PetscInt                   :: nconn
    PetscInt         , pointer :: conn_id_up(:)
    PetscInt         , pointer :: conn_id_dn(:)
    PetscReal        , pointer :: conn_dist_up(:)
    PetscReal        , pointer :: conn_dist_dn(:)
    PetscReal        , pointer :: conn_area(:)
    PetscInt         , pointer :: conn_type(:)
    PetscErrorCode             :: ierr

    call ComputeXYZCentroids(nx, ny, nz, dx, dy, dz, &
         x_min, y_min, z_min, &
         xc_1d, yc_1d, zc_1d)

    ncells = nx*ny*nz
    
    allocate(filter(nx*ny*nz))
    allocate(dx_1d (nx*ny*nz))
    allocate(dy_1d (nx*ny*nz))
    allocate(dz_1d (nx*ny*nz))
    allocate(area  (nx*ny*nz))
    allocate(volume(nx*ny*nz))

    filter (:) = 1
    dx_1d  (:) = dx
    dy_1d  (:) = dy
    dz_1d  (:) = dz
    area   (:) = dy*dz
    volume (:) = dx*dy*dz

    allocate(mesh)

    call mesh%Init()
    call mesh%SetName        (name)
    call mesh%SetOrientation (MESH_AGAINST_GRAVITY)
    call mesh%SetID          (mesh_id)
    call mesh%SetDimensions  (ncells, 0, nz)

    call mesh%SetGridCellFilter      (filter)
    call mesh%SetGeometricAttributes (VAR_XC   , xc_1d)
    call mesh%SetGeometricAttributes (VAR_YC   , yc_1d)
    call mesh%SetGeometricAttributes (VAR_ZC   , zc_1d)
    call mesh%SetGeometricAttributes (VAR_DX   , dx_1d)
    call mesh%SetGeometricAttributes (VAR_DY   , dy_1d)
    call mesh%SetGeometricAttributes (VAR_DZ   , dz_1d)
    call mesh%SetGeometricAttributes (VAR_AREA , area)
    call mesh%SetGeometricAttributes (VAR_VOLUME , volume)

    call ComputeInternalConnections(nx, ny, nz, dx, dy, dz, CONN_IN_XYZ_DIR, &
         nconn, conn_id_up, conn_id_dn, conn_dist_up, conn_dist_dn, &
         conn_area, conn_type)

    call mesh%CreateAndAddConnectionSet(CONN_SET_INTERNAL, &
         nconn,  conn_id_up, conn_id_dn,          &
         conn_dist_up, conn_dist_dn,  conn_area,  &
         conn_type)

    call vsfm_mpp%AddMesh(imesh, mesh)

    call mesh%Clean()
    
    deallocate(filter       )
    deallocate(dx_1d        )
    deallocate(dy_1d        )
    deallocate(dz_1d        )
    deallocate(area         )
    deallocate(volume       )
    deallocate(conn_id_up   )
    deallocate(conn_id_dn   )
    deallocate(conn_dist_up )
    deallocate(conn_dist_dn )
    deallocate(conn_area    )
    deallocate(conn_type    )

  end subroutine AddMesh

  !------------------------------------------------------------------------

  subroutine add_multiple_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM     , only : vsfm_mpp
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_ROOT_COL
    use MultiPhysicsProbConstants, only : MESH_SPAC_XYLEM_COL
    !
    ! !ARGUMENTS
    implicit none

    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Soil' , MESH_CLM_SOIL_COL  )
    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Root' , MESH_SPAC_ROOT_COL )
    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE for Xylem', MESH_SPAC_XYLEM_COL)

    call vsfm_mpp%SetMeshesOfGoveqns()

  end subroutine add_multiple_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:

    call add_conditions_to_goveqns_for_mms()
    call add_conditions_to_goveqns_for_coupling()

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns_for_mms()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn, ii, jj, kk
    PetscInt                            :: count, nconn
    PetscInt                  , pointer :: id_up    (:)
    PetscInt                  , pointer :: id_dn    (:)
    PetscReal                 , pointer :: dist_up  (:)
    PetscReal                 , pointer :: dist_dn  (:)
    PetscReal                 , pointer :: area     (:)
    PetscInt                  , pointer :: itype    (:)
    PetscReal                 , pointer :: unit_vec (:,:)
    class(connection_set_type), pointer :: conn_set

    ! Add source-sinks in all equations for MMS
    ieqn = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn  , &
         ss_or_bc_type = COND_SS                 , &
         name          = 'Source term for MMS'   , &
         unit          = 'kg/m^3'                , &
         cond_type     = COND_MASS_RATE          , &
         region_type   = ALL_CELLS)
    
    ieqn = 2
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn  , &
         ss_or_bc_type = COND_SS                 , &
         name          = 'Source term for MMS'   , &
         unit          = 'kg/m^3'                , &
         cond_type     = COND_MASS_RATE          , &
         region_type   = ALL_CELLS)

    ieqn = 3
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn  , &
         ss_or_bc_type = COND_SS                 , &
         name          = 'Source term for MMS'   , &
         unit          = 'kg/m^3'                , &
         cond_type     = COND_MASS_RATE          , &
         region_type   = ALL_CELLS)

    ieqn = 3
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn  , &
         ss_or_bc_type = COND_SS                 , &
         name          = 'Potential mass flux'   , &
         unit          = 'kg/m^3'                , &
         cond_type     = COND_DOWNREG_MASS_RATE_FETCH2, &
         region_type   = ALL_CELLS)

    ! Add boundary conditions in all equations for MMS
    ieqn = 1
    call ComputeBoundaryDomainConnection(num_soil, 1, 1, dx_soil, dy_soil, dz_soil, &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    call MeshCreateConnectionSet(vsfm_mpp%meshes(ieqn), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddConditionInGovEqn( ieqn , &
         ss_or_bc_type = COND_BC                 , &
         name          = 'Pressure BC for MMS'   , &
         unit          = 'Pa'                    , &
         cond_type     = COND_DIRICHLET          , &
         conn_set      = conn_set)

    nullify   (conn_set )
    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

    ieqn = 2
    call ComputeLeftBoundaryDomainConnection(num_root, 1, 1, dx_root, dy_root, dz_root, &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    call MeshCreateConnectionSet(vsfm_mpp%meshes(ieqn), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddConditionInGovEqn( ieqn , &
         ss_or_bc_type = COND_BC                 , &
         name          = 'Pressure BC for MMS'   , &
         unit          = 'Pa'                    , &
         cond_type     = COND_DIRICHLET          , &
         conn_set      = conn_set)

    nullify   (conn_set )
    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

    ieqn = 3
    call ComputeRightBoundaryDomainConnection(num_xylm, 1, 1, dx_xylm, dy_xylm, dz_xylm, &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    call MeshCreateConnectionSet(vsfm_mpp%meshes(ieqn), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddConditionInGovEqn( ieqn , &
         ss_or_bc_type = COND_BC                 , &
         name          = 'Pressure BC for MMS'   , &
         unit          = 'Pa'                    , &
         cond_type     = COND_DIRICHLET          , &
         conn_set      = conn_set)

    nullify   (conn_set )
    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

  end subroutine add_conditions_to_goveqns_for_mms

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns_for_coupling()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn, ii, jj, kk
    PetscInt                            :: num_other_goveqs
    PetscInt                  , pointer :: ieqn_others(:)
    PetscInt                            :: count, nconn
    PetscInt                  , pointer :: id_up    (:)
    PetscInt                  , pointer :: id_dn    (:)
    PetscReal                 , pointer :: dist_up  (:)
    PetscReal                 , pointer :: dist_dn  (:)
    PetscReal                 , pointer :: area     (:)
    PetscInt                  , pointer :: itype    (:)
    PetscReal                 , pointer :: unit_vec (:,:)
    class(connection_set_type), pointer :: conn_set

    allocate(ieqn_others(1))

    ! Soil mass equation --> root mass equation
    ieqn             = 1
    num_other_goveqs = 1
    ieqn_others(1)   = 2

    nconn = num_root
    allocate(id_up    (nconn))
    allocate(id_dn    (nconn))
    allocate(dist_up  (nconn))
    allocate(dist_dn  (nconn))
    allocate(area     (nconn))
    allocate(itype    (nconn))
    allocate(unit_vec (nconn,3))

    do kk = 1,num_root
       id_up(kk)      = 0
       id_dn(kk)      = kk
       dist_dn(kk)    = root_area/2.d0
       dist_up(kk)    = root_area/2.d0
       area(kk)       = root_area
       unit_vec(kk,1) = 1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_HORIZONTAL
    end do       

    call MeshCreateConnectionSet(vsfm_mpp%meshes(ieqn), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( ieqn   , &
         name               = 'Root BC in soil equation'   , &
         unit               = 'Pa'                         , &
         num_other_goveqs   = num_other_goveqs             , &
         id_of_other_goveqs = ieqn_others                  , &
         conn_set           = conn_set)

    nullify   (conn_set )
    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

    ! Root mass equation --> Soil mass equation
    ieqn             = 2
    num_other_goveqs = 1
    ieqn_others(1)   = 1

    nconn = num_root
    allocate(id_up    (nconn))
    allocate(id_dn    (nconn))
    allocate(dist_up  (nconn))
    allocate(dist_dn  (nconn))
    allocate(area     (nconn))
    allocate(itype    (nconn))
    allocate(unit_vec (nconn,3))

    do kk = 1,num_root
       id_up(kk)      = 0
       id_dn(kk)      = kk
       dist_dn(kk)    = root_area/2.d0
       dist_up(kk)    = root_area/2.d0
       area(kk)       = root_area
       unit_vec(kk,1) = 1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_HORIZONTAL
    end do       

    call MeshCreateConnectionSet(vsfm_mpp%meshes(ieqn), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( ieqn   , &
         name               = 'Soil BC in soil equation'   , &
         unit               = 'Pa'                         , &
         num_other_goveqs   = num_other_goveqs             , &
         id_of_other_goveqs = ieqn_others                  , &
         conn_set           = conn_set)

    nullify   (conn_set )
    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

    ! Root mass equation --> Xylem mass equation
    ieqn             = 2
    num_other_goveqs = 1
    ieqn_others(1)   = 3

    call ComputeRightBoundaryDomainConnection(num_root, 1, 1, dx_root, dy_root, dz_root, &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    call MeshCreateConnectionSet(vsfm_mpp%meshes(ieqn), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( ieqn   , &
         name               = 'Xylem BC in root equation'   , &
         unit               = 'Pa'                         , &
         num_other_goveqs   = num_other_goveqs             , &
         id_of_other_goveqs = ieqn_others                  , &
         conn_set           = conn_set)

    nullify   (conn_set )
    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

    ! Xylem mass equation --> Root mass equation
    ieqn             = 3
    num_other_goveqs = 1
    ieqn_others(1)   = 2

    call ComputeLeftBoundaryDomainConnection(num_xylm, 1, 1, dx_xylm, dy_xylm, dz_xylm, &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    call MeshCreateConnectionSet(vsfm_mpp%meshes(ieqn), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn( ieqn   , &
         name               = 'Root BC in xylem equation'   , &
         unit               = 'Pa'                         , &
         num_other_goveqs   = num_other_goveqs             , &
         id_of_other_goveqs = ieqn_others                  , &
         conn_set           = conn_set)

    nullify   (conn_set )
    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

  end subroutine add_conditions_to_goveqns_for_coupling

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    !
    implicit none
    !
    integer           :: ieqn
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)
    integer, pointer  :: is_bc(:)

    !
    ! Allocate auxvars
    !

    ! SOIL <---> ROOT
    ieqn               = 1
    nvars_for_coupling = 1

    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    var_ids_for_coupling    (1) = VAR_PRESSURE
    goveqn_ids_for_coupling (1) = 2

    call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

    ! ROOT <---> SOIL
    ! ROOT <---> XYLEM
    ieqn               = 2
    nvars_for_coupling = 2

    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    var_ids_for_coupling    (1) = VAR_PRESSURE
    var_ids_for_coupling    (2) = VAR_PRESSURE
    goveqn_ids_for_coupling (1) = 1
    goveqn_ids_for_coupling (2) = 3

    call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

    ! XYLEM <---> ROOT
    ieqn               = 3
    nvars_for_coupling = 1

    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    var_ids_for_coupling    (1) = VAR_PRESSURE
    goveqn_ids_for_coupling (1) = 2

    call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

    call vsfm_mpp%AllocateAuxVars()
    

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties()
    !
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbVSFM , only : VSFMMPPSetDensityType
    use EOSWaterMod          , only : DENSITY_TGDPB01
    !
    PetscInt           :: eqn_id

    eqn_id = 1; call VSFMMPPSetDensityType(vsfm_mpp, eqn_id, DENSITY_TGDPB01)
    eqn_id = 2; call VSFMMPPSetDensityType(vsfm_mpp, eqn_id, DENSITY_TGDPB01)
    eqn_id = 3; call VSFMMPPSetDensityType(vsfm_mpp, eqn_id, DENSITY_TGDPB01)

    call set_material_properties_for_soil()
    call set_material_properties_for_root()
    call set_material_properties_for_xylm()

    call set_material_properties_for_soil_bc()
    call set_material_properties_for_root_bc()

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_soil()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM , only : VSFMMPPSetRelativePermeability
    use SaturationFunction   , only : SAT_FUNC_VAN_GENUCHTEN
    use EOSWaterMod          , only : DENSITY_TGDPB01
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt           :: eqn_id
    PetscInt           :: ncells
    PetscReal, pointer :: por(:)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: lambda(:)
    PetscReal, pointer :: residual_sat(:)
    PetscReal, pointer :: perm(:)
    PetscInt , pointer :: satfunc_type(:)

    eqn_id = 1
    ncells = num_soil
    
    allocate(por             (ncells))
    allocate(perm            (ncells))
    allocate(alpha           (ncells))
    allocate(lambda          (ncells))
    allocate(residual_sat    (ncells))
    allocate(satfunc_type    (ncells))
    
    
    satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN

    call set_variable_for_problem(SOIL_POROSITY      , por          )
    call set_variable_for_problem(SOIL_PERMEABILITY  , perm         )
    call set_variable_for_problem(SOIL_SATFUNC_ALPHA , alpha        )
    call set_variable_for_problem(SOIL_SATFUNC_LAMBDA, lambda       )
    call set_variable_for_problem(SOIL_RES_SAT       , residual_sat )

    call VSFMMPPSetSoilPorosity(vsfm_mpp, eqn_id, por)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, eqn_id, perm, perm, perm)

    call VSFMMPPSetSaturationFunction(vsfm_mpp, eqn_id, satfunc_type, &
         alpha, lambda, residual_sat)

    deallocate(por             )
    deallocate(alpha           )
    deallocate(lambda          )
    deallocate(residual_sat    )
    deallocate(satfunc_type    )

  end subroutine set_material_properties_for_soil

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_soil_bc()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunctionAuxVarConn
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnIntValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnRealValue
    use SaturationFunction        , only : SAT_FUNC_VAN_GENUCHTEN
    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
    use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_MANOLI_TYPE
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use MultiPhysicsProbConstants , only : VAR_FLUX_TYPE
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_TYPE
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_UP
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_DN
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt           :: eqn_id
    PetscInt           :: ncells, nconns
    PetscInt , pointer :: soil_satfunc_type  (:)
    PetscReal, pointer :: soil_alpha         (:)
    PetscReal, pointer :: soil_lambda        (:)
    PetscReal, pointer :: soil_residual_sat  (:)
    PetscReal, pointer :: soil_conductance   (:)

    PetscInt , pointer :: flux_type          (:)
    PetscInt , pointer :: conductance_type   (:)
    PetscBool, pointer :: set_upwind_auxvar  (:)
    PetscInt , pointer :: dn_satfunc_type  (:)
    PetscReal, pointer :: dn_alpha         (:)
    PetscReal, pointer :: dn_lambda        (:)
    PetscReal, pointer :: dn_residual_sat  (:)
    PetscReal, pointer :: dn_conductance   (:)
    
    eqn_id = 1
    ncells = num_soil
    nconns = num_soil + 2
    
    allocate(soil_alpha         (ncells))
    allocate(soil_lambda        (ncells))
    allocate(soil_residual_sat  (ncells))
    allocate(soil_conductance   (ncells))
    allocate(soil_satfunc_type  (ncells))

    allocate(flux_type          (nconns))
    allocate(conductance_type   (nconns))
    allocate(set_upwind_auxvar  (nconns))

    allocate(dn_alpha           (nconns))
    allocate(dn_lambda          (nconns))
    allocate(dn_residual_sat    (nconns))
    allocate(dn_conductance     (nconns))
    allocate(dn_satfunc_type    (nconns))

    ! Set values for soil and root
    soil_satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN

    call set_variable_for_problem(SOIL_SATFUNC_ALPHA   , soil_alpha         )
    call set_variable_for_problem(SOIL_SATFUNC_LAMBDA  , soil_lambda        )
    call set_variable_for_problem(SOIL_RES_SAT         , soil_residual_sat  )
    call set_variable_for_problem(SOIL_CONDUCTANCE_VAL , soil_conductance   )
    
    ! Set values for all BCs (i.e. BC for MMS + BC for soil-root coupling)
    flux_type        (1:2     ) = DARCY_FLUX_TYPE
    flux_type        (3:nconns) = CONDUCTANCE_FLUX_TYPE

    conductance_type(1:2) = 0    ; conductance_type(3:nconns) = CONDUCTANCE_MANOLI_TYPE

    dn_alpha        (1:2) = 0.d0 ; dn_alpha        (3:nconns) = soil_alpha        (1:ncells)
    dn_lambda       (1:2) = 0.d0 ; dn_lambda       (3:nconns) = soil_lambda       (1:ncells)
    dn_residual_sat (1:2) = 0.d0 ; dn_residual_sat (3:nconns) = soil_residual_sat (1:ncells)
    dn_satfunc_type (1:2) = 0    ; dn_satfunc_type (3:nconns) = soil_satfunc_type (1:ncells)
    dn_conductance  (1:2) = 0    ; dn_conductance  (3:nconns) = soil_conductance  (1:ncells)
    
    ! Set downwind values
    set_upwind_auxvar(:) = PETSC_FALSE
    call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp , &
         eqn_id, AUXVAR_CONN_BC                    , &
         set_upwind_auxvar, dn_satfunc_type, dn_alpha    , &
         dn_lambda, dn_residual_sat)
        
    ! Set connection flux type
    call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_FLUX_TYPE, flux_type)

    ! Set conductance type
    call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_CONDUCTANCE_TYPE, conductance_type)

    call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_CONDUCTANCE_DN, dn_conductance)

    deallocate(soil_alpha         )
    deallocate(soil_lambda        )
    deallocate(soil_residual_sat  )
    deallocate(soil_conductance   )
    deallocate(soil_satfunc_type  )

    deallocate(flux_type          )
    deallocate(conductance_type   )
    deallocate(set_upwind_auxvar  )

    deallocate(dn_alpha           )
    deallocate(dn_lambda          )
    deallocate(dn_residual_sat    )
    deallocate(dn_conductance     )
    deallocate(dn_satfunc_type    )

  end subroutine set_material_properties_for_soil_bc

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_root()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
    use SaturationFunction        , only : SAT_FUNC_VAN_GENUCHTEN
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt           :: eqn_id
    PetscInt           :: ncells
    PetscReal, pointer :: por(:)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: lambda(:)
    PetscReal, pointer :: residual_sat(:)
    PetscReal, pointer :: perm(:)
    PetscInt , pointer :: satfunc_type(:)
    PetscInt , pointer :: relperm_type(:)
    PetscReal, pointer :: weibull_d(:)
    PetscReal, pointer :: weibull_c(:)

    eqn_id = 2
    ncells = num_root

    allocate (por          (ncells))
    allocate (perm         (ncells))
    allocate (alpha        (ncells))
    allocate (lambda       (ncells))
    allocate (residual_sat (ncells))
    allocate (satfunc_type (ncells))
    allocate (relperm_type (ncells))
    allocate (weibull_d    (ncells))
    allocate (weibull_c    (ncells))    
    
    call set_variable_for_problem(ROOT_POROSITY      , por           )
    call set_variable_for_problem(ROOT_PERMEABILITY  , perm          )
    call set_variable_for_problem(ROOT_SATFUNC_ALPHA , alpha         )
    call set_variable_for_problem(ROOT_SATFUNC_LAMBDA, lambda        )
    call set_variable_for_problem(ROOT_RES_SAT       , residual_sat  )

    call VSFMMPPSetSoilPorosity(vsfm_mpp, eqn_id, por)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, eqn_id, perm, perm, perm)

#ifdef USE_VG
    satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN
    call VSFMMPPSetSaturationFunction(vsfm_mpp, eqn_id, satfunc_type, &
         alpha, lambda, residual_sat)
#else
    satfunc_type(:) = SAT_FUNC_FETCH2
    call VSFMMPPSetSaturationFunction(vsfm_mpp, eqn_id, satfunc_type, &
         alpha, lambda, residual_sat)
    relperm_type(:) = RELPERM_FUNC_WEIBULL
    call set_variable_for_problem(ROOT_WEIBULL_D     , weibull_d     )
    call set_variable_for_problem(ROOT_WEIBULL_C     , weibull_c     )
    call VSFMMPPSetRelativePermeability(vsfm_mpp, eqn_id, relperm_type, weibull_d, weibull_c)
#endif

    deallocate(por          )
    deallocate(perm         )
    deallocate(alpha        )
    deallocate(lambda       )
    deallocate(residual_sat )
    deallocate(satfunc_type )
    deallocate(weibull_d    )
    deallocate(weibull_c    )

  end subroutine set_material_properties_for_root

  !------------------------------------------------------------------------

  subroutine set_material_properties_for_root_bc()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunctionAuxVarConn
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeabilityAuxVarConn
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnIntValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnRealValue
    use SaturationFunction        , only : SAT_FUNC_VAN_GENUCHTEN
    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
    use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_MANOLI_TYPE
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use MultiPhysicsProbConstants , only : VAR_FLUX_TYPE
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_TYPE
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_UP
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_DN
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt           :: eqn_id
    PetscInt           :: ncells, nconns

    PetscInt , pointer :: root_satfunc_type  (:)
    PetscReal, pointer :: root_alpha         (:)
    PetscReal, pointer :: root_lambda        (:)
    PetscReal, pointer :: root_residual_sat  (:)
    PetscInt , pointer :: root_relperm_type  (:)
    PetscReal, pointer :: root_weibull_c_val (:)
    PetscReal, pointer :: root_weibull_d_val (:)
    PetscReal, pointer :: root_conductance   (:)

    PetscInt , pointer :: flux_type          (:)
    PetscInt , pointer :: conductance_type   (:)
    PetscBool, pointer :: set_upwind_auxvar  (:)

    PetscInt , pointer :: dn_satfunc_type  (:)
    PetscReal, pointer :: dn_alpha         (:)
    PetscReal, pointer :: dn_lambda        (:)
    PetscReal, pointer :: dn_residual_sat  (:)
    PetscInt , pointer :: dn_relperm_type  (:)
    PetscReal, pointer :: dn_weibull_c_val (:)
    PetscReal, pointer :: dn_weibull_d_val (:)
    PetscReal, pointer :: dn_conductance   (:)
    
    eqn_id = 2
    
    ncells = num_soil
    allocate(root_alpha         (ncells))
    allocate(root_lambda        (ncells))
    allocate(root_residual_sat  (ncells))
    allocate(root_weibull_c_val (ncells))
    allocate(root_weibull_d_val (ncells))
    allocate(root_conductance   (ncells))
    allocate(root_satfunc_type  (ncells))
    allocate(root_relperm_type  (ncells))

    nconns = num_soil + 2
    allocate(flux_type          (nconns))
    allocate(conductance_type   (nconns))
    allocate(set_upwind_auxvar  (nconns))

    allocate(dn_alpha           (nconns))
    allocate(dn_lambda          (nconns))
    allocate(dn_residual_sat    (nconns))
    allocate(dn_weibull_c_val   (nconns))
    allocate(dn_weibull_d_val   (nconns))
    allocate(dn_conductance     (nconns))
    allocate(dn_satfunc_type    (nconns))
    allocate(dn_relperm_type    (nconns))

    ! Set values for soil and root
    root_satfunc_type(:) = SAT_FUNC_FETCH2
    root_relperm_type(:) = RELPERM_FUNC_WEIBULL

    call set_variable_for_problem(ROOT_SATFUNC_ALPHA   , root_alpha         )
    call set_variable_for_problem(ROOT_SATFUNC_LAMBDA  , root_lambda        )
    call set_variable_for_problem(ROOT_RES_SAT         , root_residual_sat  )
    call set_variable_for_problem(ROOT_WEIBULL_D       , root_weibull_d_val )
    call set_variable_for_problem(ROOT_WEIBULL_C       , root_weibull_c_val )
    call set_variable_for_problem(ROOT_CONDUCTANCE_VAL , root_conductance   )

    ! Set values for all BCs (i.e. BC for MMS + BC for soil-root coupling)
    flux_type        (1:1       ) = DARCY_FLUX_TYPE
    flux_type        (2:nconns-1) = CONDUCTANCE_FLUX_TYPE
    flux_type        (nconns    ) = DARCY_FLUX_TYPE

    conductance_type(1          ) = 0
    conductance_type(2:nconns-1 ) = CONDUCTANCE_MANOLI_TYPE
    conductance_type(nconns     ) = 0

    dn_alpha        (:) = 0.d0 ; dn_alpha        (2:nconns-1) = root_alpha        (1:ncells)
    dn_lambda       (:) = 0.d0 ; dn_lambda       (2:nconns-1) = root_lambda       (1:ncells)
    dn_residual_sat (:) = 0.d0 ; dn_residual_sat (2:nconns-1) = root_residual_sat (1:ncells)
    dn_satfunc_type (:) = 0    ; dn_satfunc_type (2:nconns-1) = root_satfunc_type (1:ncells)
    dn_weibull_d_val(:) = 0    ; dn_weibull_d_val(2:nconns-1) = root_weibull_d_val(1:ncells)
    dn_weibull_c_val(:) = 0    ; dn_weibull_c_val(2:nconns-1) = root_weibull_c_val(1:ncells)
    dn_conductance  (:) = 0    ; dn_conductance  (2:nconns-1) = root_conductance  (1:ncells)
    dn_relperm_type (:) = 0    ; dn_relperm_type (2:nconns-1) = root_relperm_type (1:ncells)

    dn_alpha        (:) = root_alpha        (1)
    dn_lambda       (:) = root_lambda       (1)
    dn_residual_sat (:) = root_residual_sat (1)
    dn_satfunc_type (:) = root_satfunc_type (1)
    dn_weibull_d_val(:) = root_weibull_d_val(1)
    dn_weibull_c_val(:) = root_weibull_c_val(1)
    dn_conductance  (:) = root_conductance  (1)
    dn_relperm_type (:) = root_relperm_type (1)

    ! Set downwind values
    set_upwind_auxvar(:) = PETSC_FALSE
#ifdef USE_VG
    dn_satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN
    call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp , &
         eqn_id, AUXVAR_CONN_BC                    , &
         set_upwind_auxvar, dn_satfunc_type, dn_alpha    , &
         dn_lambda, dn_residual_sat)
#else
    call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp , &
         eqn_id, AUXVAR_CONN_BC                    , &
         set_upwind_auxvar, dn_satfunc_type, dn_alpha    , &
         dn_lambda, dn_residual_sat)

    call VSFMMPPSetRelativePermeabilityAuxVarConn(vsfm_mpp , &
         eqn_id, AUXVAR_CONN_BC                    , &
         set_upwind_auxvar, dn_relperm_type, dn_weibull_d_val, &
         dn_weibull_c_val)
#endif
    
    ! Set connection flux type
    call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_FLUX_TYPE, flux_type)

    ! Set conductance type
    call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_CONDUCTANCE_TYPE, conductance_type)

    call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_CONDUCTANCE_DN, dn_conductance)

    deallocate(root_alpha         )
    deallocate(root_lambda        )
    deallocate(root_residual_sat  )
    deallocate(root_weibull_c_val )
    deallocate(root_weibull_d_val )
    deallocate(root_conductance   )
    deallocate(root_satfunc_type  )
    deallocate(root_relperm_type  )

    deallocate(flux_type          )
    deallocate(conductance_type   )
    deallocate(set_upwind_auxvar  )

    deallocate(dn_alpha           )
    deallocate(dn_lambda          )
    deallocate(dn_residual_sat    )
    deallocate(dn_weibull_c_val   )
    deallocate(dn_weibull_d_val   )
    deallocate(dn_conductance     )
    deallocate(dn_satfunc_type    )
    deallocate(dn_relperm_type    )

  end subroutine set_material_properties_for_root_bc

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_xylm()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSourceSinkAuxVarRealValue
    use SaturationFunction        , only : SAT_FUNC_VAN_GENUCHTEN
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_EXPONENT
    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt           :: eqn_id
    PetscInt           :: ncells
    PetscReal, pointer :: por(:)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: lambda(:)
    PetscReal, pointer :: residual_sat(:)
    PetscReal, pointer :: perm(:)
    PetscInt , pointer :: satfunc_type(:)
    PetscInt , pointer :: relperm_type(:)
    PetscReal, pointer :: weibull_d(:)
    PetscReal, pointer :: weibull_c(:)
    PetscReal, pointer :: c3(:)
    PetscReal, pointer :: phis50(:)

    eqn_id = 3
    ncells = num_xylm

    allocate (por          (ncells))
    allocate (perm         (ncells))
    allocate (alpha        (ncells))
    allocate (lambda       (ncells))
    allocate (residual_sat (ncells))
    allocate (satfunc_type (ncells))
    allocate (relperm_type (ncells))
    allocate (weibull_d    (ncells))
    allocate (weibull_c    (ncells))    
    allocate (c3           (2*ncells))
    allocate (phis50       (2*ncells))
    
    call set_variable_for_problem(XYLM_POROSITY       , por          )
    call set_variable_for_problem(XYLM_PERMEABILITY   , perm         )
    call set_variable_for_problem(XYLM_SATFUNC_ALPHA  , alpha        )
    call set_variable_for_problem(XYLM_SATFUNC_LAMBDA , lambda       )
    call set_variable_for_problem(XYLM_RES_SAT        , residual_sat )
    call set_variable_for_problem(XYLM_SINK_EXPONENT  , c3           )
    call set_variable_for_problem(XYLM_SINK_PRESSURE  , phis50       )

    call VSFMMPPSetSoilPorosity(vsfm_mpp, eqn_id, por)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, eqn_id, perm, perm, perm)

#ifdef USE_VG
    satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN
    call VSFMMPPSetSaturationFunction(vsfm_mpp, eqn_id, satfunc_type, &
         alpha, lambda, residual_sat)
#else
    satfunc_type(:) = SAT_FUNC_FETCH2
    call VSFMMPPSetSaturationFunction(vsfm_mpp, eqn_id, satfunc_type, &
         alpha, lambda, residual_sat)
    relperm_type(:) = RELPERM_FUNC_WEIBULL
    call set_variable_for_problem(XYLM_WEIBULL_D      , weibull_d    )
    call set_variable_for_problem(XYLM_WEIBULL_C      , weibull_c    )
    call VSFMMPPSetRelativePermeability(vsfm_mpp, eqn_id, relperm_type, weibull_d, weibull_c)
#endif

    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, eqn_id, &
         VAR_POT_MASS_SINK_EXPONENT, c3)

    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, eqn_id, &
         VAR_POT_MASS_SINK_PRESSURE, phis50)

    deallocate(por          )
    deallocate(perm         )
    deallocate(alpha        )
    deallocate(lambda       )
    deallocate(residual_sat )
    deallocate(satfunc_type )
    deallocate(weibull_d    )
    deallocate(weibull_c    )
    deallocate(c3           )
    deallocate(phis50       )

  end subroutine set_material_properties_for_xylm

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    !
    implicit none
    !
    PetscReal          :: theta
    PetscReal          :: Se
    PetscInt           :: ii
    PetscReal, pointer :: pressure_ic(:)
    PetscReal, pointer :: soil_ic(:),root_ic(:),xylm_ic(:)
    PetscErrorCode     :: ierr

    allocate(soil_ic(num_soil))
    allocate(root_ic(num_root))
    allocate(xylm_ic(num_xylm))

    call set_variable_for_problem(SOIL_INITIAL_PRESSURE, soil_ic)
    call set_variable_for_problem(ROOT_INITIAL_PRESSURE, root_ic)
    call set_variable_for_problem(XYLM_INITIAL_PRESSURE, xylm_ic)

    call VecGetArrayF90(vsfm_mpp%soe%solver%soln, pressure_ic, ierr); CHKERRQ(ierr)

    pressure_ic(1                   :num_soil                   ) = soil_ic(:)
    pressure_ic(1+num_soil          :num_soil+num_root          ) = root_ic(:)
    pressure_ic(1+num_soil+num_root :num_soil+num_root+num_xylm ) = xylm_ic(:)

    call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, pressure_ic, ierr); CHKERRQ(ierr)

    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)
    
    deallocate(soil_ic)
    deallocate(root_ic)
    deallocate(xylm_ic)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_boundary_conditions()
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsVSFMType           , only : sysofeqns_vsfm_type
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use MultiPhysicsProbConstants           , only : AUXVAR_BC
    use MultiPhysicsProbConstants           , only : VAR_BC_SS_CONDITION
    use GoverningEquationBaseType           , only : goveqn_base_type
    implicit none
    !
    PetscReal                  , pointer :: soil_bc(:)
    PetscReal                  , pointer :: root_bc(:)
    PetscReal                  , pointer :: xylm_bc(:)
    class(sysofeqns_base_type) , pointer :: base_soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: nconn
    PetscInt                             :: soe_auxvar_id

    nconn = 0

    allocate (soil_bc(2))
    allocate (root_bc(1))
    allocate (xylm_bc(1))

    call set_variable_for_problem(SOIL_BC_PRESSURE, soil_bc)
    call set_variable_for_problem(ROOT_BC_PRESSURE, root_bc)
    call set_variable_for_problem(XYLM_BC_PRESSURE, xylm_bc)

    base_soe => vsfm_mpp%soe

    select type(base_soe)
    class is (sysofeqns_vsfm_type)

       call base_soe%SetDataFromCLM(AUXVAR_BC, VAR_BC_SS_CONDITION, 1, soil_bc)
       call base_soe%SetDataFromCLM(AUXVAR_BC, VAR_BC_SS_CONDITION, 2, root_bc)
       call base_soe%SetDataFromCLM(AUXVAR_BC, VAR_BC_SS_CONDITION, 3, xylm_bc)

    end select

    deallocate (soil_bc)
    deallocate (root_bc)
    deallocate (xylm_bc)

  end subroutine set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_source_sink_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsVSFMType , only : sysofeqns_vsfm_type
    !
    PetscInt                             :: soe_auxvar_id
    PetscReal                  , pointer :: soil_ss(:)
    PetscReal                  , pointer :: root_ss(:)
    PetscReal                  , pointer :: xylm_ss(:)
    class(sysofeqns_base_type) , pointer :: base_soe
    PetscErrorCode                       :: ierr

    allocate(soil_ss(num_soil))
    allocate(root_ss(num_root))
    allocate(xylm_ss(num_xylm))

    call set_variable_for_problem(SOIL_MASS_SOURCE, soil_ss)
    call set_variable_for_problem(ROOT_MASS_SOURCE, root_ss)
    call set_variable_for_problem(XYLM_MASS_SOURCE, xylm_ss)

    base_soe => vsfm_mpp%soe

    select type(base_soe)
    class is (sysofeqns_vsfm_type)
       
       call base_soe%SetDataFromCLM(AUXVAR_SS, VAR_BC_SS_CONDITION, 1, soil_ss)
       call base_soe%SetDataFromCLM(AUXVAR_SS, VAR_BC_SS_CONDITION, 2, root_ss)
       call base_soe%SetDataFromCLM(AUXVAR_SS, VAR_BC_SS_CONDITION, 3, xylm_ss)

       call set_variable_for_problem(XYLM_PET, xylm_ss)
       call base_soe%SetDataFromCLM(AUXVAR_SS, VAR_BC_SS_CONDITION, 4, xylm_ss)

    end select

    deallocate(soil_ss)
    deallocate(root_ss)
    deallocate(xylm_ss)
    
  end subroutine set_source_sink_conditions

  !------------------------------------------------------------------------
   subroutine compute_soil_pressure_or_deriv(x, val, dval_dx, d2val_dx2)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     PetscReal, intent(out), optional :: d2val_dx2
     !
     PetscReal, parameter             :: a0 =  1000.d0
     PetscReal, parameter             :: a1 = -20000.d0
     PetscReal                        :: num, den

     num     = x          - x_soil_min
     den     = x_soil_max - x_soil_min

     if (present(val      )) val       =  a0                 *sin(num/den*PI) + a1 + PRESSURE_REF
     if (present(dval_dx  )) dval_dx   =  a0*PI/den          *cos(num/den*PI)
     if (present(d2val_dx2)) d2val_dx2 = -a0*((PI/den)**2.d0)*sin(num/den*PI)
     
   end subroutine compute_soil_pressure_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_soil_permeability_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: p0 = 1.d-11 ! [m^2]
     PetscReal            :: num, den

     num     = x          - x_soil_min
     den     = x_soil_max - x_soil_min

     if (present(val    )) val     = p0*2.d0
     if (present(dval_dx)) dval_dx = 0.d0

   end subroutine compute_soil_permeability_or_deriv

   !------------------------------------------------------------------------
   subroutine compute_soil_alpha_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 1.d0/4000.d0 ! [Pa^{-1}]

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_soil_alpha_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_soil_lambda_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 0.5d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_soil_lambda_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_soil_residualsat_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 0.0d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_soil_residualsat_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_root_pressure_or_deriv(x, val, dval_dx, d2val_dx2)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     PetscReal, intent(out), optional :: d2val_dx2
     !
     PetscReal, parameter             :: a0 =  2000.d0
     PetscReal, parameter             :: a1 = -25000.d0
     PetscReal                        :: num, den

     num     = x          - x_xylm_min
     den     = x_xylm_max - x_xylm_min

     if (present(val      )) val       =  a0                      *sin(num/den*PI*2.d0) + a1 + PRESSURE_REF
     if (present(dval_dx  )) dval_dx   =  a0*(2.d0*PI/den)        *cos(num/den*PI*2.d0)
     if (present(d2val_dx2)) d2val_dx2 = -a0*((2.d0*PI/den)**2.d0)*sin(num/den*PI*2.d0)

     num     = x          - x_root_min
     den     = x_root_max - x_root_min
     if (present(val      )) val       =  -a0                 *sin(num/den*PI) + a1 + PRESSURE_REF
     if (present(dval_dx  )) dval_dx   =  -a0*(PI/den)        *cos(num/den*PI)
     if (present(d2val_dx2)) d2val_dx2 = a0*((PI/den)**2.d0)*sin(num/den*PI)

   end subroutine compute_root_pressure_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_root_alpha_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  root_phi88
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_root_alpha_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_root_lambda_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  root_phi50
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_root_lambda_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_root_weibullD_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  root_c1
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_root_weibullD_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_root_weibullC_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  root_c2
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_root_weibullC_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_root_permeability_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  root_kmax*vis/1000.d0*1.125d0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_root_permeability_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_root_residualsat_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  0.d0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_root_residualsat_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_pressure_or_deriv(x, val, dval_dx, d2val_dx2)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     PetscReal, intent(out), optional :: d2val_dx2
     !
     PetscReal, parameter             :: a0 =  2000.d0
     PetscReal, parameter             :: a1 = -25000.d0
     PetscReal                        :: num, den

     num     = x          - x_xylm_min
     den     = x_xylm_max - x_xylm_min

     if (present(val      )) val       =  a0                      *sin(num/den*PI*2.d0) + a1 + PRESSURE_REF
     if (present(dval_dx  )) dval_dx   =  a0*(2.d0*PI/den)        *cos(num/den*PI*2.d0)
     if (present(d2val_dx2)) d2val_dx2 = -a0*((2.d0*PI/den)**2.d0)*sin(num/den*PI*2.d0)
     
     !if (present(val      )) val       =  a0                 *sin(num/den*PI) + a1 + PRESSURE_REF
     !if (present(dval_dx  )) dval_dx   =  a0*(PI/den)        *cos(num/den*PI)
     !if (present(d2val_dx2)) d2val_dx2 =  -a0*((PI/den)**2.d0)*sin(num/den*PI)

   end subroutine compute_xylm_pressure_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_alpha_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  xylm_phi88
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_xylm_alpha_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_lambda_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  xylm_phi50
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_xylm_lambda_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_weibullD_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  xylm_c1
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_xylm_weibullD_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_weibullC_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  xylm_c2
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_xylm_weibullC_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_permeability_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  xylm_kmax*vis/1000.d0*1.125d0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_xylm_permeability_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_residualsat_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  0.d0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_xylm_residualsat_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_sink_exponent_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  xylm_c3
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_xylm_sink_exponent_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_xylm_sink_pressure_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx

     if (present(val    )) val     =  xylm_phis50
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_xylm_sink_pressure_or_deriv

   !------------------------------------------------------------------------

   subroutine set_variable_for_problem(data_type, data_1D)
    !
    ! !DESCRIPTION:
    !
    use EOSWaterMod              , only : DENSITY_TGDPB01,DENSITY_CONSTANT
    use EOSWaterMod              , only : Density, Viscosity
    use MultiPhysicsProbConstants, only : FMWH2O
    use SaturationFunction       , only : saturation_params_type
    use SaturationFunction       , only : SatFunc_Set_VG
    use SaturationFunction       , only : SatFunc_Set_FETCH2
    use SaturationFunction       , only : SatFunc_Set_Weibull_RelPerm
    use SaturationFunction       , only : SatFunc_PressToSat
    use SaturationFunction       , only : SatFunc_PressToRelPerm
    use RichardsODEPressureConnAuxType, only : rich_ode_pres_conn_auxvar_type
    use MultiPhysicsProbConstants , only : CONDUCTANCE_MANOLI_TYPE
    !
    implicit none
    !
    PetscInt                     :: data_type
    PetscReal , pointer          :: data_1D(:)
    PetscInt                     :: ii
    PetscReal                    :: xx
    PetscReal                    :: num, den
    PetscReal                    :: a0, a1
    PetscReal                    :: P     , dP_dx, d2P_dx2
    PetscReal                    :: p0    , dp0_dx
    PetscReal                    :: k     , dk_dx
    PetscReal                    :: ds0   , ds0_dx
    PetscReal                    :: dsr   , dsr_dx
    PetscReal                    :: kr    , dkr_dx, dkr_dse
    PetscReal                    :: rho   , drho_dP, d2rho_dP2, drho_dx, d2rho_dx2
    PetscReal                    :: mu    , dmu_dP
    PetscReal                    :: dse_dP, dse_dp0
    PetscReal                    :: m
    PetscReal , parameter        :: g = GRAVITY_CONSTANT
    PetscReal                    :: dden_dT, dmu_dT
    PetscReal                    :: sat_res, se
    PetscReal                    :: dkr_dP
    PetscReal                    :: P_bc, p0_bc, m_bc, sat_res_bc, cond_bc, rho_bc
    PetscReal                    :: cond,area
    PetscReal                    :: weibull_c, weibull_d
    type(saturation_params_type) :: satParam
    type (rich_ode_pres_conn_auxvar_type) :: auxvar_conn

    select case(data_type)
    case (SOIL_POROSITY)
       data_1D(:) = 0.d0

    case (SOIL_PERMEABILITY)
       do ii = 1, num_soil
          xx = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil
          call compute_soil_permeability_or_deriv(xx, val=data_1D(ii))
       end do

    case (SOIL_SATFUNC_ALPHA)
       do ii = 1, num_soil
          xx = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil
          call compute_soil_alpha_or_deriv(xx, val=data_1D(ii))
       end do

    case (SOIL_SATFUNC_LAMBDA)
       do ii = 1, num_soil
          xx = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil
          call compute_soil_lambda_or_deriv(xx, val=data_1D(ii))
       end do

    case (SOIL_RES_SAT)
       do ii = 1, num_soil
          xx = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil
          call compute_soil_residualsat_or_deriv(xx, val=data_1D(ii))
       end do

    case (SOIL_CONDUCTANCE_VAL)
       data_1D(:) = 1.d-11

    case (SOIL_INITIAL_PRESSURE)
       P = 0.d0
       do ii = 1, num_soil
          xx = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil
          call compute_soil_pressure_or_deriv(xx, val=data_1D(ii))
          P = P + 1.d0/num_soil*data_1D(ii)
       end do
       data_1D(:) = P

    case (SOIL_PRESSURE)
       do ii = 1, num_soil
          xx = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil
          call compute_soil_pressure_or_deriv(xx, val=data_1D(ii))
       end do

    case (SOIL_BC_PRESSURE)
       ii = 1
       xx = x_soil_min
       call compute_soil_pressure_or_deriv(xx, val=data_1D(ii))

       ii = 2
       xx = x_soil_max
       call compute_soil_pressure_or_deriv(xx, val=data_1D(ii))

    case (SOIL_MASS_SOURCE)
       allocate(soil_root_flux(num_soil))
       do ii = 1, num_soil

          xx      = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil

          call compute_soil_permeability_or_deriv (xx, val=k , dval_dx=dk_dx                   )
          call compute_soil_alpha_or_deriv        (xx, val=p0, dval_dx=dp0_dx                  )
          call compute_soil_lambda_or_deriv       (xx, val=m                                   )
          call compute_soil_residualsat_or_deriv  (xx, val=sat_res                             )
          call compute_soil_pressure_or_deriv     (xx, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2 )
          call compute_soil_residualsat_or_deriv  (xx, sat_res)
          
          call Viscosity(P, 298.15d0, mu, dmu_dP, dmu_dT)
          call Density(P, 298.15d0, DENSITY_TGDPB01,  rho, drho_dP, dden_dT)

          rho       = rho     * FMWH2O
          drho_dP   = drho_dP * FMWH2O
          d2rho_dP2 = 0.d0

          call SatFunc_Set_VG(satParam, sat_res, p0, m)
          call SatFunc_PressToSat(satParam, P, se, dse_dP)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)

          call compute_root_pressure_or_deriv     (xx, val=P_bc)
#ifdef USE_VG
          p0_bc=p0; m_bc=m; sat_res_bc=sat_res;
#else
          call compute_root_alpha_or_deriv        (xx, val=p0_bc)
          call compute_root_lambda_or_deriv       (xx, val=m_bc)
          call compute_root_residualsat_or_deriv  (xx, val=sat_res_bc)
#endif
          cond = 1.d-11
          cond_bc = 2.d-11

          call auxvar_conn%Init()

          auxvar_conn%conductance_type = CONDUCTANCE_MANOLI_TYPE
          auxvar_conn%pressure_dn      = P
          auxvar_conn%pressure_up      = P_bc
          auxvar_conn%conductance_dn   = cond
          auxvar_conn%conductance_up   = cond_bc

#ifdef USE_VG
          call SatFunc_Set_VG(auxvar_conn%satParams_dn, sat_res   , p0   , m)
          call SatFunc_Set_VG(auxvar_conn%satParams_up, sat_res_bc, p0_bc, m_bc)
#else
          call SatFunc_Set_VG    (auxvar_conn%satParams_dn, sat_res   , p0   , m)
          call SatFunc_Set_FETCH2(auxvar_conn%satParams_up, p0_bc, m_bc)
          call compute_root_weibullD_or_deriv (xx, weibull_d)
          call compute_root_weibullC_or_deriv (xx, weibull_c)
          call SatFunc_Set_Weibull_RelPerm(auxvar_conn%satParams_up, weibull_d, weibull_c)
#endif

          call auxvar_conn%AuxVarCompute()
          call Density(P, 298.15d0, DENSITY_TGDPB01,  rho_bc, drho_dP, dden_dT)
          rho_bc = rho_bc*FMWH2O
          area = 1.d0
          soil_root_flux(ii) = -(0.5d0*rho + 0.5d0*rho_bc)*auxvar_conn%krg*(P_bc - P)*area
          
          call Density(P, 298.15d0, DENSITY_TGDPB01,  rho, drho_dP, dden_dT)
          rho = rho*FMWH2O

          dkr_dse   = &
               0.5d0 * se**(-0.5d0) *( 1.d0 - (1.d0 - se **(1.d0/m))**m)**2.d0 + &
               se**(0.5d0) * 2.d0   *( 1.d0 - (1.d0 - se **(1.d0/m))**m) * (1.d0 - se**(1.d0/m)) * se**(1.d0/m - 1.d0)
          dse_dp0   = 0.d0
          
          dkr_dx    = dkr_dP * dP_dx + dkr_dse * dse_dp0 * dp0_dx
          drho_dx   = drho_dP * dP_dx
          d2rho_dx2 = d2rho_dP2 * dP_dx + drho_dP * d2P_dx2

          data_1D(ii) = &
               -((k*kr/mu)*drho_dx + (rho*kr/mu)*dk_dx + (rho*k/mu)*dkr_dx)*(dP_dx) &
               -(rho*k*kr/mu)*(d2P_dx2)
          data_1D(ii) = data_1D(ii)*dx_soil
          data_1D(ii) = data_1D(ii) + soil_root_flux(ii)
       end do

    case (SOIL_LIQ_SAT)
       do ii = 1, num_soil
          xx      = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil

          call compute_soil_alpha_or_deriv        (xx, val=p0)
          call compute_soil_lambda_or_deriv       (xx, val=m )
          call compute_soil_residualsat_or_deriv  (xx, val=sat_res)
          call compute_soil_pressure_or_deriv     (xx, val=P)

          call SatFunc_Set_VG(satParam, sat_res, p0, m)
          call SatFunc_PressToSat(satParam, P, se, dse_dP)
          data_1D(ii) =  se
       end do

    case (SOIL_REL_PERM)
       do ii = 1, num_soil

          xx      = x_soil_min + dx_soil/2.d0 + (ii-1)*dx_soil

          call compute_soil_alpha_or_deriv        (xx, val=p0)
          call compute_soil_lambda_or_deriv       (xx, val=m )
          call compute_soil_residualsat_or_deriv  (xx, val=sat_res)
          call compute_soil_pressure_or_deriv     (xx, val=P)

          call SatFunc_Set_VG(satParam, sat_res, p0, m)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)
          data_1D(ii) =  kr
       end do

    case (ROOT_POROSITY)
       data_1D(:) = 0.d0

    case (ROOT_PERMEABILITY)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
#ifdef USE_VG
          call compute_soil_permeability_or_deriv (xx, data_1D(ii))
#else
          call compute_root_permeability_or_deriv (xx, data_1D(ii))
#endif
       end do

    case (ROOT_SATFUNC_ALPHA)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
#ifdef USE_VG
          call compute_soil_alpha_or_deriv (xx, data_1D(ii))
#else
          call compute_root_alpha_or_deriv (xx, data_1D(ii))
#endif
       end do

    case (ROOT_SATFUNC_LAMBDA)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
#ifdef USE_VG
          call compute_soil_lambda_or_deriv (xx, data_1D(ii))
#else
          call compute_root_lambda_or_deriv (xx, data_1D(ii))
#endif
       end do

    case (ROOT_RES_SAT)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
#ifdef USE_VG
          call compute_soil_residualsat_or_deriv (xx, data_1D(ii))
#else
          call compute_root_residualsat_or_deriv (xx, data_1D(ii))
#endif
       end do

    case (ROOT_WEIBULL_D)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
          call compute_root_weibullD_or_deriv (xx, data_1D(ii))
       end do

    case (ROOT_WEIBULL_C)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
          call compute_root_weibullC_or_deriv (xx, data_1D(ii))
       end do

    case (ROOT_CONDUCTANCE_VAL)
       data_1D(:) = 2.d-11

    case (ROOT_INITIAL_PRESSURE)
       P = 0.d0
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
          call compute_root_pressure_or_deriv (xx, data_1D(ii))
          P = P + 1.d0/num_root*data_1D(ii)
       end do
       data_1D(:) = P

    case (ROOT_PRESSURE)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
          call compute_root_pressure_or_deriv (xx, data_1D(ii))
       end do

    case (ROOT_BC_PRESSURE)
       xx = x_root_min
       ii = 1
       call compute_root_pressure_or_deriv(xx, data_1D(ii))

    case (ROOT_MASS_SOURCE)
       do ii = 1, num_root

          xx      = x_root_min + dx_root/2.d0 + (ii-1)*dx_root
          
          call compute_root_permeability_or_deriv (xx, val=k , dval_dx=dk_dx                   )
          call compute_root_alpha_or_deriv        (xx, val=p0, dval_dx=dp0_dx                  )
          call compute_root_lambda_or_deriv       (xx, val=m                                   )
          call compute_root_residualsat_or_deriv  (xx, val=sat_res                             )
          call compute_root_pressure_or_deriv     (xx, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2 )
          call compute_root_residualsat_or_deriv  (xx, sat_res)

#ifdef USE_VG
          call compute_soil_permeability_or_deriv (xx, val=k , dval_dx=dk_dx                   )
          call compute_soil_alpha_or_deriv        (xx, val=p0, dval_dx=dp0_dx                  )
          call compute_soil_lambda_or_deriv       (xx, val=m                                   )
          call compute_soil_residualsat_or_deriv  (xx, val=sat_res                             )
          !call compute_soil_pressure_or_deriv     (xx, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2 )
          call compute_soil_residualsat_or_deriv  (xx, sat_res)
#endif

          call Viscosity(P, 298.15d0, mu, dmu_dP, dmu_dT)
          call Density(P, 298.15d0, DENSITY_TGDPB01,  rho, drho_dP, dden_dT)

          rho     = rho     * FMWH2O
          drho_dP = drho_dP * FMWH2O
          d2rho_dP2 = 0.d0

#ifdef USE_VG
          call SatFunc_Set_VG(satParam, sat_res, p0, m)
#else
          call SatFunc_Set_FETCH2(satParam, p0, m)
          call compute_root_weibullD_or_deriv(xx, weibull_d)
          call compute_root_weibullC_or_deriv(xx, weibull_c)
          call SatFunc_Set_Weibull_RelPerm(satParam, weibull_d, weibull_c)
#endif

          call SatFunc_PressToSat(    satParam, P, se, dse_dP)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)

          dkr_dse   = &
               0.5d0 * se**(-0.5d0) *( 1.d0 - (1.d0 - se **(1.d0/m))**m)**2.d0 + &
               se**(0.5d0) * 2.d0   *( 1.d0 - (1.d0 - se **(1.d0/m))**m) * (1.d0 - se**(1.d0/m)) * se**(1.d0/m - 1.d0)
          dse_dp0   = 0.d0
          
          dkr_dx    = dkr_dP * dP_dx !+ dkr_dse * dse_dp0 * dp0_dx
          dkr_dx    = dkr_dP * dP_dx !+ dkr_dse * dse_dp0 * dp0_dx
          drho_dx   = drho_dP * dP_dx
          d2rho_dx2 = d2rho_dP2 * dP_dx + drho_dP * d2P_dx2

          data_1D(ii) = &
               -((k*kr/mu)*drho_dx + (rho*kr/mu)*dk_dx + (rho*k/mu)*dkr_dx)*(dP_dx) &
               -(rho*k*kr/mu)*(d2P_dx2)
          data_1D(ii) = data_1D(ii)*dx_root
          data_1D(ii) = data_1D(ii) - soil_root_flux(ii)

       end do

    case (ROOT_LIQ_SAT)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root

          call compute_root_alpha_or_deriv        (xx, val=p0)
          call compute_root_lambda_or_deriv       (xx, val=m )
          call compute_root_pressure_or_deriv     (xx, val=P )

          call SatFunc_Set_FETCH2(satParam, p0, m)
          call compute_root_weibullD_or_deriv(xx, weibull_d)
          call compute_root_weibullC_or_deriv(xx, weibull_c)
          call SatFunc_Set_Weibull_RelPerm(satParam, weibull_d, weibull_c)
          call SatFunc_PressToSat(satParam, P, se, dse_dP)

          data_1D(ii) = se
       end do

    case (ROOT_REL_PERM)
       do ii = 1, num_root
          xx = x_root_min + dx_root/2.d0 + (ii-1)*dx_root

          call compute_root_alpha_or_deriv        (xx, val=p0)
          call compute_root_lambda_or_deriv       (xx, val=m )
          call compute_root_pressure_or_deriv     (xx, val=P )

          call SatFunc_Set_FETCH2(satParam, p0, m)
          call compute_root_weibullD_or_deriv(xx, weibull_d)
          call compute_root_weibullC_or_deriv(xx, weibull_c)
          call SatFunc_Set_Weibull_RelPerm(satParam, weibull_d, weibull_c)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)

          data_1D(ii) = kr
       end do

    case (XYLM_POROSITY)
       data_1D(:) = 0.d0

    case (XYLM_PERMEABILITY)
       do ii = 1, num_xylm
          xx = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
#ifdef USE_VG
          call compute_soil_permeability_or_deriv (xx, data_1D(ii))
#else
          call compute_xylm_permeability_or_deriv (xx, data_1D(ii))
#endif
       end do

    case (XYLM_SATFUNC_ALPHA)
       data_1D(:) = xylm_phi88
       do ii = 1, num_xylm
          xx = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
#ifdef USE_VG
          call compute_soil_alpha_or_deriv (xx, data_1D(ii))
#else
          call compute_xylm_alpha_or_deriv (xx, data_1D(ii))
#endif
       end do

    case (XYLM_SATFUNC_LAMBDA)
       do ii = 1, num_xylm
          xx = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
#ifdef USE_VG
          call compute_soil_lambda_or_deriv (xx, data_1D(ii))
#else
          call compute_xylm_lambda_or_deriv (xx, data_1D(ii))
#endif
       end do

    case (XYLM_RES_SAT)
       do ii = 1, num_xylm
          xx = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
#ifdef USE_VG
          call compute_soil_residualsat_or_deriv (xx, data_1D(ii))
#else
          call compute_xylm_residualsat_or_deriv (xx, data_1D(ii))
#endif
       end do

    case (XYLM_WEIBULL_D)
       data_1D(:) = xylm_c1
       do ii = 1, num_xylm
          xx = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
          call compute_xylm_weibullD_or_deriv (xx, data_1D(ii))
       end do

    case (XYLM_WEIBULL_C)
       data_1D(:) = xylm_c2
       do ii = 1, num_xylm
          xx = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
          call compute_xylm_weibullC_or_deriv (xx, data_1D(ii))
       end do

    case (XYLM_SINK_EXPONENT)
       data_1D(1         :  num_xylm) = 0.d0
       data_1D(num_xylm+1:2*num_xylm) = xylm_c3

    case (XYLM_SINK_PRESSURE)
       data_1D(1         :  num_xylm) = 0.d0
       data_1D(num_xylm+1:2*num_xylm) = xylm_phis50

    case (XYLM_CONDUCTANCE_VAL)
       data_1D(:) = 3.d-11

    case (XYLM_INITIAL_PRESSURE)
       P = 0.d0
       do ii = 1, num_xylm
          xx          = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
          call compute_xylm_pressure_or_deriv(xx, P)
          data_1D(ii) = P - 0.1d4
       end do

    case (XYLM_PRESSURE)
       do ii = 1, num_xylm
          xx          = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
          call compute_xylm_pressure_or_deriv(xx, data_1D(ii))
       end do

    case (XYLM_BC_PRESSURE)
       xx = x_xylm_max
       ii = 1
       call compute_xylm_pressure_or_deriv(xx, data_1D(ii))

    case (XYLM_MASS_SOURCE)
       do ii = 1, num_xylm

          xx      = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
          
          call compute_xylm_permeability_or_deriv (xx, val=k , dval_dx=dk_dx                   )
          call compute_xylm_alpha_or_deriv        (xx, val=p0, dval_dx=dp0_dx                  )
          call compute_xylm_lambda_or_deriv       (xx, val=m                                   )
          call compute_xylm_residualsat_or_deriv  (xx, val=sat_res                             )
          call compute_xylm_pressure_or_deriv     (xx, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2 )
          call compute_xylm_residualsat_or_deriv  (xx, sat_res)

#ifdef USE_VG
          call compute_soil_permeability_or_deriv (xx, val=k , dval_dx=dk_dx                   )
          call compute_soil_alpha_or_deriv        (xx, val=p0, dval_dx=dp0_dx                  )
          call compute_soil_lambda_or_deriv       (xx, val=m                                   )
          call compute_soil_residualsat_or_deriv  (xx, val=sat_res                             )
          !call compute_soil_pressure_or_deriv     (xx, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2 )
          call compute_soil_residualsat_or_deriv  (xx, sat_res)
#endif

          call Viscosity(P, 298.15d0, mu, dmu_dP, dmu_dT)
          call Density(P, 298.15d0, DENSITY_TGDPB01,  rho, drho_dP, dden_dT)

          rho     = rho     * FMWH2O
          drho_dP = drho_dP * FMWH2O
          d2rho_dP2 = 0.d0

#ifdef USE_VG
          call SatFunc_Set_VG(satParam, sat_res, p0, m)
#else
          call SatFunc_Set_FETCH2(satParam, p0, m)
          call compute_xylm_weibullD_or_deriv(xx, weibull_d)
          call compute_xylm_weibullC_or_deriv(xx, weibull_c)
          call SatFunc_Set_Weibull_RelPerm(satParam, weibull_d, weibull_c)
#endif

          call SatFunc_PressToSat(satParam, P, se, dse_dP)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)
          
          dkr_dse   = &
               0.5d0 * se**(-0.5d0) *( 1.d0 - (1.d0 - se **(1.d0/m))**m)**2.d0 + &
               se**(0.5d0) * 2.d0   *( 1.d0 - (1.d0 - se **(1.d0/m))**m) * (1.d0 - se**(1.d0/m)) * se**(1.d0/m - 1.d0)
          dse_dp0   = 0.d0
          
          dkr_dx    = dkr_dP * dP_dx !+ dkr_dse * dse_dp0 * dp0_dx
          dkr_dx    = dkr_dP * dP_dx !+ dkr_dse * dse_dp0 * dp0_dx
          drho_dx   = drho_dP * dP_dx
          d2rho_dx2 = d2rho_dP2 * dP_dx + drho_dP * d2P_dx2

          data_1D(ii) = &
               -((k*kr/mu)*drho_dx + (rho*kr/mu)*dk_dx + (rho*k/mu)*dkr_dx)*(dP_dx) &
               -(rho*k*kr/mu)*(d2P_dx2)
          data_1D(ii) = data_1D(ii)*dx_xylm - max_pet*exp(-((P-PRESSURE_REF)/xylm_phis50)**xylm_c3)

       end do

    case (XYLM_PET)
       do ii = 1, num_xylm
          xx = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm
          data_1D(ii) = max_pet
       end do

    case (XYLM_LIQ_SAT)
       do ii = 1, num_xylm
          xx      = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm

          call compute_xylm_alpha_or_deriv        (xx, val=p0)
          call compute_xylm_lambda_or_deriv       (xx, val=m )
          call compute_xylm_pressure_or_deriv     (xx, val=P )

          call SatFunc_Set_FETCH2(satParam, p0, m)
          call compute_xylm_weibullD_or_deriv(xx, weibull_d)
          call compute_xylm_weibullC_or_deriv(xx, weibull_c)
          call SatFunc_Set_Weibull_RelPerm(satParam, weibull_d, weibull_c)
          call SatFunc_PressToSat(satParam, P, se, dse_dP)

          data_1D(ii) = se
       end do

    case (XYLM_REL_PERM)
       do ii = 1, num_xylm
          xx      = x_xylm_min + dx_xylm/2.d0 + (ii-1)*dx_xylm

          call compute_xylm_alpha_or_deriv        (xx, val=p0)
          call compute_xylm_lambda_or_deriv       (xx, val=m )
          call compute_xylm_pressure_or_deriv     (xx, val=P )

          call SatFunc_Set_FETCH2(satParam, p0, m)
          call compute_xylm_weibullD_or_deriv(xx, weibull_d)
          call compute_xylm_weibullC_or_deriv(xx, weibull_c)
          call SatFunc_Set_Weibull_RelPerm(satParam, weibull_d, weibull_c)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)

          data_1D(ii) = kr
       end do

    case default
       write(*,*)data_type
       write(*,*)'Unknown data_type'
       stop
    end select

  end subroutine set_variable_for_problem

  !------------------------------------------------------------------------
  subroutine save_problem_variable(vsfm_mpp, true_soln_filename, &
              soil_variable_id, root_variable_id, xylm_variable_id)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants    , only : AUXVAR_SS
    use MultiPhysicsProbConstants    , only : VAR_BC_SS_CONDITION
    use SystemOfEquationsBaseType    , only : sysofeqns_base_type
    use SystemOfEquationsThermalType , only : sysofeqns_thermal_type
    use petscvec
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    character(len=256)                   :: true_soln_filename
    PetscInt                             :: soil_variable_id, root_variable_id, xylm_variable_id
    !
    PetscInt                             :: soe_auxvar_id
    PetscReal                  , pointer :: val(:)
    Vec                                  :: true_soln
    PetscViewer                          :: viewer
    class(sysofeqns_base_type) , pointer :: base_soe
    character(len=256)                   :: string
    PetscReal, pointer                   :: soil_p(:),root_p(:),xylm_p(:)
    PetscErrorCode                       :: ierr

    base_soe => vsfm_mpp%soe

    allocate(soil_p(num_soil))
    allocate(root_p(num_root))
    allocate(xylm_p(num_xylm))

    call set_variable_for_problem(soil_variable_id, soil_p)
    call set_variable_for_problem(root_variable_id, root_p)
    call set_variable_for_problem(xylm_variable_id, xylm_p)

    string = trim(true_soln_filename)

    call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)

    call VecDuplicate(vsfm_mpp%soe%solver%soln, true_soln, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(true_soln, val, ierr); CHKERRQ(ierr)

    val(1                   :num_soil                   ) = soil_p(:)
    val(1+num_soil          :num_soil+num_root          ) = root_p(:)
    val(1+num_soil+num_root :num_soil+num_root+num_xylm ) = xylm_p(:)

    call VecRestoreArrayF90(true_soln, val, ierr); CHKERRQ(ierr)

    call VecView(true_soln,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

    deallocate(soil_p)
    deallocate(root_p)
    deallocate(xylm_p)

  end subroutine save_problem_variable

end module vsfm_spac_mms_problem
