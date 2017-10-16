module vsfm_spac_campbell_problem

  implicit none
  
#include <petsc/finclude/petsc.h>
  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscInt              :: nz_xylem
  PetscInt              :: nz_root
  PetscInt              :: nz_soil
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost
  PetscReal             :: Campbell_n
  PetscReal             :: Campbell_b
  PetscReal             :: Campbell_he
  PetscReal             :: theta_s
  PetscReal             :: VG_alpha
  PetscReal             :: VG_n
  PetscReal , parameter :: PI  = 4 * atan (1.0_8)

  public :: run_vsfm_spac_campbell_problem
  public :: output_regression_vsfm_spac_campbell_problem
  
contains

  subroutine run_vsfm_spac_campbell_problem()
    !
#include <petsc/finclude/petsc.h>
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use mpp_varpar           , only : mpp_varpar_init
    use petscsys
    use petscvec
    use petscmat
    use petscts
    use petscsnes
    use petscdm
    use petscdmda
    !    
    implicit none
    
    !
    PetscBool          :: converged
    PetscInt           :: converged_reason
    PetscErrorCode     :: ierr
    PetscReal          :: dtime
    PetscReal          :: time
    PetscInt           :: istep, nstep
    PetscBool          :: flg
    PetscBool          :: save_initial_soln, save_final_soln
    character(len=256) :: string
    character(len=256) :: output_suffix
    PetscViewer        :: viewer

    ! Set default settings
    nz_xylem          = 2
    nz_root           = 28
    nz_soil           = 50
    dtime             = 3600.d0
    nstep             = 24
    save_initial_soln = PETSC_FALSE
    save_final_soln   = PETSC_FALSE
    output_suffix     = ''

    ! Get some command line options

    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dtime,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_initial_soln',save_initial_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_final_soln',save_final_soln,flg,ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-output_suffix',output_suffix,flg,ierr)

    ! Initialize the problem
    call Init()  

    time = 0.d0

    do istep = 1, nstep

       call set_bondary_conditions(time)
       time = time + dtime

       ! Run the model
       call vsfm_mpp%sysofeqns%StepDT(dtime, istep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)

    end do

  end subroutine run_vsfm_spac_campbell_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    !
    implicit none
    !

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp()

    ! 2. Add all meshes needed for the MPP
    call add_meshes()

    ! 3. Add all governing equations
    call add_goveqns()

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns()

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars()

    ! 6. Setup the MPP
    call vsfm_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties() 

    ! 9. Set flux type
    call set_conn_flux_type()
    
    ! 8. Set initial conditions
    call set_initial_conditions()

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp()
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    !
#include <petsc/finclude/petsc.h>
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use MultiPhysicsProbVSFM , only : vsfm_mpp
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
  subroutine add_meshes()
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
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
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    !
    implicit none
    !
#include <petsc/finclude/petsc.h>
    !
    PetscReal :: dx, dy, dz
    PetscInt :: imesh, kk
    PetscInt :: nlev, nz, nconn
    PetscInt :: iconn, vert_nconn
    PetscReal, pointer :: soil_xc(:)           ! x-position of grid cell [m]
    PetscReal, pointer :: soil_yc(:)           ! y-position of grid cell [m]
    PetscReal, pointer :: soil_zc(:)           ! z-position of grid cell [m]
    PetscReal, pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
    PetscReal, pointer :: soil_area(:)         ! area of grid cell [m^2]
    PetscReal, pointer :: soil_vol(:)          ! volume of grid cell [m^3]
    PetscInt , pointer :: soil_filter(:)       ! 

    PetscInt, pointer  :: vert_conn_id_up(:)   !
    PetscInt, pointer  :: vert_conn_id_dn(:)   !
    PetscReal, pointer :: vert_conn_dist_up(:) !
    PetscReal, pointer :: vert_conn_dist_dn(:) !
    PetscReal, pointer :: vert_conn_area(:)    !
    PetscInt , pointer :: vert_conn_type(:)    !

    PetscErrorCode :: ierr

    dx = x_column/nx
    dy = y_column/ny
    dz = z_column/nz_soil

    imesh        = 1
    nlev         = nz_xylem + nz_root + nz_soil
    nz           = nlev
    ncells_local = nx*ny*nlev
    ncells_ghost = 0

    allocate(soil_xc     (nz))
    allocate(soil_yc     (nz))
    allocate(soil_zc     (nz))
    allocate(soil_dx     (nz))
    allocate(soil_dy     (nz))
    allocate(soil_dz     (nz))
    allocate(soil_area   (nz))
    allocate(soil_filter (nz))
    allocate(soil_vol    (nz))

    soil_filter (:) = 1
    soil_area   (:) = 1.d0
    soil_dx     (:) = 1.d0
    soil_dy     (:) = 1.d0
    soil_dz     (:) = 1.d0/50.d0
    soil_xc     (:) = 1.d0/2.d0
    soil_yc     (:) = 1.d0/2.d0
    soil_vol    (:) = 1.d0/50.d0

    soil_zc(1) = 0.d0
    soil_zc(2) = 0.d0
    soil_vol(31) = soil_vol(1) / 2.d0
    
    do kk = 3,nz_xylem + nz_root
       soil_zc(kk) = -(dz/2.d0 + dz * (kk - 1))
    enddo

    do kk = nz_xylem + nz_root + 1, nz_xylem + nz_root + nz_soil
       soil_zc(kk) = -(dz/2.d0 + dz * (kk - nz_xylem - nz_root - 1))
    end do

    call mpp_varpar_set_nlevsoi(nz)
    call mpp_varpar_set_nlevgrnd(nz)

    !
    ! Set up the meshes
    !    
    call vsfm_mpp%SetNumMeshes(1)

    call vsfm_mpp%MeshSetName        (imesh, 'Soil mesh')
    call vsfm_mpp%MeshSetOrientation (imesh, MESH_AGAINST_GRAVITY)
    call vsfm_mpp%MeshSetID          (imesh, MESH_CLM_SOIL_COL)
    call vsfm_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

    call vsfm_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME , soil_vol)
    !call vsfm_mpp%MeshComputeVolume          (imesh)

    nconn = 1 + nz_root * 2 + nz_soil - 1

    allocate (vert_conn_id_up   (nconn))
    allocate (vert_conn_id_dn   (nconn))
    allocate (vert_conn_dist_up (nconn))
    allocate (vert_conn_dist_dn (nconn))
    allocate (vert_conn_area    (nconn))
    allocate (vert_conn_type    (nconn))

    iconn = 0
    iconn = iconn + 1
    vert_conn_id_up(iconn)   = 1
    vert_conn_id_dn(iconn)   = 2
    vert_conn_dist_up(iconn) = 0.5d0*dz
    vert_conn_dist_dn(iconn) = 0.5d0*dz
    vert_conn_area(iconn)    = soil_area(1)
    vert_conn_type(iconn)    = CONN_VERTICAL

    do kk = 2, nz_xylem + nz_root -1

       iconn = iconn + 1
       vert_conn_id_up(iconn)   = 2
       vert_conn_id_dn(iconn)   = kk + 1
       vert_conn_dist_up(iconn) = 0.5d0*dz
       vert_conn_dist_dn(iconn) = 0.5d0*dz
       vert_conn_area(iconn)    = soil_area(kk)
       vert_conn_type(iconn)    = CONN_VERTICAL

    end do

    do kk = 2, nz_xylem + nz_root -1

       iconn = iconn + 1
       vert_conn_id_up(iconn)   = kk + 1
       vert_conn_id_dn(iconn)   = kk + 1 + nz_xylem + nz_root
       vert_conn_dist_up(iconn) = 0.5d0*dz
       vert_conn_dist_dn(iconn) = 0.5d0*dz
       vert_conn_area(iconn)    = soil_area(kk)
       vert_conn_type(iconn)    = CONN_VERTICAL

    end do

    do kk = 1, nz_soil - 1
       iconn = iconn + 1
       vert_conn_id_up(iconn)   = kk + nz_xylem + nz_root
       vert_conn_id_dn(iconn)   = kk + nz_xylem + nz_root + 1
       vert_conn_dist_up(iconn) = 0.5d0*dz
       vert_conn_dist_dn(iconn) = 0.5d0*dz
       vert_conn_area(iconn)    = soil_area(kk)
       vert_conn_type(iconn)    = CONN_VERTICAL
    end do
    
    vert_nconn = iconn

    call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
         vert_nconn,  vert_conn_id_up, vert_conn_id_dn,          &
         vert_conn_dist_up, vert_conn_dist_dn,  vert_conn_area,  &
         vert_conn_type)

    deallocate(soil_xc)
    deallocate(soil_yc)
    deallocate(soil_zc)
    deallocate(soil_dx)
    deallocate(soil_dy)
    deallocate(soil_dz)
    deallocate(soil_area)
    deallocate(soil_filter)  

  end subroutine add_meshes

  !------------------------------------------------------------------------

  subroutine add_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : GE_RE
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    !
    ! !ARGUMENTS
    implicit none

    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE', MESH_CLM_SOIL_COL)

    call vsfm_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_DOWNREGULATE_POT_MASS_RATE
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use ConnectionSetType         , only : ConnectionSetDestroy
    use MeshType                  , only : MeshCreateConnectionSet
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn

    ieqn       = 1
    call vsfm_mpp%GovEqnAddCondition(ieqn, COND_SS,   &
         'Potential Mass_Flux', 'kg/s', COND_DOWNREGULATE_POT_MASS_RATE, &
         SOIL_BOTTOM_CELLS)

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    !
    implicit none

    !
    ! Allocate auxvars
    !
    call vsfm_mpp%AllocateAuxVars()
    

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSourceSinkAuxVarRealValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_EXPONENT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use mpp_varcon                , only : denh2o
    use mpp_varcon                , only : grav
    use SaturationFunction        , only : SAT_FUNC_VAN_GENUCHTEN, SAT_FUNC_BROOKS_COREY
    !
    implicit none
    !
    PetscReal , pointer   :: vsfm_watsat(:,:)
    PetscReal , pointer   :: vsfm_hksat(:,:)
    PetscReal , pointer   :: vsfm_bsw(:,:)
    PetscReal , pointer   :: vsfm_sucsat(:,:)
    PetscReal , pointer   :: vsfm_eff_porosity(:,:)
    PetscReal , pointer   :: vsfm_residual_sat(:,:)
    PetscReal , pointer   :: por(:)
    PetscReal , pointer   :: alpha(:)
    PetscReal , pointer   :: lambda(:)
    PetscReal , pointer   :: sat_res(:)
    PetscReal , pointer   :: perm(:)
    PetscReal , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscReal             :: Ks
    PetscInt              :: begc , endc
    integer   , pointer   :: vsfm_filter(:)
    PetscInt  , pointer   :: satfunc_type(:)
    PetscReal , pointer   :: ss_auxvar_value(:)
    PetscInt              :: nz
    !-----------------------------------------------------------------------

    Ks          = 0.001d0 ! [kg s m^{-3}]
    theta_s     = 0.46d0  !
    Campbell_b  = 4.58d0  ! [-]
    Campbell_he = -4.2d0  ! [J kg^{-1}]

    VG_n        = 1.35d0  ! [-]
    VG_alpha    = 0.15d0  ! [kg J^{-1}]

    Campbell_n  = 2.d0 + 3.d0/Campbell_b

    begc = 1
    endc = 1

    nz = nz_xylem + nz_root + nz_soil

    allocate (por          (nz))
    allocate (alpha        (nz))
    allocate (lambda       (nz))
    allocate (sat_res      (nz))
    allocate (satfunc_type (nz))
    allocate (perm         (nz))

    por          (1                      : nz_xylem + nz_root           ) = 0.d0
    sat_res      (1                      : nz_xylem + nz_root           ) = 0.d0
    lambda       (1                      : nz_xylem + nz_root           ) = 1.d0/Campbell_b
    alpha        (1                      : nz_xylem + nz_root           ) = 1.d-3/(-Campbell_he)
    perm         (1                      : nz_xylem + nz_root           ) = Ks/1.d6*8.904156d-4
    satfunc_type (1                      : nz_xylem + nz_root           ) = SAT_FUNC_BROOKS_COREY

    por          (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = theta_s
    sat_res      (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = 0.01d0
    lambda       (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = 1.d0 - 1.d0/VG_n
    alpha        (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = VG_alpha * 1.d-3
    perm         (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = Ks/1.d6*8.904156d-4
    satfunc_type (nz_xylem + nz_root + 1 : nz_xylem + nz_root + nz_soil ) = SAT_FUNC_VAN_GENUCHTEN

        
    call VSFMMPPSetSoilPorosity(vsfm_mpp, 1, por)

    call VSFMMPPSetSaturationFunction(vsfm_mpp, 1, satfunc_type, &
         alpha, lambda, sat_res)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, 1, perm, perm, perm)

    deallocate(por     )
    deallocate(sat_res )
    deallocate(lambda  )
    deallocate(alpha   )

    allocate(ss_auxvar_value(1))
    
    ss_auxvar_value(:) = 10.d0
    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_EXPONENT, ss_auxvar_value)

    ss_auxvar_value(:) = -1500000.d0
    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_PRESSURE, ss_auxvar_value)

    deallocate(ss_auxvar_value)

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use petscsys
    use petscvec
    use petscmat
    use petscts
    use petscsnes
    use petscdm
    use petscdmda
    !
    implicit none
    !
    PetscReal :: theta
    PetscReal :: Se
    PetscInt  :: ii
    PetscReal, pointer :: press_ic(:), v_x(:)
    Vec :: X
    PetscViewer :: viewer
    PetscErrorCode :: ierr

    allocate(press_ic(nz_xylem + nz_root + nz_soil))

    theta = 0.20d0

    Se = theta/theta_s
    do ii = 1, nz_xylem + nz_root + nz_soil
       press_ic(ii) = (Campbell_he * Se**(-Campbell_b))* 1.d3 + 101325.d0
    enddo

    call vsfm_mpp%Restart(press_ic)

    deallocate(press_ic)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_bondary_conditions(time)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC, VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    !
    implicit none
    !
    PetscReal          :: time
    !
    PetscReal, pointer :: ss_value   (:)
    PetscInt           :: soe_auxvar_id
    PetscReal          :: TimeOfDay
    PetscReal          :: fi
    PetscReal          :: ETp
    PetscReal          :: tp

    TimeOfDay = mod(time,(3600.d0*24.d0))/3600.d0
    fi        = 0.9d0
    ETp       = 5.55555555556d-05
    
    tp = fi * ETp * 2.3d0 * (0.05d0 + sin(0.0175d0 * 7.5d0 * TimeOfDay))** 4.d0
    
    allocate(ss_value(1))
    ss_value(:) = -tp 

    soe_auxvar_id = 1
    call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

    deallocate(ss_value   )

  end subroutine set_bondary_conditions

  !------------------------------------------------------------------------
  subroutine set_conn_flux_type()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnRealValue 
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnIntValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunctionAuxVarConn
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use MultiPhysicsProbConstants , only : VAR_FLUX_TYPE
    use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE
    use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants , only : VAR_CAMPBELL_HE
    use MultiPhysicsProbConstants , only : VAR_CAMPBELL_N
    use SaturationFunction        , only : RELPERM_FUNC_CAMPBELL
    use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
    use petscsys
    !
    implicit none
    !
    PetscReal , pointer   :: cond_conn_in(:)
    PetscInt  , pointer   :: flux_type_conn_in(:)
    PetscReal , pointer   :: campbell_he_conn_in(:)
    PetscReal , pointer   :: campbell_n_conn_in(:)
    PetscInt  , pointer   :: satfunc_itype_conn_in(:)
    PetscInt              :: nconn_in
    PetscInt              :: kk
    PetscReal             :: rootDepth
    PetscReal             :: rootMin
    PetscReal             :: rw
    PetscReal             :: r1
    PetscReal             :: RL
    PetscReal             :: L
    PetscInt              :: nz_loc
    PetscReal             :: dz_loc
    PetscReal             :: Ks
    PetscReal             :: theta,psi_soil,Se,K
    PetscReal , pointer   :: z_int(:)
    PetscReal , pointer   :: Rr(:)
    PetscReal , pointer   :: bz(:)


    nconn_in = nz_xylem - 1 + nz_root * 2 + nz_soil - 1

    allocate (cond_conn_in          (nconn_in))
    allocate (flux_type_conn_in     (nconn_in))
    allocate (campbell_he_conn_in   (nconn_in))
    allocate (campbell_n_conn_in    (nconn_in))
    allocate (satfunc_itype_conn_in (nconn_in))

    satfunc_itype_conn_in(:) = 0

    flux_type_conn_in(1                             :nz_xylem - 1 + nz_root * 2) = CONDUCTANCE_FLUX_TYPE
    flux_type_conn_in(nz_xylem - 1 + nz_root * 2 + 1:nconn_in                  ) = DARCY_FLUX_TYPE
    
    call VSFMMPPSetAuxVarConnIntValue( vsfm_mpp, 1, AUXVAR_CONN_INTERNAL, VAR_FLUX_TYPE, flux_type_conn_in)
    
    nz_loc = 50
    dz_loc = 1.d0/nz_loc
    allocate(z_int(nz_loc+1))
    allocate(Rr(nz_loc))
    allocate(bz(nz_loc))

    Ks          = 0.001d0 !
    rootDepth   = 0.6d0
    rootMin     = 0.02d0
    rw          = 25000000000.d0
    r1          = 0.001d0

    !pc          = -1500.d0
    RL          = 1/(3.d6 * 1.d6)

    do kk = 1, nz_loc + 1
       z_int(kk) = (kk-1)*dz_loc
    end do

    cond_conn_in(1) = RL

    do kk = 1, nz_loc
       if (z_int(kk) > rootMin .and. z_int(kk) < rootDepth) then
          L = 40000.d0 * (rootDepth - z_int(kk)) / rootDepth
          Rr(kk) = 2.d0 * rw / (L * (z_int(kk+1) - z_int(kk-1)))
          bz(kk) = ((1.d0 - Campbell_n) * log(PI * r1 * r1 * L) / (2 * PI * L * (z_int(kk+1) - z_int(kk-1))))
       else
          Rr(kk) = 0.d0
          bz(kk) = 0.d0
       endif


       if (kk >= 3 .and. kk <= 30) then
          theta    = 0.1d0

          Se       = theta/theta_s
          psi_soil = Campbell_he * Se**(-Campbell_b)
          K        = Ks !* (Campbell_he / psi_soil)**Campbell_n;

          cond_conn_in          (kk-1     ) = 1.d-6/Rr(kk)
          cond_conn_in          (kk-2 + 29) = 1.d-6/(bz(kk)/K)
          campbell_he_conn_in   (kk-2 + 29) = -campbell_he * 1.d3
          campbell_n_conn_in    (kk-2 + 29) = campbell_n
          satfunc_itype_conn_in (kk-2 + 29) = RELPERM_FUNC_CAMPBELL

       endif
    enddo

    call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, 1, AUXVAR_CONN_INTERNAL, VAR_CONDUCTANCE, cond_conn_in)

    call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, 1, AUXVAR_CONN_INTERNAL, PETSC_FALSE, &
         satfunc_itype_conn_in, campbell_he_conn_in, campbell_n_conn_in)
        
  end subroutine set_conn_flux_type


  !------------------------------------------------------------------------
  subroutine output_regression_vsfm_spac_campbell_problem(filename_base, num_cells)
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use regression_mod            , only : regression_type
    !
    implicit none
    !
    character(len=256)    :: filename_base
    character(len=512)    :: filename
    PetscInt, intent(in)  :: num_cells
    !
    PetscInt              :: output
    character(len=64)     :: category
    character(len=64)     :: name
    PetscReal, pointer    :: data(:)
    type(regression_type) :: regression

    allocate(data(ncells_local))

    call regression%Init(filename_base, num_cells)
    call regression%OpenOutput()

    name = 'liquid_pressure'
    category = 'pressure'
    call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL,  &
         VAR_PRESSURE, -1, data)
    call regression%WriteData(name, category, data)

    name = 'liquid_saturation'
    category = 'general'
    call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL,  &
         VAR_LIQ_SAT, -1, data)
    call regression%WriteData(name, category, data)

    call regression%CloseOutput()
    
    deallocate(data)

  end subroutine output_regression_vsfm_spac_campbell_problem
  
end module vsfm_spac_campbell_problem
