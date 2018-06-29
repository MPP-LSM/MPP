module vsfm_spac_fetch2_problem

#include <petsc/finclude/petsc.h>

  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda

  implicit none

  PetscReal , parameter :: vis      = 8.904156d-4   ! Pa s
  PetscReal , parameter :: rho      = 1000.d0       ! kg m^{-3}
  PetscReal , parameter :: grav     = 9.81          ! m s^{-2}

  PetscInt  , parameter :: nx       = 1             ! -
  PetscInt  , parameter :: ny       = 1             ! -
  PetscReal , parameter :: dx       = 1.d0          ! m
  PetscReal , parameter :: dy       = 1.d0          ! m
  PetscReal , parameter :: dz       = 0.2d0         ! m

  PetscReal , parameter :: porosity = 0.55d0        ! [-]

  PetscInt  , parameter :: oak_nz       = 60        ! -
  PetscReal , parameter :: oak_Asapwood = 0.0088d0  ! m^2
  PetscReal , parameter :: oak_Acrown   = 28.8d0    ! m^2
  PetscReal , parameter :: oak_phis50   = -0.91d6   ! Pa
  PetscReal , parameter :: oak_phi50    = -2.5d6    ! Pa
  PetscReal , parameter :: oak_phi88    = -0.5d6    ! Pa
  PetscReal , parameter :: oak_c1       = 1.7d6     ! Pa
  PetscReal , parameter :: oak_c2       = 3.0d0     ! -
  PetscReal , parameter :: oak_c3       = 12.3d0    ! -
  PetscReal , parameter :: oak_kmax     = 1.6d-6    ! s

  PetscInt  , parameter :: pine_nz       = 85       ! -
  PetscReal , parameter :: pine_Asapwood = 0.0616d0 ! m^2
  PetscReal , parameter :: pine_Acrown   = 46.1d0   ! m^2
  !PetscReal , parameter :: pine_phis50   = -0.18d6  ! Pa
  PetscReal , parameter :: pine_phis50   = -0.91d6  ! Pa
  PetscReal , parameter :: pine_phi50    = -2.2d6   ! Pa
  PetscReal , parameter :: pine_phi88    = -0.5d6   ! Pa
  PetscReal , parameter :: pine_c1       = 1.2d6    ! Pa
  PetscReal , parameter :: pine_c2       = 2.8d0    ! -
  PetscReal , parameter :: pine_c3       = 10.3d0   ! -
  PetscReal , parameter :: pine_kmax     = 1.2d-6   ! s

  ! Parameters for root length density = length-of-root/volume-of-soil  [m_root/m^3_soil]
  PetscInt , parameter :: oak_root_nz        = 15
  PetscReal, parameter :: oak_root_qz        = 3.d0         ! [-]
  PetscReal, parameter :: oaK_root_d         = 3.d0         ! [m]
  PetscReal, parameter :: oak_root_rld0      = 4.d4         ! [m/m^3]
  PetscReal, parameter :: oak_root_radius    = 2.d-3        ! [m]

  PetscReal , parameter :: oak_rad_root_kmax = 1.6d-6   !
  PetscReal , parameter :: oak_rad_root_phi50= -2.5d6   !
  PetscReal , parameter :: oak_rad_root_phi88= -0.5d6   !
  PetscReal , parameter :: oak_rad_root_c1   = 1.7d6    !
  PetscReal , parameter :: oak_rad_root_c2   = 3.0d0    ! -

  PetscReal , parameter :: oak_axi_root_kmax = 1.6d-6   !
  PetscReal , parameter :: oak_axi_root_phi50= -2.5d6   !
  PetscReal , parameter :: oak_axi_root_phi88= -0.5d6   !
  PetscReal , parameter :: oak_axi_root_c1   = 1.7d6    !
  PetscReal , parameter :: oak_axi_root_c2   = 3.0d0    ! -

  PetscInt , parameter :: pine_root_nz        = 25
  PetscReal, parameter :: pine_root_qz        = 3.d0         ! [-]
  PetscReal, parameter :: pine_root_d         = 5.d0         ! [m]
  PetscReal, parameter :: pine_root_rld0      = 4.d4         ! [m/m^3]
  PetscReal, parameter :: pine_root_radius    = 2.d-3        ! [m]

  PetscReal , parameter :: pine_rad_root_kmax = 1.6d-6   !
  PetscReal , parameter :: pine_rad_root_phi50= -2.5d6   !
  PetscReal , parameter :: pine_rad_root_phi88= -0.5d6   !
  PetscReal , parameter :: pine_rad_root_c1   = 1.7d6    !
  PetscReal , parameter :: pine_rad_root_c2   = 3.0d0    ! -

  PetscReal , parameter :: pine_axi_root_kmax = 1.6d-6   !
  PetscReal , parameter :: pine_axi_root_phi50= -2.5d6   !
  PetscReal , parameter :: pine_axi_root_phi88= -0.5d6   !
  PetscReal , parameter :: pine_axi_root_c1   = 1.7d6    !
  PetscReal , parameter :: pine_axi_root_c2   = 3.0d0    ! -

  PetscInt , parameter :: soil_nz         = 25
  PetscReal, parameter :: soil_perm       = 6.83d-11      ! [m^2]
  PetscReal, parameter :: soil_sat_res    = 0.06d0        ! [-]
  PetscReal, parameter :: soil_alpha      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: soil_vg_m       = 0.33d0        ! [-]
  PetscReal, parameter :: soil_por        = 0.5d0         ! [-]

  PetscInt :: neqns

  character(len=256) :: problem_type

  character(len=256) :: ic_filename
  PetscBool          :: ic_file_specified

  PetscInt              :: ncells_local

  PetscReal             :: theta_s
  PetscReal             :: VG_alpha
  PetscReal             :: VG_n
  PetscReal , parameter :: PI  = 4 * atan (1.0_8)

  public :: run_vsfm_spac_fetch2_problem
  
contains

  subroutine run_vsfm_spac_fetch2_problem()
    !
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use mpp_varpar           , only : mpp_varpar_init
    !
    implicit none
    !
    PetscBool          :: converged
    PetscInt           :: converged_reason
    PetscInt           :: nET
    PetscErrorCode     :: ierr
    PetscReal          :: dtime
    PetscReal          :: time
    PetscInt           :: istep, nstep
    PetscBool          :: flg
    PetscBool          :: save_initial_soln, save_final_soln
    PetscBool          :: print_actual_et
    character(len=256) :: string
    character(len=256) :: pet_file
    character(len=256) :: soil_bc_file
    character(len=256) :: actual_et_file
    Vec                :: ET
    Vec                :: SoilBC
    Vec                :: Actual_ET
    PetscInt           :: vec_size
    PetscInt           :: ntimes
    PetscReal, pointer :: act_et_p(:)
    PetscViewer        :: viewer
    PetscBool          :: error_in_cmd_options

    ! Set default settings
    dtime                  = 180.d0
    nstep                  = 24
    save_initial_soln      = PETSC_FALSE
    save_final_soln        = PETSC_FALSE
    print_actual_et        = PETSC_FALSE
    error_in_cmd_options   = PETSC_FALSE
    problem_type           = 'oak'
    ic_file_specified      = PETSC_FALSE

    ! Get some command line options

    call PetscOptionsGetInt  (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-nstep', nstep,flg, ierr)

    call PetscOptionsGetBool (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-save_initial_soln', save_initial_soln, flg, ierr)

    call PetscOptionsGetBool (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-save_final_soln', save_final_soln, flg,ierr)

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-pet_file', pet_file, flg, ierr)
    if (.not.flg) then
       write(*,*)'ERROR: Specify the potential ET filename via -pet_file <filename>'
       error_in_cmd_options = PETSC_TRUE
    end if

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-soil_bc_file',soil_bc_file,flg,ierr)
    if (.not.flg) then
       write(*,*)'ERROR: Specify the soil boundary condition filename via -soil_bc_file <filename>'
       error_in_cmd_options = PETSC_TRUE
    end if
    
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-actual_et_file', actual_et_file, flg, ierr)
    if (flg) then
       print_actual_et = PETSC_TRUE
    end if

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-ic_file', ic_filename, flg, ierr)
    if (flg) ic_file_specified = PETSC_TRUE

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-problem_type', problem_type, flg, ierr)
    if (flg) then
       select case(trim(problem_type))
       case ('oak')
       case ('pine')
       case ('oak_and_pine')
       case ('oak_spac')
       case ('pine_spac')
       case ('oak_pine_spac')
       case default
          write(*,*)'ERROR: Unknown option specified via -problem_type <value>.'
          write(*,*)'       Valid value for -problem_type are:'
          write(*,*)'         oak'
          write(*,*)'         pine'
          write(*,*)'         oak_and_pine'
          write(*,*)'         oak_spac'
          write(*,*)'         pine_spac'
          write(*,*)'         oak_pine_spac'
          error_in_cmd_options = PETSC_TRUE
       end select
    end if

    if (error_in_cmd_options) stop

    select case(trim(problem_type))
    case ('oak')
       nET       = oak_nz
    case ('pine')
       nET       = pine_nz
    case ('oak_and_pine')
       nET      = oak_nz + pine_nz
    case ('oak_spac')
       nET       = oak_nz
    case ('pine_spac')
       nET       = pine_nz
    case ('oak_pine_spac')
       nET      = oak_nz + pine_nz
    case default
       write(*,*)'Unable to determine nET for problem_type = ',trim(problem_type)
       stop
    end select

    ! Read the PET file
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, pet_file, FILE_MODE_READ, viewer, ierr); CHKERRQ(ierr)
    call VecCreate(PETSC_COMM_WORLD, ET, ierr); CHKERRQ(ierr)
    call VecLoad(ET, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
    call VecGetSize(ET, vec_size, ierr); CHKERRQ(ierr)

    ntimes = vec_size/nET
    if (nstep > ntimes) then
       write(*,*)'ERROR: Value for nstep exceeds the number of timesteps ' // &
            'for which potential ET data is available.'
       write(*,*)'nstep = ',nstep
       write(*,*)'ntimes= ',ntimes
       stop
    end if

    ! Read the soil bc file
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, soil_bc_file, FILE_MODE_READ, viewer, ierr); CHKERRQ(ierr)
    call VecCreate(PETSC_COMM_WORLD, SoilBC, ierr); CHKERRQ(ierr)
    call VecLoad(SoilBC, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
    call VecGetSize(SoilBC, vec_size, ierr); CHKERRQ(ierr)

    ntimes = vec_size
    if (nstep > ntimes) then
       write(*,*)'ERROR: Value for nstep exceeds the number of timesteps ' // &
            'for which soil boundary condition data is available.'
       write(*,*)'nstep = ',nstep
       write(*,*)'ntimes= ',ntimes
       stop
    end if

    ! Initialize the problem
    call Init()

    time = 0.d0
    call VecCreateSeq(PETSC_COMM_SELF, nstep*nET, Actual_ET, ierr)
    call VecGetArrayF90(Actual_ET, act_et_p, ierr)

    call PetscViewerBinaryOpen(PETSC_COMM_SELF,'init.bin',FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
    call VecView(vsfm_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

    do istep = 1, nstep

       select case(trim(problem_type))
       case ('oak')
          call set_boundary_conditions_for_single_tree(nET, istep, ET, SoilBC)
       case ('pine')
          call set_boundary_conditions_for_single_tree(nET, istep, ET, SoilBC)
       case ('oak_and_pine')
          call set_boundary_conditions_for_two_trees(nET, istep, ET, SoilBC)
       case default
          write(*,*)'Unable to set BC for problem_type = ' // trim(problem_type)
          stop
       end select

       time = time + dtime

       ! Run the model
       call vsfm_mpp%soe%StepDT(dtime, istep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)

       if (.not.converged) then
          write(*,*)'Model failed to converge'
          stop
       end if

       if (print_actual_et) then
          select case(trim(problem_type))
          case ('oak')
             call diagnose_actual_sink_for_single_tree(nET, istep, act_et_p)
          case ('pine')
             call diagnose_actual_sink_for_single_tree(nET, istep, act_et_p)
          case ('oak_and_pine')
             call diagnose_actual_sink_for_single_tree(nET, istep, act_et_p)
          case default
             write(*,*)'Unable to diagnose actual ET for problem_type = ' // trim(problem_type)
             stop
          end select
       end if

    end do
    call VecRestoreArrayF90(Actual_ET, act_et_p, ierr)

    call PetscViewerBinaryOpen(PETSC_COMM_SELF,'final.bin',FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
    call VecView(vsfm_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

    if (print_actual_et) then
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,actual_et_file,FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(Actual_ET,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    end if

    call VecDestroy(Actual_ET, ierr); CHKERRQ(ierr)

  end subroutine run_vsfm_spac_fetch2_problem
  
  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use SystemOfEquationsVSFMType , only : VSFMSOEUpdateConnections
    !
    implicit none

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp()

    ! 2. Add all meshes needed for the MPP
    call add_meshes()

    ! 3. Add all governing equations
    call add_goveqn(neqns)

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns()

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars()

    ! 6. Setup the MPP
    call vsfm_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties()

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
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
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
    call vsfm_mpp%SetName    ('Tree level hydrodynamics')
    call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes()
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp

    implicit none

    select case(trim(problem_type))
    case ('oak')
       neqns = 1
       call vsfm_mpp%SetNumMeshes(neqns)
       call add_xylem_mesh_for_single_tree(1, oak_nz, oak_Asapwood)

    case ('pine')
       neqns = 1
       call vsfm_mpp%SetNumMeshes(neqns)
       call add_xylem_mesh_for_single_tree(1, pine_nz, pine_Asapwood)

    case ('oak_and_pine')
       neqns = 1
       call vsfm_mpp%SetNumMeshes(neqns)
       call add_xylem_mesh_for_two_trees()

    case ('oak_spac')
       neqns = 3
       call vsfm_mpp%SetNumMeshes(neqns)

       call add_xylem_mesh_for_single_tree(1, oak_nz, oak_Asapwood)
       call add_root_mesh_for_single_tree( 2, 'oak')
       call add_soil_mesh(3)

    case ('pine_spac')
       neqns = 3
       call vsfm_mpp%SetNumMeshes(neqns)

       call add_xylem_mesh_for_single_tree(1, pine_nz, pine_Asapwood)
       call add_root_mesh_for_single_tree( 2, 'pine')
       call add_soil_mesh(3)

    case ('oak_pine_spac')
       neqns = 5
       call vsfm_mpp%SetNumMeshes(neqns)

       call add_xylem_mesh_for_single_tree(1, oak_nz, oak_Asapwood)
       call add_root_mesh_for_single_tree( 2, 'oak')
       call add_xylem_mesh_for_single_tree(3, pine_nz, pine_Asapwood)
       call add_root_mesh_for_single_tree( 4, 'pine')
       call add_soil_mesh(5)

    case default
       write(*,*)'Error while adding mesh. Unsupported problem_type : ',trim(problem_type)
       stop

    end select

  end subroutine add_meshes

  !------------------------------------------------------------------------
  subroutine add_xylem_mesh_for_single_tree(imesh, local_nz, local_Asapwood)
    !
    use MeshType                  , only : mesh_type
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    !
    implicit none
    !
    PetscInt, intent(in):: imesh
    PetscInt            :: local_nz
    PetscReal           :: local_Asapwood

    !
    class(mesh_type) , pointer :: mesh
    PetscInt            :: kk
    PetscInt            :: nlev
    PetscInt            :: iconn, vert_nconn

    PetscInt            :: ncells_xylem

    PetscInt            :: nconn_xylem

    PetscReal , pointer :: xc_xylem(:)                ! x-position of grid cell [m]
    PetscReal , pointer :: yc_xylem(:)                ! y-position of grid cell [m]
    PetscReal , pointer :: zc_xylem(:)                ! z-position of grid cell [m]
    PetscReal , pointer :: dx_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: dy_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: dz_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: area_xylem(:)              ! area of grid cell [m^2]
    PetscReal , pointer :: vol_xylem(:)               ! volume of grid cell [m^3]
    PetscInt  , pointer :: filter_xylem(:)            ! 
    
    PetscInt  , pointer :: vert_conn_id_up_xylem(:)   !
    PetscInt  , pointer :: vert_conn_id_dn_xylem(:)   !
    PetscReal , pointer :: vert_conn_dist_up_xylem(:) !
    PetscReal , pointer :: vert_conn_dist_dn_xylem(:) !
    PetscReal , pointer :: vert_conn_area_xylem(:)    !
    PetscInt  , pointer :: vert_conn_type_xylem(:)    !

    PetscErrorCode      :: ierr

    ncells_xylem =local_nz

    allocate(xc_xylem     (ncells_xylem ))
    allocate(yc_xylem     (ncells_xylem ))
    allocate(zc_xylem     (ncells_xylem ))
    allocate(dx_xylem     (ncells_xylem ))
    allocate(dy_xylem     (ncells_xylem ))
    allocate(dz_xylem     (ncells_xylem ))
    allocate(area_xylem   (ncells_xylem ))
    allocate(filter_xylem (ncells_xylem ))
    allocate(vol_xylem    (ncells_xylem ))
    
    call set_xylem_geometric_attributes( &
         1, local_nz, local_Asapwood,   &
         filter_xylem, area_xylem, vol_xylem, &
         dx_xylem, dy_xylem, dz_xylem, &
         xc_xylem, yc_xylem, zc_xylem)

    call mpp_varpar_set_nlevsoi(local_nz)
    call mpp_varpar_set_nlevgrnd(local_nz)

    ! Xylem Mesh
    allocate(mesh)

    call SetMeshPreliminaries(mesh, 'Xylem mesh', ncells_xylem, local_nz)

    call SetMeshGeometricAttributes( &
         mesh, filter_xylem, &
         xc_xylem, yc_xylem, zc_xylem, &
         dx_xylem, dy_xylem, dz_xylem, &
         area_xylem, vol_xylem)
         
    nconn_xylem = local_nz - 1

    allocate (vert_conn_id_up_xylem   (nconn_xylem))
    allocate (vert_conn_id_dn_xylem   (nconn_xylem))
    allocate (vert_conn_dist_up_xylem (nconn_xylem))
    allocate (vert_conn_dist_dn_xylem (nconn_xylem))
    allocate (vert_conn_area_xylem    (nconn_xylem))
    allocate (vert_conn_type_xylem    (nconn_xylem))

    iconn = 0
    do kk = 1, local_nz - 1
       iconn = iconn + 1
       vert_conn_id_up_xylem(iconn)   = kk 
       vert_conn_id_dn_xylem(iconn)   = kk + 1
       vert_conn_dist_up_xylem(iconn) = 0.5d0*dz
       vert_conn_dist_dn_xylem(iconn) = 0.5d0*dz
       vert_conn_area_xylem(iconn)    = area_xylem(kk)
       vert_conn_type_xylem(iconn)    = CONN_VERTICAL
    end do

    call mesh%CreateAndAddConnectionSet( &
         CONN_SET_INTERNAL,              &
         nconn_xylem,  vert_conn_id_up_xylem, vert_conn_id_dn_xylem,          &
         vert_conn_dist_up_xylem, vert_conn_dist_dn_xylem,  vert_conn_area_xylem,  &
         vert_conn_type_xylem)

    call vsfm_mpp%AddMesh(imesh, mesh)

    call mesh%Clean()

    deallocate (xc_xylem                )
    deallocate (yc_xylem                )
    deallocate (zc_xylem                )
    deallocate (dx_xylem                )
    deallocate (dy_xylem                )
    deallocate (dz_xylem                )
    deallocate (area_xylem              )
    deallocate (filter_xylem            )
    deallocate (vol_xylem               )

    deallocate (vert_conn_id_up_xylem   )
    deallocate (vert_conn_id_dn_xylem   )
    deallocate (vert_conn_dist_up_xylem )
    deallocate (vert_conn_dist_dn_xylem )
    deallocate (vert_conn_area_xylem    )
    deallocate (vert_conn_type_xylem    )

  end subroutine add_xylem_mesh_for_single_tree

  !------------------------------------------------------------------------

  subroutine add_xylem_mesh_for_two_trees()
    !
    use MeshType                  , only : mesh_type
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    !
    implicit none
    !

    !
    PetscInt            :: imesh, kk
    PetscInt            :: nlev
    PetscInt            :: iconn, vert_nconn

    PetscInt            :: ncells_xylem

    PetscInt            :: nconn_xylem

    class(mesh_type) , pointer :: mesh
    PetscReal , pointer :: xc_xylem(:)                ! x-position of grid cell [m]
    PetscReal , pointer :: yc_xylem(:)                ! y-position of grid cell [m]
    PetscReal , pointer :: zc_xylem(:)                ! z-position of grid cell [m]
    PetscReal , pointer :: dx_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: dy_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: dz_xylem(:)                ! layer thickness of grid cell [m]
    PetscReal , pointer :: area_xylem(:)              ! area of grid cell [m^2]
    PetscReal , pointer :: vol_xylem(:)               ! volume of grid cell [m^3]
    PetscInt  , pointer :: filter_xylem(:)            ! 
    
    PetscInt  , pointer :: vert_conn_id_up_xylem(:)   !
    PetscInt  , pointer :: vert_conn_id_dn_xylem(:)   !
    PetscReal , pointer :: vert_conn_dist_up_xylem(:) !
    PetscReal , pointer :: vert_conn_dist_dn_xylem(:) !
    PetscReal , pointer :: vert_conn_area_xylem(:)    !
    PetscInt  , pointer :: vert_conn_type_xylem(:)    !

    PetscInt :: ibeg, iend
    PetscReal :: Asapwood
    PetscErrorCode      :: ierr

    ncells_xylem = oak_nz + pine_nz

    allocate(xc_xylem     (ncells_xylem ))
    allocate(yc_xylem     (ncells_xylem ))
    allocate(zc_xylem     (ncells_xylem ))
    allocate(dx_xylem     (ncells_xylem ))
    allocate(dy_xylem     (ncells_xylem ))
    allocate(dz_xylem     (ncells_xylem ))
    allocate(area_xylem   (ncells_xylem ))
    allocate(filter_xylem (ncells_xylem ))
    allocate(vol_xylem    (ncells_xylem ))

    ibeg = 1; iend = oak_nz; Asapwood = oak_Asapwood
    call set_xylem_geometric_attributes( &
         ibeg, iend, Asapwood, &
         filter_xylem, area_xylem, vol_xylem, &
         dx_xylem, dy_xylem, dz_xylem, &
         xc_xylem, yc_xylem, zc_xylem)

    ibeg = 1 + oak_nz; iend = oak_nz + pine_nz; Asapwood = pine_Asapwood
    call set_xylem_geometric_attributes( &
         ibeg, iend, Asapwood, &
         filter_xylem, area_xylem, vol_xylem, &
         dx_xylem, dy_xylem, dz_xylem, &
         xc_xylem, yc_xylem, zc_xylem)

    call mpp_varpar_set_nlevsoi(oak_nz + pine_nz)
    call mpp_varpar_set_nlevgrnd(oak_nz+ pine_nz)

    ! Xylem Mesh
    imesh        = 1

    allocate(mesh)

    call SetMeshPreliminaries(mesh, 'Xylem mesh', ncells_xylem, oak_nz + pine_nz)

    call SetMeshGeometricAttributes( &
         mesh, filter_xylem, &
         xc_xylem, yc_xylem, zc_xylem, &
         dx_xylem, dy_xylem, dz_xylem, &
         area_xylem, vol_xylem)

    nconn_xylem = oak_nz - 1 + pine_nz - 1

    allocate (vert_conn_id_up_xylem   (nconn_xylem))
    allocate (vert_conn_id_dn_xylem   (nconn_xylem))
    allocate (vert_conn_dist_up_xylem (nconn_xylem))
    allocate (vert_conn_dist_dn_xylem (nconn_xylem))
    allocate (vert_conn_area_xylem    (nconn_xylem))
    allocate (vert_conn_type_xylem    (nconn_xylem))

    iconn = 0
    do kk = 1, oak_nz - 1
       iconn = iconn + 1
       vert_conn_id_up_xylem(iconn)   = kk 
       vert_conn_id_dn_xylem(iconn)   = kk + 1
       vert_conn_dist_up_xylem(iconn) = 0.5d0*dz
       vert_conn_dist_dn_xylem(iconn) = 0.5d0*dz
       vert_conn_area_xylem(iconn)    = area_xylem(kk)
       vert_conn_type_xylem(iconn)    = CONN_VERTICAL
    end do

    do kk = oak_nz + 1, oak_nz + pine_nz - 1
       iconn = iconn + 1
       vert_conn_id_up_xylem(iconn)   = kk 
       vert_conn_id_dn_xylem(iconn)   = kk + 1
       vert_conn_dist_up_xylem(iconn) = 0.5d0*dz
       vert_conn_dist_dn_xylem(iconn) = 0.5d0*dz
       vert_conn_area_xylem(iconn)    = area_xylem(kk)
       vert_conn_type_xylem(iconn)    = CONN_VERTICAL
    end do

    call mesh%CreateAndAddConnectionSet( &
         CONN_SET_INTERNAL,              &
         nconn_xylem,  vert_conn_id_up_xylem, vert_conn_id_dn_xylem,          &
         vert_conn_dist_up_xylem, vert_conn_dist_dn_xylem,  vert_conn_area_xylem,  &
         vert_conn_type_xylem)

    call vsfm_mpp%AddMesh(imesh, mesh)

    call mesh%Clean()

    deallocate (xc_xylem                )
    deallocate (yc_xylem                )
    deallocate (zc_xylem                )
    deallocate (dx_xylem                )
    deallocate (dy_xylem                )
    deallocate (dz_xylem                )
    deallocate (area_xylem              )
    deallocate (filter_xylem            )
    deallocate (vol_xylem               )

    deallocate (vert_conn_id_up_xylem   )
    deallocate (vert_conn_id_dn_xylem   )
    deallocate (vert_conn_dist_up_xylem )
    deallocate (vert_conn_dist_dn_xylem )
    deallocate (vert_conn_area_xylem    )
    deallocate (vert_conn_type_xylem    )

  end subroutine add_xylem_mesh_for_two_trees

  !------------------------------------------------------------------------
  subroutine set_xylem_geometric_attributes(ibeg, iend, Asapwood, &
       filter_xylem, area_xylem, vol_xylem, &
       dx_xylem, dy_xylem, dz_xylem, &
       xc_xylem, yc_xylem, zc_xylem)

    implicit none

    PetscInt            :: ibeg, iend
    PetscReal           :: Asapwood
    PetscReal , pointer :: xc_xylem(:)
    PetscReal , pointer :: yc_xylem(:)
    PetscReal , pointer :: zc_xylem(:)
    PetscReal , pointer :: dx_xylem(:)
    PetscReal , pointer :: dy_xylem(:)
    PetscReal , pointer :: dz_xylem(:)
    PetscReal , pointer :: area_xylem(:)
    PetscReal , pointer :: vol_xylem(:)
    PetscInt  , pointer :: filter_xylem(:)

    PetscInt            :: kk, nz

    filter_xylem (ibeg:iend) = 1
    area_xylem   (ibeg:iend) = Asapwood
    dx_xylem     (ibeg:iend) = sqrt(Asapwood)
    dy_xylem     (ibeg:iend) = sqrt(Asapwood)
    dz_xylem     (ibeg:iend) = dz
    xc_xylem     (ibeg:iend) = sqrt(Asapwood)/2.d0
    yc_xylem     (ibeg:iend) = sqrt(Asapwood)/2.d0

    nz = iend - ibeg + 1

    zc_xylem(ibeg)  = -0.17d0 + nz * dz
    vol_xylem(ibeg) = area_xylem(ibeg) * dz
    do kk = ibeg+1, iend
       zc_xylem(kk)  = -(dz/2.d0 + dz * (kk - ibeg)) + nz * dz
       vol_xylem(kk) = area_xylem(kk) * dz
    enddo

  end subroutine set_xylem_geometric_attributes

  !------------------------------------------------------------------------

  subroutine add_root_mesh_for_single_tree( imesh, tree_name)

    use MeshType                  , only : mesh_type
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_SET_CONDITIONS
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL

    implicit none

    PetscInt                   :: imesh
    character(len=*)           :: tree_name

    class(mesh_type) , pointer :: mesh
    PetscInt                   :: kk
    PetscInt                   :: nz
    PetscReal                  :: qz
    PetscReal                  :: dd
    PetscReal                  :: rld0
    PetscReal                  :: radius
    PetscReal                  :: root_len_den, root_len
    PetscReal        , pointer :: xc_1d(:), yc_1d(:), zc_1d(:)
    PetscReal        , pointer :: dx_1d(:), dy_1d(:), dz_1d(:)
    PetscReal        , pointer :: area(:), vol(:)
    PetscInt         , pointer :: filter(:)
    PetscInt                   :: nconn
    PetscInt         , pointer :: conn_id_up(:)
    PetscInt         , pointer :: conn_id_dn(:)
    PetscReal        , pointer :: conn_dist_up(:)
    PetscReal        , pointer :: conn_dist_dn(:)
    PetscReal        , pointer :: conn_area(:)
    PetscInt         , pointer :: conn_type(:)
    PetscReal        , pointer :: conn_unit_vec(:,:)

    select case(trim(tree_name))
    case ('oak')
       nz     = oak_root_nz
       qz     = oak_root_qz
       dd     = oak_root_d
       rld0   = oak_root_rld0
       radius = oak_root_radius

    case ('pine')
       nz     = pine_root_nz
       qz     = pine_root_qz
       dd     = pine_root_d
       rld0   = pine_root_rld0
       radius = pine_root_radius

    case default
       write(*,*)'Unable to set root mesh for tree_name = '//trim(tree_name)
       stop
    end select

    allocate (xc_1d  (nz))
    allocate (yc_1d  (nz))
    allocate (zc_1d  (nz))
    allocate (dx_1d  (nz))
    allocate (dy_1d  (nz))
    allocate (dz_1d  (nz))
    allocate (vol    (nz))
    allocate (area   (nz))
    allocate (filter (nz))

    do kk = 1, nz
       filter(kk)   = 1

       zc_1d(kk)    = -(kk-1)*dz - dz/2.0

       root_len_den = rld0 * (1.d0 - abs(zc_1d(kk))/dd)*exp(-qz*abs(zc_1d(kk))/dd)
       root_len     = root_len_den * & ! [m/m^3]
                      (dx*dy*dz)       ! [m^3]

       vol(kk)      = PI*(radius**2.d0)*root_len
       area(kk)     = 2.d0*PI*radius*root_len

       xc_1d(kk)    = radius
       yc_1d(kk)    = root_len

       dx_1d(kk)    = dx
       dy_1d(kk)    = dy
       dz_1d(kk)    = dz

    end do

    allocate(mesh)
    call SetMeshPreliminaries(mesh, trim(tree_name) // ' root', nz, nz)

    call SetMeshGeometricAttributes( &
         mesh, filter, &
         xc_1d, yc_1d, zc_1d, &
         dx_1d, dy_1d, dz_1d, &
         area, vol)

    ! Add connection set for root-soil BC
    nconn = nz
    allocate (conn_id_up   (nconn   ))
    allocate (conn_id_dn   (nconn   ))
    allocate (conn_dist_up (nconn   ))
    allocate (conn_dist_dn (nconn   ))
    allocate (conn_area    (nconn   ))
    allocate (conn_type    (nconn   ))
    allocate (conn_unit_vec(nconn,3 ))

    conn_unit_vec(:,1) = -1.d0
    conn_unit_vec(:,2) = 0.d0
    conn_unit_vec(:,3) = 0.d0

    do kk = 1, nz
       root_len_den     = rld0 * (1.d0 - abs(zc_1d(kk))/dd)*exp(-qz*abs(zc_1d(kk))/dd)
       conn_id_up(kk)   = 0
       conn_id_dn(kk)   = kk
       conn_dist_up(kk) = 0.d0
       conn_dist_dn(kk) = (PI*root_len_den)**(-0.5d0)
       conn_area(kk)    = area(kk)
    end do

    call mesh%CreateAndAddConnectionSet( &
         CONN_SET_CONDITIONS,            &
         nconn,  conn_id_up, conn_id_dn, &
         conn_dist_up, conn_dist_dn,     &
         conn_area, conn_type, conn_unit_vec)

    ! Add connection set for root-xylem BC
    nconn = 1
    do kk = 1, nconn
       root_len_den     = rld0 * (1.d0 - abs(zc_1d(kk))/dd)*exp(-qz*abs(zc_1d(kk))/dd)
       root_len         = root_len_den * (dx*dy*dz)
       conn_id_up(kk)   = 0
       conn_id_dn(kk)   = kk
       conn_dist_up(kk) = 0.d0
       conn_dist_dn(kk) = dz/2.d0
       conn_area(kk)    = PI*(root_len**2.d0)
    end do

    call mesh%CreateAndAddConnectionSet( &
         CONN_SET_CONDITIONS,            &
         nconn,  conn_id_up, conn_id_dn, &
         conn_dist_up, conn_dist_dn,     &
         conn_area, conn_type, conn_unit_vec)

    ! Add internal connection set
    nconn = nz-1
    do kk = 1, nz-1
       root_len_den = &
            (rld0 * (1.d0 - abs(zc_1d(kk  ))/dd)*exp(-qz*abs(zc_1d(kk  ))/dd) + &
             rld0 * (1.d0 - abs(zc_1d(kk+1))/dd)*exp(-qz*abs(zc_1d(kk+1))/dd))/2.d0

       root_len = root_len_den * (dx*dy*dz)

       conn_id_up(kk)   = kk
       conn_id_dn(kk)   = kk+1

       conn_dist_up(kk) = dz/2.d0
       conn_dist_dn(kk) = dz/2.d0
       conn_area(kk)    = PI*(root_len**2.d0)
    end do

    call mesh%CreateAndAddConnectionSet( &
         CONN_SET_INTERNAL,              &
         nconn,  conn_id_up, conn_id_dn, &
         conn_dist_up, conn_dist_dn,     &
         conn_area, conn_type)

    call vsfm_mpp%AddMesh(imesh, mesh)

    call mesh%Clean()

    deallocate (xc_1d        )
    deallocate (yc_1d        )
    deallocate (zc_1d        )
    deallocate (dx_1d        )
    deallocate (dy_1d        )
    deallocate (dz_1d        )
    deallocate (vol          )
    deallocate (area         )
    deallocate (filter       )
    deallocate (conn_id_up   )
    deallocate (conn_id_dn   )
    deallocate (conn_dist_up )
    deallocate (conn_dist_dn )
    deallocate (conn_area    )
    deallocate (conn_type    )
    deallocate (conn_unit_vec)

  end subroutine add_root_mesh_for_single_tree

  !------------------------------------------------------------------------
  subroutine add_soil_mesh(imesh)

    use MeshType                  , only : mesh_type
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
    use MultiPhysicsProbConstants , only : CONN_SET_CONDITIONS
    use MultiPhysicsProbConstants , only : CONN_SET_LATERAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use MultiPhysicsProbConstants , only : CONN_IN_Z_DIR
    use mpp_mesh_utils

    implicit none

    PetscInt :: imesh

    class(mesh_type) , pointer :: mesh
    PetscInt :: ncells,nz
    PetscReal                  :: x_min, y_min, z_min
    PetscReal        , pointer :: xc_1d(:), yc_1d(:), zc_1d(:)
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

    x_min = 0.d0
    y_min = 0.d0
    z_min = 0.d0
    nz    = soil_nz

    call ComputeXYZCentroids( nx, ny, nz, dx, dy, dz, &
         x_min, y_min, z_min, xc_1d, yc_1d, zc_1d)

    ncells = nx*ny*nz
    
    allocate(filter(ncells))
    allocate(dx_1d (ncells))
    allocate(dy_1d (ncells))
    allocate(dz_1d (ncells))
    allocate(area  (ncells))
    allocate(volume(ncells))

    filter (:) = 1
    dx_1d  (:) = dx
    dy_1d  (:) = dy
    dz_1d  (:) = dz
    area   (:) = dy*dz
    volume (:) = dx*dy*dz

    allocate(mesh)

    call SetMeshPreliminaries(mesh, 'Soil', ncells, nz)

    call SetMeshGeometricAttributes( &
         mesh, filter, &
         xc_1d, yc_1d, zc_1d, &
         dx_1d, dy_1d, dz_1d, &
         area, volume)

    call ComputeInternalConnections(nx, ny, nz, dx, dy, dz, CONN_IN_Z_DIR, &
         nconn, conn_id_up, conn_id_dn, conn_dist_up, conn_dist_dn, &
         conn_area, conn_type)

    call mesh%CreateAndAddConnectionSet( &
         CONN_SET_INTERNAL,              &
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

  end subroutine add_soil_mesh

  !------------------------------------------------------------------------
  subroutine SetMeshPreliminaries(mesh, name, ncells, nz)

    use MeshType                  , only : mesh_type
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL

    implicit none

    class(mesh_type), pointer   :: mesh
    character(len=*)            :: name
    PetscInt                    :: ncells, nz

    PetscInt        , parameter :: ncells_ghost = 0

    call mesh%Init()
    call mesh%SetName                (name)
    call mesh%SetOrientation         (MESH_ALONG_GRAVITY)
    call mesh%SetID                  (MESH_CLM_SOIL_COL)
    call mesh%SetDimensions          (ncells, ncells_ghost, nz)

  end subroutine SetMeshPreliminaries

  !------------------------------------------------------------------------
  subroutine SetMeshGeometricAttributes(mesh, filter, xc, yc, zc, &
       dx, dy, dz, area, vol)

    use MeshType                  , only : mesh_type
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : VAR_VOLUME

    implicit none

    class(mesh_type), pointer :: mesh
    PetscInt        , pointer :: filter(:)
    PetscReal       , pointer :: xc(:)
    PetscReal       , pointer :: yc(:)
    PetscReal       , pointer :: zc(:)
    PetscReal       , pointer :: dx(:)
    PetscReal       , pointer :: dy(:)
    PetscReal       , pointer :: dz(:)
    PetscReal       , pointer :: area(:)
    PetscReal       , pointer :: vol(:)

    call mesh%SetGridCellFilter      (filter)
    call mesh%SetGeometricAttributes (VAR_XC     , xc)
    call mesh%SetGeometricAttributes (VAR_YC     , yc)
    call mesh%SetGeometricAttributes (VAR_ZC     , zc)
    call mesh%SetGeometricAttributes (VAR_DX     , dx)
    call mesh%SetGeometricAttributes (VAR_DY     , dy)
    call mesh%SetGeometricAttributes (VAR_DZ     , dz)
    call mesh%SetGeometricAttributes (VAR_AREA   , area)
    call mesh%SetGeometricAttributes (VAR_VOLUME , vol)

  end subroutine SetMeshGeometricAttributes

  !------------------------------------------------------------------------

  subroutine add_goveqn(neqns)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : GE_RE
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    !
    ! !ARGUMENTS
    implicit none

    PetscInt, intent(in) :: neqns

    PetscInt :: irank

    do irank = 1, neqns
       call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, 'Richards Equation ODE', irank)
    end do

    call vsfm_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqn

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()

    implicit none

    select case(trim(problem_type))
    case ('oak')
       call add_conditions_to_goveqns_for_single_tree()

    case ('pine')
       call add_conditions_to_goveqns_for_single_tree()

    case ('oak_and_pine')
       call add_conditions_to_goveqns_for_two_trees()

    case default
       write(*,*)'Error while adding condition. Unsupported problem_type : ',trim(problem_type)
       stop

    end select

  end subroutine add_conditions_to_goveqns
  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns_for_single_tree()
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
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn

    ieqn       = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Potential Mass_Flux', 'kg/s', COND_DOWNREG_MASS_RATE_FETCH2, &
         ALL_CELLS)

    ieqn       = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Bottom BC', 'Pa', COND_DIRICHLET, SOIL_BOTTOM_CELLS)

  end subroutine add_conditions_to_goveqns_for_single_tree

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns_for_two_trees()
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
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn
    PetscInt                            :: nconn
    PetscInt                            :: iconn
    PetscInt                  , pointer :: conn_id_up(:)
    PetscInt                  , pointer :: conn_id_dn(:)
    PetscReal                 , pointer :: conn_dist_up(:)
    PetscReal                 , pointer :: conn_dist_dn(:)
    PetscReal                 , pointer :: conn_area(:)
    PetscInt                  , pointer :: conn_type(:)
    PetscReal                 , pointer :: conn_unitvec(:,:)
    class(connection_set_type), pointer :: conn_set

    ieqn       = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Potential Mass_Flux', 'kg/s', COND_DOWNREG_MASS_RATE_FETCH2, &
         ALL_CELLS)

    nconn = 2

    allocate (conn_id_up   (nconn))
    allocate (conn_id_dn   (nconn))
    allocate (conn_dist_up (nconn))
    allocate (conn_dist_dn (nconn))
    allocate (conn_area    (nconn))
    allocate (conn_type    (nconn))
    allocate (conn_unitvec (nconn,3))

    iconn = 1
    conn_id_up(iconn)     = 0
    conn_id_dn(iconn)     = oak_nz
    conn_dist_up(iconn)   = 0.d0
    conn_dist_dn(iconn)   = 0.5d0*dz
    conn_area(iconn)      = oak_Asapwood
    conn_type(iconn)      = CONN_VERTICAL
    conn_unitvec(iconn,1) = 0.d0
    conn_unitvec(iconn,2) = 0.d0
    conn_unitvec(iconn,3) = 1.d0

    iconn = 2
    conn_id_up(iconn)     = 0
    conn_id_dn(iconn)     = oak_nz+pine_nz
    conn_dist_up(iconn)   = 0.d0
    conn_dist_dn(iconn)   = 0.5d0*dz
    conn_area(iconn)      = pine_Asapwood
    conn_type(iconn)      = CONN_VERTICAL
    conn_unitvec(iconn,1) = 0.d0
    conn_unitvec(iconn,2) = 0.d0
    conn_unitvec(iconn,3) = 1.d0

    call MeshCreateConnectionSet(vsfm_mpp%meshes(1), &
         nconn, conn_id_up, conn_id_dn, &
         conn_dist_up, conn_dist_dn, conn_area, &
         conn_type, conn_unitvec, conn_set)

    ieqn       = 1
    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Bottom BC', 'Pa', COND_DIRICHLET, conn_set=conn_set)

    deallocate (conn_id_up   )
    deallocate (conn_id_dn   )
    deallocate (conn_dist_up )
    deallocate (conn_dist_dn )
    deallocate (conn_area    )
    deallocate (conn_type    )
    deallocate (conn_unitvec )

  end subroutine add_conditions_to_goveqns_for_two_trees

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

    implicit none

    select case(trim(problem_type))
    case ('oak')
       call set_material_properties_for_single_tree(   &
            oak_nz, porosity, oak_phi88, oak_phi50,    &
            oak_kmax, oak_c1, oak_c2, oak_c3, oak_phis50)

    case ('pine')
       call set_material_properties_for_single_tree(   &
            pine_nz, porosity, pine_phi88, pine_phi50, &
            pine_kmax, pine_c1, pine_c2, pine_c3, pine_phis50)

    case ('oak_and_pine')
       call set_material_properties_for_two_trees()

    case default
       write(*,*)'Unable to set material prop for problem_type = ' // trim(problem_type)
       stop

    end select

  end subroutine set_material_properties
  !------------------------------------------------------------------------
  subroutine set_material_properties_for_single_tree(local_nz, local_porosity, &
       local_phi88, local_phi50, local_kmax, local_c1, local_c2, local_c3, local_phis50)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSourceSinkAuxVarRealValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_EXPONENT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    !
    implicit none
    !
    PetscInt  :: local_nz
    PetscReal :: local_porosity
    PetscReal :: local_phi88
    PetscReal :: local_phi50
    PetscReal :: local_kmax
    PetscReal :: local_c1
    PetscReal :: local_c2
    PetscReal :: local_c3
    PetscReal :: local_phis50
    !
    PetscReal , pointer   :: por(:)
    PetscReal , pointer   :: alpha(:)
    PetscReal , pointer   :: lambda(:)
    PetscReal , pointer   :: sat_res(:)
    PetscReal , pointer   :: perm(:)
    PetscInt  , pointer   :: vsfm_filter(:)
    PetscInt  , pointer   :: satfunc_type(:)
    PetscInt  , pointer   :: relperm_type(:)
    PetscReal , pointer   :: ss_exponent(:)
    PetscReal , pointer   :: ss_pressure(:)
    PetscReal , pointer   :: weibull_d(:)
    PetscReal , pointer   :: weibull_c(:)
    PetscInt              :: ieqn
    !-----------------------------------------------------------------------

    allocate (por            (local_nz))
    allocate (alpha          (local_nz))
    allocate (lambda         (local_nz))
    allocate (sat_res        (local_nz))
    allocate (satfunc_type   (local_nz))
    allocate (perm           (local_nz))
    allocate (relperm_type   (local_nz))
    allocate (weibull_d      (local_nz))
    allocate (weibull_c      (local_nz))
    allocate (ss_exponent    (local_nz))
    allocate (ss_pressure    (local_nz))

    call set_xylem_material_properties(1, local_nz,       &
         local_porosity, local_phi88, local_phi50,        &
         local_kmax, local_c1, local_c2, local_c3,        &
         local_phis50,                                    &
         por, satfunc_type, alpha, lambda, sat_res, perm, &
         relperm_type, weibull_d, weibull_c, ss_exponent, &
         ss_pressure)

    call VSFMMPPSetSoilPorosity(vsfm_mpp, 1, por)
    call VSFMMPPSetSaturationFunction(vsfm_mpp, 1, satfunc_type, &
         alpha, lambda, sat_res)
    call VSFMMPPSetSoilPermeability(vsfm_mpp, 1, perm, perm, perm)
    call VSFMMPPSetRelativePermeability(vsfm_mpp, 1, relperm_type, weibull_d, weibull_c)

    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_EXPONENT, ss_exponent)

    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_PRESSURE, ss_pressure)

    deallocate(por          )
    deallocate(sat_res      )
    deallocate(lambda       )
    deallocate(alpha        )
    deallocate(perm         )
    deallocate(satfunc_type )
    deallocate(relperm_type )
    deallocate(weibull_c    )
    deallocate(weibull_d    )
    deallocate(ss_exponent  )
    deallocate(ss_pressure  )

  end subroutine set_material_properties_for_single_tree

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_two_trees()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSourceSinkAuxVarRealValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_EXPONENT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
    !
    implicit none
    !
    PetscInt              :: nz
    PetscReal , pointer   :: por(:)
    PetscReal , pointer   :: alpha(:)
    PetscReal , pointer   :: lambda(:)
    PetscReal , pointer   :: sat_res(:)
    PetscReal , pointer   :: perm(:)
    PetscInt  , pointer   :: vsfm_filter(:)
    PetscInt  , pointer   :: satfunc_type(:)
    PetscInt  , pointer   :: relperm_type(:)
    PetscReal , pointer   :: ss_exponent(:)
    PetscReal , pointer   :: ss_pressure(:)
    PetscReal , pointer   :: weibull_d(:)
    PetscReal , pointer   :: weibull_c(:)
    PetscInt              :: ieqn
    !-----------------------------------------------------------------------

    nz = oak_nz + pine_nz
    
    allocate (por            (nz))
    allocate (alpha          (nz))
    allocate (lambda         (nz))
    allocate (sat_res        (nz))
    allocate (satfunc_type   (nz))
    allocate (perm           (nz))
    allocate (relperm_type   (nz))
    allocate (weibull_d      (nz))
    allocate (weibull_c      (nz))
    allocate (ss_exponent    (nz))
    allocate (ss_pressure    (nz))

    call set_xylem_material_properties(1, oak_nz,         &
         porosity, oak_phi88, oak_phi50,                  &
         oak_kmax, oak_c1, oak_c2, oak_c3,                &
         oak_phis50,                                      &
         por, satfunc_type, alpha, lambda, sat_res, perm, &
         relperm_type, weibull_d, weibull_c, ss_exponent, &
         ss_pressure)

    call set_xylem_material_properties(oak_nz+1, nz,      &
         porosity, pine_phi88, pine_phi50,                &
         pine_kmax, pine_c1, pine_c2, pine_c3,            &
         pine_phis50,                                     &
         por, satfunc_type, alpha, lambda, sat_res, perm, &
         relperm_type, weibull_d, weibull_c, ss_exponent, &
         ss_pressure)

    call VSFMMPPSetSoilPorosity(vsfm_mpp, 1, por)

    call VSFMMPPSetSaturationFunction(vsfm_mpp, 1, satfunc_type, &
         alpha, lambda, sat_res)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, 1, perm, perm, perm)
    
    call VSFMMPPSetRelativePermeability(vsfm_mpp, 1, relperm_type, weibull_d, weibull_c)

    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_EXPONENT, ss_exponent)

    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, 1, &
         VAR_POT_MASS_SINK_PRESSURE, ss_pressure)

    deallocate(por          )
    deallocate(sat_res      )
    deallocate(lambda       )
    deallocate(alpha        )
    deallocate(perm         )
    deallocate(satfunc_type )
    deallocate(relperm_type )
    deallocate(weibull_c    )
    deallocate(weibull_d    )
    deallocate(ss_exponent  )
    deallocate(ss_pressure  )

  end subroutine set_material_properties_for_two_trees

  !------------------------------------------------------------------------
  subroutine set_xylem_material_properties(ibeg, iend,     &
         local_porosity, local_phi88, local_phi50,         &
         local_kmax, local_c1, local_c2, local_c3,         &
         local_phis50,                                     &
         por, satfunc_type, alpha, lambda, sat_res, perm,  &
         relperm_type, weibull_d, weibull_c, ss_exponent,  &
         ss_pressure)

    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL

    implicit none

    PetscInt            :: ibeg, iend
    PetscReal           :: local_porosity
    PetscReal           :: local_phi88
    PetscReal           :: local_phi50
    PetscReal           :: local_kmax
    PetscReal           :: local_c1
    PetscReal           :: local_c2
    PetscReal           :: local_c3
    PetscReal           :: local_phis50
    PetscReal , pointer :: por(:)
    PetscReal , pointer :: alpha(:)
    PetscReal , pointer :: lambda(:)
    PetscReal , pointer :: sat_res(:)
    PetscReal , pointer :: perm(:)
    PetscInt  , pointer :: vsfm_filter(:)
    PetscInt  , pointer :: satfunc_type(:)
    PetscInt  , pointer :: relperm_type(:)
    PetscReal , pointer :: ss_exponent(:)
    PetscReal , pointer :: ss_pressure(:)
    PetscReal , pointer :: weibull_d(:)
    PetscReal , pointer :: weibull_c(:)

    por          (ibeg:iend ) = local_porosity

    satfunc_type (ibeg:iend ) = SAT_FUNC_FETCH2
    alpha        (ibeg:iend ) = local_phi88
    lambda       (ibeg:iend ) = local_phi50
    sat_res      (ibeg:iend ) = 0.d0

    perm         (ibeg:iend ) = local_kmax  * vis / rho! * 1.125d0
    relperm_type (ibeg:iend ) = RELPERM_FUNC_WEIBULL
    weibull_d    (ibeg:iend ) = local_c1
    weibull_c    (ibeg:iend ) = local_c2
    ss_exponent  (ibeg:iend ) = local_c3
    ss_pressure  (ibeg:iend ) = local_phis50

  end subroutine set_xylem_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()

    implicit none

    select case(trim(problem_type))
    case ('oak')
       call set_initial_conditions_for_single_tree(oak_nz)

    case ('pine')
       call set_initial_conditions_for_single_tree(pine_nz)

    case ('oak_and_pine')
       call set_initial_conditions_for_two_trees()

    case default
       write(*,*)'Unable to set IC for problem_type = ' // trim(problem_type)
       stop

    end select
  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_initial_conditions_for_single_tree(nz)
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
    PetscInt           :: nz
    !
    PetscReal          :: theta
    PetscReal          :: phi_root_mean_times_beta_s
    PetscInt           :: ii
    PetscReal, pointer :: press_ic(:)
    PetscErrorCode     :: ierr

    Vec                :: p_init
    PetscViewer        :: viewer
    PetscReal, pointer :: p_init_ptr(:)

    phi_root_mean_times_beta_s = 5.831916333333334e+03

    if (ic_file_specified) then
       call PetscViewerBinaryOpen(PETSC_COMM_WORLD, ic_filename, FILE_MODE_READ, &
            viewer, ierr); CHKERRQ(ierr)
       ! Load the data
       call VecCreate(PETSC_COMM_WORLD, p_init, ierr); CHKERRQ(ierr)
       call VecLoad(p_init, viewer, ierr); CHKERRQ(ierr)
       call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

       call VecCopy(p_init, vsfm_mpp%soe%solver%soln, ierr); CHKERRQ(ierr)
       call VecGetArrayF90(p_init, p_init_ptr, ierr); CHKERRQ(ierr)
       call vsfm_mpp%Restart(p_init_ptr)
       call VecRestoreArrayF90(p_init, p_init_ptr, ierr); CHKERRQ(ierr)
    else
       call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       do ii = 1, nz
          press_ic(ii) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (nz - ii)*dz) + 101325.d0
       enddo
       call vsfm_mpp%Restart(press_ic)
       call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
    end if

    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)
    
  end subroutine set_initial_conditions_for_single_tree

  !------------------------------------------------------------------------
  subroutine set_initial_conditions_for_two_trees()
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
    PetscReal          :: theta
    PetscReal          :: phi_root_mean_times_beta_s
    PetscInt           :: ii
    PetscReal, pointer :: press_ic(:)
    PetscErrorCode     :: ierr

    Vec                :: p_init
    PetscViewer        :: viewer
    PetscReal, pointer :: p_init_ptr(:)

    phi_root_mean_times_beta_s = 5.831916333333334e+03
    
    if (ic_file_specified) then
       call PetscViewerBinaryOpen(PETSC_COMM_WORLD, ic_filename, FILE_MODE_READ, &
            viewer, ierr); CHKERRQ(ierr)
       ! Load the data
       call VecCreate(PETSC_COMM_WORLD, p_init, ierr); CHKERRQ(ierr)
       call VecLoad(p_init, viewer, ierr); CHKERRQ(ierr)
       call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

       call VecCopy(p_init, vsfm_mpp%soe%solver%soln, ierr); CHKERRQ(ierr)
       call VecGetArrayF90(p_init, p_init_ptr, ierr); CHKERRQ(ierr)
       call vsfm_mpp%Restart(p_init_ptr)
       call VecRestoreArrayF90(p_init, p_init_ptr, ierr); CHKERRQ(ierr)
    else
       call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       do ii = 1, oak_nz
          press_ic(ii         ) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (oak_nz - ii)*dz) + 101325.d0
       enddo

       do ii = 1, pine_nz
          press_ic(ii + oak_nz) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (pine_nz - ii)*dz) + 101325.d0
       enddo
       call vsfm_mpp%Restart(press_ic)
       call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
    end if

    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

#if 0
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'initial_pressure.bin', FILE_MODE_READ, viewer, ierr); CHKERRQ(ierr)
    call VecCreate(PETSC_COMM_WORLD, p_init, ierr); CHKERRQ(ierr)
    call VecLoad(p_init, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(p_init, p_init_ptr, ierr); CHKERRQ(ierr)

    do ii = 1, oak_nz
       press_ic(ii         ) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (oak_nz - ii)*dz) + 101325.d0
    enddo

    do ii = 1, pine_nz
       press_ic(ii + oak_nz) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (pine_nz - ii)*dz) + 101325.d0
    enddo

    call vsfm_mpp%Restart(press_ic)
    call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)

    call VecRestoreArrayF90(p_init, p_init_ptr, ierr); CHKERRQ(ierr)

    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)
#endif
    
  end subroutine set_initial_conditions_for_two_trees

  !------------------------------------------------------------------------
  subroutine set_boundary_conditions_for_single_tree(nz, nstep, ET, SoilBC)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC, VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nz
    PetscInt           :: nstep
    !
    PetscReal, pointer :: ss_value(:)
    PetscReal, pointer :: bc_value(:)
    PetscViewer        :: viewer
    Vec                :: ET
    Vec                :: SoilBC
    PetscReal, pointer :: et_p(:)
    PetscReal, pointer :: soilbc_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk

    allocate(ss_value(nz))

    call VecGetArrayF90(ET, et_p, ierr)
    do kk = 1, nz
       ss_value(kk) = -et_p(nz*(nstep-1) + kk) * dz
    end do
    call VecRestoreArrayF90(ET, et_p, ierr)

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

    deallocate(ss_value)

    call VecGetArrayF90(SoilBC, soilbc_p, ierr)
    allocate(bc_value(1))
    bc_value(1) = soilbc_p(nstep)
    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bc_value)
    deallocate(bc_value)
    call VecRestoreArrayF90(SoilBC, soilbc_p, ierr)

  end subroutine set_boundary_conditions_for_single_tree

  !------------------------------------------------------------------------
  subroutine set_boundary_conditions_for_two_trees(nz, nstep, ET, SoilBC)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC, VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nz
    PetscInt           :: nstep
    !
    PetscReal, pointer :: ss_value(:)
    PetscReal, pointer :: bc_value(:)
    PetscViewer        :: viewer
    Vec                :: ET
    Vec                :: SoilBC
    PetscReal, pointer :: et_p(:)
    PetscReal, pointer :: soilbc_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk

    allocate(ss_value(nz))

    ! Set source sink
    call VecGetArrayF90(ET, et_p, ierr)
    
    do kk = 1, nz
       ss_value(kk) = -et_p(nz*(nstep-1) + kk)*dz
    end do
    call VecRestoreArrayF90(ET, et_p, ierr)

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

    deallocate(ss_value)

    ! Set BC
    call VecGetArrayF90(SoilBC, soilbc_p, ierr)
    allocate(bc_value(2))
    bc_value(1:2) = soilbc_p(nstep)
        soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bc_value)
    deallocate(bc_value)
    call VecRestoreArrayF90(SoilBC, soilbc_p, ierr)

  end subroutine set_boundary_conditions_for_two_trees

  !------------------------------------------------------------------------
  subroutine diagnose_actual_sink_for_single_tree(nz, nstep, act_et_p)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS_FLUX
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nz
    PetscInt           :: nstep
    PetscReal, pointer :: act_et_p(:)
    !
    PetscReal, pointer :: ss_value(:)
    PetscReal, pointer :: press(:)
    PetscInt           :: kk

    allocate(ss_value(nz))
    allocate(press(nz))

    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_SS, VAR_MASS_FLUX, 1, ss_value)
    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,VAR_PRESSURE, -1, press)
    do kk = 1, nz
       act_et_p((nstep-1)*nz + kk) = ss_value(kk)
    end do
    
    deallocate(ss_value)
    deallocate(press)

  end subroutine diagnose_actual_sink_for_single_tree

end module vsfm_spac_fetch2_problem
