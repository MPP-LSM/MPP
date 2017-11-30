module spac_component
  
#include <petsc/finclude/petsc.h>
  
  implicit none

  type, public :: spac_component_mesh_type
     PetscInt  , pointer :: id(:,:,:)
     PetscReal , pointer :: xc(:)                          !
     PetscReal , pointer :: yc(:)                          !
     PetscReal , pointer :: zc(:)                          !
     PetscReal , pointer :: area(:)                        !
     PetscReal , pointer :: vol(:)                         !
     PetscInt  , pointer :: filter(:)                      !
  contains
     procedure, public :: Copy => MeshCopy
  end type spac_component_mesh_type

  type, public :: spac_component_pp_type
     PetscReal , pointer :: por(:)
     PetscReal , pointer :: perm(:)
     PetscReal , pointer :: lambda(:)
     PetscReal , pointer :: alpha(:)
     PetscReal , pointer :: eff_porosity(:)
     PetscReal , pointer :: residual_sat(:)
     PetscInt  , pointer :: satfunc_type(:)
     PetscInt  , pointer :: relperm_type(:)
     PetscReal , pointer :: relperm_param_1(:)
     PetscReal , pointer :: relperm_param_2(:)
   contains
     procedure, public :: Copy => PPCopy
  end type spac_component_pp_type

contains

  subroutine MeshCopy(this, idx_beg, idx_end, xc, yc, zc, area, vol, filter)
    !
    implicit none
    !
    class (spac_component_mesh_type) :: this
    PetscInt                         :: idx_beg
    PetscInt                         :: idx_end
    PetscReal , pointer              :: xc(:)
    PetscReal , pointer              :: yc(:)
    PetscReal , pointer              :: zc(:)
    PetscReal , pointer              :: area(:)
    PetscReal , pointer              :: vol(:)
    PetscInt  , pointer              :: filter(:)

    xc     (idx_beg:idx_end) = this%xc     (:)
    yc     (idx_beg:idx_end) = this%yc     (:)
    zc     (idx_beg:idx_end) = this%zc     (:)
    area   (idx_beg:idx_end) = this%area   (:)
    vol    (idx_beg:idx_end) = this%vol    (:)
    filter (idx_beg:idx_end) = this%filter (:)

  end subroutine MeshCopy

  subroutine PPCopy(this, idx_beg, idx_end, por, perm, alpha, lambda, &
       relperm_type, relperm_param_1, relperm_param_2, residual_sat, &
       satfunc_type)
    !
    implicit none
    !
    class (spac_component_pp_type) :: this
    PetscInt                       :: idx_beg
    PetscInt                       :: idx_end
    PetscReal , pointer            :: por(:)
    PetscReal , pointer            :: perm(:)
    PetscReal , pointer            :: alpha(:)
    PetscReal , pointer            :: lambda(:)
    PetscInt  , pointer            :: relperm_type(:)
    PetscReal , pointer            :: relperm_param_1(:)
    PetscReal , pointer            :: relperm_param_2(:)
    PetscReal , pointer            :: residual_sat(:)
    PetscInt  , pointer            :: satfunc_type(:)

    por             (idx_beg:idx_end) = this%por             (:)
    perm            (idx_beg:idx_end) = this%perm            (:)
    alpha           (idx_beg:idx_end) = this%alpha           (:)
    lambda          (idx_beg:idx_end) = this%lambda          (:)
    relperm_type    (idx_beg:idx_end) = this%relperm_type    (:)
    relperm_param_1 (idx_beg:idx_end) = this%relperm_param_1 (:)
    relperm_param_2 (idx_beg:idx_end) = this%relperm_param_2 (:)
    residual_sat    (idx_beg:idx_end) = this%residual_sat    (:)
    satfunc_type    (idx_beg:idx_end) = this%satfunc_type    (:)

  end subroutine PPCopy

end module spac_component

module soil_parameters

#include <petsc/finclude/petsc.h>
  
  implicit none

  PetscInt  , parameter :: soil_nx                          = 2          !
  PetscInt  , parameter :: soil_ny                          = 1          !
  PetscInt  , parameter :: soil_nz                          = 20         !
  PetscReal , parameter :: soil_dx                          = 10.d0      ! [m]
  PetscReal , parameter :: soil_dy                          = 10.d0      ! [m]
  PetscReal , parameter :: soil_dz                          = 0.25d0     ! [m]
  PetscReal             :: slope                                         ! [m m^{-1}]

  PetscReal, parameter  :: perm_xy_top                       = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: perm_z_top                        = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: sat_res_top                       = 0.06d0    ! [-]
  PetscReal, parameter  :: alpha_top                         = 0.00005d0 ! [Pa^{-1}]
  PetscReal, parameter  :: vg_m_top                          = 0.33d0    ! [-]
  PetscReal, parameter  :: por_top                           = 0.5d0     ! [-]

  PetscReal, parameter  :: perm_xy_mid                       = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: perm_z_mid                        = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: sat_res_mid                       = 0.06d0    ! [-]
  PetscReal, parameter  :: alpha_mid                         = 0.00005d0 ! [Pa^{-1}]
  PetscReal, parameter  :: vg_m_mid                          = 0.33d0    ! [-]
  PetscReal, parameter  :: por_mid                           = 0.5d0     ! [-]

  PetscReal, parameter  :: perm_xy_bot                       = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: perm_z_bot                        = 6.83d-11  ! [m^2]
  PetscReal, parameter  :: sat_res_bot                       = 0.06d0    ! [-]
  PetscReal, parameter  :: alpha_bot                         = 0.00005d0 ! [Pa^{-1}]
  PetscReal, parameter  :: vg_m_bot                          = 0.33d0    ! [-]
  PetscReal, parameter  :: por_bot                           = 0.5d0     ! [-]
  
  PetscInt              :: soil_ncells                                   ! Total number of active grid cells
  PetscInt              :: soil_mesh_nconn                               ! Total number of all (vertical + horizontal) connections

  PetscBool             :: is_soil_horizontally_disconnected             ! Are soil grid cells horizontally disconnected?

  PetscReal , pointer   :: soil_xc(:)                                    ! x-position of grid cell [m]
  PetscReal , pointer   :: soil_yc(:)                                    ! y-position of grid cell [m]
  PetscReal , pointer   :: soil_zc(:)                                    ! y-position of grid cell [m]
  PetscReal , pointer   :: soil_soil_dx(:)                               ! x-position of grid cell [m]
  PetscReal , pointer   :: soil_soil_dy(:)                               ! y-position of grid cell [m
  PetscReal , pointer   :: soil_soil_dz(:)                               ! y-position of grid cell [m]
  PetscReal , pointer   :: soil_area(:)                                  ! area of grid cell [m^2]
  PetscReal , pointer   :: soil_vol(:)                                   ! volume of grid cell [m^3]
  PetscInt  , pointer   :: soil_filter(:)                                ! 

  PetscInt  , pointer   :: soil_id(:,:,:)                                ! [-] Id > 0 indicates active grid cell, otherwise inactive
  PetscReal , pointer   :: soil_xc3d(:,:,:)                              ! [m]
  PetscReal , pointer   :: soil_yc3d(:,:,:)                              ! [m]
  PetscReal , pointer   :: soil_zc3d(:,:,:)                              ! [m]

  PetscInt  , pointer   :: top_active_layer_kk_index(:,:)                ! [-] Index of the first active grid cell in z-dimension
  PetscReal , pointer   :: elevation(:,:)                                ! [m]

  PetscReal , pointer :: soil_por(:)
  PetscReal , pointer :: soil_perm(:)
  PetscReal , pointer :: soil_lambda(:)
  PetscReal , pointer :: soil_alpha(:)
  PetscReal , pointer :: soil_eff_porosity(:)
  PetscReal , pointer :: soil_residual_sat(:)
  PetscInt  , pointer :: soil_satfunc_type(:)
  PetscInt  , pointer :: soil_relperm_type(:)
  PetscReal , pointer :: soil_relperm_param_1(:)
  PetscReal , pointer :: soil_relperm_param_2(:)

end module soil_parameters

module overstory_parameters
  
#include <petsc/finclude/petsc.h>

  use spac_component
  
  implicit none

  PetscReal , parameter :: overstory_root_depth             = 2.0d0     ! [m]
  PetscReal , parameter :: overstory_canopy_height          = 17.d0     ! [m]
  PetscReal , parameter :: overstory_area_canopy            = 12.25d0   ! [m^2]
  PetscReal , parameter :: overstory_area_sapwood           = 0.013d0   ! [m^2]
  PetscReal , parameter :: overstory_LAI                    = 5.52d0    ! [m^2 m^{-2}]
  PetscReal , parameter :: overstory_area_tappering_coeff   = 0.75d0    ! [-]
  PetscReal , parameter :: overstory_branch_length          = 1.75d0    ! [m]
  PetscReal , parameter :: overstory_branch_area_ratio      = 0.15d0    ! [m^2 branch m^{-2} sapwood]
  PetscReal , parameter :: overstory_root_radius            = 2.9d-4    ! [m]
  PetscReal , parameter :: overstory_root_conductance       = 3.d-11    ! [s^{-1}]
  PetscReal , parameter :: overstory_xylem_Kmax             = 2.5d-5    ! [m s^{-1}]
  PetscReal , parameter :: overstory_xylem_vulnerability_c  = 3.5d0     ! [-]
  PetscReal , parameter :: overstory_xylem_vulnerability_d  = 480.d0    ! [m]
  PetscReal , parameter :: overstory_xylem_phi0             = -2.87d6   ! [Pa]
  PetscReal , parameter :: overstory_xylem_p                = 100.0d0   ! [-]
  PetscReal , parameter :: overstory_xylem_porosity         = 0.57d0    ! [-]
  PetscReal , parameter :: overstory_lad_profile(68)        = (/ &
       0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  &
       0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  &
       0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  &
       0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,  &
       0.00d0, 0.00d0, 0.01d0, 0.03d0, 0.04d0,  &
       0.05d0, 0.06d0, 0.07d0, 0.08d0, 0.09d0,  &
       0.13d0, 0.21d0, 0.30d0, 0.38d0, 0.46d0,  &
       0.52d0, 0.59d0, 0.65d0, 0.71d0, 0.74d0,  &
       0.78d0, 0.81d0, 0.84d0, 0.85d0, 0.85d0,  &
       0.84d0, 0.84d0, 0.83d0, 0.81d0, 0.79d0,  &
       0.77d0, 0.74d0, 0.72d0, 0.69d0, 0.67d0,  &
       0.64d0, 0.61d0, 0.57d0, 0.54d0, 0.51d0,  &
       0.46d0, 0.42d0, 0.37d0, 0.32d0, 0.27d0,  &
       0.20d0, 0.13d0, 0.05d0  &
       /)
  PetscReal , parameter :: overstory_B_profile(8) = (/ &
       45.73d0, 42.82d0, 43.02d0, 39.23d0, 36.80d0,  &
       36.52d0, 21.94d0, 22.83d0  &
       /)

  PetscInt :: overstory_xylem_nz                                        !
  PetscInt :: overstory_leaf_nz                                         !
  PetscInt :: overstory_root_nz                                         !
  PetscInt :: overstory_nconn

  PetscReal , pointer :: overstory_xylem_area_profile(:)                !
  PetscReal , pointer :: overstory_branch_length_profile(:)             !
  PetscBool , pointer :: overstory_xylem_has_branches(:)                !
  PetscReal , pointer :: overstory_root_length_profile(:)
  PetscReal , pointer :: overstory_root_area_profile(:)
  PetscReal , pointer :: overstory_root_vol_profile(:)
  PetscInt  , pointer :: overstory_branch_2_xylem_index(:)

  type (spac_component_mesh_type) :: overstory_root_mesh
  type (spac_component_mesh_type) :: overstory_xylem_mesh
  type (spac_component_mesh_type) :: overstory_leaf_mesh

  type (spac_component_pp_type) :: overstory_root_pp
  type (spac_component_pp_type) :: overstory_xylem_pp
  type (spac_component_pp_type) :: overstory_leaf_pp

end module overstory_parameters

module understory_parameters
  !
#include <petsc/finclude/petsc.h>
  !
  use spac_component
  !
  implicit none
  !
  PetscReal , parameter :: understory_root_depth            = 0.5d0     ! [m]
  PetscReal , parameter :: understory_canopy_height         = 1.d0      ! [m]
  PetscReal , parameter :: understory_area_canopy           = 9.0d0     ! [m^2]
  PetscReal , parameter :: understory_area_sapwood          = 0.010d0   ! [m^2]
  PetscReal , parameter :: understory_LAI                   = 3.19d0    ! [m^2 m^{-2}]
  PetscReal , parameter :: understory_area_tappering_coeff  = 0.75d0    ! [-]
  PetscReal , parameter :: understory_branch_length         = 1.50d0    ! [m]
  PetscReal , parameter :: understory_branch_area_ratio     = 0.15d0    ! [m^2 branch m^{-2} sapwood]
  PetscReal , parameter :: understory_root_radius           = 2.9d-4    ! [m]
  PetscReal , parameter :: understory_root_conductance      = 3.d-11    ! [s^{-1}]
  PetscReal , parameter :: understory_xylem_Kmax            = 2.5d-5    ! [m s^{-1}]
  PetscReal , parameter :: understory_xylem_vulnerability_c = 3.5d0     ! [-]
  PetscReal , parameter :: understory_xylem_vulnerability_d = 480.d0    ! [m]
  PetscReal , parameter :: understory_xylem_phi0            = -2.87d6   ! [Pa]
  PetscReal , parameter :: understory_xylem_p               = 100.0d0   ! [-]
  PetscReal , parameter :: understory_xylem_porosity        = 0.57d0    ! [-]
  PetscReal , parameter :: understory_lad_profile(24) = (/ &
       0.00d0, 0.07d0, 0.21d0, 0.35d0, 0.49d0,  &
       0.54d0, 0.57d0, 0.61d0, 0.64d0, 0.66d0,  &
       0.67d0, 0.69d0, 0.70d0, 0.70d0, 0.69d0,  &
       0.68d0, 0.66d0, 0.65d0, 0.61d0, 0.58d0,  &
       0.54d0, 0.50d0, 0.39d0, 0.28d0  &
       /)
  PetscReal , parameter :: understory_B_profile(2) = (/ &
       0.76d0, 0.16d0  &
       /)

  PetscInt :: understory_xylem_nz                                       !
  PetscInt :: understory_leaf_nz                                        !
  PetscInt :: understory_root_nz                                        !
  PetscInt :: understory_nconn
  
  PetscReal , pointer :: understory_xylem_area_profile(:)               !
  PetscReal , pointer :: understory_branch_length_profile(:)            !
  PetscBool , pointer :: understory_xylem_has_branches(:)               !
  PetscReal , pointer :: understory_root_length_profile(:)              !
  PetscReal , pointer :: understory_root_area_profile(:)                !
  PetscReal , pointer :: understory_root_vol_profile(:)                 !
  PetscInt  , pointer :: understory_branch_2_xylem_index(:)

  type (spac_component_mesh_type) :: understory_root_mesh
  type (spac_component_mesh_type) :: understory_xylem_mesh
  type (spac_component_mesh_type) :: understory_leaf_mesh

  type (spac_component_pp_type) :: understory_root_pp
  type (spac_component_pp_type) :: understory_xylem_pp
  type (spac_component_pp_type) :: understory_leaf_pp

end module understory_parameters

module problem_parameters

  use soil_parameters
  use overstory_parameters
  use understory_parameters

  implicit none

  PetscReal , parameter :: PI       = 4 * atan (1.0_8) ! [-]
  PetscReal , parameter :: vish2o   = 0.001002d0       ! [N s/m^2] @ 20 degC

  PetscReal , parameter :: init_wtd = 3.d0             ! Initial water table depth below surface [m]
  
  PetscInt              :: ncells

  PetscBool             :: multi_goveqns_formulation

  PetscReal , pointer   :: vsfm_por(:)
  PetscReal , pointer   :: vsfm_perm(:)
  PetscReal , pointer   :: vsfm_lambda(:)
  PetscReal , pointer   :: vsfm_alpha(:)
  PetscReal , pointer   :: vsfm_eff_porosity(:)
  PetscReal , pointer   :: vsfm_residual_sat(:)
  PetscInt  , pointer   :: vsfm_satfunc_type(:)
  PetscInt  , pointer   :: vsfm_relperm_type(:)
  PetscReal , pointer   :: vsfm_relperm_param_1(:)
  PetscReal , pointer   :: vsfm_relperm_param_2(:)

  PetscInt  , pointer   :: conn_flux_type(:)    !
  PetscInt  , pointer   :: conn_cond_type(:)    !
  PetscReal , pointer   :: conn_cond     (:)    !
  PetscReal , pointer   :: conn_cond_up  (:)    !
  PetscReal , pointer   :: conn_cond_dn  (:)    !
  
  PetscInt  , pointer   :: conn_satparam_up_itype(:)
  PetscReal , pointer   :: conn_satparam_up_param_1(:)
  PetscReal , pointer   :: conn_satparam_up_param_2(:)
  PetscReal , pointer   :: conn_satparam_up_param_3(:)
  PetscInt  , pointer   :: conn_relperm_up_itype(:)
  PetscReal , pointer   :: conn_relperm_up_param_1(:)
  PetscReal , pointer   :: conn_relperm_up_param_2(:)

  PetscInt  , pointer   :: conn_satparam_dn_itype(:)
  PetscReal , pointer   :: conn_satparam_dn_param_1(:)
  PetscReal , pointer   :: conn_satparam_dn_param_2(:)
  PetscReal , pointer   :: conn_satparam_dn_param_3(:)
  PetscInt  , pointer   :: conn_relperm_dn_itype(:)
  PetscReal , pointer   :: conn_relperm_dn_param_1(:)
  PetscReal , pointer   :: conn_relperm_dn_param_2(:)
  PetscReal , pointer   :: conn_relperm_dn_param_3(:)

end module problem_parameters

!------------------------------------------------------------------------
program vsfm_spac_on_hillslope

#include <petsc/finclude/petsc.h>

  use petscsys
  !
  implicit none
  !
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call run_vsfm_spac_on_hillslope()

  call PetscFinalize(ierr)

end program vsfm_spac_on_hillslope

!------------------------------------------------------------------------
subroutine run_vsfm_spac_on_hillslope()
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
  use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
  use problem_parameters
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
  PetscBool      :: converged
  PetscInt       :: converged_reason
  PetscErrorCode :: ierr
  PetscReal      :: time
  PetscReal      :: dtime
  PetscInt       :: nstep
  PetscInt       :: istep
  PetscViewer    :: viewer
  Vec            :: sat
  PetscReal, pointer :: sat_p(:)
  PetscReal, pointer :: data(:)
  PetscBool          :: flg

  ! Initialize data
  time  = 0.d0

  ! Set default settings
  dtime                             = 180.d0
  nstep                             = 1
  slope                             = 0.05d0
  is_soil_horizontally_disconnected = PETSC_FALSE
  multi_goveqns_formulation         = PETSC_FALSE
  
  call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-slope',slope,flg,ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-soil_horizontally_disconnected',is_soil_horizontally_disconnected,flg,ierr)

  call initialize_problem()

  call vsfm_mpp%soe%SetDtime(1.d0)
  call vsfm_mpp%soe%PreSolve()
  call vsfm_mpp%soe%PostSolve()

  allocate(data(soil_ncells))
  call VecDuplicate(vsfm_mpp%soe%solver%soln, sat, ierr); CHKERRQ(ierr)

  call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL, VAR_LIQ_SAT, -1, data)

  call VecGetArrayF90(sat, sat_p, ierr); CHKERRQ(ierr)
  sat_p(:) = data(:)
  call VecRestoreArrayF90(sat, sat_p, ierr); CHKERRQ(ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_SELF, 'sat_init.bin', FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
  call VecView(sat, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  do istep = 1, nstep

     !call set_bondary_conditions(istep)
     time = time + dtime

     ! Run the model
     call vsfm_mpp%soe%StepDT(dtime, istep, &
          converged, converged_reason, ierr); CHKERRQ(ierr)

  end do

  call PetscViewerBinaryOpen(PETSC_COMM_SELF, 'pressure_final.bin', FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
  call VecView(vsfm_mpp%soe%solver%soln, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL, VAR_LIQ_SAT, -1, data)

  call VecGetArrayF90(sat, sat_p, ierr); CHKERRQ(ierr)
  sat_p(:) = data(:)
  call VecRestoreArrayF90(sat, sat_p, ierr); CHKERRQ(ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_SELF, 'sat_final.bin', FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
  call VecView(sat, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  deallocate(data)
  call VecDestroy(sat, ierr)

end subroutine run_vsfm_spac_on_hillslope

!------------------------------------------------------------------------
subroutine initialize_problem()
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  !
  implicit none

  ! 1. Initialize the multi-physics-problem (MPP)
  call initialize_mpp()
  
  ! 2. Add all meshes needed for the MPP
  call add_mesh()

  ! 3. Add all governing equations
  call add_goveqn()

  ! 4. Add boundary and source-sink conditions to all governing equations
  call add_conditions_to_goveqns()

  ! 5. Allocate memory to hold auxvars
  call allocate_auxvars()

  ! 6. Setup the MPP
  call vsfm_mpp%SetupProblem()

  ! 7. Add material properities associated with all governing equations
  call set_material_properties() 

  ! 8. Set connection flux type
  call set_conn_flux_type()

  ! 9. Set initial conditions
  call set_initial_conditions()

end subroutine initialize_problem

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
  call vsfm_mpp%SetName    ('Tree hydrosoil_dynamics on a hillslope')
  call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
  call vsfm_mpp%SetMPIRank (iam)

end subroutine initialize_mpp

!------------------------------------------------------------------------
subroutine add_mesh()
  !
  use problem_parameters
  !
  implicit none
  !

  call setup_soil_mesh()
  call setup_overstory_mesh()
  call setup_understory_mesh()

  if (.not. multi_goveqns_formulation) then
     call add_single_mesh()
  else
     !call add_mulitple_meshes()
  end if

end subroutine add_mesh

!------------------------------------------------------------------------
subroutine add_single_mesh()
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnRealValue 
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnIntValue
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunctionAuxVarConn
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
  use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
  use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
  use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
  use MultiPhysicsProbConstants , only : VAR_FLUX_TYPE
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
  use MultiPhysicsProbConstants , only : CONDUCTANCE_CAMPBELL_TYPE    
  use MultiPhysicsProbConstants , only : CONDUCTANCE_MANOLI_TYPE    
  use SaturationFunction        , only : SAT_FUNC_VAN_GENUCHTEN
  use SaturationFunction        , only : SAT_FUNC_CHUANG
  use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
  use SaturationFunction        , only : RELPERM_FUNC_MUALEM
  use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
  use mpp_varcon                , only : denh2o
  use mpp_varcon                , only : grav
  use problem_parameters
  !
  implicit none
  !
#include <petsc/finclude/petsc.h>
  !
  PetscInt            :: imesh
  PetscInt            :: ii
  PetscInt            :: jj
  PetscInt            :: kk
  PetscInt            :: idx
  PetscInt            :: count
  PetscInt            :: nlev
  PetscInt            :: iconn
  PetscInt            :: soil_ncells_ghost

  PetscInt  , pointer :: conn_id_up(:)   !
  PetscInt  , pointer :: conn_id_dn(:)   !
  PetscReal , pointer :: conn_dist_up(:) !
  PetscReal , pointer :: conn_dist_dn(:) !
  PetscReal , pointer :: conn_area(:)    !
  PetscInt  , pointer :: conn_type(:)    !

  PetscReal , pointer :: xc(:)                ! x-position of grid cell [m]
  PetscReal , pointer :: yc(:)                ! y-position of grid cell [m]
  PetscReal , pointer :: zc(:)                ! y-position of grid cell [m]
  PetscReal , pointer :: dx(:)                ! [m]
  PetscReal , pointer :: dy(:)                ! [m]
  PetscReal , pointer :: dz(:)                ! [m]
  PetscReal , pointer :: area(:)              ! area of grid cell [m^2]
  PetscReal , pointer :: vol(:)               ! volume of grid cell [m^3]
  PetscInt  , pointer :: filter(:)            ! 

  PetscInt :: nconn
  PetscInt :: idx_beg, idx_end
  
  PetscErrorCode      :: ierr

  ncells = soil_ncells + &
       soil_nx * soil_ny *( overstory_root_nz  + overstory_xylem_nz  + overstory_leaf_nz) + &
       soil_nx * soil_ny *( understory_root_nz + understory_xylem_nz + understory_leaf_nz)

  allocate(xc     (ncells))
  allocate(yc     (ncells))
  allocate(zc     (ncells))
  allocate(dx     (ncells))
  allocate(dy     (ncells))
  allocate(dz     (ncells))
  allocate(area   (ncells))
  allocate(vol    (ncells))
  allocate(filter (ncells))

  allocate(vsfm_por             (ncells)); vsfm_por             (:) = 0.d0
  allocate(vsfm_perm            (ncells)); vsfm_perm            (:) = 0.d0
  allocate(vsfm_lambda          (ncells)); vsfm_lambda          (:) = 0.d0
  allocate(vsfm_alpha           (ncells)); vsfm_alpha           (:) = 0.d0
  allocate(vsfm_eff_porosity    (ncells)); vsfm_eff_porosity    (:) = 0.d0
  allocate(vsfm_residual_sat    (ncells)); vsfm_residual_sat    (:) = 0.d0
  allocate(vsfm_satfunc_type    (ncells)); vsfm_satfunc_type    (:) = 0
  allocate(vsfm_relperm_type    (ncells)); vsfm_relperm_type    (:) = 0
  allocate(vsfm_relperm_param_1 (ncells)); vsfm_relperm_param_1 (:) = 0.d0
  allocate(vsfm_relperm_param_2 (ncells)); vsfm_relperm_param_2 (:) = 0.d0

  ! Soil properties
  vsfm_por          (1:soil_ncells) = por_top
  vsfm_perm         (1:soil_ncells) = perm_z_top
  vsfm_lambda       (1:soil_ncells) = vg_m_top
  vsfm_alpha        (1:soil_ncells) = alpha_top
  vsfm_residual_sat (1:soil_ncells) = sat_res_top
  vsfm_satfunc_type (1:soil_ncells) = SAT_FUNC_VAN_GENUCHTEN
  vsfm_relperm_type (1:soil_ncells) = RELPERM_FUNC_MUALEM

  !
  ! Add soil grid cells
  !
  xc(1:soil_ncells)     = soil_xc(1:soil_ncells)
  yc(1:soil_ncells)     = soil_yc(1:soil_ncells)
  zc(1:soil_ncells)     = soil_zc(1:soil_ncells)
  dx(1:soil_ncells)     = soil_soil_dx(1:soil_ncells)
  dy(1:soil_ncells)     = soil_soil_dy(1:soil_ncells)
  dz(1:soil_ncells)     = soil_soil_dz(1:soil_ncells)
  area(1:soil_ncells)   = soil_area(1:soil_ncells)
  vol(1:soil_ncells)    = soil_vol(1:soil_ncells)
  filter(1:soil_ncells) = soil_filter(1:soil_ncells)

  ncells = soil_ncells

  !
  ! Overstory
  !

  ! Add root grid cells

  idx_beg = ncells + 1
  idx_end = ncells + soil_nx*soil_ny*overstory_root_nz
  ncells  = idx_end

  call overstory_root_mesh%copy(idx_beg, idx_end, xc, yc, zc, area, vol, filter)

  call overstory_root_pp%copy  (idx_beg, idx_end, vsfm_por, vsfm_perm, vsfm_alpha, &
       vsfm_lambda, vsfm_relperm_type, vsfm_relperm_param_1, vsfm_relperm_param_2, &
       vsfm_residual_sat, vsfm_satfunc_type)

  ! Add xylem grid cells
  idx_beg = ncells + 1
  idx_end = ncells + soil_nx*soil_ny*overstory_xylem_nz
  ncells  = idx_end

  call overstory_xylem_mesh%copy(idx_beg, idx_end, xc, yc, zc, area, vol, filter)

  call overstory_xylem_pp%copy  (idx_beg, idx_end, vsfm_por, vsfm_perm, vsfm_alpha, &
       vsfm_lambda, vsfm_relperm_type, vsfm_relperm_param_1, vsfm_relperm_param_2, &
       vsfm_residual_sat, vsfm_satfunc_type)

  ! Add leaf grid cells
  idx_beg = ncells + 1
  idx_end = ncells + soil_nx*soil_ny*overstory_leaf_nz
  ncells  = idx_end

  call overstory_leaf_mesh%copy(idx_beg, idx_end, xc, yc, zc, area, vol, filter)

  call overstory_leaf_pp%copy  (idx_beg, idx_end, vsfm_por, vsfm_perm, vsfm_alpha, &
       vsfm_lambda, vsfm_relperm_type, vsfm_relperm_param_1, vsfm_relperm_param_2, &
       vsfm_residual_sat, vsfm_satfunc_type)

  !
  ! Understory
  !

  ! Add root grid cells

  idx_beg = ncells + 1
  idx_end = ncells + soil_nx*soil_ny*understory_root_nz
  ncells  = idx_end

  call understory_root_mesh%copy(idx_beg, idx_end, xc, yc, zc, area, vol, filter)

  call understory_root_pp%copy  (idx_beg, idx_end, vsfm_por, vsfm_perm, vsfm_alpha, &
       vsfm_lambda, vsfm_relperm_type, vsfm_relperm_param_1, vsfm_relperm_param_2, &
       vsfm_residual_sat, vsfm_satfunc_type)


  ! Add xylem grid cells
  idx_beg = ncells + 1
  idx_end = ncells + soil_nx*soil_ny*understory_xylem_nz
  ncells  = idx_end

  call understory_xylem_mesh%copy(idx_beg, idx_end, xc, yc, zc, area, vol, filter)

  call understory_xylem_pp%copy  (idx_beg, idx_end, vsfm_por, vsfm_perm, vsfm_alpha, &
       vsfm_lambda, vsfm_relperm_type, vsfm_relperm_param_1, vsfm_relperm_param_2, &
       vsfm_residual_sat, vsfm_satfunc_type)

  ! Add leaf grid cells
  idx_beg = ncells + 1
  idx_end = ncells + soil_nx*soil_ny*understory_leaf_nz
  ncells  = idx_end

  call understory_leaf_mesh%copy(idx_beg, idx_end, xc, yc, zc, area, vol, filter)

  call understory_leaf_pp%copy  (idx_beg, idx_end, vsfm_por, vsfm_perm, vsfm_alpha, &
       vsfm_lambda, vsfm_relperm_type, vsfm_relperm_param_1, vsfm_relperm_param_2, &
       vsfm_residual_sat, vsfm_satfunc_type)

  !
  ! Set up the meshes
  !    
  call vsfm_mpp%SetNumMeshes(1)

  ! Soil Mesh
  imesh             = 1
  soil_ncells_ghost = 0

  call vsfm_mpp%MeshSetName                (imesh, 'Soil mesh'                             )
  call vsfm_mpp%MeshSetOrientation         (imesh, MESH_ALONG_GRAVITY                      )
  call vsfm_mpp%MeshSetID                  (imesh, MESH_CLM_SOIL_COL                       )
  call vsfm_mpp%MeshSetDimensions          (imesh, ncells, soil_ncells_ghost, soil_nz )

  call vsfm_mpp%MeshSetGridCellFilter      (imesh, filter                             )
  call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC     , xc                    )
  call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC     , yc                    )
  call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC     , zc                    )
  call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX     , dx               )
  call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY     , dy               )
  call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ     , dz               )
  call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA   , area                  )
  call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_VOLUME , vol                   )


  !
  ! Add connections
  !

  overstory_nconn   = (2*overstory_root_nz    +   & ! root-to-root and root-to-soil
                       overstory_xylem_nz - 1 +   & ! xylem-to-xylem
                       overstory_leaf_nz      ) * & ! xylem-to-leaf
                       soil_nx * soil_ny
  
  understory_nconn  = (2*understory_root_nz    +   & ! root-to-root and root-to-soil
                       understory_xylem_nz - 1 +   & ! xylem-to-xylem
                       understory_leaf_nz      ) * & ! xylem-to-leaf
                       soil_nx * soil_ny

  nconn = soil_mesh_nconn + overstory_nconn + understory_nconn
  
  allocate (conn_id_up               (nconn))
  allocate (conn_id_dn               (nconn))
  allocate (conn_dist_up             (nconn))
  allocate (conn_dist_dn             (nconn))
  allocate (conn_area                (nconn))
  allocate (conn_type                (nconn))
  allocate (conn_flux_type           (nconn)); conn_flux_type(:) = 0
  allocate (conn_cond_type           (nconn)); conn_cond_type(:) = 0
  allocate (conn_cond                (nconn)); conn_cond     (:) = 0.d0
  allocate (conn_cond_up             (nconn)); conn_cond_up  (:) = 0.d0
  allocate (conn_cond_dn             (nconn)); conn_cond_dn  (:) = 0.d0

  allocate (conn_satparam_up_itype   (nconn)); conn_satparam_up_itype   (:) = 0
  allocate (conn_satparam_up_param_1 (nconn)); conn_satparam_up_param_1 (:) = 0.d0
  allocate (conn_satparam_up_param_2 (nconn)); conn_satparam_up_param_2 (:) = 0.d0
  allocate (conn_satparam_up_param_3 (nconn)); conn_satparam_up_param_3 (:) = 0.d0
  allocate (conn_relperm_up_itype    (nconn)); conn_relperm_up_itype    (:) = 0
  allocate (conn_relperm_up_param_1  (nconn)); conn_relperm_up_param_1  (:) = 0.d0
  allocate (conn_relperm_up_param_2  (nconn)); conn_relperm_up_param_2  (:) = 0.d0

  allocate (conn_satparam_dn_itype   (nconn)); conn_satparam_dn_itype   (:) = 0
  allocate (conn_satparam_dn_param_1 (nconn)); conn_satparam_dn_param_1 (:) = 0.d0
  allocate (conn_satparam_dn_param_2 (nconn)); conn_satparam_dn_param_2 (:) = 0.d0
  allocate (conn_satparam_dn_param_3 (nconn)); conn_satparam_dn_param_3 (:) = 0.d0
  allocate (conn_relperm_dn_itype    (nconn)); conn_relperm_dn_itype    (:) = 0
  allocate (conn_relperm_dn_param_1  (nconn)); conn_relperm_dn_param_1  (:) = 0.d0
  allocate (conn_relperm_dn_param_2  (nconn)); conn_relperm_dn_param_2  (:) = 0.d0
  allocate (conn_relperm_dn_param_3  (nconn)); conn_relperm_dn_param_3  (:) = 0.d0
  
  iconn = 0

  !  
  ! Soil-to-Soil Connections: Add vertical connections
  !
  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, soil_nz-1
           if (soil_id(ii,jj,kk) > 0 .and. soil_id(ii,jj,kk+1) > 0) then
              iconn                 = iconn + 1
              conn_id_up(iconn)     = soil_id(ii,jj,kk  )
              conn_id_dn(iconn)     = soil_id(ii,jj,kk+1)
              conn_dist_up(iconn)   = 0.5d0*soil_dz
              conn_dist_dn(iconn)   = 0.5d0*soil_dz
              conn_area(iconn)      = soil_dx*soil_dy
              conn_type(iconn)      = CONN_VERTICAL
              conn_flux_type(iconn) = DARCY_FLUX_TYPE
           end if
        end do
     end do
  end do
  
  !  
  ! Soil-to-Soil Connections: Add horizontal connections in x-direction
  !
  if (.not.is_soil_horizontally_disconnected) then
     do ii = 1, soil_nx-1
        do jj = 1, soil_ny
           do kk = 1, soil_nz
              if (soil_id(ii,jj,kk) > 0 .and. soil_id(ii+1,jj,kk) > 0) then
                 iconn                 = iconn + 1
                 conn_id_up(iconn)     = soil_id(ii  ,jj,kk)
                 conn_id_dn(iconn)     = soil_id(ii+1,jj,kk)
                 conn_dist_up(iconn)   = 0.5d0*soil_dx
                 conn_dist_dn(iconn)   = 0.5d0*soil_dx
                 conn_area(iconn)      = soil_dy*soil_dz
                 conn_type(iconn)      = CONN_HORIZONTAL
                 conn_flux_type(iconn) = DARCY_FLUX_TYPE
              end if
           end do
        end do
     end do

     if (soil_ny > 1) then
        write(*,*)'Add code to for horizontal connections in y-direciton'
        stop
     endif
  endif

  !  
  ! Overstory Connections
  !

  do ii = 1, soil_nx
     do jj = 1, soil_ny

        ! Root-to-Soil: Condunctance flux type
        do kk = 1, overstory_root_nz
           iconn                 = iconn + 1
           conn_id_up(iconn)     = overstory_root_mesh%id(ii,jj,kk)
           conn_id_dn(iconn)     = soil_id          (ii,jj,kk-1+top_active_layer_kk_index(ii,jj))

           conn_dist_up(iconn)   = 0.d0
           conn_dist_dn(iconn)   = overstory_root_length_profile(kk)
           conn_area(iconn)      = overstory_root_area_profile  (kk)
           conn_type(iconn)      = CONN_HORIZONTAL
           conn_flux_type(iconn) = CONDUCTANCE_FLUX_TYPE
           conn_cond_type(iconn) = CONDUCTANCE_MANOLI_TYPE
           conn_cond_up(iconn)   = overstory_root_conductance
           conn_cond_dn(iconn)   = perm_z_top /vish2o * (denh2o * grav) / & ! [m/s]
                                   overstory_root_length_profile(kk)        ! [m]

           ! Saturation up: CHUANG
           conn_satparam_up_itype   (iconn) = vsfm_satfunc_type    (conn_id_up(iconn) )
           conn_satparam_up_param_1 (iconn) = vsfm_alpha           (conn_id_up(iconn) )
           conn_satparam_up_param_2 (iconn) = vsfm_lambda          (conn_id_up(iconn) )
           conn_satparam_up_param_3 (iconn) = vsfm_residual_sat    (conn_id_up(iconn) )
           ! Relative Perm. up: WEIBULL
           conn_relperm_up_itype    (iconn) = vsfm_relperm_type    (conn_id_up(iconn) )
           conn_relperm_up_param_1  (iconn) = vsfm_relperm_param_1 (conn_id_up(iconn) )
           conn_relperm_up_param_2  (iconn) = vsfm_relperm_param_2 (conn_id_up(iconn) )

           ! Saturation dn: CHUANG
           conn_satparam_dn_itype   (iconn) = vsfm_satfunc_type    (conn_id_dn(iconn) )
           conn_satparam_dn_param_1 (iconn) = vsfm_alpha           (conn_id_dn(iconn) )
           conn_satparam_dn_param_2 (iconn) = vsfm_lambda          (conn_id_dn(iconn) )
           conn_satparam_dn_param_3 (iconn) = vsfm_residual_sat    (conn_id_dn(iconn) )
           ! Relative Perm. dn: MUALEM
           conn_relperm_dn_itype    (iconn) = vsfm_relperm_type    (conn_id_dn(iconn) )
           conn_relperm_dn_param_1  (iconn) = vsfm_alpha           (conn_id_dn(iconn) )
           conn_relperm_dn_param_2  (iconn) = vsfm_lambda          (conn_id_dn(iconn) )
           conn_relperm_dn_param_3  (iconn) = vsfm_residual_sat    (conn_id_dn(iconn) )

        end do

        ! Xylem-to-Root: Conductance flux type
        do kk = 1, overstory_root_nz
           iconn                 = iconn + 1
           conn_id_up(iconn)     = overstory_xylem_mesh%id(ii,jj,1)
           conn_id_dn(iconn)     = overstory_root_mesh%id (ii,jj,kk)
           conn_dist_up(iconn)   = 0.1d0
           conn_dist_dn(iconn)   = 0.1d0
           conn_area(iconn)      = overstory_area_sapwood
           conn_type(iconn)      = CONN_VERTICAL
           conn_flux_type(iconn) = CONDUCTANCE_FLUX_TYPE
           conn_cond_type(iconn) = CONDUCTANCE_CAMPBELL_TYPE
           conn_cond(iconn)      = overstory_root_conductance
           conn_cond_up(iconn)   = overstory_root_conductance
           conn_cond_dn(iconn)   = overstory_root_conductance

           ! Saturation up: CHUANG
           conn_satparam_up_itype   (iconn) = vsfm_satfunc_type    (conn_id_up(iconn) )
           conn_satparam_up_param_1 (iconn) = vsfm_alpha           (conn_id_up(iconn) )
           conn_satparam_up_param_2 (iconn) = vsfm_lambda          (conn_id_up(iconn) )
           conn_satparam_up_param_3 (iconn) = vsfm_residual_sat    (conn_id_up(iconn) )
           ! Relative Perm. up: WEIBULL
           conn_relperm_up_itype    (iconn) = vsfm_relperm_type    (conn_id_up(iconn) )
           conn_relperm_up_param_1  (iconn) = vsfm_relperm_param_1 (conn_id_up(iconn) )
           conn_relperm_up_param_2  (iconn) = vsfm_relperm_param_2 (conn_id_up(iconn) )

           ! Saturation dn: CHUANG
           conn_satparam_dn_itype   (iconn) = vsfm_satfunc_type    (conn_id_dn(iconn) )
           conn_satparam_dn_param_1 (iconn) = vsfm_alpha           (conn_id_dn(iconn) )
           conn_satparam_dn_param_2 (iconn) = vsfm_lambda          (conn_id_dn(iconn) )
           conn_satparam_dn_param_3 (iconn) = vsfm_residual_sat    (conn_id_dn(iconn) )
           ! Relative Perm. dn: WEIBULL
           conn_relperm_dn_itype    (iconn) = vsfm_relperm_type    (conn_id_dn(iconn) )
           conn_relperm_dn_param_1  (iconn) = vsfm_relperm_param_1 (conn_id_dn(iconn) )
           conn_relperm_dn_param_2  (iconn) = vsfm_relperm_param_2 (conn_id_dn(iconn) )
        end do

        ! Xylem-to-Xylem: Darcy flux
        do kk = 1, overstory_xylem_nz-1
           iconn                 = iconn + 1
           conn_id_up(iconn)     = overstory_xylem_mesh%id(ii,jj,kk  )
           conn_id_dn(iconn)     = overstory_xylem_mesh%id(ii,jj,kk+1)
           conn_dist_up(iconn)   = 0.5d0*soil_dz
           conn_dist_dn(iconn)   = 0.5d0*soil_dz
           conn_area(iconn)      = overstory_area_sapwood
           conn_type(iconn)      = CONN_VERTICAL
           conn_flux_type(iconn) = DARCY_FLUX_TYPE
        end do

        ! Xylem-to-Leaf
        do kk = 1, overstory_leaf_nz
           idx                   = overstory_branch_2_xylem_index(kk)
           
           iconn                 = iconn + 1
           conn_id_up(iconn)     = overstory_xylem_mesh%id(ii,jj,idx)
           conn_id_dn(iconn)     = overstory_leaf_mesh%id (ii,jj,kk )
           conn_dist_up(iconn)   = 0.5d0*overstory_branch_length_profile(idx)
           conn_dist_dn(iconn)   = 0.5d0*overstory_branch_length_profile(idx)
           conn_area(iconn)      = overstory_xylem_area_profile(idx)*overstory_branch_area_ratio
           conn_type(iconn)      = CONN_HORIZONTAL
           conn_flux_type(iconn) = DARCY_FLUX_TYPE
        end do

     end do
  end do

  !  
  ! Understory Connections
  !

  do ii = 1, soil_nx
     do jj = 1, soil_ny

        ! Root-to-Soil: Condunctance flux type
        do kk = 1, understory_root_nz
           iconn                 = iconn + 1
           conn_id_up(iconn)     = understory_root_mesh%id(ii,jj,kk)
           conn_id_dn(iconn)     = soil_id           (ii,jj,kk-1+top_active_layer_kk_index(ii,jj))
           conn_dist_up(iconn)   = 0.d0
           conn_dist_dn(iconn)   = understory_root_length_profile(kk)
           conn_area(iconn)      = understory_root_area_profile  (kk)
           conn_type(iconn)      = CONN_HORIZONTAL
           conn_flux_type(iconn) = CONDUCTANCE_FLUX_TYPE
           conn_cond_type(iconn) = CONDUCTANCE_MANOLI_TYPE
           conn_cond_up(iconn)   = understory_root_conductance
           conn_cond_dn(iconn)   = perm_z_top /vish2o * (denh2o * grav) / & ! [m/s]
                understory_root_length_profile(kk)       ! [m]

           ! Saturation up: CHUANG
           conn_satparam_up_itype   (iconn) = vsfm_satfunc_type    (conn_id_up(iconn) )
           conn_satparam_up_param_1 (iconn) = vsfm_alpha           (conn_id_up(iconn) )
           conn_satparam_up_param_2 (iconn) = vsfm_lambda          (conn_id_up(iconn) )
           conn_satparam_up_param_3 (iconn) = vsfm_residual_sat    (conn_id_up(iconn) )
           ! Relative Perm. up: MUALEM
           conn_relperm_up_itype    (iconn) = vsfm_relperm_type    (conn_id_up(iconn) )
           conn_relperm_up_param_1  (iconn) = vsfm_relperm_param_1 (conn_id_up(iconn) )
           conn_relperm_up_param_2  (iconn) = vsfm_relperm_param_2 (conn_id_up(iconn) )

           ! Saturation dn: CHUANG
           conn_satparam_dn_itype   (iconn) = vsfm_satfunc_type    (conn_id_dn(iconn) )
           conn_satparam_dn_param_1 (iconn) = vsfm_alpha           (conn_id_dn(iconn) )
           conn_satparam_dn_param_2 (iconn) = vsfm_lambda          (conn_id_dn(iconn) )
           conn_satparam_dn_param_3 (iconn) = vsfm_residual_sat    (conn_id_dn(iconn) )
           ! Relative Perm. dn: MUALEM
           conn_relperm_dn_itype    (iconn) = vsfm_relperm_type    (conn_id_dn(iconn) )
           conn_relperm_dn_param_1  (iconn) = vsfm_alpha           (conn_id_dn(iconn) )
           conn_relperm_dn_param_2  (iconn) = vsfm_lambda          (conn_id_dn(iconn) )
           conn_relperm_dn_param_3  (iconn) = vsfm_residual_sat    (conn_id_dn(iconn) )
        end do

        ! Xylem-to-Root: Conductance flux type
        do kk                    = 1, understory_root_nz
           iconn                 = iconn + 1
           conn_id_up(iconn)     = understory_xylem_mesh%id(ii,jj,1)
           conn_id_dn(iconn)     = understory_root_mesh%id (ii,jj,kk)
           conn_dist_up(iconn)   = 0.1d0
           conn_dist_dn(iconn)   = 0.1d0
           conn_area(iconn)      = understory_area_sapwood
           conn_type(iconn)      = CONN_VERTICAL
           conn_flux_type(iconn) = CONDUCTANCE_FLUX_TYPE
           conn_cond_type(iconn) = CONDUCTANCE_CAMPBELL_TYPE
           conn_cond_up          = understory_root_conductance
           conn_cond_up(iconn)   = understory_root_conductance
           conn_cond_dn(iconn)   = understory_root_conductance

           ! Saturation up: CHUANG
           conn_satparam_up_itype   (iconn) = vsfm_satfunc_type    (conn_id_up(iconn) )
           conn_satparam_up_param_1 (iconn) = vsfm_alpha           (conn_id_up(iconn) )
           conn_satparam_up_param_2 (iconn) = vsfm_lambda          (conn_id_up(iconn) )
           conn_satparam_up_param_3 (iconn) = vsfm_residual_sat    (conn_id_up(iconn) )
           ! Relative Perm. up: WEIBULL
           conn_relperm_up_itype    (iconn) = vsfm_relperm_type    (conn_id_up(iconn) )
           conn_relperm_up_param_1  (iconn) = vsfm_relperm_param_1 (conn_id_up(iconn) )
           conn_relperm_up_param_2  (iconn) = vsfm_relperm_param_2 (conn_id_up(iconn) )

           ! Saturation dn: CHUANG
           conn_satparam_dn_itype   (iconn) = vsfm_satfunc_type    (conn_id_dn(iconn) )
           conn_satparam_dn_param_1 (iconn) = vsfm_alpha           (conn_id_dn(iconn) )
           conn_satparam_dn_param_2 (iconn) = vsfm_lambda          (conn_id_dn(iconn) )
           conn_satparam_dn_param_3 (iconn) = vsfm_residual_sat    (conn_id_dn(iconn) )
           ! Relative Perm. dn: WEIBULL
           conn_relperm_dn_itype    (iconn) = vsfm_relperm_type    (conn_id_dn(iconn) )
           conn_relperm_dn_param_1  (iconn) = vsfm_relperm_param_1 (conn_id_dn(iconn) )
           conn_relperm_dn_param_2  (iconn) = vsfm_relperm_param_2 (conn_id_dn(iconn) )
        end do

        ! Xylem-to-Xylem: Darcy flux
        do kk                    = 1, understory_xylem_nz-1
           iconn                 = iconn + 1
           conn_id_up(iconn)     = understory_xylem_mesh%id(ii,jj,kk  )
           conn_id_dn(iconn)     = understory_xylem_mesh%id(ii,jj,kk+1)
           conn_dist_up(iconn)   = 0.5d0*soil_dz
           conn_dist_dn(iconn)   = 0.5d0*soil_dz
           conn_area(iconn)      = understory_area_sapwood
           conn_type(iconn)      = CONN_VERTICAL
           conn_flux_type(iconn) = DARCY_FLUX_TYPE
        end do

        ! Xylem-to-Leaf
        do kk = 1, understory_leaf_nz
           idx                   = understory_branch_2_xylem_index(kk)
           
           iconn                 = iconn + 1
           conn_id_up(iconn)     = understory_xylem_mesh%id(ii,jj,idx)
           conn_id_dn(iconn)     = understory_leaf_mesh%id (ii,jj,kk )
           conn_dist_up(iconn)   = 0.5d0*understory_branch_length_profile(idx)
           conn_dist_dn(iconn)   = 0.5d0*understory_branch_length_profile(idx)
           conn_area(iconn)      = understory_xylem_area_profile(idx)*understory_branch_area_ratio
           conn_type(iconn)      = CONN_HORIZONTAL
           conn_flux_type(iconn) = DARCY_FLUX_TYPE
        end do

     end do
  end do

  imesh = 1
  call vsfm_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
       iconn,  conn_id_up, conn_id_dn,          &
       conn_dist_up, conn_dist_dn,  &
       conn_area,  conn_type)

  deallocate (conn_id_up   )
  deallocate (conn_id_dn   )
  deallocate (conn_dist_up )
  deallocate (conn_dist_dn )
  deallocate (conn_area    )
  deallocate (conn_type    )

end subroutine add_single_mesh

!------------------------------------------------------------------------
subroutine setup_soil_mesh()
  !
  use problem_parameters
  !
  implicit none
  !
#include <petsc/finclude/petsc.h>
  !
  PetscInt            :: count
  PetscInt            :: ii
  PetscInt            :: jj
  PetscInt            :: kk
  
  allocate(soil_id                   (soil_nx,soil_ny,soil_nz ))
  allocate(soil_xc3d                 (soil_nx,soil_ny,soil_nz ))
  allocate(soil_yc3d                 (soil_nx,soil_ny,soil_nz ))
  allocate(soil_zc3d                 (soil_nx,soil_ny,soil_nz ))
  allocate(top_active_layer_kk_index (soil_nx,soil_ny         )); top_active_layer_kk_index(:,:) = 0
  allocate(elevation                 (soil_nx,soil_ny         ))

  soil_id(:,:,:) = 0

  count = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, soil_nz
           soil_xc3d(ii,jj,kk) =  soil_dx/2 + soil_dx*(ii-1)
           soil_yc3d(ii,jj,kk) =  soil_dy/2 + soil_dy*(jj-1)
           soil_zc3d(ii,jj,kk) = -soil_dz/2 - soil_dz*(kk-1)
           if (soil_zc3d(ii,jj,kk) <= soil_zc3d(1,1,1) -slope * soil_dx*(ii-1) ) then
              count = count + 1
              soil_id(ii,jj,kk) = count
           endif
        end do
     end do
  end do

  soil_ncells       = count

  allocate (soil_xc      (soil_ncells))
  allocate (soil_yc      (soil_ncells))
  allocate (soil_zc      (soil_ncells))
  allocate (soil_soil_dx (soil_ncells))
  allocate (soil_soil_dy (soil_ncells))
  allocate (soil_soil_dz (soil_ncells))
  allocate (soil_area    (soil_ncells))
  allocate (soil_filter  (soil_ncells))
  allocate (soil_vol     (soil_ncells))

  count           = 0
  soil_mesh_nconn = 0

  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, soil_nz

           ! is the grid cell active?
           if (soil_id(ii,jj,kk) > 0) then

              count = count + 1

              soil_xc(count)      = soil_xc3d(ii,jj,kk)
              soil_yc(count)      = soil_yc3d(ii,jj,kk)
              soil_zc(count)      = soil_zc3d(ii,jj,kk)
              soil_soil_dx(count) = soil_dx
              soil_soil_dy(count) = soil_dy
              soil_soil_dz(count) = soil_dz
              soil_area(count)    = soil_dx*soil_dy
              soil_vol(count)     = soil_dx*soil_dy*soil_dz
              soil_filter(count)  = 1

              ! add vertical connection except if the grid cell is last in the z-dir
              if (kk < soil_nz) soil_mesh_nconn = soil_mesh_nconn + 1

              ! should horizontal connections be added?
              if (.not.is_soil_horizontally_disconnected) then
                 ! add horizontal connection except if the grid cell is last in x-dir
                 if (ii < soil_nx .and. soil_id(ii+1,jj,kk) > 0) soil_mesh_nconn = soil_mesh_nconn + 1
              end if

              ! is the index of first active grid cell in z-dir already set?
              if (top_active_layer_kk_index(ii,jj) == 0) then
                 top_active_layer_kk_index(ii,jj) = kk
                 elevation                (ii,jj) = soil_zc3d(ii,jj,kk) + soil_dz/2.d0
              end if

           end if
        end do
     end do
  end do

end subroutine setup_soil_mesh

!------------------------------------------------------------------------
subroutine setup_overstory_mesh()
  !
#include <petsc/finclude/petsc.h>
  !
  use mpp_varcon         , only : grav
  use mpp_varcon         , only : denh2o
  use SaturationFunction , only : RELPERM_FUNC_WEIBULL
  use SaturationFunction , only : SAT_FUNC_CHUANG
  use petscsys
  use problem_parameters
  !
  implicit none
  !
  PetscInt  :: ii,jj,kk
  PetscInt  :: idx
  PetscInt  :: id_value
  PetscInt  :: count
  PetscReal :: zz
  PetscReal :: soil_volume
  PetscReal :: root_length

  overstory_leaf_nz  = 0
  overstory_xylem_nz = int(overstory_canopy_height/soil_dz)
  overstory_root_nz  = int(overstory_root_depth   /soil_dz)

  allocate (overstory_xylem_area_profile    (overstory_xylem_nz))
  allocate (overstory_xylem_has_branches    (overstory_xylem_nz))
  allocate (overstory_branch_length_profile (overstory_xylem_nz))
  allocate (overstory_branch_2_xylem_index  (overstory_xylem_nz))
  allocate (overstory_root_length_profile   (overstory_root_nz ))
  allocate (overstory_root_area_profile     (overstory_root_nz ))
  allocate (overstory_root_vol_profile      (overstory_root_nz ))
  
  overstory_branch_length_profile(:) = 0.d0
  count = 0
  
  do kk = 1, overstory_xylem_nz
     zz = (kk-1)*soil_dz + soil_dz/2.d0
     overstory_xylem_area_profile(kk) = overstory_area_sapwood * &
          ( 1.d0 - overstory_area_tappering_coeff * zz/overstory_canopy_height) ** 2.d0

     if (overstory_lad_profile(kk) > 0.d0) then
        count                               = count + 1
        overstory_leaf_nz                   = overstory_leaf_nz + 1;
        overstory_xylem_has_branches(kk)    = PETSC_TRUE
        overstory_branch_length_profile(kk) = overstory_xylem_area_profile(kk) * &
                                              overstory_branch_area_ratio
        overstory_branch_2_xylem_index(count) = kk
     else
        overstory_xylem_has_branches(kk) = PETSC_FALSE
     end if
  end do

  soil_volume = soil_dx * soil_dy * soil_dz
  do kk = 1, overstory_root_nz
     root_length                       = overstory_B_profile(kk) * soil_volume
     overstory_root_length_profile(kk) = root_length
     overstory_root_area_profile(kk)   = 2*PI*overstory_root_radius*root_length
     overstory_root_vol_profile(kk)    = PI*(overstory_root_radius**2.d0)*root_length
  end do

  allocate(overstory_root_pp%por              (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%perm             (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%lambda           (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%alpha            (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%eff_porosity     (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%residual_sat     (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%satfunc_type     (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%relperm_type     (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%relperm_param_1  (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_pp%relperm_param_2  (soil_nx*soil_ny*overstory_root_nz))

  allocate(overstory_root_mesh%xc               (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_mesh%yc               (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_mesh%zc               (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_mesh%area             (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_mesh%vol              (soil_nx*soil_ny*overstory_root_nz))
  allocate(overstory_root_mesh%filter           (soil_nx*soil_ny*overstory_root_nz))

  allocate(overstory_xylem_pp%por             (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%perm            (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%lambda          (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%alpha           (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%eff_porosity    (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%residual_sat    (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%satfunc_type    (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%relperm_type    (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%relperm_param_1 (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_pp%relperm_param_2 (soil_nx*soil_ny*overstory_xylem_nz))

  allocate(overstory_xylem_mesh%xc              (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_mesh%yc              (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_mesh%zc              (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_mesh%area            (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_mesh%vol             (soil_nx*soil_ny*overstory_xylem_nz))
  allocate(overstory_xylem_mesh%filter          (soil_nx*soil_ny*overstory_xylem_nz))

  allocate(overstory_leaf_pp%por              (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%perm             (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%lambda           (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%alpha            (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%eff_porosity     (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%residual_sat     (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%satfunc_type     (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%relperm_type     (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%relperm_param_1  (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_pp%relperm_param_2  (soil_nx*soil_ny*overstory_leaf_nz))

  allocate(overstory_leaf_mesh%xc               (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_mesh%yc               (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_mesh%zc               (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_mesh%area             (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_mesh%vol              (soil_nx*soil_ny*overstory_leaf_nz))
  allocate(overstory_leaf_mesh%filter           (soil_nx*soil_ny*overstory_leaf_nz))

  allocate(overstory_root_mesh%id               (soil_nx, soil_ny, overstory_root_nz ))
  allocate(overstory_xylem_mesh%id              (soil_nx, soil_ny, overstory_xylem_nz))
  allocate(overstory_leaf_mesh%id               (soil_nx, soil_ny, overstory_leaf_nz ))

  if (multi_goveqns_formulation) then
     id_value = 0
  else
     id_value = soil_ncells
  end if
  
  ! Add root grid cells
  count = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, overstory_root_nz

           id_value                               = id_value + 1
           overstory_root_mesh%id           (ii,jj,kk) = id_value

           count                                  = count + 1

           overstory_root_mesh%xc              (count) = soil_xc3d(ii,jj,1)
           overstory_root_mesh%yc              (count) = soil_yc3d(ii,jj,1)
           overstory_root_mesh%zc              (count) = elevation(ii,jj) - soil_dz/2.d0 - soil_dz*(kk-1)
           overstory_root_mesh%area            (count) = overstory_root_area_profile(kk)
           overstory_root_mesh%vol             (count) = overstory_root_vol_profile(kk)
           overstory_root_mesh%filter          (count) = 1

           overstory_root_pp%por             (count) = 0.d0
           overstory_root_pp%perm            (count) = 0.d0
           overstory_root_pp%alpha           (count) = overstory_xylem_phi0
           overstory_root_pp%lambda          (count) = overstory_xylem_p
           overstory_root_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           overstory_root_pp%relperm_param_1 (count) = overstory_xylem_vulnerability_d * grav * denh2o
           overstory_root_pp%relperm_param_2 (count) = overstory_xylem_vulnerability_c
           overstory_root_pp%residual_sat    (count) = 0.d0
           overstory_root_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  ! Add xylem grid cells
  count = 0
  if (multi_goveqns_formulation) id_value = 0

  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, overstory_xylem_nz
           id_value                                = id_value + 1
           overstory_xylem_mesh%id           (ii,jj,kk) = id_value

           count                                   = count + 1
           overstory_xylem_mesh%xc              (count) = soil_xc3d(ii,jj,1)
           overstory_xylem_mesh%yc              (count) = soil_yc3d(ii,jj,1)
           overstory_xylem_mesh%zc              (count) = elevation(ii,jj) + soil_dz/2.d0 + (kk-1)*soil_dz
           overstory_xylem_mesh%area            (count) = overstory_xylem_area_profile(kk)
           overstory_xylem_mesh%vol             (count) = overstory_xylem_area_profile(kk) * soil_dz
           overstory_xylem_mesh%filter          (count) = 1

           overstory_xylem_pp%por             (count) = overstory_xylem_porosity
           overstory_xylem_pp%perm            (count) = overstory_xylem_Kmax * vish2o / (denh2o * grav)
           overstory_xylem_pp%alpha           (count) = overstory_xylem_phi0
           overstory_xylem_pp%lambda          (count) = overstory_xylem_p
           overstory_xylem_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           overstory_xylem_pp%relperm_param_1 (count) = overstory_xylem_vulnerability_d * grav * denh2o
           overstory_xylem_pp%relperm_param_2 (count) = overstory_xylem_vulnerability_c
           overstory_xylem_pp%residual_sat    (count) = 0.d0
           overstory_xylem_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  ! Add leaf grid cells
  count = 0
  if (multi_goveqns_formulation) id_value = 0

  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, overstory_leaf_nz
           id_value                               = id_value + 1
           overstory_leaf_mesh%id           (ii,jj,kk) = id_value

           count                                  = count + 1
           idx                                    = overstory_branch_2_xylem_index(kk)

           overstory_leaf_mesh%xc              (count) = soil_xc3d(ii,jj,1)-overstory_branch_length_profile(idx)
           overstory_leaf_mesh%yc              (count) = soil_yc3d(ii,jj,1)
           overstory_leaf_mesh%zc              (count) = elevation(ii,jj) + soil_dz/2.d0 + (kk-1)*soil_dz + (overstory_xylem_nz-overstory_leaf_nz)*soil_dz
           overstory_leaf_mesh%area            (count) = overstory_xylem_area_profile(kk) * overstory_branch_area_ratio

           overstory_leaf_mesh%vol             (count) = overstory_leaf_mesh%area(count) * overstory_branch_length_profile(idx)
           overstory_leaf_mesh%filter          (count) = 1

           overstory_leaf_pp%por             (count) = 0.d0
           overstory_leaf_pp%perm            (count) = overstory_xylem_Kmax * vish2o / (denh2o * grav)
           overstory_leaf_pp%alpha           (count) = overstory_xylem_phi0
           overstory_leaf_pp%lambda          (count) = overstory_xylem_p
           overstory_leaf_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           overstory_leaf_pp%relperm_param_1 (count) = overstory_xylem_vulnerability_d * grav * denh2o
           overstory_leaf_pp%relperm_param_2 (count) = overstory_xylem_vulnerability_c
           overstory_leaf_pp%residual_sat    (count) = 0.d0
           overstory_leaf_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

end subroutine setup_overstory_mesh

!------------------------------------------------------------------------
subroutine setup_understory_mesh()
  !
#include <petsc/finclude/petsc.h>
  !
  use mpp_varcon         , only : grav
  use mpp_varcon         , only : denh2o
  use SaturationFunction , only : RELPERM_FUNC_WEIBULL
  use SaturationFunction , only : SAT_FUNC_CHUANG
  use petscsys
  use problem_parameters
  !
  implicit none
  !
  PetscInt  :: ii,jj,kk
  PetscInt  :: idx
  PetscInt  :: id_value
  PetscInt  :: count
  PetscReal :: zz
  PetscReal :: soil_volume
  PetscReal :: root_length

  understory_leaf_nz  = 0
  understory_xylem_nz = int(understory_canopy_height/soil_dz)
  understory_root_nz  = int(understory_root_depth   /soil_dz)

  allocate (understory_xylem_area_profile    (understory_xylem_nz))
  allocate (understory_xylem_has_branches    (understory_xylem_nz))
  allocate (understory_branch_length_profile (understory_xylem_nz))
  allocate (understory_branch_2_xylem_index  (understory_xylem_nz))
  allocate (understory_root_length_profile   (understory_root_nz ))
  allocate (understory_root_area_profile     (understory_root_nz ))
  allocate (understory_root_vol_profile      (understory_root_nz ))

  understory_branch_length_profile(:) = 0.d0
  count = 0

  do kk = 1, understory_xylem_nz
     zz = (kk-1)*soil_dz + soil_dz/2.d0
     understory_xylem_area_profile(kk) = understory_area_sapwood * &
          ( 1.d0 - understory_area_tappering_coeff * zz/understory_canopy_height) ** 2.d0

     if (understory_lad_profile(kk) > 0.d0) then
        count                                = count + 1
        understory_leaf_nz                   = understory_leaf_nz + 1;
        understory_xylem_has_branches(kk)    = PETSC_TRUE
        understory_branch_length_profile(kk) = understory_xylem_area_profile(kk) * &
                                               understory_branch_area_ratio
        understory_branch_2_xylem_index(count) = kk
     else
        understory_xylem_has_branches(kk) = PETSC_FALSE
     end if
     
  end do

  soil_volume = soil_dx * soil_dy * soil_dz
  do kk = 1, understory_root_nz
     root_length                        = understory_B_profile(kk) * soil_volume
     understory_root_length_profile(kk) = root_length
     understory_root_area_profile(kk)   = 2*PI*understory_root_radius*root_length
     understory_root_vol_profile(kk)    = PI*(overstory_root_radius**2.d0)*root_length
  end do

  allocate(understory_root_pp%por              (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%perm             (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%lambda           (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%alpha            (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%eff_porosity     (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%residual_sat     (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%satfunc_type     (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%relperm_type     (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%relperm_param_1  (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_pp%relperm_param_2  (soil_nx*soil_ny*understory_root_nz))

  allocate(understory_root_mesh%xc               (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_mesh%yc               (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_mesh%zc               (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_mesh%area             (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_mesh%vol              (soil_nx*soil_ny*understory_root_nz))
  allocate(understory_root_mesh%filter           (soil_nx*soil_ny*understory_root_nz))

  allocate(understory_xylem_pp%por             (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%perm            (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%lambda          (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%alpha           (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%eff_porosity    (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%residual_sat    (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%satfunc_type    (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%relperm_type    (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%relperm_param_1 (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_pp%relperm_param_2 (soil_nx*soil_ny*understory_xylem_nz))

  allocate(understory_xylem_mesh%xc              (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_mesh%yc              (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_mesh%zc              (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_mesh%area            (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_mesh%vol             (soil_nx*soil_ny*understory_xylem_nz))
  allocate(understory_xylem_mesh%filter          (soil_nx*soil_ny*understory_xylem_nz))

  allocate(understory_leaf_pp%por              (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%perm             (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%lambda           (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%alpha            (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%eff_porosity     (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%residual_sat     (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%satfunc_type     (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%relperm_type     (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%relperm_param_1  (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_pp%relperm_param_2  (soil_nx*soil_ny*understory_leaf_nz))

  allocate(understory_leaf_mesh%xc               (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_mesh%yc               (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_mesh%zc               (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_mesh%area             (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_mesh%vol              (soil_nx*soil_ny*understory_leaf_nz))
  allocate(understory_leaf_mesh%filter           (soil_nx*soil_ny*understory_leaf_nz))

  allocate(understory_root_mesh%id  (soil_nx, soil_ny, understory_root_nz  ))
  allocate(understory_xylem_mesh%id (soil_nx, soil_ny, understory_xylem_nz ))
  allocate(understory_leaf_mesh%id  (soil_nx, soil_ny, understory_leaf_nz  ))
  
  if (multi_goveqns_formulation) then
     id_value = 0
  else
     id_value = soil_ncells + soil_nx * soil_ny *( overstory_root_nz  + overstory_xylem_nz  + overstory_leaf_nz)
  end if

  ! Add root grid cells
  count = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny

        do kk = 1, understory_root_nz
           id_value                                = id_value + 1
           understory_root_mesh%id           (ii,jj,kk) = id_value

           count                                   = count + 1
           understory_root_mesh%xc              (count) = soil_xc3d(ii,jj,1)
           understory_root_mesh%yc              (count) = soil_yc3d(ii,jj,1)
           understory_root_mesh%zc              (count) = elevation(ii,jj) - soil_dz/2.d0 - soil_dz*(kk-1)
           understory_root_mesh%area            (count) = understory_root_area_profile(kk)
           understory_root_mesh%vol             (count) = understory_root_vol_profile(kk)
           understory_root_mesh%filter          (count) = 1

           understory_root_pp%por             (count) = 0.d0
           understory_root_pp%perm            (count) = 0.d0
           understory_root_pp%alpha           (count) = understory_xylem_phi0
           understory_root_pp%lambda          (count) = understory_xylem_p
           understory_root_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           understory_root_pp%relperm_param_1 (count) = understory_xylem_vulnerability_d * grav * denh2o
           understory_root_pp%relperm_param_2 (count) = understory_xylem_vulnerability_c
           understory_root_pp%residual_sat    (count) = 0.d0
           understory_root_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  ! Add xylem grid cells
  count = 0
  if (multi_goveqns_formulation) id_value = 0

  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, understory_xylem_nz
           id_value                                 = id_value + 1
           understory_xylem_mesh%id           (ii,jj,kk) = id_value

           count                                         = count + 1
           understory_xylem_mesh%xc              (count) = soil_xc3d(ii,jj,1)
           understory_xylem_mesh%yc              (count) = soil_yc3d(ii,jj,1)
           understory_xylem_mesh%zc              (count) = elevation(ii,jj) + soil_dz/2.d0 + (kk-1)*soil_dz
           understory_xylem_mesh%area            (count) = understory_xylem_area_profile(kk)
           understory_xylem_mesh%vol             (count) = understory_xylem_area_profile(kk) * soil_dz
           understory_xylem_mesh%filter          (count) = 1

           understory_xylem_pp%por             (count) = understory_xylem_porosity
           understory_xylem_pp%perm            (count) = understory_xylem_Kmax * vish2o / (denh2o * grav)
           understory_xylem_pp%alpha           (count) = understory_xylem_phi0
           understory_xylem_pp%lambda          (count) = understory_xylem_p
           understory_xylem_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           understory_xylem_pp%relperm_param_1 (count) = understory_xylem_vulnerability_d * grav * denh2o
           understory_xylem_pp%relperm_param_2 (count) = understory_xylem_vulnerability_c
           understory_xylem_pp%residual_sat    (count) = 0.d0
           understory_xylem_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

  ! Add branch grid cells
  count = 0
  if (multi_goveqns_formulation) id_value = 0

  do ii = 1, soil_nx
     do jj = 1, soil_ny
        do kk = 1, understory_leaf_nz
           id_value                                = id_value + 1
           understory_leaf_mesh%id           (ii,jj,kk) = id_value

           count                                   = count + 1
           idx                                     = understory_branch_2_xylem_index(kk)

           understory_leaf_mesh%xc              (count) = soil_xc3d(ii,jj,1)-understory_branch_length_profile(idx)
           understory_leaf_mesh%yc              (count) = soil_yc3d(ii,jj,1)
           understory_leaf_mesh%zc              (count) = elevation(ii,jj) + soil_dz/2.d0 + (kk-1)*soil_dz + (understory_xylem_nz-understory_leaf_nz)*soil_dz
           understory_leaf_mesh%area            (count) = understory_xylem_area_profile(kk) * understory_branch_area_ratio

           understory_leaf_mesh%vol             (count) = understory_leaf_mesh%area(count) * understory_branch_length_profile(idx)
           understory_leaf_mesh%filter          (count) = 1

           understory_leaf_pp%por             (count) = 0.d0
           understory_leaf_pp%perm            (count) = understory_xylem_Kmax * vish2o / (denh2o * grav)
           understory_leaf_pp%alpha           (count) = understory_xylem_phi0
           understory_leaf_pp%lambda          (count) = understory_xylem_p
           understory_leaf_pp%relperm_type    (count) = RELPERM_FUNC_WEIBULL
           understory_leaf_pp%relperm_param_1 (count) = understory_xylem_vulnerability_d * grav * denh2o
           understory_leaf_pp%relperm_param_2 (count) = understory_xylem_vulnerability_c
           understory_leaf_pp%residual_sat    (count) = 0.d0
           understory_leaf_pp%satfunc_type    (count) = SAT_FUNC_CHUANG
        end do
     end do
  end do

end subroutine setup_understory_mesh

!------------------------------------------------------------------------
subroutine add_goveqn()
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

end subroutine add_goveqn

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
  use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
  use MultiPhysicsProbConstants , only : CONN_VERTICAL
  use ConnectionSetType         , only : connection_set_type
  use ConnectionSetType         , only : ConnectionSetDestroy
  use MeshType                  , only : MeshCreateConnectionSet
  use petscsys
  !
  ! !ARGUMENTS
  implicit none
  !
  PetscInt                            :: ieqn

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
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoils
  use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
  use EOSWaterMod               , only : DENSITY_TGDPB01
  use EOSWaterMod               , only : DENSITY_IFC67
  use mpp_varcon                , only : denh2o
  use mpp_varcon                , only : grav
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetDensityType
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunctionAuxVarConn
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
  use problem_parameters
  use petscsys
  !
#include <petsc/finclude/petsc.h>
  !
  implicit none
  !
  PetscInt                     :: igoveqn
  PetscBool, pointer           :: set_upwind_auxvar(:)
  !-----------------------------------------------------------------------

  igoveqn = 1

  call VSFMMPPSetSoilPorosity(vsfm_mpp, igoveqn, vsfm_por)

  call VSFMMPPSetSaturationFunction(vsfm_mpp, igoveqn, vsfm_satfunc_type, &
         vsfm_alpha, vsfm_lambda, vsfm_residual_sat)

  call VSFMMPPSetSoilPermeability(vsfm_mpp, igoveqn, vsfm_perm, vsfm_perm, vsfm_perm)
  call VSFMMPPSetDensityType(vsfm_mpp, 1, DENSITY_IFC67)

  allocate(set_upwind_auxvar(soil_mesh_nconn + overstory_nconn + understory_nconn))

  set_upwind_auxvar(:) = PETSC_TRUE
  call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       set_upwind_auxvar, conn_satparam_up_itype, conn_satparam_up_param_1, &
       conn_satparam_up_param_2, conn_satparam_up_param_3)

  call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       set_upwind_auxvar, conn_relperm_up_itype, conn_relperm_up_param_1, &
       conn_relperm_up_param_2, conn_relperm_up_param_2)

  set_upwind_auxvar(:) = PETSC_FALSE
  call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       set_upwind_auxvar, conn_satparam_dn_itype, conn_satparam_dn_param_1, &
       conn_satparam_dn_param_2, conn_satparam_dn_param_3)

  call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       set_upwind_auxvar, conn_relperm_dn_itype, conn_relperm_dn_param_1, &
       conn_relperm_dn_param_2, conn_relperm_dn_param_3)

  deallocate(vsfm_por          )
  deallocate(vsfm_perm         )
  deallocate(vsfm_lambda       )
  deallocate(vsfm_alpha        )
  deallocate(vsfm_eff_porosity )
  deallocate(vsfm_residual_sat )
  deallocate(vsfm_satfunc_type      )

end subroutine set_material_properties

!------------------------------------------------------------------------
subroutine set_initial_conditions()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use problem_parameters
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
  PetscInt           :: ii,jj,kk
  PetscInt           :: count
  PetscReal          :: wtd_for_column
  PetscBool          :: wtd_for_column_unset
  PetscReal, pointer :: press_ic(:)
  PetscErrorCode     :: ierr

  call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)

  press_ic(:) = 91325.d0
  count = 0
  do ii = 1, soil_nx
     do jj = 1, soil_ny
        wtd_for_column_unset = PETSC_TRUE
        do kk = 1, soil_nz
           if (wtd_for_column_unset .and. soil_id(ii,jj,kk) > 0 ) then
              wtd_for_column = -init_wtd - soil_dz*(kk-1)
              wtd_for_column_unset = PETSC_FALSE
           end if
           if (soil_id(ii,jj,kk) > 0) then
              count = count + 1
              press_ic(count) = 101325 + (wtd_for_column - soil_zc3d(ii,jj,kk))* 1000.d0 * 9.81d0
           endif
        end do
     end do
  end do
  call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
  call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
  call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

end subroutine set_initial_conditions

!------------------------------------------------------------------------
subroutine set_conn_flux_type()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbVSFM      , only : vsfm_mpp
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability 
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnIntValue
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnRealValue
  use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunctionAuxVarConn
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_INTERNAL
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
  use MultiPhysicsProbConstants , only : VAR_FLUX_TYPE
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_TYPE
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_UP
  use MultiPhysicsProbConstants , only : VAR_CONDUCTANCE_DN
  use MultiPhysicsProbConstants , only : CONDUCTANCE_FLUX_TYPE
  use MultiPhysicsProbConstants , only : DARCY_FLUX_TYPE
  use MultiPhysicsProbConstants , only : VAR_CAMPBELL_HE
  use MultiPhysicsProbConstants , only : VAR_CAMPBELL_N
  use SaturationFunction        , only : RELPERM_FUNC_CAMPBELL
  use MultiPhysicsProbConstants , only : AUXVAR_CONN_BC
  use petscsys
  use problem_parameters
  !
  implicit none
  !
  PetscInt :: igoveqn

  igoveqn = 1

  call VSFMMPPSetRelativePermeability(vsfm_mpp, igoveqn, vsfm_relperm_type, &
       vsfm_relperm_param_1, vsfm_relperm_param_2)

  ! Set connection flux type
  call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       VAR_FLUX_TYPE, conn_flux_type)

  ! Set conductance type
  call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       VAR_CONDUCTANCE_TYPE, conn_cond_type)

  call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       VAR_CONDUCTANCE, conn_cond)

  call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       VAR_CONDUCTANCE_UP, conn_cond_up)

  call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, igoveqn, AUXVAR_CONN_INTERNAL, &
       VAR_CONDUCTANCE_DN, conn_cond_dn)

end subroutine set_conn_flux_type
