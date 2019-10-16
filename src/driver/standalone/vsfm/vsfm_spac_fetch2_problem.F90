module vsfm_spac_fetch2_problem

#include <petsc/finclude/petsc.h>

  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  use MultiPhysicsProbVSFM, only : vsfm_mpp

  implicit none

  PetscReal , parameter :: vis      = 8.904156d-4   ! Pa s
  PetscReal , parameter :: rho      = 1000.d0       ! kg m^{-3}
  PetscReal , parameter :: grav     = 9.81          ! m s^{-2}

  PetscInt  , parameter :: nx       = 1             ! -
  PetscInt  , parameter :: ny       = 1             ! -
  PetscReal , parameter :: dx       = 1.d0          ! m
  PetscReal , parameter :: dy       = 1.d0          ! m
  PetscReal , parameter :: dz_xylem = 0.2d0         ! m
  PetscReal , parameter :: dz_soil  = 0.1d0         ! m

  PetscReal , parameter :: porosity = 0.45d0        ! [-]

  PetscInt  , parameter :: oak_nz       = 59        ! -
  !PetscReal , parameter :: oak_Asapwood = 0.0099d0  ! m^2
  PetscReal , parameter :: oak_Asapwood = 0.0255d0  ! m^2
  !PetscReal , parameter :: oak_Acrown   = 28.8d0    ! m^2
  PetscReal , parameter :: oak_phis50   = -0.35d6   ! Pa
  PetscReal , parameter :: oak_phi50    = -2.5d6    ! Pa
  PetscReal , parameter :: oak_phi88    = -0.5d6    ! Pa
  PetscReal , parameter :: oak_c1       = 1.7d6     ! Pa
  PetscReal , parameter :: oak_c2       = 3.0d0     ! -
  PetscReal , parameter :: oak_c3       = 12.3d0    ! -
  PetscReal , parameter :: oak_kmax     = 1.6d-6    ! s

  PetscInt  , parameter :: pine_nz       = 85       ! -
  !PetscReal , parameter :: pine_Asapwood = 0.0616d0 ! m^2
  PetscReal , parameter :: pine_Asapwood = 0.0142d0 ! m^2
  !PetscReal , parameter :: pine_Acrown   = 46.1d0   ! m^2
  !PetscReal , parameter :: pine_phis50   = -0.18d6  ! Pa 
  PetscReal , parameter :: pine_phis50   = -0.18d6  ! Pa
  PetscReal , parameter :: pine_phi50    = -2.2d6   ! Pa
  PetscReal , parameter :: pine_phi88    = -0.5d6   ! Pa
  PetscReal , parameter :: pine_c1       = 1.2d6    ! Pa
  PetscReal , parameter :: pine_c2       = 5.0d0    ! -
  PetscReal , parameter :: pine_c3       = 10.3d0   ! -
  PetscReal , parameter :: pine_kmax     = 1.2d-6   ! s

  PetscInt  , parameter :: maple_nz       = 85       ! -
  PetscReal , parameter :: maple_Asapwood = 0.0188d0 ! m^2
  PetscReal , parameter :: maple_phis50   = -0.25d6  ! Pa
  PetscReal , parameter :: maple_phi50    = -2.2d6   ! Pa
  PetscReal , parameter :: maple_phi88    = -0.5d6   ! Pa
  PetscReal , parameter :: maple_c1       = 1.2d6    ! Pa
  PetscReal , parameter :: maple_c2       = 5.0d0    ! -
  PetscReal , parameter :: maple_c3       = 10.3d0   ! -
  PetscReal , parameter :: maple_kmax     = 1.2d-6   ! s

  PetscInt  , parameter :: es_nz       = 85       ! -
  PetscReal , parameter :: es_Asapwood = 0.04d0 ! m^2
  PetscReal , parameter :: es_phis50   = -0.3d6  ! Pa
  PetscReal , parameter :: es_phi50    = -2.2d6   ! Pa
  PetscReal , parameter :: es_phi88    = -0.5d6   ! Pa
  PetscReal , parameter :: es_c1       = 1.2d6    ! Pa
  PetscReal , parameter :: es_c2       = 5.0d0    ! -
  PetscReal , parameter :: es_c3       = 10.3d0   ! -
  PetscReal , parameter :: es_kmax     = 1.2d-6   ! s

  ! Parameters for root length density = length-of-root/volume-of-soil  [m_root/m^3_soil]
  PetscInt , parameter :: oak_root_nz        = 60
  PetscReal, parameter :: oak_root_qz        = 3.d0         ! [-]
  !PetscReal, parameter :: oaK_root_d         = 3.d0         ! [m]
  PetscReal, parameter :: oaK_root_d         = 7.d0         ! [m]
  PetscReal, parameter :: oak_root_rld0      = 4.d4         ! [m/m^3]
  PetscReal, parameter :: oak_root_radius    = 2.d-3        ! [m]

  PetscReal , parameter :: oak_rad_root_kmax = 1.2d-6   !
  PetscReal , parameter :: oak_rad_root_phi50= -2.2d6   !
  PetscReal , parameter :: oak_rad_root_phi88= -0.5d6   !
  PetscReal , parameter :: oak_rad_root_c1   = 1.2d6    !
  PetscReal , parameter :: oak_rad_root_c2   = 5.0d0    ! -

  PetscReal , parameter :: oak_axi_root_kmax = 1.2d-6   !
  PetscReal , parameter :: oak_axi_root_phi50= -2.2d6   !
  PetscReal , parameter :: oak_axi_root_phi88= -0.5d6   !
  PetscReal , parameter :: oak_axi_root_c1   = 1.2d6    !
  PetscReal , parameter :: oak_axi_root_c2   = 5.0d0    ! -

  PetscInt , parameter :: pine_root_nz        = 60
  PetscReal, parameter :: pine_root_qz        = 3.d0         ! [-]
  PetscReal, parameter :: pine_root_d         = 7.d0         ! [m]
  PetscReal, parameter :: pine_root_rld0      = 4.d4         ! [m/m^3]
  PetscReal, parameter :: pine_root_radius    = 2.d-2        ! [m]

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

  PetscInt , parameter :: maple_root_nz        = 30
  PetscReal, parameter :: maple_root_qz        = 3.d0         ! [-]
  PetscReal, parameter :: maple_root_d         = 7.d0         ! [m]
  PetscReal, parameter :: maple_root_rld0      = 4.d4         ! [m/m^3]
  PetscReal, parameter :: maple_root_radius    = 2.d-2        ! [m]

  PetscReal , parameter :: maple_rad_root_kmax = 1.2d-6   !
  PetscReal , parameter :: maple_rad_root_phi50= -2.2d6   !
  PetscReal , parameter :: maple_rad_root_phi88= -0.5d6   !
  PetscReal , parameter :: maple_rad_root_c1   = 1.2d6    !
  PetscReal , parameter :: maple_rad_root_c2   = 5.0d0    ! -

  PetscReal , parameter :: maple_axi_root_kmax = 1.2d-6   !
  PetscReal , parameter :: maple_axi_root_phi50= -2.2d6   !
  PetscReal , parameter :: maple_axi_root_phi88= -0.5d6   !
  PetscReal , parameter :: maple_axi_root_c1   = 1.2d6    !
  PetscReal , parameter :: maple_axi_root_c2   = 5.0d0    ! -

  PetscInt , parameter :: es_root_nz        = 60
  PetscReal, parameter :: es_root_qz        = 3.d0         ! [-]
  PetscReal, parameter :: es_root_d         = 7.d0         ! [m]
  PetscReal, parameter :: es_root_rld0      = 4.d4         ! [m/m^3]
  PetscReal, parameter :: es_root_radius    = 2.d-2        ! [m]

  PetscReal , parameter :: es_rad_root_kmax = 1.2d-6   !
  PetscReal , parameter :: es_rad_root_phi50= -2.2d6   !
  PetscReal , parameter :: es_rad_root_phi88= -0.5d6   !
  PetscReal , parameter :: es_rad_root_c1   = 1.2d6    !
  PetscReal , parameter :: es_rad_root_c2   = 5.0d0    ! -

  PetscReal , parameter :: es_axi_root_kmax = 1.2d-6   !
  PetscReal , parameter :: es_axi_root_phi50= -2.2d6   !
  PetscReal , parameter :: es_axi_root_phi88= -0.5d6   !
  PetscReal , parameter :: es_axi_root_c1   = 1.2d6    !
  PetscReal , parameter :: es_axi_root_c2   = 5.0d0    ! -

  PetscInt , parameter :: soil_nz         = 60
  PetscReal, parameter :: soil_perm       = 6.83d-08      ! [m^2]
  PetscReal, parameter :: soil_perm_bot   = 6.83d-08      ! [m^2]
  PetscReal, parameter :: soil_sat_res    = 0.02d0        ! [-]
  !PetscReal, parameter :: soil_alpha      = 0.00005d0     ! [Pa^{-1}]
  PetscReal, parameter :: soil_alpha      = 0.00035d0     ! [Pa^{-1}]
  !PetscReal, parameter :: soil_alpha      = 0.001d0     ! [Pa^{-1}]
  !PetscReal, parameter :: soil_alpha      = 0.0052d0     ! [Pa^{-1}]
  PetscReal, parameter :: soil_vg_m       = 0.40d0        ! [-]
  PetscReal, parameter :: soil_por        = 0.50d0         ! [-]
  PetscReal, parameter :: soil_por_bot    = 0.50d0         ! [-]

  ! Index for various trees
  PetscInt   :: E_IDX ! Early Successional
  PetscInt   :: M_IDX ! Maple
  PetscInt   :: O_IDX ! Oak
  PetscInt   :: P_IDX ! Pine
  
  PetscInt   :: nz(4)
  PetscReal  :: Asapwood(4)
  PetscReal  :: xylem_porosity(4)
  PetscReal  :: root_porosity(4)
  PetscReal  :: phis50(4)
  PetscReal  :: phi50(4)
  PetscReal  :: phi88(4)
  PetscReal  :: c1(4)
  PetscReal  :: c2(4)
  PetscReal  :: c3(4)
  PetscReal  :: kmax(4)

  PetscInt  :: root_nz(4)
  PetscReal :: root_qz(4)
  PetscReal :: root_d(4)
  PetscReal :: root_rld0(4)
  PetscReal :: root_radius(4)

  PetscReal  :: rad_root_kmax(4)
  PetscReal  :: rad_root_phi50(4)
  PetscReal  :: rad_root_phi88(4)
  PetscReal  :: rad_root_c1(4)
  PetscReal  :: rad_root_c2(4)

  PetscReal  :: axi_root_kmax(4)
  PetscReal  :: axi_root_phi50(4)
  PetscReal  :: axi_root_phi88(4)
  PetscReal  :: axi_root_c1(4)
  PetscReal  :: axi_root_c2(4)

  PetscInt :: neqns

  character(len=256) :: problem_type
  character(len=256) :: goveqn_name(9)

  character(len=256) :: ic_filename
  PetscBool          :: ic_file_specified

  character(len=256) :: sm_ic_filename
  PetscBool          :: sm_ic_file_specified

  PetscInt              :: ncells_local

  PetscReal             :: theta_s
  PetscReal             :: VG_alpha
  PetscReal             :: VG_n
  PetscReal , parameter :: PI  = 4 * atan (1.0_8)

  PetscInt :: GE_oak_xylm
  PetscInt :: GE_oak_root
  PetscInt :: GE_oak_pine_xylm
  PetscInt :: GE_pine_xylm
  PetscInt :: GE_pine_root
  PetscInt :: GE_e_xylm
  PetscInt :: GE_e_root
  PetscInt :: GE_m_xylm
  PetscInt :: GE_m_root
  PetscInt :: GE_o_xylm
  PetscInt :: GE_o_root
  PetscInt :: GE_p_xylm
  PetscInt :: GE_p_root
  PetscInt :: GE_soil

  PetscBool          :: soil_bc_specified
  PetscBool          :: sm_bc_specified
  PetscBool          :: soil_ss_specified

  public :: run_vsfm_spac_fetch2_problem
  
contains

  !------------------------------------------------------------------------
  subroutine SetUpTreeProperties()
    !
    implicit none
    !

    E_IDX = 1; M_IDX = 2; O_IDX = 3; P_IDX = 4

    nz             (E_IDX) = 110               ; nz             (M_IDX) = 110                  ; nz             (O_IDX) = 110                ; nz             (P_IDX) = 110                 ;
    Asapwood       (E_IDX) = es_Asapwood       ; Asapwood       (M_IDX) = maple_Asapwood       ; Asapwood       (O_IDX) = oak_Asapwood       ; Asapwood       (P_IDX) = pine_Asapwood       ;
    !xylem_porosity (E_IDX) = porosity          ; xylem_porosity (M_IDX) = porosity             ; xylem_porosity (O_IDX) = porosity           ; xylem_porosity (P_IDX) = porosity            ;
    !root_porosity  (E_IDX) = porosity          ; root_porosity  (M_IDX) = porosity             ; root_porosity  (O_IDX) = porosity           ; root_porosity  (P_IDX) = porosity            ;
    !xylem_porosity (E_IDX) = 0.d0              ; xylem_porosity (M_IDX) = 0.d0                 ; xylem_porosity (O_IDX) = 0.d0               ; xylem_porosity (P_IDX) = 0.d0            ;
    !root_porosity  (E_IDX) = 0.d0              ; root_porosity  (M_IDX) = 0.d0                 ; root_porosity  (O_IDX) = 0.d0               ; root_porosity  (P_IDX) = 0.d0            ;
    xylem_porosity (E_IDX) = 1.d0              ; xylem_porosity (M_IDX) = 1.d0                 ; xylem_porosity (O_IDX) = 1.d0               ; xylem_porosity (P_IDX) = 1.d0            ;
    root_porosity  (E_IDX) = 1.d0              ; root_porosity  (M_IDX) = 1.d0                 ; root_porosity  (O_IDX) = 1.d0               ; root_porosity  (P_IDX) = 1.d0            ;
    phis50         (E_IDX) = es_phis50         ; phis50         (M_IDX) = maple_phis50         ; phis50         (O_IDX) = oak_phis50         ; phis50         (P_IDX) = pine_phis50         ;
    phi50          (E_IDX) = es_phi50          ; phi50          (M_IDX) = maple_phi50          ; phi50          (O_IDX) = oak_phi50          ; phi50          (P_IDX) = pine_phi50          ;
    phi88          (E_IDX) = es_phi88          ; phi88          (M_IDX) = maple_phi88          ; phi88          (O_IDX) = oak_phi88          ; phi88          (P_IDX) = pine_phi88          ;
    c1             (E_IDX) = es_c1             ; c1             (M_IDX) = maple_c1             ; c1             (O_IDX) = oak_c1             ; c1             (P_IDX) = pine_c1             ;
    c2             (E_IDX) = es_c2             ; c2             (M_IDX) = maple_c2             ; c2             (O_IDX) = oak_c2             ; c2             (P_IDX) = pine_c2             ;
    c3             (E_IDX) = es_c3             ; c3             (M_IDX) = maple_c3             ; c3             (O_IDX) = oak_c3             ; c3             (P_IDX) = pine_c3             ;
    kmax           (E_IDX) = es_kmax           ; kmax           (M_IDX) = maple_kmax           ; kmax           (O_IDX) = oak_kmax           ; kmax           (P_IDX) = pine_kmax           ;
    root_nz        (E_IDX) = es_root_nz        ; root_nz        (M_IDX) = maple_root_nz        ; root_nz        (O_IDX) = oak_root_nz        ; root_nz        (P_IDX) = pine_root_nz        ;
    root_qz        (E_IDX) = es_root_qz        ; root_qz        (M_IDX) = maple_root_qz        ; root_qz        (O_IDX) = oak_root_qz        ; root_qz        (P_IDX) = pine_root_qz        ;
    root_d         (E_IDX) = es_root_d         ; root_d         (M_IDX) = maple_root_d         ; root_d         (O_IDX) = oak_root_d         ; root_d         (P_IDX) = pine_root_d         ;
    root_rld0      (E_IDX) = es_root_rld0      ; root_rld0      (M_IDX) = maple_root_rld0      ; root_rld0      (O_IDX) = oak_root_rld0      ; root_rld0      (P_IDX) = pine_root_rld0      ;
    root_radius    (E_IDX) = es_root_radius    ; root_radius    (M_IDX) = maple_root_radius    ; root_radius    (O_IDX) = oak_root_radius    ; root_radius    (P_IDX) = pine_root_radius    ;
    rad_root_kmax  (E_IDX) = es_rad_root_kmax  ; rad_root_kmax  (M_IDX) = maple_rad_root_kmax  ; rad_root_kmax  (O_IDX) = oak_rad_root_kmax  ; rad_root_kmax  (P_IDX) = pine_rad_root_kmax  ;
    rad_root_phi50 (E_IDX) = es_rad_root_phi50 ; rad_root_phi50 (M_IDX) = maple_rad_root_phi50 ; rad_root_phi50 (O_IDX) = oak_rad_root_phi50 ; rad_root_phi50 (P_IDX) = pine_rad_root_phi50 ;
    rad_root_phi88 (E_IDX) = es_rad_root_phi88 ; rad_root_phi88 (M_IDX) = maple_rad_root_phi88 ; rad_root_phi88 (O_IDX) = oak_rad_root_phi88 ; rad_root_phi88 (P_IDX) = pine_rad_root_phi88 ;
    rad_root_c1    (E_IDX) = es_rad_root_c1    ; rad_root_c1    (M_IDX) = maple_rad_root_c1    ; rad_root_c1    (O_IDX) = oak_rad_root_c1    ; rad_root_c1    (P_IDX) = pine_rad_root_c1    ;
    rad_root_c2    (E_IDX) = es_rad_root_c2    ; rad_root_c2    (M_IDX) = maple_rad_root_c2    ; rad_root_c2    (O_IDX) = oak_rad_root_c2    ; rad_root_c2    (P_IDX) = pine_rad_root_c2    ;
    axi_root_kmax  (E_IDX) = es_axi_root_kmax  ; axi_root_kmax  (M_IDX) = maple_axi_root_kmax  ; axi_root_kmax  (O_IDX) = oak_axi_root_kmax  ; axi_root_kmax  (P_IDX) = pine_axi_root_kmax  ;
    axi_root_phi50 (E_IDX) = es_axi_root_phi50 ; axi_root_phi50 (M_IDX) = maple_axi_root_phi50 ; axi_root_phi50 (O_IDX) = oak_axi_root_phi50 ; axi_root_phi50 (P_IDX) = pine_axi_root_phi50 ;
    axi_root_phi88 (E_IDX) = es_axi_root_phi88 ; axi_root_phi88 (M_IDX) = maple_axi_root_phi88 ; axi_root_phi88 (O_IDX) = oak_axi_root_phi88 ; axi_root_phi88 (P_IDX) = pine_axi_root_phi88 ;
    axi_root_c1    (E_IDX) = es_axi_root_c1    ; axi_root_c1    (M_IDX) = maple_axi_root_c1    ; axi_root_c1    (O_IDX) = oak_axi_root_c1    ; axi_root_c1    (P_IDX) = pine_axi_root_c1    ;
    axi_root_c2    (E_IDX) = es_axi_root_c2    ; axi_root_c2    (M_IDX) = maple_axi_root_c2    ; axi_root_c2    (O_IDX) = oak_axi_root_c2    ; axi_root_c2    (P_IDX) = pine_axi_root_c2    ;

end subroutine SetUpTreeProperties

  !------------------------------------------------------------------------
  subroutine run_vsfm_spac_fetch2_problem()
    !
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use mpp_varpar           , only : mpp_varpar_init
    use SystemOfEquationsVSFMType, only : sysofeqns_vsfm_type
    use SystemOfEquationsBaseType, only : sysofeqns_base_type
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
    PetscBool          :: save_initial_soln
    PetscBool          :: save_final_soln
    PetscBool          :: save_actual_et
    PetscBool          :: save_sat
    PetscBool          :: save_mass
    PetscBool          :: save_pressure
    PetscBool          :: save_internal_mass_flux
    PetscBool          :: save_boundary_mass_flux
    PetscBool          :: save_coupling_mass_flux
    character(len=256) :: string
    character(len=256) :: pet_file
    character(len=256) :: soil_bc_file
    character(len=256) :: soil_ss_file
    character(len=256) :: sm_bc_file
    character(len=256) :: actual_et_file
    character(len=256) :: sat_file
    character(len=256) :: mass_file
    character(len=256) :: pressure_file
    character(len=256) :: internal_mass_flux_file
    character(len=256) :: boundary_mass_flux_file
    character(len=256) :: coupling_mass_flux_file
    Vec                :: ET
    Vec                :: SoilBC
    Vec                :: SoilMoistureBC
    Vec                :: SoilSS
    Vec                :: Actual_ET
    Vec                :: Sat
    Vec                :: Mass
    Vec                :: Pressure
    Vec                :: Internal_Mass_Flux
    Vec                :: Boundary_Mass_Flux
    Vec                :: Coupling_Mass_Flux
    PetscInt           :: nsat
    PetscInt           :: num_internal_mass_flux
    PetscInt           :: num_boundary_mass_flux
    PetscInt           :: num_coupling_mass_flux
    PetscInt           :: vec_size
    PetscInt           :: ntimes
    PetscReal, pointer :: act_et_p(:)
    PetscReal, pointer :: sat_p(:)
    PetscReal, pointer :: mass_p(:)
    PetscReal, pointer :: pressure_p(:)
    PetscReal, pointer :: internal_mass_flux_p(:)
    PetscReal, pointer :: boundary_mass_flux_p(:)
    PetscReal, pointer :: coupling_mass_flux_p(:)
    PetscViewer        :: viewer
    PetscBool          :: error_in_cmd_options
    class(sysofeqns_base_type), pointer :: soe

    ! Set default settings
    dtime                  = 1800.d0    ! [s]
    nstep                  = 24
    save_initial_soln      = PETSC_FALSE
    save_final_soln        = PETSC_FALSE
    save_actual_et         = PETSC_FALSE
    save_sat               = PETSC_FALSE
    save_mass              = PETSC_FALSE
    save_pressure          = PETSC_FALSE
    save_internal_mass_flux= PETSC_FALSE
    save_boundary_mass_flux= PETSC_FALSE
    save_coupling_mass_flux= PETSC_FALSE
    error_in_cmd_options   = PETSC_FALSE
    problem_type           = 'oak'
    ic_file_specified      = PETSC_FALSE
    sm_ic_file_specified   = PETSC_FALSE
    soil_bc_specified      = PETSC_FALSE
    soil_ss_specified      = PETSC_FALSE
    sm_bc_specified        = PETSC_FALSE

    !
    call SetUpTreeProperties()
    
    ! Get some command line options

    call PetscOptionsGetInt  (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-nstep', nstep,flg, ierr)

    call PetscOptionsGetBool (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-save_initial_soln', save_initial_soln, flg, ierr)

    call PetscOptionsGetBool (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-save_final_soln', save_final_soln, flg,ierr)

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-actual_et_file', actual_et_file, flg, ierr)
    if (flg) then
       save_actual_et = PETSC_TRUE
    end if

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-sat_file', sat_file, flg, ierr)
    if (flg) save_sat = PETSC_TRUE

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-mass_file', mass_file, flg, ierr)
    if (flg) save_mass = PETSC_TRUE

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-pressure_file', pressure_file, flg, ierr)
    if (flg) save_pressure = PETSC_TRUE

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-internal_mass_flux_file', internal_mass_flux_file, flg, ierr)
    if (flg) then
       save_internal_mass_flux = PETSC_TRUE
    end if

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-boundary_mass_flux_file', boundary_mass_flux_file, flg, ierr)
    if (flg) then
       save_boundary_mass_flux = PETSC_TRUE
    end if

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-coupling_mass_flux_file', coupling_mass_flux_file, flg, ierr)
    if (flg) then
       save_coupling_mass_flux = PETSC_TRUE
    end if

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-ic_file', ic_filename, flg, ierr)
    if (flg) ic_file_specified = PETSC_TRUE

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-sm_ic_file', sm_ic_filename, flg, ierr)
    if (flg) sm_ic_file_specified = PETSC_TRUE

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
       case ('emop_spac')
       case ('e_spac')
       case ('m_spac')
       case ('o_spac')
       case ('p_spac')
       case ('e')
       case ('m')
       case ('o')
       case ('p')
       case ('soil')
       case default
          write(*,*)'ERROR: Unknown option specified via -problem_type <value>.'
          write(*,*)'       Valid value for -problem_type are:'
          write(*,*)'         oak'
          write(*,*)'         pine'
          write(*,*)'         oak_and_pine'
          write(*,*)'         oak_spac'
          write(*,*)'         pine_spac'
          write(*,*)'         oak_pine_spac'
          write(*,*)'         soil'
          error_in_cmd_options = PETSC_TRUE
       end select
    end if

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-soil_bc_file',soil_bc_file,flg,ierr)
    if (flg) soil_bc_specified = PETSC_TRUE
    
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-sm_bc_file',sm_bc_file,flg,ierr)
    if (flg) sm_bc_specified = PETSC_TRUE

    if ( trim(problem_type) == 'oak'          .or. &
         trim(problem_type) == 'pine'         .or. &
         trim(problem_type) == 'oak_and_pine' .or. &
         trim(problem_type) == 'e'            .or. &
         trim(problem_type) == 'm'            .or. &
         trim(problem_type) == 'o'            .or. &
         trim(problem_type) == 'p'                 &
         ) then
       if ((.not.soil_bc_specified) .and. (.not.sm_bc_specified)) then
          write(*,*)'ERROR: Specify the soil boundary condition filename via either options::'
          write(*,*)'       -soil_bc_file <filename> For pressure BC'
          write(*,*)'       -sm_bc_file   <filename> For soil moisture BC'
          error_in_cmd_options = PETSC_TRUE
       end if
    end if

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-soil_ss_file',soil_ss_file,flg,ierr)
    if (flg) soil_ss_specified = PETSC_TRUE

    if ( trim(problem_type) == 'oak_spac'      .or. &
         trim(problem_type) == 'pine_spac'     .or. &
         trim(problem_type) == 'oak_pine_spac' .or. &
         trim(problem_type) == 'soil'               &
         !trim(problem_type) == 'emop_spac'           &
         !trim(problem_type) == 'e_spac'        .or. &
         !trim(problem_type) == 'm_spac'        .or. &
         !trim(problem_type) == 'o_spac'        .or. &
         !trim(problem_type) == 'p_spac'             &
         ) then
       if (.not.flg) then
          write(*,*)'ERROR: Specify the soil source-sink condition filename via -soil_ss_file <filename>'
          error_in_cmd_options = PETSC_TRUE
       end if
    end if

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-pet_file', pet_file, flg, ierr)
    select case(trim(problem_type))
    case ('soil')
    case default
       if (.not.flg) then
          write(*,*)'ERROR: Specify the potential ET filename via -pet_file <filename>'
          error_in_cmd_options = PETSC_TRUE
       end if
    end select

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
    case ('emop_spac')
       nET      = nz(E_IDX) + nz(M_IDX) + nz(O_IDX) + nz(P_IDX)
    case ('e_spac')
       nET      = nz(E_IDX)
    case ('m_spac')
       nET      = nz(M_IDX)
    case ('o_spac')
       nET      = nz(O_IDX)
    case ('p_spac')
       nET      = nz(P_IDX)
    case ('e')
       nET      = nz(E_IDX)
    case ('m')
       nET      = nz(M_IDX)
    case ('o')
       nET      = nz(O_IDX)
    case ('p')
       nET      = nz(P_IDX)
    case ('soil')
       nET      = 0
    case default
       write(*,*)'Unable to determine nET for problem_type = ',trim(problem_type)
       stop
    end select

    ! Read the PET file
    if (nET>0) then
       write(*,*)'pet : ',trim(pet_file)
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
    endif

    ! Read the soil bc file
    if (soil_bc_specified) then
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
    end if

    write(*,*)'sm_bc_specified: ',sm_bc_specified
    write(*,*)trim(sm_bc_file)
    if (sm_bc_specified) then
       call PetscViewerBinaryOpen(PETSC_COMM_WORLD, sm_bc_file, FILE_MODE_READ, viewer, ierr); CHKERRQ(ierr)
       call VecCreate(PETSC_COMM_WORLD, SoilMoistureBC, ierr); CHKERRQ(ierr)
       call VecLoad(SoilMoistureBC, viewer, ierr); CHKERRQ(ierr)
       call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
       call VecGetSize(SoilMoistureBC, vec_size, ierr); CHKERRQ(ierr)

       ntimes = vec_size
       if (nstep > ntimes) then
          write(*,*)'ERROR: Value for nstep exceeds the number of timesteps ' // &
               'for which soil moisture boundary condition data is available.'
          write(*,*)'nstep = ',nstep
          write(*,*)'ntimes= ',ntimes
          stop
       end if
    end if

    ! Read the soil ss file
    if (soil_ss_specified) then
       call PetscViewerBinaryOpen(PETSC_COMM_WORLD, soil_ss_file, FILE_MODE_READ, viewer, ierr); CHKERRQ(ierr)
       call VecCreate(PETSC_COMM_WORLD, SoilSS, ierr); CHKERRQ(ierr)
       call VecLoad(SoilSS, viewer, ierr); CHKERRQ(ierr)
       call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
       call VecGetSize(SoilSS, vec_size, ierr); CHKERRQ(ierr)

       ntimes = vec_size
       if (nstep > ntimes) then
          write(*,*)'ERROR: Value for nstep exceeds the number of timesteps ' // &
               'for which soil source-sink condition data is available.'
          write(*,*)'nstep = ',nstep
          write(*,*)'ntimes= ',ntimes
          stop
       end if
    end if

    ! Initialize the problem
    call Init()

    if (save_initial_soln) then
       string = 'initial_soln.bin'
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(vsfm_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    time = 0.d0
    if (save_actual_et) then
       call VecCreateSeq(PETSC_COMM_SELF, nstep*nET, Actual_ET, ierr)
       call VecGetArrayF90(Actual_ET, act_et_p, ierr)
    endif

    if (save_sat) then
       call VecGetSize(vsfm_mpp%soe%solver%soln, nsat, ierr)
       call VecCreateSeq(PETSC_COMM_SELF, nstep*nsat, Sat, ierr)
       call VecGetArrayF90(Sat, sat_p, ierr)
    endif

    if (save_mass) then
       call VecGetSize(vsfm_mpp%soe%solver%soln, nsat, ierr)
       call VecCreateSeq(PETSC_COMM_SELF, nstep*nsat, Mass, ierr)
       call VecGetArrayF90(Mass, mass_p, ierr)
    endif

    if (save_pressure) then
       call VecGetSize(vsfm_mpp%soe%solver%soln, nsat, ierr)
       call VecCreateSeq(PETSC_COMM_SELF, nstep*nsat, Pressure, ierr)
       call VecGetArrayF90(Pressure, pressure_p, ierr)
    endif

    if (save_internal_mass_flux) then
       soe => vsfm_mpp%soe
       num_internal_mass_flux = 0
       select type(soe)
       class is (sysofeqns_vsfm_type)
          num_internal_mass_flux = soe%num_auxvars_conn_in
       end select
       call VecCreateSeq(PETSC_COMM_SELF, nstep*num_internal_mass_flux, Internal_Mass_Flux, ierr)
       call VecGetArrayF90(Internal_Mass_Flux, internal_mass_flux_p, ierr)
    endif

    if (save_boundary_mass_flux) then
       soe => vsfm_mpp%soe
       num_boundary_mass_flux = 0
       select type(soe)
       class is (sysofeqns_vsfm_type)
          num_boundary_mass_flux = soe%num_auxvars_bc
       end select
       call VecCreateSeq(PETSC_COMM_SELF, nstep*num_boundary_mass_flux, Boundary_Mass_Flux, ierr)
       call VecGetArrayF90(Boundary_Mass_Flux, boundary_mass_flux_p, ierr)
    endif

    if (save_coupling_mass_flux) then
       soe => vsfm_mpp%soe
       num_coupling_mass_flux = 0
       select type(soe)
       class is (sysofeqns_vsfm_type)
          num_coupling_mass_flux = soe%num_auxvars_bc_otr_geqn
       end select
       call VecCreateSeq(PETSC_COMM_SELF, nstep*num_coupling_mass_flux, Coupling_Mass_Flux, ierr)
       call VecGetArrayF90(Coupling_Mass_Flux, coupling_mass_flux_p, ierr)
    endif

    do istep = 1, nstep

       write(*,*)'istep = ',istep
       select case(trim(problem_type))
       case ('oak','pine','e','m','o','p')
          call set_source_sink_for_single_tree(nET, istep, ET)
          call set_boundary_conditions_for_single_tree(istep, SoilBC)

       case ('oak_and_pine')
          call set_source_sink_for_single_tree(nET, istep, ET)
          call set_boundary_conditions_for_two_trees(istep, SoilBC)

       case ('oak_spac')
          call set_source_sink_for_single_tree(nET, istep, ET)
          call set_source_sink_for_soil(istep, SoilSS)

       case ('pine_spac')
          call set_source_sink_for_single_tree(nET, istep, ET)
          call set_source_sink_for_soil(istep, SoilSS)

       case ('oak_pine_spac')
          call set_source_sink_for_oak_and_pine(istep, ET)
          call set_source_sink_for_soil(istep, SoilSS)

       case ('soil')
          if (soil_ss_specified) call set_source_sink_for_soil(istep, SoilSS)
          if (soil_bc_specified) call set_boundary_conditions_for_single_tree(istep, SoilBC)

       case ('emop_spac','e_spac','m_spac','o_spac','p_spac')
          call set_source_sink_for_emop(nET, istep, ET)
          if (soil_ss_specified) call set_source_sink_for_soil(istep, SoilSS)
          if (soil_bc_specified) call set_boundary_conditions_for_single_tree(istep, SoilBC)
          if (sm_bc_specified) call set_soil_moisture_boundary_conditions_for_single_tree(istep, SoilMoistureBC)

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

       if (save_actual_et) then
          select case(trim(problem_type))
          case ('oak','pine','oak_and_pine','oak_spac','pine_spac','e','m','o','p')
             call diagnose_actual_sink_for_single_tree(nET, istep, act_et_p)
          case ('oak_pine_spac')
             call diagnose_actual_sink_for_oak_and_pine(istep, act_et_p)
          case ('soil')
          case ('emop_spac','e_spac','m_spac','o_spac','p_spac')
             call diagnose_actual_sink_for_emop(nET, istep, act_et_p)
          case default
             write(*,*)'Unable to diagnose actual ET for problem_type = ' // trim(problem_type)
             stop
          end select
       end if

       if (save_sat) then
          call save_saturation(nsat, istep, sat_p)
       end if

       if (save_mass) then
          call save_liquid_mass(nsat, istep, mass_p)
       end if

       if (save_pressure) then
          call save_pressure_data(nsat, istep, pressure_p)
       end if

       if (save_internal_mass_flux) then
          call save_internal_mass_flux_data(vsfm_mpp%soe, num_internal_mass_flux, istep, internal_mass_flux_p)
       end if

       if (save_boundary_mass_flux) then
          call save_boundary_mass_flux_data(vsfm_mpp%soe, num_boundary_mass_flux, istep, boundary_mass_flux_p)
       end if

       if (save_coupling_mass_flux) then
          call save_coupling_mass_flux_data(vsfm_mpp%soe, num_coupling_mass_flux, istep, coupling_mass_flux_p)
       end if

    end do

    if (save_actual_et) then
       call VecRestoreArrayF90(Actual_ET, act_et_p, ierr)
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,actual_et_file,FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(Actual_ET,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
       call VecDestroy(Actual_ET, ierr); CHKERRQ(ierr)
    end if

    if (save_sat) then
       call VecRestoreArrayF90(Sat, sat_p, ierr)
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,sat_file,FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(Sat,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
       call VecDestroy(Sat, ierr); CHKERRQ(ierr)
    end if

    if (save_mass) then
       call VecRestoreArrayF90(Mass, mass_p, ierr)
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,mass_file,FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(Mass,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
       call VecDestroy(Mass, ierr); CHKERRQ(ierr)
    end if

    if (save_pressure) then
       call VecRestoreArrayF90(Pressure, pressure_p, ierr)
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,pressure_file,FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(Pressure,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
       call VecDestroy(Pressure, ierr); CHKERRQ(ierr)
    end if

    if (save_final_soln) then
       string = 'final_soln.bin'
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(vsfm_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    if (save_internal_mass_flux) then
       call VecRestoreArrayF90(Internal_Mass_Flux, internal_mass_flux_p, ierr)
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,internal_mass_flux_file,FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(Internal_Mass_Flux,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
       call VecDestroy(Internal_Mass_Flux, ierr); CHKERRQ(ierr)
    end if

    if (save_boundary_mass_flux) then
       call VecRestoreArrayF90(Boundary_Mass_Flux, boundary_mass_flux_p, ierr)
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,boundary_mass_flux_file,FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(Boundary_Mass_Flux,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
       call VecDestroy(Boundary_Mass_Flux, ierr); CHKERRQ(ierr)
    end if

    if (save_coupling_mass_flux) then
       call VecRestoreArrayF90(Coupling_Mass_Flux, coupling_mass_flux_p, ierr)
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,coupling_mass_flux_file,FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(Coupling_Mass_Flux,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
       call VecDestroy(Coupling_Mass_Flux, ierr); CHKERRQ(ierr)
    end if

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
    call vsfm_mpp%soe%PrintInfo()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties()

    call VSFMSOEUpdateConnections(vsfm_mpp%soe, MPP_VSFM_SNES_CLM)

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

    GE_oak_xylm      = 0
    GE_oak_root      = 0
    GE_oak_pine_xylm = 0
    GE_pine_xylm     = 0
    GE_pine_root     = 0
    GE_soil          = 0

    select case(trim(problem_type))
    case ('oak')
       neqns       = 1
       GE_oak_xylm = 1; goveqn_name(1) = 'Oak xylem:';

    case ('pine')
       neqns        = 1
       GE_pine_xylm = 1; goveqn_name(1) = 'Pine xylem:';

    case ('oak_and_pine')
       neqns            = 1
       GE_oak_pine_xylm = 1; goveqn_name(1) = 'Oak_Pine xylem:';

    case ('oak_spac')
       neqns       = 3
       GE_oak_xylm = 1; goveqn_name(1) = 'Oak xylem:';
       GE_oak_root = 2; goveqn_name(2) = 'Oak root :';
       GE_soil     = 3; goveqn_name(3) = 'Soil     :';

    case ('pine_spac')
       neqns        = 3
       GE_pine_xylm = 1; goveqn_name(1) = 'Pine xylem:';
       GE_pine_root = 2; goveqn_name(2) = 'Pine root :';
       GE_soil      = 3; goveqn_name(3) = 'Soil      :';

    case ('oak_pine_spac')
       neqns        = 5
       GE_oak_xylm  = 1; goveqn_name(1) = 'Oak xylem:';
       GE_oak_root  = 2; goveqn_name(2) = 'Oak root  :';
       GE_pine_xylm = 3; goveqn_name(3) = 'Pine xylem:';
       GE_pine_root = 4; goveqn_name(4) = 'Pine root :';
       GE_soil      = 5; goveqn_name(5) = 'Soil      :';

    case ('emop_spac')
       neqns      = 9
       GE_e_xylm  = 1; goveqn_name(1) = 'E xylem:';
       GE_e_root  = 2; goveqn_name(2) = 'E root :';
       GE_m_xylm  = 3; goveqn_name(3) = 'M xylem:';
       GE_m_root  = 4; goveqn_name(4) = 'M root :';
       GE_o_xylm  = 5; goveqn_name(5) = 'O xylem:';
       GE_o_root  = 6; goveqn_name(6) = 'O root :';
       GE_p_xylm  = 7; goveqn_name(7) = 'P xylem:';
       GE_p_root  = 8; goveqn_name(8) = 'P root :';
       GE_soil    = 9; goveqn_name(9) = 'Soil   :';

    case ('e_spac')
       neqns      = 3
       GE_e_xylm  = 1; goveqn_name(1) = 'E xylem:';
       GE_e_root  = 2; goveqn_name(2) = 'E root :';
       GE_soil    = 3; goveqn_name(3) = 'Soil   :';

    case ('m_spac')
       neqns      = 3
       GE_m_xylm  = 1; goveqn_name(1) = 'M xylem:';
       GE_m_root  = 2; goveqn_name(2) = 'M root :';
       GE_soil    = 3; goveqn_name(3) = 'Soil   :';

    case ('o_spac')
       neqns      = 3
       GE_o_xylm  = 1; goveqn_name(1) = 'O xylem:';
       GE_o_root  = 2; goveqn_name(2) = 'O root :';
       GE_soil    = 3; goveqn_name(3) = 'Soil   :';

    case ('p_spac')
       neqns      = 3
       GE_p_xylm  = 1; goveqn_name(1) = 'P xylem:';
       GE_p_root  = 2; goveqn_name(2) = 'P root :';
       GE_soil    = 3; goveqn_name(3) = 'Soil   :';

    case ('soil')
       neqns       = 1
       GE_soil     = 1; goveqn_name(1) = 'Soil:';

    case ('e')
       neqns      = 1
       GE_e_xylm  = 1; goveqn_name(1) = 'E xylem:';

    case ('m')
       neqns      = 1
       GE_m_xylm  = 1; goveqn_name(1) = 'M xylem:';

    case ('o')
       neqns      = 1
       GE_o_xylm  = 1; goveqn_name(1) = 'O xylem:';

    case ('p')
       neqns      = 1
       GE_p_xylm  = 1; goveqn_name(1) = 'P xylem:';

    case default
       write(*,*)'Error while adding mesh. Unsupported problem_type : ',trim(problem_type)
       stop

    end select

    call vsfm_mpp%SetNumMeshes(neqns)

    if (GE_oak_xylm      > 0) call add_xylem_mesh_for_single_tree(GE_oak_xylm, oak_nz, oak_Asapwood)
    if (GE_oak_root      > 0) call add_root_mesh_for_single_tree( GE_oak_root, 'oak')
    if (GE_pine_xylm     > 0) call add_xylem_mesh_for_single_tree(GE_pine_xylm, pine_nz, pine_Asapwood)
    if (GE_pine_root     > 0) call add_root_mesh_for_single_tree( GE_pine_root, 'pine')

    if (GE_e_xylm        > 0) call add_xylem_mesh_for_single_tree(GE_e_xylm, nz(E_IDX), Asapwood(E_IDX))
    if (GE_e_root        > 0) call add_root_mesh_for_single_tree( GE_e_root, 'e')
    if (GE_m_xylm        > 0) call add_xylem_mesh_for_single_tree(GE_m_xylm, nz(M_IDX), Asapwood(M_IDX))
    if (GE_m_root        > 0) call add_root_mesh_for_single_tree( GE_m_root, 'm')
    if (GE_o_xylm        > 0) call add_xylem_mesh_for_single_tree(GE_o_xylm, nz(O_IDX), Asapwood(O_IDX))
    if (GE_o_root        > 0) call add_root_mesh_for_single_tree( GE_o_root, 'o')
    if (GE_p_xylm        > 0) call add_xylem_mesh_for_single_tree(GE_p_xylm, nz(P_IDX), Asapwood(P_IDX))
    if (GE_p_root        > 0) call add_root_mesh_for_single_tree( GE_p_root, 'p')

    if (GE_soil          > 0) call add_soil_mesh(GE_soil)
    if (GE_oak_pine_xylm > 0) call add_xylem_mesh_for_two_trees()

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
    PetscReal , pointer :: dz(:)                ! layer thickness of grid cell [m]
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
    allocate(dz     (ncells_xylem ))
    allocate(area_xylem   (ncells_xylem ))
    allocate(filter_xylem (ncells_xylem ))
    allocate(vol_xylem    (ncells_xylem ))
    
    call set_xylem_geometric_attributes( &
         1, local_nz, local_Asapwood,   &
         filter_xylem, area_xylem, vol_xylem, &
         dx_xylem, dy_xylem, dz, &
         xc_xylem, yc_xylem, zc_xylem)

    call mpp_varpar_set_nlevsoi(local_nz)
    call mpp_varpar_set_nlevgrnd(local_nz)

    ! Xylem Mesh
    allocate(mesh)

    call SetMeshPreliminaries(mesh, 'Xylem mesh', ncells_xylem, local_nz)

    call SetMeshGeometricAttributes( &
         mesh, filter_xylem, &
         xc_xylem, yc_xylem, zc_xylem, &
         dx_xylem, dy_xylem, dz, &
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
       vert_conn_dist_up_xylem(iconn) = 0.5d0*dz_xylem
       vert_conn_dist_dn_xylem(iconn) = 0.5d0*dz_xylem
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
    deallocate (dz                )
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
    PetscReal , pointer :: dz(:)                ! layer thickness of grid cell [m]
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
    allocate(dz     (ncells_xylem ))
    allocate(area_xylem   (ncells_xylem ))
    allocate(filter_xylem (ncells_xylem ))
    allocate(vol_xylem    (ncells_xylem ))

    ibeg = 1; iend = oak_nz; Asapwood = oak_Asapwood
    call set_xylem_geometric_attributes( &
         ibeg, iend, Asapwood, &
         filter_xylem, area_xylem, vol_xylem, &
         dx_xylem, dy_xylem, dz, &
         xc_xylem, yc_xylem, zc_xylem)

    ibeg = 1 + oak_nz; iend = oak_nz + pine_nz; Asapwood = pine_Asapwood
    call set_xylem_geometric_attributes( &
         ibeg, iend, Asapwood, &
         filter_xylem, area_xylem, vol_xylem, &
         dx_xylem, dy_xylem, dz, &
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
         dx_xylem, dy_xylem, dz, &
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
       vert_conn_dist_up_xylem(iconn) = 0.5d0*dz_xylem
       vert_conn_dist_dn_xylem(iconn) = 0.5d0*dz_xylem
       vert_conn_area_xylem(iconn)    = area_xylem(kk)
       vert_conn_type_xylem(iconn)    = CONN_VERTICAL
    end do

    do kk = oak_nz + 1, oak_nz + pine_nz - 1
       iconn = iconn + 1
       vert_conn_id_up_xylem(iconn)   = kk 
       vert_conn_id_dn_xylem(iconn)   = kk + 1
       vert_conn_dist_up_xylem(iconn) = 0.5d0*dz_xylem
       vert_conn_dist_dn_xylem(iconn) = 0.5d0*dz_xylem
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
    deallocate (dz                )
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
       dx_xylem, dy_xylem, dz, &
       xc_xylem, yc_xylem, zc_xylem)

    implicit none

    PetscInt            :: ibeg, iend
    PetscReal           :: Asapwood
    PetscReal , pointer :: xc_xylem(:)
    PetscReal , pointer :: yc_xylem(:)
    PetscReal , pointer :: zc_xylem(:)
    PetscReal , pointer :: dx_xylem(:)
    PetscReal , pointer :: dy_xylem(:)
    PetscReal , pointer :: dz(:)
    PetscReal , pointer :: area_xylem(:)
    PetscReal , pointer :: vol_xylem(:)
    PetscInt  , pointer :: filter_xylem(:)

    PetscInt            :: kk, nz

    filter_xylem (ibeg:iend) = 1
    area_xylem   (ibeg:iend) = Asapwood
    dx_xylem     (ibeg:iend) = sqrt(Asapwood)
    dy_xylem     (ibeg:iend) = sqrt(Asapwood)
    dz           (ibeg:iend) = dz_xylem
    xc_xylem     (ibeg:iend) = 0.d0!sqrt(Asapwood)/2.d0
    yc_xylem     (ibeg:iend) = 0.d0!sqrt(Asapwood)/2.d0

    nz = iend - ibeg + 1

    zc_xylem(ibeg)  = -0.17d0*0.d0 + nz * dz_xylem
    vol_xylem(ibeg) = area_xylem(ibeg) * dz_xylem
    do kk = ibeg+1, iend
       zc_xylem(kk)  = -(dz_xylem/2.d0 + dz_xylem * (kk - ibeg)) + nz * dz_xylem
       vol_xylem(kk) = area_xylem(kk) * dz_xylem
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
    PetscInt                   :: IDX
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
    PetscReal                  :: local_Asapwood

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

    case ('e')
       IDX    = E_IDX
       nz     = root_nz(IDX)
       qz     = root_qz(IDX)
       dd     = root_d(IDX)
       rld0   = root_rld0(IDX)
       radius = root_radius(IDX)
       local_Asapwood = Asapwood(IDX)

    case ('m')
       IDX    = M_IDX
       nz     = root_nz(IDX)
       qz     = root_qz(IDX)
       dd     = root_d(IDX)
       rld0   = root_rld0(IDX)
       radius = root_radius(IDX)
       local_Asapwood = Asapwood(IDX)

    case ('o')
       IDX    = O_IDX
       nz     = root_nz(IDX)
       qz     = root_qz(IDX)
       dd     = root_d(IDX)
       rld0   = root_rld0(IDX)
       radius = root_radius(IDX)
       local_Asapwood = Asapwood(IDX)

    case ('p')
       IDX    = P_IDX
       nz     = root_nz(IDX)
       qz     = root_qz(IDX)
       dd     = root_d(IDX)
       rld0   = root_rld0(IDX)
       radius = root_radius(IDX)
       local_Asapwood = Asapwood(IDX)

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

       zc_1d(kk)    = -(kk-1)*dz_soil - dz_soil/2.0

       root_len_den = rld0 * (1.d0 - abs(zc_1d(kk))/dd)*exp(-qz*abs(zc_1d(kk))/dd)
       root_len_den = rld0 * exp(-qz*abs(zc_1d(kk))/dd)
       root_len     = root_len_den * & ! [m/m^3]
                      (dx*dy*dz_soil)       ! [m^3]

       vol(kk)      = PI*(radius**2.d0)*root_len
       area(kk)     = 2.d0*PI*radius*root_len

       xc_1d(kk)    = 0.d0!radius
       yc_1d(kk)    = 0.d0!root_len

       dx_1d(kk)    = dx
       dy_1d(kk)    = dy
       dz_1d(kk)    = dz_soil

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
       root_len_den     = rld0 * exp(-qz*abs(zc_1d(kk))/dd)
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
       conn_id_up(kk)   = 0
       conn_id_dn(kk)   = kk
       conn_dist_up(kk) = 0.d0
       conn_dist_dn(kk) = dz_soil/2.d0
       conn_area(kk)    = local_Asapwood!PI*(radius**2.d0)
    end do

    call mesh%CreateAndAddConnectionSet( &
         CONN_SET_CONDITIONS,            &
         nconn,  conn_id_up, conn_id_dn, &
         conn_dist_up, conn_dist_dn,     &
         conn_area, conn_type, conn_unit_vec)

    ! Add internal connection set
    nconn = nz-1
    do kk = 1, nz-1
       conn_id_up(kk)   = kk
       conn_id_dn(kk)   = kk+1

       conn_dist_up(kk) = dz_soil/2.d0
       conn_dist_dn(kk) = dz_soil/2.d0
       conn_area(kk)    = local_Asapwood!PI*(radius**2.d0)
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
    PetscInt :: kk, ncells,nz
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
    z_min = -soil_nz*dz_soil
    nz    = soil_nz

    call ComputeXYZCentroids( nx, ny, nz, dx, dy, dz_soil, &
         x_min, y_min, z_min, xc_1d, yc_1d, zc_1d)

    do kk = 1, nz
       zc_1d(kk)    = -(kk-1)*dz_soil - dz_soil/2.0
    end do

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
    dz_1d  (:) = dz_soil
    area   (:) = dy*dz_soil
    volume (:) = dx*dy*dz_soil

    allocate(mesh)

    call SetMeshPreliminaries(mesh, 'Soil', ncells, nz)

    call SetMeshGeometricAttributes( &
         mesh, filter, &
         xc_1d, yc_1d, zc_1d, &
         dx_1d, dy_1d, dz_1d, &
         area, volume)

    call ComputeInternalConnections(nx, ny, nz, dx, dy, dz_soil, CONN_IN_Z_DIR, &
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
       call vsfm_mpp%AddGovEqnWithMeshRank(GE_RE, trim(goveqn_name(irank))//'Richards Equation ODE', irank)
    end do

    call vsfm_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqn

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()

    implicit none

    select case(trim(problem_type))
    case ('oak')
       call add_bc_to_goveqns_for_single_tree(GE_oak_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_oak_xylm)

    case ('pine')
       call add_bc_to_goveqns_for_single_tree(GE_pine_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_pine_xylm)

    case ('oak_and_pine')
       call add_bc_to_goveqns_for_two_trees(GE_oak_pine_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_oak_pine_xylm)

    case ('oak_spac')
       call add_ss_to_goveqns_for_single_tree(GE_oak_xylm)
       call add_ss_for_soil(GE_soil)

    case ('pine_spac')
       call add_ss_to_goveqns_for_single_tree(GE_pine_xylm)
       call add_ss_for_soil(GE_soil)

    case ('oak_pine_spac')
       call add_ss_to_goveqns_for_single_tree(GE_oak_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_pine_xylm)
       call add_ss_for_soil(GE_soil)

    case ('soil')
       if (soil_ss_specified) call add_ss_for_soil(GE_soil)
       if (soil_bc_specified) call add_top_bc_to_soil_goveqn(GE_soil)

    case ('emop_spac')
       call add_ss_to_goveqns_for_single_tree(GE_e_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_m_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_o_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_p_xylm)
       if (soil_ss_specified) call add_ss_for_soil(GE_soil)
       if (soil_bc_specified .or. sm_bc_specified) call add_top_bc_to_soil_goveqn(GE_soil)
       call add_bottom_bc_to_soil_goveqn(GE_SOIL)

    case ('e_spac')
       call add_ss_to_goveqns_for_single_tree(GE_e_xylm)
       if (soil_ss_specified) call add_ss_for_soil(GE_soil)
       if (soil_bc_specified .or. sm_bc_specified) call add_top_bc_to_soil_goveqn(GE_soil)
       call add_bottom_bc_to_soil_goveqn(GE_SOIL)

    case ('m_spac')
       call add_ss_to_goveqns_for_single_tree(GE_m_xylm)
       if (soil_ss_specified) call add_ss_for_soil(GE_soil)
       if (soil_bc_specified .or. sm_bc_specified) call add_top_bc_to_soil_goveqn(GE_soil)
       call add_bottom_bc_to_soil_goveqn(GE_SOIL)

    case ('o_spac')
       call add_ss_to_goveqns_for_single_tree(GE_o_xylm)
       if (soil_ss_specified) call add_ss_for_soil(GE_soil)
       if (soil_bc_specified .or. sm_bc_specified) call add_top_bc_to_soil_goveqn(GE_soil)
       call add_bottom_bc_to_soil_goveqn(GE_SOIL)

    case ('p_spac')
       call add_ss_to_goveqns_for_single_tree(GE_p_xylm)
       if (soil_ss_specified) call add_ss_for_soil(GE_soil)
       if (soil_bc_specified .or. sm_bc_specified) call add_top_bc_to_soil_goveqn(GE_soil)
       call add_bottom_bc_to_soil_goveqn(GE_SOIL)

    case ('e')
       call add_bc_to_goveqns_for_single_tree(GE_e_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_e_xylm)

    case ('m')
       call add_bc_to_goveqns_for_single_tree(GE_m_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_m_xylm)

    case ('o')
       call add_bc_to_goveqns_for_single_tree(GE_o_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_o_xylm)

    case ('p')
       call add_bc_to_goveqns_for_single_tree(GE_p_xylm)
       call add_ss_to_goveqns_for_single_tree(GE_p_xylm)

    case default
       write(*,*)'Error while adding condition. Unsupported problem_type : ',trim(problem_type)
       stop

    end select

    if (GE_oak_xylm  > 0 .and. GE_oak_root  > 0) call add_xylm2root_coupling_bc('oak')
    if (GE_pine_xylm > 0 .and. GE_pine_root > 0) call add_xylm2root_coupling_bc('pine')
    if (GE_e_xylm    > 0 .and. GE_e_root    > 0) call add_xylm2root_coupling_bc('e')
    if (GE_m_xylm    > 0 .and. GE_m_root    > 0) call add_xylm2root_coupling_bc('m')
    if (GE_o_xylm    > 0 .and. GE_o_root    > 0) call add_xylm2root_coupling_bc('o')
    if (GE_p_xylm    > 0 .and. GE_p_root    > 0) call add_xylm2root_coupling_bc('p')
    if (GE_oak_root  > 0 .and. GE_soil      > 0) call add_root2soil_coupling_bc('oak')
    if (GE_pine_root > 0 .and. GE_soil      > 0) call add_root2soil_coupling_bc('pine')
    if (GE_e_root    > 0 .and. GE_soil      > 0) call add_root2soil_coupling_bc('e')
    if (GE_m_root    > 0 .and. GE_soil      > 0) call add_root2soil_coupling_bc('m')
    if (GE_o_root    > 0 .and. GE_soil      > 0) call add_root2soil_coupling_bc('o')
    if (GE_p_root    > 0 .and. GE_soil      > 0) call add_root2soil_coupling_bc('p')

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine add_bc_to_goveqns_for_single_tree(ieqn)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt, intent(in) :: ieqn

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Bottom BC', 'Pa', COND_DIRICHLET, SOIL_BOTTOM_CELLS)

  end subroutine add_bc_to_goveqns_for_single_tree

  !------------------------------------------------------------------------
  subroutine add_top_bc_to_soil_goveqn(ieqn)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt, intent(in) :: ieqn

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Top BC', 'Pa', COND_DIRICHLET, SOIL_TOP_CELLS)

  end subroutine add_top_bc_to_soil_goveqn

  !------------------------------------------------------------------------
  subroutine add_bottom_bc_to_soil_goveqn(ieqn)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt, intent(in) :: ieqn

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Bottom BC', 'Pa', COND_DIRICHLET, SOIL_BOTTOM_CELLS)

  end subroutine add_bottom_bc_to_soil_goveqn

  !------------------------------------------------------------------------
  subroutine add_ss_to_goveqns_for_single_tree(ieqn)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DOWNREG_MASS_RATE_FETCH2
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt, intent(in) :: ieqn

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Potential Mass_Flux', 'kg/s', COND_DOWNREG_MASS_RATE_FETCH2, &
         ALL_CELLS)

  end subroutine add_ss_to_goveqns_for_single_tree

  !------------------------------------------------------------------------
  subroutine add_bc_to_goveqns_for_two_trees(ieqn)
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
    conn_dist_dn(iconn)   = 0.5d0*dz_soil
    conn_area(iconn)      = oak_Asapwood
    conn_type(iconn)      = CONN_VERTICAL
    conn_unitvec(iconn,1) = 0.d0
    conn_unitvec(iconn,2) = 0.d0
    conn_unitvec(iconn,3) = 1.d0

    iconn = 2
    conn_id_up(iconn)     = 0
    conn_id_dn(iconn)     = oak_nz+pine_nz
    conn_dist_up(iconn)   = 0.d0
    conn_dist_dn(iconn)   = 0.5d0*dz_soil
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

  end subroutine add_bc_to_goveqns_for_two_trees

  !------------------------------------------------------------------------
  subroutine add_ss_for_soil(ieqn)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use petscsys
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt, intent(in) :: ieqn

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS,   &
         'Infiltration', 'kg/s', COND_MASS_RATE, &
         SOIL_TOP_CELLS)

  end subroutine add_ss_for_soil

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

    !
    ! Allocate auxvars
    !
    if (GE_oak_root > 0 .and. GE_soil > 0) then
       call SetCouplingAuxVarsInXylemAndRootEqns(GE_oak_xylm, GE_oak_root)
    end if

    if (GE_pine_root > 0 .and. GE_soil > 0) then
       call SetCouplingAuxVarsInXylemAndRootEqns(GE_pine_xylm, GE_pine_root)
    end if

    if (GE_e_root > 0 .and. GE_soil > 0) then
       call SetCouplingAuxVarsInXylemAndRootEqns(GE_e_xylm, GE_e_root)
    end if

    if (GE_m_root > 0 .and. GE_soil > 0) then
       call SetCouplingAuxVarsInXylemAndRootEqns(GE_m_xylm, GE_m_root)
    end if

    if (GE_o_root > 0 .and. GE_soil > 0) then
       call SetCouplingAuxVarsInXylemAndRootEqns(GE_o_xylm, GE_o_root)
    end if

    if (GE_p_root > 0 .and. GE_soil > 0) then
       call SetCouplingAuxVarsInXylemAndRootEqns(GE_p_xylm, GE_p_root)
    end if

    if (GE_soil > 0) then

       ! Allocate maximum possible memory
       allocate (var_ids_for_coupling    (neqns))
       allocate (goveqn_ids_for_coupling (neqns))

       var_ids_for_coupling (:) = VAR_PRESSURE

       nvars_for_coupling = 0
       if (GE_oak_root > 0) then
          nvars_for_coupling = nvars_for_coupling + 1
          goveqn_ids_for_coupling (nvars_for_coupling) = GE_oak_root
       end if

       if (GE_pine_root > 0) then
          nvars_for_coupling = nvars_for_coupling + 1
          goveqn_ids_for_coupling (nvars_for_coupling) = GE_pine_root
       end if

       if (GE_e_root > 0) then
          nvars_for_coupling = nvars_for_coupling + 1
          goveqn_ids_for_coupling (nvars_for_coupling) = GE_e_root
       end if

       if (GE_m_root > 0) then
          nvars_for_coupling = nvars_for_coupling + 1
          goveqn_ids_for_coupling (nvars_for_coupling) = GE_m_root
       end if

       if (GE_o_root > 0) then
          nvars_for_coupling = nvars_for_coupling + 1
          goveqn_ids_for_coupling (nvars_for_coupling) = GE_o_root
       end if

       if (GE_p_root > 0) then
          nvars_for_coupling = nvars_for_coupling + 1
          goveqn_ids_for_coupling (nvars_for_coupling) = GE_p_root
       end if

       ieqn = GE_soil

       call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
            var_ids_for_coupling, goveqn_ids_for_coupling)

       deallocate(var_ids_for_coupling   )
       deallocate(goveqn_ids_for_coupling)
       
    end if

    call vsfm_mpp%AllocateAuxVars()

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------

  subroutine SetCouplingAuxVarsInXylemAndRootEqns(GE_xylm, GE_root)

    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants, only : VAR_PRESSURE

    implicit none

    integer           :: GE_xylm, GE_root

    integer           :: ieqn
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)

    nvars_for_coupling = 1
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    ieqn = GE_xylm
    var_ids_for_coupling    (1) = VAR_PRESSURE
    goveqn_ids_for_coupling (1) = GE_root

    call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

    nvars_for_coupling = 2
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    ieqn = GE_root
    var_ids_for_coupling    (1) = VAR_PRESSURE
    var_ids_for_coupling    (2) = VAR_PRESSURE
    goveqn_ids_for_coupling (1) = GE_xylm
    goveqn_ids_for_coupling (2) = GE_soil

    call vsfm_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

  end subroutine SetCouplingAuxVarsInXylemAndRootEqns

  !------------------------------------------------------------------------
  subroutine set_material_properties()

    implicit none

    PetscInt :: ii, ieqn, IDX
    
    select case(trim(problem_type))
    case ('oak')
       call set_material_properties_for_single_tree(   &
            GE_oak_xylm, &
            oak_nz, porosity, oak_phi88, oak_phi50,    &
            oak_kmax, oak_c1, oak_c2, oak_c3, oak_phis50)

    case ('pine')
       call set_material_properties_for_single_tree(   &
            GE_pine_xylm, &
            pine_nz, porosity, pine_phi88, pine_phi50, &
            pine_kmax, pine_c1, pine_c2, pine_c3, pine_phis50)

    case ('oak_and_pine')
       call set_material_properties_for_two_trees(GE_oak_pine_xylm)

    case ('oak_spac')
       call set_material_properties_for_single_tree(   &
            GE_oak_xylm, &
            oak_nz, porosity, oak_phi88, oak_phi50,    &
            oak_kmax, oak_c1, oak_c2, oak_c3, oak_phis50)

       call set_material_properties_for_single_tree(   &
            GE_oak_root, oak_root_nz, porosity, &
            oak_axi_root_phi88, oak_axi_root_phi50, oak_axi_root_kmax, &
            oak_axi_root_c1, oak_axi_root_c2, oak_c3, &
            oak_phis50, set_ss_values=.false.)

       call set_material_properties_for_soil(GE_soil)

       call set_material_properties_for_soil_bc()
       call set_material_properties_for_root_bc('oak')

    case ('pine_spac')
       call set_material_properties_for_single_tree(   &
            GE_pine_xylm, &
            pine_nz, porosity, pine_phi88, pine_phi50, &
            pine_kmax, pine_c1, pine_c2, pine_c3, pine_phis50)

       call set_material_properties_for_single_tree(   &
            GE_pine_root, pine_root_nz, porosity, &
            pine_axi_root_phi88, pine_axi_root_phi50, pine_axi_root_kmax, &
            pine_axi_root_c1, pine_axi_root_c2, pine_c3, &
            pine_phis50, set_ss_values=.false.)

       call set_material_properties_for_soil(GE_soil)

       call set_material_properties_for_soil_bc()
       call set_material_properties_for_root_bc('pine')

    case ('oak_pine_spac')
       call set_material_properties_for_single_tree(   &
            GE_oak_xylm, &
            oak_nz, porosity, oak_phi88, oak_phi50,    &
            oak_kmax, oak_c1, oak_c2, oak_c3, oak_phis50)

       call set_material_properties_for_single_tree(   &
            GE_oak_root, oak_root_nz, porosity, &
            oak_axi_root_phi88, oak_axi_root_phi50, oak_axi_root_kmax, &
            oak_axi_root_c1, oak_axi_root_c2, oak_c3, &
            oak_phis50, set_ss_values=.false.)

       call set_material_properties_for_single_tree(   &
            GE_pine_xylm, &
            pine_nz, porosity, pine_phi88, pine_phi50, &
            pine_kmax, pine_c1, pine_c2, pine_c3, pine_phis50)

       call set_material_properties_for_single_tree(   &
            GE_pine_root, pine_root_nz, porosity, &
            pine_axi_root_phi88, pine_axi_root_phi50, pine_axi_root_kmax, &
            pine_axi_root_c1, pine_axi_root_c2, pine_c3, &
            pine_phis50, set_ss_values=.false.)

       call set_material_properties_for_soil(GE_soil)

       call set_material_properties_for_soil_bc()
       call set_material_properties_for_root_bc('oak')
       call set_material_properties_for_root_bc('pine')

    case ('soil')
       call set_material_properties_for_soil(GE_soil)
       if (soil_bc_specified) call set_material_properties_for_soil_bc()

    case ('emop_spac','e_spac','m_spac','o_spac','p_spac','e','m','o','p')
       do ii = 1,4
          select case(ii)
          case (1)
             IDX = E_IDX; ieqn = GE_e_xylm
          case (2)
             IDX = M_IDX; ieqn = GE_m_xylm
          case (3)
             IDX = O_IDX; ieqn = GE_o_xylm
          case (4)
             IDX = P_IDX; ieqn = GE_p_xylm
          end select
             
          if (ieqn == 0) cycle

          call set_material_properties_for_single_tree(   &
               ieqn, &
               nz(IDX), &
               xylem_porosity(IDX), &
               phi88(IDX), &
               phi50(IDX),    &
               kmax(IDX), &
               !axi_root_kmax(IDX), &
               c1(IDX), &
               c2(IDX), &
               c3(IDX), &
               phis50(IDX) &
               )
       end do

       do ii = 1,4
          select case(ii)
          case (1)
             IDX = E_IDX; ieqn = GE_e_root
          case (2)
             IDX = M_IDX; ieqn = GE_m_root
          case (3)
             IDX = O_IDX; ieqn = GE_o_root
          case (4)
             IDX = P_IDX; ieqn = GE_p_root
          end select
             
          if (ieqn == 0) cycle

          call set_material_properties_for_single_tree(   &
               ieqn, &
               root_nz(IDX), &
               root_porosity(IDX), &
               axi_root_phi88(IDX), &
               axi_root_phi50(IDX), &
               axi_root_kmax(IDX), &
               axi_root_c1(IDX), &
               axi_root_c2(IDX), &
               c3(IDX), &
               phis50(IDX), &
               set_ss_values=.false.)
       end do

       if (GE_soil > 0) call set_material_properties_for_soil(GE_soil)
       if (GE_soil > 0) then
          call set_material_properties_for_soil_bc()
       end if

       if (GE_e_root > 0) call set_material_properties_for_root_bc('e')
       if (GE_m_root > 0) call set_material_properties_for_root_bc('m')
       if (GE_o_root > 0) call set_material_properties_for_root_bc('o')
       if (GE_p_root > 0) call set_material_properties_for_root_bc('p')

    case default
       write(*,*)'Unable to set material prop for problem_type = ' // trim(problem_type)
       stop

    end select

  end subroutine set_material_properties
  !------------------------------------------------------------------------
  subroutine set_material_properties_for_single_tree(igoveqn, local_nz, local_porosity, &
       local_phi88, local_phi50, local_kmax, local_c1, local_c2, local_c3, local_phis50, &
       set_ss_values)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetDensityType
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSourceSinkAuxVarRealValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_EXPONENT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use EOSWaterMod          , only : DENSITY_CONSTANT
    !
    implicit none
    !
    PetscInt  :: igoveqn
    PetscInt  :: local_nz
    PetscReal :: local_porosity
    PetscReal :: local_phi88
    PetscReal :: local_phi50
    PetscReal :: local_kmax
    PetscReal :: local_c1
    PetscReal :: local_c2
    PetscReal :: local_c3
    PetscReal :: local_phis50
    PetscBool,optional :: set_ss_values
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
    PetscBool             :: set_ss
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

    set_ss = .true.
    if (present(set_ss_values)) set_ss = set_ss_values

    call VSFMMPPSetDensityType(vsfm_mpp, igoveqn, DENSITY_TGDPB01)
    !call VSFMMPPSetDensityType(vsfm_mpp, igoveqn, DENSITY_CONSTANT)

    call set_xylem_material_properties(1, local_nz,       &
         local_porosity, local_phi88, local_phi50,        &
         local_kmax, local_c1, local_c2, local_c3,        &
         local_phis50,                                    &
         por, satfunc_type, alpha, lambda, sat_res, perm, &
         relperm_type, weibull_d, weibull_c, ss_exponent, &
         ss_pressure)

    call VSFMMPPSetSoilPorosity(vsfm_mpp, igoveqn, por)

    call VSFMMPPSetSaturationFunction(vsfm_mpp, igoveqn, satfunc_type, &
         alpha, lambda, sat_res)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, igoveqn, perm, perm, perm)

    call VSFMMPPSetRelativePermeability(vsfm_mpp, igoveqn, &
         relperm_type, weibull_d, weibull_c)

    if (set_ss) then
       call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, igoveqn, &
            VAR_POT_MASS_SINK_EXPONENT, ss_exponent)

       call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, igoveqn, &
            VAR_POT_MASS_SINK_PRESSURE, ss_pressure)
    end if

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
  subroutine set_material_properties_for_two_trees(igoveqn)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetDensityType
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSourceSinkAuxVarRealValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeability
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_POT_MASS_SINK_EXPONENT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use EOSWaterMod          , only : DENSITY_CONSTANT
    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL
    !
    implicit none
    !
    PetscInt              :: igoveqn
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

    call VSFMMPPSetDensityType(vsfm_mpp, igoveqn, DENSITY_TGDPB01)
    !call VSFMMPPSetDensityType(vsfm_mpp, igoveqn, DENSITY_CONSTANT)

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

    call VSFMMPPSetSoilPorosity(vsfm_mpp, igoveqn, por)

    call VSFMMPPSetSaturationFunction(vsfm_mpp, igoveqn, satfunc_type, &
         alpha, lambda, sat_res)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, igoveqn, perm, perm, perm)
    
    call VSFMMPPSetRelativePermeability(vsfm_mpp, igoveqn, &
         relperm_type, weibull_d, weibull_c)

    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, igoveqn, &
         VAR_POT_MASS_SINK_EXPONENT, ss_exponent)

    call VSFMMPPSetSourceSinkAuxVarRealValue(vsfm_mpp, igoveqn, &
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
  subroutine set_material_properties_for_soil(igoveqn)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbVSFM , only : VSFMMPPSetDensityType
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM , only : VSFMMPPSetRelativePermeability
    use SaturationFunction   , only : SAT_FUNC_VAN_GENUCHTEN
    use EOSWaterMod          , only : DENSITY_TGDPB01
    use EOSWaterMod          , only : DENSITY_CONSTANT
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt           :: igoveqn
    !
    PetscInt           :: ncells
    PetscReal, pointer :: por(:)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: lambda(:)
    PetscReal, pointer :: residual_sat(:)
    PetscReal, pointer :: perm(:)
    PetscInt , pointer :: satfunc_type(:)
    PetscInt :: k

    call VSFMMPPSetDensityType(vsfm_mpp, igoveqn, DENSITY_TGDPB01)
    !call VSFMMPPSetDensityType(vsfm_mpp, igoveqn, DENSITY_CONSTANT)

    ncells = soil_nz
    
    allocate(por             (ncells))
    allocate(perm            (ncells))
    allocate(alpha           (ncells))
    allocate(lambda          (ncells))
    allocate(residual_sat    (ncells))
    allocate(satfunc_type    (ncells))
    
    por(:)          = soil_por
    perm(:)         = soil_perm
    do k=10,ncells
       por(k)   = soil_por_bot
       perm(k)  = soil_perm_bot
    end do
    alpha(:)        = soil_alpha
    lambda(:)       = soil_vg_m
    residual_sat(:) = soil_sat_res
    satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN

    call VSFMMPPSetSoilPorosity(vsfm_mpp, igoveqn, por)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, igoveqn, perm, perm, perm)

    call VSFMMPPSetSaturationFunction(vsfm_mpp, igoveqn, satfunc_type, &
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
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetRelativePermeabilityAuxVarConn
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnIntValue
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetAuxVarConnRealValue
    use SaturationFunction        , only : SAT_FUNC_VAN_GENUCHTEN
    use SaturationFunction        , only : SAT_FUNC_FETCH2
    use SaturationFunction        , only : RELPERM_FUNC_WEIBULL, RELPERM_FUNC_MUALEM
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
    PetscInt           :: nconns
    PetscInt , pointer :: flux_type          (:)
    PetscInt , pointer :: conductance_type   (:)
    PetscBool, pointer :: set_upwind_auxvar  (:)
    PetscInt , pointer :: dn_satfunc_type  (:)
    PetscReal, pointer :: dn_alpha         (:)
    PetscReal, pointer :: dn_lambda        (:)
    PetscReal, pointer :: dn_residual_sat  (:)
    PetscReal, pointer :: dn_conductance   (:)
    PetscReal, pointer :: dn_weibull_c_val (:)
    PetscReal, pointer :: dn_weibull_d_val (:)
    PetscInt , pointer :: dn_relperm_type  (:)
    PetscInt :: kk, nz, GE, iroot
    PetscReal :: zc_1d, qz, dd, rld0, root_len_den, droot
    PetscInt :: nroot, IDX
    PetscInt, pointer :: GE_root_id(:)
    PetscReal :: c1, c2, kmax, phi88, phi50
    
    eqn_id = GE_soil

    allocate(GE_root_id(neqns))

    nconns = 0
    nroot  = 0
    if (GE_oak_root  > 0) then
       nconns            = nconns + oak_root_nz
       nroot             = nroot + 1
       GE_root_id(nroot) = GE_oak_root
    end if
    if (GE_pine_root > 0) then
       nconns            = nconns + pine_root_nz
       nroot             = nroot + 1
       GE_root_id(nroot) = GE_pine_root
    end if
    if (GE_e_root    > 0) then
       nconns            = nconns + root_nz(E_IDX)
       nroot             = nroot + 1
       GE_root_id(nroot) = GE_e_root
    end if
    if (GE_m_root    > 0) then
       nconns            = nconns + root_nz(M_IDX)
       nroot             = nroot + 1
       GE_root_id(nroot) = GE_m_root
    end if
    if (GE_o_root    > 0) then
       nconns            = nconns + root_nz(O_IDX)
       nroot             = nroot + 1
       GE_root_id(nroot) = GE_o_root
    end if
    if (GE_p_root    > 0) then
       nconns            = nconns + root_nz(P_IDX)
       nroot             = nroot + 1
       GE_root_id(nroot) = GE_p_root
    end if

    if (soil_bc_specified .or. sm_bc_specified) nconns = nconns + 1
    nconns = nconns + 1

    allocate(flux_type          (nconns))
    allocate(conductance_type   (nconns))
    allocate(set_upwind_auxvar  (nconns))
    allocate(dn_alpha           (nconns))
    allocate(dn_lambda          (nconns))
    allocate(dn_residual_sat    (nconns))
    allocate(dn_conductance     (nconns))
    allocate(dn_satfunc_type    (nconns))
    allocate(dn_weibull_c_val   (nconns))
    allocate(dn_weibull_d_val   (nconns))
    allocate(dn_relperm_type    (nconns))

    ! Set values for all BCs
    flux_type        (1:nconns) = CONDUCTANCE_FLUX_TYPE

    conductance_type(1:nconns) = CONDUCTANCE_MANOLI_TYPE
    dn_alpha        (1:nconns) = soil_alpha
    dn_lambda       (1:nconns) = soil_vg_m
    dn_residual_sat (1:nconns) = soil_sat_res
    dn_satfunc_type (1:nconns) = SAT_FUNC_VAN_GENUCHTEN

    nconns = 0

    if (soil_bc_specified .or. sm_bc_specified) then
       nconns = nconns + 1;
       flux_type(nconns) = DARCY_FLUX_TYPE;
       conductance_type(nconns) = 0;
       dn_relperm_type(nconns) = RELPERM_FUNC_MUALEM

       nconns = nconns + 1;
       flux_type(nconns) = DARCY_FLUX_TYPE;
       conductance_type(nconns) = 0;
       dn_relperm_type(nconns) = RELPERM_FUNC_MUALEM
    end if

    do iroot = 1, nroot
       if (GE_root_id(iroot) == GE_oak_root) then
          nz   = oak_root_nz
          qz   = oak_root_qz
          rld0 = oak_root_rld0
          dd   = oak_root_d
          GE   = GE_oak_root
          c1   = oak_rad_root_c1
          c2   = oak_rad_root_c2
          kmax = oak_rad_root_kmax
          phi88= oak_rad_root_phi88
          phi50= oak_rad_root_phi50

       elseif (GE_root_id(iroot) == GE_pine_root) then
          nz   = pine_root_nz
          qz   = pine_root_qz
          rld0 = pine_root_rld0
          dd   = pine_root_d
          GE   = GE_pine_root
          c1   = pine_rad_root_c1
          c2   = pine_rad_root_c2
          kmax = pine_rad_root_kmax
          phi88= pine_rad_root_phi88
          phi50= pine_rad_root_phi50

       elseif (GE_root_id(iroot) == GE_e_root) then
          IDX  = E_IDX
          nz   = root_nz(IDX)
          qz   = root_qz(IDX)
          rld0 = root_rld0(IDX)
          dd   = root_d(IDX)
          GE   = GE_root_id(iroot)
          c1   = rad_root_c1(IDX)
          c2   = rad_root_c2(IDX)
          kmax = rad_root_kmax(IDX)
          phi88= rad_root_phi88(IDX)
          phi50= rad_root_phi50(IDX)

       elseif (GE_root_id(iroot) == GE_m_root) then
          IDX  = M_IDX
          nz   = root_nz(IDX)
          qz   = root_qz(IDX)
          rld0 = root_rld0(IDX)
          dd   = root_d(IDX)
          GE   = GE_root_id(iroot)
          c1   = rad_root_c1(IDX)
          c2   = rad_root_c2(IDX)
          kmax = rad_root_kmax(IDX)
          phi88= rad_root_phi88(IDX)
          phi50= rad_root_phi50(IDX)

       elseif (GE_root_id(iroot) == GE_o_root) then
          IDX  = O_IDX
          nz   = root_nz(IDX)
          qz   = root_qz(IDX)
          rld0 = root_rld0(IDX)
          dd   = root_d(IDX)
          GE   = GE_root_id(iroot)
          c1   = rad_root_c1(IDX)
          c2   = rad_root_c2(IDX)
          kmax = rad_root_kmax(IDX)
          phi88= rad_root_phi88(IDX)
          phi50= rad_root_phi50(IDX)

       elseif (GE_root_id(iroot) == GE_p_root) then
          IDX  = P_IDX
          nz   = root_nz(IDX)
          qz   = root_qz(IDX)
          rld0 = root_rld0(IDX)
          dd   = root_d(IDX)
          GE   = GE_root_id(iroot)
          c1   = rad_root_c1(IDX)
          c2   = rad_root_c2(IDX)
          kmax = rad_root_kmax(IDX)
          phi88= rad_root_phi88(IDX)
          phi50= rad_root_phi50(IDX)
       end if

       
       if (GE > 0) then
          do kk = 1, nz
             nconns = nconns + 1

             zc_1d = -(kk-1)*dz_soil - dz_soil/2.0

             root_len_den = rld0 * (1.d0 - abs(zc_1d)/dd)*exp(-qz*abs(zc_1d)/dd)
             root_len_den = rld0 * exp(-qz*abs(zc_1d)/dd)

             droot = (PI*root_len_den)**(-0.5d0)

             if (kk<5) then
                dn_conductance(nconns) = soil_perm/droot
             else
                dn_conductance(nconns) = soil_perm_bot/droot
             end if
             dn_conductance(nconns) = (kmax*vis/rho)/droot

             dn_alpha        (nconns) = phi88
             dn_lambda       (nconns) = phi50
             dn_residual_sat (nconns) = 0.d0
             dn_satfunc_type (nconns) = SAT_FUNC_FETCH2
             dn_weibull_d_val(nconns) = c1
             dn_weibull_c_val(nconns) = c2
             dn_relperm_type (nconns) = RELPERM_FUNC_WEIBULL
          end do
       end if
    end do

    ! Set downwind values
    set_upwind_auxvar(:) = PETSC_FALSE
    call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp , &
         eqn_id, AUXVAR_CONN_BC                    , &
         set_upwind_auxvar, dn_satfunc_type, dn_alpha    , &
         dn_lambda, dn_residual_sat)
        
    call VSFMMPPSetRelativePermeabilityAuxVarConn(vsfm_mpp , &
         eqn_id, AUXVAR_CONN_BC                    , &
         set_upwind_auxvar, dn_relperm_type, dn_weibull_d_val, &
         dn_weibull_c_val)

    ! Set connection flux type
    call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_FLUX_TYPE, flux_type)

    ! Set conductance type
    call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_CONDUCTANCE_TYPE, conductance_type)

    call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_CONDUCTANCE_DN, dn_conductance)

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

  subroutine set_material_properties_for_root_bc(tree_name)
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
    character(len=*)   :: tree_name
    !
    PetscInt           :: eqn_id
    PetscInt           :: ncells, nconns

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
    PetscReal          :: phi88, phi50, c1, c2, kmax
    PetscInt :: kk
    PetscReal :: zc_1d, qz, dd, rld0, root_len_den, droot
    PetscInt :: IDX

    select case(trim(tree_name))
    case ('oak')
       eqn_id = GE_oak_root
       ncells = oak_root_nz
       phi88  = oak_rad_root_phi88
       phi50  = oak_rad_root_phi50
       c1     = oak_rad_root_c1
       c2     = oak_rad_root_c2
       kmax   = oak_rad_root_kmax
       qz     = oak_root_qz
       rld0   = oak_root_rld0
       dd     = oak_root_d
    case ('pine')
       eqn_id = GE_pine_root
       ncells = pine_root_nz
       phi88  = pine_rad_root_phi88
       phi50  = pine_rad_root_phi50
       c1     = pine_rad_root_c1
       c2     = pine_rad_root_c2
       kmax   = pine_rad_root_kmax
       qz     = pine_root_qz
       rld0   = pine_root_rld0
       dd     = pine_root_d
    case ('e')
       IDX    = E_IDX
       eqn_id = GE_e_root
       ncells = root_nz(IDX)
       phi88  = rad_root_phi88(IDX)
       phi50  = rad_root_phi50(IDX)
       c1     = rad_root_c1(IDX)
       c2     = rad_root_c2(IDX)
       kmax   = rad_root_kmax(IDX)
       qz     = root_qz(IDX)
       rld0   = root_rld0(IDX)
       dd     = root_d(IDX)
    case ('m')
       IDX    = M_IDX
       eqn_id = GE_m_root
       ncells = root_nz(IDX)
       phi88  = rad_root_phi88(IDX)
       phi50  = rad_root_phi50(IDX)
       c1     = rad_root_c1(IDX)
       c2     = rad_root_c2(IDX)
       kmax   = rad_root_kmax(IDX)
       qz     = root_qz(IDX)
       rld0   = root_rld0(IDX)
       dd     = root_d(IDX)
    case ('o')
       IDX    = O_IDX
       eqn_id = GE_o_root
       ncells = root_nz(IDX)
       phi88  = rad_root_phi88(IDX)
       phi50  = rad_root_phi50(IDX)
       c1     = rad_root_c1(IDX)
       c2     = rad_root_c2(IDX)
       kmax   = rad_root_kmax(IDX)
       qz     = root_qz(IDX)
       rld0   = root_rld0(IDX)
       dd     = root_d(IDX)
    case ('p')
       IDX    = P_IDX
       eqn_id = GE_p_root
       ncells = root_nz(IDX)
       phi88  = rad_root_phi88(IDX)
       phi50  = rad_root_phi50(IDX)
       c1     = rad_root_c1(IDX)
       c2     = rad_root_c2(IDX)
       kmax   = rad_root_kmax(IDX)
       qz     = root_qz(IDX)
       rld0   = root_rld0(IDX)
       dd     = root_d(IDX)
    case default
       write(*,*)'Unknown tree_name: ' //trim(tree_name)
    end select
    
    nconns = ncells + 1
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

    ! Set values for all BCs
    flux_type(1)        = DARCY_FLUX_TYPE
    flux_type(2:nconns) = CONDUCTANCE_FLUX_TYPE
    conductance_type(1) = 0
    conductance_type(2:nconns) = CONDUCTANCE_MANOLI_TYPE

    dn_alpha        (2:nconns) = phi88
    dn_lambda       (2:nconns) = phi50
    dn_residual_sat (2:nconns) = 0.d0
    dn_satfunc_type (2:nconns) = SAT_FUNC_FETCH2
    dn_weibull_d_val(2:nconns) = c1
    dn_weibull_c_val(2:nconns) = c2
    dn_relperm_type (2:nconns) = RELPERM_FUNC_WEIBULL

    do kk = 1, ncells
       zc_1d = -(kk-1)*dz_soil - dz_soil/2.0

       root_len_den = rld0 * (1.d0 - abs(zc_1d)/dd)*exp(-qz*abs(zc_1d)/dd)
       root_len_den = rld0 * exp(-qz*abs(zc_1d)/dd)

       droot = (PI*root_len_den)**(-0.5d0)

       dn_conductance(kk+1) = (kmax*vis/rho)/droot
    end do

    ! Set downwind values
    set_upwind_auxvar(:) = PETSC_FALSE
    call VSFMMPPSetSaturationFunctionAuxVarConn(vsfm_mpp , &
         eqn_id, AUXVAR_CONN_BC                    , &
         set_upwind_auxvar, dn_satfunc_type, dn_alpha    , &
         dn_lambda, dn_residual_sat)

    call VSFMMPPSetRelativePermeabilityAuxVarConn(vsfm_mpp , &
         eqn_id, AUXVAR_CONN_BC                    , &
         set_upwind_auxvar, dn_relperm_type, dn_weibull_d_val, &
         dn_weibull_c_val)
    
    ! Set connection flux type
    call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_FLUX_TYPE, flux_type)

    ! Set conductance type
    call VSFMMPPSetAuxVarConnIntValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_CONDUCTANCE_TYPE, conductance_type)

    call VSFMMPPSetAuxVarConnRealValue(vsfm_mpp, eqn_id, AUXVAR_CONN_BC, &
         VAR_CONDUCTANCE_DN, dn_conductance)

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
  subroutine set_initial_conditions()

    implicit none
    !
    PetscInt           :: ii, jj
    PetscReal          :: zc
    PetscReal, pointer :: press_ic(:)
    Vec                :: p_init, sm_init
    PetscReal, pointer :: sm_init_p(:)
    PetscReal          :: Se, n, Psoil_ic, Pc
    PetscInt           :: soil_idx_offset, offset, GE, nlyrs
    PetscInt           :: e_xylm_offset, m_xylm_offset, p_xylm_offset, o_xylm_offset
    PetscInt           :: e_root_offset, m_root_offset, p_root_offset, o_root_offset
    PetscViewer        :: viewer
    PetscErrorCode     :: ierr

    select case(trim(problem_type))
    case ('oak')
       call set_initial_conditions_for_single_tree(oak_nz)

    case ('pine')
       call set_initial_conditions_for_single_tree(pine_nz)

    case ('oak_and_pine')
       call set_initial_conditions_for_two_trees()
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    case ('oak_spac')

       call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       jj = 0
       do ii = 1, oak_nz
          jj = jj + 1
          if (ii == 1) then
             zc = -0.17d0 + oak_nz*dz_soil
          else
             zc = -dz_xylem/2.d0 - dz_xylem*(ii-1) + oak_nz*dz_xylem
          end if
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)
       enddo
       do ii = 1, oak_root_nz
          jj = jj + 1
          zc = -dz_soil/2.d0 - dz_soil*(ii-1)
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)
       enddo
       do ii = 1, soil_nz
          jj = jj + 1
          zc = -dz_soil/2.d0 - dz_soil*(ii-1)
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)
       enddo
       call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)

    case ('pine_spac')
       call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       jj = 0
       do ii = 1, pine_nz
          jj = jj + 1
          if (ii == 1) then
             zc = -0.17d0 + pine_nz*dz_xylem
          else
             zc = -dz_xylem/2.d0 - dz_xylem*(ii-1) + pine_nz*dz_xylem
          end if
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)
       enddo
       do ii = 1, pine_root_nz
          jj = jj + 1
          zc = -dz_soil/2.d0 - dz_soil*(ii-1)
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)
       enddo
       do ii = 1, soil_nz
          jj = jj + 1
          zc = -dz_soil/2.d0 - dz_soil*(ii-1)
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)
       enddo
       call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    case ('oak_pine_spac')
       call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       jj = 0
       do ii = 1, oak_nz
          jj = jj + 1
          if (ii == 1) then
             zc = -0.17d0 + oak_nz*dz_xylem
          else
             zc = -dz_xylem/2.d0 - dz_xylem*(ii-1) + oak_nz*dz_xylem
          end if
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)/10.d0
       enddo
       do ii = 1, oak_root_nz
          jj = jj + 1
          zc = -dz_soil/2.d0 - dz_soil*(ii-1)
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)/10.d0
       enddo
       do ii = 1, pine_nz
          jj = jj + 1
          if (ii == 1) then
             zc = -0.17d0 + pine_nz*dz_xylem
          else
             zc = -dz_xylem/2.d0 - dz_xylem*(ii-1) + pine_nz*dz_xylem
          end if
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)/10.d0
       enddo
       do ii = 1, pine_root_nz
          jj = jj + 1
          zc = -dz_soil/2.d0 - dz_soil*(ii-1)
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)/10.d0
       enddo
       do ii = 1, soil_nz
          jj = jj + 1
          zc = -dz_soil/2.d0 - dz_soil*(ii-1)
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)/10.d0
       enddo
       call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    case ('soil')
       call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       jj = 0
       do ii = 1, soil_nz
          jj = jj + 1
          zc = -dz_soil/2.d0 - dz_soil*(ii-1)
          press_ic(jj) = 101325.d0 - rho*grav*(zc + 6.d0)/10.d0
       enddo
       press_ic(01) = -124112.449286d0
       press_ic(02) = -33619.504565d0
       press_ic(03) = -39515.138408d0
       press_ic(04) = -57259.991628d0
       press_ic(05) = -86854.064226d0
       press_ic(06) = -102152.206157d0
       press_ic(07) = -103154.417420d0
       press_ic(08) = -104156.628683d0
       press_ic(09) = -105158.839946d0
       press_ic(10) = -106161.051210d0
       press_ic(11) = -120217.157742d0
       press_ic(12) = -147327.159544d0
       press_ic(13) = -174437.161346d0
       press_ic(14) = -201547.163148d0
       press_ic(15) = -228657.164950d0
       press_ic(16) = -228657.164950d0
       press_ic(17) = -228657.164950d0
       press_ic(18) = -228657.164950d0
       press_ic(19) = -228657.164950d0
       press_ic(20) = -228657.164950d0
       press_ic(21) = -228657.164950d0
       press_ic(22) = -228657.164950d0
       press_ic(23) = -228657.164950d0
       press_ic(24) = -228657.164950d0
       press_ic(25) = -228657.164950d0
       press_ic(26) = -228657.164950d0
       press_ic(27) = -228657.164950d0
       press_ic(28) = -228657.164950d0
       press_ic(29) = -228657.164950d0
       press_ic(30) = -228657.164950d0
       call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    case ('emop_spac','m_spac','o_spac','p_spac','e','m','o','p')
       call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       press_ic(:) = 90000.d0
       call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)

       if (ic_file_specified) then
         call PetscViewerBinaryOpen(PETSC_COMM_WORLD, ic_filename, FILE_MODE_READ, &
              viewer, ierr); CHKERRQ(ierr)
         call VecCreate(PETSC_COMM_WORLD, p_init, ierr); CHKERRQ(ierr)
         call VecLoad(p_init, viewer, ierr); CHKERRQ(ierr)
         call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
         call VecCopy(p_init, vsfm_mpp%soe%solver%soln, ierr); CHKERRQ(ierr)
         call VecDestroy(p_init, ierr); CHKERRQ(ierr)
       end if

       if (sm_ic_file_specified) then
          call PetscViewerBinaryOpen(PETSC_COMM_WORLD, sm_ic_filename, FILE_MODE_READ, &
               viewer, ierr); CHKERRQ(ierr)
          call VecCreate(PETSC_COMM_WORLD, sm_init, ierr); CHKERRQ(ierr)
          call VecLoad(sm_init, viewer, ierr); CHKERRQ(ierr)
          call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
          call VecGetArrayF90(sm_init, sm_init_p, ierr); CHKERRQ(ierr)
          call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr);
          n = 1.d0/(1.d0-soil_vg_m)
          soil_idx_offset = 0
          e_xylm_offset = soil_idx_offset; if (GE_e_xylm>0) soil_idx_offset = soil_idx_offset + nz(E_IDX)
          m_xylm_offset = soil_idx_offset; if (GE_m_xylm>0) soil_idx_offset = soil_idx_offset + nz(M_IDX)
          o_xylm_offset = soil_idx_offset; if (GE_o_xylm>0) soil_idx_offset = soil_idx_offset + nz(O_IDX)
          p_xylm_offset = soil_idx_offset; if (GE_p_xylm>0) soil_idx_offset = soil_idx_offset + nz(P_IDX)
          e_root_offset = soil_idx_offset; if (GE_e_root>0) soil_idx_offset = soil_idx_offset + root_nz(E_IDX)
          m_root_offset = soil_idx_offset; if (GE_m_root>0) soil_idx_offset = soil_idx_offset + root_nz(M_IDX)
          o_root_offset = soil_idx_offset; if (GE_o_root>0) soil_idx_offset = soil_idx_offset + root_nz(O_IDX)
          p_root_offset = soil_idx_offset; if (GE_p_root>0) soil_idx_offset = soil_idx_offset + root_nz(P_IDX)
          do jj = 1,soil_nz
             Se = (sm_init_p(jj) - soil_sat_res)/(1.d0 - soil_sat_res)
             Pc = -((Se**(-1.d0/soil_vg_m) - 1.d0)**(1.d0/n))/soil_alpha;
             Psoil_ic = 101325 + Pc
             press_ic(soil_idx_offset+jj) = Psoil_ic
          end do

          ! Set root pressure
          do ii = 1, 4
             select case (ii)
             case (1)
                GE     = GE_e_root
                offset = e_root_offset
                nlyrs   = root_nz(E_IDX)
             case (2)
                GE     = GE_m_root
                offset = m_root_offset
                nlyrs   = root_nz(M_IDX)
             case (3)
                GE     = GE_o_root
                offset = o_root_offset
                nlyrs   = root_nz(O_IDX)
             case (4)
                GE     = GE_p_root
                offset = p_root_offset
                nlyrs   = root_nz(P_IDX)
             end select
             if (GE>0) then
                do jj = 1,nlyrs
                   press_ic(offset+jj) = press_ic(soil_idx_offset+jj)
                end do
             endif
          enddo

          ! Set xylem pressure
          do ii = 1, 4
             select case (ii)
             case (1)
                GE     = GE_e_xylm
                offset = e_xylm_offset
                nlyrs   = nz(E_IDX)
             case (2)
                GE     = GE_m_xylm
                offset = m_xylm_offset
                nlyrs   = nz(M_IDX)
             case (3)
                GE     = GE_o_xylm
                offset = o_xylm_offset
                nlyrs   = nz(O_IDX)
             case (4)
                GE     = GE_p_xylm
                offset = p_xylm_offset
                nlyrs   = nz(P_IDX)
             end select
             if (GE>0) then
                do jj = 1, nlyrs
                   press_ic(offset+jj) = press_ic(soil_idx_offset+1)-1000.d0*9.8d0*nlyrs*dz_xylem/10.d0+1000.d0*9.8d0*(jj)*dz_xylem/10.d0
                end do
             endif
          enddo
          call VecRestoreArrayF90(sm_init, sm_init_p, ierr); CHKERRQ(ierr)
          call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr);
       end if

       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    case ('e_spac')
       call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)
       press_ic(:) = 90000.d0
       call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr)

       if (ic_file_specified) then
         call PetscViewerBinaryOpen(PETSC_COMM_WORLD, ic_filename, FILE_MODE_READ, &
              viewer, ierr); CHKERRQ(ierr)
         call VecCreate(PETSC_COMM_WORLD, p_init, ierr); CHKERRQ(ierr)
         call VecLoad(p_init, viewer, ierr); CHKERRQ(ierr)
         call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
         call VecCopy(p_init, vsfm_mpp%soe%solver%soln, ierr); CHKERRQ(ierr)
         call VecDestroy(p_init, ierr); CHKERRQ(ierr)
       end if

       if (sm_ic_file_specified) then
          call PetscViewerBinaryOpen(PETSC_COMM_WORLD, sm_ic_filename, FILE_MODE_READ, &
               viewer, ierr); CHKERRQ(ierr)
          call VecCreate(PETSC_COMM_WORLD, sm_init, ierr); CHKERRQ(ierr)
          call VecLoad(sm_init, viewer, ierr); CHKERRQ(ierr)
          call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
          call VecGetArrayF90(sm_init, sm_init_p, ierr); CHKERRQ(ierr)
          call VecGetArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr);
          n = 1.d0/(1.d0-soil_vg_m)
          do jj = 1,soil_nz
             Se = (sm_init_p(jj) - soil_sat_res)/(1.d0 - soil_sat_res)
             Pc = -((Se**(-1.d0/soil_vg_m) - 1.d0)**(1.d0/n))/soil_alpha;
             Psoil_ic = 101325 + Pc
             press_ic(nz(E_IDX)+root_nz(E_IDX)+jj) = Psoil_ic
          end do
          do jj = 1,root_nz(E_IDX)
             press_ic(nz(E_IDX)+jj) = press_ic(nz(E_IDX)+root_nz(E_IDX)+jj)
          end do
          do jj = 1,nz(E_IDX)
             press_ic(jj) = press_ic(nz(E_IDX)+1)-1000.d0*9.8d0*nz(E_IDX)*dz_soil/10.d0+1000.d0*9.8d0*(jj)*dz_soil/10.d0
          end do
          call VecRestoreArrayF90(sm_init, sm_init_p, ierr); CHKERRQ(ierr)
          call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, press_ic, ierr); CHKERRQ(ierr);
       end if

       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
       call VecCopy(vsfm_mpp%soe%solver%soln, vsfm_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

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
          press_ic(ii) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (nz - ii)*dz_soil) + 101325.d0
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
          press_ic(ii         ) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (oak_nz - ii)*dz_xylem) + 101325.d0
       enddo

       do ii = 1, pine_nz
          press_ic(ii + oak_nz) = -phi_root_mean_times_beta_s - rho * grav * (0.17d0 + (pine_nz - ii)*dz_xylem) + 101325.d0
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
  subroutine set_source_sink_for_single_tree(nz, istep, ET)
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
    PetscInt           :: istep
    Vec                :: ET
    !
    PetscReal, pointer :: ss_value(:)
    PetscReal, pointer :: et_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk

    allocate(ss_value(nz))

    call VecGetArrayF90(ET, et_p, ierr)
    do kk = 1, nz
       ss_value(kk) = -et_p(nz*(istep-1) + kk) * dz_xylem
    end do
    call VecRestoreArrayF90(ET, et_p, ierr)

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

    deallocate(ss_value)

  end subroutine set_source_sink_for_single_tree

  !------------------------------------------------------------------------
  subroutine set_source_sink_for_soil(istep, SoilSS)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: istep
    Vec                :: SoilSS
    !
    PetscReal, pointer :: ss_value(:)
    PetscReal, pointer :: ss_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk

    allocate(ss_value(1))

    call VecGetArrayF90(SoilSS, ss_p, ierr)
    ss_value(1) = ss_p(istep)/1000.d0
    call VecRestoreArrayF90(SoilSS, ss_p, ierr)

    soe_auxvar_id = 2
    if (GE_oak_xylm > 0 .and. GE_pine_xylm > 1) soe_auxvar_id = 3
    if (GE_soil == 1) soe_auxvar_id = 1
    if (GE_soil == 9) soe_auxvar_id = 5

    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

    deallocate(ss_value)

  end subroutine set_source_sink_for_soil

  !------------------------------------------------------------------------
  subroutine set_boundary_conditions_for_single_tree(istep, SoilBC)
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
    PetscInt           :: istep
    Vec                :: SoilBC
    !
    PetscReal, pointer :: bc_value(:)
    PetscReal, pointer :: et_p(:)
    PetscReal, pointer :: soilbc_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk

    call VecGetArrayF90(SoilBC, soilbc_p, ierr)
    allocate(bc_value(1))
    bc_value(1) = soilbc_p(istep)

    soe_auxvar_id = 1
    !if (GE_e_xylm > 0) soe_auxvar_id = soe_auxvar_id + 1
    !if (GE_m_xylm > 0) soe_auxvar_id = soe_auxvar_id + 1
    !if (GE_o_xylm > 0) soe_auxvar_id = soe_auxvar_id + 1
    !if (GE_p_xylm > 0) soe_auxvar_id = soe_auxvar_id + 1

    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bc_value)
    deallocate(bc_value)
    call VecRestoreArrayF90(SoilBC, soilbc_p, ierr)

  end subroutine set_boundary_conditions_for_single_tree

  !------------------------------------------------------------------------
  subroutine set_soil_moisture_boundary_conditions_for_single_tree(istep, SoilMoistureBC)
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
    PetscInt           :: istep
    Vec                :: SoilMoistureBC
    !
    PetscReal, pointer :: bc_value(:)
    PetscReal, pointer :: smbc_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscReal          :: Se, n, Psoil, Pc

    call VecGetArrayF90(SoilMoistureBC, smbc_p, ierr)
    allocate(bc_value(1))

    n = 1.d0/(1.d0-soil_vg_m)
    Se    = (smbc_p(istep) - soil_sat_res)/(1.d0 - soil_sat_res)
    Pc    = -((Se**(-1.d0/soil_vg_m) - 1.d0)**(1.d0/n))/soil_alpha;
    Psoil = 101325 + Pc
    bc_value(1) = Psoil

    soe_auxvar_id = 1

    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bc_value)

    deallocate(bc_value)
    call VecRestoreArrayF90(SoilMoistureBC, smbc_p, ierr)

  end subroutine set_soil_moisture_boundary_conditions_for_single_tree

  !------------------------------------------------------------------------
  subroutine set_boundary_conditions_for_two_trees(istep, SoilBC)
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
    PetscInt           :: istep
    Vec                :: SoilBC
    !
    PetscReal, pointer :: bc_value(:)
    PetscReal, pointer :: et_p(:)
    PetscReal, pointer :: soilbc_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk

    ! Set BC
    call VecGetArrayF90(SoilBC, soilbc_p, ierr)
    allocate(bc_value(2))
    bc_value(1:2) = soilbc_p(istep)
        soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bc_value)
    deallocate(bc_value)
    call VecRestoreArrayF90(SoilBC, soilbc_p, ierr)

  end subroutine set_boundary_conditions_for_two_trees

  !------------------------------------------------------------------------
  subroutine set_source_sink_for_oak_and_pine(istep, ET)
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
    PetscInt           :: istep
    !
    PetscReal, pointer :: oak_ss_value(:)
    PetscReal, pointer :: pine_ss_value(:)
    PetscViewer        :: viewer
    Vec                :: ET
    PetscReal, pointer :: et_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk, idx

    allocate(oak_ss_value(oak_nz))
    allocate(pine_ss_value(pine_nz))

    ! Set source sink
    call VecGetArrayF90(ET, et_p, ierr)

    idx = (oak_nz+pine_nz)*(istep-1)

    ! Save PET data for oak
    do kk = 1, oak_nz
       idx = idx + 1
       oak_ss_value(kk) = -et_p(idx)*dz_xylem
    end do

    ! Save PET data for pine
    do kk = 1, pine_nz
       idx = idx + 1
       pine_ss_value(kk) = -et_p(idx)*dz_xylem
    end do
    call VecRestoreArrayF90(ET, et_p, ierr)

    ! Set PET for oak
    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, oak_ss_value)

    ! Set PET for pine
    soe_auxvar_id = 2
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, pine_ss_value)

    deallocate(oak_ss_value)
    deallocate(pine_ss_value)

  end subroutine set_source_sink_for_oak_and_pine

  !------------------------------------------------------------------------
  subroutine set_source_sink_for_emop(nET, istep, ET)
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
    PetscInt           :: nET
    PetscInt           :: istep
    Vec                :: ET
    !
    PetscInt           :: ii
    PetscReal, pointer :: ss_value(:)
    PetscViewer        :: viewer
    PetscReal, pointer :: et_p(:)
    PetscErrorCode     :: ierr
    PetscInt           :: soe_auxvar_id
    PetscInt           :: kk, ss_idx, IDX, ieqn

    ! Set source sink
    call VecGetArrayF90(ET, et_p, ierr)

    ss_idx        = nET*(istep-1)
    soe_auxvar_id = 0

    do ii = 1,4
       select case(ii)
       case (1)
          IDX = E_IDX; ieqn = GE_e_xylm;
       case (2)
          IDX = M_IDX; ieqn = GE_m_xylm;
       case (3)
          IDX = O_IDX; ieqn = GE_o_xylm;
       case (4)
          IDX = P_IDX; ieqn = GE_p_xylm;
       end select

       if (ieqn==0) cycle

       allocate(ss_value(nz(IDX)))

       ! Save PET data for xylem
       do kk = 1, nz(IDX)
          ss_idx = ss_idx + 1
          ss_value(kk) = -et_p(ss_idx)*dz_xylem
       end do

       ! Set PET for xylem
       soe_auxvar_id = soe_auxvar_id + 1
       call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_SS,  &
            VAR_BC_SS_CONDITION, soe_auxvar_id, ss_value)

       deallocate(ss_value)
    end do
    call VecRestoreArrayF90(ET, et_p, ierr)


  end subroutine set_source_sink_for_emop

  !------------------------------------------------------------------------
  subroutine add_xylm2root_coupling_bc(tree_name)
    !
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    !
    implicit none
    !
    character(len=*)                     :: tree_name
    !
    PetscInt                             :: mesh_id, region_id
    PetscInt                             :: GE_xylm, GE_root
    PetscInt                             :: ieqn
    PetscInt                             :: num_other_goveqns
    PetscInt                             :: nz_local
    PetscInt                   , pointer :: ieqn_others(:)
    class(connection_set_type) , pointer :: conn_set

    num_other_goveqns = 1 
    allocate(ieqn_others(num_other_goveqns))

    select case(trim(tree_name))
    case ('oak')
       mesh_id = GE_oak_root
       GE_xylm = GE_oak_xylm
       GE_root = GE_oak_root
       nz_local= oak_nz

    case ('pine')
       mesh_id = GE_pine_root
       GE_xylm = GE_pine_xylm
       GE_root = GE_pine_root
       nz_local= pine_nz

    case ('e')
       mesh_id = GE_e_root
       GE_xylm = GE_e_xylm
       GE_root = GE_e_root
       nz_local= nz(E_IDX)

    case ('m')
       mesh_id = GE_m_root
       GE_xylm = GE_m_xylm
       GE_root = GE_m_root
       nz_local= nz(M_IDX)

    case ('o')
       mesh_id = GE_o_root
       GE_xylm = GE_o_xylm
       GE_root = GE_o_root
       nz_local= nz(O_IDX)

    case ('p')
       mesh_id = GE_p_root
       GE_xylm = GE_p_xylm
       GE_root = GE_p_root
       nz_local= nz(P_IDX)

    case default
       write(*,*)'Unable to set root mesh for tree_name = '//trim(tree_name)
       stop
    end select
    
    region_id = 2 ! 

    ! Add xylem BC in root equation
    
    allocate(conn_set)

    ieqn           = GE_root
    ieqn_others(1) = GE_xylm

    call vsfm_mpp%meshes(mesh_id)%conditions_conn_set_list%GetConnectionSet(region_id, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Xylem BC in root equation',         &
         unit='Pa',                                &
         num_other_goveqs=num_other_goveqns,       &
         id_of_other_goveqs=ieqn_others,           &
         conn_set=conn_set)

    ! Add root BC in xylem equation

    allocate(conn_set)

    ieqn           = GE_xylm
    ieqn_others(1) = GE_root
    call vsfm_mpp%meshes(mesh_id)%conditions_conn_set_list%GetConnectionSet(region_id, conn_set)

    ! Change the downward id to be the last control volume in the column
    call conn_set%conn(1)%SetIDDn(nz_local)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Root BC in xylem equation',         &
         unit='Pa',                                &
         num_other_goveqs=num_other_goveqns,       &
         id_of_other_goveqs=ieqn_others,           &
         conn_set=conn_set)

    deallocate(ieqn_others)

  end subroutine add_xylm2root_coupling_bc

  !------------------------------------------------------------------------
  subroutine add_root2soil_coupling_bc(tree_name)
    !
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    !
    implicit none
    !
    character(len=*)                     :: tree_name
    !
    PetscInt                             :: mesh_id, region_id
    PetscInt                             :: GE_root
    PetscInt                             :: ieqn
    PetscInt                             :: num_other_goveqns
    PetscInt                             :: nz
    PetscInt                   , pointer :: ieqn_others(:)
    class(connection_set_type) , pointer :: conn_set

    num_other_goveqns = 1 
    allocate(ieqn_others(num_other_goveqns))

    select case(trim(tree_name))
    case ('oak')
       mesh_id = GE_oak_root
       GE_root = GE_oak_root

    case ('pine')
       mesh_id = GE_pine_root
       GE_root = GE_pine_root

    case ('e')
       mesh_id = GE_e_root
       GE_root = GE_e_root

    case ('m')
       mesh_id = GE_m_root
       GE_root = GE_m_root

    case ('o')
       mesh_id = GE_o_root
       GE_root = GE_o_root

    case ('p')
       mesh_id = GE_p_root
       GE_root = GE_p_root

    case default
       write(*,*)'Unable to set root mesh for tree_name = '//trim(tree_name)
       stop
    end select
    
    region_id = 1 ! 

    ! Add soil BC in root equation
    
    allocate(conn_set)

    ieqn           = GE_root
    ieqn_others(1) = GE_soil

    call vsfm_mpp%meshes(mesh_id)%conditions_conn_set_list%GetConnectionSet(region_id, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Soil BC in root equation',         &
         unit='Pa',                                &
         num_other_goveqs=num_other_goveqns,       &
         id_of_other_goveqs=ieqn_others,           &
         conn_set=conn_set)

    ! Add root BC in soil equation

    allocate(conn_set)

    ieqn           = GE_soil
    ieqn_others(1) = GE_root
    call vsfm_mpp%meshes(mesh_id)%conditions_conn_set_list%GetConnectionSet(region_id, conn_set)

    call vsfm_mpp%soe%AddCouplingBCsInGovEqn(ieqn, &
         name='Root BC in soil equation',         &
         unit='Pa',                                &
         num_other_goveqs=num_other_goveqns,       &
         id_of_other_goveqs=ieqn_others,           &
         conn_set=conn_set)

    deallocate(ieqn_others)

  end subroutine add_root2soil_coupling_bc

  !------------------------------------------------------------------------
  subroutine diagnose_actual_sink_for_single_tree(nz, istep, act_et_p)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use MultiPhysicsProbConstants , only : VAR_MASS_FLUX
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nz
    PetscInt           :: istep
    PetscReal, pointer :: act_et_p(:)
    !
    PetscReal, pointer :: ss_value(:)
    PetscInt           :: kk

    allocate(ss_value(nz))

    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_SS, VAR_MASS_FLUX, 1, ss_value)
    do kk = 1, nz
       act_et_p((istep-1)*nz + kk) = ss_value(kk)
    end do
    
    deallocate(ss_value)

  end subroutine diagnose_actual_sink_for_single_tree

  !------------------------------------------------------------------------
  subroutine diagnose_actual_sink_for_oak_and_pine(istep, act_et_p)
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
    PetscInt           :: istep
    PetscReal, pointer :: act_et_p(:)
    !
    PetscReal, pointer :: oak_ss_value(:)
    PetscReal, pointer :: pine_ss_value(:)
    PetscInt           :: kk, idx

    allocate(oak_ss_value(oak_nz))
    allocate(pine_ss_value(pine_nz))

    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_SS, VAR_MASS_FLUX, 1, oak_ss_value)
    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_SS, VAR_MASS_FLUX, 2, pine_ss_value)

    idx = (istep-1)*(oak_nz+pine_nz)

    do kk = 1, oak_nz
       idx = idx + 1
       act_et_p(idx) = oak_ss_value(kk)
    end do
    
    do kk = 1, pine_nz
       idx = idx + 1
       act_et_p(idx) = pine_ss_value(kk)
    end do
    
    deallocate(oak_ss_value)
    deallocate(pine_ss_value)

  end subroutine diagnose_actual_sink_for_oak_and_pine

  !------------------------------------------------------------------------
  subroutine diagnose_actual_sink_for_emop(nET, istep, act_et_p)
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
    PetscInt           :: nET, istep
    PetscReal, pointer :: act_et_p(:)
    !
    PetscReal, pointer :: ss_value(:)
    PetscInt           :: kk, et_idx, IDX, ii, soe_auxvar_id, ieqn

    et_idx        = (istep-1)*(nET)
    soe_auxvar_id = 0

    do ii = 1, 4

       select case(ii)
       case (1)
          IDX = E_IDX; ieqn = GE_e_xylm
       case (2)
          IDX = M_IDX; ieqn = GE_m_xylm
       case (3)
          IDX = O_IDX; ieqn = GE_o_xylm
       case (4)
          IDX = P_IDX; ieqn = GE_p_xylm
       end select

       if (ieqn == 0) cycle

       soe_auxvar_id = soe_auxvar_id + 1
       allocate(ss_value(nz(IDX)))

       call vsfm_mpp%soe%GetDataForCLM(AUXVAR_SS, VAR_MASS_FLUX, soe_auxvar_id, ss_value)

       do kk = 1, nz(IDX)
          et_idx = et_idx + 1
          act_et_p(et_idx) = ss_value(kk)
       end do
    
       deallocate(ss_value)

    end do

  end subroutine diagnose_actual_sink_for_emop

  !------------------------------------------------------------------------
  subroutine save_internal_mass_flux_data(soe, num_internal_mass_flux, istep, mass_flux_p)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : VAR_MASS_FLUX
    use SystemOfEquationsVSFMType, only : sysofeqns_vsfm_type
    use SystemOfEquationsBaseType, only : sysofeqns_base_type
    use petscsys
    use petscvec
    !
    implicit none
    !
    class(sysofeqns_base_type) :: soe
    PetscInt           :: num_internal_mass_flux
    PetscInt           :: istep
    PetscReal, pointer :: mass_flux_p(:)
    !
    PetscReal, pointer :: ss_value(:)
    PetscInt           :: iauxvar

    select type(soe)
    class is (sysofeqns_vsfm_type)
       do iauxvar = 1, soe%num_auxvars_conn_in
          mass_flux_p((istep-1)*num_internal_mass_flux + iauxvar) = soe%aux_vars_conn_in(iauxvar)%mass_flux
       end do
    end select

  end subroutine save_internal_mass_flux_data

  !------------------------------------------------------------------------
  subroutine save_boundary_mass_flux_data(soe, num_boundary_mass_flux, istep, mass_flux_p)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : VAR_MASS_FLUX
    use SystemOfEquationsVSFMType, only : sysofeqns_vsfm_type
    use SystemOfEquationsBaseType, only : sysofeqns_base_type
    use petscsys
    use petscvec
    !
    implicit none
    !
    class(sysofeqns_base_type) :: soe
    PetscInt           :: num_boundary_mass_flux
    PetscInt           :: istep
    PetscReal, pointer :: mass_flux_p(:)
    !
    PetscReal, pointer :: ss_value(:)
    PetscInt           :: iauxvar

    select type(soe)
    class is (sysofeqns_vsfm_type)
       do iauxvar = 1, soe%num_auxvars_bc
          mass_flux_p((istep-1)*num_boundary_mass_flux + iauxvar) = soe%aux_vars_bc(iauxvar)%mass_flux
       end do
    end select

  end subroutine save_boundary_mass_flux_data

  !------------------------------------------------------------------------
  subroutine save_coupling_mass_flux_data(soe, num_coupling_mass_flux, istep, mass_flux_p)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : VAR_MASS_FLUX
    use SystemOfEquationsVSFMType, only : sysofeqns_vsfm_type
    use SystemOfEquationsBaseType, only : sysofeqns_base_type
    use petscsys
    use petscvec
    !
    implicit none
    !
    class(sysofeqns_base_type) :: soe
    PetscInt           :: num_coupling_mass_flux
    PetscInt           :: istep
    PetscReal, pointer :: mass_flux_p(:)
    !
    PetscReal, pointer :: ss_value(:)
    PetscInt           :: iauxvar

    select type(soe)
    class is (sysofeqns_vsfm_type)
       do iauxvar = 1, soe%num_auxvars_bc_otr_geqn
          mass_flux_p((istep-1)*num_coupling_mass_flux + iauxvar) = soe%aux_vars_bc_otr_geqn(iauxvar)%mass_flux
       end do
    end select

  end subroutine save_coupling_mass_flux_data

  !------------------------------------------------------------------------
  subroutine save_saturation(nz, istep, liq_sat_p)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nz
    PetscInt           :: istep
    PetscReal, pointer :: liq_sat_p(:)
    !
    PetscReal, pointer :: liq_sat(:)
    PetscInt           :: kk, idx

    allocate(liq_sat(nz))

    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL, VAR_LIQ_SAT, -1, liq_sat)

    idx = (istep-1)*nz

    do kk = 1, nz
       idx = idx + 1
       liq_sat_p(idx) = liq_sat(kk)
    end do
    
    deallocate(liq_sat)

  end subroutine save_saturation

  !------------------------------------------------------------------------
  subroutine save_liquid_mass(nz, istep, liq_mass_p)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nz
    PetscInt           :: istep
    PetscReal, pointer :: liq_mass_p(:)
    !
    PetscReal, pointer :: liq_mass(:)
    PetscInt           :: kk, idx

    allocate(liq_mass(nz))

    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL, VAR_MASS, -1, liq_mass)

    idx = (istep-1)*nz

    do kk = 1, nz
       idx = idx + 1
       liq_mass_p(idx) = liq_mass(kk)
    end do
    
    deallocate(liq_mass)

  end subroutine save_liquid_mass

  !------------------------------------------------------------------------
  subroutine save_pressure_data(nz, istep, pressure_p)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use petscsys
    use petscvec
    !
    implicit none
    !
    PetscInt           :: nz
    PetscInt           :: istep
    PetscReal, pointer :: pressure_p(:)
    !
    PetscReal, pointer :: soln_p(:)
    PetscInt           :: kk, idx
    PetscErrorCode     :: ierr

    call VecGetArrayF90(vsfm_mpp%soe%solver%soln, soln_p, ierr)

    idx = (istep-1)*nz

    do kk = 1, nz
       idx = idx + 1
       pressure_p(idx) = soln_p(kk)
    end do
    
    call VecRestoreArrayF90(vsfm_mpp%soe%solver%soln, soln_p, ierr)

  end subroutine save_pressure_data

end module vsfm_spac_fetch2_problem
