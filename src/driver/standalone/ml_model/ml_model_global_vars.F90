module ml_model_global_vars

  implicit none

#include <petsc/finclude/petsc.h>

  type, public :: condition_type
     PetscInt :: ndata
     PetscReal, pointer :: data(:)
  end type condition_type

  type, public :: boundary_condition_type
     type(condition_type) :: iskyb_vis
     type(condition_type) :: iskyd_vis
     type(condition_type) :: iskyb_nir
     type(condition_type) :: iskyd_nir
     type(condition_type) :: Irsky
     type(condition_type) :: tref
     type(condition_type) :: qref
     type(condition_type) :: pref
     type(condition_type) :: uref
     type(condition_type) :: co2ref
     type(condition_type) :: o2ref
     type(condition_type) :: albsoib_vis
     type(condition_type) :: albsoib_nir
     type(condition_type) :: albsoid_vis
     type(condition_type) :: albsoid_nir
     type(condition_type) :: tg
     type(condition_type) :: soil_t
     type(condition_type) :: sza
     type(condition_type) :: rhg
     type(condition_type) :: soilres
     type(condition_type) :: soil_tk
     type(condition_type) :: pref_prev
     type(condition_type) :: h2osoi_vol
     type(condition_type) :: fwet
     type(condition_type) :: fdry
  end type boundary_condition_type

  type, public :: internal_condition_type

     ! shortwave model
     type(condition_type) :: ileaf_sun_vis
     type(condition_type) :: ileaf_shd_vis
     type(condition_type) :: ileaf_sun_nir
     type(condition_type) :: ileaf_shd_nir
     type(condition_type) :: isoil_vis
     type(condition_type) :: isoil_nir

     ! longwave model
     type(condition_type) :: Labs_leaf_sun
     type(condition_type) :: Labs_leaf_shd
     type(condition_type) :: Labs_soil

     ! leaf boundary layer model
     type(condition_type) :: gbh
     type(condition_type) :: gbv
     type(condition_type) :: gbc

     ! photosynthesis model
     type(condition_type) :: gs_sun
     type(condition_type) :: gs_shd

     ! multi-layer canopy model
     type(condition_type) :: Tcan
     type(condition_type) :: Tleaf_sun
     type(condition_type) :: Tleaf_shd
     type(condition_type) :: Tair
     type(condition_type) :: qair
     type(condition_type) :: wind

  end type internal_condition_type

  type, public :: canopy_level_accumulator
     type(condition_type) :: ustar       ! friction velocity [m/s]
     type(condition_type) :: lup         ! upward longwave radiation above canopy [W/m^2]

     type(condition_type) :: labs_soi    ! absorbed longwave radiation [W/m^2]
     type(condition_type) :: rnabs_soi   ! net radiation [W/m^2]
     type(condition_type) :: sh_soi      ! sensible heat flux [W/m^2]
     type(condition_type) :: lh_soi      ! latent heat flux [W/m^2]
     type(condition_type) :: et_soi      ! water vapor flux [mol H2O/m^2/s]
     type(condition_type) :: g_soi       ! ground heat flux [W/m^2]
     type(condition_type) :: gac0_soi    ! aerodynamic conductance for soil [mol/m^2/s]
  end type canopy_level_accumulator

  type, public :: vertical_level_accumulator
     type(condition_type) :: sh_air      ! sensible heat flux [W/m^2]
     type(condition_type) :: et_air      ! water vapor flux [mol H2O/m^2/s]
     type(condition_type) :: st_air      ! heat storage [W/m^2]
     type(condition_type) :: gac_air     ! aerodynamic conductance [mol/m^2/s]

     type(condition_type) :: labs_leaf_sun   ! absorbed longwave radiation [W/m^2]
     type(condition_type) :: labs_leaf_shd   ! absorbed longwave radiation [W/m^2]
     type(condition_type) :: rn_leaf_sun     ! net radiation [W/m^2]
     type(condition_type) :: rn_leaf_shd     ! net radiation [W/m^2]
     type(condition_type) :: sh_leaf_sun     ! sensible heat flux [W/m^2]
     type(condition_type) :: sh_leaf_shd     ! sensible heat flux [W/m^2]
     type(condition_type) :: lh_leaf_sun     ! latent heat flux [W/m^2]
     type(condition_type) :: lh_leaf_shd     ! latent heat flux [W/m^2]
     type(condition_type) :: tr_leaf_sun     ! transpiration flux [mol H2O/m^2_leaf/s]
     type(condition_type) :: tr_leaf_shd     ! transpiration flux [mol H2O/m^2_leaf/s]
     type(condition_type) :: st_leaf_sun     ! heat storage [W/m^2]
     type(condition_type) :: st_leaf_shd     ! heat storage [W/m^2]
     type(condition_type) :: anet_leaf_sun   ! net photosynthesis [umol CO2/m^2_leaf/s]
     type(condition_type) :: anet_leaf_shd   ! net photosynthesis [umol CO2/m^2_leaf/s]
     type(condition_type) :: agross_leaf_sun ! gross photosynthesis [umo CO2/m^2_leaf/s]
     type(condition_type) :: agross_leaf_shd ! gross photosynthesis [umo CO2/m^2_leaf/s]
     type(condition_type) :: gs_leaf_sun     ! stomatal conductance [mol H2O/m^2_leaf/s]
     type(condition_type) :: gs_leaf_shd     ! stomatal conductance [mol H2O/m^2_leaf/s]
  end type vertical_level_accumulator


  character(len=1024) :: bc_file
  character(len=1024) :: ic_file

  ! Mesh attributes
  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x        = 1.d0
  PetscReal , parameter :: y        = 1.d0
  PetscReal , parameter :: dz_cair  = 0.5d0
  PetscReal , parameter :: z_cair   = 46.d0
  PetscInt  , parameter :: nz_cair  = 92
  PetscReal , parameter :: hc       = 21.d0
  PetscInt  , parameter :: nveg     = 42
  PetscInt              :: nbot, ntop
  PetscBool             :: output_data
  PetscBool             :: checkpoint_data
  PetscBool             :: use_ic

  ! Vegetation parameters
  PetscReal, pointer :: dlai(:), dsai(:), dpai(:), sumpai(:), cumpai(:), fssh(:), leaf_td(:)

  ! Problem parameters
  PetscInt  :: ncair
  PetscInt  :: ntree
  PetscInt  :: beg_step, end_step, nsubstep

  ! Boundary conditions
  type(boundary_condition_type)    :: bnd_cond
  type(internal_condition_type)    :: int_cond
  type(canopy_level_accumulator)   :: canp_lev_vars
  type(vertical_level_accumulator) :: vert_lev_vars


  ! Shortwave model
  PetscInt :: SHORTWAVE_MESH
  PetscInt :: SHORTWAVE_GE

  ! Longwave model
  PetscInt :: LONGWAVE_MESH
  PetscInt :: LONGWAVE_GE

  ! Leaf boundary layer model
  PetscInt :: LBL_MESH
  PetscInt :: LBL_GE

  ! Photosynthesis model
  PetscInt :: PHOTOSYNTHESIS_MESH
  PetscInt :: PHOTOSYNTHESIS_GE
  PetscInt :: c3psn, gstype

  ! Multi-layer canopy model
  PetscInt :: CAIR_MESH
  PetscInt :: CLEF_MESH
  PetscInt :: SOIL_MESH

  PetscInt :: CAIR_TEMP_GE
  PetscInt :: CAIR_VAPR_GE
  PetscInt :: CLEF_TEMP_SUN_GE
  PetscInt :: CLEF_TEMP_SHD_GE

end module ml_model_global_vars

