module ml_model_global_vars

  implicit none

#include <petsc/finclude/petsc.h>

  type, public :: condition_type
     PetscInt :: ndata
     PetscReal, pointer :: data(:)
  end type condition_type
  
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

  ! Vegetation parameters
  PetscReal, pointer :: dlai(:), dsai(:), dpai(:), sumpai(:), cumpai(:), fssh(:), leaf_td(:)

  ! Problem parameters
  PetscInt  :: ncair
  PetscInt  :: ntree

  ! Boundary conditions
  type(condition_type) :: Iskyb_vis,  Iskyd_vis, Iskyb_nir,  Iskyd_nir
  type(condition_type) :: Irsky
  type(condition_type) :: Ileaf_sun_vis, Ileaf_shd_vis
  type(condition_type) :: Ileaf_sun_nir, Ileaf_shd_nir
  type(condition_type) :: Isoil_vis, Isoil_nir
  type(condition_type) :: Labs_leaf_sun, Labs_leaf_shd, Labs_soil
  type(condition_type) :: Pref, Uref, Tref, Rhref, Qref, eair
  type(condition_type) :: Tleaf_sun, Tleaf_shd ! dimension = nveg
  type(condition_type) :: Tcan
  type(condition_type) :: Tair, qair, wind ! dimension = ncan
  type(condition_type) :: gs_sun, gs_shd
  type(condition_type) :: gbh, gbv, gbc
  type(condition_type) :: co2ref, o2ref
  type(condition_type) :: Albsoib_vis, Albsoib_nir
  type(condition_type) :: Albsoid_vis, Albsoid_nir
  type(condition_type) :: tg, soil_t
  type(condition_type) :: sza

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

