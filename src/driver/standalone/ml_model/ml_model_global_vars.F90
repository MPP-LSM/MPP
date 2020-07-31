module ml_model_global_vars

  implicit none

#include <petsc/finclude/petsc.h>

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
  PetscInt :: nbot, ntop

  ! Vegetation parameters
  PetscReal, pointer :: dpai(:), cumlai(:), fssh(:)

  ! Problem parameters
  PetscInt  :: ncair
  PetscInt  :: ntree

  ! Boundary conditions
  PetscReal, pointer :: Iskyb_vis(:),  Iskyd_vis(:), Iskyb_nir(:),  Iskyd_nir(:)
  PetscReal, pointer :: Irsky(:)
  PetscReal, pointer :: Pref(:), Uref(:), Tref(:), Rhref(:)

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

