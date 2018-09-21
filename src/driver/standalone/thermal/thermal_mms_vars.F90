module thermal_mms_vars

  implicit none

#include <petsc/finclude/petsc.h>

  PetscInt              :: nx,ny,nz
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost
  PetscReal             :: dx, dy, dz
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscReal , pointer   :: soil_xc_3d(:,:,:)
  PetscReal , pointer   :: soil_yc_3d(:,:,:)
  PetscReal , pointer   :: soil_zc_3d(:,:,:)
  PetscInt  , pointer   :: soil_id(:,:,:)
  PetscReal , parameter :: PI = 4 * atan(1.0d0)

  PetscInt  , parameter :: STEADY_STATE_1D = 1
  PetscInt  , parameter :: STEADY_STATE_2D = 2
  PetscInt  , parameter :: STEADY_STATE_3D = 3

  PetscInt  , parameter :: DATA_XC                          = 1
  PetscInt  , parameter :: DATA_YC                          = 2
  PetscInt  , parameter :: DATA_ZC                          = 3
  PetscInt  , parameter :: DATA_THERMAL_CONDUCTIVITY        = 4
  PetscInt  , parameter :: DATA_TEMPERATURE_BC              = 5
  PetscInt  , parameter :: DATA_HEAT_SOURCE                 = 6
  PetscInt  , parameter :: DATA_TEMPERATURE                 = 7
  PetscInt  , parameter :: DATA_THERMAL_CONDUCTIVITY_1D_VEC = 8

end module thermal_mms_vars
