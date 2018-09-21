module vsfm_mms_vars

  implicit none

#include <petsc/finclude/petsc.h>

  PetscReal , parameter :: PI = 4 * atan(1.0d0)

  PetscInt              :: problem_type
  
  PetscInt              :: nx,ny,nz
  PetscInt              :: ncells
  PetscInt              :: ncells_ghost

  PetscReal             :: dx, dy, dz
  PetscReal             :: xlim, ylim, zlim
  PetscReal             :: x_min , x_max
  PetscReal             :: y_min , y_max
  PetscReal             :: z_min , z_max
  
  PetscReal , pointer   :: soil_xc_3d(:,:,:)
  PetscReal , pointer   :: soil_yc_3d(:,:,:)
  PetscReal , pointer   :: soil_zc_3d(:,:,:)
  PetscInt  , pointer   :: soil_id_3d(:,:,:)

  PetscBool             :: fully_saturated
  
  PetscInt  , parameter :: STEADY_STATE_SOIL_ONLY_1D = 1

  PetscInt  , parameter :: DATA_XC                   = 1
  PetscInt  , parameter :: DATA_YC                   = 2
  PetscInt  , parameter :: DATA_ZC                   = 3
  PetscInt  , parameter :: DATA_PERMEABILITY         = 4
  PetscInt  , parameter :: DATA_PRESSURE_BC          = 5
  PetscInt  , parameter :: DATA_MASS_SOURCE          = 6
  PetscInt  , parameter :: DATA_PRESSURE             = 7
  PetscInt  , parameter :: DATA_POROSITY             = 8
  PetscInt  , parameter :: DATA_SATFUNC_TYPE         = 9
  PetscInt  , parameter :: DATA_SATFUNC_ALPHA        = 10
  PetscInt  , parameter :: DATA_SATFUNC_LAMBDA       = 11
  PetscInt  , parameter :: DATA_RES_SAT              = 12
  PetscInt  , parameter :: DATA_INITIAL_PRESSURE     = 13
  PetscInt  , parameter :: DATA_LIQUID_SATURATION    = 14

end module vsfm_mms_vars
