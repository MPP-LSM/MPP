module ml_model_global_vars

  implicit none

#include <petsc/finclude/petsc.h>

  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x        = 1.d0
  PetscReal , parameter :: y        = 1.d0
  PetscReal , parameter :: dz_cair  = 0.5d0
  PetscReal , parameter :: z_cair   = 46.d0
  PetscInt  , parameter :: nz_cair  = 92

  PetscReal , parameter :: hc       = 21.d0
  PetscInt  , parameter :: nveg     = 42

  PetscInt :: nbot
  PetscInt :: ntop

  PetscReal, pointer :: dpai(:), cumlai(:), fssh(:)

  PetscInt  :: ncair
  PetscInt  :: ntree

  PetscInt :: SHORTWAVE_MESH
  PetscInt :: SHORTWAVE_GE

  PetscInt :: LONGWAVE_MESH
  PetscInt :: LONGWAVE_GE

  PetscInt :: LBL_MESH
  PetscInt :: LBL_GE

  PetscInt :: CAIR_MESH
  PetscInt :: CLEF_MESH
  PetscInt :: SOIL_MESH

  PetscInt :: CAIR_TEMP_GE
  PetscInt :: CAIR_VAPR_GE
  PetscInt :: CLEF_TEMP_SUN_GE
  PetscInt :: CLEF_TEMP_SHD_GE
  PetscInt :: SOIL_TEMP_GE
  PetscInt :: SOIL_VAPR_GE

  PetscInt :: CLEF_REGION_IN_CAIR_MESH
  PetscInt :: CAIR_REGION_IN_CLEF_MESH

end module ml_model_global_vars

