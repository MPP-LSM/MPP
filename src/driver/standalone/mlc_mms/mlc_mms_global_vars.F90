module mlc_mms_global_vars

  implicit none

#include <petsc/finclude/petsc.h>

  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x        = 1.d0
  PetscReal , parameter :: y        = 1.d0

  PetscReal :: hc
  PetscReal :: z_cair
  !PetscReal :: z_cleaf

  PetscInt :: nz_cair, nz_cleaf
  PetscInt :: ncair
  PetscInt :: ntree

  PetscReal :: zmax_l, zmin_l, z_l
  PetscReal :: dz_cair

  PetscInt :: nbot, ntop

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

end module mlc_mms_global_vars

