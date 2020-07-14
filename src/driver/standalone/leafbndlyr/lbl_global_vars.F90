module lbl_global_vars

  implicit none

#include <petsc/finclude/petsc.h>

  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x        = 1.d0
  PetscReal , parameter :: y        = 1.d0
  PetscReal , parameter :: z_cair   = 46.d0
  !PetscReal , parameter :: z_cleaf  = 21.d0
  PetscReal , parameter :: z_cleaf  = 46.d0
  PetscInt  , parameter :: nz_cair  = 92
  !PetscInt  , parameter :: nz_cleaf = 42
  PetscInt  , parameter :: nz_cleaf = 92
  PetscInt  :: ncair
  PetscInt  :: ntree
  PetscReal , parameter :: hc       = 21.d0

  PetscInt, parameter :: nbot = 6
  PetscInt, parameter :: ntop = 42

  PetscInt :: LBL_MESH
  PetscInt :: LBL_GE

end module lbl_global_vars

