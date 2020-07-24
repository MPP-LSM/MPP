module shortwave_global_vars
  implicit none

#include <petsc/finclude/petsc.h>

  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x        = 1.d0
  PetscReal , parameter :: y        = 1.d0
  PetscReal , parameter :: z_cair   = 46.d0
  PetscReal , parameter :: z_cleaf  = 46.d0
  PetscInt  , parameter :: nz_cair  = 60
  PetscInt  , parameter :: nz_cleaf = 60
  PetscInt  :: ncair
  PetscInt  :: ntree
  PetscReal , parameter :: hc       = 21.d0

  PetscInt, parameter :: nbot = 6
  PetscInt, parameter :: ntop = 42

  PetscInt :: SHORTWAVE_GE
  PetscInt :: SHORTWAVE_MESH

end module shortwave_global_vars
