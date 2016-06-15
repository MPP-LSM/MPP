program main

  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use mpp_varpar           , only : mpp_varpar_init
  !
  implicit none
  !
#include "finclude/petscsys.h"
  !
  PetscInt       :: nlevsoi, nlevgrnd, nlevsno
  PetscErrorCode :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD

  call mpp_varpar_init()

  call PetscFinalize(ierr)

end program main
