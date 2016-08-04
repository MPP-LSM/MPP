program vsfm_vchannel
  !
  use vsfm_vchannel_problem      , only : run_vsfm_vchannel_problem
  !
  implicit none
  !
#include "finclude/petscsys.h"
  !
  PetscErrorCode     :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_vsfm_vchannel_problem()

  call PetscFinalize(ierr)

end program vsfm_vchannel
