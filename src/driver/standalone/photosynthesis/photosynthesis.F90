program photosynthesis
  !
#include <petsc/finclude/petsc.h>
  !
  use photosynthesis_problem, only : run_photosynthesis_problem
  use petscsys
  !
  implicit none
  !
  PetscErrorCode     :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_photosynthesis_problem()

  call PetscFinalize(ierr)

end program photosynthesis
  
