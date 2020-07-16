program longwave
  !
#include <petsc/finclude/petsc.h>
  !
  use longwave_problem, only : run_longwave_problem
  use petscsys
  !
  implicit none
  !
  PetscErrorCode     :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_longwave_problem()

  call PetscFinalize(ierr)

end program longwave
  
