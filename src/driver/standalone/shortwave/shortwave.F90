program shortwave
  !
#include <petsc/finclude/petsc.h>
  !
  use shortwave_problem, only : run_shortwave_problem
  use petscsys
  !
  implicit none
  !
  PetscErrorCode     :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_shortwave_problem()

  call PetscFinalize(ierr)

end program shortwave
  
