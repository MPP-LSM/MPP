program mlc
  !
#include <petsc/finclude/petsc.h>
  !
  use mlc_problem, only : run_mlc_problem
  use petscsys
  !
  implicit none
  !
  PetscErrorCode     :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_mlc_problem()

  call PetscFinalize(ierr)

end program mlc
