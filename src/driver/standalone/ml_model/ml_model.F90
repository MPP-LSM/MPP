program ml_model
  !
#include <petsc/finclude/petsc.h>
  !
  !use ml_model_problem, only : run_ml_model_problem
  use petscsys
  !
  implicit none
  !
  PetscErrorCode     :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  !call run_ml_model_problem()

  call PetscFinalize(ierr)

end program ml_model
