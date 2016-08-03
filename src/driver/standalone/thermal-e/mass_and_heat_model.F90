program mass_and_heat_model
  !
  use mass_and_heat_model_problem , only : run_mass_and_heat_model_problem
  !
  implicit none
  !
#include "finclude/petscsys.h"
  !
  PetscErrorCode     :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_mass_and_heat_model_problem()

  call PetscFinalize(ierr)

end program mass_and_heat_model
