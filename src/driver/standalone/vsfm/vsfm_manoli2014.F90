program vsfm_manoli2014
  !
  use vsfm_manoli2014_problem      , only : run_vsfm_manoli2014_problem
  !
  implicit none
  !
#include "finclude/petscsys.h"
  !
  PetscErrorCode     :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  write(*,*)'call run_vsfm_manoli2014_problem()'
  call run_vsfm_manoli2014_problem()

  call PetscFinalize(ierr)

end program vsfm_manoli2014
