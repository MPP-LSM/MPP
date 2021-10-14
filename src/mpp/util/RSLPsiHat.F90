module RSLPsiHat
#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
  use petscsys
  use petscvec

  implicit none

  PetscInt, parameter :: nZ = 276, nL = 41 ! Dimensions of RSL psihat look-up tables
  PetscReal :: zdtgridM(nZ,1)              ! Grid of zdt on which psihat is given for momentum
  PetscReal :: dtLgridM(1,nL)              ! Grid of dtL on which psihat is given for momentum
  PetscReal :: psigridM(nZ,nL)             ! Grid of psihat values for momentum
  PetscReal :: zdtgridH(nZ,1)              ! Grid of zdt on which psihat is given for heat
  PetscReal :: dtLgridH(1,nL)              ! Grid of dtL on which psihat is given for heat
  PetscReal :: psigridH(nZ,nL)             ! Grid of psihat values for heat

  public :: InitializePsiHat

contains

  !------------------------------------------------------------------------
  subroutine InitializePsiHat()
    !
    implicit none
    !
    character(len=1024) , parameter :: psihat_filename = '../../share/rsl_psihat.bin'
    PetscInt                        :: size , ii, jj, count
    Vec                             :: psihat_vec
    PetscReal           , pointer   :: psihat_p(:)
    PetscViewer                     :: viewer
    PetscErrorCode                  :: ierr

    ! Open the binary file
    call PetscViewerBinaryOpen(PETSC_COMM_SELF, psihat_filename, FILE_MODE_READ, viewer, ierr); CHKERRQ(ierr)
    call VecCreate(PETSC_COMM_SELF, psihat_vec, ierr); CHKERRQ(ierr)
    call VecLoad(psihat_vec, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(psihat_vec, psihat_p, ierr); CHKERRA(ierr)

    call VecGetSize(psihat_vec, size, ierr); CHKERRA(ierr)

    if (size /= (nZ + nL + nZ*nL)*2) then
       write(iulog,*)'ERROR: The size of PETSc vector for PsiHat is not of correct size'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    count = 0

    ! The data in the PETSc vector is arranged in the following order
    ! 1. dtLgridM
    ! 2. zdtgridM
    ! 3. psigridM
    ! 4. dtLgridH
    ! 5. zdtgridH
    ! 6. psigridH
    !
    do jj = 1, nL
       count = count + 1;
       dtLgridM(1,jj) = psihat_p(count)
    end do

    do ii = 1, nZ
       count = count + 1
       zdtgridM(ii,1) = psihat_p(count)
    end do

    do ii = 1, nZ
       do jj = 1, nL
          count = count + 1
          psigridM(ii,jj) = psihat_p(count)
       end do
    end do

    do jj = 1, nL
       count = count + 1;
       dtLgridH(1,jj) = psihat_p(count)
    end do

    do ii = 1, nZ
       count = count + 1
       zdtgridH(ii,1) = psihat_p(count)
    end do

    do ii = 1, nZ
       do jj = 1, nL
          count = count + 1
          psigridH(ii,jj) = psihat_p(count)
       end do
    end do

    call VecRestoreArrayF90(psihat_vec, psihat_p, ierr); CHKERRA(ierr)
    call VecDestroy(psihat_vec, ierr); CHKERRA(ierr)

  end subroutine InitializePsiHat

#endif

end module RSLPsiHat
