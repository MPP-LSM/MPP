module MathUtilsMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Math tools
  !
  ! !USES:
  use mpp_varctl, only : iulog
  use mpp_abortutils, only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: quadratic          ! Solve a quadratic equation for its two roots

contains

  !-----------------------------------------------------------------------
  subroutine quadratic (a, b, c, r1, r2)
    !
    ! !DESCRIPTION:
    ! Solve a quadratic equation for its two roots
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in)  :: a,b,c       ! Terms for quadratic equation
    PetscReal, intent(out) :: r1,r2       ! Roots of quadratic equation
    !
    ! !LOCAL VARIABLES:
    PetscReal :: q                        ! Temporary term for quadratic solution
    !---------------------------------------------------------------------

    if (a == 0.d0) then
       write (iulog,*) 'Quadratic solution error: a = ',a
       call endrun()
    end if

    if (b >= 0.d0) then
       q = -0.5d0 * (b + sqrt(b*b - 4.d0*a*c))
    else
       q = -0.5d0 * (b - sqrt(b*b - 4.d0*a*c))
    end if

    r1 = q / a
    if (q /= 0.d0) then
       r2 = c / q
    else
       r2 = 1.e36
    end if

  end subroutine quadratic

#endif

end module MathUtilsMod
