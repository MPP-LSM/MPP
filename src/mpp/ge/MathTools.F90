module MathToolsMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Math tools
  !
  ! !USES:
  use mpp_varctl, only : iulog
  use mpp_abortutils, only : endrun
  use CanopyTurbulenceAuxType , only : canopy_turbulence_auxvar_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: hybrid             ! Solve for the root of a function using secant and Brent's methods
  public :: zbrent             ! Use Brent's method to find the root of a function
  public :: tridiag            ! Solve a tridiagonal system of equations
  public :: beta_function      ! Evaluate the beta function at p and q: B(p,q)
  public :: log_gamma_function ! Evaluate the log natural of the gamma function at x: ln(G(x))

  interface
    function xfunc (icair, cturb, x) result(f)
    use CanopyTurbulenceAuxType, only : canopy_turbulence_auxvar_type
    integer :: icair
    PetscReal :: x, f
    class(canopy_turbulence_auxvar_type) :: cturb
    end function xfunc
  end interface

contains

  !-----------------------------------------------------------------------
  function hybrid (msg, icair, cturb, func, xa, xb, tol) result(root)
    !
    ! !DESCRIPTION:
    ! Solve for the root of a function, given initial estimates xa and xb.
    ! The root is updated until its accuracy is tol.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*)                    :: msg        ! String to be printed
    integer             , intent(in)    :: icair      ! index of canopy air
    procedure (xfunc)                   :: func       ! Function to solve
    PetscReal           , intent(in)    :: xa, xb     ! Initial estimates of root
    PetscReal           , intent(in)    :: tol        ! Error tolerance
    type(canopy_turbulence_auxvar_type) , intent(inout) :: cturb
    !
    ! !LOCAL VARIABLES:
    PetscReal                           :: x0, x1     ! Estimates of root
    PetscReal                           :: f0, f1     ! Function value for x0 and x1
    PetscReal                           :: minx       ! x0 or x1 that gives smallest function value
    PetscReal                           :: minf       ! Smallest function value using x0 or x1
    PetscReal                           :: root       ! Root

    integer                             :: iter       ! Iteration loop index
    PetscReal                           :: dx         ! Change in root
    PetscReal                           :: x          ! Updated root
    integer             , parameter     :: itmax = 40 ! Maximum number of iterations
    !---------------------------------------------------------------------

    x0 = xa
    f0 = func (icair, cturb, x0)
    if (f0 == 0.d0) then
       root = x0
       return
    end if

    x1 = xb
    f1 = func (icair, cturb, x1)
    if (f1 == 0.d0) then
       root = x1
       return
    end if

    if (f1 < f0) then
       minx = x1
       minf = f1
    else
       minx = x0
       minf = f0
    end if

    ! First use the secant method, and then use the brent method as a backup

    iter = 0
    do
       iter = iter + 1
       dx = -f1 * (x1 - x0) / (f1 - f0)
       x = x1 + dx
       if (abs(dx) < tol) then
          x0 = x
          exit
       end if
       x0 = x1
       f0 = f1
       x1 = x
       f1 = func (icair, cturb, x1)
       if (f1 < minf) then
          minx = x1
          minf = f1
       end if

       ! If a root zone is found, use the brent method for a robust backup strategy

       if (f1 * f0 < 0.d0) then
          x = zbrent (msg, icair, cturb, func, x0, x1, tol)
          x0 = x
          exit
       end if

       ! In case of failing to converge within itmax iterations stop at the minimum function

       if (iter > itmax) then
          f1 = func (icair, cturb, minx)
          x0 = minx
          exit
       end if

    end do

    root = x0

  end function hybrid

  !-----------------------------------------------------------------------
  function zbrent (msg, icair, cturb, func, xa, xb, tol) result(root)
    !
    ! !DESCRIPTION:
    ! Use Brent's method to find the root of a function, which is known to exist
    ! between xa and xb. The root is updated until its accuracy is tol.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*)                    :: msg            ! String to be printed
    integer             , intent(in)    :: icair          ! index of canopy air
    procedure (xfunc)                   :: func           ! Function to solve
    PetscReal           , intent(in)    :: xa, xb         ! Minimum and maximum of the variable domain to search
    PetscReal           , intent(in)    :: tol            ! Error tolerance
    type(canopy_turbulence_auxvar_type) , intent(inout) :: cturb
    !
    ! !LOCAL VARIABLES:
    integer, parameter                  :: itmax = 50     ! Maximum number of iterations
    PetscReal           , parameter     :: eps = 1.e-08 ! Relative error tolerance
    integer                             :: iter           ! Iteration loop index
    PetscReal                           :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    PetscReal                           :: root
    !---------------------------------------------------------------------

    a = xa
    b = xb
    fa = func (icair, cturb, a)
    fb = func (icair, cturb, b)

    if ((fa > 0.d0 .and. fb > 0.d0) .or. (fa < 0.d0 .and. fb < 0.d0)) then
       write (iulog,*) 'zbrent: Root must be bracketed'
       write (iulog,*) 'called from: ',msg
       write (iulog,*) xa, fa
       write (iulog,*) xb, fb
       call endrun()
    end if
    c = b
    fc = fb
    iter = 0
    do
       if (iter == itmax) exit
       iter = iter + 1
       if ((fb > 0.d0 .and. fc > 0.d0) .or. (fb < 0.d0 .and. fc < 0.d0)) then
          c = a
          fc = fa
          d = b - a
          e = d
       end if
       if (abs(fc) < abs(fb)) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
       end if
       tol1 = 2.d0 * eps * abs(b) + 0.5d0 * tol
       xm = 0.5d0 * (c - b)
       if (abs(xm) <= tol1 .or. fb == 0.d0) exit
       if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s = fb / fa
          if (a == c) then
             p = 2.d0 * xm * s
             q = 1.d0 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * (2.d0 * xm * q * (q - r) - (b-a) * (r - 1.d0))
             q = (q - 1.d0) * (r - 1.d0) * (s - 1.d0)
          end if
          if (p > 0.d0) q = -q
          p = abs(p)
          if (2.d0*p < min(3.d0*xm*q-abs(tol1*q),abs(e*q))) then
             e = d
             d = p / q
          else
             d = xm
             e = d
          end if
       else
          d = xm
          e = d
       end if
       a = b
       fa = fb
       if (abs(d) > tol1) then
          b = b + d
       else
          b = b + sign(tol1,xm)
       end if
       fb = func (icair, cturb, b)
       if (fb == 0.d0) exit
    end do
    root = b

    if (iter == itmax) then
       write (iulog,*) 'zbrent: Maximum number of interations exceeded'
       write (iulog,*) 'called from: ',msg
       call endrun()
    end if

  end function zbrent

  !-----------------------------------------------------------------------
  subroutine tridiag (a, b, c, r, u, n)
    !
    ! !DESCRIPTION:
    ! Solve a tridiagonal system of equations
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    integer,  intent(in)  :: n        ! Number of soil layers
    PetscReal, intent(in)  :: a(n)     ! A vector for tridiagonal solution
    PetscReal, intent(in)  :: b(n)     ! B vector for tridiagonal solution
    PetscReal, intent(in)  :: c(n)     ! C vector for tridiagonal solution
    PetscReal, intent(in)  :: r(n)     ! R vector for tridiagonal solution
    PetscReal, intent(out) :: u(n)     ! U vector for tridiagonal solution
    !
    ! !LOCAL VARIABLES:
    PetscReal :: gam(n)                ! Temporary calculation
    PetscReal :: bet                   ! Temporary calculation
    integer  :: j                     ! Soil layer index
    !---------------------------------------------------------------------

    ! Tridiagonal solution:
    !
    ! Solve for U given the set of equations F x U = R, where U is a vector
    ! of length N, R is a vector of length N, and F is an N x N tridiagonal
    ! matrix defined by the vectors A, B, C (each of length N). A(1) and
    ! C(N) are undefined and are not referenced by the subroutine.
    !
    !    | b(1) c(1)   0  ...                      |   | u(1)   |   | r(1)   |
    !    | a(2) b(2) c(2) ...                      |   | u(2)   |   | r(2)   |
    !    |                ...                      | x | ...    | = | ...    |
    !    |                ... a(n-1) b(n-1) c(n-1) |   | u(n-1) |   | r(n-1) |
    !    |                ...   0    a(n)   b(n)   |   | u(n)   |   | r(n)   |
    !

    bet = b(1)
    u(1) = r(1) / bet
    do j = 2, n
       gam(j) = c(j-1) / bet
       bet = b(j) - a(j)*gam(j)
       u(j) = (r(j) - a(j)*u(j-1)) / bet
    end do
    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do

  end subroutine tridiag

  !-----------------------------------------------------------------------
  function beta_function (p, q) result(beta)
    !
    ! !DESCRIPTION:
    ! Return the value of the beta function evaluated at p and q: B(p,q)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in) :: p      ! Input argument
    PetscReal, intent(in) :: q      ! Input argument
    PetscReal :: beta               ! Beta function: B(p,q)

    beta = exp(log_gamma_function(p) + log_gamma_function(q) - log_gamma_function(p+q))
  
  end function beta_function

  !-----------------------------------------------------------------------
  function log_gamma_function (x) result(gammaln)
    !
    ! !DESCRIPTION:
    ! Return the value of the log natural of the gamma function evaluated at x: ln(G(x)) 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in) :: x     ! Input argument
    PetscReal :: gammaln           ! ln(G(x))
    PetscReal :: y
    PetscReal :: tmp
    PetscReal :: ser
    integer  :: j

    PetscReal, parameter :: coef(6) = (/ 76.18009172947146e0, -86.50532032941677e0, &
         24.01409824083091e0, -1.231739572450155e0, 0.1208650973866179e-02, -0.5395239384953e-05 /)
    PetscReal, parameter :: stp = 2.5066282746310005d0

    y = x
    tmp = x + 5.5d0
    tmp = (x + 0.5d0) * log(tmp) - tmp
    ser = 1.000000000190015d0
    do j = 1, 6
       y = y + 1.d0
       ser = ser + coef(j) / y
    end do
    gammaln = tmp + log(stp*ser/x)

  end function log_gamma_function

#endif

end module MathToolsMod
