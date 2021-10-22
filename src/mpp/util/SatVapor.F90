module WaterVaporMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate saturation vapor pressure and latent heat of vaporization
  !
  ! !USES:
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SatVap     ! Saturation vapor pressure and derivative
  public :: LatVap     ! Latent heat of vaporization
  !
  ! This code is from Gordan Bonan's (NCAR) repository at
  ! https://github.com/gbonan/CLM-ml_v0
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SatVap (t, es, desdt)
    !
    ! !DESCRIPTION:
    ! Compute saturation vapor pressure and change in saturation vapor pressure
    ! with respect to temperature. Polynomial approximations are from:
    ! Flatau et al (1992) Polynomial fits to saturation vapor pressure.
    ! Journal of Applied Meteorology 31:1507-1513
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : TFRZ
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in)  :: t        ! Temperature (K)
    PetscReal, intent(out) :: es       ! Vapor pressure (Pa)
    PetscReal, intent(out) :: desdt    ! d(es)/d(t) (Pa/K)
    !
    ! !LOCAL VARIABLES:
    PetscReal :: tc                    ! Temperature (C)
    !---------------------------------------------------------------------

    ! For water vapor (temperature range is 0C to 100C)
 
    PetscReal, parameter :: a0 =  6.11213476d0
    PetscReal, parameter :: a1 =  0.444007856d0
    PetscReal, parameter :: a2 =  0.143064234d-01
    PetscReal, parameter :: a3 =  0.264461437d-03
    PetscReal, parameter :: a4 =  0.305903558d-05
    PetscReal, parameter :: a5 =  0.196237241d-07
    PetscReal, parameter :: a6 =  0.892344772d-10
    PetscReal, parameter :: a7 = -0.373208410d-12
    PetscReal, parameter :: a8 =  0.209339997d-15
 
    ! and for derivative
 
    PetscReal, parameter :: b0 =  0.444017302d0
    PetscReal, parameter :: b1 =  0.286064092d-01
    PetscReal, parameter :: b2 =  0.794683137d-03
    PetscReal, parameter :: b3 =  0.121211669d-04
    PetscReal, parameter :: b4 =  0.103354611d-06
    PetscReal, parameter :: b5 =  0.404125005d-09
    PetscReal, parameter :: b6 = -0.788037859d-12
    PetscReal, parameter :: b7 = -0.114596802d-13
    PetscReal, parameter :: b8 =  0.381294516d-16
 
    ! For ice (temperature range is -75C to 0C)
 
    PetscReal, parameter :: c0 =  6.11123516d0
    PetscReal, parameter :: c1 =  0.503109514d0
    PetscReal, parameter :: c2 =  0.188369801d-01
    PetscReal, parameter :: c3 =  0.420547422d-03
    PetscReal, parameter :: c4 =  0.614396778d-05
    PetscReal, parameter :: c5 =  0.602780717d-07
    PetscReal, parameter :: c6 =  0.387940929d-09
    PetscReal, parameter :: c7 =  0.149436277d-11
    PetscReal, parameter :: c8 =  0.262655803d-14
 
    ! and for derivative
 
    PetscReal, parameter :: d0 =  0.503277922d0
    PetscReal, parameter :: d1 =  0.377289173d-01
    PetscReal, parameter :: d2 =  0.126801703d-02
    PetscReal, parameter :: d3 =  0.249468427d-04
    PetscReal, parameter :: d4 =  0.313703411d-06
    PetscReal, parameter :: d5 =  0.257180651d-08
    PetscReal, parameter :: d6 =  0.133268878d-10
    PetscReal, parameter :: d7 =  0.394116744d-13
    PetscReal, parameter :: d8 =  0.498070196d-16

    tc = t - tfrz
    if (tc > 100.0d0) tc = 100.0d0
    if (tc < -75.0d0) tc = -75.0d0

    if (tc >= 0.0d0) then
       es    = a0 + tc*(a1 + tc*(a2 + tc*(a3 + tc*(a4 &
             + tc*(a5 + tc*(a6 + tc*(a7 + tc*a8)))))))
       desdt = b0 + tc*(b1 + tc*(b2 + tc*(b3 + tc*(b4 &
             + tc*(b5 + tc*(b6 + tc*(b7 + tc*b8)))))))
    else
       es    = c0 + tc*(c1 + tc*(c2 + tc*(c3 + tc*(c4 &
             + tc*(c5 + tc*(c6 + tc*(c7 + tc*c8)))))))
       desdt = d0 + tc*(d1 + tc*(d2 + tc*(d3 + tc*(d4 &
             + tc*(d5 + tc*(d6 + tc*(d7 + tc*d8)))))))
    end if

    es    = es    * 100.d0            ! Convert from mb to Pa
    desdt = desdt * 100.d0            ! Convert from mb to Pa

  end subroutine SatVap

  !-----------------------------------------------------------------------
  function LatVap (t) result(lambda)
    !
    ! !DESCRIPTION:
    ! Latent heat of vaporization in relation to air temperature
    !
    ! !USES:
    use MultiPhysicsProbConstants, only: TFRZ, HVAP, HSUB, MM_H2O
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in)  :: t     ! Temperature (K)
    !
    ! !LOCAL VARIABLES:
    PetscReal :: lambda             ! Molar latent heat of vaporization (J/mol)
    !---------------------------------------------------------------------

    if (t > TFRZ) then
       lambda = HVAP
    else
       lambda = HSUB
    endif
    lambda = lambda * MM_H2O ! Molar latent heat of vaporization (J/mol)

  end function LatVap

#endif

end module WaterVaporMod
