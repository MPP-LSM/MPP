
module EOSWaterMod

#ifdef USE_PETSC_LIB

  ! !USES:
  use mpp_varctl         , only : iulog
  use mpp_abortutils         , only : endrun
  use mpp_shr_log_mod        , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  PetscInt, parameter, public :: DENSITY_CONSTANT  = 1
  PetscInt, parameter, public :: DENSITY_TGDPB01   = 2

  public :: Density
  public :: Viscosity
  public :: InternalEnergyAndEnthalpy

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine Density(p, t_K, density_itype, den, dden_dp, dden_dT)
    !
    ! !DESCRIPTION:
    ! Given pressure, temperature and type of density formulation, compute:
    ! - density of water, and
    ! - first derivative of density w.r.t pressure
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : FMWH2O
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)   :: p
    PetscReal, intent(in)   :: t_K
    PetscInt, intent(in)    :: density_itype
    PetscReal, intent(out)  :: den
    PetscReal, intent(out)  :: dden_dp
    PetscReal, intent(out)  :: dden_dT

    select case(density_itype)
    case (DENSITY_CONSTANT)
       call DensityConstant(p, t_K, den, dden_dp, dden_dT)
    case (DENSITY_TGDPB01)
       call DensityTGDPB01 (p , t_K, den, dden_dp,dden_dT )
    case default
       write(iulog,*)'Density: Unknown denity_itype. '
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine Density

  !------------------------------------------------------------------------
  subroutine DensityConstant(p, t_K, den, dden_dp, dden_dT)
    !
    ! !DESCRIPTION:
    ! Return constant density of water
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : FMWH2O
    use mpp_varcon                , only : denh2o
    !
    implicit none
    !
    ! !ARGUMENTS    
    PetscReal, intent(in)  :: p         ! [Pa]
    PetscReal, intent(in)  :: t_K       ! [K]
    PetscReal, intent(out) :: den       ! [kmol m^{-3}]
    PetscReal, intent(out) :: dden_dp   ! [kmol m^{-3} Pa^{-1}]
    PetscReal, intent(out) :: dden_dT   ! [kmol m^{-3} K^{-1}]

    den     = denh2o/FMWH2O ! [kmol m^{-3}]
    dden_dp = 0.d0
    dden_dT = 0.d0

  end subroutine DensityConstant

  !------------------------------------------------------------------------
  subroutine DensityTGDPB01(p, t_K, den, dden_dp, dden_dT)
    !
    ! !DESCRIPTION:
    ! Return density and deriv. w.r.t. pressure based on Tanaka et al. (2001)
    !
    ! Reference:
    ! Tanaka M. , G. Girard, R. Davis, A. Peuto, and N. Bignell. 2001.
    ! Recommended table for the density of water between 0 °C
    ! and 40 °C based on recent experimental reports. Metrologia,
    ! 38:301-309 [doi:10.1088/0026-1394/38/4/3].
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : FMWH2O
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: p                ! [Pa]
    PetscReal, intent(in)  :: t_K              ! [K]
    PetscReal, intent(out) :: den              ! [kmol m^{-3}]
    PetscReal, intent(out) :: dden_dp          ! [kmol m^{-3} Pa^{-1}]
    PetscReal, intent(out) :: dden_dT          ! [kmol m^{-3} K^{-1}]

    !
    PetscReal,parameter    :: a1 = -3.983035d0     ! [degC]
    PetscReal,parameter    :: a2 = 301.797d0       ! [degC]
    PetscReal,parameter    :: a3 = 522528.9d0      ! [degC^{2}]
    PetscReal,parameter    :: a4 = 69.34881d0      ! [degC]
    PetscReal,parameter    :: a5 = 999.974950d0    ! [kg m^{-3}]
    PetscReal,parameter    :: k0 = 50.74d-11       ! [Pa^{-1}]
    PetscReal,parameter    :: k1 = -0.326d-11      ! [Pa^{-1} degC^{-1}]
    PetscReal,parameter    :: k2 = 0.00416d-11     ! [Pa^{-1} degC^{-2}]
    PetscReal,parameter    :: p0 = 101325.d0       ! [Pa]
    PetscReal              :: t_c
    PetscReal              :: dent
    PetscReal              :: kappa
    PetscReal              :: ddent_dt
    PetscReal              :: ddent_dt_1
    PetscReal              :: ddent_dt_2
    PetscReal              :: ddent_dt_3
    PetscReal              :: ddent_dp
    PetscReal              :: dkappa_dp
    PetscReal              :: dkappa_dt

    t_c = t_K - 273.15d0

    ! Density of water as function of t_K
    dent = a5*(1.d0 - ((t_c + a1)**2.d0)*(t_c + a2)/a3/(t_c + a4))

    ! Compressibility of water
    if (p > p0) then
       kappa = (1.d0 + (k0 + k1*t_c + k2*t_c**2.d0)*(p - p0))
    else
       kappa = 1.d0
    endif

    ! Density of water
    den = dent*kappa/FMWH2O ! [kmol m^{-3}]

    ! Derivative
    ddent_dp   = 0.d0
    ddent_dt_1 = -((t_c + a1)**2.d0)/a3/(t_c + a4)
    ddent_dt_2 = -2.d0*(t_c + a1)*(t_c + a2)/a3/(t_c + a4)
    ddent_dt_3 =  ((t_c + a1)**2.d0)*(t_c + a2)/a3/((t_c + a4)**2.d0)
    ddent_dt   = a5*(ddent_dt_1 + ddent_dt_2 + ddent_dt_3)

    if (p > p0) then
       dkappa_dp  = (k0 + k1*t_c + k2*t_c**2.d0)
       dkappa_dt  = (k1 + 2.d0*k2*t_c)*(p - p0)
    else
       dkappa_dp  = 0.d0
       dkappa_dt  = 0.d0
    endif

    dden_dT    = (ddent_dt*kappa + dent*dkappa_dt)/FMWH2O
    dden_dp    = (ddent_dp*kappa + dent*dkappa_dp)/FMWH2O

  end subroutine DensityTGDPB01

  !------------------------------------------------------------------------
  subroutine Viscosity(p, t_K, vis, dvis_dp)
    !
    ! !DESCRIPTION:
    ! Return viscosity of water
    !
    implicit none
    !
    ! !ARGUMENTS    
    PetscReal, intent(in)  :: p
    PetscReal, intent(in)  :: t_K
    PetscReal, intent(out) :: vis
    PetscReal, intent(out) :: dvis_dp

    vis     = 8.904156d-4 ! [Pa s]
    dvis_dp = 0.d0

  end subroutine Viscosity

  !------------------------------------------------------------------------
  subroutine InternalEnergyAndEnthalpy(p, t_K, den, dden_dT, &
       u, h, du_dT, dh_dT)
    !
    ! !DESCRIPTION:
    ! Return internal energy of water
    !
    use MultiPhysicsProbConstants, only : FMWH2O
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: p        ! [Pa]
    PetscReal, intent(in)  :: t_K      ! [K]
    PetscReal, intent(in)  :: den      ! [kg m^{-3}]
    PetscReal, intent(in)  :: dden_dT  ! [kg m^{-3} K^{-1}]
    PetscReal, intent(out) :: u        ! [J kmol^{-1}]
    PetscReal, intent(out) :: h        ! [J kmol^{-1}]
    PetscReal, intent(out) :: du_dT    ! [J kmol^{-1} K^{-1}]
    PetscReal, intent(out) :: dh_dT    ! [J kmol^{-1} K^{-1}]
    !
    PetscReal, parameter :: u0 = 4.217 * 1.d3 ! [J kg^{-1} K^{-1}]

    u     = u0 * (t_K - 273.15d0)         ! [J kg^{-1}]
    du_dT = u0                            ! [J kg^{-1} K^{-1}]

    h     = u + P/den                     ! [J kg^{-1}]
    dh_dT = du_dT - P/(den**2.d0)*dden_dT ! [J kg^{-1} K^{-1}]

    u     = u * FMWH2O                    ! [J kmol^{-1}]
    h     = h * FMWH2O                    ! [J kmol^{-1}]
    du_dT = du_dT * FMWH2O                ! [J kmol^{-1} K^{-1}]
    dh_dT = dh_dT * FMWH2O                ! [J kmol^{-1} K^{-1}]

  end subroutine InternalEnergyAndEnthalpy

#endif

end module EOSWaterMod
