module CanopyTurbulence

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl              , only : iulog
  use mpp_abortutils          , only : endrun
  use mpp_shr_log_mod         , only : errMsg => shr_log_errMsg
  use CanopyTurbulenceAuxType , only : canopy_turbulence_auxvar_type
  use MathToolsMod            , only : hybrid

  implicit none
  private

  public :: ObukhovLength
  public :: WindProfile
  public :: AerodynamicConductances

  PetscReal, parameter :: cd = 0.25d0               ! RSL - leaf drag coefficient (dimensionless)
  PetscReal, parameter :: beta_neutral_max = 0.35d0 ! RSL - maximum value for beta in neutral conditions
  PetscReal, parameter :: cr = 0.3d0                ! RSL - parameter to calculate beta_neutral
  PetscReal, parameter :: c2 = 0.5d0                ! RSL - depth scale multiplier
  PetscReal, parameter :: Pr0 = 0.5d0               ! RSL - neutral value for Pr (Sc)
  PetscReal, parameter :: Pr1 = 0.3d0               ! RSL - magnitude of variation of Pr (Sc) with stability
  PetscReal, parameter :: Pr2 = 2.0d0               ! RSL - scale of variation of Pr (Sc) with stability
  PetscReal, parameter :: z0mg = 0.01d0             ! RSL - roughness length of ground (m)
  PetscReal, parameter :: wind_forc_min = 1.0d0     ! Minimum wind speed at forcing height (m/s)
  PetscReal, parameter :: eta_max = 20.d0           ! Maximum value for "eta" parameter (used to constrain lm/beta)
  PetscReal, parameter :: zeta_min = -2.0d0         ! Minimum value for Monin-Obukhov zeta parameter
  PetscReal, parameter :: zeta_max = 1.0d0          ! Maximum value for Monin-Obukhov zeta parameter
  PetscReal, parameter :: beta_min = 0.2d0          ! Minimum value for H&F beta parameter
  PetscReal, parameter :: beta_max = 0.5d0          ! Maximum value for H&F beta parameter
  PetscReal, parameter :: wind_min = 0.1d0          ! Minimum wind speed in canopy (m/s)
  PetscReal, parameter :: ra_max = 500.d0           ! Maximum aerodynamic resistance (s/m)

contains

  !------------------------------------------------------------------------
  subroutine ObukhovLength (cturb)
    !
    implicit none
    !
    class(canopy_turbulence_auxvar_type) :: cturb
    !
    PetscReal :: obu0, obu1, tol, dummy
    PetscInt  :: icair

    do icair = 1, cturb%ncair

       cturb%Lc(icair) = cturb%hc(icair)/(cd * cturb%pai(icair))

       ! Calculate Obukhov 
       obu0 = 100.d0          ! Initial estimate for Obukhov length (m)
       obu1 = -100.d0         ! Initial estimate for Obukhov length (m)
       tol  = 0.01d0           ! Accuracy tolerance for Obukhov length (m)
       dummy = hybrid ('ObukhovFunction', icair, cturb, Obukhov, obu0, obu1, tol)
       cturb%obu(icair) = cturb%obu_ustar(icair)

    end do

  end subroutine ObukhovLength

  !------------------------------------------------------------------------
  subroutine WindProfile (cturb)
    !
    use MultiPhysicsProbConstants, only : VKC, GRAVITY_CONSTANT
    !
    implicit none
    !
    class(canopy_turbulence_auxvar_type) :: cturb
    !
    PetscReal :: lm, lm_over_beta
    PetscInt  :: icair, k
    PetscReal :: dp, z_minus_d, h_minus_d
    PetscReal :: psi_m_hc, psi_m_rsl_hc, psi_m_zs, psi_m_rsl_zs, psim

    do icair = 1, cturb%ncair

       ! Above-canopy
       h_minus_d = cturb%hc(icair) - cturb%disp(icair)

       psi_m_hc     = psim_monin_obukhov (h_minus_d / cturb%obu(icair) );
       psi_m_rsl_hc = 0.d0

       do k = cturb%ntop(icair)+1, cturb%ncair
          z_minus_d = cturb%zs(icair, k) - cturb%disp(icair)

          psi_m_zs     = psim_monin_obukhov (z_minus_d / cturb%obu(icair) );
          psi_m_rsl_zs = 0.d0

          psim = -psi_m_zs + psi_m_hc + psi_m_rsl_zs - psi_m_rsl_hc + VKC/cturb%beta(icair)
          cturb%wind(icair,k) = cturb%ustar(icair) / VKC * (log(z_minus_d/h_minus_d) + psim)
       end do

       ! At top of the canopy
       cturb%ucan(icair) = cturb%ustar(icair) / cturb%beta(icair)

       ! Wind profile within canopy
       lm = 2.d0 * cturb%beta(icair)**3.d0 * cturb%Lc(icair)
       lm_over_beta = lm/cturb%beta(icair)

       do k = 2, cturb%ntop(icair)
          cturb%wind(icair,k) = cturb%ucan(icair) * exp((cturb%zs(icair,k) - cturb%hc(icair))/lm_over_beta)
          cturb%wind(icair,k) = max(cturb%wind(icair,k), 0.1d0)
       end do

       ! At ground
       cturb%wind(icair,1) = 0.d0
    end do

  end subroutine WindProfile

  !------------------------------------------------------------------------
  subroutine AerodynamicConductances (cturb)
    !
    use MultiPhysicsProbConstants, only : VKC, GRAVITY_CONSTANT
    !
    implicit none
    !
    class(canopy_turbulence_auxvar_type) :: cturb
    !
    PetscReal :: lm, lm_over_beta, sumres
    PetscInt  :: icair, k
    PetscReal :: z_minus_d, h_minus_d, zlog
    PetscReal :: zl, zu, res
    PetscReal :: psi_c_hc, psi_c_rsl_hc, psi_c_zs, psi_c_rsl_zs
    PetscReal :: psic, psic1, psic2
    PetscReal :: ga_below_hc, ga_above_hc
    PetscReal :: zoc_g, zlog_m, zlog_c, ustar_g

    do icair = 1, cturb%ncair

       ! Above-canopy
       h_minus_d = cturb%hc(icair) - cturb%disp(icair)

       psi_c_hc     = psic_monin_obukhov (h_minus_d / cturb%obu(icair) );
       psi_c_rsl_hc = 0.d0

       do k = cturb%ntop(icair)+1, cturb%ncan_lev-1

          ! Lower height
          z_minus_d = cturb%zs(icair, k) - cturb%disp(icair)
          psi_c_zs     = psic_monin_obukhov (z_minus_d / cturb%obu(icair) );
          psi_c_rsl_zs = 0.d0
          psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

          ! Upper height
          z_minus_d = cturb%zs(icair, k+1) - cturb%disp(icair)
          psi_c_zs     = psic_monin_obukhov (z_minus_d / cturb%obu(icair) );
          psi_c_rsl_zs = 0.d0
          psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

          psic = psic2 - psic1
          zlog = log((cturb%zs(icair,k+1) - cturb%disp(icair))/ (cturb%zs(icair,k) - cturb%disp(icair)))
          cturb%ga_prof(icair,k) = cturb%rhomol(icair) * VKC * cturb%ustar(icair) / (zlog + psic)
       end do

       ! At top to the reference height
       k = cturb%ncan_lev
       z_minus_d = cturb%zs(icair, k) - cturb%disp(icair)
       psi_c_zs     = psic_monin_obukhov (z_minus_d / cturb%obu(icair) );
       psi_c_rsl_zs = 0.d0
       psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

       z_minus_d = cturb%zref(icair) - cturb%disp(icair)
       psi_c_zs     = psic_monin_obukhov (z_minus_d / cturb%obu(icair) );
       psi_c_rsl_zs = 0.d0
       psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

       psic = psic2 - psic1
       zlog = log((cturb%zref(icair) - cturb%disp(icair))/ (cturb%zs(icair,k) - cturb%disp(icair)))
       cturb%ga_prof(icair,k) = cturb%rhomol(icair) * VKC * cturb%ustar(icair) / (zlog + psic)

       ! Within canopy
       lm = 2.d0 * cturb%beta(icair)**3.d0 * cturb%Lc(icair)
       lm_over_beta = lm/cturb%beta(icair)

       do k = 2, cturb%ntop(icair)-1
          zl = cturb%zs(icair,k  ) - cturb%hc(icair)
          zu = cturb%zs(icair,k+1) - cturb%hc(icair)
          res = cturb%PrSc(icair) / ( cturb%beta(icair) * cturb%ustar(icair)) * &
               (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
          cturb%ga_prof(icair,k) = cturb%rhomol(icair)/res;
       end do

       ! Top of canopy layer
       k = cturb%ntop(icair)
       zl = cturb%zs(icair,k) - cturb%hc(icair)
       zu = cturb%hc(icair)   - cturb%hc(icair)
       res = cturb%PrSc(icair) / ( cturb%beta(icair) * cturb%ustar(icair)) * &
            (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
       ga_below_hc = cturb%rhomol(icair)/res;

       z_minus_d = cturb%hc(icair) - cturb%disp(icair)
       psi_c_zs     = psic_monin_obukhov (z_minus_d / cturb%obu(icair) );
       psi_c_rsl_zs = 0.d0
       psic1 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

       z_minus_d = cturb%zs(icair, k+1) - cturb%disp(icair)
       psi_c_zs     = psic_monin_obukhov (z_minus_d / cturb%obu(icair) );
       psi_c_rsl_zs = 0.d0
       psic2 = -psi_c_zs + psi_c_hc + psi_c_rsl_zs - psi_c_rsl_hc

       psic = psic2 - psic1
       zlog = log((cturb%zs(icair,k+1) - cturb%disp(icair))/ (cturb%hc(icair) - cturb%disp(icair)))
       ga_above_hc = cturb%rhomol(icair) * VKC * cturb%ustar(icair) / (zlog + psic)

       cturb%ga_prof(icair,k) = 1.d0 / (1.d0/ga_below_hc + 1.d0/ga_above_hc)

       ! Check aerodynamic resistance sums to 1/gac
       sumres = 1.d0/ga_above_hc
       do k = cturb%ntop(icair)+1, cturb%ncan_lev
          sumres = sumres + 1.d0/cturb%ga_prof(icair,k)
       end do

       if (abs(1.d0/sumres - cturb%gac(icair)) > 1e-06) then
          write(iulog,*)'ERROR: *** Above canopy conductances does not sums to 1/gac ***'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       ! At ground
       k = 1
       zoc_g = 0.1d0 * z0mg
       zlog_m = log(cturb%zs(icair,k+1)/z0mg)
       zlog_c = log(cturb%zs(icair,k+1)/zoc_g)
       ustar_g = cturb%wind(icair,k+1) * VKC /zlog_m;
       ustar_g = max(ustar_g, 0.01d0)
       res = zlog_c/(VKC * ustar_g)
       cturb%ga_prof(icair,k) = cturb%rhomol(icair)/res

       !gac0(p) = rhomol(p) * vkc * ustar_g / zlog_c   !!! CLMml v0 CODE !!!

       res = min (cturb%rhomol(icair)/cturb%ga_prof(icair,k), ra_max)
       cturb%ga_prof(icair,k) = cturb%rhomol(icair) / res

       ! Limit resistance to < 500 s/m
       k = 1
       do k = 1+1, cturb%ncan_lev
          res = min(cturb%rhomol(icair)/cturb%ga_prof(icair,k), 500.d0)
          cturb%ga_prof(icair,k) = cturb%rhomol(icair)/res
       end do
    end do

  end subroutine AerodynamicConductances

  !------------------------------------------------------------------------
   function Obukhov (icair, cturb, obu_val) result(obu_dif)
     !
     use MultiPhysicsProbConstants, only : VKC, GRAVITY_CONSTANT
     use MultiPhysicsProbConstants, only : MM_H2O, MM_DRY_AIR
     !
     implicit none
     !
     class(canopy_turbulence_auxvar_type) :: cturb
     PetscInt :: icair
     PetscReal :: obu_val, obu_dif
     !
     PetscReal :: beta_neutral, LcL, beta
     PetscReal :: a, b, c, d, q, r
     PetscReal :: dp, z_minus_d, h_minus_d
     PetscReal :: phi_m_hc
     PetscReal :: phi_c_hc
     PetscReal :: psim, psi_m_zref, psi_m_hc, psi_m_rsl_zref, psi_m_rsl_hc
     PetscReal :: psic, psi_c_zref, psi_c_hc, psi_c_rsl_zref, psi_c_rsl_hc
     PetscReal :: zlog, tvstar
     PetscReal :: c1

     ! Limit the value of Obukhov length
     if (abs(obu_val) < 0.1d0) obu_val = 0.1d0

     ! neutral value of beta
     c1 = (VKC / log((cturb%hc(icair) + z0mg)/z0mg))**2.d0
     beta_neutral = min(sqrt(c1 + cr*cturb%pai(icair)), beta_neutral_max)

     LcL = cturb%Lc(icair)/obu_val

     if (LcL <= 0.d0) then

        b = 16.d0 * LcL * beta_neutral**4.d0
        beta = sqrt( 0.5d0 * (-b + sqrt( b**2.d0 + 4.d0 * beta_neutral**4.d0)) )
     else
        a = 5.d0 * LcL
        b = 0.d0
        c = 1.d0
        d = -beta_neutral
        q = (2.d0*b**3.d0 - 9.d0*a*b*c + 27.d0*(a**2)*d)**2.d0 - 4.d0*(b**2.d0 - 3.d0*a*c)**3.d0
        q = sqrt(q)
        r = 0.5d0 * (q + 2.d0*b**3.d0 - 9.d0*a*b*c + 27.d0*(a**2)*d)
        r = r**(1.d0/3.d0)
        beta = -(b+r)/(3.d0*a) - (b**2.d0 - 3.d0*a*c)/(3.d0*a*r)
     end if

     beta = min(beta_max, max(beta, beta_min))
     cturb%beta(icair) = beta

     dp = beta**2.d0 * cturb%Lc(icair)
     dp = dp * (1.d0 - exp(-0.25d0*cturb%pai(icair)/beta**2.d0))
     dp = min(cturb%hc(icair), dp)
     cturb%disp(icair) = max(cturb%hc(icair) - dp, 0.d0)

     z_minus_d = cturb%zref(icair) - cturb%disp(icair) ! zdisp
     h_minus_d = cturb%hc(icair)   - cturb%disp(icair)

     ! Prandtl (Schmidt) number at canopy top
     cturb%PrSc(icair) = Pr0 + Pr1 * tanh(Pr2*cturb%Lc(icair)/obu_val);
     cturb%PrSc(icair) = (1.d0 - beta_neutral/beta_neutral_max) * 1.d0 + (beta_neutral/beta_neutral_max) * cturb%PrSc(icair)

     ! Compute Monin-Obokhov phi functions for momentum and scalar
     phi_m_hc = phim_monin_obukhov(h_minus_d / obu_val)
     phi_c_hc = phic_monin_obukhov(h_minus_d / obu_val)

     cturb%c2(icair)  = c2;
     cturb%c1m(icair) = (1.d0 - VKC/ (2.d0 * beta * phi_m_hc)) * exp(0.5d0 * cturb%c2(icair))

     ! Compute psi_hat
     psi_m_rsl_zref = 0.d0
     psi_m_rsl_hc   = 0.d0

     psi_c_rsl_zref = 0.d0
     psi_c_rsl_hc   = 0.d0

     ! Compute the Monin-Obukhov psi functions for momentum and scalar
     ! at the reference height and the canopy top
     psi_m_zref = psim_monin_obukhov (z_minus_d / obu_val);
     psi_m_hc   = psim_monin_obukhov (h_minus_d / obu_val);

     psi_c_zref = psic_monin_obukhov (z_minus_d / obu_val);
     psi_c_hc   = psic_monin_obukhov (h_minus_d / obu_val);

     ! Calculate u*, T*, q* and Tv*
     zlog = log(z_minus_d/h_minus_d);
     psim = -psi_m_zref + psi_m_hc + psi_m_rsl_zref - psi_m_rsl_hc + VKC/beta;
     psic = -psi_c_zref + psi_c_hc + psi_c_rsl_zref - psi_c_rsl_hc;

     cturb%ustar(icair)     =  cturb%uref(icair)                       * VKC / (zlog + psim)
     cturb%tstar(icair)     = (cturb%thref(icair) - cturb%tcan(icair)) * VKC / (zlog + psic)
     cturb%vstar(icair)     = (cturb%vref(icair)  - cturb%vcan(icair)) * VKC / (zlog + psic)
     cturb%obu_ustar(icair) = obu_val

     ! Aerodynamic conductance
     cturb%gac(icair) = cturb%rhomol(icair) * VKC * cturb%ustar(icair) / (zlog + psic)

     ! New Obukhov length
     tvstar = cturb%tstar(icair) + 0.61d0 * cturb%thref(icair) * cturb%vstar(icair) * MM_H2O/MM_DRY_AIR
     cturb%obu(icair) = cturb%ustar(icair)**2.d0 * cturb%thvref(icair)/ (VKC * GRAVITY_CONSTANT * tvstar)

     obu_dif = cturb%obu(icair) - obu_val

   end function Obukhov

  !-----------------------------------------------------------------------
  function phim_monin_obukhov (zeta) result(phi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov phi stability function for momentum
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    PetscReal :: phi               ! phi for momentum
    !---------------------------------------------------------------------

    if (zeta < 0.d0) then           ! --- unstable
       phi = 1.d0 / sqrt(sqrt(1.d0 - 16.d0 * zeta))
    else                             ! --- stable
       phi = 1.d0 + 5.d0 * zeta
    end if

  end function phim_monin_obukhov

  !-----------------------------------------------------------------------
  function phic_monin_obukhov (zeta) result(phi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov phi stability function for scalars
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    PetscReal :: phi               ! phi for scalars
    !---------------------------------------------------------------------

    if (zeta < 0.d0) then           ! --- unstable
       phi = 1.d0 / sqrt(1.d0 - 16.d0 * zeta)
    else                             ! --- stable
       phi = 1.d0 + 5.d0 * zeta
    end if

  end function phic_monin_obukhov

  !-----------------------------------------------------------------------
  function psim_monin_obukhov (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov psi stability function for momentum
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    PetscReal :: x                 ! (1 - 16*zeta)**1/4
    PetscReal :: psi               ! psi for momentum
    PetscReal, parameter :: pi  = 4 * atan (1.0d0)
    !---------------------------------------------------------------------

    if (zeta < 0.d0) then           ! --- unstable
       x = sqrt(sqrt(1.d0 - 16.d0 * zeta))
       psi = 2.d0 * log((1.d0+x)/2.d0) + log((1.d0+x*x)/2.d0) - 2.d0*atan(x) + pi * 0.5d0
    else                             ! --- stable
       psi = -5.d0 * zeta
    end if

  end function psim_monin_obukhov

  !-----------------------------------------------------------------------
  function psic_monin_obukhov (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov psi stability function for scalars
    !
    ! !ARGUMENTS:
    implicit none
    PetscReal, intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    PetscReal :: x                 ! (1 - 16*zeta)**1/4
    PetscReal :: psi               ! psi for scalars
    !---------------------------------------------------------------------

    if (zeta < 0.d0) then           ! --- unstable
       x = sqrt(sqrt(1.d0 - 16.d0 * zeta))
       psi = 2.d0 * log((1.d0+x*x)/2.d0)
    else                             ! --- stable
       psi = -5.d0 * zeta
    end if

  end function psic_monin_obukhov

#endif

end module CanopyTurbulence
