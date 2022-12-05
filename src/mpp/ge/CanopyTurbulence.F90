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
    PetscReal :: psim, psic

    do icair = 1, cturb%ncair

       ! Above-canopy
       h_minus_d = cturb%hc(icair) - cturb%disp(icair)

       do k = cturb%ntop(icair)+1, cturb%ncan_lev
          z_minus_d = cturb%zs(icair, k) - cturb%disp(icair)

          call ComputePsiRSL (cturb%zs(icair,k), cturb%hc(icair), cturb%disp(icair), &
               cturb%obu(icair), cturb%beta(icair), cturb%PrSc(icair), psim, psic)

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
    PetscReal :: psic, psic1, psic2, psim1, psim2
    PetscReal :: ga_below_hc, ga_above_hc
    PetscReal :: zoc_g, zlog_m, zlog_c, ustar_g

    associate ( &
         hc       => cturb%hc       , &
         zs       => cturb%zs       , &
         obu      => cturb%obu      , &
         ntop     => cturb%ntop     , &
         disp     => cturb%disp     , &
         beta     => cturb%beta     , &
         ga_prof  => cturb%ga_prof  , &
         zref     => cturb%zref     , &
         Lc       => cturb%Lc       , &
         rhomol   => cturb%rhomol   , &
         ustar    => cturb%ustar    , &
         wind     => cturb%wind     , &
         PrSc     => cturb%PrSc       &
         )

      do icair = 1, cturb%ncair

         ! Above-canopy
         h_minus_d = hc(icair) - disp(icair)

         do k = ntop(icair)+1, cturb%ncan_lev-1

            ! Lower height
            call ComputePsiRSL (zs(icair,k  ), hc(icair), disp(icair), obu(icair), beta(icair), PrSc(icair), psim1, psic1)
            call ComputePsiRSL (zs(icair,k+1), hc(icair), disp(icair), obu(icair), beta(icair), PrSc(icair), psim2, psic2)

            psic = psic2 - psic1
            zlog = log((zs(icair,k+1) - disp(icair))/ (zs(icair,k) - disp(icair)))
            ga_prof(icair,k) = rhomol(icair) * VKC * ustar(icair) / (zlog + psic)
         end do

         ! At top to the reference height
         k = cturb%ncan_lev
         call ComputePsiRSL (zs(icair,k), hc(icair), disp(icair), obu(icair), beta(icair), PrSc(icair), psim1, psic1)
         call ComputePsiRSL (zref(icair), hc(icair), disp(icair), obu(icair), beta(icair), PrSc(icair), psim2, psic2)

         psic = psic2 - psic1
         zlog = log((zref(icair) - disp(icair))/ (zs(icair,k) - disp(icair)))
         ga_prof(icair,k) = rhomol(icair) * VKC * ustar(icair) / (zlog + psic)

         ! Within canopy
         lm = 2.d0 * beta(icair)**3.d0 * Lc(icair)
         lm_over_beta = lm/beta(icair)

         do k = 2, ntop(icair)-1
            zl = zs(icair,k  ) - hc(icair)
            zu = zs(icair,k+1) - hc(icair)
            res = PrSc(icair) / ( beta(icair) * ustar(icair)) * &
                 (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
            ga_prof(icair,k) = rhomol(icair)/res;
         end do

         ! Top of canopy layer
         k = ntop(icair)
         zl = zs(icair,k) - hc(icair)
         zu = hc(icair)   - hc(icair)
         res = PrSc(icair) / ( beta(icair) * ustar(icair)) * &
              (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
         ga_below_hc = rhomol(icair)/res;

         call ComputePsiRSL (hc(icair    ), hc(icair), disp(icair), obu(icair), beta(icair), PrSc(icair), psim1, psic1)
         call ComputePsiRSL (zs(icair,k+1), hc(icair), disp(icair), obu(icair), beta(icair), PrSc(icair), psim2, psic2)

         psic = psic2 - psic1
         zlog = log((zs(icair,k+1) - disp(icair))/ (hc(icair) - disp(icair)))
         ga_above_hc = rhomol(icair) * VKC * ustar(icair) / (zlog + psic)

         ga_prof(icair,k) = 1.d0 / (1.d0/ga_below_hc + 1.d0/ga_above_hc)

         ! Check aerodynamic resistance sums to 1/gac
         sumres = 1.d0/ga_above_hc
         do k = ntop(icair)+1, cturb%ncan_lev
            sumres = sumres + 1.d0/ga_prof(icair,k)
         end do

         if (abs(1.d0/sumres - cturb%gac(icair)) > 1e-06) then
            write(iulog,*)'ERROR: *** Above canopy conductances does not sums to 1/gac ***'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

         ! At ground
         k = 1
         zoc_g = 0.1d0 * z0mg
         zlog_m = log(zs(icair,k+1)/z0mg)
         zlog_c = log(zs(icair,k+1)/zoc_g)
         ustar_g = wind(icair,k+1) * VKC /zlog_m;
         ustar_g = max(ustar_g, 0.01d0)
         res = zlog_c/(VKC * ustar_g)
         ga_prof(icair,k) = rhomol(icair)/res

         !gac0(p) = rhomol(p) * vkc * ustar_g / zlog_c   !!! CLMml v0 CODE !!!

         res = min (rhomol(icair)/ga_prof(icair,k), ra_max)
         ga_prof(icair,k) = rhomol(icair) / res

         ! Limit resistance to < 500 s/m
         k = 1
         do k = 1+1, cturb%ncan_lev
            res = min(rhomol(icair)/ga_prof(icair,k), 500.d0)
            ga_prof(icair,k) = rhomol(icair)/res
         end do
      end do

    end associate

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
     PetscReal :: psim, psic
     PetscReal :: zlog, tvstar
     PetscReal :: c1, zeta, obu_cur

     obu_cur = obu_val
     ! Limit the value of Obukhov length
     if (abs(obu_cur) < 0.1d0) obu_cur = 0.1d0

     ! neutral value of beta
     c1 = (VKC / log((cturb%hc(icair) + z0mg)/z0mg))**2.d0
     beta_neutral = min(sqrt(c1 + cr*cturb%pai(icair)), beta_neutral_max)

     LcL = cturb%Lc(icair)/obu_cur

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
     cturb%PrSc(icair) = Pr0 + Pr1 * tanh(Pr2*cturb%Lc(icair)/obu_cur);
     cturb%PrSc(icair) = (1.d0 - beta_neutral/beta_neutral_max) * 1.d0 + (beta_neutral/beta_neutral_max) * cturb%PrSc(icair)

     zeta = (cturb%zref(icair) - cturb%disp(icair)) / obu_cur
     if (zeta >= 0.d0) then
        zeta = min(zeta_max, max(zeta,0.01d0))
     else
        zeta = max(zeta_min, min(zeta,-0.01d0))
     end if
     obu_cur = (cturb%zref(icair) - cturb%disp(icair)) / zeta

     ! Compute Monin-Obokhov phi functions for momentum and scalar
     phi_m_hc = phim_monin_obukhov(h_minus_d / obu_cur)
     phi_c_hc = phic_monin_obukhov(h_minus_d / obu_cur)

     cturb%c2(icair)  = c2;
     cturb%c1m(icair) = (1.d0 - VKC/ (2.d0 * beta * phi_m_hc)) * exp(0.5d0 * cturb%c2(icair))

     ! Compute psi terms in eqn. 6.78 and 6.101 in Bonan's book
     call ComputePsiRSL(cturb%zref(icair), cturb%hc(icair), cturb%disp(icair), obu_cur, &
          cturb%beta(icair), cturb%PrSc(icair), psim, psic)

     ! Calculate u*, T*, q* and Tv*
     zlog = log(z_minus_d/h_minus_d);
     cturb%ustar(icair)     =  cturb%uref(icair)                       * VKC / (zlog + psim)
     cturb%tstar(icair)     = (cturb%thref(icair) - cturb%tcan(icair)) * VKC / (zlog + psic)
     cturb%qstar(icair)     = (cturb%qref(icair)  - cturb%qcan(icair)) * VKC / (zlog + psic)
     cturb%obu_ustar(icair) = obu_cur

     ! Aerodynamic conductance
     cturb%gac(icair) = cturb%rhomol(icair) * VKC * cturb%ustar(icair) / (zlog + psic)

     ! New Obukhov length
     tvstar = cturb%tstar(icair) + 0.61d0 * cturb%thref(icair) * cturb%qstar(icair)
     cturb%obu(icair) = cturb%ustar(icair)**2.d0 * cturb%thvref(icair)/ (VKC * GRAVITY_CONSTANT * tvstar)

     obu_dif = cturb%obu(icair) - obu_val

   end function Obukhov

   !-----------------------------------------------------------------------
   subroutine ComputePsiRSL (za, hc, disp, obu, beta, PrSc, psim, psic)
     !
     use MultiPhysicsProbConstants, only : VKC, GRAVITY_CONSTANT
     use RSLPsiHat
     !
     implicit none
     !
     PetscReal, intent(in) :: za, hc, disp, obu, beta, PrSc
     PetscReal, intent(out) :: psim, psic
     !
     PetscReal :: c1, phim, phic
     PetscReal :: psihat1, psihat2
     PetscReal :: z_minus_d, h_minus_d
     PetscReal :: psi_m_zref     , psi_m_hc
     PetscReal :: psi_m_rsl_zref , psi_m_rsl_hc
     PetscReal :: psi_c_zref     , psi_c_hc
     PetscReal :: psi_c_rsl_zref , psi_c_rsl_hc

     z_minus_d = za - disp
     h_minus_d = hc - disp

     !
     ! Computation for momentum
     !

     ! Compute the Monin-Obukhov psi for momentum
     ! at the reference height and the canopy top
     phim = phim_monin_obukhov(h_minus_d/obu)
     c1 = (1.d0 - VKC/(2.d0 * beta * phim)) * exp(0.5d0*c2)

     psi_m_zref = psim_monin_obukhov (z_minus_d / obu);
     psi_m_hc   = psim_monin_obukhov (h_minus_d / obu);

     ! Compute the RSL psi_hat for momentum
     ! at the reference height and the canopy top
     call LookupPsihat((za - hc)/h_minus_d, h_minus_d/obu, zdtgridM, dtLgridM, psigridM, psihat1)
     call LookupPsihat((hc - hc)/h_minus_d, h_minus_d/obu, zdtgridM, dtLgridM, psigridM, psihat2)

     psi_m_rsl_zref = psihat1 * c1
     psi_m_rsl_hc   = psihat2 * c1

     psim = -psi_m_zref + psi_m_hc + psi_m_rsl_zref - psi_m_rsl_hc + VKC/beta;

     !
     ! Computation for scalar
     !

     ! Compute the Monin-Obukhov psi for scalar
     ! at the reference height and the canopy top
     phic = phic_monin_obukhov(h_minus_d/obu)
     c1 = (1.d0 - PrSc*VKC/(2.d0 * beta * phic)) * exp(0.5d0*c2)

     psi_c_zref = psic_monin_obukhov (z_minus_d / obu);
     psi_c_hc   = psic_monin_obukhov (h_minus_d / obu);

     ! Compute the RSL psi_hat for momentum
     ! at the reference height and the canopy top
     call LookupPsihat((za - hc)/h_minus_d, h_minus_d/obu, zdtgridH, dtLgridH, psigridH, psihat1)
     call LookupPsihat((hc - hc)/h_minus_d, h_minus_d/obu, zdtgridH, dtLgridH, psigridH, psihat2)
     psi_c_rsl_zref = psihat1 * c1
     psi_c_rsl_hc   = psihat2 * c1

     psic = -psi_c_zref + psi_c_hc + psi_c_rsl_zref - psi_c_rsl_hc;

   end subroutine ComputePsiRSL

   !-----------------------------------------------------------------------
   subroutine LookupPsihat (zdt, dtL, zdtgrid, dtLgrid, psigrid, psihat)
     !
     use RSLPsiHat, only : nZ, nL
     !
     ! !ARGUMENTS:
     implicit none
     !
     PetscReal, intent(in) :: zdt, dtL
     PetscReal, intent(in) :: zdtgrid(nZ,1), dtLgrid(1,nL), psigrid(nZ,nL)
     PetscReal, intent(out) :: psihat
     !
     !LOCAL VARIABLES
     PetscInt  :: ii, jj
     PetscInt  :: L1, L2, Z1, Z2
     PetscReal :: wL1, wL2, wZ1, wZ2
     !---------------------------------------------------------------------


     !  dtLgrid
     !  /|\
     !   |
     !   |
     !   |
     !   |
     !   |
     !   |
     !   |
     !   |
     !   ---------------------->
     !                      zdtgrid

     L1 = 0 ; L2 = 0
     if (dtL <= dtLgrid(1,1)) then
        L1 = 1
        L2 = 1
        wL1 = 0.5d0
        wL2 = 0.5d0
     else if (dtL > dtLgrid(1,nL)) then
        L1 = nL
        L2 = nL
        wL1 = 0.5d0
        wL2 = 0.5d0
     else
        do jj = 1, nL-1
           if ((dtL <= dtLgrid(1,jj+1)) .and. (dtL > dtLgrid(1,jj))) then
              L1 = jj
              L2 = jj + 1
              wL1 = (dtLgrid(1,L2) - dtL) / (dtLgrid(1,L2) - dtLgrid(1,L1))
              wL2 = 1.d0 - wL1
           end if
        end do
     end if

     if (L1 == 0 .or. L2 == 0) then
        call endrun (msg=' ERROR: LookupPsihat error: indices L1 and L2 not found')
     end if

     ! Find indices and weights for zdt values which bracket the specified zdt

     Z1 = 0 ; Z2 = 0
     if (zdt > zdtgrid(1,1)) then
        Z1 = 1
        Z2 = 1
        wZ1 = 0.5d0
        wZ2 = 0.5d0
     else if (zdt < zdtgrid(nZ,1)) then
        Z1 = nZ
        Z2 = nZ
        wZ1 = 0.5d0
        wZ2 = 0.5d0
     else
        do ii = 1, nZ-1
           if ((zdt >= zdtgrid(ii+1,1)) .and. (zdt < zdtgrid(ii,1))) then
              Z1 = ii
              Z2 = ii + 1
              wZ1 = (zdt - zdtgrid(ii+1,1)) / (zdtgrid(ii,1) - zdtgrid(ii+1,1))
              wZ2 = 1.d0 - wZ1
           end if
        end do
     end if

     if (Z1 == 0 .or. Z2 == 0) then
        call endrun (msg=' ERROR: LookupPsihat error: indices Z1 and Z2 not found')
     end if

     ! Calculate psihat as a weighted average of the values of psihat on the grid

     psihat = wZ1 * wL1 * psigrid(Z1,L1) + wZ2 * wL1 * psigrid(Z2,L1) &
          + wZ1 * wL2 * psigrid(Z1,L2) + wZ2 * wL2 * psigrid(Z2,L2)

   end subroutine LookupPsihat

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
