module PhotosynthesisAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C3
  use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C4
  use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_MEDLYN
  use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_BBERRY

  implicit none

  private

  type, public :: photosynthesis_auxvar_type

     ! Primary indepedent variables
     PetscReal :: ci        ! leaf intercellular CO2 (umol/mol)

     PetscReal :: tleaf     ! leaf temperature (K)
     PetscReal :: gbv       ! leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
     PetscReal :: gbc       ! leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
     PetscReal :: eair      ! vapor pressure profile (Pa)

     PetscReal :: cair      ! atmospheric CO2 profile (umol/mol)
     PetscReal :: o2ref     ! atmospheric O2 at reference height (mmol/mol)
     PetscReal :: apar      ! leaf absorbed PAR (umol photon/m2 leaf/s)

     PetscInt  :: colim     ! Photosynthesis co-limitation: 0 = no. 1 = yes
     PetscInt  :: c3psn     ! photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
     PetscInt  :: gstype    ! Stomatal conductance: 0 = Medlyn. 1 = Ball-Berry. 2 = WUE optimization

     PetscReal :: ceair     ! vapor pressure of air, constrained
     PetscReal :: esat      ! saturation vapor pressure (Pa)
     PetscReal :: desat     ! derivative of saturation vapor pressure w.r.t. leaf temperature (Pa/K)

     PetscReal :: g0opt     ! Ball-Berry minimum leaf conductance, unstressed (mol H2O/m2/s)
     PetscReal :: g1opt     ! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
     PetscReal :: g0        ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
     PetscReal :: g1        ! Ball-Berry slope of conductance-photosynthesis relationship
     PetscReal :: dpai      ! layer plant area index (m2/m2)
     PetscReal :: btran     ! Ball-Berry soil wetness factor (-)

     PetscReal :: vcmax25   ! leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
     PetscReal :: jmax25    ! C3 - maximum electron transport rate at 25C for canopy layer (umol/m2/s)
     PetscReal :: rd25      ! leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)
     PetscReal :: kp25      ! C4 - initial slope of CO2 response curve at 25C for canopy layer (mol/m2/s)

     PetscReal :: kc25      ! Michaelis-Menten constant for CO2 at 25C (umol/mol)
     PetscReal :: ko25      ! Michaelis-Menten constant for O2 at 25C (mmol/mol)
     PetscReal :: cp25      ! CO2 compensation point at 25C (umol/mol)

     PetscReal :: rdha      ! activation energy for rd (J/mol)
     PetscReal :: kcha      ! activation energy for kc (J/mol)
     PetscReal :: koha      ! activation energy for ko (J/mol)
     PetscReal :: cpha      ! activation energy for cp (J/mol)
     PetscReal :: vcmaxha   ! activation energy for vcmax (J/mol)
     PetscReal :: jmaxha    ! activation energy for jmax (J/mol)

     PetscReal :: rdhd      ! Deactivation energy for rd (J/mol)
     PetscReal :: vcmaxhd   ! Deactivation energy for vcmax (J/mol)
     PetscReal :: jmaxhd    ! Deactivation energy for jmax (J/mol)

     PetscReal :: rdse      ! entropy term for rd (J/mol/K)
     PetscReal :: vcmaxse   ! entropy term for vcmax (J/mol/K)
     PetscReal :: jmaxse    ! entropy term for jmax (J/mol/K)

     PetscReal :: vcmaxc    ! scaling factor for high temperature inhibition (25 C = 1.0)
     PetscReal :: jmaxc     ! scaling factor for high temperature inhibition (25 C = 1.0)
     PetscReal :: rdc       ! scaling factor for high temperature inhibition (25 C = 1.0)

     PetscReal :: colim_c3  ! Empirical curvature parameter for C3 co-limitation
     PetscReal :: colim_c4a ! Empirical curvature parameter for C4 co-limitation
     PetscReal :: colim_c4b ! Empirical curvature parameter for C4 co-limitation

     PetscReal :: qe_c4     ! Quantum yield, used only for C4 (mol CO2 / mol photons)
     PetscReal :: phi_psii  ! Quantum yield of PS II
     PetscReal :: theta_j   ! Empirical curvature parameter for electron transport rate

     PetscReal :: vcmax     ! maximum carboxylation rate (umol/m2/s)
     PetscReal :: jmax      ! maximum electron transport rate (umol/m2/s)
     PetscReal :: je        ! Electron transport rate (umol/m2/s)

     PetscReal :: kc        ! Michaelis-Menten constant for CO2 (umol/mol)
     PetscReal :: ko        ! Michaelis-Menten constant for O2 (mmol/mol)     
     PetscReal :: cp        ! CO2 compensation point (umol/mol)

     PetscReal :: rd        ! leaf respiration rate (umol CO2/m2 leaf/s)
     PetscReal :: kp        ! C4: Initial slope of CO2 response curve at 25C (mol/m2/s)

     PetscReal :: hs        ! leaf fractional humidity at leaf surface (dimensionless)
     PetscReal :: vpd       ! leaf vapor pressure deficit (Pa)
     PetscReal :: ac        ! leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal :: aj        ! leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal :: ap        ! leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal :: ag        ! leaf gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal :: an        ! leaf net photosynthesis (umol CO2/m2 leaf/s)
     PetscReal :: cs        ! leaf surface CO2 (umol/mol)
     PetscReal :: gs        ! leaf stomatal conductance (mol H2O/m2 leaf/s)

     PetscInt  :: pathway_and_stomatal_params_defined
   contains

     procedure, public :: Init          => PhotosynthesisInit
     procedure, public :: AuxVarCompute => PhotosynthesisAuxVarCompute

  end type photosynthesis_auxvar_type

  private :: ft    ! photosynthesis temperature response
  private :: fth   ! photosynthesis temperature inhibition
  private :: fth25 ! scaling factor for photosynthesis temperature inhibition

contains


  !------------------------------------------------------------------------
  function ft (tl, ha) result(ans)
    ! !DESCRIPTION:
    ! Photosynthesis temperature response
    !
    use MultiPhysicsProbConstants , only : RGAS, TFRZ
    !
    implicit none
    !
    ! !ARGUMENTS:
    PetscReal :: tl                   ! Leaf temperature (K)
    PetscReal :: ha                   ! Activation energy (J/mol)
    !
    ! !LOCAL VARIABLES:
    PetscReal :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = exp( ha / (RGAS * (TFRZ + 25.d0)) * (1.d0 - (TFRZ + 25.d0) / tl) )

  end function ft

  !-----------------------------------------------------------------------
  function fth (tl, hd, se, c) result(ans)
    !
    ! !DESCRIPTION:
    ! Photosynthesis temperature inhibition
    !
    use MultiPhysicsProbConstants , only : RGAS
    !
    implicit none
    !
    ! !ARGUMENTS:
    PetscReal :: tl                   ! Leaf temperature (K)
    PetscReal :: hd                   ! Deactivation energy (J/mol)
    PetscReal :: se                   ! Entropy term (J/mol/K)
    PetscReal :: c                    ! Scaling factor for high temperature inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    PetscReal :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = c / ( 1.d0 + exp( (-hd + se*tl) / (RGAS*tl) ) )

  end function fth

  !-----------------------------------------------------------------------
  function fth25 (hd, se) result(ans)
    !
    ! !DESCRIPTION:
    ! Scaling factor for photosynthesis temperature inhibition
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : RGAS, TFRZ
    !
    implicit none
    !
    ! !ARGUMENTS:
    PetscReal, intent(in) :: hd       ! Deactivation energy (J/mol)
    PetscReal, intent(in) :: se       ! Entropy term (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    PetscReal :: ans                  ! Temperature function value
    !---------------------------------------------------------------------

    ans = 1.d0 + exp( (-hd + se * (TFRZ + 25.d0)) / (RGAS * (TFRZ + 25.d0)) )

  end function fth25

  !------------------------------------------------------------------------
  subroutine PhotosynthesisInit(this)
    !
    implicit none
    !
    class(photosynthesis_auxvar_type) :: this
    PetscInt                          :: c3psn
    PetscInt                          :: gstype

    this%ci      = 0.d0

    this%tleaf   = 0.d0
    this%gbv     = 0.d0
    this%gbc     = 0.d0

    this%o2ref   = 0.d0
    this%ceair   = 0.d0
    this%esat    = 0.d0

    this%pathway_and_stomatal_params_defined = 0
    this%colim  =  1
    this%c3psn  = -1
    this%gstype = -1

    ! --- Kc, Ko, Cp at 25C
    this%kc25 = 404.9d0;
    this%ko25 = 278.4d0;
    this%cp25 = 42.75d0;

    ! --- Activation energy
    this%kcha    = 79430.d0;
    this%koha    = 36380.d0;
    this%cpha    = 37830.d0;
    this%rdha    = 46390.d0;
    this%vcmaxha = 65330.d0;
    this%jmaxha  = 43540.d0;

    ! --- High temperature deactivation
    ! Deactivation energy (J/mol)
    this%rdhd    = 150000.d0;
    this%vcmaxhd = 150000.d0;
    this%jmaxhd  = 150000.d0;

    ! Entropy term (J/mol/K)
    this%rdse    = 490.d0;
    this%vcmaxse = 490.d0;
    this%jmaxse  = 490.d0;

    ! Scaling factors for high temperature inhibition (25 C = 1.0).
    ! The factor "c" scales the deactivation to a value of 1.0 at 25C.
    this%vcmaxc = fth25 (this%vcmaxhd , this%vcmaxse );
    this%jmaxc  = fth25 (this%jmaxhd  , this%jmaxse  );
    this%rdc    = fth25 (this%rdhd    , this%rdse    );

    this%phi_psii = 0.85d0;
    this%theta_j  = 0.90d0;
    this%colim_c3 = 0.98d0;

    this%colim_c4a = 0.80d0;
    this%colim_c4b = 0.95d0;
    this%qe_c4     = 0.05d0;

    this%dpai    = 0.d0
    this%btran   = 1.d0
    this%eair    = 0.d0
    this%cair    = 0.d0
    this%apar    = 0.d0

    this%rd      = 0.d0
    this%hs      = 0.d0
    this%vpd     = 0.d0
    this%ac      = 0.d0
    this%aj      = 0.d0
    this%ap      = 0.d0
    this%ag      = 0.d0
    this%an      = 0.d0
    this%cs      = 0.d0
    this%gs      = 0.d0

  end subroutine PhotosynthesisInit

  !------------------------------------------------------------------------
  subroutine SetPathwayParameters(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : TFRZ
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    select case (this%c3psn)
    case (VAR_PHOTOSYNTHETIC_PATHWAY_C4) ! C4

       this%vcmax25  = 40.d0;
       this%jmax25   = 0.d0;
       this%kp25     = 0.02d0 * this%vcmax25;
       this%rd25     = 0.025d0 * this%vcmax25;

    case (VAR_PHOTOSYNTHETIC_PATHWAY_C3) ! C3

       this%vcmax25 = 57.7d0;
       this%jmax25  = 1.67d0 * this%vcmax25;
       this%kp25    = 0.d0;
       this%rd25    = 0.015d0 * this%vcmax25;

    case default
       write(iulog,*)'Unsupported photosynthesis pathway: ',this%c3psn
       call exit(0)
    end select

  end subroutine SetPathwayParameters

  !------------------------------------------------------------------------
  subroutine SetStomatalConductanceParameters(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : TFRZ
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this

    select case (this%c3psn)
    case (VAR_PHOTOSYNTHETIC_PATHWAY_C4) ! C4

       if (this%gstype == VAR_STOMATAL_CONDUCTANCE_BBERRY) then

          this%g0opt = 0.04d0;       ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
          this%g1opt = 4.0d0;        ! Ball-Berry slope of conductance-photosynthesis relationship

       elseif (this%gstype == VAR_STOMATAL_CONDUCTANCE_MEDLYN) then

          this%g0opt = 0.0d0;        ! Medlyn minimum leaf conductance (mol H2O/m2/s)
          this%g1opt = 1.62d0;       ! Medlyn slope of conductance-photosynthesis relationship

       else
          write(iulog,*)'Unsupported stomatal conductance: ',this%gstype
          call exit(0)
       end if

    case (VAR_PHOTOSYNTHETIC_PATHWAY_C3) ! C3

       if (this%gstype == VAR_STOMATAL_CONDUCTANCE_BBERRY) then

          this%g0opt = 0.01d0;       ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
          this%g1opt = 9.0d0;        ! Ball-Berry slope of conductance-photosynthesis relationship

       elseif (this%gstype == VAR_STOMATAL_CONDUCTANCE_MEDLYN) then

          this%g0opt = 0.0d0;        ! Medlyn minimum leaf conductance (mol H2O/m2/s)
          this%g1opt = 4.45d0;       ! Medlyn slope of conductance-photosynthesis relationship

       else
          write(iulog,*)'Unsupported stomatal conductance: ',this%gstype
          call exit(0)
       end if

    end select
  end subroutine SetStomatalConductanceParameters

  !------------------------------------------------------------------------
  subroutine PhotosynthesisAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    ! Computes various secondary quantities (sat, den, etc) based on
    ! the primary quantity (pressure).
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : TFRZ
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal :: t1, t2, t3

    if (this%pathway_and_stomatal_params_defined  == 0) then
       call SetPathwayParameters(this)
       call SetStomatalConductanceParameters(this)
       this%pathway_and_stomatal_params_defined = 1
    end if

    if (this%dpai > 0.d0) then

       select case(this%c3psn)
       case (VAR_PHOTOSYNTHETIC_PATHWAY_C4)

          call C4_Temperature_Response(this)
          call C4_Metabolic_Photosynthesis_Rate(this)
          call C4_Net_Assimilation(this)

       case (VAR_PHOTOSYNTHETIC_PATHWAY_C3)

          call C3_Temperature_Response(this)
          call Compute_Electron_Transport_Rate(this)
          call C3_Metabolic_Photosynthesis_Rate(this)
          call C3_Net_Assimilation(this)

       end select

       ! CO2 at leaf surface
       this%cs = max(this%cair - this%an/this%gbc, 1.d0)

       ! Saturation vapor pressure at leaf temperature
       call SatVap (this%tleaf, this%esat, this%desat)

       ! Constrain eair >= 0.05*esat[tleaf] so that solution does not blow up. This 
       ! ensures that hs does not go to zero. Also eair <= esat[tleaf] so that hs <= 1. 
       this%ceair = min( max(this%eair, 0.20d0*this%esat), this%esat )

       select case(this%gstype)
       case (VAR_STOMATAL_CONDUCTANCE_BBERRY)

          call GsBallBerry(this)

       case (VAR_STOMATAL_CONDUCTANCE_MEDLYN)

          call GsMedlyn(this)

       end select

    else
       write(iulog,*)'PhotosynthesisAuxVarCompute: Add code when dapi = 0.d0'
       call exit(0)
    end if

  end subroutine PhotosynthesisAuxVarCompute

  !------------------------------------------------------------------------

  subroutine C4_Temperature_Response(this)
    ! !USES:
    use MultiPhysicsProbConstants , only : TFRZ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal :: t1, t2, t3

    t1         = 2.0**( (this%tleaf - (TFRZ + 25.d0)) / 10.d0 )
    t2         = 1.d0 + exp( 0.2d0*((TFRZ + 15.d0) - this%tleaf))
    t3         = 1.d0 + exp( 0.3d0*(this%tleaf - (TFRZ + 40.d0)))
    this%vcmax = this%vcmax25 * t1 / (t2 * t3)

    t3         = 1.d0 + exp( 1.3d0*(this%tleaf - (TFRZ + 55.d0)))
    this%rd    = this%rd25 * t1 / t3

    this%kp    = this%kp25 * t1

    this%vcmax = this%vcmax * this%btran

  end subroutine C4_Temperature_Response

  !------------------------------------------------------------------------

  subroutine C3_Temperature_Response(this)
    ! !USES:
    use MultiPhysicsProbConstants , only : TFRZ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this

    this%kc    = this%kc25    * ft(this%tleaf, this%kcha    )
    this%ko    = this%ko25    * ft(this%tleaf, this%koha    )
    this%cp    = this%cp25    * ft(this%tleaf, this%cpha    )

    this%vcmax = this%vcmax25 * ft(this%tleaf, this%vcmaxha ) * fth(this%tleaf, this%vcmaxhd, this%vcmaxse, this%vcmaxc) 
    this%jmax  = this%jmax25  * ft(this%tleaf, this%jmaxha  ) * fth(this%tleaf, this%jmaxhd , this%jmaxse , this%jmaxc )
    this%rd    = this%rd25    * ft(this%tleaf, this%rdha    ) * fth(this%tleaf, this%rdhd   , this%rdse   , this%rdc   )

    this%kp = 0.d0

    this%vcmax = this%vcmax * this%btran

  end subroutine C3_Temperature_Response

  !------------------------------------------------------------------------

  subroutine Compute_Electron_Transport_Rate(this)
    ! !USES:
    !
    use MathUtilsMod, only : quadratic
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal :: qabs
    PetscReal :: aquad, bquad, cquad
    PetscReal :: root1, root2

    qabs = 0.5d0 * this%phi_psii * this%apar

    aquad = this%theta_j
    bquad = -(qabs + this%jmax)
    cquad = qabs * this%jmax

    call quadratic(aquad, bquad, cquad, root1, root2)

    this%je = min(root1, root2)

  end subroutine Compute_Electron_Transport_Rate

  !------------------------------------------------------------------------

  subroutine C4_Metabolic_Photosynthesis_Rate(this)
    ! !USES:
    use MultiPhysicsProbConstants , only : TFRZ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this

    ! Rubisco-limited photosynthesis
    this%ac = this%vcmax

    ! RuBP-limited photosynthesis
    this%aj = this%qe_c4 * this%apar

    ! PEP carboxylase-limited (CO2-limited)
    this%ap = this%kp * max(this%ci, 0.d0)

  end subroutine C4_Metabolic_Photosynthesis_Rate

  !------------------------------------------------------------------------

  subroutine C3_Metabolic_Photosynthesis_Rate(this)
    ! !USES:
    use MultiPhysicsProbConstants , only : TFRZ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this

    ! Rubisco-limited photosynthesis
    this%ac = this%vcmax * max(this%ci - this%cp, 0.d0)/(this%ci + this%kc*(1.d0 + this%o2ref/this%ko))

    ! RuBP-limited photosynthesis
    this%aj = this%je    * max(this%ci - this%cp, 0.d0)/(4.d0*this%ci + 8.d0*this%cp)

    ! Product-limited photosynthesis
    this%ap = 0.d0

  end subroutine C3_Metabolic_Photosynthesis_Rate

  !------------------------------------------------------------------------

  subroutine C4_Net_Assimilation(this)
    ! !USES:
    !
   use MathUtilsMod, only : quadratic
   !
   implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal :: aquad, bquad, cquad
    PetscReal :: root1, root2
    PetscReal :: ai

    select case(this%colim)
    case (1)
       aquad = this%colim_c4a
       bquad = -(this%ac + this%aj)
       cquad =   this%ac * this%aj

       call quadratic(aquad, bquad, cquad, root1, root2)
       ai = min(root1, root2)

       aquad = this%colim_c4b
       bquad = -(ai + this%ap)
       cquad =   ai * this%ap

       call quadratic(aquad, bquad, cquad, root1, root2)
       this%ag = min(root1, root2)

    case (0)

       this%ag = min(this%ac, this%aj, this%ap)

    end select
    
    this%ac = max(this%ac, 0.d0)
    this%aj = max(this%aj, 0.d0)
    this%ap = max(this%ap, 0.d0)
    this%ag = max(this%ag, 0.d0)

    this%an = this%ag - this%rd

  end subroutine C4_Net_Assimilation

  !------------------------------------------------------------------------

  subroutine C3_Net_Assimilation(this)
    !
    use MathUtilsMod, only : quadratic
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal :: aquad, bquad, cquad
    PetscReal :: root1, root2
    
    select case(this%colim)
    case (1)
       aquad = this%colim_c3
       bquad = -(this%ac + this%aj)
       cquad =   this%ac * this%aj

       call quadratic(aquad, bquad, cquad, root1, root2)

       this%ag = min(root1, root2)

    case (0)

       this%ag = min(this%ac, this%aj)

    end select

    this%ac = max(this%ac, 0.d0)
    this%aj = max(this%aj, 0.d0)
    this%ap = max(this%ap, 0.d0)
    this%ag = max(this%ag, 0.d0)

    this%an = this%ag - this%rd

  end subroutine C3_Net_Assimilation

  !------------------------------------------------------------------------
  subroutine GsBallBerry(this)
    ! !USES:
    !
    use MathUtilsMod, only : quadratic
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal :: aquad, bquad, cquad
    PetscReal :: root1, root2

    this%g0 = max( this%g0opt * this%btran, 1.d-06 )
    this%g1 = this%g1opt

    if (this%an > 0.d0) then

       aquad = this%cs
       bquad = this%cs   * (this%gbv - this%g0) - this%g1*this%an
       cquad = -this%gbv * (this%cs * this%g0 + this%g1 * this%an * this%ceair/this%esat)

       call quadratic (aquad, bquad, cquad, root1, root2)
       this%gs = max(root1, root2)
       
    else
       this%gs = this%g0
    end if


  end subroutine GsBallBerry

  !------------------------------------------------------------------------
  subroutine GsMedlyn(this)
    ! !USES:
    !
    use MathUtilsMod, only : quadratic
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal, parameter :: vpd_min = 100.d0
    PetscReal :: aquad, bquad, cquad
    PetscReal :: root1, root2
    PetscReal :: vpd_term, term

    this%g0 = this%g0opt
    this%g1 = this%g1opt

    if (this%an > 0.d0) then

       vpd_term = max((this%esat - this%ceair), vpd_min) * 0.001d0
       term = 1.6d0 * this%an / this%cs

       aquad = 1.d0
       bquad = -(2.d0 * (this%g0 + term) + (this%g1 * term)**2.d0 / (this%gbv * vpd_term))
       cquad = this%g0 * this%g0 + (2.d0 * this%g0 + term * (1.d0 - this%g1 * this%g1 / vpd_term)) * term

       call quadratic (aquad, bquad, cquad, root1, root2)
       this%gs = max(root1,root2)

    else

       this%gs = this%g0

    end if

  end subroutine GsMedlyn

#endif

end module PhotosynthesisAuxType

