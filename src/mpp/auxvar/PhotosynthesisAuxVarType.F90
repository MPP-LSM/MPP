module PhotosynthesisAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C3
  use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C4
  use MultiPhysicsProbConstants , only : VAR_SCM_MEDLYN
  use MultiPhysicsProbConstants , only : VAR_SCM_BBERRY
  use MultiPhysicsProbConstants , only : VAR_SCM_BONAN14
  use MultiPhysicsProbConstants , only : VAR_SCM_MODIFIED_BONAN14
  use MultiPhysicsProbConstants , only : VAR_SCM_MANZONI11
  use MultiPhysicsProbConstants , only : VAR_SCM_WUE
  use MultiPhysicsProbConstants , only : VAR_SCM_OSMWANG
  use petscsys

  implicit none

  private

  type, public :: root_auxvar_type
     PetscReal :: biomass ! fine root biomass (g biomass/m^2)
     PetscReal :: radius  ! fine root radius (m)
     PetscReal :: density ! fine root density (g biomass/m^3 root)
     PetscReal :: resist  ! hydraulic resistance of root (MPa.s.g/mmol H2O)
  end type root_auxvar_type

  type, public :: soil_auxvar_type
     PetscInt           :: nlevsoi       ! number of soil layers

     PetscReal, pointer :: h2osoi_vol(:) ! volumetric water content (m^3/m^3)
     PetscReal, pointer :: watsat(:)     ! volumetric content at saturation (i.e. porosity) (m^3/m^3)
     PetscReal, pointer :: psi_sat(:)    ! saturated matric potential
     PetscReal, pointer :: psi(:)        ! matric potential (mm)
     PetscReal, pointer :: hksat(:)      ! hydraulic conductivity at saturation (mm H2O/s)
     PetscReal, pointer :: bsw(:)        ! Capp and Hornberger 'b' parameter
     PetscReal, pointer :: rootfr(:)     ! fraction of roots
     PetscReal, pointer :: dz(:)         ! layer thickness (m)
   contains

     procedure, public :: AllocateMemory => SoilAuxVarAllocateMemory

  end type soil_auxvar_type

  type, public :: plant_auxvar_type
     PetscInt  :: nleaf                   ! number of leaves: 1 (single leaf) and 2 (sunlit and shaded leaves)

     PetscReal, pointer :: leaf_psi(:)    ! water potential from previous time step (MPa)
     PetscReal, pointer :: leaf_height(:) ! height (m)
     PetscReal, pointer :: leaf_capc(:)   ! leaf capacitance (mmol H2O/ m^2 leaf/MPa)
     PetscReal, pointer :: leaf_lsc(:)    ! leaf-specific conducntance (mmol H2O/m^2/MPa)
     PetscReal, pointer :: leaf_minlwp(:) ! minimum leaf water potential (MPa)
     PetscReal, pointer :: leaf_lai(:)    ! leaf area index (m^2/m^2)
     PetscReal, pointer :: k_stem2leaf(:) ! hydraulic conductance (mmol H2O/m^2 leaf area/MPa)

     PetscReal, pointer :: resist_soil(:) ! hydraulic resistance (MPa.s.m2/mmol H2O)
     PetscReal, pointer :: psi_soil(:)    ! weighted water potential (MPa)
     PetscReal, pointer :: dpsi_leaf(:)   ! change in weighted water potential (MPa)

     PetscReal, pointer :: weibull_c(:)   ! Weibull parameter [-]
     PetscReal, pointer :: weibull_b(:)   ! Weibull parameter [MPa]
     PetscReal, pointer :: leaf_psi_crt(:)! critical leaf pressure [MPa]

     PetscReal          :: dtime          ! Time step for plant hydraulics
   contains

     procedure, public :: AllocateMemory => PlantAuxVarAllocateMemory

  end type plant_auxvar_type

  type, public :: photosynthesis_auxvar_type

     ! Primary indepedent variables
     PetscReal , pointer :: ci(:)      ! leaf intercellular CO2 (umol/mol)

     PetscInt            :: ndof       ! no. of unknowns. 1 for BB, Medlyn, and WUE; 2 for Bonan14

     PetscReal           :: tleaf      ! leaf temperature (K)
     PetscReal           :: tleaf_prev ! leaf temperature at previous time step (K)

     PetscReal           :: gbv        ! leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
     PetscReal           :: gbc        ! leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
     PetscReal           :: eair       ! vapor pressure profile (Pa)
     PetscReal           :: pref       ! atmospheric pressure at reference height (Pa)

     PetscReal           :: cair       ! atmospheric CO2 profile (umol/mol)
     PetscReal           :: o2ref      ! atmospheric O2 at reference height (mmol/mol)
     PetscReal           :: apar       ! leaf absorbed PAR (umol photon/m2 leaf/s)

     PetscInt            :: colim      ! Photosynthesis co-limitation: 0 = no. 1 = yes
     PetscInt            :: c3psn      ! photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
     PetscInt            :: gstype     ! Stomatal conductance: 0 = Medlyn. 1 = Ball-Berry. 2 = WUE optimization

     PetscReal           :: ceair      ! vapor pressure of air, constrained
     PetscReal           :: esat       ! saturation vapor pressure (Pa)
     PetscReal           :: desat      ! derivative of saturation vapor pressure w.r.t. leaf temperature (Pa/K)

     PetscReal           :: g0opt      ! Ball-Berry minimum leaf conductance, unstressed (mol H2O/m2/s)
     PetscReal           :: g1opt      ! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
     PetscReal           :: g0         ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
     PetscReal           :: g1         ! Ball-Berry slope of conductance-photosynthesis relationship
     PetscReal           :: dpai       ! layer plant area index (m2/m2)
     PetscReal           :: btran      ! Ball-Berry soil wetness factor (-)

     PetscReal           :: vcmax25    ! leaf maximum carboxylation rate at 25C for canopy layer (umol/m2/s)
     PetscReal           :: jmax25     ! C3 - maximum electron transport rate at 25C for canopy layer (umol/m2/s)
     PetscReal           :: rd25       ! leaf respiration rate at 25C for canopy layer (umol CO2/m2/s)
     PetscReal           :: kp25       ! C4 - initial slope of CO2 response curve at 25C for canopy layer (mol/m2/s)

     PetscReal           :: kc25       ! Michaelis-Menten constant for CO2 at 25C (umol/mol)
     PetscReal           :: ko25       ! Michaelis-Menten constant for O2 at 25C (mmol/mol)
     PetscReal           :: cp25       ! CO2 compensation point at 25C (umol/mol)

     PetscReal           :: rdha       ! activation energy for rd (J/mol)
     PetscReal           :: kcha       ! activation energy for kc (J/mol)
     PetscReal           :: koha       ! activation energy for ko (J/mol)
     PetscReal           :: cpha       ! activation energy for cp (J/mol)
     PetscReal           :: vcmaxha    ! activation energy for vcmax (J/mol)
     PetscReal           :: jmaxha     ! activation energy for jmax (J/mol)

     PetscReal           :: rdhd       ! Deactivation energy for rd (J/mol)
     PetscReal           :: vcmaxhd    ! Deactivation energy for vcmax (J/mol)
     PetscReal           :: jmaxhd     ! Deactivation energy for jmax (J/mol)

     PetscReal           :: rdse       ! entropy term for rd (J/mol/K)
     PetscReal           :: vcmaxse    ! entropy term for vcmax (J/mol/K)
     PetscReal           :: jmaxse     ! entropy term for jmax (J/mol/K)

     PetscReal           :: vcmaxc     ! scaling factor for high temperature inhibition (25 C = 1.0)
     PetscReal           :: jmaxc      ! scaling factor for high temperature inhibition (25 C = 1.0)
     PetscReal           :: rdc        ! scaling factor for high temperature inhibition (25 C = 1.0)

     PetscReal           :: colim_c3   ! Empirical curvature parameter for C3 co-limitation
     PetscReal           :: colim_c4a  ! Empirical curvature parameter for C4 co-limitation
     PetscReal           :: colim_c4b  ! Empirical curvature parameter for C4 co-limitation

     PetscReal           :: qe_c4      ! Quantum yield, used only for C4 (mol CO2 / mol photons)
     PetscReal           :: phi_psii   ! Quantum yield of PS II
     PetscReal           :: theta_j    ! Empirical curvature parameter for electron transport rate

     PetscReal           :: vcmax      ! maximum carboxylation rate (umol/m2/s)
     PetscReal           :: jmax       ! maximum electron transport rate (umol/m2/s)
     PetscReal           :: je         ! Electron transport rate (umol/m2/s)

     PetscReal           :: kc         ! Michaelis-Menten constant for CO2 (umol/mol)
     PetscReal           :: ko         ! Michaelis-Menten constant for O2 (mmol/mol)
     PetscReal           :: cp         ! CO2 compensation point (umol/mol)

     PetscReal           :: rd         ! leaf respiration rate (umol CO2/m2 leaf/s)
     PetscReal           :: kp         ! C4: Initial slope of CO2 response curve at 25C (mol/m2/s)

     PetscReal           :: hs         ! leaf fractional humidity at leaf surface (dimensionless)
     PetscReal           :: vpd        ! leaf vapor pressure deficit (Pa)

     PetscReal , pointer :: ac(:)      ! leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal , pointer :: aj(:)      ! leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal , pointer :: ap(:)      ! leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal , pointer :: ag(:)      ! leaf gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal , pointer :: an(:)      ! leaf net photosynthesis (umol CO2/m2 leaf/s)
     PetscReal           :: cs         ! leaf surface CO2 (umol/mol)
     PetscReal , pointer :: gs(:)      ! leaf stomatal conductance (mol H2O/m2 leaf/s)
     PetscReal , pointer :: gleaf_c(:) ! leaf-to-air conductance to CO2 (mol CO2/m2 leaf/s)
     PetscReal , pointer :: gleaf_w(:) ! leaf-to-air conductance to H2O (mol H2O/m2 leaf/s)

     PetscReal           :: ac_soln      ! leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal           :: aj_soln      ! leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal           :: ap_soln      ! leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal           :: ag_soln      ! leaf gross photosynthesis (umol CO2/m2 leaf/s)
     PetscReal           :: an_soln      ! leaf net photosynthesis (umol CO2/m2 leaf/s)
     PetscReal           :: gs_soln      ! leaf stomatal conductance (mol H2O/m2 leaf/s)
     PetscReal           :: gleaf_c_soln ! leaf-to-air conductance to CO2 (mol CO2/m2 leaf/s)
     PetscReal           :: gleaf_w_soln ! leaf-to-air conductance to H2O (mol H2O/m2 leaf/s)

     PetscReal           :: iota               ! stomatal efficiency (umol CO2/ mol H2O)
     PetscReal , pointer :: residual_wue(:)    ! residual for WUE equation (umol CO2/ mol H2O)
     PetscReal , pointer :: residual_hyd(:)    ! residual for plant hydraulics equation (Pa)
     PetscBool , pointer :: soln_is_bounded(:) ! Is 0 if solution is not bounded between gs_min_wue and gs_max_wue otherwise 1
     PetscReal           :: manzoni11_beta

     PetscReal           :: fdry       ! Fraction of plant area index that is green and dry
     PetscReal           :: fwet       ! Fraction of plant area index that is wet

     PetscReal , pointer :: dac_dci(:) ! deriavate of leaf rubisco-limited gross photosynthesis wrt leaf CO2 (mol/m2 leaf/s)
     PetscReal , pointer :: daj_dci(:) ! deriavate of leaf RuBP-limited gross photosynthesis wrt leaf CO2 (mol/m2 leaf/s)
     PetscReal , pointer :: dap_dci(:) ! deriavate of leaf product-limited (C3), CO2-limited (C4) gross photosynthesis wrt leaf CO2 (mol/m2 leaf/s)
     PetscReal , pointer :: dag_dci(:) ! deriavate of leaf gross photosynthesis wrt leaf CO2 (mol/m2 leaf/s)
     PetscReal , pointer :: dan_dci(:) ! deriavate of leaf net photosynthesis wrt leaf CO2 (mol/m2 leaf/s)

     PetscInt            :: pathway_and_stomatal_params_defined

     type(root_auxvar_type) , pointer :: root
     type(soil_auxvar_type) , pointer :: soil
     type(plant_auxvar_type), pointer :: plant

   contains

     procedure, public :: Init                         => PhotosynthesisInit
     procedure, public :: AuxVarCompute                => PhotosynthesisAuxVarCompute
     procedure, public :: DetermineIfSolutionIsBounded => PhotosynthesisDetermineIfSolutionIsBounded
     procedure, public :: PreSolve                     => PhotosynthesisPreSolve
     procedure, public :: PostSolve                    => PhotosynthesisPostSolve
     procedure, public :: SetDefaultParameters         => PhotosynthesisAuxVarSetDefaultParameters

  end type photosynthesis_auxvar_type

  private :: ft    ! photosynthesis temperature response
  private :: fth   ! photosynthesis temperature inhibition
  private :: fth25 ! scaling factor for photosynthesis temperature inhibition

  PetscReal, parameter :: gs_min = 1.d-6
  PetscReal, parameter :: gs_min_wue = 0.005d0
  PetscReal, parameter :: gs_max_wue = 2.0d0
  PetscReal, parameter :: gs_delta_wue = 1.d-5
  PetscReal, parameter :: gs_min_bonan14 = 0.002d0
  PetscReal, parameter :: gs_delta_bonan14 = 0.001d0

contains

  !------------------------------------------------------------------------
  subroutine SoilAuxVarAllocateMemory(this)
    !
    implicit none
    !
    class(soil_auxvar_type) :: this

    if (this%nlevsoi == 0) then
       call endrun(msg=' ERROR: Number of soil layers is zero '//&
            errMsg(__FILE__, __LINE__))
    end if

    allocate(this%h2osoi_vol( this%nlevsoi))
    allocate(this%watsat(     this%nlevsoi))
    allocate(this%psi_sat(    this%nlevsoi))
    allocate(this%psi(        this%nlevsoi))
    allocate(this%hksat(      this%nlevsoi))
    allocate(this%bsw(        this%nlevsoi))
    allocate(this%rootfr(     this%nlevsoi))
    allocate(this%dz(         this%nlevsoi))

    this%h2osoi_vol (:) = 0.d0
    this%watsat     (:) = 0.d0
    this%psi        (:) = 0.d0
    this%hksat      (:) = 0.d0
    this%bsw        (:) = 0.d0
    this%rootfr     (:) = 0.d0
    this%dz         (:) = 0.d0

  end subroutine SoilAuxVarAllocateMemory

  !------------------------------------------------------------------------
  subroutine PlantAuxVarAllocateMemory(this)
    !
    implicit none
    !
    class(plant_auxvar_type) :: this

    if (this%nleaf == 0) then
       call endrun(msg=' ERROR: Number of leaves is zero '//&
            errMsg(__FILE__, __LINE__))
    end if

    allocate(this%leaf_psi    (this%nleaf))
    allocate(this%leaf_height (this%nleaf))
    allocate(this%leaf_capc   (this%nleaf))
    allocate(this%leaf_lsc    (this%nleaf))
    allocate(this%leaf_minlwp (this%nleaf))
    allocate(this%leaf_lai    (this%nleaf))
    allocate(this%k_stem2leaf (this%nleaf))
    allocate(this%resist_soil (this%nleaf))
    allocate(this%psi_soil    (this%nleaf))
    allocate(this%dpsi_leaf   (this%nleaf))

    this%leaf_psi    (:) = 0.d0
    this%leaf_height (:) = 0.d0
    this%leaf_capc   (:) = 0.d0
    this%leaf_lsc    (:) = 0.d0
    this%leaf_minlwp (:) = 0.d0
    this%leaf_lai    (:) = 0.d0
    this%k_stem2leaf (:) = 0.d0
    this%resist_soil (:) = 0.d0
    this%psi_soil    (:) = 0.d0
    this%dpsi_leaf   (:) = 0.d0

  end subroutine PlantAuxVarAllocateMemory

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
  subroutine PhotosynthesisInit(this, ndof)
    !
    implicit none
    !
    class(photosynthesis_auxvar_type) :: this
    PetscInt                          :: ndof


    this%ndof = ndof

    allocate(this%ci(ndof))
    this%ci(:) = 0.d0

    this%tleaf     = 0.d0
    this%tleaf_prev= 0.d0
    this%gbv       = 0.d0
    this%gbc       = 0.d0

    this%o2ref   = 0.d0
    this%ceair   = 0.d0
    this%esat    = 0.d0
    this%pref    = 101325.d0

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

    allocate(this%ac(ndof))
    allocate(this%aj(ndof))
    allocate(this%ap(ndof))
    allocate(this%ag(ndof))
    allocate(this%an(ndof))

    this%ac(:)   = 0.d0
    this%aj(:)   = 0.d0
    this%ap(:)   = 0.d0
    this%ag(:)   = 0.d0
    this%an(:)   = 0.d0

    allocate(this%gs(ndof))
    allocate(this%gleaf_c(ndof))
    allocate(this%gleaf_w(ndof))
    this%gs(:)      = 0.d0
    this%gleaf_c(:) = 0.d0
    this%gleaf_w(:) = 0.d0

    this%cs      = 0.d0

    this%iota    = 750.d0
    this%manzoni11_beta = -1.5d0 ! Liu et al. (2020) Nature Climate Change
    allocate(this%residual_wue(ndof))
    allocate(this%residual_hyd(ndof))
    allocate(this%soln_is_bounded(ndof))
    this%soln_is_bounded(:) = PETSC_TRUE

    this%fdry = 0.d0
    this%fwet = 0.d0

    allocate(this%dac_dci(ndof))
    allocate(this%daj_dci(ndof))
    allocate(this%dap_dci(ndof))
    allocate(this%dag_dci(ndof))
    allocate(this%dan_dci(ndof))

    this%dac_dci(:) = 0.d0
    this%daj_dci(:) = 0.d0
    this%dap_dci(:) = 0.d0
    this%dag_dci(:) = 0.d0
    this%dan_dci(:) = 0.d0

    allocate(this%root)
    allocate(this%soil)
    allocate(this%plant)

    this%root%biomass = 0.d0
    this%root%radius  = 0.d0
    this%root%density = 0.d0
    this%root%resist  = 0.d0

    this%soil%nlevsoi = 0

    allocate(this%plant%weibull_c(ndof))
    allocate(this%plant%weibull_b(ndof))
    allocate(this%plant%leaf_psi_crt(ndof))

    this%plant%nleaf  = 0
    this%plant%dtime  = 300.d0

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
       call endrun(msg=errMsg(__FILE__, __LINE__))

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

       select case(this%gstype)
       case (VAR_SCM_BBERRY)

          this%g0opt = 0.04d0;       ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
          this%g1opt = 4.0d0;        ! Ball-Berry slope of conductance-photosynthesis relationship

       case (VAR_SCM_MEDLYN)

          this%g0opt = 1.0d-4;        ! Medlyn minimum leaf conductance (mol H2O/m2/s)
          this%g1opt = 1.62d0;       ! Medlyn slope of conductance-photosynthesis relationship

       case(VAR_SCM_WUE, VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14, VAR_SCM_MANZONI11, VAR_SCM_OSMWANG)

       case default
          write(iulog,*)'Unsupported stomatal conductance: ',this%gstype
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

    case (VAR_PHOTOSYNTHETIC_PATHWAY_C3) ! C3

       select case(this%gstype)
          case(VAR_SCM_BBERRY)

          this%g0opt = 0.01d0;       ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
          this%g1opt = 9.0d0;        ! Ball-Berry slope of conductance-photosynthesis relationship

       case (VAR_SCM_MEDLYN)

          this%g0opt = 1.0d-4;        ! Medlyn minimum leaf conductance (mol H2O/m2/s)
          this%g1opt = 4.45d0;       ! Medlyn slope of conductance-photosynthesis relationship

       case (VAR_SCM_WUE, VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14, VAR_SCM_MANZONI11, VAR_SCM_OSMWANG)

       case default
          write(iulog,*)'Unsupported stomatal conductance: ',this%gstype
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

    end select
  end subroutine SetStomatalConductanceParameters

  !------------------------------------------------------------------------
  subroutine SetPlantParameters(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    implicit none
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type)    :: this

    select case(this%gstype)
    case (VAR_SCM_BBERRY, VAR_SCM_MEDLYN, VAR_SCM_WUE, VAR_SCM_BONAN14, VAR_SCM_MANZONI11)

    case(VAR_SCM_MODIFIED_BONAN14)
       this%plant%weibull_b(:) = 2.d0
       this%plant%weibull_c(:) = 5.d0

    case (VAR_SCM_OSMWANG)
       this%plant%weibull_b(:) = 1.2d0
       this%plant%weibull_c(:) = 5.d0

       this%plant%leaf_psi_crt(:) = - 1.2d0 * (log(1.d3))**(1.d0/5.d0)

    case default
       write(iulog,*)'Unsupported stomatal conductance: ',this%gstype
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SetPlantParameters
  !------------------------------------------------------------------------
  subroutine ComputeSoilResistance(this)
    !
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT, MM_H2O
    !
    implicit none
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type)    :: this
    !
    class(soil_auxvar_type)  , pointer   :: soil
    class(plant_auxvar_type) , pointer   :: plant
    class(root_auxvar_type)  , pointer   :: root
    PetscInt                             :: ileaf, j
    PetscReal                            :: s, hk, head, totevap
    PetscReal                            :: root_biomass_density, root_length_density, root_dist, root_cross_sec_area
    PetscReal                            :: soilr, soilr1, soilr2, blw_grnd_conductance, vwc
    PetscReal                , pointer   :: psi_mpa(:), evap(:)
    PetscReal                , parameter :: g = 9.80665d0
    PetscReal                , parameter :: denh2o = 1000.d0
    !PetscReal                , parameter :: mmh2o = 18.02d-3 ! molecular mass of water (kg/mol)

    soil  => this%soil
    plant => this%plant
    root  => this%root

    allocate(psi_mpa (soil%nlevsoi))
    allocate(evap    (soil%nlevsoi))

    head                = g * denh2o * 1.d-6 ! MPa/m
    root_cross_sec_area = PETSC_PI * (root%radius**2.d0);

    do ileaf = 1, plant%nleaf

       blw_grnd_conductance = 0.d0

       do j = 1, soil%nlevsoi

          vwc = max(soil%h2osoi_vol(j),1.0d-6)/(soil%dz(j)*denh2o)

          s = max(min(vwc/soil%watsat(j), 1.d0), 0.01d0);

          hk = soil%hksat(j) * s**(2.d0 * soil%bsw(j) + 3.d0); ! mm/s
          hk = hk * 1.d-03 / head;                              ! mm/s -> m/s -> m2/s/MPa
          hk = hk * denh2o / MM_H2O * 1000.d0;                  ! m2/s/MPa -> mmol/m/s/MPa

          ! Matric potential for each layer (MPa)
          soil%psi(j) = soil%psi_sat(j)* s**(-soil%bsw(j))
          psi_mpa(j) = soil%psi(j) * 1.d-03 * head;          ! mm -> m -> MPa

          ! Root biomass density (g biomass / m3 soil)
          root_biomass_density = root%biomass * soil%rootfr(j) / soil%dz(j);
          root_biomass_density = max(root_biomass_density, 1.d-10);

          ! Root length density (m root per m3 soil)
          root_length_density = root_biomass_density / (root%density * root_cross_sec_area);

          ! Distance between roots (m)
          root_dist = sqrt(1.d0 / (root_length_density * PETSC_PI));

          ! Soil-to-root resistance (MPa.s.m2/mmol H2O)
          soilr1 = log(root_dist/root%radius) / (2.d0 * PETSC_PI * root_length_density * soil%dz(j) * hk);

          ! Root-to-stem resistance (MPa.s.m2/mmol H2O)
          soilr2 = root%resist / (root_biomass_density * soil%dz(j));

          ! Belowground resistance (MPa.s.m2/mmol H2O)
          soilr = soilr1 + soilr2;

          ! Total belowground resistance. First sum the conductances (1/soilr)
          ! for each soil layer and then convert back to a resistance after the
          ! summation

          blw_grnd_conductance = blw_grnd_conductance + 1.d0 / soilr;

          ! Maximum transpiration for each layer (mmol H2O/m2/s)

          evap(j) = (psi_mpa(j) - plant%leaf_minlwp(ileaf)) / soilr;
          evap(j) = max (evap(j), 0.d0);
       end do

       ! Belowground resistance: resistance = 1 / conductance
       plant%resist_soil(ileaf) = plant%leaf_lai(ileaf)/blw_grnd_conductance

       ! Weighted soil water potential (MPa)
       plant%psi_soil(ileaf) = 0.d0;
       totevap               = 0.d0

       do j = 1, soil%nlevsoi
          plant%psi_soil(ileaf) = plant%psi_soil(ileaf) + psi_mpa(j) * evap(j);
          totevap = totevap + evap(j)
       end do

       if (totevap > 0.d0) then
          plant%psi_soil(ileaf) = plant%psi_soil(ileaf)/totevap
       else
          plant%psi_soil(ileaf) = plant%leaf_minlwp(ileaf)
       end if

       plant%leaf_lsc(ileaf) = 1.d0/(1.d0/ plant%k_stem2leaf(ileaf) + plant%resist_soil(ileaf))

    enddo

    deallocate(psi_mpa)
    deallocate(evap   )

  end subroutine ComputeSoilResistance

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
    class(photosynthesis_auxvar_type)    :: this
    !
    class(plant_auxvar_type) , pointer   :: plant
    PetscReal                            :: t1, t2, t3, etflx
    PetscReal                            :: an, gbc, delta_c
    PetscReal                            :: gs_val_wue, gs_val_hyd, an_low, an_high, check, psi_new, cost_low, cost_high
    PetscReal                            :: etflx_low, etflx_high, etflx_c
    PetscReal                            :: factor
    PetscInt                 , parameter :: ileaf = 1 ! Currently only one leaf is supported
    PetscReal                , parameter :: g = 9.80665d0
    PetscReal                , parameter :: denh2o = 1000.d0
    PetscInt , parameter :: idof_wue = 1
    PetscInt , parameter :: idof_hyd = 2
    PetscInt :: ii

    plant => this%plant
    select case (this%gstype)
    case (VAR_SCM_BBERRY)
       call PhotosynthesisAuxVarCompute_SemiEmpirical(this)
       call ComputeSoilResistance(this)

    case (VAR_SCM_MEDLYN)
       call PhotosynthesisAuxVarCompute_SemiEmpirical(this)
       call ComputeSoilResistance(this)

    case (VAR_SCM_WUE, VAR_SCM_MANZONI11)
       gs_val_wue = this%gs(idof_wue)

       this%gs(idof_wue) = gs_val_wue - gs_delta_wue
       call PhotosynthesisAuxVarCompute_WUE(this)
       an_low = this%an(idof_wue)

       this%gs(idof_wue) = gs_val_wue
       call PhotosynthesisAuxVarCompute_WUE(this)
       an_high = this%an(idof_wue)

       call ComputeSoilResistance(this)

       factor = 1.d0
       if (this%gstype == VAR_SCM_MANZONI11) then
          factor = exp(this%manzoni11_beta * plant%leaf_psi(ileaf))
       end if
       this%residual_wue(idof_wue) = (an_high - an_low) - this%iota * factor * gs_delta_wue * this%vpd

    case (VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14)
       gs_val_wue = this%gs(idof_wue)
       gs_val_hyd = this%gs(idof_hyd)

       this%gs(idof_wue) = gs_val_wue - gs_delta_bonan14
       this%gs(idof_hyd) = gs_val_hyd - gs_delta_bonan14
       call PhotosynthesisAuxVarCompute_WUE(this)
       an_low = this%an(idof_wue)

       this%gs(idof_wue) = gs_val_wue
       this%gs(idof_hyd) = gs_val_hyd
       call PhotosynthesisAuxVarCompute_WUE(this)
       an_high = this%an(idof_wue)

       this%residual_wue(idof_wue) = (an_high - an_low) - this%iota * gs_delta_bonan14 * this%vpd

       !
       ! Residual for hydraulics equation =  psi_{t} + dpsi_{t+1} - leaf_minlwp
       !
       call ComputeSoilResistance(this)

       etflx = (this%esat - this%eair)/this%pref * this%gleaf_w(idof_hyd) * this%fdry

       plant => this%plant
       if (this%gstype == VAR_SCM_BONAN14) then
          call ComputeChangeInPsi(plant, etflx)
          this%residual_hyd(idof_hyd) =  plant%leaf_psi(ileaf) + plant%dpsi_leaf(ileaf) - plant%leaf_minlwp(ileaf)
       else
          call ComputePsi_ModifiedBonan14(plant, etflx, psi_new)
          plant%dpsi_leaf(ileaf) = psi_new - plant%leaf_psi(ileaf)
          this%residual_hyd(idof_hyd) = psi_new - plant%leaf_minlwp(ileaf)
          !this%residual_hyd(idof_hyd) = 0.d0
       end if

    case (VAR_SCM_OSMWANG)

       gs_val_wue = this%gs(idof_wue)

       ! compute value at 'gs - delta_gs'
       this%gs(idof_wue) = gs_val_wue - gs_delta_wue
       call PhotosynthesisAuxVarCompute_WUE(this)

       an_low = this%an(idof_wue)
       etflx_low  = (this%esat - this%eair)/this%pref * this%gleaf_w(idof_wue) * this%fdry

       ! compute value at 'gs'
       this%gs(idof_wue) = gs_val_wue
       call PhotosynthesisAuxVarCompute_WUE(this)

       an_high = this%an(idof_wue)
       etflx_high  = (this%esat - this%eair)/this%pref * this%gleaf_w(idof_wue) * this%fdry

       call ComputeSoilResistance(this)

       call ComputeCriticalFlow(plant, etflx_c)

       cost_low  = an_low *etflx_low /etflx_c
       cost_high = an_high*etflx_high/etflx_c

       this%residual_wue(idof_wue) = (an_high - an_low)/(etflx_high - etflx_low) - an_high/(etflx_c - etflx_high)!, 1.d-7)

       this%residual_wue(idof_wue) = (an_high - an_low)/(etflx_high - etflx_low) - (cost_high - cost_low)/(etflx_high - etflx_low)

    case default
       write(iulog,*)'Unsupported stomatal conductance: ',this%gstype
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select
  end subroutine PhotosynthesisAuxVarCompute

  !------------------------------------------------------------------------
  subroutine ComputeCriticalFlow(plant, etflx_crt)
    !
    type(plant_auxvar_type) , pointer       :: plant
    PetscReal               , intent(inout) :: etflx_crt        ! [mol H2O/m2/s]
    !
    PetscInt                , parameter     :: ileaf = 1
    PetscReal               , parameter     :: g = 9.80665d0    ! [m/s^2]
    PetscReal               , parameter     :: denh2o = 1000.d0 ! [kg/m^3]
    PetscReal                               :: Kmax             ! [mmol/m^2/MPa]
    PetscReal                               :: kl               ! [mmol/m^2/MPa]
    PetscReal                               :: w_c              ! [-]
    PetscReal                               :: w_b              ! [MPa]
    PetscReal                               :: psi_crt          ! [MPa]
    PetscReal                               :: psi_soil         ! [MPa]
    PetscReal                               :: head             ! [MPa/m]
    PetscReal                               :: h                ! [m]

    head     = g * denh2o * 1.d-6

    Kmax     = plant%leaf_lsc     (ileaf)
    w_b      = plant%weibull_b    (ileaf)
    w_c      = plant%weibull_c    (ileaf)
    psi_crt  = plant%leaf_psi_crt (ileaf)
    psi_soil = plant%psi_soil     (ileaf)
    h        = plant%leaf_height  (ileaf)

    kl   = Kmax*exp(-(-psi_soil/w_b)**w_c)
    etflx_crt = max(Kl * (psi_soil - psi_crt - head * h)/1.d3, 1.d-7)

  end subroutine ComputeCriticalFlow

  !------------------------------------------------------------------------
  subroutine ComputeChangeInPsi(plant, etflx)
    !
    type(plant_auxvar_type) , pointer   :: plant
    PetscReal                            :: etflx, dtime
    !
    PetscReal                , parameter :: g = 9.80665d0
    PetscReal                , parameter :: denh2o = 1000.d0
    PetscReal                            :: head, a, b
    PetscInt                 , parameter :: ileaf = 1 ! Currently only one leaf is supported

    head = g * denh2o * 1.d-6 ! MPa/m

    a  = plant%psi_soil(ileaf) - head * plant%leaf_height(ileaf) - 1.d3 * etflx/plant%leaf_lsc(ileaf)
    b  = plant%leaf_capc(ileaf) / plant%leaf_lsc(ileaf)

    plant%dpsi_leaf(ileaf) = (a - plant%leaf_psi(ileaf)) * (1.d0 - exp(-plant%dtime/b));

  end subroutine ComputeChangeInPsi

  !------------------------------------------------------------------------
  subroutine ComputePsi_ModifiedBonan14(plant, etflx, psi_new)
    !
    use MathUtilsMod, only : quadratic
    !
    type(plant_auxvar_type) , pointer     :: plant
    PetscReal                             :: etflx
    PetscReal               , intent(out) :: psi_new
    !
    PetscInt                , parameter :: ileaf = 1 ! Currently only one leaf is supported
    PetscReal               , parameter :: g = 9.80665d0
    PetscReal               , parameter :: denh2o = 1000.d0
    PetscReal               , parameter :: w_b = 2 ! [MPa]
    PetscReal               , parameter :: w_c = 5 ! [-]
    PetscReal               , parameter :: pert = 1.d-2 ! [MPa]
    PetscReal                           :: dK_dpsi, Cp, dt, K, Kmax, psi_old, psi_soil, h, E, head
    PetscReal                           :: Kl, Kl_pert, dKl_dpsi
    PetscReal                           :: aquad, bquad, cquad
    PetscReal                           :: root1, root2
    PetscReal                           :: k_s2l, resist_s, factor

    head    = g * denh2o * 1.d-6 ! MPa/m
    Cp      = plant%leaf_capc(ileaf)
    dt      = plant%dtime
    psi_old = plant%leaf_psi(ileaf)
    psi_soil= plant%psi_soil(ileaf)
    h       = plant%leaf_height(ileaf)
    E       = 1.d3*etflx

    ! Compute leaf-specific hydraulic conductance while accounting
    ! for leaf water potential on stem-to-leaf conductance through
    ! a weibull function
    factor  = exp(-(-psi_old/w_b)**w_c)
    k_s2l   = plant%k_stem2leaf(ileaf)
    resist_s= plant%resist_soil(ileaf)
    Kl      = 1.d0/( 1.d0/(k_s2l * factor) + resist_s)

    psi_new = ( (Cp/dt) * psi_old + Kl*(psi_soil - head * h) - E )/ ( Cp/dt + Kl)

  end subroutine ComputePsi_ModifiedBonan14

  !------------------------------------------------------------------------
  subroutine PhotosynthesisAuxVarSetDefaultParameters(this)
    !
    ! !DESCRIPTION:
    ! Set default parameters depending on photosynthetic pathway (i.e., C3/C4)
    ! and stomatal conductance model.
    !
    ! !USES:
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !

    if (this%pathway_and_stomatal_params_defined  == 0) then
       call SetPathwayParameters(this)
       call SetStomatalConductanceParameters(this)
       call SetPlantParameters(this)
       this%pathway_and_stomatal_params_defined = 1
    end if

end subroutine PhotosynthesisAuxVarSetDefaultParameters


  !------------------------------------------------------------------------
  subroutine PhotosynthesisAuxVarCompute_SemiEmpirical(this)
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
    PetscReal                            :: t1, t2, t3, a, b, dy, head
    PetscReal                            :: an, gbc, delta_c
    PetscInt                             :: ileaf, idof
    class(plant_auxvar_type) , pointer   :: plant
    PetscReal                , parameter :: g = 9.80665d0
    PetscReal                , parameter :: denh2o = 1000.d0

    call this%SetDefaultParameters()

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
       do idof = 1, this%ndof
          this%cs = max(this%cair - this%an(idof)/this%gbc, 1.d0)
       end do

       ! Saturation vapor pressure at leaf temperature
       call SatVap (this%tleaf, this%esat, this%desat)

       ! Constrain eair >= 0.05*esat[tleaf] so that solution does not blow up. This 
       ! ensures that hs does not go to zero. Also eair <= esat[tleaf] so that hs <= 1. 
       this%ceair = min( max(this%eair, 0.20d0*this%esat), this%esat )

       ! Constrain the vapor pressure to be less than equal to leaf esat
       this%ceair = min( this%eair, this%esat )

       select case(this%gstype)
       case (VAR_SCM_BBERRY)

          call GsBallBerry(this)

          do idof = 1, this%ndof
             if (this%gs(idof) > 0.d0) then
                this%gleaf_c(idof) = 1.d0/(1.0d0/this%gbc + 1.6d0/this%gs(idof))
                this%gleaf_w(idof) = 1.d0/(1.0d0/this%gbv + 1.0d0/this%gs(idof))
             else
                this%gleaf_c(idof) = 0.d0
                this%gleaf_w(idof) = 0.d0
             end if
          enddo

      case (VAR_SCM_MEDLYN)

          call GsMedlyn(this)

          do idof = 1, this%ndof
             if (this%gs(idof) > 0.d0) then
                this%gleaf_c(idof) = 1.d0/(1.0d0/this%gbc + 1.6d0/this%gs(idof))
                this%gleaf_w(idof) = 1.d0/(1.0d0/this%gbv + 1.0d0/this%gs(idof))
             else
                this%gleaf_c(idof) = 0.d0
                this%gleaf_w(idof) = 0.d0
             end if
          enddo

       end select

    else
       write(iulog,*)'PhotosynthesisAuxVarCompute: Add code when dapi = 0.d0'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine PhotosynthesisAuxVarCompute_SemiEmpirical

    !------------------------------------------------------------------------
  subroutine PhotosynthesisAuxVarCompute_WUE(this)
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
    PetscReal                            :: t1, t2, t3, a, b, dy, head
    PetscReal                            :: an, gbc, delta_c
    PetscInt                             :: ileaf, idof
    class(plant_auxvar_type) , pointer   :: plant
    PetscReal                , parameter :: g = 9.80665d0
    PetscReal                , parameter :: denh2o = 1000.d0

    if (this%pathway_and_stomatal_params_defined  == 0) then
       call SetPathwayParameters(this)
       call SetStomatalConductanceParameters(this)
       call SetPlantParameters(this)
       this%pathway_and_stomatal_params_defined = 1
    end if

    if (this%dpai > 0.d0) then

       select case(this%c3psn)
       case (VAR_PHOTOSYNTHETIC_PATHWAY_C4)
          write(*,*)'PhotosynthesisAuxVarCompute2 not implemented for C4'
          call endrun(msg=errMsg(__FILE__, __LINE__))

       case (VAR_PHOTOSYNTHETIC_PATHWAY_C3)

          call C3_Temperature_Response(this)
          call Compute_Electron_Transport_Rate(this)
          call C3_Net_Assimilation_From_Gs(this)

       end select

       ! CO2 at leaf surface
       do idof = 1, this%ndof
          this%cs = max(this%cair - this%an(idof)/this%gbc, 1.d0)
       end do

       ! Saturation vapor pressure at leaf temperature
       call SatVap (this%tleaf, this%esat, this%desat)

       ! Constrain eair >= 0.05*esat[tleaf] so that solution does not blow up. This 
       ! ensures that hs does not go to zero. Also eair <= esat[tleaf] so that hs <= 1. 
       this%ceair = min( max(this%eair, 0.20d0*this%esat), this%esat )

       ! Constrain the vapor pressure to be less than equal to leaf esat
       this%ceair = min( this%eair, this%esat )

       select case(this%gstype)
       case (VAR_SCM_WUE, VAR_SCM_MANZONI11)
          do idof = 1, this%ndof
             this%hs = (this%gbv * this%eair + this%gs(idof) * this%esat)/((this%gbv + this%gs(idof)) * this%esat)
             this%vpd = max((this%esat - this%hs * this%esat), 0.1d0)/this%pref
          end do

       case (VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14)

          idof = 1
          this%hs = (this%gbv * this%eair + this%gs(idof) * this%esat)/((this%gbv + this%gs(idof)) * this%esat)
          this%vpd = max((this%esat - this%hs * this%esat), 0.1d0)/this%pref

       end select

    else
       write(iulog,*)'PhotosynthesisAuxVarCompute: Add code when dapi = 0.d0'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

  end subroutine PhotosynthesisAuxVarCompute_WUE

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
    !
    PetscInt :: idof

    do idof = 1, this%ndof
       ! Rubisco-limited photosynthesis
       this%ac(idof) = this%vcmax
       this%dac_dci(idof) = 0.d0

       ! RuBP-limited photosynthesis
       this%aj(idof) = this%qe_c4 * this%apar
       this%daj_dci(idof) = 0.d0

       ! PEP carboxylase-limited (CO2-limited)
       if (this%ci(idof) > 0.d0) then
          this%ap(idof)      = this%kp * this%ci(idof)
          this%dap_dci(idof) = this%kp
       else
          this%ap(idof) = 0.d0
          this%dap_dci(idof) = 0.d0
       end if
    enddo

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
    !
    PetscReal :: a, b
    PetscInt  :: idof

    do idof = 1, this%ndof
       if (this%ci(idof) - this%cp > 0.d0) then

          ! Rubisco-limited photosynthesis
          a = this%vcmax
          b = this%kc*(1.d0 + this%o2ref/this%ko)

          this%ac(idof)      = a*(this%ci(idof) - this%cp)/(this%ci(idof) + b)
          this%dac_dci(idof) = a*(b + this%cp)/((this%ci(idof) + b)**2.d0)

          ! RuBP-limited photosynthesis
          a = this%je/4.d0
          b = 2.d0*this%cp

          this%aj(idof)      = a*(this%ci(idof) - this%cp)/(this%ci(idof) + b)
          this%daj_dci(idof) = a*(b + this%cp)/((this%ci(idof) + b)**2.d0)

       else
          this%ac(idof)      = 0.d0
          this%dac_dci(idof) = 0.d0
          this%aj(idof)      = 0.d0
          this%daj_dci(idof) = 0.d0
       end if

       ! Product-limited photosynthesis
       this%ap(idof)      = 0.d0
       this%dap_dci(idof) = 0.d0
    enddo

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
    PetscReal :: ai, dai_dci, denom
    PetscInt  :: idof

    do idof = 1, this%ndof
       select case(this%colim)
       case (1)
          aquad = this%colim_c4a
          bquad = -(this%ac(idof) + this%aj(idof))
          cquad =   this%ac(idof) * this%aj(idof)

          call quadratic(aquad, bquad, cquad, root1, root2)
          ai = min(root1, root2)
          denom   = this%ac(idof) + this%aj(idof) - 2.d0 * this%colim_c4a * ai
          dai_dci = (this%dac_dci(idof) * (this%aj(idof) - ai) + this%daj_dci(idof) * (this%ac(idof) - ai)) / denom

          aquad = this%colim_c4b
          bquad = -(ai + this%ap(idof))
          cquad =   ai * this%ap(idof)

          call quadratic(aquad, bquad, cquad, root1, root2)
          this%ag(idof) = min(root1, root2)

          if (this%ag(idof) > 0.d0) then
             denom        = ai + this%ap(idof) - 2.d0 * this%colim_c4b * this%ag(idof)
             this%dag_dci(idof) = (dai_dci * (this%ap(idof) - this%ag(idof)) + this%dap_dci(idof) * (ai - this%ag(idof)))/denom
          else
             this%dag_dci(idof) = 0.d0
          end if

       case (0)

          this%ag(idof) = min(this%ac(idof), this%aj(idof), this%ap(idof))

          if (this%ac(idof) < this%aj(idof) .and. this%ac(idof) < this%ap(idof)) then
             this%dag_dci(idof) = this%dac_dci(idof)

          elseif (this%aj(idof) < this%ac(idof) .and. this%aj(idof) < this%ap(idof)) then
             this%dag_dci(idof) = this%daj_dci(idof)

          else
             this%dag_dci(idof) = this%dap_dci(idof)
          end if

       end select

       if (this%ac(idof) < 0.d0) then
          this%ac(idof)      = 0.d0
          this%dac_dci(idof) = 0.d0
       endif

       if (this%aj(idof) < 0.d0) then
          this%aj(idof)      = 0.d0
          this%daj_dci(idof) = 0.d0
       endif

       if (this%ap(idof) < 0.d0) then
          this%ap(idof)      = 0.d0
          this%dap_dci(idof) = 0.d0
       endif

       if (this%ag(idof) < 0.d0) then
          this%ag(idof)      = 0.d0
          this%dag_dci(idof) = 0.d0
       end if

       this%an(idof)      = this%ag(idof) - this%rd
       this%dan_dci(idof) = this%dag_dci(idof)
    enddo

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
    PetscReal :: root1, root2, denom
    PetscInt  :: idof
    
    do idof = 1, this%ndof
       select case(this%colim)
       case (1)
          aquad = this%colim_c3
          bquad = -(this%ac(idof) + this%aj(idof))
          cquad =   this%ac(idof) * this%aj(idof)

          call quadratic(aquad, bquad, cquad, root1, root2)

          this%ag(idof) = min(root1, root2)
          if (this%ag(idof) > 0.d0) then
             denom        = this%ac(idof) + this%aj(idof) - 2.d0 * this%colim_c3 * this%ag(idof)
             this%dag_dci(idof) = (this%dac_dci(idof) * (this%aj(idof) - this%ag(idof)) + this%daj_dci(idof) * (this%ac(idof) - this%ag(idof))) / denom
          else
             this%dag_dci(idof) = 0.d0
          endif

       case (0)

          this%ag(idof) = min(this%ac(idof), this%aj(idof))
          if (this%ac(idof) < this%aj(idof)) then
             this%dag_dci(idof) = this%dac_dci(idof)
          else
             this%dag_dci(idof) = this%daj_dci(idof)
          end if

       end select

       if (this%ac(idof) < 0.d0) then
          this%ac(idof)      = 0.d0
          this%dac_dci(idof) = 0.d0
       endif

       if (this%aj(idof) < 0.d0) then
          this%aj(idof)      = 0.d0
          this%daj_dci(idof) = 0.d0
       endif

       if (this%ap(idof) < 0.d0) then
          this%ap(idof)      = 0.d0
          this%dap_dci(idof) = 0.d0
       endif

       if (this%ag(idof) < 0.d0) then
          this%ag(idof)      = 0.d0
          this%dag_dci(idof) = 0.d0
       end if

       this%an(idof)      = this%ag(idof) - this%rd
       this%dan_dci(idof) = this%dag_dci(idof)
    enddo

  end subroutine C3_Net_Assimilation

  !------------------------------------------------------------------------

  subroutine C3_Net_Assimilation_From_Gs(this)
    !
    use MathUtilsMod, only : quadratic
    !
    implicit none
    !
    ! !ARGUMENTS
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal                         :: a, b
    PetscReal                         :: aquad, bquad, cquad
    PetscReal                         :: root1, root2, denom
    PetscInt                          :: idof
    PetscReal, parameter              :: ci_min = 1.d0

    do idof = 1, this%ndof

       this%gleaf_c(idof) = 1.d0/(1.0d0/this%gbc + 1.6d0/this%gs(idof))
       this%gleaf_w(idof) = 1.d0/(1.0d0/this%gbv + 1.0d0/this%gs(idof))

       ! Rubisco-limited photosynthesis
       a = this%vcmax
       b = this%kc*(1.d0 + this%o2ref/this%ko)

       aquad = 1.d0/this%gleaf_c(idof)
       bquad = -(this%cair + b) - (a - this%rd)/this%gleaf_c(idof)
       cquad = a * (this%cair - this%cp) - this%rd * (this%cair + b)

       call quadratic(aquad, bquad, cquad, root1, root2)

       this%ac(idof) = min(root1, root2) + this%rd

       ! RuBP-limited photosynthesis
       a = this%je/4.d0
       b = 2.d0*this%cp

       aquad = 1.d0/this%gleaf_c(idof)
       bquad = -(this%cair + b) - (a - this%rd)/this%gleaf_c(idof)
       cquad = a * (this%cair - this%cp) - this%rd * (this%cair + b)

       call quadratic(aquad, bquad, cquad, root1, root2)

       this%aj(idof) = min(root1, root2) + this%rd

       select case(this%colim)
       case (1)
          aquad = this%colim_c3
          bquad = -(this%ac(idof) + this%aj(idof))
          cquad =   this%ac(idof) * this%aj(idof)

          call quadratic(aquad, bquad, cquad, root1, root2)

          this%ag(idof) = min(root1, root2)

       case (0)

          this%ag(idof) = min(this%ac(idof), this%aj(idof))
       end select

       if (this%ac(idof) < 0.d0) then
          this%ac(idof) = 0.d0
       endif

       if (this%aj(idof) < 0.d0) then
          this%aj(idof) = 0.d0
       endif

       if (this%ap(idof) < 0.d0) then
          this%ap(idof) = 0.d0
       endif

       if (this%ag(idof) < 0.d0) then
          this%ag(idof) = 0.d0
       end if

       ! Net photosynthesis
       this%an(idof) = this%ag(idof) - this%rd

       this%ci(idof) = this%cair - this%an(idof)/this%gleaf_c(idof)
       this%ci(idof) = max(this%ci(idof), ci_min)

    enddo

  end subroutine C3_Net_Assimilation_From_Gs

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
    PetscInt  :: idof

    this%g0 = max( this%g0opt * this%btran, gs_min )
    this%g1 = this%g1opt

    do idof = 1, this%ndof

       if (this%an(idof) > 0.d0) then

          aquad = this%cs
          bquad = this%cs   * (this%gbv - this%g0) - this%g1*this%an(idof)
          cquad = -this%gbv * (this%cs * this%g0 + this%g1 * this%an(idof) * this%ceair/this%esat)

          call quadratic (aquad, bquad, cquad, root1, root2)
          this%gs(idof) = max(root1, root2)

       else
          this%gs(idof) = this%g0
       end if
    enddo

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
    PetscInt  :: idof

    this%g0 = this%g0opt
    this%g1 = this%g1opt

    do idof = 1, this%ndof
    if (this%an(idof) > 0.d0) then

       vpd_term = max((this%esat - this%ceair), vpd_min) * 0.001d0
       term = 1.6d0 * this%an(idof) / this%cs

       aquad = 1.d0
       bquad = -(2.d0 * (this%g0 + term) + (this%g1 * term)**2.d0 / (this%gbv * vpd_term))
       cquad = this%g0 * this%g0 + (2.d0 * this%g0 + term * (1.d0 - this%g1 * this%g1 / vpd_term)) * term

       call quadratic (aquad, bquad, cquad, root1, root2)
       this%gs(idof) = max(root1,root2)

    else

       this%gs(idof) = this%g0

    end if
    enddo

  end subroutine GsMedlyn

  !------------------------------------------------------------------------
  subroutine PhotosynthesisDetermineIfSolutionIsBounded(this)
    !
    implicit none
    !
    class(photosynthesis_auxvar_type) :: this
    PetscBool                         :: bounded
    !
    PetscInt, parameter               :: idof_wue = 1
    PetscInt, parameter               :: idof_hyd = 2
    PetscReal                         :: res_wue_1, res_wue_2
    PetscReal                         :: res_hyd_1, res_hyd_2
    PetscReal                         :: gs_1, gs_2

    select case (this%gstype)
    case (VAR_SCM_WUE, VAR_SCM_MANZONI11, VAR_SCM_OSMWANG)
       this%gs(idof_wue) = gs_min_wue
       call this%AuxVarCompute()
       res_wue_1 = this%residual_wue(idof_wue)

       this%gs(idof_wue) = gs_max_wue
       call this%AuxVarCompute()
       res_wue_2 = this%residual_wue(idof_wue)

       if (res_wue_1 * res_wue_2 > 0.d0) then
          this%soln_is_bounded(idof_wue) = PETSC_FALSE
       else
          this%soln_is_bounded(idof_wue) = PETSC_TRUE
       end if

    case (VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14)
    
          ! Residual at minimum gs
          this%gs(idof_wue) = gs_min_bonan14
          this%gs(idof_hyd) = gs_min_bonan14

          call this%AuxVarCompute()

          res_wue_1 = this%residual_wue(idof_wue)
          res_hyd_1 = this%residual_hyd(idof_hyd)

          ! Residual at maximum gs
          this%gs(idof_wue) = gs_max_wue
          this%gs(idof_hyd) = gs_max_wue

          call this%AuxVarCompute()

          res_wue_2 = this%residual_wue(idof_wue)
          res_hyd_2 = this%residual_hyd(idof_hyd)

          if ( min(res_wue_1,res_hyd_1) * min(res_wue_2,res_hyd_2) < 0.d0) then
             if ( res_wue_1 * res_wue_2 < 0.d0) then
                this%soln_is_bounded(idof_wue) = PETSC_TRUE
             else
                this%soln_is_bounded(idof_wue) = PETSC_FALSE
             end if
             if ( res_hyd_1 * res_hyd_2 < 0.d0) then
                this%soln_is_bounded(idof_hyd) = PETSC_TRUE
             else
                this%soln_is_bounded(idof_hyd) = PETSC_FALSE
             end if
          else
             this%soln_is_bounded(idof_wue) = PETSC_FALSE
             this%soln_is_bounded(idof_hyd) = PETSC_FALSE
          end if

       case (VAR_SCM_BBERRY, VAR_SCM_MEDLYN)

         if (this%c3psn == VAR_PHOTOSYNTHETIC_PATHWAY_C3) then
            this%ci(idof_wue) = 0.7d0 * this%cair
         else
            this%ci(idof_wue) = 0.4d0 * this%cair
         endif
         call this%AuxVarCompute()
         gs_1 = this%gs(1)
         res_wue_1 = this%an(idof_wue)

         res_wue_2 = -999.d0
         gs_2 = -999.d0

         if (this%an(idof_wue) < 0.d0) then
            this%soln_is_bounded(idof_wue) = PETSC_FALSE
         else
            this%ci(idof_wue) = 0.99d0 * this%cair
            call this%AuxVarCompute()
            res_wue_2 = this%an(idof_wue)
            gs_2 = this%gs(1)

            if (this%an(idof_wue) < 0.d0) then
               this%soln_is_bounded(idof_wue) = PETSC_FALSE
            else
               this%soln_is_bounded(idof_wue) = PETSC_TRUE
            endif
         endif

       case default
          write(*,*)'Unknown stomatal model'
          call endrun(msg=errMsg(__FILE__, __LINE__))

       end select

  end subroutine PhotosynthesisDetermineIfSolutionIsBounded

  !------------------------------------------------------------------------
  subroutine PhotosynthesisPreSolve(this)
    !
    use MultiPhysicsProbConstants , only : VAR_SCM_BONAN14
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    class(photosynthesis_auxvar_type) :: this
    !
    PetscReal                         :: etflx, esat, desat
    PetscInt, parameter               :: ileaf = 1

    call SatVap (this%tleaf_prev, esat, desat)

    etflx = (esat + desat * (this%tleaf - this%tleaf_prev) - this%eair)/this%pref * this%gleaf_w_soln * this%fdry
    if (this%gstype == VAR_SCM_MODIFIED_BONAN14) then
       call ComputePsi_ModifiedBonan14(this%plant, etflx, this%plant%leaf_psi(ileaf))
    else
       call ComputeChangeInPsi(this%plant, etflx)
       this%plant%leaf_psi(ileaf) = this%plant%leaf_psi(ileaf) + this%plant%dpsi_leaf(ileaf)
    end if

  end subroutine PhotosynthesisPreSolve

  !------------------------------------------------------------------------
  subroutine PhotosynthesisPostSolve(this)
    !
   use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    class(photosynthesis_auxvar_type)    :: this
    !
    class(plant_auxvar_type) , pointer   :: plant
    PetscInt                 , parameter :: idof_wue = 1
    PetscInt                 , parameter :: idof_hyd = 2
    PetscInt                             :: ileaf, idof
    PetscReal                            :: a, b, head
    PetscReal                , parameter :: g = 9.80665d0
    PetscReal                , parameter :: denh2o = 1000.d0

    select case(this%gstype)
    case (VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14)

       ! Determine solution for which DOF (wue or hydrualics) is the solution
       ! based on co-optimization
       if ( this%soln_is_bounded(idof_wue) .and. &
            this%soln_is_bounded(idof_hyd)) then

          if (this%gs(idof_wue) < this%gs(idof_hyd)) then
             idof = idof_wue
          else
             idof = idof_hyd
          end if

       elseif (      this%soln_is_bounded(idof_wue) .and. &
            (.not.this%soln_is_bounded(idof_hyd))) then

          idof = idof_wue

       elseif ((.not.this%soln_is_bounded(idof_wue)) .and. &
            (     this%soln_is_bounded(idof_hyd))) then

          idof = idof_wue

       else
          idof = idof_wue
       end if
    case (VAR_SCM_BBERRY, VAR_SCM_MANZONI11, VAR_SCM_MEDLYN, VAR_SCM_WUE, VAR_SCM_OSMWANG)
       idof = idof_wue
    case default
        write(*,*)'Unknown stomatal conductance model'
        call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    this%ac_soln      = this%ac(idof)
    this%aj_soln      = this%aj(idof)
    this%ap_soln      = this%ap(idof)
    this%ag_soln      = this%ag(idof)
    this%an_soln      = this%an(idof)
    this%gs_soln      = this%gs(idof)
    this%gleaf_c_soln = this%gleaf_c(idof)
    this%gleaf_w_soln = this%gleaf_w(idof)
    if (this%gs_soln > 0.d0) then
       this%gleaf_c_soln = 1.d0/(1.0d0/this%gbc + 1.6d0/this%gs_soln)
       this%gleaf_w_soln = 1.d0/(1.0d0/this%gbv + 1.0d0/this%gs_soln)
    end if

  end subroutine PhotosynthesisPostSolve
#endif

end module PhotosynthesisAuxType

