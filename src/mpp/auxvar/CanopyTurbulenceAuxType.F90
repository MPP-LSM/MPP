module CanopyTurbulenceAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg

  implicit none
  private

  type, public :: canopy_turbulence_auxvar_type

     PetscInt :: ncair
     PetscInt :: ncan_lev

     PetscInt , pointer :: ntop(:)      ! Index of top canopy layer
     PetscReal, pointer :: hc(:)        ! Canopy height (m)

     PetscReal, pointer :: zref(:)      ! Reference height (m)
     PetscReal, pointer :: pref(:)      ! Atmospheric pressure (Pa)
     PetscReal, pointer :: uref(:)      ! Wind speed at reference height (m/s)
     PetscReal, pointer :: vref(:)      ! Water vapor at reference height (mol/mol)
     PetscReal, pointer :: tref(:)      ! Air temperature at reference height (K)
     PetscReal, pointer :: rhref(:)     ! Relative humidity at reference height (%)
     !PetscReal, pointer :: shref(:)     ! Specific humidity at reference height (kg/kg)

     PetscReal, pointer :: ucan(:)      ! Wind speed at canopy top (m/s)
     PetscReal, pointer :: vcan(:)      ! Water vapor at canopy top (mol/mol)
     PetscReal, pointer :: tcan(:)      ! Air temperature at canopy top (K)

     PetscReal, pointer :: rhomol(:)    ! Molar density (mol/m3)
     PetscReal, pointer :: rhoair(:)    ! Air density at reference height (kg/m3)
     PetscReal, pointer :: cpair(: )    ! Specific heat of air at constant pressure, at reference height (J/mol/K)
     PetscReal, pointer :: mmair(:)     ! Molecular mass of air at reference height (kg/mol)
     PetscReal, pointer :: thref(:)     ! Potential temperature at reference height (K)
     PetscReal, pointer :: thvref(:)    ! Virtual potential temperature at reference height (K)

     PetscReal, pointer :: Lc(:)        ! Canopy density length scale (m)

     PetscReal, pointer :: c1m(:)       ! Roughness sublayer c1 parameter for momentum (dimensionless)
     PetscReal, pointer :: c1c(:)       ! Roughness sublayer c1 parameter for scalars (dimensionless)
     PetscReal, pointer :: c2(:)        ! Roughness sublayer depth scale multiplier (dimensionless)
     PetscReal, pointer :: disp(:)      ! Displacement height (m)
     PetscReal, pointer :: beta(:)      ! u* / u(hc)
     PetscReal, pointer :: PrSc(:)      ! Prandtl (Schmidt) number at canopy top
     PetscReal, pointer :: ustar(:)     ! Friction velocity (m/s)
     PetscReal, pointer :: tstar(:)     ! Temperature scale (K)
     PetscReal, pointer :: vstar(:)     ! Water vapor scale (mol/mol)
     PetscReal, pointer :: gac(:)       ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
     PetscReal, pointer :: obu_ustar(:) ! Obukhov length used for u* (m)
     PetscReal, pointer :: obu(:)       ! Value for Obukhov length (m)

     PetscReal, pointer :: pai(:)       ! Canopy plant area index (m2/m2): LAI + SAI

     PetscReal, pointer :: tksoi(:)     ! Soil thermal conductivity (W/m/K)
     PetscReal, pointer :: dzsoi(:)     ! Soil layer depth (m)
     PetscReal, pointer :: tsoi(:)      ! Soil temperature (K)
     PetscReal, pointer :: ressoi(:)    ! Soil evaporative resistance (s/m)
     PetscReal, pointer :: rhgsoi(:)    ! Relative humidity of airspace at soil surface (fraction)
     PetscReal, pointer :: rnsoi(:)     ! Net radiation at ground (W/m2)

     PetscReal, pointer :: zs(:,:)
     PetscReal, pointer :: wind(:,:)
     PetscReal, pointer :: ga_prof(:,:)
   contains
     procedure, public :: Init                    => CAirTurbInit
     procedure, public :: ComputeDerivedAtmInputs => CAirTurbComputeDerivedAtmInputs
  end type canopy_turbulence_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine CAirTurbInit(this, ncair)
    !
    implicit none
    !
    class (canopy_turbulence_auxvar_type) :: this
    PetscInt :: ncair
    !
    PetscInt :: icair, k

    this%ncair = ncair
    this%ncan_lev = 93

    allocate(this%ntop      (ncair)); this%ntop(:) = 43;
    allocate(this%hc        (ncair))

    allocate(this%zref      (ncair))

    allocate(this%pref      (ncair))
    allocate(this%uref      (ncair))
    allocate(this%vref      (ncair))
    allocate(this%tref      (ncair))
    allocate(this%rhref     (ncair))
    !allocate(this%shref     (ncair))

    allocate(this%ucan      (ncair))
    allocate(this%vcan      (ncair))
    allocate(this%tcan      (ncair))

    allocate(this%rhomol    (ncair))
    allocate(this%rhoair    (ncair))
    allocate(this%mmair     (ncair))
    allocate(this%cpair     (ncair))
    allocate(this%thref     (ncair))
    allocate(this%thvref    (ncair))

    allocate(this%Lc        (ncair))

    allocate(this%c1m       (ncair))
    allocate(this%c1c       (ncair))
    allocate(this%c2        (ncair))
    allocate(this%disp      (ncair))
    allocate(this%beta      (ncair))
    allocate(this%PrSc      (ncair))
    allocate(this%ustar     (ncair))
    allocate(this%tstar     (ncair))
    allocate(this%vstar     (ncair))
    allocate(this%gac       (ncair))
    allocate(this%obu_ustar (ncair))
    allocate(this%obu       (ncair))

    allocate(this%pai       (ncair))

    allocate(this%tksoi     (ncair))
    allocate(this%dzsoi     (ncair))
    allocate(this%tsoi      (ncair))
    allocate(this%ressoi    (ncair))
    allocate(this%rhgsoi    (ncair))
    allocate(this%rnsoi     (ncair))

    allocate(this%zs(ncair, this%ncan_lev))
    allocate(this%wind(ncair,this%ncan_lev))
    allocate(this%ga_prof(ncair,this%ncan_lev))

    do icair = 1,ncair
       this%zs(icair,1) = 0.d0
       do k = 2,this%ncan_lev
          this%zs(icair,k) = 0.25d0 + (k-2)*0.5d0
       end do
    end do

  end subroutine CAirTurbInit

  !------------------------------------------------------------------------
  subroutine CAirTurbComputeDerivedAtmInputs(this,icair)
    !
    use WaterVaporMod             , only : SatVap
    use MultiPhysicsProbConstants , only : MM_H2O, MM_DRY_AIR, CPD, CPW, RGAS
    !
    implicit none
    !
    class (canopy_turbulence_auxvar_type) :: this
    !
    PetscReal :: esat, desatdt, eref, vref_kg_kg
    PetscInt :: icair

    call satvap (this%tref(icair), esat, desatdt);

    eref = (this%rhref(icair) / 100.d0) * esat
    this%vref(icair) = eref / this%pref(icair)

    this%rhomol(icair) = this%pref(icair) / (RGAS * this%tref(icair))
    this%rhoair(icair) = this%rhomol(icair) * MM_DRY_AIR * (1.d0 - (1.d0 - MM_H2O/MM_DRY_AIR) * eref / this%pref(icair))
    this%mmair(icair)  = this%rhoair(icair) / this%rhomol(icair)
    this%thref(icair)  = this%tref(icair) + 0.0098d0 * this%zref(icair)

    vref_kg_kg = MM_H2O/MM_DRY_AIR * eref / (this%pref(icair) - (1.d0 - MM_H2O/MM_DRY_AIR) * eref)

    this%cpair(icair)  = CPD * (1.d0 + (CPW/CPD - 1.d0) * vref_kg_kg) * this%mmair(icair)
    this%thvref(icair) = this%thref(icair) * (1.d0 + 0.61d0 * vref_kg_kg)

  end subroutine CAirTurbComputeDerivedAtmInputs

#endif

end module CanopyTurbulenceAuxType
