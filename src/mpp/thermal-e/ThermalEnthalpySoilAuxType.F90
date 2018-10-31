module ThermalEnthalpySoilAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use PorosityFunctionMod , only : porosity_params_type
  use SaturationFunction  , only : saturation_params_type
  use RichardsODEPressureAuxType, only : rich_ode_pres_auxvar_type
  !
  implicit none
  !
  private

  type, public, extends(rich_ode_pres_auxvar_type) :: therm_enthalpy_soil_auxvar_type

     private

     PetscReal                    :: ul                       ! internal energy liquid [J kmol^{-1}]
     PetscReal                    :: hl                       ! enthalpy liquid [J kmol^{-3}]

     PetscReal                    :: dul_dP                   ! [J kmol^{-1} Pa^{-1}]
     PetscReal                    :: dhl_dP                   ! [J kmol^{-3} Pa^{-1}]

     PetscReal                    :: dul_dT                   ! [J kmol^{-1} K^{-1}]
     PetscReal                    :: dhl_dT                   ! [J kmol^{-3} K^{-1}]
     PetscReal                    :: dkr_dT                   ! [K^{-1}]
     PetscReal                    :: dsat_dT                  ! [K^{-1}]

     PetscReal                    :: Kel                      ! Kersten number liquid [-]
     PetscReal                    :: therm_cond_wet           ! wet thermal conductivity [J s^{-1} m^{-1} K^{-1}]
     PetscReal                    :: therm_cond_dry           ! dry thermal conductivity [J s^{-1} m^{-1} K^{-1}]
     PetscReal                    :: therm_cond               ! thermal conductivity [J s^{-1} m^{-1} K^{-1}]
     PetscReal                    :: dtherm_cond_dP
     PetscReal                    :: dKel_dp
     PetscReal                    :: therm_alpha

     PetscReal                    :: den_soil                 ! [kg m^{-3}]
     PetscReal                    :: heat_cap_soil            ! [J kg^{-1} K^{-1}]

     PetscInt                     :: int_energy_enthalpy_type ! [-]

   contains
     procedure, public :: Init                        => ThermEnthalpyAuxVarInit
     procedure, public :: Set                         => ThermEnthalpyAuxVarSet
     procedure, public :: Get                         => ThermEnthalpyAuxVarGet
     procedure, public :: Compute                     => ThermEnthalpyAuxVarCompute     

     procedure, public :: SetIntEnergyLiquid          => ThermEnthalpyAuxVarSetIntEnergyLiquid
     procedure, public :: SetEnthalpyLiquid           => ThermEnthalpyAuxVarSetEnthalpyLiquid
     procedure, public :: SetDUlDP                    => ThermEnthalpyAuxVarSetDUlDP
     procedure, public :: SetDHlDP                    => ThermEnthalpyAuxVarSetDHlDP
     procedure, public :: SetDUlDT                    => ThermEnthalpyAuxVarSetDUlDT
     procedure, public :: SetDHlDT                    => ThermEnthalpyAuxVarSetDHlDT
     procedure, public :: SetDKrdT                    => ThermEnthalpyAuxVarSetDKrdT
     procedure, public :: SetDSatdT                   => ThermEnthalpyAuxVarSetDSatdT
     procedure, public :: SetKerstenNumberLiquid      => ThermEnthalpyAuxVarSetKerstenNumberLiquid
     procedure, public :: SetThermalCondWet           => ThermEnthalpyAuxVarSetThermalCondWet
     procedure, public :: SetThermalCondDry           => ThermEnthalpyAuxVarSetThermalCondDry
     procedure, public :: SetThermalCond              => ThermEnthalpyAuxVarSetThermalCond
     procedure, public :: SetDThermalCondDP           => ThermEnthalpyAuxVarSetDThermalCondDP
     procedure, public :: SetDKelDp                   => ThermEnthalpyAuxVarSetDKelDp
     procedure, public :: SetKerstenNumberLiquidCoeff => ThermEnthalpyAuxVarSetKerstenNumberLiquidCoeff
     procedure, public :: SetDensitySoil              => ThermEnthalpyAuxVarSetDensitySoil
     procedure, public :: SetHeatCapSoil              => ThermEnthalpyAuxVarSetHeatCapSoil
     procedure, public :: SetIntEnergyType            => ThermEnthalpyAuxVarSetIntEnergyType

     procedure, public :: GetIntEnergyLiquid          => ThermEnthalpyAuxVarGetIntEnergyLiquid
     procedure, public :: GetEnthalpyLiquid           => ThermEnthalpyAuxVarGetEnthalpyLiquid
     procedure, public :: GetDUlDP                    => ThermEnthalpyAuxVarGetDUlDP
     procedure, public :: GetDHlDP                    => ThermEnthalpyAuxVarGetDHlDP
     procedure, public :: GetDUlDT                    => ThermEnthalpyAuxVarGetDUlDT
     procedure, public :: GetDHlDT                    => ThermEnthalpyAuxVarGetDHlDT
     procedure, public :: GetDKrdT                    => ThermEnthalpyAuxVarGetDKrdT
     procedure, public :: GetDSatdT                   => ThermEnthalpyAuxVarGetDSatdT
     procedure, public :: GetKerstenNumberLiquid      => ThermEnthalpyAuxVarGetKerstenNumberLiquid
     procedure, public :: GetThermalCondWet           => ThermEnthalpyAuxVarGetThermalCondWet
     procedure, public :: GetThermalCondDry           => ThermEnthalpyAuxVarGetThermalCondDry
     procedure, public :: GetThermalCond              => ThermEnthalpyAuxVarGetThermalCond
     procedure, public :: GetDThermalCondDP           => ThermEnthalpyAuxVarGetDThermalCondDP
     procedure, public :: GetDKelDp                   => ThermEnthalpyAuxVarGetDKelDp
     procedure, public :: GetKerstenNumberLiquidCoeff => ThermEnthalpyAuxVarGetKerstenNumberLiquidCoeff
     procedure, public :: GetDensitySoil              => ThermEnthalpyAuxVarGetDensitySoil
     procedure, public :: GetHeatCapSoil              => ThermEnthalpyAuxVarGetHeatCapSoil

  end type therm_enthalpy_soil_auxvar_type

  public :: ThermEnthalpyAuxVarCopy

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use EOSWaterMod               , only : INT_ENERGY_ENTHALPY_CONSTANT
    use EOSWaterMod               , only : INT_ENERGY_ENTHALPY_IFC67
    use RichardsODEPressureAuxType, only : RichODEPressureAuxVarInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) :: this

    call RichODEPressureAuxVarInit(this)
    
    this%ul                       = 0.d0
    this%hl                       = 0.d0

    this%dul_dT                   = 0.d0
    this%dhl_dT                   = 0.d0
    this%dkr_dT                   = 0.d0
    this%dsat_dT                  = 0.d0

    this%dul_dP                   = 0.d0
    this%dhl_dP                   = 0.d0

    this%Kel                      = 0.d0
    this%therm_cond_wet           = 0.d0
    this%therm_cond_dry           = 0.d0
    this%therm_cond               = 0.d0
    this%dtherm_cond_dP           = 0.d0
    this%dKel_dp                  = 0.d0
    this%therm_alpha              = 0.d0

    call this%SetPermeability(8.3913d-12)

    this%den_soil                 = 0.d0
    this%heat_cap_soil            = 0.d0

    this%int_energy_enthalpy_type = INT_ENERGY_ENTHALPY_CONSTANT

  end subroutine ThermEnthalpyAuxVarInit

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarCopy(this, auxvar)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) :: this
    class(therm_enthalpy_soil_auxvar_type) :: auxvar

    call this%SetTemperature(auxvar%GetTemperature())
    call this%SetPressure(auxvar%GetPressure())

    this%ul                       = auxvar%ul
    this%hl                       = auxvar%hl
    call this%SetDensity(auxvar%GetDensity())
    call this%SetViscosity(auxvar%GetViscosity())
    call this%SetLiquidSaturation(auxvar%GetLiquidSaturation())
    call this%SetRelativePermeability(auxvar%GetRelativePermeability())

    this%dul_dT                   = auxvar%dul_dT
    this%dhl_dT                   = auxvar%dhl_dT
    call this%SetDDenDT(auxvar%GetDDenDT())
    call this%SetDVisDT(auxvar%GetDVisDT())
    this%dsat_dT                  = auxvar%dsat_dT
    this%dkr_dT                   = auxvar%dkr_dT

    this%dul_dP                   = auxvar%dul_dP
    this%dhl_dP                   = auxvar%dhl_dP
    call this%SetDDenDP(auxvar%GetdDenDP())
    call this%SetDVisdP(auxvar%GetDVisdP())
    call this%SetDSatdP(auxvar%GetDSatdP())
    call this%SetDKrdP(auxvar%GetDKrdP())
    call this%SetDPordP(auxvar%GetDPordP())

    this%Kel                      = auxvar%Kel
    this%therm_cond_wet           = auxvar%therm_cond_wet
    this%therm_cond_dry           = auxvar%therm_cond_dry
    this%therm_cond               = auxvar%therm_cond
    this%dtherm_cond_dP           = auxvar%dtherm_cond_dP
    this%dKel_dp                  = auxvar%dKel_dp
    call this%SetPorosity(auxvar%GetPorosity())
    call this%SetPermeabilityXYZ(auxvar%GetPermeabilityXYZ())
    this%therm_alpha              = auxvar%therm_alpha

    this%den_soil                 = auxvar%den_soil
    this%heat_cap_soil            = auxvar%heat_cap_soil

    call this%SetDensityType(auxvar%GetDensityType())
    this%int_energy_enthalpy_type = auxvar%int_energy_enthalpy_type

    call this%SetConditionValue(auxvar%GetConditionValue())

    call this%porParams%Copy(auxvar%porParams)
    call this%satParams%Copy(auxvar%satParams)

  end subroutine ThermEnthalpyAuxVarCopy

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSet(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) , intent(inout) :: this
    PetscInt                               , intent(in)    :: var_type
    PetscReal                              , intent(in)    :: variable_value

    select case(var_type)
    case (VAR_TEMPERATURE)
       call this%SetTemperature(variable_value)

    case (VAR_BC_SS_CONDITION)
       call this%SetConditionValue(variable_value)

    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermEnthalpyAuxVarSet

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarGet(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) , intent(in)  :: this
    PetscInt                               , intent(in)  :: var_type
    PetscReal                              , intent(out) :: variable_value

    select case(var_type)
    case (VAR_TEMPERATURE)
       variable_value = this%GetTemperature()

    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine ThermEnthalpyAuxVarGet

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    !
    use EOSWaterMod               , only : Density
    use EOSWaterMod               , only : Viscosity
    use PorosityFunctionMod       , only : PorosityFunctionComputation
    use SaturationFunction        , only : SatFunc_PressToSat
    use SaturationFunction        , only : SatFunc_PressToRelPerm
    use EOSWaterMod               , only : InternalEnergyAndEnthalpy
    use MultiPhysicsProbConstants , only : FMWH2O, PRESSURE_REF
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal            :: pressure, temperature
    PetscReal            :: sat, dsat_dP
    PetscReal            :: kr, dkr_dP
    PetscReal            :: por, dpor_dP
    PetscReal            :: den, dden_dP, dden_dT
    PetscReal            :: vis, dvis_dP, dvis_dT
    PetscReal, parameter :: frac_liq_sat = 1.d0

    pressure = this%GetPressure()
    temperature = this%GetTemperature()

    ! Compute saturation
    call SatFunc_PressToSat(this%satParams , pressure, sat, dsat_dP)
    call this%SetLiquidSaturation(sat)
    call this%SetDSatDP(dsat_dP)

    ! Compute relative permeability
    call SatFunc_PressToRelPerm(this%satParams, pressure, frac_liq_sat, &
         kr, dkr_dP)
    call this%SetRelativePermeability(kr)
    call this%SetDKrDP(dkr_dP)

    ! Compute porosity
    call PorosityFunctionComputation(this%porParams, pressure, por, dpor_dP)
    call this%SetPorosity(por)
    call this%SetDPorDP(dpor_dP)

    if (pressure < PRESSURE_REF)  pressure = PRESSURE_REF

    ! Compute density
    call Density(pressure, temperature, this%GetDensityType(), &
         den, dden_dP, dden_dT)
    call this%setDensity(den)
    call this%SetDDenDP(dden_dP)
    call this%SetDDenDT(dden_dT)

    ! Compute viscosity
    call Viscosity(pressure, temperature, vis, dvis_dP, dvis_dT)
    call this%SetViscosity(vis)
    call this%SetDVisDP(dvis_dP)
    call this%SetDVisDT(dvis_dT)

    ! Compute internal energy and enthalpy
    call InternalEnergyAndEnthalpy(pressure, temperature, &
         this%int_energy_enthalpy_type, den*FMWH2O, &
         this%GetDDenDT()*FMWH2O, dden_dP*FMWH2O, &
         this%ul, this%hl, this%dul_dT, this%dhl_dT, &
         this%dul_dP, this%dhl_dP)

    this%Kel        = (this%GetLiquidSaturation() + 1.d-6 )**(this%therm_alpha)
    this%dKel_dp    = this%therm_alpha*(this%GetLiquidSaturation() + 1.d-6)**(this%therm_alpha - 1.d0)* &
                      dsat_dp

    this%therm_cond     = this%therm_cond_wet*this%Kel + &
                          this%therm_cond_dry*(1.d0 - this%Kel)
    this%dtherm_cond_dP = (this%therm_cond_wet - this%therm_cond_dry)*this%dKel_dp


  end subroutine ThermEnthalpyAuxVarCompute

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetIntEnergyLiquid(this, val)
    !
    ! !DESCRIPTION:
    ! Set internal energy liquid
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%ul = val
 
  end subroutine ThermEnthalpyAuxVarSetIntEnergyLiquid
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetEnthalpyLiquid(this, val)
    !
    ! !DESCRIPTION:
    ! Set enthalpy liquid
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%hl = val
 
  end subroutine ThermEnthalpyAuxVarSetEnthalpyLiquid
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDUlDP(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of internal energy liquid w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dul_dP = val
 
  end subroutine ThermEnthalpyAuxVarSetDUlDP
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDHlDP(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of enthalpy liquid w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dhl_dP = val
 
  end subroutine ThermEnthalpyAuxVarSetDHlDP
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDUlDT(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of internal energy liquid w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dul_dT = val
 
  end subroutine ThermEnthalpyAuxVarSetDUlDT
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDHlDT(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of enthalpy liquid w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dhl_dT = val
 
  end subroutine ThermEnthalpyAuxVarSetDHlDT
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDKrdT(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of relative permeability w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dkr_dT = val
 
  end subroutine ThermEnthalpyAuxVarSetDKrdT
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDSatdT(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of saturation w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dsat_dT = val
 
  end subroutine ThermEnthalpyAuxVarSetDSatdT
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetKerstenNumberLiquid(this, val)
    !
    ! !DESCRIPTION:
    ! Set Kersten number liquid
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%Kel = val
 
  end subroutine ThermEnthalpyAuxVarSetKerstenNumberLiquid
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetThermalCondWet(this, val)
    !
    ! !DESCRIPTION:
    ! Set wet thermal conductivity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%therm_cond_wet = val
 
  end subroutine ThermEnthalpyAuxVarSetThermalCondWet
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetThermalCondDry(this, val)
    !
    ! !DESCRIPTION:
    ! Set dry thermal conductivity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%therm_cond_dry = val
 
  end subroutine ThermEnthalpyAuxVarSetThermalCondDry
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetThermalCond(this, val)
    !
    ! !DESCRIPTION:
    ! Set effective thermal conductivity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%therm_cond = val
 
  end subroutine ThermEnthalpyAuxVarSetThermalCond
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDThermalCondDP(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of effective thermal conductivity w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dtherm_cond_dP = val
 
  end subroutine ThermEnthalpyAuxVarSetDThermalCondDP
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDKelDp(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of Kersten number liquid w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dKel_dp = val
 
  end subroutine ThermEnthalpyAuxVarSetDKelDp
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetKerstenNumberLiquidCoeff(this, val)
    !
    ! !DESCRIPTION:
    ! Set exponent value for computing Kersten number parameterization
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%therm_alpha = val
 
  end subroutine ThermEnthalpyAuxVarSetKerstenNumberLiquidCoeff
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetDensitySoil(this, val)
    !
    ! !DESCRIPTION:
    ! Set density of soil
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%den_soil = val
 
  end subroutine ThermEnthalpyAuxVarSetDensitySoil
 
  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetHeatCapSoil(this, val)
    !
    ! !DESCRIPTION:
    ! Set heat capacity of soil
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%heat_cap_soil = val
 
  end subroutine ThermEnthalpyAuxVarSetHeatCapSoil

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSetIntEnergyType(this, val)
    !
    ! !DESCRIPTION:
    ! Set internal energy type
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscInt                         :: val
 
    this%int_energy_enthalpy_type = val
 
  end subroutine ThermEnthalpyAuxVarSetIntEnergyType

  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetIntEnergyLiquid(this)
    !
    ! !DESCRIPTION:
    ! Get internal energy liquid
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetIntEnergyLiquid
 
    ThermEnthalpyAuxVarGetIntEnergyLiquid = this%ul
 
  end function ThermEnthalpyAuxVarGetIntEnergyLiquid
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetEnthalpyLiquid(this)
    !
    ! !DESCRIPTION:
    ! Get enthalpy liquid
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetEnthalpyLiquid
 
    ThermEnthalpyAuxVarGetEnthalpyLiquid = this%hl
 
  end function ThermEnthalpyAuxVarGetEnthalpyLiquid
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDUlDP(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of internal energy liquid w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDUlDP
 
    ThermEnthalpyAuxVarGetDUlDP = this%dul_dP
 
  end function ThermEnthalpyAuxVarGetDUlDP
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDHlDP(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of enthalpy liquid w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDHlDP
 
    ThermEnthalpyAuxVarGetDHlDP = this%dhl_dP
 
  end function ThermEnthalpyAuxVarGetDHlDP
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDUlDT(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of internal energy liquid w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDUlDT
 
    ThermEnthalpyAuxVarGetDUlDT = this%dul_dT
 
  end function ThermEnthalpyAuxVarGetDUlDT
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDHlDT(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of enthalpy liquid w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDHlDT
 
    ThermEnthalpyAuxVarGetDHlDT = this%dhl_dT
 
  end function ThermEnthalpyAuxVarGetDHlDT
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDKrdT(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of relative permeability w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDKrdT
 
    ThermEnthalpyAuxVarGetDKrdT = this%dkr_dT
 
  end function ThermEnthalpyAuxVarGetDKrdT
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDSatdT(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of saturation w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDSatdT
 
    ThermEnthalpyAuxVarGetDSatdT = this%dsat_dT
 
  end function ThermEnthalpyAuxVarGetDSatdT
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetKerstenNumberLiquid(this)
    !
    ! !DESCRIPTION:
    ! Get Kersten number liquid
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetKerstenNumberLiquid
 
    ThermEnthalpyAuxVarGetKerstenNumberLiquid = this%Kel
 
  end function ThermEnthalpyAuxVarGetKerstenNumberLiquid
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetThermalCondWet(this)
    !
    ! !DESCRIPTION:
    ! Get wet thermal conductivity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetThermalCondWet
 
    ThermEnthalpyAuxVarGetThermalCondWet = this%therm_cond_wet
 
  end function ThermEnthalpyAuxVarGetThermalCondWet
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetThermalCondDry(this)
    !
    ! !DESCRIPTION:
    ! Get dry thermal conductivity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetThermalCondDry
 
    ThermEnthalpyAuxVarGetThermalCondDry = this%therm_cond_dry
 
  end function ThermEnthalpyAuxVarGetThermalCondDry
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetThermalCond(this)
    !
    ! !DESCRIPTION:
    ! Get effective thermal conductivity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetThermalCond
 
    ThermEnthalpyAuxVarGetThermalCond = this%therm_cond
 
  end function ThermEnthalpyAuxVarGetThermalCond
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDThermalCondDP(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of effective thermal conductivity w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDThermalCondDP
 
    ThermEnthalpyAuxVarGetDThermalCondDP = this%dtherm_cond_dP
 
  end function ThermEnthalpyAuxVarGetDThermalCondDP
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDKelDp(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of Kersten number liquid w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDKelDp
 
    ThermEnthalpyAuxVarGetDKelDp = this%dKel_dp
 
  end function ThermEnthalpyAuxVarGetDKelDp
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetKerstenNumberLiquidCoeff(this)
    !
    ! !DESCRIPTION:
    ! Get exponent value for computing Kersten number parameterization
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetKerstenNumberLiquidCoeff
 
    ThermEnthalpyAuxVarGetKerstenNumberLiquidCoeff = this%therm_alpha
 
  end function ThermEnthalpyAuxVarGetKerstenNumberLiquidCoeff
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetDensitySoil(this)
    !
    ! !DESCRIPTION:
    ! Get density of soil
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetDensitySoil
 
    ThermEnthalpyAuxVarGetDensitySoil = this%den_soil
 
  end function ThermEnthalpyAuxVarGetDensitySoil
 
  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetHeatCapSoil(this)
    !
    ! !DESCRIPTION:
    ! Get heat capacity of soil
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermEnthalpyAuxVarGetHeatCapSoil
 
    ThermEnthalpyAuxVarGetHeatCapSoil = this%heat_cap_soil
 
  end function ThermEnthalpyAuxVarGetHeatCapSoil

  !------------------------------------------------------------------------
  function ThermEnthalpyAuxVarGetIntEnergyType(this)
    !
    ! !DESCRIPTION:
    ! Get internal energy type
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_enthalpy_soil_auxvar_type) :: this
    PetscInt                         :: ThermEnthalpyAuxVarGetIntEnergyType
 
    ThermEnthalpyAuxVarGetIntEnergyType = this%int_energy_enthalpy_type
 
  end function ThermEnthalpyAuxVarGetIntEnergyType


#endif

end module ThermalEnthalpySoilAuxType
