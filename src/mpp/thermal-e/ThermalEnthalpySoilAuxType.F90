module ThermalEnthalpySoilAuxType

#ifdef USE_PETSC_LIB

  ! !USES:
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use PorosityFunctionMod , only : porosity_params_type
  use SaturationFunction  , only : saturation_params_type
  !
  implicit none
  !
  private

#include "finclude/petscsys.h"

  type, public :: therm_enthalpy_soil_auxvar_type

     PetscReal                    :: temperature     ! [Kelvins]

     PetscReal                    :: pressure        ! [Pa]

     PetscReal                    :: ul              ! internal energy liquid [J kmol^{-1}]
     PetscReal                    :: hl              ! enthalpy liquid [J kmol^{-3}]
     PetscReal                    :: denl            ! density liquid [kmol m^{-3}]
     PetscReal                    :: vis             ! visocity liquid [Pa s]
     PetscReal                    :: satl            ! saturation liquid [-]
     PetscReal                    :: krl             ! relative permeability liquid [-]
     
     PetscReal                    :: dul_dT          ! [J kmol^{-1} K^{-1}]
     PetscReal                    :: dhl_dT          ! [J kmol^{-3} K^{-1}]
     PetscReal                    :: ddenl_dT        ! [kmol m^{-3} K^{-1}]
     PetscReal                    :: dvis_dT         ! [Pa s K^{-1}]
     PetscReal                    :: dsatl_dT        ! [K^{-1}]
     PetscReal                    :: dkrl_dT         ! [K^{-1}]

     PetscReal                    :: Kel             ! Kersten number liquid [-]
     PetscReal                    :: therm_cond_wet      ! wet thermal conductivity [J s^{-1} m^{-3} K^{-1}]
     PetscReal                    :: therm_cond_dry      ! dry thermal conductivity [J s^{-1} m^{-3} K^{-1}]
     PetscReal                    :: therm_cond          ! thermal conductivity [J s^{-1} m^{-3} K^{-1}]
     PetscReal                    :: por             ! [-]
     PetscReal                    :: therm_alpha

     PetscReal                    :: den_soil        ! [kg m^{-3}]
     PetscReal                    :: heat_cap_soil   ! [J kg^{-1} K^{-1}]

     PetscInt                     :: density_type    ! [-]

     ! If the auxvar corresponds to boundary condition
     ! or source sink, the value is stored in this variable
     PetscReal                    :: condition_value !

     type(porosity_params_type)   :: porParams
     type(saturation_params_type) :: satParams

   contains
     procedure, public :: Init    => ThermEnthalpyAuxVarInit
     procedure, public :: Copy    => ThermEnthalpyAuxVarCopy
     procedure, public :: Set     => ThermEnthalpyAuxVarSet
     procedure, public :: Get     => ThermEnthalpyAuxVarGet
     procedure, public :: Compute => ThermEnthalpyAuxVarCompute     
  end type therm_enthalpy_soil_auxvar_type
  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use EOSWaterMod         , only : DENSITY_CONSTANT
    use PorosityFunctionMod , only : PorosityFunctionInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) :: this

    this%temperature     = 0.d0
    this%pressure        = 0.d0

    this%ul              = 0.d0
    this%hl              = 0.d0
    this%denl            = 0.d0
    this%vis             = 0.d0
    this%satl            = 0.d0
    this%krl             = 0.d0

    this%dul_dT          = 0.d0
    this%dhl_dT          = 0.d0
    this%ddenl_dT        = 0.d0
    this%dvis_dT         = 0.d0
    this%dsatl_dT        = 0.d0
    this%dkrl_dT         = 0.d0

    this%Kel             = 0.d0
    this%therm_cond_wet      = 0.d0
    this%therm_cond_dry      = 0.d0
    this%therm_cond          = 0.d0
    this%therm_alpha           = 0.d0

    this%den_soil        = 0.d0
    this%heat_cap_soil   = 0.d0

    this%density_type    = DENSITY_CONSTANT
    this%condition_value = 0.d0

    call PorosityFunctionInit(this%porParams)

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

    this%temperature     = auxvar%temperature
    this%pressure        = auxvar%pressure

    this%ul              = auxvar%ul
    this%hl              = auxvar%hl
    this%denl            = auxvar%denl
    this%vis             = auxvar%vis
    this%satl            = auxvar%satl
    this%krl             = auxvar%krl

    this%dul_dT          = auxvar%dul_dT
    this%dhl_dT          = auxvar%dhl_dT
    this%ddenl_dT        = auxvar%ddenl_dT
    this%dvis_dT         = auxvar%dvis_dT
    this%dsatl_dT        = auxvar%dsatl_dT
    this%dkrl_dT         = auxvar%dkrl_dT

    this%Kel             = auxvar%Kel
    this%therm_cond_wet      = auxvar%therm_cond_wet
    this%therm_cond_dry      = auxvar%therm_cond_dry
    this%therm_cond          = auxvar%therm_cond
    this%therm_alpha           = auxvar%therm_alpha

    this%den_soil        = auxvar%den_soil
    this%heat_cap_soil   = auxvar%heat_cap_soil

    this%density_type    = auxvar%density_type
    this%condition_value = auxvar%condition_value

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
       this%temperature     = variable_value

    case (VAR_BC_SS_CONDITION)
       this%condition_value = variable_value

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
       variable_value = this%temperature

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
    use MultiPhysicsProbConstants , only : FMWH2O
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal :: dsat_dP
    PetscReal :: dkr_dP
    PetscReal :: dden_dP
    PetscReal :: dvis_dP
    PetscReal :: dpor_dP
    PetscReal, parameter :: frac_liq_sat = 1.d0

    ! Compute saturation
    call SatFunc_PressToSat(this%satParams , this%pressure, this%satl, dsat_dP)

    ! Compute relative permeability
    call SatFunc_PressToRelPerm(this%satParams, this%pressure, frac_liq_sat, &
         this%krl, dkr_dP)

    ! Compute density
    call Density(this%pressure, this%temperature, this%density_type, &
         this%denl, dden_dP, this%ddenl_dT)

    ! Compute viscosity
    call Viscosity(this%pressure, this%temperature, this%vis, dvis_dP)

    ! Compute porosity
    call PorosityFunctionComputation(this%porParams, this%pressure, this%por, dpor_dP)

    ! Compute internal energy and enthalpy
    call InternalEnergyAndEnthalpy(this%pressure, this%temperature, this%denl*FMWH2O, &
         this%ddenl_dT*FMWH2O, this%ul, this%hl, this%dul_dT, this%dhl_dT)

    this%Kel        = (this%satl)**(this%therm_alpha)
    this%therm_cond = this%therm_cond_wet*this%Kel + &
                      this%therm_cond_dry*(1.d0 - this%Kel)

  end subroutine ThermEnthalpyAuxVarCompute

#endif

end module ThermalEnthalpySoilAuxType
