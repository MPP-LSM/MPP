module RichardsODEPressureAuxType

#ifdef USE_PETSC_LIB

  ! !USES:
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use PorosityFunctionMod , only : porosity_params_type
  use SaturationFunction  , only : saturation_params_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include <petsc/finclude/petsc.h>

  type, public :: rich_ode_pres_auxvar_type

     ! primary unknown independent variable
     PetscReal, private :: pressure               ! [Pa]

     PetscReal, private :: pressure_prev          ! [Pa]

     ! independent variable available from:
     !  - another governing equation, or
     !  - another system-of-equation
     PetscReal, private :: temperature            ! [K]
     PetscReal, private :: frac_liq_sat           ! [-]

     ! If the auxvar corresponds to boundary condition
     ! or source sink, the value is stored in this
     ! variable
     PetscReal, private :: condition_value        ! Depends

     ! parameters
     PetscReal, private :: perm(3)                ! [m^2]
     PetscReal, private :: por                    ! [-]
     PetscInt, private  :: density_type           ! [-]
     PetscReal, private :: pot_mass_sink_pressure ! [Pa]
     PetscReal, private :: pot_mass_sink_exponent ! [-]

     ! derived quantities = f(state_variables, parameters)
     PetscReal, private :: vis                    ! [Pa s]
     PetscReal, private :: kr                     ! [-]
     PetscReal, private :: sat                    ! [-]
     PetscReal, private :: den                    ! [kg m^{-3}]

     PetscReal, private :: dpor_dP                ! [Pa^{-1}]
     PetscReal, private :: dvis_dP                ! [s]
     PetscReal, private :: dkr_dP                 ! [Pa^{-1}]
     PetscReal, private :: dsat_dP                ! [Pa^{-1}]
     PetscReal, private :: dden_dP                ! [kg m^{-3} Pa^{-1}]

     PetscReal, private :: dvis_dT                ! [Pa s K^{-1}]
     PetscReal, private :: dden_dT                ! [kmol m^{-3} K^{-1}]

     type(porosity_params_type)   :: porParams
     type(saturation_params_type) :: satParams

   contains
     procedure, public :: Init                    => RichODEPressureAuxVarInit
     procedure, public :: SetValue                => RichODEPressureAuxVarSetValue
     procedure, public :: GetValue                => RichODEPressureAuxVarGetValue
     procedure, public :: AuxVarCompute           => RichODEPressureAuxVarCompute
     procedure, public :: SetPressure             => RichODEPressureAuxVarSetPressure
     procedure, public :: SetPressurePrev         => RichODEPressureAuxVarSetPressurePrev
     procedure, public :: SetTemperature          => RichODEPressureAuxVarSetTemperature
     procedure, public :: SetFracLiqSat           => RichODEPressureAuxVarSetFracLiqSat
     procedure, public :: SetConditionValue       => RichODEPressureAuxVarSetConditionValue
     procedure, public :: SetPermeability         => RichODEPressureAuxVarSetPermeability
     procedure, public :: SetPermeabilityX        => RichODEPressureAuxVarSetPermeabilityX
     procedure, public :: SetPermeabilityY        => RichODEPressureAuxVarSetPermeabilityY
     procedure, public :: SetPermeabilityZ        => RichODEPressureAuxVarSetPermeabilityZ
     procedure, public :: SetPermeabilityXYZ      => RichODEPressureAuxVarSetPermeabilityXYZ
     procedure, public :: SetPorosity             => RichODEPressureAuxVarSetPorosity
     procedure, public :: SetDensityType          => RichODEPressureAuxVarSetDensityType
     procedure, public :: SetPotMassSinkPressure  => RichODEPressureAuxVarSetPotMassSinkPressure
     procedure, public :: SetPotMassSinkExponent  => RichODEPressureAuxVarSetPotMassSinkExponent
     procedure, public :: SetViscosity            => RichODEPressureAuxVarSetViscosity
     procedure, public :: SetRelativePermeability => RichODEPressureAuxVarSetRelativePermeability
     procedure, public :: SetLiquidSaturation     => RichODEPressureAuxVarSetLiquidSaturation
     procedure, public :: SetDensity              => RichODEPressureAuxVarSetDensity
     procedure, public :: SetDPorDP               => RichODEPressureAuxVarSetDPorDP
     procedure, public :: SetDVisDP               => RichODEPressureAuxVarSetDVisDP
     procedure, public :: SetDKrDP                => RichODEPressureAuxVarSetDKrDP
     procedure, public :: SetDSatDP               => RichODEPressureAuxVarSetDSatDP
     procedure, public :: SetDDenDP               => RichODEPressureAuxVarSetDDenDP
     procedure, public :: SetDVisDT               => RichODEPressureAuxVarSetDVisDT
     procedure, public :: SetDDenDT               => RichODEPressureAuxVarSetDDenDT
     procedure, public :: GetPressure             => RichODEPressureAuxVarGetPressure
     procedure, public :: GetPressurePrev         => RichODEPressureAuxVarGetPressurePrev
     procedure, public :: GetTemperature          => RichODEPressureAuxVarGetTemperature
     procedure, public :: GetFracLiqSat           => RichODEPressureAuxVarGetFracLiqSat
     procedure, public :: GetConditionValue       => RichODEPressureAuxVarGetConditionValue
     procedure, public :: GetPermeabilityX        => RichODEPressureAuxVarGetPermeabilityX
     procedure, public :: GetPermeabilityY        => RichODEPressureAuxVarGetPermeabilityY
     procedure, public :: GetPermeabilityZ        => RichODEPressureAuxVarGetPermeabilityZ
     procedure, public :: GetPermeabilityXYZ      => RichODEPressureAuxVarGetPermeabilityXYZ
     procedure, public :: GetPorosity             => RichODEPressureAuxVarGetPorosity
     procedure, public :: GetDensityType          => RichODEPressureAuxVarGetDensityType
     procedure, public :: GetPotMassSinkPressure  => RichODEPressureAuxVarGetPotMassSinkPressure
     procedure, public :: GetPotMassSinkExponent  => RichODEPressureAuxVarGetPotMassSinkExponent
     procedure, public :: GetViscosity            => RichODEPressureAuxVarGetViscosity
     procedure, public :: GetRelativePermeability => RichODEPressureAuxVarGetRelativePermeability
     procedure, public :: GetLiquidSaturation     => RichODEPressureAuxVarGetLiquidSaturation
     procedure, public :: GetDensity              => RichODEPressureAuxVarGetDensity
     procedure, public :: GetDPorDP               => RichODEPressureAuxVarGetDPorDP
     procedure, public :: GetDVisDP               => RichODEPressureAuxVarGetDVisDP
     procedure, public :: GetDKrDP                => RichODEPressureAuxVarGetDKrDP
     procedure, public :: GetDSatDP               => RichODEPressureAuxVarGetDSatDP
     procedure, public :: GetDDenDP               => RichODEPressureAuxVarGetDDenDP
     procedure, public :: GetDVisDT               => RichODEPressureAuxVarGetDVisDT
     procedure, public :: GetDDenDT               => RichODEPressureAuxVarGetDDenDT
  end type rich_ode_pres_auxvar_type

  public :: RichODEPressureAuxVarInit
  public :: RichODEPressureAuxVarCopy

  !------------------------------------------------------------------------
contains


  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize an auxiliary variable
    !
    ! !USES:
    use EOSWaterMod         , only : DENSITY_CONSTANT
    use PorosityFunctionMod , only : PorosityFunctionInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type) :: this

    this%pressure                = 0.d0
    this%temperature             = 273.15d0 + 25.d0
    this%frac_liq_sat            = 1.d0
    this%condition_value         = 0.d0
    this%pressure_prev           = 3.5355d3
    this%perm(:)                 = 0.d0
    this%por                     = 0.d0
    this%pot_mass_sink_pressure  = 0.d0
    this%pot_mass_sink_exponent  = 0.d0

    this%satParams%sat_func_type = 0
    this%satParams%sat_res       = 0.d0
    this%satParams%alpha         = 0.d0
    this%satParams%bc_lambda     = 0.d0
    this%satParams%vg_m          = 0.d0
    this%satParams%vg_n          = 0.d0

    call PorosityFunctionInit(this%porParams)

    this%vis                     = 0.d0
    this%dvis_dP                 = 0.d0
    this%dvis_dT                 = 0.d0
    this%dpor_dP                 = 0.d0
    this%dkr_dP                  = 0.d0
    this%dsat_dP                 = 0.d0
    this%dden_dP                 = 0.d0
    this%dden_dT                 = 0.d0

    this%density_type            = DENSITY_CONSTANT

  end subroutine RichODEPressureAuxVarInit

  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarCopy(this, auxvar)
    !
    ! !DESCRIPTION:
    ! Copies content of `auxvar`
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type) :: this
    class(rich_ode_pres_auxvar_type) :: auxvar

    this%pressure                = auxvar%pressure
    this%temperature             = auxvar%temperature
    this%frac_liq_sat            = auxvar%frac_liq_sat
    this%condition_value         = auxvar%condition_value
    this%pressure_prev           = auxvar%pressure_prev
    this%perm(:)                 = auxvar%perm(:)
    this%por                     = auxvar%por
    this%pot_mass_sink_pressure  = auxvar%pot_mass_sink_pressure
    this%pot_mass_sink_exponent  = auxvar%pot_mass_sink_exponent

    this%satParams%sat_func_type = auxvar%satParams%sat_func_type
    this%satParams%sat_res       = auxvar%satParams%sat_res
    this%satParams%alpha         = auxvar%satParams%alpha
    this%satParams%bc_lambda     = auxvar%satParams%bc_lambda
    this%satParams%vg_m          = auxvar%satParams%vg_m
    this%satParams%vg_n          = auxvar%satParams%vg_n

    call this%porParams%Copy(auxvar%porParams)
    call this%satParams%Copy(auxvar%satParams)

    this%vis                     = auxvar%vis
    this%dvis_dP                 = auxvar%dvis_dP
    this%dpor_dP                 = auxvar%dpor_dP
    this%dkr_dP                  = auxvar%dkr_dP
    this%dsat_dP                 = auxvar%dsat_dP
    this%dden_dP                 = auxvar%dden_dP
    this%density_type            = auxvar%density_type

  end subroutine RichODEPressureAuxVarCopy

  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetValue(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_PRESSURE_PREV
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type), intent(inout) :: this
    PetscInt, intent(in)                            :: var_type
    PetscReal, intent(in)                           :: variable_value

    select case(var_type)
    case (VAR_PRESSURE)
       this%pressure        = variable_value
    case (VAR_TEMPERATURE)
       this%temperature     = variable_value
    case (VAR_PRESSURE_PREV)
       this%pressure_prev   = variable_value
    case (VAR_BC_SS_CONDITION)
       this%condition_value = variable_value
    case default
       write(iulog,*) 'In RichODEPressureAuxVarSetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichODEPressureAuxVarSetValue


  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarGetValue(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants, only : VAR_LIQ_SAT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type) :: this
    PetscInt, intent(in)             :: var_type
    PetscReal, intent(out)           :: variable_value

    select case(var_type)
    case (VAR_PRESSURE)
       variable_value = this%pressure
    case (VAR_TEMPERATURE)
       variable_value = this%temperature
    case (VAR_BC_SS_CONDITION)
       variable_value = this%condition_value
    case (VAR_LIQ_SAT)
       variable_value = this%sat
    case default
       write(iulog,*) 'In RichODEPressureAuxVarGetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichODEPressureAuxVarGetValue


  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    ! Computes various secondary quantities (sat, den, etc) based on
    ! the primary quantity (pressure).
    !
    ! !USES:
    use EOSWaterMod           , only : Density
    use EOSWaterMod           , only : Viscosity
    use PorosityFunctionMod   , only : PorosityFunctionComputation
    use SaturationFunction    , only : SatFunc_PressToSat
    use SaturationFunction    , only : SatFunc_PressToRelPerm
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type)   :: this

    ! Compute saturation
    call SatFunc_PressToSat(this%satParams , &
         this%pressure                     , &
         this%sat                          , &
         this%dsat_dP                        &
         )

    ! Compute relative permeability
    call SatFunc_PressToRelPerm(this%satParams , &
         this%pressure                         , &
         this%frac_liq_sat                     , &
         this%kr                               , &
         this%dkr_dP                             &
         )

    ! Compute density
    call Density(this%pressure                      , &
         this%temperature                           , &
         this%density_type                          , &
         this%den                                   , &
         this%dden_dP                               , &
         this%dden_dT                                 &
         )

    ! Compute viscosity
    call Viscosity(this%pressure                    , &
         this%temperature                           , &
         this%vis                                   , &
         this%dvis_dP                               , &
         this%dvis_dT                                 &
         )

    ! Compute porosity
    call PorosityFunctionComputation(this%porParams , &
         this%pressure                              , &
         this%por                                   , &
         this%dpor_dP                                 &
         )

  end subroutine RichODEPressureAuxVarCompute

  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPressure(this, val)
    !
    ! !DESCRIPTION:
    ! Set pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%pressure = val
 
  end subroutine RichODEPressureAuxVarSetPressure
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPressurePrev(this, val)
    !
    ! !DESCRIPTION:
    ! Set pressure at previous timestep
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%pressure_prev = val
 
  end subroutine RichODEPressureAuxVarSetPressurePrev
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetTemperature(this, val)
    !
    ! !DESCRIPTION:
    ! Set temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%temperature = val
 
  end subroutine RichODEPressureAuxVarSetTemperature
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetFracLiqSat(this, val)
    !
    ! !DESCRIPTION:
    ! Set fraction of liquid saturation
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%frac_liq_sat = val
 
  end subroutine RichODEPressureAuxVarSetFracLiqSat
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetConditionValue(this, val)
    !
    ! !DESCRIPTION:
    ! Set boundary or source-sink value
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%condition_value = val
 
  end subroutine RichODEPressureAuxVarSetConditionValue
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPermeability(this, val)
    !
    ! !DESCRIPTION:
    ! Set permeability in x, y, and z-dir to a homogeneous value
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%perm(:) = val
 
  end subroutine RichODEPressureAuxVarSetPermeability
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPermeabilityX(this, val)
    !
    ! !DESCRIPTION:
    ! Set permeability in x-dir
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%perm(1) = val
 
  end subroutine RichODEPressureAuxVarSetPermeabilityX
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPermeabilityY(this, val)
    !
    ! !DESCRIPTION:
    ! Set permeability in y-dir
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%perm(2) = val
 
  end subroutine RichODEPressureAuxVarSetPermeabilityY
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPermeabilityZ(this, val)
    !
    ! !DESCRIPTION:
    ! Set permeability in z-dir
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%perm(3) = val
 
  end subroutine RichODEPressureAuxVarSetPermeabilityz
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPermeabilityXYZ(this, val)
    !
    ! !DESCRIPTION:
    ! Set porosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val(3)
 
    this%perm(:) = val(:)
 
  end subroutine RichODEPressureAuxVarSetPermeabilityXYZ
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPorosity(this, val)
    !
    ! !DESCRIPTION:
    ! Set porosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%por = val
 
  end subroutine RichODEPressureAuxVarSetPorosity
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDensityType(this, val)
    !
    ! !DESCRIPTION:
    ! Set porosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscInt                         :: val
 
    this%density_type = val
 
  end subroutine RichODEPressureAuxVarSetDensityType
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPotMassSinkPressure(this, val)
    !
    ! !DESCRIPTION:
    ! Set pressure value for the parameterization downregulating potential mass sink
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%pot_mass_sink_pressure = val
 
  end subroutine RichODEPressureAuxVarSetPotMassSinkPressure
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetPotMassSinkExponent(this, val)
    !
    ! !DESCRIPTION:
    ! Set exponent value for the parameterization downregulating potential mass sink
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%pot_mass_sink_exponent = val
 
  end subroutine RichODEPressureAuxVarSetPotMassSinkExponent
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetViscosity(this, val)
    !
    ! !DESCRIPTION:
    ! Set viscosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%vis = val
 
  end subroutine RichODEPressureAuxVarSetViscosity
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetRelativePermeability(this, val)
    !
    ! !DESCRIPTION:
    ! Set relative permeability
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%kr = val
 
  end subroutine RichODEPressureAuxVarSetRelativePermeability
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetLiquidSaturation(this, val)
    !
    ! !DESCRIPTION:
    ! Set liquid saturation
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%sat = val
 
  end subroutine RichODEPressureAuxVarSetLiquidSaturation
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDensity(this, val)
    !
    ! !DESCRIPTION:
    ! Set density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%den = val
 
  end subroutine RichODEPressureAuxVarSetDensity
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDPorDP(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of porosity w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dpor_dP = val
 
  end subroutine RichODEPressureAuxVarSetDPorDP
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDVisDP(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of viscosity w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dvis_dP = val
 
  end subroutine RichODEPressureAuxVarSetDVisDP
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDKrDP(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of relative permeability w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dkr_dP = val
 
  end subroutine RichODEPressureAuxVarSetDKrDP
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDSatDP(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of saturation w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dsat_dP = val
 
  end subroutine RichODEPressureAuxVarSetDSatDP
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDDenDP(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of density w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dden_dP = val
 
  end subroutine RichODEPressureAuxVarSetDDenDP
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDVisDT(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of viscosity w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dvis_dT = val
 
  end subroutine RichODEPressureAuxVarSetDVisDT

  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetDDenDT(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of density w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dden_dT = val
 
  end subroutine RichODEPressureAuxVarSetDDenDT

  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPressure(this)
    !
    ! !DESCRIPTION:
    ! Get pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPressure
 
    RichODEPressureAuxVarGetPressure = this%pressure
 
  end function RichODEPressureAuxVarGetPressure
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPressurePrev(this)
    !
    ! !DESCRIPTION:
    ! Get pressure at previous timestep
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPressurePrev
 
    RichODEPressureAuxVarGetPressurePrev = this%pressure_prev
 
  end function RichODEPressureAuxVarGetPressurePrev
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetTemperature(this)
    !
    ! !DESCRIPTION:
    ! Get temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetTemperature
 
    RichODEPressureAuxVarGetTemperature = this%temperature
 
  end function RichODEPressureAuxVarGetTemperature
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetFracLiqSat(this)
    !
    ! !DESCRIPTION:
    ! Get fraction of liquid saturation
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetFracLiqSat
 
    RichODEPressureAuxVarGetFracLiqSat = this%frac_liq_sat
 
  end function RichODEPressureAuxVarGetFracLiqSat
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetConditionValue(this)
    !
    ! !DESCRIPTION:
    ! Get boundary or source-sink value
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetConditionValue
 
    RichODEPressureAuxVarGetConditionValue = this%condition_value
 
  end function RichODEPressureAuxVarGetConditionValue
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPermeabilityX(this)
    !
    ! !DESCRIPTION:
    ! Get  permeability in x-dir
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPermeabilityX
 
    RichODEPressureAuxVarGetPermeabilityX = this%perm(1)
 
  end function RichODEPressureAuxVarGetPermeabilityX
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPermeabilityY(this)
    !
    ! !DESCRIPTION:
    ! Get permeability in y-dir
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPermeabilityY
 
    RichODEPressureAuxVarGetPermeabilityY = this%perm(2)
 
  end function RichODEPressureAuxVarGetPermeabilityY
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPermeabilityZ(this)
    !
    ! !DESCRIPTION:
    ! Get permeability in z-dir
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPermeabilityZ
 
    RichODEPressureAuxVarGetPermeabilityZ = this%perm(3)
 
  end function RichODEPressureAuxVarGetPermeabilityZ
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPermeabilityXYZ(this)
    !
    ! !DESCRIPTION:
    ! Get permeability in x, y, and z direction
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPermeabilityXYZ(3)
 
    RichODEPressureAuxVarGetPermeabilityXYZ(:) = this%perm(:)
 
  end function RichODEPressureAuxVarGetPermeabilityXYZ
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPorosity(this)
    !
    ! !DESCRIPTION:
    ! Get porosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPorosity
 
    RichODEPressureAuxVarGetPorosity = this%por
 
  end function RichODEPressureAuxVarGetPorosity
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDensityType(this)
    !
    ! !DESCRIPTION:
    ! Get porosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscInt                         :: RichODEPressureAuxVarGetDensityType
 
    RichODEPressureAuxVarGetDensityType = this%density_type
 
  end function RichODEPressureAuxVarGetDensityType
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPotMassSinkPressure(this)
    !
    ! !DESCRIPTION:
    ! Get pressure value for the parameterization downregulating potential mass sink
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPotMassSinkPressure
 
    RichODEPressureAuxVarGetPotMassSinkPressure = this%pot_mass_sink_pressure
 
  end function RichODEPressureAuxVarGetPotMassSinkPressure
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetPotMassSinkExponent(this)
    !
    ! !DESCRIPTION:
    ! Get exponent value for the parameterization downregulating potential mass sink
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetPotMassSinkExponent
 
    RichODEPressureAuxVarGetPotMassSinkExponent = this%pot_mass_sink_exponent
 
  end function RichODEPressureAuxVarGetPotMassSinkExponent
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetViscosity(this)
    !
    ! !DESCRIPTION:
    ! Get viscosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetViscosity
 
    RichODEPressureAuxVarGetViscosity = this%vis
 
  end function RichODEPressureAuxVarGetViscosity
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetRelativePermeability(this)
    !
    ! !DESCRIPTION:
    ! Get relative permeability
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetRelativePermeability
 
    RichODEPressureAuxVarGetRelativePermeability = this%kr
 
  end function RichODEPressureAuxVarGetRelativePermeability
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetLiquidSaturation(this)
    !
    ! !DESCRIPTION:
    ! Get liquid saturation
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetLiquidSaturation
 
    RichODEPressureAuxVarGetLiquidSaturation = this%sat
 
  end function RichODEPressureAuxVarGetLiquidSaturation
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDensity(this)
    !
    ! !DESCRIPTION:
    ! Get density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetDensity
 
    RichODEPressureAuxVarGetDensity = this%den
 
  end function RichODEPressureAuxVarGetDensity
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDPorDP(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of porosity w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetDPorDP
 
    RichODEPressureAuxVarGetDPorDP = this%dpor_dP
 
  end function RichODEPressureAuxVarGetDPorDP
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDVisDP(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of viscosity w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetDVisDP
 
    RichODEPressureAuxVarGetDVisDP = this%dvis_dP
 
  end function RichODEPressureAuxVarGetDVisDP
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDKrDP(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of relative permeability w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetDKrDP
 
    RichODEPressureAuxVarGetDKrDP = this%dkr_dP
 
  end function RichODEPressureAuxVarGetDKrDP
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDSatDP(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of saturation w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetDSatDP
 
    RichODEPressureAuxVarGetDSatDP = this%dsat_dP
 
  end function RichODEPressureAuxVarGetDSatDP
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDDenDP(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of density w.r.t. pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetDDenDP
 
    RichODEPressureAuxVarGetDDenDP = this%dden_dP
 
  end function RichODEPressureAuxVarGetDDenDP
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDVisDT(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of viscosity w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetDVisDT
 
    RichODEPressureAuxVarGetDVisDT = this%dvis_dT
 
  end function RichODEPressureAuxVarGetDVisDT
 
  !------------------------------------------------------------------------
  function RichODEPressureAuxVarGetDDenDT(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of density w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureAuxVarGetDDenDT
 
    RichODEPressureAuxVarGetDDenDT = this%dden_dT
 
  end function RichODEPressureAuxVarGetDDenDT

#endif

end module RichardsODEPressureAuxType
