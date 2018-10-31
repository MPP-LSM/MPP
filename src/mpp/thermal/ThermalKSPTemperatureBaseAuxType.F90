module ThermalKSPTemperatureBaseAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !
  ! !PUBLIC TYPES:
  use petscsys

  implicit none

  private

  type, public :: therm_ksp_temp_base_auxvar_type

     private
     
     ! primary unknown independent variable
     PetscReal :: temperature     ! [K]

     PetscBool :: is_active       !

     ! independent variable available from:
     !  - another governing equation, or
     !  - another system-of-equation
     PetscReal :: frac            ! [-]
     PetscReal :: dhsdT           ! [W/m^2/K]
     PetscReal :: dist_up         ! [m]
     PetscReal :: dist_dn         ! [m]

     ! derived quantities = f(state_variables, parameters)
     PetscReal :: therm_cond      ! [W/(m K)]
     PetscReal :: heat_cap_pva    ! [J/(m3 K)]

     ! If the auxvar corresponds to boundary condition
     ! or source sink, the value is stored in this
     ! variable
     PetscReal :: condition_value ! Depends

   contains
     procedure, public :: BaseInit                  => ThermKSPTempBaseAuxVarInit

     procedure, public :: SetTemperature            => ThermKSPTempBaseAuxVarSetTemperature
     procedure, public :: SetActive                 => ThermKSPTempBaseAuxVarSetActive
     procedure, public :: SetInactive               => ThermKSPTempBaseAuxVarSetInActive
     procedure, public :: SetArealFraction          => ThermKSPTempBaseAuxVarSetArealFraction
     procedure, public :: SetDHsDT                  => ThermKSPTempBaseAuxVarSetDHsDT
     procedure, public :: SetDistanceUpwind         => ThermKSPTempBaseAuxVarSetDistanceUpwind
     procedure, public :: SetDistanceDownwind       => ThermKSPTempBaseAuxVarSetDistanceDownwind
     procedure, public :: SetThermalConductivity    => ThermKSPTempBaseAuxVarSetThermalConductivity
     procedure, public :: SetVolumetricHeatCapacity => ThermKSPTempBaseAuxVarSetVolumetricHeatCapacity
     procedure, public :: SetConditionValue         => ThermKSPTempBaseAuxVarSetConditionValue

     procedure, public :: GetTemperature            => ThermKSPTempBaseAuxVarGetTemperature
     procedure, public :: GetArealFraction          => ThermKSPTempBaseAuxVarGetArealFraction
     procedure, public :: GetDHsDT                  => ThermKSPTempBaseAuxVarGetDHsDT
     procedure, public :: GetDistanceUpwind         => ThermKSPTempBaseAuxVarGetDistanceUpwind
     procedure, public :: GetDistanceDownwind       => ThermKSPTempBaseAuxVarGetDistanceDownwind
     procedure, public :: GetThermalConductivity    => ThermKSPTempBaseAuxVarGetThermalConductivity
     procedure, public :: GetVolumetricHeatCapacity => ThermKSPTempBaseAuxVarGetVolumetricHeatCapacity
     procedure, public :: GetConditionValue         => ThermKSPTempBaseAuxVarGetConditionValue

     procedure, public :: IsActive                  => ThermKSPTempBaseAuxVarIsActive
  end type therm_ksp_temp_base_auxvar_type

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_base_auxvar_type)   :: this

    this%temperature     = 273.15d0
    this%is_active       = PETSC_FALSE

    this%dhsdT           = 0.d0
    this%frac            = 0.d0
    this%dist_up         = 0.d0
    this%dist_dn         = 0.d0

    this%therm_cond      = 0.d0
    this%heat_cap_pva    = 0.d0
    this%condition_value = 0.d0

  end subroutine ThermKSPTempBaseAuxVarInit

  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetTemperature(this, val)
    !
    ! !DESCRIPTION:
    ! Set temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    PetscReal                        :: val
 
    this%temperature = val
 
  end subroutine ThermKSPTempBaseAuxVarSetTemperature
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetActive(this)
    !
    ! !DESCRIPTION:
    ! Set temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
 
    this%is_active = PETSC_TRUE
 
  end subroutine ThermKSPTempBaseAuxVarSetActive
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetInactive(this)
    !
    ! !DESCRIPTION:
    ! Set temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
 
    this%is_active = PETSC_FALSE
 
  end subroutine ThermKSPTempBaseAuxVarSetInactive
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetArealFraction(this, val)
    !
    ! !DESCRIPTION:
    ! Set areal fraction of the ground occupied
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    PetscReal                        :: val
 
    this%frac = val
 
  end subroutine ThermKSPTempBaseAuxVarSetArealFraction
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetDHsDT(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of ground heat flux w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dhsdT = val
 
  end subroutine ThermKSPTempBaseAuxVarSetDHsDT
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetDistanceUpwind(this, val)
    !
    ! !DESCRIPTION:
    ! Set upwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dist_up = val
 
  end subroutine ThermKSPTempBaseAuxVarSetDistanceUpwind
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetDistanceDownwind(this, val)
    !
    ! !DESCRIPTION:
    ! Set downwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dist_dn = val
 
  end subroutine ThermKSPTempBaseAuxVarSetDistanceDownwind
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetThermalConductivity(this, val)
    !
    ! !DESCRIPTION:
    ! Set thermal conductivity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    PetscReal                        :: val
 
    this%therm_cond = val
 
  end subroutine ThermKSPTempBaseAuxVarSetThermalConductivity
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetVolumetricHeatCapacity(this, val)
    !
    ! !DESCRIPTION:
    ! Set volumetric heat capacity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    PetscReal                        :: val
 
    this%heat_cap_pva = val
 
  end subroutine ThermKSPTempBaseAuxVarSetVolumetricHeatCapacity
 
  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarSetConditionValue(this, val)
    !
    ! !DESCRIPTION:
    ! Set value of boundary or source-sink condition
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    PetscReal                        :: val
 
    this%condition_value = val
 
  end subroutine ThermKSPTempBaseAuxVarSetConditionValue

  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarGetTemperature(this)
    !
    ! !DESCRIPTION:
    ! Get temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempBaseAuxVarGetTemperature
 
    ThermKSPTempBaseAuxVarGetTemperature = this%temperature
 
  end function ThermKSPTempBaseAuxVarGetTemperature
 
  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarGetArealFraction(this)
    !
    ! !DESCRIPTION:
    ! Get areal fraction of the ground occupied
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempBaseAuxVarGetArealFraction
 
    ThermKSPTempBaseAuxVarGetArealFraction = this%frac
 
  end function ThermKSPTempBaseAuxVarGetArealFraction
 
  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarGetDHsDT(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of ground heat flux w.r.t. temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempBaseAuxVarGetDHsDT
 
    ThermKSPTempBaseAuxVarGetDHsDT = this%dhsdT
 
  end function ThermKSPTempBaseAuxVarGetDHsDT
 
  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarGetDistanceUpwind(this)
    !
    ! !DESCRIPTION:
    ! Get upwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempBaseAuxVarGetDistanceUpwind
 
    ThermKSPTempBaseAuxVarGetDistanceUpwind = this%dist_up
 
  end function ThermKSPTempBaseAuxVarGetDistanceUpwind
 
  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarGetDistanceDownwind(this)
    !
    ! !DESCRIPTION:
    ! Get downwind distance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempBaseAuxVarGetDistanceDownwind
 
    ThermKSPTempBaseAuxVarGetDistanceDownwind = this%dist_dn
 
  end function ThermKSPTempBaseAuxVarGetDistanceDownwind
 
  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarGetThermalConductivity(this)
    !
    ! !DESCRIPTION:
    ! Get thermal conductivity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempBaseAuxVarGetThermalConductivity
 
    ThermKSPTempBaseAuxVarGetThermalConductivity = this%therm_cond
 
  end function ThermKSPTempBaseAuxVarGetThermalConductivity
 
  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarGetVolumetricHeatCapacity(this)
    !
    ! !DESCRIPTION:
    ! Get volumetric heat capacity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempBaseAuxVarGetVolumetricHeatCapacity
 
    ThermKSPTempBaseAuxVarGetVolumetricHeatCapacity = this%heat_cap_pva
 
  end function ThermKSPTempBaseAuxVarGetVolumetricHeatCapacity
 
  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarGetConditionValue(this)
    !
    ! !DESCRIPTION:
    ! Get value of boundary or source-sink condition
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempBaseAuxVarGetConditionValue
 
    ThermKSPTempBaseAuxVarGetConditionValue = this%condition_value
 
  end function ThermKSPTempBaseAuxVarGetConditionValue

  !------------------------------------------------------------------------
  function ThermKSPTempBaseAuxVarIsActive(this)
    !
    ! !DESCRIPTION:
    ! Query active flag
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_base_auxvar_type) :: this
    !
    PetscBool :: ThermKSPTempBaseAuxVarIsActive
 
    ThermKSPTempBaseAuxVarIsActive = this%is_active
 
  end function ThermKSPTempBaseAuxVarIsActive

#endif

end module ThermalKSPTemperatureBaseAuxType

