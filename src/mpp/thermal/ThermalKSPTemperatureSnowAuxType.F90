module  ThermalKSPTemperatureSnowAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !
  ! !USES:
  use ThermalKSPTemperatureBaseAuxType
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none

  private

  type, public, extends(therm_ksp_temp_base_auxvar_type)  :: therm_ksp_temp_snow_auxvar_type

     private

     PetscReal :: liq_areal_den       ! [kg/m^2] = h2osoi_liq
     PetscReal :: ice_areal_den       ! [kg/m^2] = h2osoi_ice
     PetscInt  :: num_snow_layer      ! number of snow layer
     PetscReal :: dz                  ! [m]
     PetscReal :: tuning_factor       !

   contains
     procedure, public :: Init                  => ThermKSPTempSnowAuxVarInit
     procedure, public :: AuxVarCompute         => ThermKSPTempSnowAuxVarCompute

     procedure, public :: SetLiqArealDensity    => ThermalKSPTemperatureSnowAuxSetLiqArealDensity
     procedure, public :: SetIceArealDensity    => ThermalKSPTemperatureSnowAuxSetIceArealDensity
     procedure, public :: SetNumberOfSnowLayers => ThermalKSPTemperatureSnowAuxSetNumberOfSnowLayers
     procedure, public :: SetTuningFactor       => ThermalKSPTemperatureSnowAuxSetTuningFactor
     procedure, public :: SetDepth              => ThermalKSPTemperatureSnowAuxSetDepth

     procedure, public :: GetLiqArealDensity    => ThermalKSPTemperatureSnowAuxGetLiqArealDensity
     procedure, public :: GetIceArealDensity    => ThermalKSPTemperatureSnowAuxGetIceArealDensity
     procedure, public :: GetNumberOfSnowLayers => ThermalKSPTemperatureSnowAuxGetNumberOfSnowLayers
     procedure, public :: GetTuningFactor       => ThermalKSPTemperatureSnowAuxGetTuningFactor
     procedure, public :: GetDepth              => ThermalKSPTemperatureSnowAuxGetDepth

  end type therm_ksp_temp_snow_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSnowAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_snow_auxvar_type)   :: this

    call this%BaseInit()

    this%liq_areal_den            = 0.d0
    this%ice_areal_den            = 0.d0
    this%num_snow_layer           = 0
    this%dz = 0.d0
    this%tuning_factor      = 1.d0

  end subroutine ThermKSPTempSnowAuxVarInit

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSnowAuxVarCompute(this, dz, vol)
    !
    ! !DESCRIPTION:
    !
    use mpp_varcon      , only : cpice,  cpliq, tkair, tkice
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_snow_auxvar_type)   :: this
    PetscReal                                :: dz
    PetscReal                                :: vol
    !
    ! LOCAL VARIABLES
    PetscReal :: bw

    if (.not.this%IsActive()) then
       return
    else
       bw = (this%ice_areal_den + this%liq_areal_den)/(this%GetArealFraction() * dz)
       call this%SetThermalConductivity(tkair + (7.75d-5*bw + 1.105d-6*bw*bw)*(tkice-tkair))
       call this%SetVolumetricHeatCapacity( (cpliq*this%liq_areal_den + cpice*this%ice_areal_den)/dz)
    endif

  end subroutine ThermKSPTempSnowAuxVarCompute

  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSnowAuxSetLiqArealDensity(this, val)
    !
    ! !DESCRIPTION:
    ! Set liquid areal density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    PetscReal                        :: val
 
    this%liq_areal_den = val
 
  end subroutine ThermalKSPTemperatureSnowAuxSetLiqArealDensity
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSnowAuxSetIceArealDensity(this, val)
    !
    ! !DESCRIPTION:
    ! Set ice areal density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    PetscReal                        :: val
 
    this%ice_areal_den = val
 
  end subroutine ThermalKSPTemperatureSnowAuxSetIceArealDensity
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSnowAuxSetNumberOfSnowLayers(this, val)
    !
    ! !DESCRIPTION:
    ! Set of number snow layers
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    PetscInt                        :: val
 
    this%num_snow_layer = val
 
  end subroutine ThermalKSPTemperatureSnowAuxSetNumberOfSnowLayers
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSnowAuxSetTuningFactor(this, val)
    !
    ! !DESCRIPTION:
    ! Set tuning factor
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    PetscReal                        :: val
 
    this%tuning_factor = val
 
  end subroutine ThermalKSPTemperatureSnowAuxSetTuningFactor
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSnowAuxSetDepth(this, val)
    !
    ! !DESCRIPTION:
    ! Set depth
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dz = val
 
  end subroutine ThermalKSPTemperatureSnowAuxSetDepth
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSnowAuxGetLiqArealDensity(this)
    !
    ! !DESCRIPTION:
    ! Get liquid areal density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSnowAuxGetLiqArealDensity
 
    ThermalKSPTemperatureSnowAuxGetLiqArealDensity = this%liq_areal_den
 
  end function ThermalKSPTemperatureSnowAuxGetLiqArealDensity
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSnowAuxGetIceArealDensity(this)
    !
    ! !DESCRIPTION:
    ! Get ice areal density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSnowAuxGetIceArealDensity
 
    ThermalKSPTemperatureSnowAuxGetIceArealDensity = this%ice_areal_den
 
  end function ThermalKSPTemperatureSnowAuxGetIceArealDensity
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSnowAuxGetNumberOfSnowLayers(this)
    !
    ! !DESCRIPTION:
    ! Get snow water
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    !
    PetscInt                         :: ThermalKSPTemperatureSnowAuxGetNumberOfSnowLayers
 
    ThermalKSPTemperatureSnowAuxGetNumberOfSnowLayers = this%num_snow_layer
 
  end function ThermalKSPTemperatureSnowAuxGetNumberOfSnowLayers
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSnowAuxGetTuningFactor(this)
    !
    ! !DESCRIPTION:
    ! Get tuning factor
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSnowAuxGetTuningFactor
 
    ThermalKSPTemperatureSnowAuxGetTuningFactor = this%tuning_factor
 
  end function ThermalKSPTemperatureSnowAuxGetTuningFactor
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSnowAuxGetDepth(this)
    !
    ! !DESCRIPTION:
    ! Get depth
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_snow_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSnowAuxGetDepth
 
    ThermalKSPTemperatureSnowAuxGetDepth = this%dz
 
  end function ThermalKSPTemperatureSnowAuxGetDepth

#endif

end module ThermalKSPTemperatureSnowAuxType
