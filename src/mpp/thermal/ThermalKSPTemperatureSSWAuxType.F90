module ThermalKSPTemperatureSSWAuxType

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

  type, public, extends(therm_ksp_temp_base_auxvar_type)  :: therm_ksp_temp_ssw_auxvar_type
     private

     PetscReal :: dz ! [m]
   contains
     procedure, public :: Init          => ThermKSPTempSSWAuxVarInit
     procedure, public :: AuxVarCompute => ThermKSPTempSSWAuxVarCompute
     procedure, public :: SetWaterDepth => ThermKSPTempSSWAuxSetWaterDepth
     procedure, public :: GetWaterDepth => ThermKSPTempSSWAuxGetWaterDepth
  end type therm_ksp_temp_ssw_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_ssw_auxvar_type)   :: this

    call this%BaseInit()

    this%dz = 0.d0
    
  end subroutine ThermKSPTempSSWAuxVarInit

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWAuxVarCompute(this, dz, vol)
    !
    ! !DESCRIPTION:
    !
    use mpp_varcon      , only : cpliq, denh2o, tkwat
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_ssw_auxvar_type)   :: this
    PetscReal                                :: dz
    PetscReal                                :: vol
    !
    ! LOCAL VARIABLES

    if (.not.this%IsActive()) then
       return
    else
       call this%SetThermalConductivity(tkwat)
       call this%SetVolumetricHeatCapacity(cpliq*denh2o)
    endif

  end subroutine ThermKSPTempSSWAuxVarCompute

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSSWAuxSetWaterDepth(this, val)
    !
    ! !DESCRIPTION:
    ! Set height of standing water
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_ssw_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dz = val
 
  end subroutine ThermKSPTempSSWAuxSetWaterDepth
 
  !------------------------------------------------------------------------
  function ThermKSPTempSSWAuxGetWaterDepth(this)
    !
    ! !DESCRIPTION:
    ! Get height of standing water
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_ssw_auxvar_type) :: this
    !
    PetscReal                        :: ThermKSPTempSSWAuxGetWaterDepth
 
    ThermKSPTempSSWAuxGetWaterDepth = this%dz
 
  end function ThermKSPTempSSWAuxGetWaterDepth

#endif

end module ThermalKSPTemperatureSSWAuxType
