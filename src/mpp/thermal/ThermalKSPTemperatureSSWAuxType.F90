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
     PetscReal :: dz ! [m]
   contains
     procedure, public :: Init                  => ThermKSPTempSSWAuxVarInit
     procedure, public :: AuxVarCompute         => ThermKSPTempSSWAuxVarCompute
  end type therm_ksp_temp_ssw_auxvar_type

  PetscReal, private, parameter :: thin_sfclayer = 1.0d-6   ! Threshold for thin surface layer

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

    if (.not.this%is_active) then
       return
    else

       this%therm_cond    = tkwat
       !this%heat_cap_pva  = cpliq*denh2o
       if ( (dz*this%frac*1.d3 > thin_sfclayer) .and. (this%frac > thin_sfclayer) ) then
          this%heat_cap_pva  = max(thin_sfclayer, cpliq*denh2o)
       else
          this%heat_cap_pva = thin_sfclayer
       endif
    endif

  end subroutine ThermKSPTempSSWAuxVarCompute

#endif

end module ThermalKSPTemperatureSSWAuxType
