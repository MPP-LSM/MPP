module CanopyAirTemperatureConnAuxType

#ifdef USE_PETSC_LIB

  ! !USES:
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include <petsc/finclude/petsc.h>

  type, public :: cair_temp_conn_auxvar_type

     PetscReal :: ga ! Aerodynamic conductance for scalars (mol/m2/s)

   contains

     procedure, public :: Init => CAirTempConnAuxVarInit

  end type cair_temp_conn_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine CAirTempConnAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(cair_temp_conn_auxvar_type) :: this

    this%ga = 0.d0

  end subroutine CAirTempConnAuxVarInit

#endif

end module CanopyAirTemperatureConnAuxType
