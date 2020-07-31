module ml_model_boundary_conditions
  !
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
  use ml_model_global_vars
  !
  implicit none
  !
#include <petsc/finclude/petsc.h>

  public :: read_boundary_conditions

contains

  !------------------------------------------------------------------------
  subroutine read_boundary_conditions(istep)
    !
    implicit none
    !
    PetscInt :: istep
    !
    PetscInt :: icair

    do icair = 1, ncair

       Iskyb_vis(icair) = 0.8d0
       Iskyd_vis(icair) = 0.2d0
       Iskyb_nir(icair) = 0.8d0
       Iskyd_nir(icair) = 0.2d0

       Irsky(icair)     = 400.d0

       Pref(icair)      = 98620.d0
       Uref(icair)      = 5.d0
       Tref(icair)      = 295.d0
       Rhref(icair)     = 80.d0

    end do

  end subroutine read_boundary_conditions

end module ml_model_boundary_conditions
