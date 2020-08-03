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
    use ml_model_utils, only : set_value_in_condition
    !
    implicit none
    !
    PetscInt :: istep
    !
    PetscInt :: icair

    do icair = 1, ncair

       call set_value_in_condition(Iskyb_vis, icair, 0.8d0)
       call set_value_in_condition(Iskyd_vis, icair, 0.2d0)
       call set_value_in_condition(Iskyb_nir, icair, 0.8d0)
       call set_value_in_condition(Iskyd_nir, icair, 0.2d0)

       call set_value_in_condition(Irsky, icair     , 400.d0)

       call set_value_in_condition(Pref, icair      , 98620.d0)
       call set_value_in_condition(Uref, icair      , 5.d0)
       call set_value_in_condition(Tref, icair     , 295.d0)
       call set_value_in_condition(Rhref, icair    , 80.d0)

    end do

  end subroutine read_boundary_conditions

end module ml_model_boundary_conditions
