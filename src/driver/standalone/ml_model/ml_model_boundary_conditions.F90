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
  public :: allocate_memory

contains

  !------------------------------------------------------------------------
  subroutine allocate_memory()
    !
    use ml_model_global_vars, only : condition_type
    use ml_model_utils      , only : allocate_memory_for_condition
    use ml_model_meshes     , only : nleaf
    !
    implicit none

    call allocate_memory_for_condition(Iskyb_vis , ncair)
    call allocate_memory_for_condition(Iskyd_vis , ncair)
    call allocate_memory_for_condition(Iskyb_nir , ncair)
    call allocate_memory_for_condition(Iskyd_nir , ncair)
    call allocate_memory_for_condition(Irsky     , ncair)

    call allocate_memory_for_condition(Tref  , ncair)
    call allocate_memory_for_condition(Qref  , ncair)
    call allocate_memory_for_condition(Pref  , ncair)
    call allocate_memory_for_condition(Rhref , ncair)

    call allocate_memory_for_condition(co2ref , ncair)
    call allocate_memory_for_condition(o2ref , ncair)

    call allocate_memory_for_condition(Uref  , ncair)

    call allocate_memory_for_condition(Albsoib_vis, ncair)
    call allocate_memory_for_condition(Albsoib_nir, ncair)
    call allocate_memory_for_condition(Albsoid_vis, ncair)
    call allocate_memory_for_condition(Albsoid_nir, ncair)

    call allocate_memory_for_condition(tg, ncair)
    call allocate_memory_for_condition(soil_t, ncair)

    call allocate_memory_for_condition(sza, ncair)

    call allocate_memory_for_condition(gbh , ncair*ntree*(ntop-nbot+1)*nleaf)
    call allocate_memory_for_condition(gbv , ncair*ntree*(ntop-nbot+1)*nleaf)
    call allocate_memory_for_condition(gbc , ncair*ntree*(ntop-nbot+1)*nleaf)

    call allocate_memory_for_condition(gs_sun, ncair*ntree*(ntop-nbot+1))
    call allocate_memory_for_condition(gs_shd, ncair*ntree*(ntop-nbot+1))

    call allocate_memory_for_condition(Tcan      , ncair)

    call allocate_memory_for_condition(Tair      , ncair*nz_cair)
    call allocate_memory_for_condition(eair      , ncair*nz_cair)

    call allocate_memory_for_condition(Tleaf_sun , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(Tleaf_shd , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(Labs_leaf_sun, ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(Labs_leaf_shd, ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(gs_shd    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(Ileaf_sun_vis    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(Ileaf_shd_vis    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(Ileaf_sun_nir    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(Ileaf_shd_nir    , ncair*(ntop-nbot+1) )
    
    call allocate_memory_for_condition(Labs_soil, ncair)
    call allocate_memory_for_condition(Isoil_vis, ncair)
    call allocate_memory_for_condition(Isoil_nir, ncair)

  end subroutine allocate_memory

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

      ! 1-2
      call set_value_in_condition(Iskyb_vis, icair, 58.634213093025487d0)
      call set_value_in_condition(Iskyb_nir, icair, 80.957228667870822d0)

      ! 3-4
       call set_value_in_condition(Iskyd_vis, icair, 56.857286906974515d0)
       call set_value_in_condition(Iskyd_nir, icair, 34.534271332129173d0)

      ! 5
       call set_value_in_condition(Irsky, icair, 329.34600000000000d0)

      ! 6-8
       call set_value_in_condition(Tref, icair, 295.93499389648440d0)
       call set_value_in_condition(Qref, icair, 9.4800086099820265d-3)
       call set_value_in_condition(Pref, icair, 98620.000000000000d0)

       ! 9-10
       call set_value_in_condition(co2ref, icair, 367.00000000000000d0)
       call set_value_in_condition(o2ref, icair, 209.00000000000000d0)

       ! 11
       call set_value_in_condition(Uref, icair, 5.1689999999999996d0)
!       call set_value_in_condition(Rhref, icair    , 80.d0)

       ! 12-15
       call set_value_in_condition(Albsoib_vis, icair, 0.13634140074253082d0)
       call set_value_in_condition(Albsoib_nir, icair, 0.20634140074253082d0)
       call set_value_in_condition(Albsoid_vis, icair, 0.13634140074253082d0)
       call set_value_in_condition(Albsoid_nir, icair, 0.20634140074253082d0)

       ! 16
       if (istep == 1) then
          call set_value_in_condition(tg, icair, 295.93499389648440d0)
       endif

       ! 17
       call set_value_in_condition(soil_t, icair, 294.84927368164062d0)

       ! 18
       call set_value_in_condition(sza, icair, 1.3473335314944674d0)

  end do

  end subroutine read_boundary_conditions

end module ml_model_boundary_conditions
