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
    use ml_model_utils      , only : allocate_memory_for_condition
    use ml_model_meshes     , only : nleaf
    !
    implicit none

    call allocate_memory_for_condition(bnd_cond%Iskyb_vis , ncair)
    call allocate_memory_for_condition(bnd_cond%Iskyd_vis , ncair)
    call allocate_memory_for_condition(bnd_cond%Iskyb_nir , ncair)
    call allocate_memory_for_condition(bnd_cond%Iskyd_nir , ncair)
    call allocate_memory_for_condition(bnd_cond%Irsky     , ncair)

    call allocate_memory_for_condition(bnd_cond%tref  , ncair)
    call allocate_memory_for_condition(bnd_cond%qref  , ncair)
    call allocate_memory_for_condition(bnd_cond%pref  , ncair)

    call allocate_memory_for_condition(bnd_cond%co2ref , ncair)
    call allocate_memory_for_condition(bnd_cond%o2ref , ncair)

    call allocate_memory_for_condition(bnd_cond%uref  , ncair)

    call allocate_memory_for_condition(bnd_cond%Albsoib_vis, ncair)
    call allocate_memory_for_condition(bnd_cond%Albsoib_nir, ncair)
    call allocate_memory_for_condition(bnd_cond%Albsoid_vis, ncair)
    call allocate_memory_for_condition(bnd_cond%Albsoid_nir, ncair)

    call allocate_memory_for_condition(bnd_cond%tg, ncair)
    call allocate_memory_for_condition(bnd_cond%soil_t, ncair)

    call allocate_memory_for_condition(bnd_cond%sza, ncair)

    call allocate_memory_for_condition(bnd_cond%rhg, ncair)
    call allocate_memory_for_condition(bnd_cond%soilres, ncair)

    call allocate_memory_for_condition(int_cond%gbh , ncair*ntree*(ntop-nbot+1)*nleaf)
    call allocate_memory_for_condition(int_cond%gbv , ncair*ntree*(ntop-nbot+1)*nleaf)
    call allocate_memory_for_condition(int_cond%gbc , ncair*ntree*(ntop-nbot+1)*nleaf)

    call allocate_memory_for_condition(int_cond%gs_sun, ncair*ntree*(ntop-nbot+1))
    call allocate_memory_for_condition(int_cond%gs_shd, ncair*ntree*(ntop-nbot+1))

    call allocate_memory_for_condition(int_cond%Tcan      , ncair)

    call allocate_memory_for_condition(int_cond%Tair      , ncair*nz_cair)
    call allocate_memory_for_condition(int_cond%Qair      , ncair*nz_cair)
    call allocate_memory_for_condition(int_cond%Wind      , ncair*nz_cair)

    call allocate_memory_for_condition(int_cond%Tleaf_sun , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Tleaf_shd , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Labs_leaf_sun, ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Labs_leaf_shd, ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%gs_shd    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Ileaf_sun_vis    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Ileaf_shd_vis    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Ileaf_sun_nir    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Ileaf_shd_nir    , ncair*(ntop-nbot+1) )
    
    call allocate_memory_for_condition(int_cond%Labs_soil, ncair)
    call allocate_memory_for_condition(int_cond%Isoil_vis, ncair)
    call allocate_memory_for_condition(int_cond%Isoil_nir, ncair)

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
      call set_value_in_condition(bnd_cond%Iskyb_vis, icair, 58.634213093025487d0)
      call set_value_in_condition(bnd_cond%Iskyb_nir, icair, 80.957228667870822d0)

      ! 3-4
       call set_value_in_condition(bnd_cond%Iskyd_vis, icair, 56.857286906974515d0)
       call set_value_in_condition(bnd_cond%Iskyd_nir, icair, 34.534271332129173d0)

      ! 5
       call set_value_in_condition(bnd_cond%Irsky, icair, 329.34600000000000d0)

      ! 6-9
       call set_value_in_condition(bnd_cond%tref, icair, 295.93499389648440d0)
       call set_value_in_condition(bnd_cond%qref, icair, 9.4800086099820265d-3)
       call set_value_in_condition(bnd_cond%pref, icair, 98620.000000000000d0)
       call set_value_in_condition(bnd_cond%uref, icair, 5.1689999999999996d0)

       ! 10-11
       call set_value_in_condition(bnd_cond%co2ref, icair, 367.00000000000000d0)
       call set_value_in_condition(bnd_cond%o2ref, icair, 209.00000000000000d0)

       ! 12-15
       call set_value_in_condition(bnd_cond%Albsoib_vis, icair, 0.13634140074253082d0)
       call set_value_in_condition(bnd_cond%Albsoib_nir, icair, 0.20634140074253082d0)
       call set_value_in_condition(bnd_cond%Albsoid_vis, icair, 0.13634140074253082d0)
       call set_value_in_condition(bnd_cond%Albsoid_nir, icair, 0.20634140074253082d0)

       ! 16
       if (istep == 1) then
          call set_value_in_condition(bnd_cond%tg, icair, 295.93499389648440d0)
       endif

       ! 17
       call set_value_in_condition(bnd_cond%soil_t, icair, 294.84927368164062d0)

       ! 18
       call set_value_in_condition(bnd_cond%sza, icair, 1.3473335314944674d0)

       ! 19-20
       call set_value_in_condition(bnd_cond%rhg    , icair, 0.9984057411945876d0)
       call set_value_in_condition(bnd_cond%soilres, icair, 3361.509423807650d0)
    end do

  end subroutine read_boundary_conditions

end module ml_model_boundary_conditions
