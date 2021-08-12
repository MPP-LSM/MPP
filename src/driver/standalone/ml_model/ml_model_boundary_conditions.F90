module ml_model_boundary_conditions
  !
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
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
    implicit none

    call allocate_memory_for_bc()
    call allocate_memory_for_ic()
    call allocate_memory_for_canopy_level_vars()
    call allocate_memory_for_vertical_level_vars()

  end subroutine allocate_memory

  !------------------------------------------------------------------------
  subroutine allocate_memory_for_bc()
    !
    use ml_model_utils      , only : allocate_memory_for_condition
    use ml_model_meshes     , only : nleaf
    use ml_model_global_vars, only : ncair, ntree, ntop, nbot, nz_cair, bnd_cond
    !
    implicit none
    !
    call allocate_memory_for_condition(bnd_cond%Iskyb_vis   , ncair)
    call allocate_memory_for_condition(bnd_cond%Iskyd_vis   , ncair)
    call allocate_memory_for_condition(bnd_cond%Iskyb_nir   , ncair)
    call allocate_memory_for_condition(bnd_cond%Iskyd_nir   , ncair)
    call allocate_memory_for_condition(bnd_cond%Irsky       , ncair)

    call allocate_memory_for_condition(bnd_cond%tref        , ncair)
    call allocate_memory_for_condition(bnd_cond%qref        , ncair)
    call allocate_memory_for_condition(bnd_cond%pref        , ncair)

    call allocate_memory_for_condition(bnd_cond%co2ref      , ncair)
    call allocate_memory_for_condition(bnd_cond%o2ref       , ncair)

    call allocate_memory_for_condition(bnd_cond%uref        , ncair)

    call allocate_memory_for_condition(bnd_cond%Albsoib_vis , ncair)
    call allocate_memory_for_condition(bnd_cond%Albsoib_nir , ncair)
    call allocate_memory_for_condition(bnd_cond%Albsoid_vis , ncair)
    call allocate_memory_for_condition(bnd_cond%Albsoid_nir , ncair)

    call allocate_memory_for_condition(bnd_cond%tg          , ncair)
    call allocate_memory_for_condition(bnd_cond%soil_t      , ncair)

    call allocate_memory_for_condition(bnd_cond%sza         , ncair)

    call allocate_memory_for_condition(bnd_cond%rhg         , ncair)
    call allocate_memory_for_condition(bnd_cond%soilres     , ncair)
    call allocate_memory_for_condition(bnd_cond%soil_tk     , ncair)
    call allocate_memory_for_condition(bnd_cond%h2osoi_vol  , ncair*10)

    call allocate_memory_for_condition(bnd_cond%pref_prev   , ncair)

    call allocate_memory_for_condition(bnd_cond%fwet        , ncair*ntree*(ntop-nbot+1)*nleaf)
    call allocate_memory_for_condition(bnd_cond%fdry        , ncair*ntree*(ntop-nbot+1)*nleaf)

  end subroutine allocate_memory_for_bc

  !------------------------------------------------------------------------
  subroutine allocate_memory_for_ic()
    !
    use ml_model_utils      , only : allocate_memory_for_condition
    use ml_model_meshes     , only : nleaf
    use ml_model_global_vars, only : ncair, ntree, ntop, nbot, nz_cair, int_cond
    !
    implicit none
    !
    call allocate_memory_for_condition(int_cond%gbh           , ncair*ntree*(ntop-nbot+1)*nleaf)
    call allocate_memory_for_condition(int_cond%gbv           , ncair*ntree*(ntop-nbot+1)*nleaf)
    call allocate_memory_for_condition(int_cond%gbc           , ncair*ntree*(ntop-nbot+1)*nleaf)

    call allocate_memory_for_condition(int_cond%gs_sun        , ncair*ntree*(ntop-nbot+1))
    call allocate_memory_for_condition(int_cond%gs_shd        , ncair*ntree*(ntop-nbot+1))

    call allocate_memory_for_condition(int_cond%Tcan          , ncair)

    call allocate_memory_for_condition(int_cond%Tair          , ncair*nz_cair)
    call allocate_memory_for_condition(int_cond%Qair          , ncair*nz_cair)
    call allocate_memory_for_condition(int_cond%Wind          , ncair*nz_cair)

    call allocate_memory_for_condition(int_cond%Tleaf_sun     , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Tleaf_shd     , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Labs_leaf_sun , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Labs_leaf_shd , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%gs_shd        , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Ileaf_sun_vis , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Ileaf_shd_vis , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Ileaf_sun_nir , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(int_cond%Ileaf_shd_nir , ncair*(ntop-nbot+1) )
    
    call allocate_memory_for_condition(int_cond%Labs_soil     , ncair)
    call allocate_memory_for_condition(int_cond%Isoil_vis     , ncair)
    call allocate_memory_for_condition(int_cond%Isoil_nir     , ncair)

  end subroutine allocate_memory_for_ic

  !------------------------------------------------------------------------
  subroutine allocate_memory_for_canopy_level_vars()
    !
    use ml_model_utils      , only : allocate_memory_for_condition
    use ml_model_global_vars, only : ncair, canp_lev_vars
    !
    implicit none
    !
    call allocate_memory_for_condition(canp_lev_vars%ustar     , ncair)
    call allocate_memory_for_condition(canp_lev_vars%lup       , ncair)

    call allocate_memory_for_condition(canp_lev_vars%labs_soi  , ncair)
    call allocate_memory_for_condition(canp_lev_vars%rnabs_soi , ncair)
    call allocate_memory_for_condition(canp_lev_vars%sh_soi    , ncair)
    call allocate_memory_for_condition(canp_lev_vars%lh_soi    , ncair)
    call allocate_memory_for_condition(canp_lev_vars%et_soi    , ncair)
    call allocate_memory_for_condition(canp_lev_vars%g_soi     , ncair)
    call allocate_memory_for_condition(canp_lev_vars%gac0_soi  , ncair)

  end subroutine allocate_memory_for_canopy_level_vars

  !------------------------------------------------------------------------
  subroutine allocate_memory_for_vertical_level_vars()
    !
    use ml_model_utils      , only : allocate_memory_for_condition
    use ml_model_global_vars, only : ncair, nbot, ntop, ntree, vert_lev_vars
    use ml_model_meshes     , only : nleaf
    !
    implicit none
    !
    PetscInt :: size

    size = ncair*(ntop-nbot+1)

    call allocate_memory_for_condition(vert_lev_vars%sh_air      , size)
    call allocate_memory_for_condition(vert_lev_vars%et_air      , size)
    call allocate_memory_for_condition(vert_lev_vars%st_air      , size)
    call allocate_memory_for_condition(vert_lev_vars%gac_air     , size)

    size = ncair*(ntop-nbot+1)*ntree*nleaf

    call allocate_memory_for_condition(vert_lev_vars%labs_leaf   , size)
    call allocate_memory_for_condition(vert_lev_vars%rn_leaf     , size)
    call allocate_memory_for_condition(vert_lev_vars%sh_leaf     , size)
    call allocate_memory_for_condition(vert_lev_vars%lh_leaf     , size)
    call allocate_memory_for_condition(vert_lev_vars%tr_leaf     , size)
    call allocate_memory_for_condition(vert_lev_vars%st_leaf     , size)
    call allocate_memory_for_condition(vert_lev_vars%anet_leaf   , size)
    call allocate_memory_for_condition(vert_lev_vars%agross_leaf , size)
    call allocate_memory_for_condition(vert_lev_vars%gs_leaf     , size)

  end subroutine allocate_memory_for_vertical_level_vars

  !------------------------------------------------------------------------
  subroutine read_boundary_conditions(istep, bc_data)
    !
    use ml_model_utils, only : set_value_in_condition, get_value_from_condition
    use petscsys
    use petscvec
    use ml_model_global_vars, only : bnd_cond, ncair
    !
    implicit none
    !
    PetscInt            :: istep
    Vec                 :: bc_data
    !
    PetscInt, parameter :: ncol = 31
    PetscInt            :: icair, offset, size, k, j
    PetscReal, pointer  :: bc_p(:)
    PetscReal           :: pref_prev
    PetscErrorCode      :: ierr
    PetscInt, parameter :: nlev = 10

    offset = (istep-1)*ncol

    call VecGetSize(bc_data, size, ierr); CHKERRQ(ierr)
    if (istep*ncol > size) then
       write(*,*)'ERROR: Time step exceeds the boundary condition dataset'
       call exit(0)
    end if

    call VecGetArrayF90(bc_data, bc_p, ierr); CHKERRQ(ierr)

    do icair = 1, ncair

       ! 1-2
       call set_value_in_condition(bnd_cond%Iskyb_vis   , icair, bc_p(offset +  1))
       call set_value_in_condition(bnd_cond%Iskyb_nir   , icair, bc_p(offset +  2))

       ! 3-4
       call set_value_in_condition(bnd_cond%Iskyd_vis   , icair, bc_p(offset +  3))
       call set_value_in_condition(bnd_cond%Iskyd_nir   , icair, bc_p(offset +  4))

       ! 5
       call set_value_in_condition(bnd_cond%Irsky       , icair, bc_p(offset +  5))

       ! 6-9
       call set_value_in_condition(bnd_cond%tref        , icair, bc_p(offset +  6))
       call set_value_in_condition(bnd_cond%qref        , icair, bc_p(offset +  7))
       call set_value_in_condition(bnd_cond%pref        , icair, bc_p(offset +  8))
       call set_value_in_condition(bnd_cond%uref        , icair, bc_p(offset +  9))

       ! 10-11
       call set_value_in_condition(bnd_cond%co2ref      , icair, bc_p(offset +  10))
       call set_value_in_condition(bnd_cond%o2ref       , icair, bc_p(offset +  11))

       ! 12-15
       call set_value_in_condition(bnd_cond%Albsoib_vis , icair, bc_p(offset +  12))
       call set_value_in_condition(bnd_cond%Albsoib_nir , icair, bc_p(offset +  13))
       call set_value_in_condition(bnd_cond%Albsoid_vis , icair, bc_p(offset +  14))
       call set_value_in_condition(bnd_cond%Albsoid_nir , icair, bc_p(offset +  15))

       ! 16
       if (istep == 1) then
          call set_value_in_condition(bnd_cond%tg       , icair, bc_p(offset +  16))
       endif

       ! 17
       call set_value_in_condition(bnd_cond%soil_t      , icair, bc_p(offset +  17))

       ! 18
       call set_value_in_condition(bnd_cond%sza         , icair, bc_p(offset +  18))

       ! 19-20
       call set_value_in_condition(bnd_cond%rhg         , icair, bc_p(offset +  19))
       call set_value_in_condition(bnd_cond%soilres     , icair, bc_p(offset +  20))

       ! 21
       call set_value_in_condition(bnd_cond%soil_tk     , icair, bc_p(offset +  21))

       ! 22-31
       do j = 1, nlev
          call set_value_in_condition(bnd_cond%h2osoi_vol, (icair-1)*nlev + j, bc_p(offset +  21 + j))
       end do
#if 0
       idx = offset + 21
       do k = 1, ncair*ntree*(ntop-nbot+1)*nleaf
          idx = idx + 1
          call set_value_in_condition(bnd_cond%fwet, k, bc_p(idx)
       end do

       do k = 1, ncair*ntree*(ntop-nbot+1)*nleaf
          idx = idx + 1
          call set_value_in_condition(bnd_cond%fdry, k, bc_p(idx)
       end do
#endif

       if (istep == 1) then
          pref_prev = get_value_from_condition(bnd_cond%pref, icair)
       else
          pref_prev = bc_p( (istep - 2)*ncol + 8)
       end if
       call set_value_in_condition(bnd_cond%pref_prev, icair, pref_prev)

    end do

    call VecRestoreArrayF90(bc_data, bc_p, ierr)

  end subroutine read_boundary_conditions

end module ml_model_boundary_conditions
