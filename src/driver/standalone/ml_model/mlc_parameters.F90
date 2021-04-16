module mlc_parameters

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use ml_model_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: mlc_set_parameters

contains

  !------------------------------------------------------------------------
  subroutine mlc_set_parameters(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp

    call set_air_temp_ge_parameters(mlc_mpp)
    call set_air_vapor_ge_parameters(mlc_mpp)
    call set_common_canopy_leaf_parameters(mlc_mpp, CLEF_TEMP_SUN_GE)
    call set_common_canopy_leaf_parameters(mlc_mpp, CLEF_TEMP_SHD_GE)
    call set_turbulence_parameters(mlc_mpp)
    call set_soil_parameters(mlc_mpp)

  end subroutine mlc_set_parameters

  !------------------------------------------------------------------------
  subroutine set_air_temp_ge_parameters(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType , only : goveqn_cair_temp_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, ileaf, icair, icell

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CAIR_TEMP_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_cair_temp_type)


       do icair = 1, ncair
          do k = 1, nz_cair
             icell = (icair-1)*(nz_cair+1) + k

             cur_goveq%aux_vars_in(icell)%leaf_dpai(:) = dpai(k)
             cur_goveq%aux_vars_in(icell)%leaf_fwet(:) = 0.d0
             cur_goveq%aux_vars_in(icell)%leaf_fdry(:) = 0.8218390792391702d0

             do ileaf = 1, ntree
                cur_goveq%aux_vars_in(icell)%leaf_fssh(        ileaf) = fssh(k)
                cur_goveq%aux_vars_in(icell)%leaf_fssh(ntree + ileaf) = 1.d0 - fssh(k)
             enddo
          end do
          cur_goveq%aux_vars_in((nz_cair+1)*(icair-1)+1)%is_soil = PETSC_TRUE
       enddo

    end select

  end subroutine set_air_temp_ge_parameters

  !------------------------------------------------------------------------
  subroutine set_air_vapor_ge_parameters(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyAirVaporType  , only : goveqn_cair_vapor_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, ileaf, icair, icell

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CAIR_VAPR_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_cair_vapor_type)

      do icair = 1, ncair
         do k = 1, nz_cair
            icell = (icair-1)*(nz_cair+1) + k
            cur_goveq%aux_vars_in(icell)%gbv  = 2.496430918408511d0

            cur_goveq%aux_vars_in(icell)%leaf_dpai(:) = dpai(k)
            cur_goveq%aux_vars_in(icell)%leaf_fwet(:) = 0.d0
            cur_goveq%aux_vars_in(icell)%leaf_fdry(:) = 0.8218390792391702d0

            do ileaf = 1, ntree
               cur_goveq%aux_vars_in(icell)%leaf_fssh(        ileaf) = fssh(k)
               cur_goveq%aux_vars_in(icell)%leaf_fssh(ntree + ileaf) = 1.d0 - fssh(k)
            end do
         end do
         cur_goveq%aux_vars_in((nz_cair+1)*(icair-1)+1)%is_soil = PETSC_TRUE
      end do
    end select

  end subroutine set_air_vapor_ge_parameters

  !------------------------------------------------------------------------
  subroutine set_common_canopy_leaf_parameters(mlc_mpp, ge_rank)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp
    PetscInt           :: ge_rank
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, i, num_int, icell, icair, itree

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(ge_rank, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_cleaf_temp_type)


       do icair = 1, ncair
          do itree = 1, ntree
             do k = 1, nz_cair + 1
                icell = (icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cair+1) + k
                cur_goveq%aux_vars_in(icell)%cp   = 744.5333333333334d0
                cur_goveq%aux_vars_in(icell)%fwet = 0.d0
                cur_goveq%aux_vars_in(icell)%fdry = 0.8218390792391702d0
                cur_goveq%aux_vars_in(icell)%dpai = dpai(k)
                if (cur_goveq%rank_in_soe_list == CLEF_TEMP_SUN_GE) then
                   cur_goveq%aux_vars_in(icell)%fssh = fssh(k)
                else
                   cur_goveq%aux_vars_in(icell)%fssh = 1.d0 - fssh(k)
                end if
             end do
          end do
       end do

    end select

  end subroutine set_common_canopy_leaf_parameters

  !------------------------------------------------------------------------
  subroutine set_turbulence_parameters(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsMLCType  , only : sysofeqns_mlc_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type)                   :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscInt                             :: p

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do p = 1, ncair

       soe%cturb%pai(p)   = 5.051612734794617d0

       soe%cturb%hc(p)    = 21.d0
       soe%cturb%zref(p)  = 46.d0

    end do

  end subroutine set_turbulence_parameters

  !------------------------------------------------------------------------
  subroutine set_soil_parameters(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsMLCType  , only : sysofeqns_mlc_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type)                   :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscInt                             :: icair

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do icair = 1, ncair

       soe%cturb%tksoi(icair)  = 1.261326601469150d0
       soe%cturb%dzsoi(icair)  = 7.1006354171935350d-003
       soe%cturb%ressoi(icair) = 3361.509423807650d0
       soe%cturb%rhgsoi(icair) = 0.9984057411945876d0

    end do

  end subroutine set_soil_parameters

end module mlc_parameters
