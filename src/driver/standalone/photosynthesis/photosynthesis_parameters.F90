module photosynthesis_parameters

  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbPhotosynthesis , only : mpp_photosynthesis_type
  use WaterVaporMod                  , only : SatVap
  use photosynthesis_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: set_parameters

contains

  !------------------------------------------------------------------------
  subroutine set_parameters(phtsyn_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use SystemOfEquationsPhotosynthesisType , only : sysofeqns_photosynthesis_type
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnPhotosynthesisType            , only : goveqn_photosynthesis_type
    use MultiPhysicsProbConstants           , only : VAR_PHOTOSYNTHETIC_PATHWAY_C3, VAR_PHOTOSYNTHETIC_PATHWAY_C4
    use MultiPhysicsProbConstants           , only : VAR_STOMATAL_CONDUCTANCE_MEDLYN, VAR_STOMATAL_CONDUCTANCE_BBERRY
    use MultiPhysicsProbConstants           , only : VAR_WUE
    use MultiPhysicsProbConstants           , only : TFRZ
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_photosynthesis_type)                         :: phtsyn_mpp

    class(sysofeqns_base_type)                , pointer   :: base_soe
    class(sysofeqns_photosynthesis_type)      , pointer   :: soe
    class(goveqn_base_type)                   , pointer   :: cur_goveq
    PetscReal                                 , parameter :: tau = 0.1d0
    PetscReal                                 , parameter :: rho = 0.1d0
    PetscReal                                             :: tair, esat_tair, desat_tair, relhum, vpd_tleaf
    PetscReal :: esat_25C, desat_25C, eair
    PetscReal :: esat_current, desat_current
    PetscInt                                              :: k , icair, icell

    call phtsyn_mpp%soe%SetPointerToIthGovEqn(PHTSYN_GE, cur_goveq)

    call SatVap (273.15d0 + 25.d0, esat_25C, desat_25C)

    tair = TFRZ + 25.d0
    relhum = 80.d0
    call SatVap (tair, esat_tair, desat_tair)
    eair = esat_tair * relhum /100.d0
    vpd_tleaf = esat_25C - eair
   
    select type(cur_goveq)
    class is (goveqn_photosynthesis_type)
       do icair = 1, ncair
          do k = 1, nz_cair+1
             icell = (icair-1)*(nz_cair+1) + k

             cur_goveq%aux_vars_in(icell)%colim  = 1
             cur_goveq%aux_vars_in(icell)%c3psn  = c3psn
             cur_goveq%aux_vars_in(icell)%gstype = gstype

             cur_goveq%aux_vars_in(icell)%cair  = 380.d0                       ! [mmol/mol]
             cur_goveq%aux_vars_in(icell)%o2ref = 0.209d0 * 1000.d0            ! [mmol/mol]
             cur_goveq%aux_vars_in(icell)%apar  = 2000.d0 * ( 1.d0 - rho - tau) ! [mmol photon/m2/s]
             
             cur_goveq%aux_vars_in(icell)%tleaf = TFRZ + 11.d0 + 0.25d0 * (k - 1)

             ! Values at tleaf = 11C and tair = 25C
             cur_goveq%aux_vars_in(icell)%gbv   = 2.224407920268566d0
             cur_goveq%aux_vars_in(icell)%gbc   = 1.637448199187622d0

             ! Values at tleaf = 25C and tair = 25C
             !cur_goveq%aux_vars_in(icell)%gbv   = 2.304789280798361d0
             !cur_goveq%aux_vars_in(icell)%gbc   = 1.694493818233533d0

             call SatVap (cur_goveq%aux_vars_in(icell)%tleaf, esat_current, desat_current)

             if (cur_goveq%aux_vars_in(icell)%gstype == VAR_STOMATAL_CONDUCTANCE_MEDLYN .or. &
                 cur_goveq%aux_vars_in(icell)%gstype == VAR_WUE) then
                cur_goveq%aux_vars_in(icell)%eair = esat_current - vpd_tleaf
             else
                cur_goveq%aux_vars_in(icell)%eair = esat_current * relhum/100.d0
             end if

             cur_goveq%aux_vars_in(icell)%btran = 1.d0
             cur_goveq%aux_vars_in(icell)%dpai  = 1.d0
          end do
       end do
    end select

  end subroutine set_parameters

end module photosynthesis_parameters
