module lbl_parameters

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbLBL , only : mpp_lbl_type
  use lbl_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: set_parameters

contains

  !------------------------------------------------------------------------
  subroutine set_parameters(lbl_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsLBLType  , only : sysofeqns_lbl_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnLeafBoundaryLayer   , only : goveqn_leaf_bnd_lyr_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_lbl_type) :: lbl_mpp

    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_lbl_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    !
    PetscInt                             :: k, icair, icell

    call lbl_mpp%soe%SetPointerToIthGovEqn(LBL_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_leaf_bnd_lyr_type)
       do icair = 1, ncair
          do k = 1, nz_cair+1
             icell = (icair-1)*(nz_cair+1) + k
             cur_goveq%aux_vars_in(icell)%patm  = 101325.d0                       ! [Pa]
             cur_goveq%aux_vars_in(icell)%wind  = 5.d0                            ! [m/s]
             cur_goveq%aux_vars_in(icell)%tair  = 273.15d0 +  25.d0               ! [K]
             cur_goveq%aux_vars_in(icell)%tleaf = 273.15d0 + 11.d0 + (k-1)*0.25d0 ! [K]
             cur_goveq%aux_vars_in(icell)%dleaf = 0.05d0                          ! [m]
          end do
       end do
    end select

  end subroutine set_parameters

end module lbl_parameters
