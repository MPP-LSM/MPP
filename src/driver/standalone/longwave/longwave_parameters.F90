module longwave_parameters

  use mpp_varctl               , only : iulog
  use mpp_abortutils           , only : endrun
  use mpp_shr_log_mod          , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbLongwave , only : mpp_longwave_type
  use WaterVaporMod            , only : SatVap
  use longwave_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: set_parameters

contains

  !------------------------------------------------------------------------
  subroutine set_parameters(longwave_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType     , only : sysofeqns_base_type
    use SystemOfEquationsLongwaveType , only : sysofeqns_longwave_type
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnLongwaveType            , only : goveqn_longwave_type
    use MultiPhysicsProbConstants     , only : TFRZ
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_longwave_type)             :: longwave_mpp
    !
    class(goveqn_base_type) , pointer   :: cur_goveq
    class(connection_set_type) , pointer  :: cur_conn_set
    class(condition_type)      , pointer  :: cur_cond
    PetscInt                            :: k, icell, ileaf, iconn, sum_conn
    PetscReal               , parameter :: emleaf = 0.98d0
    PetscReal               , parameter :: emgrnd = 1.00d0
    PetscReal               , parameter :: Irsky  = 400.d0
    PetscReal               , parameter :: td = 0.915d0

    call longwave_mpp%soe%SetPointerToIthGovEqn(LONGWAVE_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_longwave_type)
       do k = 1, nz_cair+1
          icell = k

          cur_goveq%aux_vars_in(icell)%trans    = td
          cur_goveq%aux_vars_in(icell)%leaf_rho = 1.d0 - emleaf
          cur_goveq%aux_vars_in(icell)%leaf_tau = 0.d0
          cur_goveq%aux_vars_in(icell)%leaf_emiss = emleaf

          if (k == 1) then
             cur_goveq%aux_vars_in(icell)%is_soil = PETSC_TRUE
             cur_goveq%aux_vars_in(icell)%soil_temperature = TFRZ + 20.d0             
             cur_goveq%aux_vars_in(icell)%soil_emiss = emgrnd             
          end if

          do ileaf = 1, cur_goveq%aux_vars_in(icell)%nleaf
             cur_goveq%aux_vars_in(icell)%leaf_temperature(ileaf) = TFRZ + 25.d0
             cur_goveq%aux_vars_in(icell)%leaf_fssh(ileaf)    = 1.0d0
             cur_goveq%aux_vars_in(icell)%leaf_dpai(ileaf)        = 0.1d0
          end do

       end do

       sum_conn = 0
       cur_cond => cur_goveq%boundary_conditions%first
       do
          if (.not.associated(cur_cond)) exit

          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             icell = cur_conn_set%conn(iconn)%GetIDDn()

             cur_goveq%aux_vars_bc(sum_conn)%Idn = Irsky

          enddo
          cur_cond => cur_cond%next
       enddo

    end select

  end subroutine set_parameters

end module longwave_parameters
