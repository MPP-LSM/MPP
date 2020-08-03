module shortwave_parameters

  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbShortwave , only : mpp_shortwave_type
  use WaterVaporMod             , only : SatVap
  use shortwave_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: set_parameters

contains

  !------------------------------------------------------------------------
  subroutine set_parameters(shortwave_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsShortwaveType , only : sysofeqns_shortwave_type
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnShortwaveType            , only : goveqn_shortwave_type
    use ConditionType                  , only : condition_type
    use ConnectionSetType              , only : connection_set_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type)               :: shortwave_mpp
    !
    class(goveqn_base_type)    , pointer   :: cur_goveq
    class(connection_set_type) , pointer   :: cur_conn_set
    class(condition_type)      , pointer   :: cur_cond
    PetscInt                               :: k , icell, ileaf, iconn, sum_conn

    PetscReal                              :: sumlai, cumlai

    PetscReal                  , parameter :: clumpfac  = 1.d0
    PetscReal                  , parameter :: lai       = 6.0d0
    PetscReal                  , parameter :: lai_inc   = 0.1d0
    PetscReal                  , parameter :: Iskyb_vis = 0.8d0
    PetscReal                  , parameter :: Iskyd_vis = 0.2d0
    PetscReal                  , parameter :: Iskyb_nir = 0.8d0
    PetscReal                  , parameter :: Iskyd_nir = 0.2d0
    PetscReal                  , parameter :: Kb        = 0.577350269189626d0
    PetscReal                  , parameter :: td        = 0.913235689378651d0

    call shortwave_mpp%soe%SetPointerToIthGovEqn(shortwave_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_shortwave_type)

       do k = 1, nz_cair+1
          icell = k

          if (k == 1) then
             cur_goveq%aux_vars_in(icell)%is_soil = PETSC_TRUE
             cur_goveq%aux_vars_in(icell)%soil_albedo_b(1) = 0.1d0 ! vis + direct
             cur_goveq%aux_vars_in(icell)%soil_albedo_b(2) = 0.2d0 ! nir + direct

             cur_goveq%aux_vars_in(icell)%soil_albedo_d(1) = 0.1d0 ! vis + diffuse
             cur_goveq%aux_vars_in(icell)%soil_albedo_d(2) = 0.2d0 ! nir + diffuse
          else
             cur_goveq%aux_vars_in(icell)%Iskyb(1) = Iskyb_vis
             cur_goveq%aux_vars_in(icell)%Iskyb(2) = Iskyb_nir
             cur_goveq%aux_vars_in(icell)%Iskyd(1) = Iskyd_vis
             cur_goveq%aux_vars_in(icell)%Iskyd(2) = Iskyd_nir

             cur_goveq%aux_vars_in(icell)%leaf_rho(1)   = 0.10d0
             cur_goveq%aux_vars_in(icell)%leaf_rho(2)   = 0.45d0
             cur_goveq%aux_vars_in(icell)%leaf_tau(1)   = 0.05d0
             cur_goveq%aux_vars_in(icell)%leaf_tau(2)   = 0.25d0
             cur_goveq%aux_vars_in(icell)%leaf_omega(1) = 0.55d0
             cur_goveq%aux_vars_in(icell)%leaf_omega(2) = 0.30d0

             sumlai = 6.d0 - (k-1)*lai_inc + lai_inc/2.d0
             cumlai = 6.d0 - (k-2)*lai_inc

             cur_goveq%aux_vars_in(icell)%leaf_dlai        = lai_inc
             cur_goveq%aux_vars_in(icell)%leaf_fraction(1) = clumpfac * exp(-Kb * sumlai * clumpfac)
             cur_goveq%aux_vars_in(icell)%leaf_fraction(2) = 1.d0 - clumpfac * exp(-Kb * sumlai * clumpfac) 

             cur_goveq%aux_vars_in(icell)%leaf_tb    = exp(-Kb * lai_inc * clumpfac)
             cur_goveq%aux_vars_in(icell)%leaf_tbcum = exp(-Kb * cumlai  * clumpfac)
             cur_goveq%aux_vars_in(icell)%leaf_td    = td
          end if
       end do

       sum_conn = 0
       cur_cond => cur_goveq%boundary_conditions%first
       do
          if (.not.associated(cur_cond)) exit

          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             icell = cur_conn_set%conn(iconn)%GetIDDn()

             cur_goveq%aux_vars_bc(sum_conn)%Iskyb(1) = Iskyb_vis
             cur_goveq%aux_vars_bc(sum_conn)%Iskyb(2) = Iskyb_nir
             cur_goveq%aux_vars_bc(sum_conn)%Iskyd(1) = Iskyd_vis
             cur_goveq%aux_vars_bc(sum_conn)%Iskyd(2) = Iskyd_nir

          enddo
          cur_cond => cur_cond%next
       enddo

    end select

  end subroutine set_parameters

end module shortwave_parameters
