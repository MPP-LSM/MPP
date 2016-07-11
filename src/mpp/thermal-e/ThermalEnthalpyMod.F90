
module ThermalEnthalpyMod

#ifdef USE_PETSC_LIB
  
  ! !USES:
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ThermalEnthalpyFlux

  !------------------------------------------------------------------------------

contains

  subroutine ThermalEnthalpyFlux( &
       T_up,                      &
       h_up,                      &
       dh_up_dT_up,               &
       den_up,                    &
       dden_up_dT_up,             &
       therm_cond_up,             &
       T_dn,                      &
       h_dn,                      &
       dh_dn_dT_dn,               &
       den_dn,                    &
       dden_dn_dT_dn,             &
       therm_cond_dn,             &
       v_darcy,                   &
       dv_darcy_dT_up,            &
       dv_darcy_dT_dn,            &
       area,                      &
       dist_up,                   &
       dist_dn,                   &
       compute_deriv,             &
       internal_conn,             &
       cond_type,                 &
       flux,                      &
       dflux_dT_up,               &
       dflux_dT_dn                &
       )
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : COND_DIRICHLET
    use MultiPhysicsProbConstants, only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: T_up
    PetscReal, intent(in)  :: h_up
    PetscReal, intent(in)  :: dh_up_dT_up
    PetscReal, intent(in)  :: den_up
    PetscReal, intent(in)  :: dden_up_dT_up
    PetscReal, intent(in)  :: therm_cond_up
    PetscReal, intent(in)  :: T_dn
    PetscReal, intent(in)  :: h_dn
    PetscReal, intent(in)  :: dh_dn_dT_dn
    PetscReal, intent(in)  :: den_dn
    PetscReal, intent(in)  :: dden_dn_dT_dn
    PetscReal, intent(in)  :: therm_cond_dn
    PetscReal, intent(in)  :: v_darcy
    PetscReal, intent(in)  :: dv_darcy_dT_up
    PetscReal, intent(in)  :: dv_darcy_dT_dn
    PetscReal, intent(in)  :: area
    PetscReal, intent(in)  :: dist_up
    PetscReal, intent(in)  :: dist_dn
    PetscBool, intent(in)  :: compute_deriv
    PetscBool, intent(in)  :: internal_conn
    PetscInt,  intent(in)  :: cond_type
    PetscReal, intent(out) :: flux
    PetscReal, intent(out) :: dflux_dT_up
    PetscReal, intent(out) :: dflux_dT_dn
    !
    ! !LOCAL VARIABLES
    PetscReal :: upweight
    PetscReal :: den_ave
    PetscReal :: dden_ave_dT_up
    PetscReal :: dden_ave_dT_dn
    PetscReal :: therm_cond_ave_over_dist
    PetscReal :: h
    PetscReal :: dh_dT_up
    PetscReal :: dh_dT_dn

    if (internal_conn) then
       upweight = dist_up/(dist_up + dist_dn)
       therm_cond_ave_over_dist = (therm_cond_up * therm_cond_dn) / (dist_up*therm_cond_dn + dist_dn*therm_cond_up)
    else

       select case(cond_type)
       case (COND_DIRICHLET)
          upweight   = 0.d0
          therm_cond_ave_over_dist = therm_cond_dn/(dist_up + dist_dn)

       case (COND_DIRICHLET_FRM_OTR_GOVEQ)
          upweight   = dist_up/(dist_up + dist_dn)
          therm_cond_ave_over_dist = (therm_cond_up * therm_cond_dn) / (dist_up*therm_cond_dn + dist_dn*therm_cond_up)

       case default
          write(iulog,*)'Unknown cond_type Add additional code.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    endif

    den_ave = upweight*den_up + (1.d0 - upweight)*den_dn

    if (v_darcy < 0.d0) then
       h = h_up
    else
       h = h_dn
    end if


    flux = v_darcy * den_ave * h * area + &
         -therm_cond_ave_over_dist*(T_up - T_dn)*area

    if (compute_deriv) then

       dden_ave_dT_up = upweight         *dden_up_dT_up
       dden_ave_dT_dn = (1.d0 - upweight)*dden_dn_dT_dn

       if (v_darcy < 0.d0) then
          dh_dT_up = dh_up_dT_up
          dh_dT_dn = 0.d0
       else
          dh_dT_up = 0.d0
          dh_dT_dn = dh_dn_dT_dn
       endif


       dflux_dT_up = dv_darcy_dT_up * den_ave        * h        * area + &
                     v_darcy        * dden_ave_dT_up * h        * area + &
                     v_darcy        * den_ave        * dh_dT_up * area + &
                     -therm_cond_ave_over_dist * area

       dflux_dT_dn = dv_darcy_dT_dn * den_ave        * h        * area + &
                     v_darcy        * dden_ave_dT_dn * h        * area + &
                     v_darcy        * den_ave        * dh_dT_dn * area + &
                     +therm_cond_ave_over_dist * area
    endif

  end subroutine ThermalEnthalpyFlux

#endif

end module ThermalEnthalpyMod
