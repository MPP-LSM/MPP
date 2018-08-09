
module ThermalEnthalpyMod

#ifdef USE_PETSC_LIB
  
#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                 , only : iulog
  use mpp_abortutils             , only : endrun
  use mpp_shr_log_mod            , only : errMsg => shr_log_errMsg
  use ThermalEnthalpySoilAuxType , only : therm_enthalpy_soil_auxvar_type
  use ConnectionSetType          , only : connection_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ThermalEnthalpyFlux
  public :: ThermalEnthalpyFluxDerivativeWrtPressure

  !------------------------------------------------------------------------------

contains

  subroutine ThermalEnthalpyFlux( &
       aux_var_up,                &
       aux_var_dn,                &
       mflux,                     &
       dmflux_dT_up,              &
       dmflux_dT_dn,              &
       conn_up2dn,                &
       compute_deriv,             &
       internal_conn,             &
       cond_type,                 &
       eflux,                     &
       deflux_dT_up,              &
       deflux_dT_dn               &
       )
    !
    ! !DESCRIPTION:
    ! Computes energy flux and it's derivative with respect to upwind/downwind
    ! temperature.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : COND_DIRICHLET
    use MultiPhysicsProbConstants, only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    type (therm_enthalpy_soil_auxvar_type) , intent(in)  :: aux_var_up
    type (therm_enthalpy_soil_auxvar_type) , intent(in)  :: aux_var_dn
    PetscReal, intent(in)  :: mflux
    PetscReal, intent(in)  :: dmflux_dT_up
    PetscReal, intent(in)  :: dmflux_dT_dn
    type(connection_type)                  , intent(in)  :: conn_up2dn
    PetscBool, intent(in)  :: compute_deriv
    PetscBool, intent(in)  :: internal_conn
    PetscInt,  intent(in)  :: cond_type
    PetscReal, intent(out) :: eflux
    PetscReal, intent(out) :: deflux_dT_up
    PetscReal, intent(out) :: deflux_dT_dn
    !
    ! !LOCAL VARIABLES
    PetscReal :: T_up
    PetscReal :: h_up
    PetscReal :: dh_up_dT_up
    PetscReal :: den_up
    PetscReal :: dden_up_dT_up
    PetscReal :: therm_cond_up
    PetscReal :: T_dn
    PetscReal :: h_dn
    PetscReal :: dh_dn_dT_dn
    PetscReal :: den_dn
    PetscReal :: dden_dn_dT_dn
    PetscReal :: therm_cond_dn
    PetscReal :: area
    PetscReal :: dist_up
    PetscReal :: dist_dn
    PetscReal :: upweight
    PetscReal :: den_ave
    PetscReal :: dden_ave_dT_up
    PetscReal :: dden_ave_dT_dn
    PetscReal :: therm_cond_ave_over_dist
    PetscReal :: h
    PetscReal :: dh_dT_up
    PetscReal :: dh_dT_dn

    T_up                 = aux_var_up%temperature
    h_up                 = aux_var_up%hl
    dh_up_dT_up          = aux_var_up%dhl_dT
    den_up               = aux_var_up%den
    dden_up_dT_up        = aux_var_up%dden_dT
    therm_cond_up        = aux_var_up%therm_cond

    T_dn                 = aux_var_dn%temperature
    h_dn                 = aux_var_dn%hl
    dh_dn_dT_dn          = aux_var_dn%dhl_dT
    den_dn               = aux_var_dn%den
    dden_dn_dT_dn        = aux_var_dn%dden_dT
    therm_cond_dn        = aux_var_dn%therm_cond

    area                 = conn_up2dn%GetArea()
    dist_up              = conn_up2dn%GetDistUp()
    dist_dn              = conn_up2dn%GetDistDn()

    if (internal_conn) then
       upweight = dist_up/(dist_up + dist_dn)
       therm_cond_ave_over_dist = (therm_cond_up * therm_cond_dn) / &
            (dist_up*therm_cond_dn + dist_dn*therm_cond_up)
    else

       select case(cond_type)
       case (COND_DIRICHLET)
          upweight   = 0.d0
          therm_cond_ave_over_dist = therm_cond_dn/(dist_up + dist_dn)

       case (COND_DIRICHLET_FRM_OTR_GOVEQ)
          upweight   = dist_up/(dist_up + dist_dn)
          therm_cond_ave_over_dist = (therm_cond_up * therm_cond_dn) / &
               (dist_up*therm_cond_dn + dist_dn*therm_cond_up)

       case default
          write(iulog,*)'Unknown cond_type Add additional code.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    endif

    den_ave = upweight*den_up + (1.d0 - upweight)*den_dn

    if (mflux <= 0.d0) then
       h = h_up
    else
       h = h_dn
    end if

    eflux = mflux * h + &
         (-therm_cond_ave_over_dist*(T_up - T_dn)*area)

    if (compute_deriv) then

       dden_ave_dT_up = upweight         *dden_up_dT_up
       dden_ave_dT_dn = (1.d0 - upweight)*dden_dn_dT_dn

       if (mflux < 0.d0) then
          dh_dT_up = dh_up_dT_up
          dh_dT_dn = 0.d0
       else
          dh_dT_up = 0.d0
          dh_dT_dn = dh_dn_dT_dn
       endif


       deflux_dT_up = dmflux_dT_up * h         + &
                      mflux        * dh_dT_up  + &
                      (-therm_cond_ave_over_dist * area)

       deflux_dT_dn = dmflux_dT_dn * h        + &
                      mflux        * dh_dT_dn + &
                      (+therm_cond_ave_over_dist * area)
    endif

  end subroutine ThermalEnthalpyFlux

  !------------------------------------------------------------------------------
  subroutine ThermalEnthalpyFluxDerivativeWrtPressure( &
       aux_var_up,                                     &
       aux_var_dn,                                     &
       mflux,                                          &
       dmflux_dP_up,                                   &
       dmflux_dP_dn,                                   &
       conn_up2dn,                                     &
       compute_deriv,                                  &
       internal_conn,                                  &
       cond_type,                                      &
       eflux,                                          &
       deflux_dP_up,                                   &
       deflux_dP_dn                                    &
       )
    !
    ! !DESCRIPTION:
    ! Computes derivative of energy flux with respect to upwind/downwind
    ! pressure.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : COND_DIRICHLET
    use MultiPhysicsProbConstants, only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    type (therm_enthalpy_soil_auxvar_type) , intent(in)  :: aux_var_up
    type (therm_enthalpy_soil_auxvar_type) , intent(in)  :: aux_var_dn
    PetscReal                              , intent(in)  :: mflux
    PetscReal                              , intent(in)  :: dmflux_dP_up
    PetscReal                              , intent(in)  :: dmflux_dP_dn
    type(connection_type)                  , intent(in)  :: conn_up2dn
    PetscBool                              , intent(in)  :: compute_deriv
    PetscBool                              , intent(in)  :: internal_conn
    PetscInt                               , intent(in)  :: cond_type
    PetscReal                              , intent(out) :: eflux
    PetscReal                              , intent(out) :: deflux_dP_up
    PetscReal                              , intent(out) :: deflux_dP_dn
    !
    ! !LOCAL VARIABLES
    PetscReal :: T_up
    PetscReal :: h_up
    PetscReal :: dh_up_dP_up
    PetscReal :: den_up
    PetscReal :: dden_up_dP_up
    PetscReal :: therm_cond_up
    PetscReal :: dtherm_cond_up_dP_up
    PetscReal :: T_dn
    PetscReal :: h_dn
    PetscReal :: dh_dn_dP_dn
    PetscReal :: den_dn
    PetscReal :: dden_dn_dP_dn
    PetscReal :: therm_cond_dn
    PetscReal :: dtherm_cond_dn_dP_dn
    PetscReal :: area
    PetscReal :: dist_up
    PetscReal :: dist_dn
    PetscReal :: upweight
    PetscReal :: den_ave
    PetscReal :: dden_ave_dP_up
    PetscReal :: dden_ave_dP_dn
    PetscReal :: dDk_dp_up
    PetscReal :: dDk_dp_dn
    PetscReal :: therm_cond_ave_over_dist
    PetscReal :: h
    PetscReal :: dh_dP_up
    PetscReal :: dh_dP_dn

    T_up                 = aux_var_up%temperature
    h_up                 = aux_var_up%hl
    dh_up_dP_up          = aux_var_up%dhl_dP
    den_up               = aux_var_up%den
    dden_up_dP_up        = aux_var_up%dden_dP
    therm_cond_up        = aux_var_up%therm_cond
    dtherm_cond_up_dP_up = aux_var_up%dtherm_cond_dP

    T_dn                 = aux_var_dn%temperature
    h_dn                 = aux_var_dn%hl
    dh_dn_dP_dn          = aux_var_dn%dhl_dP
    den_dn               = aux_var_dn%den
    dden_dn_dP_dn        = aux_var_dn%dden_dP
    therm_cond_dn        = aux_var_dn%therm_cond
    dtherm_cond_dn_dP_dn = aux_var_dn%dtherm_cond_dP

    area                 = conn_up2dn%GetArea()
    dist_up              = conn_up2dn%GetDistUp()
    dist_dn              = conn_up2dn%GetDistDn()

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

    if (mflux <= 0.d0) then
       h = h_up
    else
       h = h_dn
    end if


    eflux = mflux * h + &
         (-therm_cond_ave_over_dist*(T_up - T_dn)*area)

    if (compute_deriv) then

       dden_ave_dP_up = upweight         *dden_up_dP_up
       dden_ave_dP_dn = (1.d0 - upweight)*dden_dn_dP_dn

       if (mflux < 0.d0) then
          dh_dP_up = dh_up_dP_up
          dh_dP_dn = 0.d0
       else
          dh_dP_up = 0.d0
          dh_dP_dn = dh_dn_dP_dn
       endif

       if (internal_conn) then
          dDk_dp_up = therm_cond_ave_over_dist**2.d0/therm_cond_up**2.d0*dist_up*dtherm_cond_up_dP_up
          dDk_dp_dn = therm_cond_ave_over_dist**2.d0/therm_cond_dn**2.d0*dist_dn*dtherm_cond_dn_dP_dn
       else
          select case(cond_type)
          case (COND_DIRICHLET)
             upweight   = 0.d0
             therm_cond_ave_over_dist = therm_cond_dn/(dist_up + dist_dn)
             dDk_dp_up = 0.d0
             dDk_dp_dn = 1.d0/(dist_up + dist_dn)*dtherm_cond_dn_dP_dn

          case (COND_DIRICHLET_FRM_OTR_GOVEQ)
             dDk_dp_up = therm_cond_ave_over_dist**2.d0/therm_cond_up**2.d0*dist_up*dtherm_cond_up_dP_up
             dDk_dp_dn = therm_cond_ave_over_dist**2.d0/therm_cond_dn**2.d0*dist_dn*dtherm_cond_dn_dP_dn

          case default
             write(iulog,*)'Unknown cond_type Add additional code.'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       endif

       deflux_dP_up = dmflux_dP_up * h        + &
                      mflux        * dh_dP_up + &
                     (-dDk_dp_up*(T_up - T_dn)*area)


       deflux_dP_dn = dmflux_dP_dn * h        + &
                      mflux        * dh_dP_dn + &
                     (-dDk_dp_dn*(T_up - T_dn)*area)
    endif

  end subroutine ThermalEnthalpyFluxDerivativeWrtPressure

#endif

end module ThermalEnthalpyMod
