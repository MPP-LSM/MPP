module RichardsODEPressureConnAuxType

#ifdef USE_PETSC_LIB

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use SaturationFunction  , only : saturation_params_type
  !
  implicit none
  private

#include <petsc/finclude/petsc.h>

  type, public :: rich_ode_pres_conn_auxvar_type

     PetscReal, private :: pressure_up               ! [Pa]
     PetscReal, private :: pressure_dn               ! [Pa]
     PetscReal, private :: kr                        ! [-] used in CONDUCTANCE_CAMPBELL_TYPE

     PetscInt, private  :: flux_type
     PetscInt, private  :: conductance_type

     PetscReal, private :: conductance               ! used in CONDUCTANCE_CAMPBELL_TYPE
     PetscReal, private :: conductance_up            ! used in CONDUCTANCE_MANOLI_TYPE
     PetscReal, private :: conductance_dn            ! used in CONDUCTANCE_MANOLI_TYPE
     PetscReal, private :: dkr_dP_up                 ! used in CONDUCTANCE_CAMPBELL_TYPE
     PetscReal, private :: dkr_dP_dn                 ! used in CONDUCTANCE_CAMPBELL_TYPE
     PetscReal, private :: krg
     PetscReal, private :: dkrg_dP_up
     PetscReal, private :: dkrg_dP_dn
     PetscReal, private :: conductance_upwind_weight ! used in CONDUCTANCE_CAMPBELL_TYPE

     type(saturation_params_type) :: satParams_up
     type(saturation_params_type) :: satParams_dn

   contains
     procedure, public :: Init                       => RichODEPressureConnAuxInit
     procedure, public :: SetRealValue               => RichODEPressureConnAuxSetRealValue
     procedure, public :: SetIntValue                => RichODEPressureConnAuxSetIntValue
     procedure, public :: AuxVarCompute              => RichODEPressureConnAuxVarCompute

     procedure, public :: SetUpwindPressure          => RichODEPressureConnAuxSetUpwindPressure
     procedure, public :: SetDownwindPressure        => RichODEPressureConnAuxSetDownwindPressure
     procedure, public :: SetRelativePermeability    => RichODEPressureConnAuxSetRelativePermeability
     procedure, public :: SetFluxType                => RichODEPressureConnAuxSetFluxType
     procedure, public :: SetConductanceType         => RichODEPressureConnAuxSetConductanceType
     procedure, public :: SetConductance             => RichODEPressureConnAuxSetConductance
     procedure, public :: SetUpwindConductance       => RichODEPressureConnAuxSetUpwindConductance
     procedure, public :: SetDownwindConductance     => RichODEPressureConnAuxSetDownwindConductance
     procedure, public :: SetDKrDPup                 => RichODEPressureConnAuxSetDKrDPup
     procedure, public :: SetDKrDPdn                 => RichODEPressureConnAuxSetDKrDPdn
     procedure, public :: SetKrg                     => RichODEPressureConnAuxSetKrg
     procedure, public :: SetDKrgDPup                => RichODEPressureConnAuxSetDKrgDPup
     procedure, public :: SetDKrgDPdn                => RichODEPressureConnAuxSetDKrgDPdn
     procedure, public :: SetUpwindConductanceWeight => RichODEPressureConnAuxSetUpwindConductanceWeight

     procedure, public :: GetUpwindPressure          => RichODEPressureConnAuxGetUpwindPressure
     procedure, public :: GetDownwindPressure        => RichODEPressureConnAuxGetDownwindPressure
     procedure, public :: GetRelativePermeability    => RichODEPressureConnAuxGetRelativePermeability
     procedure, public :: GetFluxType                => RichODEPressureConnAuxGetFluxType
     procedure, public :: GetConductanceType         => RichODEPressureConnAuxGetConductanceType
     procedure, public :: GetConductance             => RichODEPressureConnAuxGetConductance
     procedure, public :: GetUpwindConductance       => RichODEPressureConnAuxGetUpwindConductance
     procedure, public :: GetDownwindConductance     => RichODEPressureConnAuxGetDownwindConductance
     procedure, public :: GetDKrDPup                 => RichODEPressureConnAuxGetDKrDPup
     procedure, public :: GetDKrDPdn                 => RichODEPressureConnAuxGetDKrDPdn
     procedure, public :: GetKrg                     => RichODEPressureConnAuxGetKrg
     procedure, public :: GetDKrgDPup                => RichODEPressureConnAuxGetDKrgDPup
     procedure, public :: GetDKrgDPdn                => RichODEPressureConnAuxGetDKrgDPdn
     procedure, public :: GetUpwindConductanceWeight => RichODEPressureConnAuxGetUpwindConductanceWeight

  end type rich_ode_pres_conn_auxvar_type

contains
  
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize the object
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_CAMPBELL_TYPE
    use petscsys
    !
    implicit none
    !
    class(rich_ode_pres_conn_auxvar_type) :: this

    this%pressure_up               = 0.d0
    this%pressure_dn               = 0.d0
    this%conductance               = 0.d0
    this%conductance_up            = 0.d0
    this%conductance_dn            = 0.d0

    this%flux_type                 = DARCY_FLUX_TYPE
    this%conductance_type          = CONDUCTANCE_CAMPBELL_TYPE
    this%dkr_dP_up                 = 0.d0
    this%dkr_dP_dn                 = 0.d0
    this%kr                        = 0.d0
    this%krg                       = 0.d0
    this%dkrg_dP_up                = 0.d0
    this%dkrg_dP_dn                = 0.d0
    this%conductance_upwind_weight = 0.d0

    call this%satParams_up%Init()
    call this%satParams_dn%Init()
    
  end subroutine RichODEPressureConnAuxInit

  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetRealValue(this, var_type, variable_value )
    !
    ! !DESCRIPTION:
    ! Set a real value
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_CONDUCTANCE
    use MultiPhysicsProbConstants, only : VAR_PRESSURE_UP
    use MultiPhysicsProbConstants, only : VAR_PRESSURE_DN
    use MultiPhysicsProbConstants, only : VAR_CONDUCTANCE_UP
    use MultiPhysicsProbConstants, only : VAR_CONDUCTANCE_DN
    use MultiPhysicsProbConstants, only : VAR_CAMPBELL_HE
    use MultiPhysicsProbConstants, only : VAR_CAMPBELL_N
    !
    implicit none
    !
    class(rich_ode_pres_conn_auxvar_type) , intent(inout) :: this
    PetscInt                              , intent(in)    :: var_type
    PetscReal                             , intent(in)    :: variable_value

    select case(var_type)
    case (VAR_CONDUCTANCE)
       this%conductance = variable_value

    case (VAR_PRESSURE_UP)
       this%pressure_up = variable_value

    case (VAR_PRESSURE_DN)
       this%pressure_dn = variable_value

    case (VAR_CONDUCTANCE_UP)
       this%conductance_up = variable_value

    case (VAR_CONDUCTANCE_DN)
       this%conductance_dn = variable_value

    case default
       write(iulog,*) 'In RichODEPressureConnAuxSetRealValue: unknown type ',var_type
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select
    
     end subroutine RichODEPressureConnAuxSetRealValue

  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetIntValue(this, var_type, variable_value )
    !
    ! !DESCRIPTION:
    ! Set an integer value
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_FLUX_TYPE
    use MultiPhysicsProbConstants, only : VAR_CONDUCTANCE_TYPE
    use MultiPhysicsProbConstants, only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_CAMPBELL_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_MANOLI_TYPE
    !
    implicit none
    !
    class(rich_ode_pres_conn_auxvar_type) , intent(inout) :: this
    PetscInt                              , intent(in)    :: var_type
    PetscInt                              , intent(in)    :: variable_value

    select case(var_type)
    case (VAR_FLUX_TYPE)
       select case(variable_value)
       case (DARCY_FLUX_TYPE)
          this%flux_type = variable_value
          
       case (CONDUCTANCE_FLUX_TYPE)
          this%flux_type = variable_value
          
       case default
          write(iulog,*)'Unknown value for VAR_FLUX_TYPE = ',variable_value
          call endrun(msg=errMsg(__FILE__,__LINE__))
       end select

    case (VAR_CONDUCTANCE_TYPE)
       select case(variable_value)
       case (CONDUCTANCE_CAMPBELL_TYPE)
          this%conductance_type = CONDUCTANCE_CAMPBELL_TYPE

       case (CONDUCTANCE_MANOLI_TYPE)
          this%conductance_type = CONDUCTANCE_MANOLI_TYPE

       case (0)
          ! Do nothing

       case default
          write(iulog,*)'Unknown value for VAR_CONDUCTANCE_TYPE = ',variable_value
          call endrun(msg=errMsg(__FILE__,__LINE__))
       end select

    case default
       write(iulog,*) 'In RichODEPressureConnAuxSetIntValue: unknown type ',var_type
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine RichODEPressureConnAuxSetIntValue


  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    ! Computes various secondary quantities (sat, den, etc) based on
    ! the primary quantity (pressure).
    !
    ! !USES:
    use SaturationFunction       , only : SatFunc_PressToRelPerm
    use MultiPhysicsProbConstants, only : CONDUCTANCE_CAMPBELL_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_MANOLI_TYPE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                             :: frac_liq_sat
    PetscReal                             :: kr_up, kr_dn
    PetscReal                             :: krg_up, krg_dn
    PetscReal                             :: denom
    PetscReal                             :: dkr_up_dP_up, dkr_dn_dP_dn

    select case(this%conductance_type)
    case(CONDUCTANCE_CAMPBELL_TYPE)
       if (this%satParams_up%relperm_func_type == 0 .and. &
            this%satParams_dn%relperm_func_type == 0  ) then

          this%kr         = 1.d0
          this%dkr_dP_up  = 0.d0
          this%dkr_dP_dn  = 0.d0
          this%krg        = this%kr        * this%conductance
          this%dkrg_dP_up = this%dkr_dP_up * this%conductance
          this%dkrg_dP_dn = this%dkr_dP_dn * this%conductance

       else

          if (this%satParams_up%relperm_func_type == 0) then

             frac_liq_sat = 1.d0

             call SatFunc_PressToRelPerm( &
                  this%satParams_dn, this%pressure_dn, &
                  frac_liq_sat, kr_dn, dkr_dn_dP_dn)

             this%kr         = kr_dn
             this%dkr_dP_up  = 0.d0
             this%dkr_dP_dn  = dkr_dn_dP_dn
             this%krg        = this%kr        * this%conductance
             this%dkrg_dP_up = this%dkr_dP_up * this%conductance
             this%dkrg_dP_dn = this%dkr_dP_dn * this%conductance

          elseif (this%satParams_dn%relperm_func_type == 0) then

             call SatFunc_PressToRelPerm( &
                  this%satParams_up, this%pressure_up, &
                  frac_liq_sat, kr_up, dkr_up_dP_up)

             this%kr         = kr_up
             this%dkr_dP_up  = dkr_up_dP_up
             this%dkr_dP_dn  = 0.d0
             this%krg        = this%kr        * this%conductance
             this%krg        = this%kr        * this%conductance
             this%dkrg_dP_up = this%dkr_dP_up * this%conductance
             this%dkrg_dP_dn = this%dkr_dP_dn * this%conductance

          else

             call SatFunc_PressToRelPerm( &
                  this%satParams_up, this%pressure_up, &
                  frac_liq_sat, kr_up, dkr_up_dP_up)

             call SatFunc_PressToRelPerm( &
                  this%satParams_dn, this%pressure_dn, &
                  frac_liq_sat, kr_dn, dkr_dn_dP_dn)

             this%kr = &
                  (       this%conductance_upwind_weight  * kr_up + &
                  (1.d0 - this%conductance_upwind_weight) * kr_dn)

             this%dkr_dP_up =  &
                  this%conductance_upwind_weight          * dkr_up_dP_up

             this%dkr_dP_dn = &
                  (1.d0 - this%conductance_upwind_weight) * dkr_dn_dP_dn

             this%krg        = this%kr        * this%conductance
             this%dkrg_dP_up = this%dkr_dP_up * this%conductance
             this%dkrg_dP_dn = this%dkr_dP_dn * this%conductance
          end if

       endif

    case (CONDUCTANCE_MANOLI_TYPE)

       call SatFunc_PressToRelPerm(              &
            this%satParams_up, this%pressure_up, &
            frac_liq_sat, kr_up, dkr_up_dP_up)

       call SatFunc_PressToRelPerm(              &
            this%satParams_dn, this%pressure_dn, &
            frac_liq_sat, kr_dn, dkr_dn_dP_dn)

       krg_up          = kr_up * this%conductance_up
       krg_dn          = kr_dn * this%conductance_dn
       denom           = krg_up + krg_dn

       this%krg        = krg_up * krg_dn / denom
       this%dkrg_dP_up = (krg_dn/denom)**2.d0 * dkr_up_dP_up * this%conductance_up
       this%dkrg_dP_dn = (krg_up/denom)**2.d0 * dkr_dn_dP_dn * this%conductance_dn

    case default
       write(iulog,*) 'Unknown conductance_type ', this%conductance_type
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine RichODEPressureConnAuxVarCompute

  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetUpwindPressure(this, val)
    !
    ! !DESCRIPTION:
    ! Set upwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%pressure_up = val
 
  end subroutine RichODEPressureConnAuxSetUpwindPressure
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetDownwindPressure(this, val)
    !
    ! !DESCRIPTION:
    ! Set downwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%pressure_dn = val
 
  end subroutine RichODEPressureConnAuxSetDownwindPressure
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetRelativePermeability(this, val)
    !
    ! !DESCRIPTION:
    ! Set relative permeability
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%kr = val
 
  end subroutine RichODEPressureConnAuxSetRelativePermeability
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetFluxType(this, val)
    !
    ! !DESCRIPTION:
    ! Set flux type
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscInt                              :: val
 
    this%flux_type = val
 
  end subroutine RichODEPressureConnAuxSetFluxType
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetConductanceType(this, val)
    !
    ! !DESCRIPTION:
    ! Set conductance type
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscInt                              :: val
 
    this%conductance_type = val
 
  end subroutine RichODEPressureConnAuxSetConductanceType
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetConductance(this, val)
    !
    ! !DESCRIPTION:
    ! Set conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%conductance = val
 
  end subroutine RichODEPressureConnAuxSetConductance
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetUpwindConductance(this, val)
    !
    ! !DESCRIPTION:
    ! Set upwind conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%conductance_up = val
 
  end subroutine RichODEPressureConnAuxSetUpwindConductance
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetDownwindConductance(this, val)
    !
    ! !DESCRIPTION:
    ! Set downwind conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%conductance_dn = val
 
  end subroutine RichODEPressureConnAuxSetDownwindConductance
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetDKrDPup(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of relative permeability w.r.t. upwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dkr_dP_up = val
 
  end subroutine RichODEPressureConnAuxSetDKrDPup
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetDKrDPdn(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of relative permeability w.r.t. downwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dkr_dP_dn = val
 
  end subroutine RichODEPressureConnAuxSetDKrDPdn
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetKrg(this, val)
    !
    ! !DESCRIPTION:
    ! Set product of relative permeability and conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%krg = val
 
  end subroutine RichODEPressureConnAuxSetKrg
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetDKrgDPup(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of the product of relative permeability and conductance w.r.t. upwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dkrg_dP_up = val
 
  end subroutine RichODEPressureConnAuxSetDKrgDPup
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetDKrgDPdn(this, val)
    !
    ! !DESCRIPTION:
    ! Set derivative of the product of relative permeability and conductance w.r.t. downwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dkrg_dP_dn = val
 
  end subroutine RichODEPressureConnAuxSetDKrgDPdn
 
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetUpwindConductanceWeight(this, val)
    !
    ! !DESCRIPTION:
    ! Set weight of upwind conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                        :: val
 
    this%conductance_upwind_weight = val
 
  end subroutine RichODEPressureConnAuxSetUpwindConductanceWeight

  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetUpwindPressure(this)
    !
    ! !DESCRIPTION:
    ! Get upwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetUpwindPressure
 
    RichODEPressureConnAuxGetUpwindPressure = this%pressure_up
 
  end function RichODEPressureConnAuxGetUpwindPressure
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetDownwindPressure(this)
    !
    ! !DESCRIPTION:
    ! Get downwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetDownwindPressure
 
    RichODEPressureConnAuxGetDownwindPressure = this%pressure_dn
 
  end function RichODEPressureConnAuxGetDownwindPressure
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetRelativePermeability(this)
    !
    ! !DESCRIPTION:
    ! Get relative permeability
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetRelativePermeability
 
    RichODEPressureConnAuxGetRelativePermeability = this%kr
 
  end function RichODEPressureConnAuxGetRelativePermeability
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetFluxType(this)
    !
    ! !DESCRIPTION:
    ! Get flux type
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscInt                         :: RichODEPressureConnAuxGetFluxType
 
    RichODEPressureConnAuxGetFluxType = this%flux_type
 
  end function RichODEPressureConnAuxGetFluxType

  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetConductanceType(this)
    !
    ! !DESCRIPTION:
    ! Get conductance type
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscInt                         :: RichODEPressureConnAuxGetConductanceType
 
    RichODEPressureConnAuxGetConductanceType = this%conductance_type
 
  end function RichODEPressureConnAuxGetConductanceType

  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetConductance(this)
    !
    ! !DESCRIPTION:
    ! Get conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetConductance
 
    RichODEPressureConnAuxGetConductance = this%conductance
 
  end function RichODEPressureConnAuxGetConductance
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetUpwindConductance(this)
    !
    ! !DESCRIPTION:
    ! Get upwind conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetUpwindConductance
 
    RichODEPressureConnAuxGetUpwindConductance = this%conductance_up
 
  end function RichODEPressureConnAuxGetUpwindConductance
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetDownwindConductance(this)
    !
    ! !DESCRIPTION:
    ! Get downwind conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetDownwindConductance
 
    RichODEPressureConnAuxGetDownwindConductance = this%conductance_dn
 
  end function RichODEPressureConnAuxGetDownwindConductance
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetDKrDPup(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of relative permeability w.r.t. upwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetDKrDPup
 
    RichODEPressureConnAuxGetDKrDPup = this%dkr_dP_up
 
  end function RichODEPressureConnAuxGetDKrDPup
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetDKrDPdn(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of relative permeability w.r.t. downwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetDKrDPdn
 
    RichODEPressureConnAuxGetDKrDPdn = this%dkr_dP_dn
 
  end function RichODEPressureConnAuxGetDKrDPdn
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetKrg(this)
    !
    ! !DESCRIPTION:
    ! Get product of relative permeability and conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetKrg
 
    RichODEPressureConnAuxGetKrg = this%krg
 
  end function RichODEPressureConnAuxGetKrg
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetDKrgDPup(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of the product of relative permeability and conductance w.r.t. upwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetDKrgDPup
 
    RichODEPressureConnAuxGetDKrgDPup = this%dkrg_dP_up
 
  end function RichODEPressureConnAuxGetDKrgDPup
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetDKrgDPdn(this)
    !
    ! !DESCRIPTION:
    ! Get derivative of the product of relative permeability and conductance w.r.t. downwind pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetDKrgDPdn
 
    RichODEPressureConnAuxGetDKrgDPdn = this%dkrg_dP_dn
 
  end function RichODEPressureConnAuxGetDKrgDPdn
 
  !------------------------------------------------------------------------
  function RichODEPressureConnAuxGetUpwindConductanceWeight(this)
    !
    ! !DESCRIPTION:
    ! Get weight of upwind conductance
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(rich_ode_pres_conn_auxvar_type) :: this
    !
    PetscReal                        :: RichODEPressureConnAuxGetUpwindConductanceWeight
 
    RichODEPressureConnAuxGetUpwindConductanceWeight = this%conductance_upwind_weight
 
  end function RichODEPressureConnAuxGetUpwindConductanceWeight

#endif

end module RichardsODEPressureConnAuxType
