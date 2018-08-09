module GoveqnThermalEnthalpySoilType

#ifdef USE_PETSC_LIB
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Govneqn data type allocation
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                              , only : iulog
  use mpp_abortutils                          , only : endrun
  use mpp_shr_log_mod                         , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType               , only : goveqn_base_type
  use ThermalEnthalpySoilAuxType              , only : therm_enthalpy_soil_auxvar_type
  use SystemOfEquationsThermalEnthalpyAuxType , only : sysofeqns_thermal_enthalpy_auxvar_type
  use petscsys
  use petscvec
  use petscmat
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_thermal_enthalpy_soil_type

     Vec :: accum_prev

     type (therm_enthalpy_soil_auxvar_type), pointer :: aux_vars_in(:)  ! Internal state.
     type (therm_enthalpy_soil_auxvar_type), pointer :: aux_vars_bc(:)  ! Boundary conditions.
     type (therm_enthalpy_soil_auxvar_type), pointer :: aux_vars_ss(:)  ! Source-sink.

     PetscInt, pointer                               :: soe_auxvars_bc_offset (:) ! SoE auxvar offset corresponding to BCs
     PetscInt, pointer                               :: soe_auxvars_ss_offset (:) ! SoE auxvar offset corresponding to SSs

   contains

     procedure, public :: Setup                   => ThermEnthalpySoilSetup
     procedure, public :: AllocateAuxVars         => ThermEnthalpySoilAllocateAuxVars

     procedure, public :: GetFromSOEAuxVarsIntrn  => ThermEnthalpySoilGetFromSOEAuxVarsIntrn
     procedure, public :: GetFromSOEAuxVarsBC     => ThermEnthalpySoilGetFromSOEAuxVarsBC
     procedure, public :: GetFromSOEAuxVarsSS     => ThermEnthalpySoilGetFromSOEAuxVarsSS
     procedure, public :: GetDataFromSOEAuxVar    => ThermEnthalpySoilGetDataFromSOEAuxVar

     procedure, public :: SetFromSOEAuxVarsIntrn  => ThermEnthalpySoilSetFromSOEAuxVarsIntrn
     procedure, public :: SetFromSOEAuxVarsBC     => ThermEnthalpySoilSetFromSOEAuxVarsBC
     procedure, public :: SetFromSOEAuxVarsSS     => ThermEnthalpySoilSetFromSOEAuxVarsSS

     procedure, public :: SetDataInSOEAuxVar      => ThermEnthalpySoilSetDataInSOEAuxVar
     procedure, public :: SetSOEAuxVarOffsets     => ThermEnthalpySoilSetSOEAuxVarOffsets
     procedure, public :: SetDensityType          => ThermEnthalpySetDensityType
     procedure, public :: SetIntEnergyEnthalpyType=> ThermEnthalpySetIntEnergyEnthalpyType

     procedure, public :: UpdateAuxVars           => ThermEnthalpySoilUpdateAuxVars
     procedure, public :: UpdateAuxVarsIntrn      => ThermEnthalpySoilUpdateAuxVarsIntrn
     procedure, public :: UpdateAuxVarsBC         => ThermEnthalpySoilUpdateAuxVarsBC
     procedure, public :: UpdateAuxVarsSS         => ThermEnthalpySoilUpdateAuxVarsSS

     procedure, public :: ComputeResidual         => ThermEnthalpySoilComputeResidual
     procedure, public :: ComputeJacobian         => ThermEnthalpySoilComputeJacobian
     procedure, public :: ComputeOffDiagJacobian  => ThermEnthalpySoilComputeOffDiagJacobian
     procedure, public :: PreSolve                => ThermEnthalpySoilPreSolve

     procedure, public :: CreateVectors           => ThermEnthalpySoilCreateVectors

     procedure, public :: SetSoilPermeability     => ThermEnthalpySetSoilPermeability
     procedure, public :: SetRelativePermeability => ThermEnthalpySetRelativePermeability
     procedure, public :: SetSoilPorosity         => ThermEnthalpySetSoilPorosity
     procedure, public :: SetSaturationFunction   => ThermEnthalpySetSaturationFunction
     procedure, public :: SetHeatCapacity         => ThermEnthalpySetHeatCapacity
     procedure, public :: SetThermalCondWet       => ThermEnthalpySetThermalCondWet
     procedure, public :: SetThermalCondDry       => ThermEnthalpySetThermalCondDry
     procedure, public :: SetThermalAlpha         => ThermEnthalpySetThermalAlpha
     procedure, public :: SetSoilDensity          => ThermEnthalpySetSoilDensity

  end type goveqn_thermal_enthalpy_soil_type

  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation for Thermal equation.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_THERM_SOIL_EBASED
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this

    call this%Create()

    this%name       = "Soil thermal equation based on enthalpy"
    this%itype      = GE_THERM_SOIL_EBASED
    this%mesh_itype = MESH_CLM_SOIL_COL

    nullify(this%aux_vars_in)
    nullify(this%aux_vars_bc)
    nullify(this%aux_vars_ss)

    nullify(this%soe_auxvars_bc_offset)
    nullify(this%soe_auxvars_ss_offset)

  end subroutine ThermEnthalpySoilSetup

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates memory for storing auxiliary variables associated with:
    !   + Internal control volumes,
    !   + Boundary condtions,
    !   + Source-sink condition.
    !
    ! !USES:
    use ConditionType, only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    !
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: ncells_cond
    PetscInt                                 :: ncond
    PetscInt                                 :: icond

    ! Allocate memory and initialize aux vars: For internal connections
    allocate(this%aux_vars_in(this%mesh%ncells_all))
    do icond = 1,this%mesh%ncells_all
       call this%aux_vars_in(icond)%Init()
    enddo

    ! Allocate memory and initialize aux vars: For boundary connections
    ncells_cond = 0
    ncond       = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       ncells_cond = ncells_cond + cur_cond%ncells
       ncond       = ncond + 1
       cur_cond => cur_cond%next
    enddo
    allocate(this%aux_vars_bc(ncells_cond))
    do icond = 1,ncells_cond
       call this%aux_vars_bc(icond)%Init()
    enddo
    allocate(this%soe_auxvars_bc_offset(ncond))
    this%soe_auxvars_bc_offset(:) = -1

    ! Allocate memory and initialize aux vars: For source sink connections
    ncells_cond = 0
    ncond       = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       ncells_cond = ncells_cond + cur_cond%ncells
       ncond       = ncond + 1
       cur_cond => cur_cond%next
    enddo
    allocate(this%aux_vars_ss(ncells_cond))
    do icond = 1,ncells_cond
       call this%aux_vars_ss(icond)%Init()
    enddo
    allocate(this%soe_auxvars_ss_offset(ncond))
    this%soe_auxvars_ss_offset(:) = -1

  end subroutine ThermEnthalpySoilAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilGetDataFromSOEAuxVar(this, soe_avar_type, soe_avars, &
       offset)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    !
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)      , intent(inout)        :: this
    PetscInt                                      , intent(in)           :: soe_avar_type
    type (sysofeqns_thermal_enthalpy_auxvar_type) , intent(in)           :: soe_avars(:)
    PetscInt                                      , intent(in), optional :: offset
    !
    ! !LOCAL VARIABLES
    PetscInt                                                             :: iauxvar_off

    select case(soe_avar_type)
    case (AUXVAR_INTERNAL)
       if (present(offset)) then
          iauxvar_off = offset
       else
          iauxvar_off = 0
       endif
       call ThermEnthalpySoilGetFromSOEAuxVarsIntrn(this, soe_avars, iauxvar_off)
    case (AUXVAR_BC)
       call ThermEnthalpySoilGetFromSOEAuxVarsBC(this, soe_avars)
    case (AUXVAR_SS)
       call ThermEnthalpySoilGetFromSOEAuxVarsSS(this, soe_avars)
    case default
       write(iulog,*) 'ThermEnthalpySoilGetDataFromSOEAuxVar: soe_avar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermEnthalpySoilGetDataFromSOEAuxVar

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilGetFromSOEAuxVarsIntrn(this, soe_avars, offset)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for internal nodes
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)     , intent(inout) :: this
    type(sysofeqns_thermal_enthalpy_auxvar_type) , intent(in)    :: soe_avars(:)
    PetscInt                                     , intent(in)    :: offset
    !
    ! LOCAL VARIABLES
    PetscInt :: iauxvar
    PetscInt :: nauxvar

    nauxvar = size(this%aux_vars_in)
    if( nauxvar > size(soe_avars) ) then
       write(iulog,*) 'size(this%aux_vars_in) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do iauxvar = 1, nauxvar
       this%aux_vars_in(iauxvar)%temperature = soe_avars(iauxvar+offset)%temperature
       this%aux_vars_in(iauxvar)%pressure    = soe_avars(iauxvar+offset)%pressure
    enddo

  end subroutine ThermEnthalpySoilGetFromSOEAuxVarsIntrn

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilSetFromSOEAuxVarsIntrn(this, var_type, ndata, data)
    !
    ! !DESCRIPTION:
    ! Saves values in the GE auxiliary variable for internal nodes.
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) , intent(inout)       :: this
    PetscInt                                 , intent(in)          :: var_type
    PetscInt                                 , intent(in)          :: ndata
    PetscReal                                , intent(in), pointer :: data(:)
    !
    ! LOCAL VARIABLES
    PetscInt :: iauxvar
    PetscInt :: nauxvar

    nauxvar = size(this%aux_vars_in)
    if( nauxvar /= ndata ) then
       write(iulog,*) 'size(this%aux_vars_in) /= ndata)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case (var_type)
    case (VAR_TEMPERATURE)
       do iauxvar = 1, nauxvar
          this%aux_vars_in(iauxvar)%temperature = data(iauxvar)
       enddo

    case (VAR_PRESSURE)
       do iauxvar = 1, nauxvar
          this%aux_vars_in(iauxvar)%pressure    = data(iauxvar)
       enddo

    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermEnthalpySoilSetFromSOEAuxVarsIntrn

  !------------------------------------------------------------------------

  subroutine ThermEnthalpySoilGetFromSOEAuxVarsBC(this, soe_avars)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for bondary conditions
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) , intent(inout)               :: this
    type(sysofeqns_thermal_enthalpy_auxvar_type), dimension(:), intent(in) :: soe_avars
    !
    ! LOCAL VARIABLES
    integer                                         :: iauxvar
    integer                                         :: iauxvar_off 
    integer                                         :: iconn
    integer                                         :: nauxVar_ge
    integer                                         :: nauxVar_soe
    integer                                         :: condition_id
    integer                                         :: sum_conn
    integer                                         :: cell_id
    type(condition_type)      , pointer             :: cur_cond
    type(connection_set_type) , pointer             :: cur_conn_set
    character(len=256)                              :: string
    type (therm_enthalpy_soil_auxvar_type), pointer :: ge_avars(:)

    ge_avars => this%aux_vars_bc

    nauxVar_ge = size(this%aux_vars_bc)
    nauxVar_soe = size(soe_avars)
    if( nauxVar_ge > nauxVar_soe ) then
       write(iulog,*) 'size(this%aux_vars_bc) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_cond => cur_cond%next
          cycle
       endif

       ! Find first soe-auxvar corresponding to goveqn-auxvar.
       iauxvar_off = -1
       do iauxvar = 1, nauxVar_soe
          if(  &
               soe_avars(iauxvar)%is_bc  &
               .and.  &
               soe_avars(iauxvar)%goveqn_id == this%rank_in_soe_list  &
               .and.  &
               soe_avars(iauxvar)%condition_id == condition_id  &
               ) then
             iauxvar_off = iauxvar - 1
             exit
          end if
       end do
       if (iauxvar_off < 0) then
          write(iulog,*) 'iauxvar_off < 0'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       cur_conn_set => cur_cond%conn_set
       do iconn = 1, cur_conn_set%num_connections

          sum_conn = sum_conn + 1
          cell_id = cur_conn_set%conn(iconn)%GetIDDn()

          select case(cur_cond%itype)
          case (COND_DIRICHLET, COND_HEAT_FLUX)
             ge_avars(sum_conn)%condition_value =  &
                  soe_avars(iconn + iauxvar_off)%condition_value

          case (COND_DIRICHLET_FRM_OTR_GOVEQ)
             ! Do nothing

          case default
             write(string,*) cur_cond%itype
             write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       enddo

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySoilGetFromSOEAuxVarsBC

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilSetFromSOEAuxVarsBC(this, var_type, ndata, data)
    !
    ! !DESCRIPTION:
    ! Saves values in the GE auxiliary variable for boundary conditions.
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) , intent(inout)       :: this
    PetscInt                                 , intent(in)          :: var_type
    PetscInt                                 , intent(in)          :: ndata
    PetscReal                                , intent(in), pointer :: data(:)
    !
    ! LOCAL VARIABLES
    PetscInt                                                       :: nauxvar
    integer                                                        :: iconn
    integer                                                        :: condition_id
    integer                                                        :: sum_conn
    integer                                                        :: data_counter
    type(condition_type)                     , pointer             :: cur_cond
    type(connection_set_type)                , pointer             :: cur_conn_set
    type (therm_enthalpy_soil_auxvar_type)   , pointer             :: ge_avars(:)

    ge_avars => this%aux_vars_bc

    nauxvar = size(this%aux_vars_bc)
    if( ndata > nauxvar ) then
       write(iulog,*) 'ndata > size(this%aux_vars_bc)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    data_counter = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_conn_set => cur_cond%conn_set
          sum_conn = sum_conn + cur_conn_set%num_connections
          cur_cond => cur_cond%next
          cycle
       endif

       cur_conn_set => cur_cond%conn_set
       do iconn = 1, cur_conn_set%num_connections

          sum_conn = sum_conn + 1
          data_counter = data_counter + 1

          select case(cur_cond%itype)
          case (COND_DIRICHLET)
             if (var_type == VAR_BC_SS_CONDITION) ge_avars(sum_conn)%condition_value = data(data_counter)
             !if (var_type == VAR_PRESSURE       ) ge_avars(sum_conn)%pressure        = data(data_counter)
             !if (var_type == VAR_TEMPERATURE    ) ge_avars(sum_conn)%condition_value = data(data_counter)
          case (COND_HEAT_FLUX)
             if (var_type == VAR_BC_SS_CONDITION) ge_avars(sum_conn)%condition_value =  data(data_counter)

          case (COND_DIRICHLET_FRM_OTR_GOVEQ)
             ! Do nothing

          case default
             write(iulog,*) 'Unknown cur_cond%itype = ', cur_cond%itype
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       enddo

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySoilSetFromSOEAuxVarsBC

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilGetFromSOEAuxVarsSS(this, soe_avars)
    !
    ! !DESCRIPTION:
    ! Copies values from SoE auxiliary variable into GE auxiliary variable
    ! for bondary conditions
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type), intent(inout)       :: this
    type(sysofeqns_thermal_enthalpy_auxvar_type), dimension(:), intent(in) :: soe_avars
    !
    ! LOCAL VARIABLES
    integer                             :: iauxvar
    integer                             :: iauxvar_off
    integer                             :: iconn
    integer                             :: nauxVar_ge
    integer                             :: nauxVar_soe
    integer                             :: condition_id
    integer                             :: sum_conn
    PetscReal                           :: var_value
    type(condition_type)      , pointer :: cur_cond
    type(connection_set_type) , pointer :: cur_conn_set
    character(len=256)                  :: string

    nauxVar_ge  = size(this%aux_vars_ss)
    nauxVar_soe = size(soe_avars)

    if( nauxVar_ge > nauxVar_soe ) then
       write(iulog,*) 'size(this%aux_vars_ss) > size(soe_avars)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       ! Find first soe-auxvar corresponding to goveqn-auxvar.
       iauxvar_off = -1
       do iauxvar = 1, nauxVar_soe
          if(  &
               soe_avars(iauxvar)%is_ss  &
               .and.  &
               soe_avars(iauxvar)%goveqn_id == this%rank_in_soe_list  &
               .and.  &
               soe_avars(iauxvar)%condition_id == condition_id  &
               ) then
             iauxvar_off = iauxvar - 1
             exit
          end if
       end do
       if (iauxvar_off < 0) then
          write(iulog,*) 'iauxvar_off < 0'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

       cur_conn_set => cur_cond%conn_set
       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          select case(cur_cond%itype)
          case (COND_HEAT_RATE)

             var_value = soe_avars(iconn + iauxvar_off)%condition_value

             this%aux_vars_ss(sum_conn)%condition_value = var_value
             cur_cond%value(iconn)                      = var_value
          case default
             write(string,*) cur_cond%itype
             write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       enddo

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySoilGetFromSOEAuxVarsSS

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilSetFromSOEAuxVarsSS(this, var_type, ndata, data)
    !
    ! !DESCRIPTION:
    ! Saves values in the GE auxiliary variable for source-sink conditions.
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) , intent(inout)       :: this
    PetscInt                                 , intent(in)          :: var_type
    PetscInt                                 , intent(in)          :: ndata
    PetscReal                                , intent(in), pointer :: data(:)
    !
    ! LOCAL VARIABLES
    PetscInt                                                       :: nauxvar
    integer                                                        :: iconn
    integer                                                        :: condition_id
    integer                                                        :: sum_conn
    integer                                                        :: data_counter
    type(condition_type)                     , pointer             :: cur_cond
    type(connection_set_type)                , pointer             :: cur_conn_set
    type (therm_enthalpy_soil_auxvar_type)   , pointer             :: ge_avars(:)

    ge_avars => this%aux_vars_ss

    nauxvar = size(this%aux_vars_ss)
    if( ndata > nauxvar ) then
       write(iulog,*) 'ndata > size(this%aux_vars_bc)'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    condition_id = 0
    sum_conn = 0
    data_counter = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       condition_id = condition_id + 1

       cur_conn_set => cur_cond%conn_set
       do iconn = 1, cur_conn_set%num_connections

          sum_conn = sum_conn + 1
          data_counter = data_counter + 1

          select case(cur_cond%itype)
          case (COND_HEAT_RATE)
             ge_avars(sum_conn)%condition_value = data(data_counter)
             cur_cond%value(iconn)              = data(data_counter)

          case default
             write(iulog,*) 'Unknown cur_cond%itype = ', cur_cond%itype
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select
       enddo

       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySoilSetFromSOEAuxVarsSS

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilSetDataInSOEAuxVar(this, soe_avar_type, soe_avars, &
       offset)
    !
    ! !DESCRIPTION:
    ! Copies data GE auxiliary variable into SoE auxiliary variable.
    ! This is done as a part of post solve.
    !
    use MultiPhysicsProbConstants    , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)                        :: this
    PetscInt                                                        :: soe_avar_type
    type (sysofeqns_thermal_enthalpy_auxvar_type), dimension(:), intent(out) :: soe_avars
    PetscInt, optional                                              :: offset
    !
    ! !LOCAL VARIABLES
    PetscInt                                                       :: iauxvar
    PetscInt                                                       :: iauxvar_off

    if (present(offset)) then
       iauxvar_off = offset
    else
       iauxvar_off = 0
    endif

    select case(soe_avar_type)
    case (AUXVAR_INTERNAL)
       this%aux_vars_in => this%aux_vars_in

       if ( size(this%aux_vars_in) > size(soe_avars)) then
          write(iulog,*) 'size(this%aux_vars_in) > size(soe_avars)'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       do iauxvar = 1, size(this%aux_vars_in)
          if (this%mesh%is_active(iauxvar)) then
             soe_avars(iauxvar+iauxvar_off)%temperature =  &
                  this%aux_vars_in(iauxvar)%temperature
          endif
       enddo

    case default
       write(iulog,*) 'ThermEnthalpySetDataInSOEAuxVar: soe_avar_type not supported'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermEnthalpySoilSetDataInSOEAuxVar

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilSetSOEAuxVarOffsets(this, bc_offset_count, bc_offsets, &
    ss_offset_count, ss_offsets)
    !
    ! !DESCRIPTION:
    !
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    PetscInt                                 :: bc_offset_count
    PetscInt, pointer                        :: bc_offsets(:)
    PetscInt                                 :: ss_offset_count
    PetscInt, pointer                        :: ss_offsets(:)
    !
    ! !LOCAL VARIABLES
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: cond_count
    PetscInt                                 :: cond_itype_to_exclude

    cond_itype_to_exclude = COND_NULL
    call this%GetNConditionsExcptCondItype(COND_BC, &
         cond_itype_to_exclude, cond_count)

    if (bc_offset_count > cond_count) then
       write(iulog,*) 'ERROR: bc_offset_count > cond_count'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cond_itype_to_exclude = COND_NULL
    call this%GetNConditionsExcptCondItype(COND_SS, &
         cond_itype_to_exclude, cond_count)

    if (ss_offset_count > cond_count) then
       write(iulog,*) 'ERROR: ss_offset_count > cond_count'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cond_count = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       cond_count = cond_count + 1
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          this%soe_auxvars_bc_offset(cond_count) = bc_offsets(cond_count)
       endif
       cur_cond => cur_cond%next
    enddo

    cond_count = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       cond_count = cond_count + 1
       this%soe_auxvars_ss_offset(cond_count) = bc_offsets(cond_count)
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySoilSetSOEAuxVarOffsets

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetDensityType(this, density_type)
    !
    ! !DESCRIPTION:
    ! Set type of density formulation to auxiliary variables associated with
    ! internal, boundary, and source-sink conditions.
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : VAR_DENSITY_TYPE
    use EOSWaterMod               , only : DENSITY_CONSTANT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use EOSWaterMod               , only : DENSITY_IFC67
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    PetscInt                                 :: density_type
    !
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: sum_conn
    PetscInt                                 :: icond

    if (density_type /= DENSITY_CONSTANT .and. &
        density_type /= DENSITY_TGDPB01  .and. &
        density_type /= DENSITY_IFC67          &
        ) then
       write(iulog,*) 'Unknown value for VAR_DENSITY_TYPE: ',density_type
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! For internal connections
    do icond = 1,this%mesh%ncells_all
       this%aux_vars_in(icond)%density_type = density_type
    enddo

    ! For boundary conditions
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       do icond = 1,cur_cond%ncells
          sum_conn = sum_conn + 1
          this%aux_vars_bc(sum_conn)%density_type = density_type
       enddo
       cur_cond => cur_cond%next
    enddo

    ! For source sink conditions
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       do icond = 1,cur_cond%ncells
          sum_conn = sum_conn + 1
          this%aux_vars_ss(sum_conn)%density_type = density_type
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySetDensityType

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetIntEnergyEnthalpyType(this, itype)
    !
    ! !DESCRIPTION:
    ! Set type of density formulation to auxiliary variables associated with
    ! internal, boundary, and source-sink conditions.
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : VAR_DENSITY_TYPE
    use EOSWaterMod               , only : INT_ENERGY_ENTHALPY_CONSTANT
    use EOSWaterMod               , only : INT_ENERGY_ENTHALPY_IFC67
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    PetscInt                                 :: itype
    !
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: sum_conn
    PetscInt                                 :: icond

    if (itype /= INT_ENERGY_ENTHALPY_CONSTANT .and. &
        itype /= INT_ENERGY_ENTHALPY_IFC67          &
        ) then
       write(iulog,*) 'Unknown value for Internal-Energy & Enthalpy type: ',itype
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! For internal connections
    do icond = 1,this%mesh%ncells_all
       this%aux_vars_in(icond)%int_energy_enthalpy_type = itype
    enddo

    ! For boundary conditions
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       do icond = 1,cur_cond%ncells
          sum_conn = sum_conn + 1
          this%aux_vars_bc(sum_conn)%int_energy_enthalpy_type = itype
       enddo
       cur_cond => cur_cond%next
    enddo

    ! For source sink conditions
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       do icond = 1,cur_cond%ncells
          sum_conn = sum_conn + 1
          this%aux_vars_ss(sum_conn)%int_energy_enthalpy_type = itype
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySetIntEnergyEnthalpyType

  !------------------------------------------------------------------------

  subroutine ThermEnthalpySoilUpdateAuxVars(this)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this

    call this%UpdateAuxVarsIntrn()
    call this%UpdateAuxVarsBC()
    call this%UpdateAuxVarsSS()

  end subroutine ThermEnthalpySoilUpdateAuxVars

  !------------------------------------------------------------------------

  subroutine ThermEnthalpySoilUpdateAuxVarsIntrn(this)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    !
    ! !LOCAL VARIABLES
    PetscInt :: ghosted_id

    ! Update aux vars for internal cells
    do ghosted_id = 1, this%mesh%ncells_local
       call this%aux_vars_in(ghosted_id)%Compute()
    enddo

  end subroutine ThermEnthalpySoilUpdateAuxVarsIntrn

  !------------------------------------------------------------------------

  subroutine ThermEnthalpySoilUpdateAuxVarsBC(this)
    !
    ! !DESCRIPTION:
    !
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ghosted_id
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscReal                                :: temperature
    type(condition_type),pointer             :: cur_cond
    class(connection_set_type), pointer      :: cur_conn_set
    character(len=256)                       :: string

    ! Update aux vars for boundary cells
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          select case(cur_cond%itype)
          case (COND_DIRICHLET)
             this%aux_vars_bc(sum_conn)%temperature = &
                  this%aux_vars_bc(sum_conn)%condition_value

          case (COND_DIRICHLET_FRM_OTR_GOVEQ)
             ! Do nothing

          case default
             write(string,*) cur_cond%itype
             write(iulog,*) 'Unknown cur_cond%itype = ' // trim(string)
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end select

          call this%aux_vars_bc(sum_conn)%Compute()
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySoilUpdateAuxVarsBC

  !------------------------------------------------------------------------

  subroutine ThermEnthalpySoilUpdateAuxVarsSS(this)
    !
    ! !DESCRIPTION:
    !
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: ghosted_id
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    type(condition_type),pointer             :: cur_cond
    class(connection_set_type), pointer      :: cur_conn_set

    ! Update aux vars for source/sink cells
    sum_conn = 0
    cur_cond => this%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()
          this%aux_vars_ss(sum_conn)%temperature =  &
              this%aux_vars_in(ghosted_id)%temperature
          call this%aux_vars_ss(sum_conn)%Compute()
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermEnthalpySoilUpdateAuxVarsSS

  !------------------------------------------------------------------------

  subroutine ThermEnthalpySoilComputeResidual(this, X, F, ierr)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    Vec                                      :: X
    Vec                                      :: F
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    PetscReal, pointer                      :: f_p(:)
    PetscReal, pointer                      :: accum_prev_p(:)
    
    call VecGetArrayF90(F, f_p, ierr); CHKERRQ(ierr);
    call VecGetArrayF90(this%accum_prev, accum_prev_p, ierr); CHKERRQ(ierr)

    f_p(:) = 0.d0

    call ThermalEnthalpySoilAccum(this, f_p)
    f_p(:) = f_p(:) - accum_prev_p(:)

    call ThermalEnthalpySoilDivergence(this, f_p)

    call VecRestoreArrayF90(this%accum_prev, accum_prev_p, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(F, F_p, ierr); CHKERRQ(ierr)

  end subroutine ThermEnthalpySoilComputeResidual

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilComputeJacobian(this, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the jacobian matrix for the discretized energy equation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    Vec                                      :: X
    Mat                                      :: A
    Mat                                      :: B
    PetscErrorCode                           :: ierr

    ! Computes the following:
    ! \sum (\rho \mathbf{q})_{i,j} \cdot \mathbf{n}_{i,j} A_{i,j}  - Q_i
    call ThermalEnthalpySoilDivergenceDeriv(this, B, ierr); CHKERRQ(ierr)

    ! \left( \frac{(\phi * s * \rho)_i V_i}{dt} \right)^{k+1}
    call ThermalEnthalpySoilAccumDeriv(this, B, ierr); CHKERRQ(ierr)

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    if ( A /= B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
       call MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    endif

  end subroutine ThermEnthalpySoilComputeJacobian

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilComputeOffDiagJacobian(this, X_1, X_2, A, B, &
       itype_of_other_goveq, rank_of_other_goveq,        &
       ierr)
    !
    ! !DESCRIPTION:
    ! Computes the off-diagonal jacobian sub matrix associated with the
    ! coupling the given governing equation with another governing equation.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : GE_THERM_SOIL_EBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    Vec                                      :: X_1
    Vec                                      :: X_2
    Mat                                      :: A
    Mat                                      :: B
    PetscInt                                 :: itype_of_other_goveq
    PetscInt                                 :: rank_of_other_goveq
    PetscErrorCode                           :: ierr
    !
    ! LOCAL VARIABLES
    character(len=256)                       :: string

    select case(itype_of_other_goveq)
    case (GE_RE)
       call OffDiagJacobian_Pressure(this, rank_of_other_goveq, &
            B, ierr)
    case (GE_THERM_SOIL_EBASED)
       call OffDiagJacobian_Temperature(this, rank_of_other_goveq, &
            B, ierr)
    case default
       write(string,*) itype_of_other_goveq
       write(iulog,*) 'Unknown id_of_other_goveq = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    if ( A /= B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
       call MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    endif

  end subroutine ThermEnthalpySoilComputeOffDiagJacobian

  !------------------------------------------------------------------------

  subroutine ThermalEnthalpySoilAccum(geq_soil, f_p)
    !
    ! !DESCRIPTION:
    !
    ! \frac{\partial \xi}{\partial dt} dV
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: geq_soil
    PetscReal, intent(out)     :: f_p(:)
    !
    ! !LOCAL VARIABLES
    PetscInt                                 :: cell_id
    PetscReal                                :: dtInv
    type (therm_enthalpy_soil_auxvar_type), pointer :: aux_vars(:)

    aux_vars => geq_soil%aux_vars_in
    
    dtInv = 1.d0/geq_soil%dtime

    ! Interior cells
    do cell_id = 1, geq_soil%mesh%ncells_local

       if (geq_soil%mesh%is_active(cell_id)) then

          ! phi*rho_liq*sat_liq*u_liq + (1-phi)*rho_soil*u_soil

          f_p(cell_id) = f_p(cell_id) +           &
               (                                  &
                aux_vars(cell_id)%por  *          &
                aux_vars(cell_id)%den *          &
                aux_vars(cell_id)%sat *          &
                aux_vars(cell_id)%ul              &
                +                                 &
                (1.d0 - aux_vars(cell_id)%por)  * &
                aux_vars(cell_id)%den_soil      * &
                aux_vars(cell_id)%heat_cap_soil * &
                (aux_vars(cell_id)%temperature-273.15d0)     &
                ) * geq_soil%mesh%vol(cell_id) * dtInv
       endif
    enddo

  end subroutine ThermalEnthalpySoilAccum

  !------------------------------------------------------------------------

  subroutine ThermalEnthalpySoilAccumDeriv(geq_soil, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! \frac{\partial \xi}{\partial dt} dV
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)        :: geq_soil
    Mat                                             :: B
    PetscErrorCode                                  :: ierr
    !
    ! !LOCAL VARIABLES
    type (therm_enthalpy_soil_auxvar_type), pointer :: aux_vars(:)
    PetscInt                                        :: cell_id
    PetscInt                                        :: row
    PetscInt                                        :: col
    PetscReal                                       :: dtInv
    PetscReal                                       :: por
    PetscReal                                       :: den
    PetscReal                                       :: sat
    PetscReal                                       :: ul
    PetscReal                                       :: dden_dT 
    PetscReal                                       :: dsat_dT
    PetscReal                                       :: dul_dT
    PetscReal                                       :: derivative
    PetscReal                                       :: den_soil    
    PetscReal                                       :: heat_cap_soil

    aux_vars => geq_soil%aux_vars_in
    
    dtInv = 1.d0/geq_soil%dtime

    ! Interior cells
    do cell_id = 1, geq_soil%mesh%ncells_local

       if (geq_soil%mesh%is_active(cell_id)) then

          por           = aux_vars(cell_id)%por
          den_soil      = aux_vars(cell_id)%den_soil
          heat_cap_soil = aux_vars(cell_id)%heat_cap_soil

          den          = aux_vars(cell_id)%den
          sat          = aux_vars(cell_id)%sat
          ul            = aux_vars(cell_id)%ul

          dden_dT      = aux_vars(cell_id)%dden_dT
          dsat_dT      = aux_vars(cell_id)%dsat_dT
          dul_dT        = aux_vars(cell_id)%dul_dT
          
          ! d/dT( phi*rho_liq*sat_liq*u_liq + (1-phi)*rho_soil*u_soil)

          derivative =                      &
               por * dden_dT * sat * ul + &
               por * den * dsat_dT * ul + &
               por * den * sat * dul_dT       
               
          derivative = (derivative + (1.d0 - por)* den_soil * heat_cap_soil)* &
               geq_soil%mesh%vol(cell_id) * dtInv                
       else
          derivative = 1.d0
       endif

       row = cell_id - 1
       col = cell_id - 1

       call MatSetValuesLocal(B, 1, row, 1, col, derivative, ADD_VALUES, ierr); CHKERRQ(ierr)
    enddo

  end subroutine ThermalEnthalpySoilAccumDeriv

  !------------------------------------------------------------------------

  subroutine ThermalEnthalpySoilDivergence(geq_soil, f_p)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : COND_NULL
    use ThermalEnthalpyMod        , only : ThermalEnthalpyFlux
    use RichardsMod               , only : RichardsFlux
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: geq_soil
    PetscReal, intent(out)                   :: f_p (:)
    !
    ! !LOCAL VARIABLES
    type(condition_type),pointer             :: cur_cond
    class(connection_set_type), pointer      :: cur_conn_set
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscReal                                :: mflux
    PetscReal                                :: dmflux_dT_up
    PetscReal                                :: dmflux_dT_dn
    PetscReal                                :: eflux
    PetscReal                                :: area
    PetscReal                                :: dummy_var1
    PetscReal                                :: dummy_var2
    PetscBool                                :: compute_deriv
    PetscBool                                :: internal_conn
    PetscBool                                :: swap_order
    PetscInt                                 :: cond_type

    compute_deriv = PETSC_FALSE

    ! Interior cells
    cur_conn_set => geq_soil%mesh%intrn_conn_set_list%first
    sum_conn   = 0
    swap_order = PETSC_FALSE
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(iconn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()

          internal_conn = PETSC_TRUE
          cond_type     = COND_NULL

          if ((.not.geq_soil%mesh%is_active(cell_id_up)) .or. &
              (.not.geq_soil%mesh%is_active(cell_id_dn))) cycle

          call RichardsFlux(                     &
               geq_soil%aux_vars_in(cell_id_up), &
               geq_soil%aux_vars_in(cell_id_dn), &
               cur_conn_set%conn(iconn),         &
               compute_deriv,                    &
               internal_conn,                    &
               swap_order,                       &
               cond_type,                        &
               mflux,                            &
               dummy_var1,                       &
               dummy_var2                        &
               )

          dmflux_dT_up = 0.d0
          dmflux_dT_dn = 0.d0

          call ThermalEnthalpyFlux(              &
               geq_soil%aux_vars_in(cell_id_up), &
               geq_soil%aux_vars_in(cell_id_dn), &
               mflux,                            &
               dmflux_dT_up,                     &
               dmflux_dT_dn,                     &
               cur_conn_set%conn(iconn),         &
               compute_deriv,                    &
               internal_conn,                    &
               cond_type,                        &
               eflux,                            &
               dummy_var1,                       &
               dummy_var2                        &
               )

          f_p(cell_id_up) = f_p(cell_id_up) - eflux
          f_p(cell_id_dn) = f_p(cell_id_dn) + eflux

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    ! Boundary cells
    cur_cond => geq_soil%boundary_conditions%first
    sum_conn   = 0
    swap_order = PETSC_FALSE
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          cell_id  = cur_conn_set%conn(iconn)%GetIDDn()
          sum_conn = sum_conn + 1

          if (.not.geq_soil%mesh%is_active(cell_id )) cycle
          
          select case(cur_cond%itype)
          case(COND_DIRICHLET_FRM_OTR_GOVEQ, COND_DIRICHLET)

             internal_conn = PETSC_FALSE
             cond_type     = cur_cond%itype

             call RichardsFlux(                   &
                  geq_soil%aux_vars_bc(sum_conn), &
                  geq_soil%aux_vars_in(cell_id ), &
                  cur_conn_set%conn(iconn),       &
                  compute_deriv,                  &
                  internal_conn,                  &
                  swap_order,                     &
                  cond_type,                      &
                  mflux,                          &
                  dummy_var1,                     &
                  dummy_var2                      &
                  )

             dmflux_dT_up = 0.d0
             dmflux_dT_dn = 0.d0

             call ThermalEnthalpyFlux(             &
                  geq_soil%aux_vars_bc(sum_conn ), &
                  geq_soil%aux_vars_in(cell_id  ), &
                  mflux,                           &
                  dmflux_dT_up,                    &
                  dmflux_dT_dn,                    &
                  cur_conn_set%conn(iconn),        &
                  compute_deriv,                   &
                  internal_conn,                   &
                  cond_type,                       &
                  eflux,                           &
                  dummy_var1,                      &
                  dummy_var2                       &
                  )

             f_p(cell_id) = f_p(cell_id) + eflux

          case (COND_HEAT_FLUX)             
             area = cur_conn_set%conn(iconn)%GetArea()

             f_p(cell_id) = f_p(cell_id) + cur_cond%value(iconn)*area

          case default
           write(iulog,*)'ThermalKSPTempSoilDivergence: Unknown boundary condition type'
           call endrun(msg=errMsg(__FILE__, __LINE__))

        end select

       enddo
       cur_cond => cur_cond%next
    enddo    

    ! Source-sink cells
    cur_cond => geq_soil%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          cell_id = cur_conn_set%conn(iconn)%GetIDDn()

          if ((.not.geq_soil%mesh%is_active(cell_id))) cycle

          select case(cur_cond%itype)
          case(COND_HEAT_RATE)
             f_p(cell_id) = f_p(cell_id) + cur_cond%value(iconn)

          case default
           write(iulog,*)'ThermalKSPTempSoilDivergence: Unknown source-sink condition type'
           call endrun(msg=errMsg(__FILE__, __LINE__))

        end select

     enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermalEnthalpySoilDivergence

  !------------------------------------------------------------------------

  subroutine ThermalEnthalpySoilDivergenceDeriv(geq_soil, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ, COND_DIRICHLET
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : COND_NULL
    use ThermalEnthalpyMod        , only : ThermalEnthalpyFlux
    use RichardsMod               , only : RichardsFlux, RichardsFluxDerivativeWrtTemperature
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: geq_soil
    Mat                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    type(condition_type),pointer             :: cur_cond
    class(connection_set_type), pointer      :: cur_conn_set
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscInt                                 :: row
    PetscInt                                 :: col
    PetscReal                                :: flux
    PetscReal                                :: area
    PetscReal                                :: Jup
    PetscReal                                :: Jdn
    PetscBool                                :: compute_deriv
    PetscBool                                :: internal_conn
    PetscBool                                :: swap_order
    PetscInt                                 :: cond_type
    PetscReal                                :: val
    PetscReal                                :: mflux
    PetscReal                                :: eflux
    PetscReal                                :: dmflux_dT_up
    PetscReal                                :: dmflux_dT_dn

    compute_deriv = PETSC_TRUE

    ! Interior cells
    cur_conn_set => geq_soil%mesh%intrn_conn_set_list%first
    sum_conn   = 0
    swap_order = PETSC_FALSE
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(iconn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()

          internal_conn = PETSC_TRUE
          cond_type     = COND_NULL
          swap_order    = PETSC_FALSE

          if ((.not.geq_soil%mesh%is_active(cell_id_up)) .or. &
              (.not.geq_soil%mesh%is_active(cell_id_dn))) cycle

          call RichardsFluxDerivativeWrtTemperature(         &
               geq_soil%aux_vars_in(cell_id_up),             &
               geq_soil%aux_vars_in(cell_id_dn),             &
               cur_conn_set%conn(iconn),                     &
               compute_deriv,                                &
               internal_conn,                                &
               swap_order,                                   &
               cond_type,                                    &
               mflux,                                        &
               dmflux_dT_up,                                 &
               dmflux_dT_dn                                  &
               )

          call ThermalEnthalpyFlux(              &
               geq_soil%aux_vars_in(cell_id_up), &
               geq_soil%aux_vars_in(cell_id_dn), &
               mflux,                            &
               dmflux_dT_up,                     &
               dmflux_dT_dn,                     &
               cur_conn_set%conn(iconn),         &
               compute_deriv,                    &
               internal_conn,                    &
               cond_type,                        &
               flux,                             &
               Jup,                              &
               Jdn                               &
               )

          row = cell_id_up - 1
          col = cell_id_up - 1
          val = -Jup
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_up - 1
          col = cell_id_dn - 1
          val = -Jdn
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_dn - 1
          col = cell_id_up - 1
          val = Jup
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_dn - 1
          col = cell_id_dn - 1
          val = Jdn
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    ! Boundary cells
    cur_cond => geq_soil%boundary_conditions%first
    sum_conn   = 0
    swap_order = PETSC_FALSE
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          cell_id  = cur_conn_set%conn(iconn)%GetIDDn()
          sum_conn = sum_conn + 1

          if (.not.geq_soil%mesh%is_active(cell_id )) cycle
          
          select case(cur_cond%itype)
          case(COND_DIRICHLET_FRM_OTR_GOVEQ, COND_DIRICHLET)

             internal_conn = PETSC_FALSE
             cond_type     = cur_cond%itype

             call RichardsFluxDerivativeWrtTemperature(        &
                  geq_soil%aux_vars_bc(sum_conn),              &
                  geq_soil%aux_vars_in(cell_id),               &
                  cur_conn_set%conn(iconn),                    &
                  compute_deriv,                               &
                  internal_conn,                               &
                  swap_order,                                  &
                  cond_type,                                   &
                  mflux,                                       &
                  dmflux_dT_up,                                &
                  dmflux_dT_dn                                 &
                  )

             call ThermalEnthalpyFlux(             &
                  geq_soil%aux_vars_bc(sum_conn ), &
                  geq_soil%aux_vars_in(cell_id  ), &
                  mflux,                           &
                  dmflux_dT_up,                    &
                  dmflux_dT_dn,                    &
                  cur_conn_set%conn(iconn),        &
                  compute_deriv,                   &
                  internal_conn,                   &
                  cond_type,                       &
                  flux,                            &
                  Jup,                             &
                  Jdn                              &
                  )

             row = cell_id - 1
             col = cell_id - 1
             val = Jdn
             call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          case (COND_HEAT_FLUX)             
             ! Do nothing
             
          case default
           write(iulog,*)'ThermalKSPTempSoilDivergence: Unknown boundary condition type'
           call endrun(msg=errMsg(__FILE__, __LINE__))

        end select

       enddo
       cur_cond => cur_cond%next
    enddo

    ! Source-sink cells
    cur_cond => geq_soil%source_sinks%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          cell_id = cur_conn_set%conn(iconn)%GetIDDn()

          if ((.not.geq_soil%mesh%is_active(cell_id))) cycle

          select case(cur_cond%itype)
          case(COND_HEAT_RATE)
             ! Do nothing

          case default
           write(iulog,*)'ThermalKSPTempSoilDivergence: Unknown source-sink condition type'
           call endrun(msg=errMsg(__FILE__, __LINE__))

        end select

     enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine ThermalEnthalpySoilDivergenceDeriv

  !------------------------------------------------------------------------
  subroutine OffDiagJacobian_Temperature(geq_soil, rank_of_other_goveq, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the derivative of energy residual equation w.r.t to pressure
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ, COND_DIRICHLET
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : COND_NULL
    use ThermalEnthalpyMod        , only : ThermalEnthalpyFlux
    use RichardsMod               , only : RichardsFlux, RichardsFluxDerivativeWrtTemperature
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: geq_soil
    PetscInt                                 :: rank_of_other_goveq
    Mat                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    type(condition_type),pointer             :: cur_cond
    class(connection_set_type), pointer      :: cur_conn_set
    PetscInt                                 :: ieqn
    PetscInt                                 :: iconn
    PetscInt                                 :: sum_conn
    PetscInt                                 :: cell_id_dn
    PetscInt                                 :: cell_id_up
    PetscInt                                 :: cell_id
    PetscInt                                 :: row
    PetscInt                                 :: col
    PetscReal                                :: flux
    PetscReal                                :: area
    PetscReal                                :: Jup
    PetscReal                                :: Jdn
    PetscBool                                :: compute_deriv
    PetscBool                                :: internal_conn
    PetscBool                                :: cur_cond_used
    PetscBool                                :: swap_order
    PetscInt                                 :: cond_type
    PetscReal                                :: val
    PetscReal                                :: mflux
    PetscReal                                :: eflux
    PetscReal                                :: dmflux_dT_up
    PetscReal                                :: dmflux_dT_dn

    compute_deriv = PETSC_TRUE

    ! Boundary cells
    cur_cond => geq_soil%boundary_conditions%first
    sum_conn   = 0
    swap_order = PETSC_FALSE
    do
       if (.not.associated(cur_cond)) exit

       cur_cond_used = PETSC_FALSE

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then

          do ieqn = 1, cur_cond%num_other_goveqs

             if (cur_cond%rank_of_other_goveqs(ieqn) == rank_of_other_goveq) then

                cur_cond_used = PETSC_TRUE

                cur_conn_set => cur_cond%conn_set

                do iconn = 1, cur_conn_set%num_connections

                   cell_id  = cur_conn_set%conn(iconn)%GetIDDn()
                   sum_conn = sum_conn + 1

                   if (.not.geq_soil%mesh%is_active(cell_id )) cycle

                   internal_conn = PETSC_FALSE
                   cond_type     = cur_cond%itype

                   call RichardsFluxDerivativeWrtTemperature(        &
                        geq_soil%aux_vars_bc(sum_conn),              &
                        geq_soil%aux_vars_in(cell_id),               &
                        cur_conn_set%conn(iconn),                    &
                        compute_deriv,                               &
                        internal_conn,                               &
                        swap_order,                                  &
                        cond_type,                                   &
                        mflux,                                       &
                        dmflux_dT_up,                                &
                        dmflux_dT_dn                                 &
                        )

                   call ThermalEnthalpyFlux(             &
                        geq_soil%aux_vars_bc(sum_conn ), &
                        geq_soil%aux_vars_in(cell_id  ), &
                        mflux,                           &
                        dmflux_dT_up,                    &
                        dmflux_dT_dn,                    &
                        cur_conn_set%conn(iconn),        &
                        compute_deriv,                   &
                        internal_conn,                   &
                        cond_type,                       &
                        flux,                            &
                        Jup,                             &
                        Jdn                              &
                        )

                   row = cell_id - 1
                   col = cur_conn_set%conn(iconn)%GetIDUp() - 1
                   val = Jup
                   call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

                enddo
             endif
          enddo

       endif

       if (.not. cur_cond_used) sum_conn = sum_conn + cur_cond%conn_set%num_connections

       cur_cond => cur_cond%next
    enddo

  end subroutine OffDiagJacobian_Temperature

    !------------------------------------------------------------------------
  subroutine OffDiagJacobian_Pressure(geq_soil, rank_of_other_goveq, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the derivative of energy residual equation w.r.t to pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)         :: geq_soil
    PetscInt                                         :: rank_of_other_goveq
    Mat                                              :: B
    PetscErrorCode                                   :: ierr
    !
    ! !LOCAL VARIABLES
    PetscBool                                :: coupling_via_BC

    call OffDiagJacobian_Pressure_ForBoundaryAuxVars( &
         geq_soil, rank_of_other_goveq, B, coupling_via_BC, ierr)

    if (coupling_via_BC) return

    call OffDiagJacobian_Pressure_ForInternalAuxVars( &
         geq_soil, rank_of_other_goveq, B, ierr)

  end subroutine OffDiagJacobian_Pressure

    !------------------------------------------------------------------------
  subroutine OffDiagJacobian_Pressure_ForBoundaryAuxVars(geq_soil, &
       rank_of_other_goveq, B, coupling_via_BC, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the derivative of energy residual equation w.r.t to pressure
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_NULL
    use ThermalEnthalpyMod        , only : ThermalEnthalpyFluxDerivativeWrtPressure
    use RichardsMod               , only : RichardsFlux
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use ThermalEnthalpyMod        , only : ThermalEnthalpyFlux
    use RichardsMod               , only : RichardsFlux
    use RichardsMod               , only : RichardsFlux_Internal
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)         :: geq_soil
    PetscInt                                         :: rank_of_other_goveq
    Mat                                              :: B
    PetscBool, intent (out)                          :: coupling_via_BC
    PetscErrorCode                                   :: ierr
    !
    ! !LOCAL VARIABLES
    type (therm_enthalpy_soil_auxvar_type) , pointer :: aux_vars(:)
    type(condition_type)                   , pointer :: cur_cond
    class(connection_set_type)             , pointer :: cur_conn_set
    PetscInt                                         :: iconn
    PetscInt                                         :: ieqn
    PetscInt                                         :: sum_conn
    PetscInt                                         :: cell_id_dn
    PetscInt                                         :: cell_id_up
    PetscInt                                         :: cell_id
    PetscInt                                         :: row
    PetscInt                                         :: col
    PetscInt                                         :: cond_type
    PetscReal                                        :: flux
    PetscReal                                        :: area
    PetscReal                                        :: Jup
    PetscReal                                        :: Jdn
    PetscBool                                        :: compute_deriv
    PetscBool                                        :: internal_conn
    PetscReal                                        :: val
    PetscReal                                        :: mflux
    PetscReal                                        :: eflux
    PetscReal                                        :: dmflux_dP_up
    PetscReal                                        :: dmflux_dP_dn
    PetscReal                                        :: dmflux_dT_up
    PetscReal                                        :: dmflux_dT_dn
    PetscReal                                        :: dtInv
    PetscReal                                        :: por
    PetscReal                                        :: den
    PetscReal                                        :: sat
    PetscReal                                        :: ul
    PetscReal                                        :: dpor_dP
    PetscReal                                        :: dden_dP
    PetscReal                                        :: dsat_dP
    PetscReal                                        :: dul_dP
    PetscReal                                        :: derivative
    PetscReal                                        :: den_soil
    PetscReal                                        :: heat_cap_soil
    PetscReal                                        :: temperature
    PetscBool                                        :: eqns_are_coupled
    PetscInt                                         :: ivar
    PetscBool                                        :: swap_order
    PetscBool                                        :: cur_cond_used

    coupling_via_BC  = PETSC_FALSE
    eqns_are_coupled = PETSC_FALSE
    compute_deriv    = PETSC_TRUE

    sum_conn = 0
    cur_cond => geq_soil%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       cur_cond_used = PETSC_FALSE

       if (cur_cond%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then

          do ieqn = 1, cur_cond%num_other_goveqs

             if (cur_cond%rank_of_other_goveqs(ieqn) == rank_of_other_goveq) then

                cur_cond_used = PETSC_TRUE

                if (.not.(cur_cond%is_the_other_GE_coupled_via_int_auxvars(ieqn))) &
                     coupling_via_BC = PETSC_TRUE

                cur_conn_set => cur_cond%conn_set

                do iconn = 1, cur_conn_set%num_connections
                   sum_conn = sum_conn + 1

                   cell_id = cur_conn_set%conn(iconn)%GetIDDn()

                   internal_conn = PETSC_FALSE
                   cond_type     = cur_cond%itype

                   swap_order = cur_cond%swap_order

                   if (.not.swap_order) then

                      call RichardsFlux(                    &
                           geq_soil%aux_vars_bc(sum_conn ), &
                           geq_soil%aux_vars_in(cell_id  ), &
                           cur_conn_set%conn(iconn),        &
                           compute_deriv,                   &
                           internal_conn,                   &
                           swap_order,                      &
                           cond_type,                       &
                           mflux,                           &
                           dmflux_dP_up,                    &
                           dmflux_dP_dn                     &
                           )

                      dmflux_dP_up = -dmflux_dP_up
                      dmflux_dP_dn = -dmflux_dP_dn

                      call ThermalEnthalpyFluxDerivativeWrtPressure( &
                           geq_soil%aux_vars_bc(sum_conn ),          &
                           geq_soil%aux_vars_in(cell_id  ),          &
                           mflux,                                    &
                           dmflux_dP_up,                             &
                           dmflux_dP_dn,                             &
                           cur_conn_set%conn(iconn),                 &
                           compute_deriv,                            &
                           internal_conn,                            &
                           cond_type,                                &
                           flux,                                     &
                           Jup,                                      &
                           Jdn                                       &
                           )

                      if (cur_cond%is_the_other_GE_coupled_via_int_auxvars(ieqn)) then
                         val = Jdn
                      else
                         val = Jup
                      endif

                   else

                      call RichardsFlux_Internal(           &
                           geq_soil%aux_vars_in(cell_id  ), &
                           geq_soil%aux_vars_bc(sum_conn ), &
                           cur_conn_set%conn(iconn),        &
                           compute_deriv,                   &
                           internal_conn,                   &
                           cond_type,                       &
                           mflux,                           &
                           dmflux_dP_up,                    &
                           dmflux_dP_dn                     &
                           )

                      dmflux_dP_up = -dmflux_dP_up
                      dmflux_dP_dn = -dmflux_dP_dn

                      call ThermalEnthalpyFluxDerivativeWrtPressure( &
                           geq_soil%aux_vars_in(cell_id  ),          &
                           geq_soil%aux_vars_bc(sum_conn ),          &
                           mflux,                                    &
                           dmflux_dP_up,                             &
                           dmflux_dP_dn,                             &
                           cur_conn_set%conn(iconn),                 &
                           compute_deriv,                            &
                           internal_conn,                            &
                           cond_type,                                &
                           flux,                                     &
                           Jup,                                      &
                           Jdn                                       &
                           )

                      if (cur_cond%is_the_other_GE_coupled_via_int_auxvars(ieqn)) then
                         val = -Jup
                      else
                         val = -Jdn
                      endif
                   endif


                   row = cell_id - 1
                   col = cur_conn_set%conn(iconn)%GetIDUp() - 1

                   if (cur_cond%is_the_other_GE_coupled_via_int_auxvars(ieqn)) then
                      col = row
                   endif

                   call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

                enddo

             endif
          enddo

       endif

       if (.not. cur_cond_used) sum_conn = sum_conn + cur_cond%conn_set%num_connections

       cur_cond => cur_cond%next
    enddo
    coupling_via_BC  = PETSC_FALSE


  end subroutine OffDiagJacobian_Pressure_ForBoundaryAuxVars

    !------------------------------------------------------------------------
  subroutine OffDiagJacobian_Pressure_ForInternalAuxVars(geq_soil, &
       rank_of_other_goveq, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the derivative of energy residual equation w.r.t to pressure
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_NULL
    use ThermalEnthalpyMod        , only : ThermalEnthalpyFluxDerivativeWrtPressure
    use RichardsMod               , only : RichardsFlux
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use ThermalEnthalpyMod        , only : ThermalEnthalpyFlux
    use RichardsMod               , only : RichardsFlux
    use RichardsMod               , only : RichardsFlux_Internal
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)         :: geq_soil
    PetscInt                                         :: rank_of_other_goveq
    Mat                                              :: B
    PetscErrorCode                                   :: ierr
    !
    ! !LOCAL VARIABLES
    type (therm_enthalpy_soil_auxvar_type) , pointer :: aux_vars(:)
    type(condition_type)                   , pointer :: cur_cond
    class(connection_set_type)             , pointer :: cur_conn_set
    PetscInt                                         :: iconn
    PetscInt                                         :: ieqn
    PetscInt                                         :: sum_conn
    PetscInt                                         :: cell_id_dn
    PetscInt                                         :: cell_id_up
    PetscInt                                         :: cell_id
    PetscInt                                         :: row
    PetscInt                                         :: col
    PetscInt                                         :: cond_type
    PetscReal                                        :: flux
    PetscReal                                        :: area
    PetscReal                                        :: Jup
    PetscReal                                        :: Jdn
    PetscBool                                        :: compute_deriv
    PetscBool                                        :: internal_conn
    PetscReal                                        :: val
    PetscReal                                        :: mflux
    PetscReal                                        :: eflux
    PetscReal                                        :: dmflux_dP_up
    PetscReal                                        :: dmflux_dP_dn
    PetscReal                                        :: dmflux_dT_up
    PetscReal                                        :: dmflux_dT_dn
    PetscReal                                        :: dtInv
    PetscReal                                        :: por
    PetscReal                                        :: den
    PetscReal                                        :: sat
    PetscReal                                        :: ul
    PetscReal                                        :: dpor_dP
    PetscReal                                        :: dden_dP
    PetscReal                                        :: dsat_dP
    PetscReal                                        :: dul_dP
    PetscReal                                        :: derivative
    PetscReal                                        :: den_soil
    PetscReal                                        :: heat_cap_soil
    PetscReal                                        :: temperature
    PetscInt                                         :: ivar
    PetscBool                                        :: swap_order
    PetscBool                                        :: cur_cond_used

    compute_deriv    = PETSC_TRUE

    aux_vars => geq_soil%aux_vars_in

    dtInv = 1.d0/geq_soil%dtime

    ! Interior cells
    do cell_id = 1, geq_soil%mesh%ncells_local

       if (geq_soil%mesh%is_active(cell_id)) then

          por           = aux_vars(cell_id)%por
          den_soil      = aux_vars(cell_id)%den_soil
          heat_cap_soil = aux_vars(cell_id)%heat_cap_soil

          den          = aux_vars(cell_id)%den
          sat          = aux_vars(cell_id)%sat
          ul            = aux_vars(cell_id)%ul

          dpor_dP       = aux_vars(cell_id)%dpor_dP
          dden_dP      = aux_vars(cell_id)%dden_dP
          dsat_dP      = aux_vars(cell_id)%dsat_dP
          dul_dP        = aux_vars(cell_id)%dul_dP

          temperature   = aux_vars(cell_id)%temperature

          derivative =                                 &
               dpor_dP * den     * sat     * ul     + &
               por     * dden_dP * sat     * ul     + &
               por     * den     * dsat_dP * ul     + &
               por     * den     * sat     * dul_dP + &
               (-dpor_dP * den_soil * heat_cap_soil *(temperature - 273.15d0))

          derivative = derivative * geq_soil%mesh%vol(cell_id) * dtInv
       else
          derivative = 1.d0
       endif

       row = cell_id - 1
       col = cell_id - 1

       call MatSetValuesLocal(B, 1, row, 1, col, derivative, ADD_VALUES, ierr); CHKERRQ(ierr)
    enddo

    ! Interior cells
    cur_conn_set => geq_soil%mesh%intrn_conn_set_list%first
    sum_conn   = 0
    swap_order = PETSC_FALSE
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(iconn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()

          internal_conn = PETSC_TRUE
          cond_type     = COND_NULL

          if ((.not.geq_soil%mesh%is_active(cell_id_up)) .or. &
               (.not.geq_soil%mesh%is_active(cell_id_dn))) cycle

          call RichardsFlux(                     &
               geq_soil%aux_vars_in(cell_id_up), &
               geq_soil%aux_vars_in(cell_id_dn), &
               cur_conn_set%conn(iconn),         &
               compute_deriv,                    &
               internal_conn,                    &
               swap_order,                       &
               cond_type,                        &
               mflux,                            &
               dmflux_dP_up,                     &
               dmflux_dP_dn                      &
               )

          dmflux_dP_up = -dmflux_dP_up
          dmflux_dP_dn = -dmflux_dP_dn

          call ThermalEnthalpyFluxDerivativeWrtPressure( &
               geq_soil%aux_vars_in(cell_id_up),         &
               geq_soil%aux_vars_in(cell_id_dn),         &
               mflux,                                    &
               dmflux_dP_up,                             &
               dmflux_dP_dn,                             &
               cur_conn_set%conn(iconn),                 &
               compute_deriv,                            &
               internal_conn,                            &
               cond_type,                                &
               flux,                                     &
               Jup,                                      &
               Jdn                                       &
               )

          row = cell_id_up - 1
          col = cell_id_up - 1
          val = -Jup
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_up - 1
          col = cell_id_dn - 1
          val = -Jdn
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_dn - 1
          col = cell_id_up - 1
          val = Jup
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = cell_id_dn - 1
          col = cell_id_dn - 1
          val = Jdn
          call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    sum_conn = 0
    cur_cond => geq_soil%boundary_conditions%first
    internal_conn = PETSC_FALSE
    do
       if (.not.associated(cur_cond)) exit

       cur_cond_used = PETSC_FALSE

       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then

          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             cell_id = cur_conn_set%conn(iconn)%GetIDDn()

             cond_type     = cur_cond%itype

             swap_order = cur_cond%swap_order

             if (.not.swap_order) then

                call RichardsFlux(                    &
                     geq_soil%aux_vars_bc(sum_conn ), &
                     geq_soil%aux_vars_in(cell_id  ), &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     swap_order,                      &
                     cond_type,                       &
                     mflux,                           &
                     dmflux_dP_up,                    &
                     dmflux_dP_dn                     &
                     )

                dmflux_dP_up = -dmflux_dP_up
                dmflux_dP_dn = -dmflux_dP_dn

                call ThermalEnthalpyFluxDerivativeWrtPressure( &
                     geq_soil%aux_vars_bc(sum_conn ),          &
                     geq_soil%aux_vars_in(cell_id  ),          &
                     mflux,                                    &
                     dmflux_dP_up,                             &
                     dmflux_dP_dn,                             &
                     cur_conn_set%conn(iconn),                 &
                     compute_deriv,                            &
                     internal_conn,                            &
                     cond_type,                                &
                     flux,                                     &
                     Jup,                                      &
                     Jdn                                       &
                     )

                val = Jdn

             else

                call RichardsFlux_Internal(           &
                     geq_soil%aux_vars_in(cell_id  ), &
                     geq_soil%aux_vars_bc(sum_conn ), &
                     cur_conn_set%conn(iconn),        &
                     compute_deriv,                   &
                     internal_conn,                   &
                     cond_type,                       &
                     mflux,                           &
                     dmflux_dP_up,                    &
                     dmflux_dP_dn                     &
                     )

                dmflux_dP_up = -dmflux_dP_up
                dmflux_dP_dn = -dmflux_dP_dn

                call ThermalEnthalpyFluxDerivativeWrtPressure( &
                     geq_soil%aux_vars_in(cell_id  ),          &
                     geq_soil%aux_vars_bc(sum_conn ),          &
                     mflux,                                    &
                     dmflux_dP_up,                             &
                     dmflux_dP_dn,                             &
                     cur_conn_set%conn(iconn),                 &
                     compute_deriv,                            &
                     internal_conn,                            &
                     cond_type,                                &
                     flux,                                     &
                     Jup,                                      &
                     Jdn                                       &
                     )

                val = -Jup
             endif

             row = cell_id - 1
             col = row

             call MatSetValuesLocal(B, 1, row, 1, col, val, ADD_VALUES, ierr); CHKERRQ(ierr)

          enddo

       else
          sum_conn = sum_conn + cur_cond%conn_set%num_connections
       endif

       cur_cond => cur_cond%next
    enddo

  end subroutine OffDiagJacobian_Pressure_ForInternalAuxVars

  !------------------------------------------------------------------------

  subroutine DiffHeatFlux(T_up, therm_cond_up, T_dn, therm_cond_dn, &
                          dist_up, dist_dn, flux)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: T_up
    PetscReal, intent(in)  :: therm_cond_up
    PetscReal, intent(in)  :: T_dn
    PetscReal, intent(in)  :: therm_cond_dn
    PetscReal, intent(in)  :: dist_up
    PetscReal, intent(in)  :: dist_dn
    PetscReal, intent(out) :: flux
    !
    ! !LOCAL VARIABLES
    PetscReal :: therm_cond

    ! Distance weighted harmonic average
    therm_cond = therm_cond_up*therm_cond_dn*(dist_up+dist_dn)/ &
                 (therm_cond_up*dist_dn + therm_cond_dn*dist_up)

    flux = -therm_cond*(T_up - T_dn)/(dist_up + dist_dn)

  end subroutine DiffHeatFlux

  !------------------------------------------------------------------------

  subroutine ThermEnthalpySoilPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Performs pre-solve computation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this

    ! !LOCAL VARIABLES
    PetscErrorCode                           :: ierr
    PetscReal, dimension(:), pointer         :: f_p

    ! Computes contribution to residual equation (accumulation term) based on
    ! the pressure values from previous time step.

    call VecGetArrayF90(this%accum_prev, f_p, ierr); CHKERRQ(ierr)
    f_p(:) = 0.d0
    call ThermalEnthalpySoilAccum(this, f_p)
    call VecRestoreArrayF90(this%accum_prev, f_p, ierr); CHKERRQ(ierr)

  end subroutine ThermEnthalpySoilPreSolve

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySoilCreateVectors(this)
    !
    ! !DESCRIPTION:
    ! Creates required PETSc vectors.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type) :: this
    !
    PetscErrorCode                           :: ierr

    call VecCreateSeq(PETSC_COMM_SELF, this%mesh%ncells_local, &
            this%accum_prev, ierr)
    CHKERRQ(ierr)

  end subroutine ThermEnthalpySoilCreateVectors

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetSoilPermeability(this, perm_x, perm_y, perm_z)
    !
    ! !DESCRIPTION:
    ! Set soil permeability values
    !
    ! !USES:
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxVarSetAbsPerm
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscReal                                 , pointer, intent(in) :: perm_x(:)
    PetscReal                                 , pointer, intent(in) :: perm_y(:)
    PetscReal                                 , pointer, intent(in) :: perm_z(:)

    if (this%mesh%ncells_local /= size(perm_x)) then
       write(iulog,*) 'No. of values for soil permeability is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(perm_x)
       write(iulog,*) 'No. of grid cells           = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxVarSetAbsPerm(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
           this%boundary_conditions, this%source_sinks, perm_x, perm_y, perm_z)

  end subroutine ThermEnthalpySetSoilPermeability

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetRelativePermeability(this, relperm_type, param_1, param_2)
    !
    ! !DESCRIPTION:
    ! Set soil permeability values
    !
    ! !USES:
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxVarSetRelPerm
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscInt                                  , pointer, intent(in) :: relperm_type(:)
    PetscReal                                 , pointer, intent(in) :: param_1(:)
    PetscReal                                 , pointer, intent(in) :: param_2(:)

    if (this%mesh%ncells_local /= size(relperm_type)) then
       write(iulog,*) 'No. of values for relative permeability is not equal to no. of grid cells.'
       write(iulog,*) 'No. of rel. perm. values = ',size(relperm_type)
       write(iulog,*) 'No. of grid cells        = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxVarSetRelPerm(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
           this%boundary_conditions, this%source_sinks, relperm_type, param_1, param_2)

  end subroutine ThermEnthalpySetRelativePermeability

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetSoilPorosity(this, por)
    !
    ! !DESCRIPTION:
    ! Set soil porosity value
    !
    ! !USES:
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxVarSetPorosity
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscReal                                 , pointer, intent(in) :: por(:)

    if (this%mesh%ncells_local /= size(por)) then
       write(iulog,*) 'No. of values for soil porosity is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(por)
       write(iulog,*) 'No. of grid cells           = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxVarSetPorosity(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
              this%boundary_conditions, this%source_sinks, por)

  end subroutine ThermEnthalpySetSoilPorosity

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetSaturationFunction(this, satfunc_type, &
       alpha, lambda, sat_res)
    !
    ! !DESCRIPTION:
    ! Set saturation function
    !
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxVarSetSatFunc
    ! !USES:
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscInt                                  , pointer, intent(in) :: satfunc_type(:)
    PetscReal                                 , pointer, intent(in) :: alpha(:)
    PetscReal                                 , pointer, intent(in) :: lambda(:)
    PetscReal                                 , pointer, intent(in) :: sat_res(:)

    if (this%mesh%ncells_local /= size(alpha)) then
       write(iulog,*) 'No. of values for saturation function is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil sat. func. values = ',size(alpha)
       write(iulog,*) 'No. of grid cells             = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxVarSetSatFunc(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
           this%boundary_conditions, this%source_sinks, satfunc_type, alpha, lambda, sat_res)

  end subroutine ThermEnthalpySetSaturationFunction

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetHeatCapacity(this, data)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxSetHeatCap
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscReal                                 , pointer, intent(in) :: data(:)

    if (this%mesh%ncells_local /= size(data)) then
       write(iulog,*) 'No. of values for heat capacity is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(data)
       write(iulog,*) 'No. of grid cells           = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxSetHeatCap(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
           this%boundary_conditions, this%source_sinks, data)

  end subroutine ThermEnthalpySetHeatCapacity

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetThermalCondWet(this, data)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxSetThermalCondWet
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscReal                                 , pointer, intent(in) :: data(:)

    if (this%mesh%ncells_local /= size(data)) then
       write(iulog,*) 'No. of values for wet thermal conductivity is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(data)
       write(iulog,*) 'No. of grid cells           = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxSetThermalCondWet(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
           this%boundary_conditions, this%source_sinks, data)

  end subroutine ThermEnthalpySetThermalCondWet

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetThermalCondDry(this, data)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxSetThermalCondDry
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscReal                                 , pointer, intent(in) :: data(:)

    if (this%mesh%ncells_local /= size(data)) then
       write(iulog,*) 'No. of values for dry thermal conductivity is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(data)
       write(iulog,*) 'No. of grid cells           = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxSetThermalCondDry(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
           this%boundary_conditions, this%source_sinks, data)

  end subroutine ThermEnthalpySetThermalCondDry

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetThermalAlpha(this, data)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxSetThermalAlpha
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscReal                                 , pointer, intent(in) :: data(:)

    if (this%mesh%ncells_local /= size(data)) then
       write(iulog,*) 'No. of values for thermal alpha is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(data)
       write(iulog,*) 'No. of grid cells           = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxSetThermalAlpha(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
           this%boundary_conditions, this%source_sinks, data)

  end subroutine ThermEnthalpySetThermalAlpha

  !------------------------------------------------------------------------
  subroutine ThermEnthalpySetSoilDensity(this, data)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ThermalEnthalpySoilAuxMod  , only : ThermalEnthalpySoilAuxSetSoilDensity
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    class(goveqn_thermal_enthalpy_soil_type)  , intent(inout)       :: this
    PetscReal                                 , pointer, intent(in) :: data(:)

    if (this%mesh%ncells_local /= size(data)) then
       write(iulog,*) 'No. of values for soil density is not equal to no. of grid cells.'
       write(iulog,*) 'No. of soil porosity values = ',size(data)
       write(iulog,*) 'No. of grid cells           = ',this%mesh%ncells_local
    endif

    call ThermalEnthalpySoilAuxSetSoilDensity(this%aux_vars_in, this%aux_vars_bc, this%aux_vars_ss, &
           this%boundary_conditions, this%source_sinks, data)

  end subroutine ThermEnthalpySetSoilDensity

#endif
  
end module GoveqnThermalEnthalpySoilType
