module GoveqnCanopyAirVaporType

#ifdef USE_PETSC_LIB
!#define USE_BONAN_FORMULATION
  
#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType , only : goveqn_base_type
  use CanopyAirVaporAuxType     , only : cair_vapor_auxvar_type
  use CanopyAirVaporConnAuxType , only : cair_vapor_conn_auxvar_type
  use petscvec
  use petscmat
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_cair_vapor_type

     type(cair_vapor_auxvar_type), pointer ::  aux_vars_in(:)
     type(cair_vapor_auxvar_type), pointer ::  aux_vars_bc(:)

     type(cair_vapor_conn_auxvar_type) , pointer :: aux_vars_conn_in(:)
     type(cair_vapor_conn_auxvar_type) , pointer :: aux_vars_conn_bc(:)

     PetscInt                                   :: nLeaf, nleafGE
     PetscInt                         , pointer :: LeafGE2CAirAuxVar_in(:,:)

   contains

     procedure, public :: Setup                     => CAirVaporSetup
     procedure, public :: AllocateAuxVars           => CAirVaporAllocateAuxVars
     procedure, public :: GetFromSoeAuxVarsCturb    => CAirVaporGetFromSoeAuxVarsCturb
     procedure, public :: SavePrimaryIndependentVar => CAirVaporSavePrmIndepVar
     procedure, public :: GetRValues                => CAirVaporGetRValues
     procedure, public :: SetRValues                => CAirVaporSetRValues
     procedure, public :: ComputeRHS                => CAirVaporComputeRHS
     procedure, public :: ComputeOperatorsDiag      => CAirVaporComputeOperatorsDiag
     procedure, public :: ComputeOperatorsOffDiag   => CAirVaporComputeOperatorsOffDiag

     procedure, public :: SetLeaf2CAirMap
     procedure, public :: SetDefaultLeaf2CAirMap

  end type goveqn_cair_vapor_type

contains

  !------------------------------------------------------------------------
  subroutine CAirVaporSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_VAPOR
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this

    call this%Create()

    this%name  = ""
    this%itype = GE_CANOPY_AIR_VAPOR

    nullify(this%aux_vars_in)
    nullify(this%aux_vars_bc)
    nullify(this%aux_vars_conn_in)
    nullify(this%aux_vars_conn_bc)

    this%nLeaf = 0
    this%nleafGE = 0
    nullify(this%LeafGE2CAirAuxVar_in)

  end subroutine CAirVaporSetup

  !------------------------------------------------------------------------
  subroutine SetLeaf2CAirMap(this, CAirIdForLeaf, nCAirIdForLeaf)

    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    PetscInt, pointer :: CAirIdForLeaf(:)
    PetscInt          :: nCAirIdForLeaf
    !
    PetscInt :: ileaf, icair
    PetscInt, pointer :: nleaf(:)

    allocate(this%LeafGE2CAirAuxVar_in(nCAirIdForLeaf,2))
    allocate(nLeaf(this%mesh%ncells_all))

    this%nleaf = nCAirIdForLeaf
    nleaf(:) = 0
    do ileaf = 1, nCAirIdForLeaf
       icair = CAirIdForLeaf(ileaf)
       nLeaf(icair) = nLeaf(icair) + 1

       this%LeafGE2CAirAuxVar_in(ileaf, 1) = icair       ! Canopy airspace id for the ileaf-th leaf
       this%LeafGE2CAirAuxVar_in(ileaf, 2) = nLeaf(icair)! The id of leaf within the canopy airspace
    enddo

  end subroutine SetLeaf2CAirMap

  !------------------------------------------------------------------------
  subroutine SetDefaultLeaf2CAirMap(this)

    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    !
    PetscInt :: ileaf, icair

    allocate(this%LeafGE2CAirAuxVar_in(this%mesh%ncells_all,2))

    this%nleaf = this%mesh%ncells_all
    do ileaf = 1, this%mesh%ncells_all
       icair = ileaf

       this%LeafGE2CAirAuxVar_in(icair, 1) = icair  ! Canopy airspace id for the ileaf-th leaf
       this%LeafGE2CAirAuxVar_in(icair, 2) = 1      ! The id of leaf within the canopy airspace
    enddo

  end subroutine SetDefaultLeaf2CAirMap

  !------------------------------------------------------------------------
  subroutine CAirVaporAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates memory for storing auxiliary variables associated with:
    !   + Internal control volumes,
    !   + Boundary condtions,
    !   + Source-sink condition.
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use CouplingVariableType      , only : coupling_variable_type
    use MultiPhysicsProbConstants , only : VAR_LEAF_TEMPERATURE
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    !
    type(condition_type)         , pointer :: cur_cond
    type(coupling_variable_type) , pointer :: cpl_var
    class(connection_set_type)   , pointer :: cur_conn_set
    PetscInt                               :: ncells_cond
    PetscInt                               :: ncond
    PetscInt                               :: icond, icell, icair, ileaf
    PetscInt                     , pointer :: nleaf(:)
    PetscInt                               :: sum_conn

    ! Allocate memory and initialize aux vars: For internal connections
    allocate(this%aux_vars_in(this%mesh%ncells_all))
    allocate(nLeaf(this%mesh%ncells_all))

    ! Determine the number of leaves that are coupled with the canopy
    ! air space
    this%nleafGE = 0
    cpl_var => this%coupling_vars%first
    do
       if (.not.associated(cpl_var)) exit
       if ( cpl_var%variable_type == VAR_LEAF_TEMPERATURE) then
          this%nleafGE = this%nleafGE + 1
       end if
       cpl_var => cpl_var%next
    end do

    if (this%nLeaf == 0) then
       call endrun(msg="Leaf2CAirMap needs be set before allocating auxvars "//errmsg(__FILE__, __LINE__))
    endif

    nLeaf(:) = 0
    do ileaf = 1, this%nLeaf
       icair = this%LeafGE2CAirAuxVar_in(ileaf,1)
       if (icair > this%mesh%ncells_local) then
          call endrun(msg="Leaf2CAirMap is incorrect as the canopy airspace id is greater than number of grid cells "//errmsg(__FILE__, __LINE__))
       endif
       nLeaf(icair) = nLeaf(icair) + 1
    enddo

    do icell = 1,this%mesh%ncells_all
       call this%aux_vars_in(icell)%Init(nleaf(icell) * this%nleafGE)
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
       call this%aux_vars_bc(icond)%Init(0)
    enddo

    ! Find number of internal connections
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    sum_conn = 0
    do
       if (.not.associated(cur_conn_set)) exit
       sum_conn = sum_conn + cur_conn_set%num_connections
       cur_conn_set => cur_conn_set%next
    enddo
    allocate(this%aux_vars_conn_in(sum_conn))
    do icond = 1,sum_conn
       call this%aux_vars_conn_in(icond)%Init()
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
    allocate(this%aux_vars_conn_bc (ncells_cond))
    do icond = 1,ncells_cond
       call this%aux_vars_conn_bc(icond)%Init()
    enddo

  end subroutine CAirVaporAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine CAirVaporGetFromSoeAuxVarsCturb(this, cturb)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use CanopyTurbulenceAuxType, only : canopy_turbulence_auxvar_type 
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type)        :: this
    class(canopy_turbulence_auxvar_type) :: cturb
    !
    PetscInt :: ncair, icair, icell, iconn, k

    ! all layers
    do icair = 1, cturb%ncair
       do k = 1, cturb%ncan_lev
          icell = (icair-1)*cturb%ncan_lev + k
          this%aux_vars_in(icell)%cpair    = cturb%cpair(icair)
          this%aux_vars_in(icell)%rhomol   = cturb%rhomol(icair)
          this%aux_vars_in(icell)%pref     = cturb%pref(icair)

          if (this%aux_vars_in(icell)%is_soil) then
             this%aux_vars_in(icell)%soil_rhg         = cturb%rhgsoi(icair)
             this%aux_vars_in(icell)%soil_rn          = cturb%rnsoi(icair)
             this%aux_vars_in(icell)%soil_tk          = cturb%tksoi(icair)
             this%aux_vars_in(icell)%soil_dz          = cturb%dzsoi(icair)
             this%aux_vars_in(icell)%soil_resis       = cturb%ressoi(icair)
             this%aux_vars_in(icell)%soil_temperature = cturb%tsoi(icair)
          endif
       end do
    end do

    ! internal connections
    do icair = 1, cturb%ncair
       do k = 1, cturb%ncan_lev-1
          iconn = (icair-1)*(cturb%ncan_lev-1) + k
          this%aux_vars_conn_in(iconn)%ga    = cturb%ga_prof(icair,k)
       end do

    end do

    ! top-layer
    do icair = 1, cturb%ncair
       iconn = icair
       this%aux_vars_conn_bc(iconn)%ga     = cturb%ga_prof(icair, k)
       this%aux_vars_bc(iconn)%water_vapor = cturb%vref(icair)
       this%aux_vars_bc(iconn)%pref        = cturb%pref(icair)
    end do

  end subroutine CAirVaporGetFromSoeAuxVarsCturb

  !------------------------------------------------------------------------
  subroutine CAirVaporSavePrmIndepVar(this, x)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use CanopyTurbulenceAuxType, only : canopy_turbulence_auxvar_type 
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    Vec                           :: x
    !
    PetscScalar, pointer :: x_p(:)
    PetscInt             :: ghosted_id, size
    PetscErrorCode       :: ierr
    
    call VecGetLocalSize(x, size, ierr); CHKERRQ(ierr)

    if (size /= this%mesh%ncells_local) then
       call endrun(msg="ERROR size of vector /= number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    end if

    call VecGetArrayF90(x, x_p, ierr); CHKERRQ(ierr)

    do ghosted_id = 1, this%mesh%ncells_local
       this%aux_vars_in(ghosted_id)%water_vapor = x_p(ghosted_id)
    end do

    call VecRestoreArrayF90(x, x_p, ierr); CHKERRQ(ierr)

  end subroutine CAirVaporSavePrmIndepVar

  !------------------------------------------------------------------------
  subroutine CAirVaporGetRValues (this, auxvar_type, var_type, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    PetscInt                      :: auxvar_type
    PetscInt                      :: var_type
    PetscInt                      :: nauxvar
    PetscReal, pointer            :: var_values(:)

    if (nauxvar > this%mesh%ncells_all) then
      call endrun(msg="ERROR nauxvar exceeds the number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    endif

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call CAirVaporGetRValuesFromAuxVars(this%aux_vars_in, var_type, nauxvar, var_values)
    case (AUXVAR_BC)
       call CAirVaporGetRValuesFromAuxVars(this%aux_vars_bc, var_type, nauxvar, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CAirVaporGetRValues

  !------------------------------------------------------------------------
  subroutine CAirVaporSetRValues (this, auxvar_type, var_type, geq_leaf_temp_rank, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    PetscInt                      :: auxvar_type
    PetscInt                      :: var_type
    PetscInt                      :: geq_leaf_temp_rank
    PetscInt                      :: nauxvar
    PetscReal, pointer            :: var_values(:)

    if (nauxvar > this%mesh%ncells_all*this%nleafGE) then
      call endrun(msg="ERROR nauxvar exceeds the number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    endif

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call CAirVaporSetRValuesFromAuxVars(this%aux_vars_in, var_type, geq_leaf_temp_rank, nauxvar, var_values, this%LeafGE2CAirAuxVar_in, this%nleafGE)
    case (AUXVAR_BC)
       call CAirVaporSetRValuesFromAuxVars(this%aux_vars_bc, var_type, geq_leaf_temp_rank, nauxvar, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CAirVaporSetRValues

  !------------------------------------------------------------------------
  subroutine CAirVaporGetRValuesFromAuxVars (aux_var, var_type, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_WATER_VAPOR
    !
    implicit none
    !
    ! !ARGUMENTS
    type(cair_vapor_auxvar_type) , pointer :: aux_var(:)
    PetscInt                               :: var_type
    PetscInt                               :: nauxvar
    PetscReal                   , pointer  :: var_values(:)
    !
    PetscInt :: iauxvar

    select case(var_type)
    case(VAR_WATER_VAPOR)
       do iauxvar = 1,nauxvar
          var_values(iauxvar) = aux_var(iauxvar)%water_vapor
       end do
    case default
       write(iulog,*) 'CAirVaporGetRValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CAirVaporGetRValuesFromAuxVars
       
  !------------------------------------------------------------------------
  subroutine CAirVaporSetRValuesFromAuxVars (aux_var, var_type, geq_leaf_temp_rank, nauxvar, var_values, LeafGE2CAirAuxVar, nleafGE)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_LEAF_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(cair_vapor_auxvar_type) , pointer :: aux_var(:)
    PetscInt                               :: var_type
    PetscInt                               :: geq_leaf_temp_rank
    PetscInt                               :: nauxvar
    PetscReal                    , pointer :: var_values(:)
    PetscInt          , optional , pointer :: LeafGE2CAirAuxVar(:,:)
    PetscInt          , optional           :: nleafGE
    !
    PetscInt :: iauxvar,  leaf_idx, cair_auxvar_idx

    select case(var_type)
    case(VAR_LEAF_TEMPERATURE)
       if (.not.present(LeafGE2CAirAuxVar) .or. .not.present(nleafGE)) then
          write(iulog,*) 'Optional arguments LeafGE2CAirAuxVar and nleafGE needed for setting leaf temperature'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
       do iauxvar = 1,nauxvar
          cair_auxvar_idx = LeafGE2CAirAuxVar(iauxvar,1)
          leaf_idx        = LeafGE2CAirAuxVar(iauxvar,2) + &
                            (geq_leaf_temp_rank-1)*aux_var(cair_auxvar_idx)%nleaf/nleafGE
          aux_var(cair_auxvar_idx)%leaf_temperature(leaf_idx) = var_values(iauxvar)
       end do
    case(VAR_TEMPERATURE)
       do iauxvar = 1,nauxvar
          aux_var(iauxvar)%temperature = var_values(iauxvar)
       end do
    case default
       write(iulog,*) 'CAirVaporGetSValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
       
  end subroutine CAirVaporSetRValuesFromAuxVars

  !------------------------------------------------------------------------
  subroutine CAirVaporComputeRHS(this, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    !
    use WaterVaporMod, only : SatVap
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    Vec                           :: B
    PetscErrorCode                :: ierr
    !
    PetscScalar, pointer          :: b_p(:)

    call VecGetArrayF90(B, b_p, ierr)

    call CAirVaporComputeRhsAccumulation(this, b_p)
    call CAirVaporRhsDivergence(this, b_p)

    call VecRestoreArrayF90(B, b_p, ierr)

  end subroutine CAirVaporComputeRHS

  !------------------------------------------------------------------------
  function SoilAirIConn(this, icell)
    !
    use ConnectionSetType, only : connection_set_type
    !
    implicit none
    !
    class(goveqn_cair_vapor_type)         :: this
    PetscInt                             :: icell, SoilAirIConn
    !
    class(connection_set_type) , pointer :: cur_conn_set
    PetscInt                             :: iconn
    PetscBool                            :: found

    cur_conn_set => this%mesh%intrn_conn_set_list%first

    ! Find the internal connection that is between the soil layer overlying
    ! air space
    found = PETSC_FALSE
    do iconn = 1, cur_conn_set%num_connections
       if (cur_conn_set%conn(iconn)%GetIDUp() == icell .or. &
            cur_conn_set%conn(iconn)%GetIDDn() == icell) then
          found = PETSC_TRUE
          exit
       end if
    end do
    if (.not. found) then
       write(iulog,*)'Internal connection not found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    SoilAirIConn = iconn

  end function SoilAirIConn

  !------------------------------------------------------------------------
  function gs0(this, icell, iconn)
    !
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    class(goveqn_cair_vapor_type) :: this
    PetscInt                               :: icell, iconn
    PetscReal                              :: gs0
    !
    class(cair_vapor_auxvar_type) , pointer :: auxvar(:)
    PetscReal                              :: qsat, s0
    PetscReal                              :: gsw, gsa

    auxvar => this%aux_vars_in

    call SatVap(auxvar(icell)%temperature, qsat, s0)
    qsat = qsat/this%aux_vars_in(icell)%pref
    s0   = s0  /this%aux_vars_in(icell)%pref

    gsw = 1.d0 / auxvar(icell)%soil_resis * auxvar(icell)%rhomol
    gsa = this%aux_vars_conn_in(iconn)%ga
    gs0 = gsw * gsa / (gsw + gsa)

  end function gs0

  !------------------------------------------------------------------------
  function alpha0(this, icell, iconn)
    !
    ! Equation 16.86 from Bonan (2019)
    !
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    class(goveqn_cair_vapor_type)           :: this
    PetscInt                               :: icell, iconn
    !
    PetscReal                              :: numer, denom, alpha0

    numer = this%aux_vars_in(icell)%cpair * this%aux_vars_conn_in(iconn)%ga
    denom = gamma0(this, icell, iconn)

    alpha0 = numer/denom

  end function alpha0

  !------------------------------------------------------------------------
  function beta0(this, icell, iconn)
    !
    ! Equation 16.87 from Bonan (2019)
    !
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    class(goveqn_cair_vapor_type)           :: this
    PetscInt                               :: icell, iconn
    PetscReal                              :: beta0
    !
    PetscReal                              :: qsat, s0
    PetscReal                              :: gs_0, lambda, numer, denom

    lambda = HVAP * MM_H2O

    call SatVap(this%aux_vars_in(icell)%temperature, qsat, s0)
    qsat = qsat/this%aux_vars_in(icell)%pref
    s0   = s0  /this%aux_vars_in(icell)%pref

    gs_0 = gs0(this, icell, iconn)

    numer = lambda * gs_0
    denom = gamma0(this, icell, iconn)

    beta0 = numer/denom

  end function beta0

  !------------------------------------------------------------------------
  function delta0(this, icell, iconn)
    !
    ! Equation 16.88 from Bonan (2019)
    !
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    class(goveqn_cair_vapor_type)           :: this
    PetscInt                                :: icell, iconn
    PetscReal                               :: delta0
    !
    class(cair_vapor_auxvar_type) , pointer :: auxvar(:)
    PetscReal                               :: qsat, s0
    PetscReal                               :: gs_0, lambda, numer, denom

    auxvar => this%aux_vars_in
    lambda = HVAP * MM_H2O

    call SatVap(this%aux_vars_in(icell)%temperature, qsat, s0)
    qsat = qsat/this%aux_vars_in(icell)%pref
    s0   = s0  /this%aux_vars_in(icell)%pref

    gs_0 = gs0(this, icell, iconn)

    numer = &
         ( auxvar(icell)%soil_rn                                                         & ! (Rn_soil
         - lambda * auxvar(icell)%soil_rhg * gs_0 * (qsat - s0*auxvar(icell)%temperature) & !  lambda * h_{s0} * g_{s0} * (qsat(T0) - s0*T0)
         + auxvar(icell)%soil_tk/auxvar(icell)%soil_dz * auxvar(icell)%soil_temperature & !  kappa_1/dz_0.5 * T_{-1})
         )

    denom = gamma0(this, icell, iconn)

    delta0 = numer/denom

  end function delta0

  !------------------------------------------------------------------------
  function gamma0(this, icell, iconn)
    !
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    class(goveqn_cair_vapor_type)          :: this
    PetscInt                               :: icell, iconn
    PetscReal                              :: gamma0
    !
    class(cair_vapor_auxvar_type), pointer :: auxvar(:)
    PetscReal                              :: qsat, s0, gs_0, lambda

    auxvar => this%aux_vars_in
    lambda = HVAP * MM_H2O

    call SatVap(auxvar(icell)%temperature, qsat, s0)
    qsat = qsat/this%aux_vars_in(icell)%pref
    s0   = s0  /this%aux_vars_in(icell)%pref

    gs_0 = gs0(this, icell, iconn)

    ! The ground temperatue contributes to the canopy air layer above the ground
    gamma0 = &
         this%aux_vars_in(icell)%cpair *this%aux_vars_conn_in(iconn)%ga + & ! Cp * ga
         lambda * this%aux_vars_in(icell)%soil_rhg * gs_0 * s0           + & ! lambda * h0 * gs0 * s0
         this%aux_vars_in(icell)%soil_tk/this%aux_vars_in(icell)%soil_dz    ! k/dz

  end function gamma0

  !------------------------------------------------------------------------
  subroutine CAirVaporComputeRhsAccumulation(this, b_p)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    !
    use WaterVaporMod             , only : SatVap
    use MultiPhysicsProbConstants , only : MM_H2O, MM_DRY_AIR
    use MultiPhysicsProbConstants , only : HVAP
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    Vec                           :: B
    PetscErrorCode                :: ierr
    !
    PetscInt                                :: icell, ileaf, iconn
    PetscScalar                   , pointer :: b_p(:)
    class(cair_vapor_auxvar_type) , pointer :: auxvar(:)
    class(connection_set_type)    , pointer :: cur_conn_set
    PetscReal                               :: qsat, s0, si, gleaf, gleaf_et, factor, gs0_value
    PetscReal                               :: delta0_value, value, lambda

    auxvar => this%aux_vars_in
    cur_conn_set => this%mesh%intrn_conn_set_list%first

    lambda = HVAP * MM_H2O

    ! Internal cells
    do icell = 1, this%mesh%ncells_local
       if (auxvar(icell)%is_soil) then

          ! Find the internal connection that is between the soil layer overlying
          ! air space
          iconn = SoilAirIConn(this, icell)

          call SatVap(auxvar(icell)%temperature, qsat, s0)
          qsat = qsat/this%aux_vars_in(icell)%pref
          s0   = s0  /this%aux_vars_in(icell)%pref

          ! The specific humidity (or partial pressure of water vapor) does not change
          ! for gound cell
          b_p(icell) = this%aux_vars_in(icell)%water_vapor/this%aux_vars_in(icell)%pref

          ! The specific humidity (or partial pressure of water vapor) at ground
          ! contributes to the canopy air cell above the ground
          delta0_value = delta0(this, icell, iconn)

          gs0_value = gs0(this, icell, iconn)

          value = gs0_value * &
               auxvar(icell)%soil_rhg*(qsat + s0*(delta0_value - auxvar(icell)%temperature))

#ifdef USE_BONAN_FORMULATION
          b_p(icell+1) = b_p(icell+1) + value
#else
          b_p(icell+1) = b_p(icell+1) + value/this%mesh%vol(icell+1)
#endif
       else

          factor = 1.d0/ (MM_H2O / MM_DRY_AIR + (1.d0 - MM_H2O / MM_DRY_AIR))
#ifdef USE_BONAN_FORMULATION
          b_p(icell) = b_p(icell) + &
               auxvar(icell)%rhomol / this%dtime * factor * auxvar(icell)%water_vapor/this%aux_vars_in(icell)%pref * this%mesh%vol(icell)
#else
          b_p(icell) = b_p(icell) + &
               auxvar(icell)%rhomol / this%dtime * factor * auxvar(icell)%water_vapor/this%aux_vars_in(icell)%pref
#endif
          do ileaf = 1,auxvar(icell)%nleaf
             call SatVap(auxvar(icell)%leaf_temperature(ileaf), qsat, si)
             qsat = qsat/this%aux_vars_in(icell)%pref
             si   = si  /this%aux_vars_in(icell)%pref

             ! Add contribution from leaf, if present in the given layer
             if (auxvar(icell)%leaf_dpai(ileaf) > 0.d0) then
                gleaf = &
                     auxvar(icell)%leaf_gs(ileaf) * auxvar(icell)%gbv/ &
                     (auxvar(icell)%leaf_gs(ileaf) + auxvar(icell)%gbv)

                gleaf_et = &
                     gleaf             * auxvar(icell)%leaf_fdry(ileaf) + &
                     auxvar(icell)%gbv * auxvar(icell)%leaf_fwet(ileaf)

                gleaf_et = gleaf_et * auxvar(icell)%leaf_fssh(ileaf) * auxvar(icell)%leaf_dpai(ileaf)

#ifdef USE_BONAN_FORMULATION
                b_p(icell) = b_p(icell) &
                     + gleaf_et * (qsat - si*auxvar(icell)%leaf_temperature(ileaf))
#else
                b_p(icell) = b_p(icell) &
                     + gleaf_et * (qsat - si*auxvar(icell)%leaf_temperature(ileaf)) /this%mesh%vol(icell)
#endif
             end if
          end do
       end if
    end do

  end subroutine CAirVaporComputeRhsAccumulation

  !------------------------------------------------------------------------
  subroutine CAirVaporRhsDivergence(this, b_p)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    use MultiPhysicsProbConstants , only : HVAP
    use WaterVaporMod             , only : SatVap
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : MM_H2O, MM_DRY_AIR
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type)           :: this
    PetscScalar                   , pointer :: b_p(:)
    !
    type(condition_type)          , pointer :: cur_cond
    class(connection_set_type)    , pointer :: cur_conn_set
    class(cair_vapor_auxvar_type) , pointer :: auxvar(:)
    PetscReal                               :: dist_up, dist_dn, dist, factor
    PetscInt                                :: iconn, sum_conn, cell_id

    ! Boundary cells
    auxvar => this%aux_vars_bc
    cur_cond => this%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          sum_conn = sum_conn + 1

          select case(cur_cond%itype)
          case(COND_DIRICHLET)

             dist_up = cur_conn_set%conn(iconn)%GetDistUp()
             dist_dn = cur_conn_set%conn(iconn)%GetDistDn()
             dist    = dist_up + dist_dn

             cell_id  = cur_conn_set%conn(iconn)%GetIDDn()

             factor = 1.d0/ (MM_H2O / MM_DRY_AIR + (1.d0 - MM_H2O / MM_DRY_AIR) * auxvar(sum_conn)%water_vapor)
#ifdef USE_BONAN_FORMULATION
             b_p(cell_id) = b_p(cell_id) + &
                  this%aux_vars_conn_bc(sum_conn)%ga * auxvar(sum_conn)%water_vapor * factor
#else
             b_p(cell_id) = b_p(cell_id) + &
                  this%aux_vars_conn_bc(sum_conn)%ga/dist * auxvar(sum_conn)%water_vapor * factor
#endif

          case default
             write(iulog,*)'Unknown boundary condition type'
             call endrun(msg=errMsg(__FILE__, __LINE__))

          end select
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine CAirVaporRhsDivergence

  !------------------------------------------------------------------------

  subroutine CAirVaporComputeOperatorsDiag(this, A, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use mpp_varcon                , only : cnfac
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    Mat                           :: A
    Mat                           :: B
    PetscErrorCode                :: ierr
    !
    ! !LOCAL VARIABLES
    class(connection_set_type) , pointer :: cur_conn_set
    type(condition_type)       , pointer :: cur_cond
    PetscInt                             :: icell, ileaf, iconn, sum_conn
    PetscInt                             :: cell_id, cell_id_up, cell_id_dn
    PetscInt                             :: row, col
    PetscReal                            :: value, ga, dist, gsw, gsa, gs0_value
    PetscReal                            :: qsat, si, gleaf, gleaf_et
    PetscReal                            :: beta0_value, lambda

    lambda = HVAP * MM_H2O

    cur_conn_set => this%mesh%intrn_conn_set_list%first

    ! For interior and top cell
    do icell = 1, this%mesh%ncells_local
       row = icell-1; col = icell-1

       if (this%aux_vars_in(icell)%is_soil) then

          ! The specific humidity (or partial pressure of water vapor) does not change
          ! for gound cell
          value = 1.d0
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

          ! The specific humidity (or partial pressure of water vapor) at ground
          ! contributes to the canopy air cell above the ground

          ! Find the internal connection that is between the soil layer overlying
          ! air space
          iconn = SoilAirIConn(this, icell)

          call SatVap(this%aux_vars_in(icell)%temperature, qsat, si)
          qsat = qsat/this%aux_vars_in(icell)%pref
          si   = si  /this%aux_vars_in(icell)%pref

          gs0_value = gs0(this, icell, iconn)

          beta0_value = beta0(this, icell, iconn)

#ifdef USE_BONAN_FORMULATION
          value = gs0_value * si * this%aux_vars_in(icell)%soil_rhg * beta0_value
#else
          value = gs0_value * si * this%aux_vars_in(icell)%soil_rhg * beta0_value/this%mesh%vol(icell+1)
#endif

          row = icell; col = icell;
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

       else

#ifdef USE_BONAN_FORMULATION
          value = this%aux_vars_in(icell)%rhomol/this%dtime*this%mesh%vol(icell)
#else
          value = this%aux_vars_in(icell)%rhomol/this%dtime
#endif
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

          do ileaf = 1, this%aux_vars_in(icell)%nleaf
             if (this%aux_vars_in(icell)%leaf_dpai(ileaf) > 0.d0) then

                call SatVap(this%aux_vars_in(icell)%leaf_temperature(ileaf), qsat, si)
                qsat = qsat/this%aux_vars_in(icell)%pref
                si   = si  /this%aux_vars_in(icell)%pref

                gleaf = &
                     this%aux_vars_in(icell)%leaf_gs(ileaf) * this%aux_vars_in(icell)%gbv/ &
                     (this%aux_vars_in(icell)%leaf_gs(ileaf) + this%aux_vars_in(icell)%gbv)

                gleaf_et = &
                     gleaf                       * this%aux_vars_in(icell)%leaf_fdry(ileaf) + &
                     this%aux_vars_in(icell)%gbv * this%aux_vars_in(icell)%leaf_fwet(ileaf)

                gleaf_et = gleaf_et * this%aux_vars_in(icell)%leaf_fssh(ileaf) * this%aux_vars_in(icell)%leaf_dpai(ileaf)

#ifdef USE_BONAN_FORMULATION
                value = gleaf_et
#else
                value = gleaf_et/this%mesh%vol(icell)
#endif

                call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
             end if
          end do
       end if

    enddo

    ! Interior cells
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections

          cell_id_up = cur_conn_set%conn(iconn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(iconn)%GetIDDn()

          dist = cur_conn_set%conn(iconn)%GetDistUp() + &
               cur_conn_set%conn(iconn)%GetDistDn()

          if (this%aux_vars_in(cell_id_up)%is_soil .or. &
               this%aux_vars_in(cell_id_dn)%is_soil ) then

             if (this%aux_vars_in(cell_id_up)%is_soil) then
                cell_id = cell_id_up
             else
                cell_id = cell_id_dn
             end if

             gs0_value = gs0(this, cell_id, iconn)

#ifdef USE_BONAN_FORMULATION
             value = gs0_value
#else
             value = gs0_value/dist
#endif
          else
             ga = this%aux_vars_conn_in(iconn)%ga
#ifdef USE_BONAN_FORMULATION
             value = ga
#else
             value = ga/dist
#endif
          end if

          if (.not.this%aux_vars_in(cell_id_up)%is_soil) then

             ! Skip if cell_id_up = canopy-air and cell_id_dn = soil
             if (.not.this%aux_vars_in(cell_id_dn)%is_soil) then
                call MatSetValuesLocal(B, 1, cell_id_up-1, 1, cell_id_dn-1, -value, ADD_VALUES, ierr); CHKERRQ(ierr)
             end if

             call MatSetValuesLocal(B, 1, cell_id_up-1, 1, cell_id_up-1,  value, ADD_VALUES, ierr); CHKERRQ(ierr)
          endif

          if (.not.this%aux_vars_in(cell_id_dn)%is_soil) then

             ! Skip if cell_id_dn = canopy-air and cell_id_up = soil
             if (.not.this%aux_vars_in(cell_id_up)%is_soil) then
                call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_up-1, -value, ADD_VALUES, ierr); CHKERRQ(ierr)
             end if

             call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_dn-1,  value, ADD_VALUES, ierr); CHKERRQ(ierr)
          end if

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    ! Boundary cells
    cur_cond => this%boundary_conditions%first
    sum_conn = 0
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections

          sum_conn = sum_conn + 1

          select case(cur_cond%itype)
          case(COND_DIRICHLET)

             dist = cur_conn_set%conn(iconn)%GetDistUp() + &
                  cur_conn_set%conn(iconn)%GetDistDn()

             cell_id  = cur_conn_set%conn(iconn)%GetIDDn()
             row = cell_id-1; col = cell_id-1;

#ifdef USE_BONAN_FORMULATION
             value = this%aux_vars_conn_bc(sum_conn)%ga
#else
             value = this%aux_vars_conn_bc(sum_conn)%ga/dist
#endif

             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

          case default
             write(iulog,*)'Unknown boundary condition type'
             call endrun(msg=errMsg(__FILE__, __LINE__))

          end select
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine CAirVaporComputeOperatorsDiag

  !------------------------------------------------------------------------

  subroutine CAirVaporComputeOperatorsOffDiag(this, A, B, &
       itype_of_other_goveq, rank_of_other_goveq, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use mpp_varcon                , only : cnfac
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : HVAP
    use WaterVaporMod             , only : SatVap
    use MultiPhysicsProbConstants , only : GE_CANOPY_AIR_TEMP
    use MultiPhysicsProbConstants , only : GE_CANOPY_LEAF_TEMP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_vapor_type) :: this
    Mat                           :: A
    Mat                           :: B
    PetscInt                      :: itype_of_other_goveq
    PetscInt                      :: rank_of_other_goveq
    PetscErrorCode                :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                      :: icell, ileaf, row, col, cair_auxvar_idx, leaf_idx, geq_leaf_temp_rank, iconn
    PetscReal                     :: value, qsat, si, gs0, gsa, gsw, alpha0_value
    PetscReal                     :: gleaf, gleaf_et

    select case (itype_of_other_goveq)
       case (GE_CANOPY_AIR_TEMP)
          ! For soil-vapror <--> air temperature
          do icell = 1, this%mesh%ncells_local
             if (this%aux_vars_in(icell)%is_soil) then
             
                ! The specific humidity (or partial pressure of water vapor) does not change
                ! for gound cell
                row = icell-1; col = icell-1
                value = 0.d0
                call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

                ! The specific humidity (or partial pressure of water vapor) at ground
                ! contributes to the canopy air cell above the ground
                call SatVap(this%aux_vars_in(icell)%temperature, qsat, si)
                qsat = qsat/this%aux_vars_in(icell)%pref
                si   = si  /this%aux_vars_in(icell)%pref

                iconn = SoilAirIConn(this, icell)
                gsw = 1.d0 / this%aux_vars_in(icell)%soil_resis * this%aux_vars_in(icell)%rhomol
                gsa = this%aux_vars_conn_in(iconn)%ga
                gs0 = gsw * gsa / (gsw + gsa)

                alpha0_value = alpha0(this, icell, iconn)

#ifdef USE_BONAN_FORMULATION
                value = - this%aux_vars_in(icell)%soil_rhg * si * gs0 * alpha0_value
#else
                value = - this%aux_vars_in(icell)%soil_rhg * si * gs0 * alpha0_value / this%mesh%vol(icell+1)
#endif
                row = icell; col = icell;
                call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
             end if
          end do

    case (GE_CANOPY_LEAF_TEMP)
       geq_leaf_temp_rank = rank_of_other_goveq
       do ileaf = 1, this%nLeaf
          cair_auxvar_idx = this%LeafGE2CAirAuxVar_in(ileaf,1)
          leaf_idx        = this%LeafGE2CAirAuxVar_in(ileaf,2) + &
                            (geq_leaf_temp_rank-1)*this%aux_vars_in(cair_auxvar_idx)%nleaf/this%nleafGE

          if (this%aux_vars_in(cair_auxvar_idx)%leaf_dpai(leaf_idx) > 0.d0) then
             row = cair_auxvar_idx-1; col = ileaf-1

             call SatVap(this%aux_vars_in(cair_auxvar_idx)%leaf_temperature(leaf_idx), qsat, si)
             qsat = qsat/this%aux_vars_in(cair_auxvar_idx)%pref
             si   = si  /this%aux_vars_in(cair_auxvar_idx)%pref

             gleaf = &
                  this%aux_vars_in(cair_auxvar_idx)%leaf_gs(leaf_idx) * this%aux_vars_in(cair_auxvar_idx)%gbv/ &
                  (this%aux_vars_in(cair_auxvar_idx)%leaf_gs(leaf_idx) + this%aux_vars_in(cair_auxvar_idx)%gbv)

             gleaf_et = &
                  gleaf                       * this%aux_vars_in(cair_auxvar_idx)%leaf_fdry(leaf_idx) + &
                  this%aux_vars_in(cair_auxvar_idx)%gbv * this%aux_vars_in(cair_auxvar_idx)%leaf_fwet(leaf_idx)

             gleaf_et = gleaf_et * this%aux_vars_in(cair_auxvar_idx)%leaf_fssh(leaf_idx) * this%aux_vars_in(cair_auxvar_idx)%leaf_dpai(leaf_idx)

#ifdef USE_BONAN_FORMULATION
             value = -si*gleaf_et
#else
             value = -si*gleaf_et/this%mesh%vol(cair_auxvar_idx)
#endif

             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
          end if

       end do
    end select
       
    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  end subroutine CAirVaporComputeOperatorsOffDiag

#endif

end module GoveqnCanopyAirVaporType
