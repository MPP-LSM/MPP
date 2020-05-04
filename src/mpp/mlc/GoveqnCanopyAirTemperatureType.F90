module GoveqnCanopyAirTemperatureType

#ifdef USE_PETSC_LIB
!#define USE_BONAN_FORMULATION

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                      , only : iulog
  use mpp_abortutils                  , only : endrun
  use mpp_shr_log_mod                 , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType       , only : goveqn_base_type
  use CanopyAirTemperatureAuxType     , only : cair_temp_auxvar_type
  use CanopyAirTemperatureConnAuxType , only : cair_temp_conn_auxvar_type
  use CanopyTurbulenceAuxType         , only : canopy_turbulence_auxvar_type
  use petscvec
  use petscmat
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_cair_temp_type

     type(cair_temp_auxvar_type)      , pointer :: aux_vars_in(:)
     type(cair_temp_auxvar_type)      , pointer :: aux_vars_bc(:)
     type(cair_temp_conn_auxvar_type) , pointer :: aux_vars_conn_in(:)
     type(cair_temp_conn_auxvar_type) , pointer :: aux_vars_conn_bc(:)

   contains

     procedure, public :: Setup                     => CAirTempSetup
     procedure, public :: AllocateAuxVars           => CAirTempAllocateAuxVars
     procedure, public :: PreSolve                  => CAirTempPreSolve
     procedure, public :: GetFromSoeAuxVarsCturb    => CAirTempGetFromSoeAuxVarsCtrub
     procedure, public :: SavePrimaryIndependentVar => CAirTempSavePrmIndepVar
     procedure, public :: GetRValues                => CAirTempGetRValues
     procedure, public :: SetRValues                => CAirTempSetRValues
     procedure, public :: ComputeRhs                => CAirTempComputeRhs
     procedure, public :: ComputeOperatorsDiag      => CAirTempComputeOperatorsDiag
     procedure, public :: ComputeOperatorsOffDiag   => CAirTempComputeOperatorsOffDiag

  end type goveqn_cair_temp_type

contains

  !------------------------------------------------------------------------
  subroutine CAirTempSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_TEMP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_temp_type) :: this

    call this%Create()

    this%name  = ""
    this%itype = GE_CANOPY_AIR_TEMP

    nullify(this%aux_vars_in)
    nullify(this%aux_vars_bc)
    nullify(this%aux_vars_conn_in)
    nullify(this%aux_vars_conn_bc)

  end subroutine CAirTempSetup

  !------------------------------------------------------------------------
  subroutine CAirTempAllocateAuxVars(this)
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
    class(goveqn_cair_temp_type) :: this
    !
    type(condition_type)         , pointer :: cur_cond
    type(coupling_variable_type) , pointer :: cpl_var
    class(connection_set_type)   , pointer :: cur_conn_set
    PetscInt                               :: ncells_cond
    PetscInt                               :: ncond
    PetscInt                               :: icond
    PetscInt                               :: nleaf
    PetscInt                               :: sum_conn

    !
    ! Determine the number of leaves that are coupled with the canopy
    ! air space
    nleaf = 0
    cpl_var => this%coupling_vars%first
    do
       if (.not.associated(cpl_var)) exit
       if ( cpl_var%variable_type == VAR_LEAF_TEMPERATURE) then
          nleaf = nleaf + 1
       end if
       cpl_var => cpl_var%next
    end do

    ! Allocate memory and initialize aux vars: For internal connections
    allocate(this%aux_vars_in(this%mesh%ncells_all))
    do icond = 1,this%mesh%ncells_all
       call this%aux_vars_in(icond)%Init(nleaf)
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
       call this%aux_vars_bc(icond)%Init(nleaf)
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

  end subroutine CAirTempAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine CAirTempPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_TEMP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_temp_type) :: this

    write(*,*)' In CAirTempPreSolve'
    call exit(0)
  end subroutine CAirTempPreSolve

  !------------------------------------------------------------------------
  subroutine CAirTempGetFromSoeAuxVarsCtrub(this, cturb)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use CanopyTurbulenceAuxType, only : canopy_turbulence_auxvar_type 
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_temp_type) :: this
    class(canopy_turbulence_auxvar_type) :: cturb
    !
    PetscInt :: ncair, icair, icell, level, iconn

    ncair = 1

    ! all layers
    do icair = 1, ncair
       do icell = 1, this%mesh%ncells_local
          level = icell
          this%aux_vars_in(icell)%cpair  = cturb%cpair(icair)
          this%aux_vars_in(icell)%rhomol = cturb%rhomol(icair)
          this%aux_vars_in(icell)%pref   = cturb%pref(icair)
       end do
    end do

    icair = 1
    do icell = 1, this%mesh%ncells_local
       level = icell
       this%aux_vars_conn_in(icell)%ga = cturb%ga_prof(icair,level)
    enddo

    ! soil-layer
    icair = 1
    icell = 1
    this%aux_vars_in(icell)%temperature = cturb%tsoi(icair)
    this%aux_vars_in(icell)%soil_rhg    = cturb%rhgsoi(icair)
    this%aux_vars_in(icell)%soil_rn     = cturb%rnsoi(icair)
    this%aux_vars_in(icell)%soil_tk     = cturb%tksoi(icair)
    this%aux_vars_in(icell)%soil_dz     = cturb%dzsoi(icair)
    this%aux_vars_in(icell)%soil_resis  = cturb%ressoi(icair)
    this%aux_vars_in(icell)%soil_temperature = cturb%tsoi(icair)

    ! top-layer
    icair = 1
    iconn = 1
    do icair = 1, ncair
       level = this%mesh%ncells_local
       this%aux_vars_conn_bc(iconn)%ga     = cturb%ga_prof(icair, level)
       this%aux_vars_bc(iconn)%temperature = cturb%thref(icair)
    end do

  end subroutine CAirTempGetFromSoeAuxVarsCtrub

  !------------------------------------------------------------------------
  subroutine CAirTempSavePrmIndepVar (this, x)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use CanopyTurbulenceAuxType, only : canopy_turbulence_auxvar_type 
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_temp_type) :: this
    Vec :: x
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
       this%aux_vars_in(ghosted_id)%temperature = x_p(ghosted_id)
    end do

    call VecRestoreArrayF90(x, x_p, ierr); CHKERRQ(ierr)

  end subroutine CAirTempSavePrmIndepVar

  !------------------------------------------------------------------------
  subroutine CAirTempGetRValues (this, auxvar_type, var_type, nauxvar, var_values)
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
    class(goveqn_cair_temp_type) :: this
    PetscInt                     :: auxvar_type
    PetscInt                     :: var_type
    PetscInt                     :: nauxvar
    PetscReal, pointer           :: var_values(:)

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call CAirTempGetRValuesFromAuxVars(this%aux_vars_in, var_type, nauxvar, var_values)
    case (AUXVAR_BC)
       call CAirTempGetRValuesFromAuxVars(this%aux_vars_bc, var_type, nauxvar, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CAirTempGetRValues

  !------------------------------------------------------------------------
  subroutine CAirTempSetRValues (this, auxvar_type, var_type, leaf_idx, nauxvar, var_values)
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
    class(goveqn_cair_temp_type) :: this
    PetscInt                     :: auxvar_type
    PetscInt                     :: var_type
    PetscInt                     :: leaf_idx
    PetscInt                     :: nauxvar
    PetscReal, pointer           :: var_values(:)

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call CAirTempSetRValuesFromAuxVars(this%aux_vars_in, var_type, leaf_idx, nauxvar, var_values)
    case (AUXVAR_BC)
       call CAirTempSetRValuesFromAuxVars(this%aux_vars_bc, var_type, leaf_idx, nauxvar, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CAirTempSetRValues

  !------------------------------------------------------------------------
  subroutine CAirTempGetRValuesFromAuxVars (aux_var, var_type, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(cair_temp_auxvar_type) , pointer :: aux_var(:)
    PetscInt                              :: var_type
    PetscInt                              :: nauxvar
    PetscReal                   , pointer :: var_values(:)
    !
    PetscInt :: iauxvar

    select case(var_type)
    case(VAR_TEMPERATURE)
       do iauxvar = 1,nauxvar
          var_values(iauxvar) = aux_var(iauxvar)%temperature
       end do
    case default
       write(iulog,*) 'CAirTempGetRValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CAirTempGetRValuesFromAuxVars

  !------------------------------------------------------------------------
  subroutine CAirTempSetRValuesFromAuxVars (aux_var, var_type, leaf_idx, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_LEAF_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_WATER_VAPOR
    !
    implicit none
    !
    ! !ARGUMENTS
    type(cair_temp_auxvar_type) , pointer :: aux_var(:)
    PetscInt                              :: var_type
    PetscInt                              :: leaf_idx
    PetscInt                              :: nauxvar
    PetscReal                   , pointer :: var_values(:)
    !
    PetscInt :: iauxvar

    select case(var_type)
    case(VAR_LEAF_TEMPERATURE)
       do iauxvar = 1,nauxvar
          aux_var(iauxvar)%leaf_temperature(leaf_idx) = var_values(iauxvar)
       end do
    case(VAR_WATER_VAPOR)
       do iauxvar = 1,nauxvar
          aux_var(iauxvar)%vcan = var_values(iauxvar)
       end do
    case default
       write(iulog,*) 'CAirTempGetSValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
       
  end subroutine CAirTempSetRValuesFromAuxVars

  !------------------------------------------------------------------------
  subroutine CAirTempComputeRhs(this, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    use MultiPhysicsProbConstants, only : HVAP
    use WaterVaporMod, only : SatVap
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_temp_type) :: this
    Vec                          :: B
    PetscErrorCode               :: ierr
    !
    PetscScalar                  , pointer :: b_p(:)

    call VecGetArrayF90(B, b_p, ierr)

    call CAirTempRhsAccumulation(this, b_p)
    call CAirTempRhsDivergence(this, b_p)

    call VecRestoreArrayF90(B, b_p, ierr)

  end subroutine CAirTempComputeRhs

  !------------------------------------------------------------------------
  subroutine CAirTempRhsAccumulation(this, b_p)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    use WaterVaporMod             , only : SatVap
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_temp_type) :: this
    PetscScalar                  , pointer :: b_p(:)
    !
    PetscInt                               :: icell
    class(cair_temp_auxvar_type) , pointer :: auxvar(:)
    PetscReal                              :: qsat, si, gsw, gsa, gs0, lambda

    auxvar => this%aux_vars_in

    lambda = HVAP * MM_H2O

    ! Soil layer
    icell = 1
    call SatVap(auxvar(icell)%temperature, qsat, si)
    qsat = qsat/this%aux_vars_in(icell)%pref
    si   = si  /this%aux_vars_in(icell)%pref

    gsw = 1.d0 / auxvar(icell)%soil_resis * auxvar(icell)%rhomol
    gsa = this%aux_vars_conn_in(1)%ga
    gs0 = gsw * gsa / (gsw + gsa)
    
    b_p(icell) = b_p(icell) &
         + ( auxvar(icell)%soil_rn &
            - lambda * auxvar(icell)%soil_rhg * gs0* (qsat - si*auxvar(icell)%temperature) &
            + auxvar(icell)%soil_tk/auxvar(icell)%soil_dz * auxvar(icell)%soil_temperature &
            )
    
    ! Internal layers
    do icell = 2, this%mesh%ncells_local-1
#ifdef USE_BONAN_FORMULATION
       b_p(icell) = b_p(icell) + auxvar(icell)%rhomol / this%dtime * auxvar(icell)%temperature * this%mesh%vol(icell)
#else
       b_p(icell) = b_p(icell) + auxvar(icell)%rhomol / this%dtime * auxvar(icell)%temperature
#endif
    end do

    ! Top layer
    icell = this%mesh%ncells_local
#ifdef USE_BONAN_FORMULATION
    b_p(icell) = auxvar(icell)%rhomol / this%dtime * auxvar(icell)%temperature * this%mesh%vol(icell)
#else
    b_p(icell) = auxvar(icell)%rhomol / this%dtime * auxvar(icell)%temperature
#endif

  end subroutine CAirTempRhsAccumulation

  !------------------------------------------------------------------------
  subroutine CAirTempRhsDivergence(this, b_p)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_temp_type)           :: this
    PetscScalar                  , pointer :: b_p(:)
    !
    class(cair_temp_auxvar_type) , pointer :: auxvar(:)
    type(condition_type)         , pointer :: cur_cond
    class(connection_set_type)   , pointer :: cur_conn_set
    PetscReal                              :: dist_up, dist_dn, dist
    PetscInt                               :: iconn, sum_conn, cell_id

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

#ifdef USE_BONAN_FORMULATION
             b_p(cell_id) = b_p(cell_id) + &
                  this%aux_vars_conn_bc(sum_conn)%ga * auxvar(sum_conn)%temperature
#else
             b_p(cell_id) = b_p(cell_id) + &
                  this%aux_vars_conn_bc(sum_conn)%ga/dist * auxvar(sum_conn)%temperature
#endif

          case default
             write(iulog,*)'Unknown boundary condition type'
             call endrun(msg=errMsg(__FILE__, __LINE__))

          end select
       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine CAirTempRhsDivergence

  !------------------------------------------------------------------------

  subroutine CAirTempComputeOperatorsDiag(this, A, B, ierr)
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
    class(goveqn_cair_temp_type) :: this
    Mat                          :: A
    Mat                          :: B
    PetscErrorCode               :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                             :: icell, ileaf, iconn, sum_conn
    PetscInt                             :: cell_id, cell_id_up, cell_id_dn
    PetscInt                             :: row, col
    PetscReal                            :: value, ga, dist, gsw, gsa, gs0
    PetscReal                            :: qsat, si, lambda
    class(connection_set_type) , pointer :: cur_conn_set
    type(condition_type)       , pointer :: cur_cond

    lambda = HVAP * MM_H2O

    ! For soil cell
    icell = 1
    row = icell-1; col = icell-1;
    call SatVap(this%aux_vars_in(icell)%temperature, qsat, si)
    qsat = qsat/this%aux_vars_in(icell)%pref
    si   = si  /this%aux_vars_in(icell)%pref

    gsw = 1.d0 / this%aux_vars_in(icell)%soil_resis * this%aux_vars_in(icell)%rhomol
    gsa = this%aux_vars_conn_in(1)%ga
    gs0 = gsw * gsa / (gsw + gsa)

    value = &
         this%aux_vars_in(icell)%cpair *this%aux_vars_conn_in(1)%ga + &
         lambda * this%aux_vars_in(icell)%soil_rhg * gs0 * si + &
         this%aux_vars_in(icell)%soil_tk/this%aux_vars_in(icell)%soil_dz

    call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

    ! For soil-cell <--> air temperature
    row = icell-1; col = icell+1-1
    value = -this%aux_vars_in(icell)%cpair * this%aux_vars_conn_in(1)%ga

    call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
    
    ! For interior and top cell
    do icell = 2, this%mesh%ncells_local
       row = icell-1; col = icell-1

#ifdef USE_BONAN_FORMULATION
       value = this%aux_vars_in(icell)%rhomol/this%dtime * this%mesh%vol(icell)
#else
       value = this%aux_vars_in(icell)%rhomol/this%dtime
#endif

       call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

       do ileaf = 1, this%aux_vars_in(icell)%nleaf
          if (this%aux_vars_in(icell)%leaf_dpai(ileaf) > 0.d0) then

#ifdef USE_BONAN_FORMULATION
             value = 2.d0 * this%aux_vars_in(icell)%gbh * &
                  this%aux_vars_in(icell)%leaf_fssh(ileaf) * &
                  this%aux_vars_in(icell)%leaf_dpai(ileaf)
#else
             value = 2.d0 * this%aux_vars_in(icell)%gbh * &
                  this%aux_vars_in(icell)%leaf_fssh(ileaf) * &
                  this%aux_vars_in(icell)%leaf_dpai(ileaf) / &
                  this%mesh%vol(icell)
#endif

             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
          end if
       end do

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

          ga = this%aux_vars_conn_in(iconn)%ga
#ifdef USE_BONAN_FORMULATION
          value = ga
#else
          value = ga/dist
#endif

          if (cell_id_up > 1) then
             call MatSetValuesLocal(B, 1, cell_id_up-1, 1, cell_id_dn-1, -value, ADD_VALUES, ierr); CHKERRQ(ierr)
             call MatSetValuesLocal(B, 1, cell_id_up-1, 1, cell_id_up-1,  value, ADD_VALUES, ierr); CHKERRQ(ierr)
          endif

          if (cell_id_dn > 1) then
             call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_up-1, -value, ADD_VALUES, ierr); CHKERRQ(ierr)
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

  end subroutine CAirTempComputeOperatorsDiag

  !------------------------------------------------------------------------

  subroutine CAirTempComputeOperatorsOffDiag(this, A, B, &
       itype_of_other_goveq, rank_of_other_goveq, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    use MultiPhysicsProbConstants , only : GE_CANOPY_AIR_VAPOR
    use MultiPhysicsProbConstants , only : GE_CANOPY_LEAF_TEMP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cair_temp_type) :: this
    Mat                          :: A
    Mat                          :: B
    PetscInt                     :: itype_of_other_goveq
    PetscInt                     :: rank_of_other_goveq
    PetscErrorCode               :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt :: icell, ileaf, row, col
    PetscReal :: value, gleaf, gleaf_et, gsw, gsa, gs0, lambda

    lambda = HVAP * MM_H2O

    select case (itype_of_other_goveq)
    case (GE_CANOPY_AIR_VAPOR)
       ! For soil-cell <--> air vapor
       icell = 1
       row = icell-1; col = icell-1

       gsw = 1.d0 / this%aux_vars_in(icell)%soil_resis * this%aux_vars_in(icell)%rhomol
       gsa = this%aux_vars_conn_in(1)%ga
       gs0 = gsw * gsa / (gsw + gsa)

       value = - lambda * gs0
       call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

    case (GE_CANOPY_LEAF_TEMP)
       ileaf = rank_of_other_goveq - 2
       do icell = 1, this%mesh%ncells_local
          row = icell-1; col = icell-1
          if (this%aux_vars_in(icell)%leaf_dpai(ileaf) > 0.d0) then

#ifdef USE_BONAN_FORMULATION
             value = -2.d0 * this%aux_vars_in(icell)%gbh * &
                  this%aux_vars_in(icell)%leaf_fssh(ileaf) * &
                  this%aux_vars_in(icell)%leaf_dpai(ileaf)
#else
             value = -2.d0 * this%aux_vars_in(icell)%gbh * &
                  this%aux_vars_in(icell)%leaf_fssh(ileaf) * &
                  this%aux_vars_in(icell)%leaf_dpai(ileaf) / &
                  this%mesh%vol(icell)
#endif

             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
          end if
       end do

    end select

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  end subroutine CAirTempComputeOperatorsOffDiag

#endif

end module GoveqnCanopyAirTemperatureType
