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

  end subroutine CAirVaporSetup

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
    PetscInt                               :: icond
    PetscInt                               :: nleaf
    PetscInt                               :: sum_conn

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
    PetscInt :: ncair, icair, icell, level, iconn

    ncair = 1

    ! all layers
    do icair = 1, ncair
       do icell = 1, this%mesh%ncells_local
          level = icell
          this%aux_vars_in(icell)%cpair    = cturb%cpair(icair)
          this%aux_vars_in(icell)%rhomol   = cturb%rhomol(icair)
          this%aux_vars_in(icell)%pref     = cturb%pref(icair)

          if (this%aux_vars_in(icell)%is_soil) then
             this%aux_vars_in(icell)%soil_resis  = cturb%ressoi(icair)
             this%aux_vars_in(icell)%soil_rhg    = cturb%rhgsoi(icair)
             this%aux_vars_in(icell)%temperature = cturb%tsoi(icair)
          endif
       end do
    end do

    ! internal connections
    do icair = 1, ncair
       do icell = 1, this%mesh%ncells_local-1
          level = icell
          this%aux_vars_conn_in(icell)%ga    = cturb%ga_prof(icair,level)
       end do

    end do

    ! top-layer
    iconn = 1
    do icair = 1, ncair
       level = this%mesh%ncells_local
       this%aux_vars_conn_bc(iconn)%ga = cturb%ga_prof(icair, level)
       this%aux_vars_bc(iconn)%water_vapor    = cturb%vref(icair)
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
  subroutine CAirVaporSetRValues (this, auxvar_type, var_type, leaf_idx, nauxvar, var_values)
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
    PetscInt                      :: leaf_idx
    PetscInt                      :: nauxvar
    PetscReal, pointer            :: var_values(:)

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call CAirVaporSetRValuesFromAuxVars(this%aux_vars_in, var_type, leaf_idx, nauxvar, var_values)
    case (AUXVAR_BC)
       call CAirVaporSetRValuesFromAuxVars(this%aux_vars_bc, var_type, leaf_idx, nauxvar, var_values)
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
  subroutine CAirVaporSetRValuesFromAuxVars (aux_var, var_type, leaf_idx, nauxvar, var_values)
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
    PetscInt                               :: leaf_idx
    PetscInt                               :: nauxvar
    PetscReal                    , pointer :: var_values(:)
    !
    PetscInt :: iauxvar

    select case(var_type)
    case(VAR_LEAF_TEMPERATURE)
       do iauxvar = 1,nauxvar
          aux_var(iauxvar)%leaf_temperature(leaf_idx) = var_values(iauxvar)
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
  subroutine CAirVaporComputeRhsAccumulation(this, b_p)
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
    PetscInt                                :: icell, ileaf
    PetscScalar                   , pointer :: b_p(:)
    class(cair_vapor_auxvar_type) , pointer :: auxvar(:)
    PetscReal                               :: qsat, si, gleaf, gleaf_et

    auxvar => this%aux_vars_in

    icell = 1
    
    ! Internal cells
    do icell = 1, this%mesh%ncells_local
       if (auxvar(icell)%is_soil) then

          call SatVap(auxvar(icell)%temperature, qsat, si)
          qsat = qsat/this%aux_vars_in(icell)%pref
          si   = si  /this%aux_vars_in(icell)%pref

          b_p(icell) = b_p(icell) + &
               auxvar(icell)%soil_rhg*(qsat - si*auxvar(icell)%temperature)
       else

#ifdef USE_BONAN_FORMULATION
          b_p(icell) = b_p(icell) + &
               auxvar(icell)%rhomol / this%dtime * auxvar(icell)%water_vapor * this%mesh%vol(icell)
#else
          b_p(icell) = b_p(icell) + &
               auxvar(icell)%rhomol / this%dtime * auxvar(icell)%water_vapor
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
    PetscReal                               :: dist_up, dist_dn, dist
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

#ifdef USE_BONAN_FORMULATION
             b_p(cell_id) = b_p(cell_id) + &
                  this%aux_vars_conn_bc(sum_conn)%ga * auxvar(sum_conn)%water_vapor
#else
             b_p(cell_id) = b_p(cell_id) + &
                  this%aux_vars_conn_bc(sum_conn)%ga/dist * auxvar(sum_conn)%water_vapor
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
    use MultiPhysicsProbConstants , only : HVAP
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
    PetscInt                             :: icell, ileaf, iconn, sum_conn
    PetscInt                             :: cell_id, cell_id_up, cell_id_dn
    PetscInt                             :: row, col
    PetscReal                            :: value, ga, dist, gsw, gsa, gs0
    PetscReal                            :: qsat, si, gleaf, gleaf_et
    class(connection_set_type) , pointer :: cur_conn_set
    type(condition_type)       , pointer :: cur_cond

    ! For soil cell
    icell = 1
    ! For interior and top cell
    do icell = 1, this%mesh%ncells_local
       row = icell-1; col = icell-1

       if (this%aux_vars_in(icell)%is_soil) then
          value = 1.d0

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
             gsw = 1.d0 / this%aux_vars_in(cell_id)%soil_resis * this%aux_vars_in(cell_id)%rhomol
             gsa = this%aux_vars_conn_in(iconn)%ga
             gs0 = gsw * gsa / (gsw + gsa)
#ifdef USE_BONAN_FORMULATION
             value = gs0
#else
             value = gs0/dist
#endif
          else
             ga = this%aux_vars_conn_in(iconn)%ga
#ifdef USE_BONAN_FORMULATION
             value = ga
#else
             value = ga/dist
#endif
          end if

          call MatSetValuesLocal(B, 1, cell_id_up-1, 1, cell_id_dn-1, -value, ADD_VALUES, ierr); CHKERRQ(ierr)
          call MatSetValuesLocal(B, 1, cell_id_up-1, 1, cell_id_up-1,  value, ADD_VALUES, ierr); CHKERRQ(ierr)

          call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_up-1, -value, ADD_VALUES, ierr); CHKERRQ(ierr)
          call MatSetValuesLocal(B, 1, cell_id_dn-1, 1, cell_id_dn-1,  value, ADD_VALUES, ierr); CHKERRQ(ierr)

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
    PetscInt :: icell, ileaf, row, col
    PetscReal :: value, qsat, si
    PetscReal :: gleaf, gleaf_et

    select case (itype_of_other_goveq)
       case (GE_CANOPY_AIR_TEMP)
          ! For soil-vapror <--> air temperature
          do icell = 1, this%mesh%ncells_local
             if (this%aux_vars_in(icell)%is_soil) then
             
                row = icell-1; col = icell-1
                call SatVap(this%aux_vars_in(icell)%temperature, qsat, si)
                qsat = qsat/this%aux_vars_in(icell)%pref
                si   = si  /this%aux_vars_in(icell)%pref

                value = - this%aux_vars_in(icell)%soil_rhg * si

                call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
             end if
          end do

    case (GE_CANOPY_LEAF_TEMP)
       ileaf = rank_of_other_goveq - 2
       do icell = 1, this%mesh%ncells_local
          row = icell-1; col = icell-1

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
             value = -si*gleaf_et
#else
             value = -si*gleaf_et/this%mesh%vol(icell)
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
