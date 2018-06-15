module GoverningEquationBaseType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Govneqn data type allocation
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
  use petscvec
  use petscmat

  ! !USES:
  use mpp_varctl         , only : iulog
  use mpp_abortutils     , only : endrun
  use mpp_shr_log_mod    , only : errMsg => shr_log_errMsg
  use MeshType           , only : mesh_type
  use ConditionType      , only : condition_list_type
  use CouplingVariableType, only : coupling_variable_list_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: goveqn_base_type

     character (len=256)             :: name                                 ! name of governing equation (GE)
     PetscInt                        :: itype                                ! identifier
     PetscInt                        :: rank_in_soe_list                     ! rank of governing equation in SoE list
     class(mesh_type),pointer        :: mesh                                 ! pointer to the mesh
     PetscInt                        :: mesh_itype                           ! type of mesh
     PetscInt                        :: mesh_rank                            ! rank of mesh within the MPP mesh list
     type(condition_list_type)       :: boundary_conditions                  ! boundary conditions to the GE
     type(condition_list_type)       :: source_sinks                         ! source/sinks in the GE

     PetscReal                       :: dtime                                ! time step [sec]

     ! Track variables supplied by other governing equations.
     type(coupling_variable_list_type) :: coupling_vars

     class(goveqn_base_type), pointer :: next
   contains
     procedure, public :: Create                            => GoveqnBaseCreate
     procedure, public :: Destroy                           => GoveqnBaseDestroy
     procedure, public :: PrintInfo                         => GoveqnBasePrintInfo
     procedure, public :: PreSolve                          => GoveqnBasePreSolve
     procedure, public :: ComputeResidual                   => GoveqnBaseComputeResidual
     procedure, public :: ComputeJacobian                   => GoveqnBaseComputeJacobian
     procedure, public :: ComputeOffDiagJacobian            => GoveqnBaseComputeOffDiagJacobian
     procedure, public :: ComputeRHS                        => GoveqnBaseComputeRHS
     procedure, public :: ComputeOperatorsDiag              => GoveqnBaseComputeOperatorsDiag
     procedure, public :: ComputeOperatorsOffDiag           => GoveqnBaseComputeOperatorsOffDiag
     procedure, public :: SetDtime                          => GoveqnBaseSetDtime
     procedure, public :: GetNConditionsExcptCondItype      => GoveqnBaseGetNConditionsExcptCondItype
     procedure, public :: GetNCellsInCondsExcptCondItype    => GoveqnBaseGetNCellsInCondsExcptCondItype
     procedure, public :: GetCondNamesExcptCondItype        => GoveqnBaseGetCondNamesExcptCondItype
     procedure, public :: AddCondition                      => GoveqnBaseAddCondition
     procedure, public :: AddCouplingBC                     => GoveqnBaseAddCouplingBC
     procedure, public :: UpdateConditionConnSet            => GoveqnBaseUpdateConditionConnSet
     procedure, public :: GetMeshIType                      => GoveqnBaseGetMeshIType
     procedure, public :: GetMeshGridCellIsActive           => GoveqnBaseGetMeshGridCellIsActive
     procedure, public :: GetConnIDDnForCondsExcptCondItype => GoveqnBaseGetConnIDDnForCondsExcptCondItype
     procedure, public :: IsCoupledToOtherEquation          => GoveqnIsCoupledToOtherEquation
  end type goveqn_base_type
  !------------------------------------------------------------------------

  public :: GoveqnBasePrintInfo

contains

  !------------------------------------------------------------------------
  subroutine GoveqnBaseCreate(this)
    !
    ! !DESCRIPTION:
    ! Initialze a GE object
    !
    ! !USES:
    use CouplingVariableType , only : CouplingVariableListCreate
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    this%name             = ""
    this%itype            = -1
    this%rank_in_soe_list = -1
    this%mesh_itype       = 0
    this%mesh_rank        = 0
    this%dtime            = 0.d0

    nullify(this%mesh)
    call this%boundary_conditions%Init()
    call this%source_sinks%Init()

    call CouplingVariableListCreate(this%coupling_vars)

    nullify(this%next                               )

  end subroutine GoveqnBaseCreate

  !------------------------------------------------------------------------
  subroutine GoveqnBaseDestroy(this)
    !
    ! !DESCRIPTION:
    ! Release allocated memory
    !
    ! !USES:
    use CouplingVariableType , only : CouplingVariableListDestroy
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    call this%boundary_conditions%Clean()
    call this%source_sinks%Clean()
    call CouplingVariableListDestroy (this%coupling_vars       )

  end subroutine GoveqnBaseDestroy

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeResidual(this, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc SNES Function evaluation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: X
    Vec                     :: F
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseResidual must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeResidual

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeOffDiagJacobian(this, X_1, X_2, A, B, &
       itype_of_other_goveq, &
       rank_of_other_goveq, &
       ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc SNES Jacobian corresponding to off-diagonal
    ! matrix coupling between two GEs
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: X_1
    Vec                     :: X_2
    Mat                     :: A
    Mat                     :: B
    PetscInt                :: itype_of_other_goveq
    PetscInt                :: rank_of_other_goveq
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseJacobianOffDiag must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeOffDiagJacobian

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeJacobian(this, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc SNES Jacobian
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: X
    Mat                     :: A
    Mat                     :: B
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseJacobian must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeJacobian

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeRHS(this, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Vec                     :: B
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseComputeRHS must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeRHS

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeOperatorsDiag(this, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc KSP Operator matrix
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Mat                     :: A
    Mat                     :: B
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseComputeOperatorsDiag must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeOperatorsDiag

  !------------------------------------------------------------------------
  subroutine GoveqnBaseComputeOperatorsOffDiag(this, A, B, &
       itype_of_other_goveq, rank_of_other_goveq, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc KSP Operator matrix
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    Mat                     :: A
    Mat                     :: B
    PetscInt                :: itype_of_other_goveq
    PetscInt                :: rank_of_other_goveq
    PetscErrorCode          :: ierr

    write(iulog,*)'GoveqnBaseComputeOperatorsOffDiag must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBaseComputeOperatorsOffDiag

  !------------------------------------------------------------------------
  subroutine GoveqnBaseSetDtime(this, dtime)
    !
    ! !DESCRIPTION:
    ! Sets time step
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    PetscReal               :: dtime

    this%dtime = dtime

  end subroutine GoveqnBaseSetDtime

  !------------------------------------------------------------------------
  subroutine GoveqnBasePrintInfo(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use CouplingVariableType , only : CouplingVariableListPrintInfo
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    write(iulog,*)'    ---------------------------------------------------------'
    write(iulog,*)'    Goveqn_name       : ',trim(this%name)
    write(iulog,*)'    Goveqn_id         : ',this%itype
    write(iulog,*)'    Goveqn_mesh_itype : ',this%mesh_itype
    write(iulog,*)' '

    write(iulog,*)'    BC'
    call this%boundary_conditions%PrintInfo()
    write(iulog,*)'    SS'
    call this%source_sinks%PrintInfo()
    write(iulog,*)'    Coupling Vars'
    call CouplingVariableListPrintInfo(this%coupling_vars)

  end subroutine GoveqnBasePrintInfo

  !------------------------------------------------------------------------
  subroutine GoveqnBasePreSolve(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PreSolve
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this

    write(iulog,*) 'In GoveqnBasePreSolve: This needs to be extended by ' // &
         'child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine GoveqnBasePreSolve

  !------------------------------------------------------------------------
  subroutine GoveqnBaseGetNConditionsExcptCondItype(this, cond_type, &
              cond_itype_to_exclude, num_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)            :: this
    PetscInt             , intent(in)  :: cond_type
    PetscInt             , intent(in)  :: cond_itype_to_exclude
    PetscInt             , intent(out) :: num_conds
    !
    ! !LOCAL VARIABLES:
    character(len=256)                 :: string

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       call this%boundary_conditions%GetNumCondsExcptCondItype( &
            cond_itype_to_exclude, num_conds)

    case (COND_SS)
       call this%source_sinks%GetNumCondsExcptCondItype( &
            cond_itype_to_exclude, num_conds)

    case default
       write(string,*) cond_type
       write(iulog,*) 'Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine GoveqnBaseGetNConditionsExcptCondItype

  !------------------------------------------------------------------------
  subroutine GoveqnBaseGetNCellsInCondsExcptCondItype(this, cond_type, &
              cond_itype_to_exclude, num_conds, ncells_for_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)                     :: this
    PetscInt             , intent(in)           :: cond_type
    PetscInt             , intent(in)           :: cond_itype_to_exclude
    PetscInt             , intent(out)          :: num_conds
    PetscInt             , intent(out), pointer :: ncells_for_conds(:)
    !
    ! !LOCAL VARIABLES:
    character(len=256)                          :: string

    call this%GetNConditionsExcptCondItype( &
         cond_type, cond_itype_to_exclude, num_conds)

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       call this%boundary_conditions%GetNumCellsForCondsExcptCondItype( &
            cond_itype_to_exclude, ncells_for_conds)

    case (COND_SS)
       call this%source_sinks%GetNumCellsForCondsExcptCondItype( &
            cond_itype_to_exclude, ncells_for_conds)

    case default
       write(string,*) cond_type
       write(iulog,*) 'Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine GoveqnBaseGetNCellsInCondsExcptCondItype

  !------------------------------------------------------------------------
  subroutine GoveqnBaseGetCondNamesExcptCondItype(this, cond_type, &
                cond_itype_to_exclude, num_conds, cond_names)
    !
    ! !DESCRIPTION:
    ! Returns the total number and names of conditions (eg. boundary condition
    ! or source-sink) present.
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)            :: this
    PetscInt             , intent(in)  :: cond_type
    PetscInt             , intent(in)  :: cond_itype_to_exclude
    PetscInt             , intent(out) :: num_conds
    character (len=256)  , pointer     :: cond_names(:)
    type(condition_type) , pointer     :: cur_cond
    !
    ! !LOCAL VARIABLES
    character(len=256)                 :: string

    call this%GetNConditionsExcptCondItype( &
         cond_type, cond_itype_to_exclude, num_conds)

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       call this%boundary_conditions%GetCondNamesExcptCondItype( &
            cond_itype_to_exclude, cond_names)

    case (COND_SS)
       call this%source_sinks%GetCondNamesExcptCondItype( &
            cond_itype_to_exclude, cond_names)

    case default
       write(string,*) cond_type
       write(iulog,*) 'Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine GoveqnBaseGetCondNamesExcptCondItype

  !------------------------------------------------------------------------
  subroutine GoveqnBaseAddCondition(this, ss_or_bc_type, name, unit,    &
       cond_type, region_type, conn_set)
    !
    ! !DESCRIPTION:
    ! Adds boundary/source-sink condition to governing equation
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConditionType             , only : ConditionNew
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)                     :: this
    PetscInt                                    :: ss_or_bc_type
    character(len =*)                           :: name
    character(len =*)                           :: unit
    PetscInt                                    :: cond_type
    PetscInt, optional                          :: region_type
    type(connection_set_type),pointer, optional :: conn_set
    !
    type(condition_type), pointer :: cond

    cond => ConditionNew()

    call ConditionMinimalSetup(cond, this%mesh, name, unit, cond_type, &
         region_type, conn_set)

    select case (ss_or_bc_type)
    case (COND_BC)
       if (cond_type == COND_DIRICHLET_FRM_OTR_GOVEQ) then
          write(iulog,*) 'Call AddCouplingBC'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif
       call this%boundary_conditions%AddCondition(cond)

    case (COND_SS)
       call this%source_sinks%AddCondition(cond)

    case default
       write(iulog,*) 'Unknown condition type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine GoveqnBaseAddCondition

  !------------------------------------------------------------------------
  subroutine GoveqnBaseAddCouplingBC(this, name, unit,    &
       num_other_goveqs, id_of_other_goveqs, &
       itype_of_other_goveqs, icoupling_of_other_goveqns, &
       region_type, conn_set)
    !
    ! !DESCRIPTION:
    ! Adds boundary condition to the governing equation that is
    ! used to couple it with another governing equation
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConditionType             , only : ConditionNew
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)                                      :: this
    character(len =*)         , intent(in)                       :: name
    character(len =*)         , intent(in)                       :: unit
    PetscInt                  , intent(in)                       :: num_other_goveqs
    PetscInt                  , intent(in), pointer              :: id_of_other_goveqs(:)
    PetscInt                  , intent(in), pointer              :: itype_of_other_goveqs(:)
    PetscBool                 , intent(in), pointer, optional    :: icoupling_of_other_goveqns(:)
    PetscInt                  , intent(in), optional             :: region_type
    type(connection_set_type) , intent(inout), pointer, optional :: conn_set
    !
    type(condition_type)      , pointer                          :: cond

    cond => ConditionNew()

    call ConditionMinimalSetup(cond, this%mesh, name, unit, COND_DIRICHLET_FRM_OTR_GOVEQ, &
         region_type, conn_set)

    cond%num_other_goveqs = num_other_goveqs

    allocate (cond%rank_of_other_goveqs                     (cond%num_other_goveqs))
    allocate (cond%itype_of_other_goveqs                    (cond%num_other_goveqs))
    allocate (cond%swap_order_of_other_goveqs               (cond%num_other_goveqs))
    allocate (cond%is_the_other_GE_coupled_via_int_auxvars  (cond%num_other_goveqs))

    cond%rank_of_other_goveqs  (:) = id_of_other_goveqs(:)
    cond%itype_of_other_goveqs (:) = itype_of_other_goveqs(:)

    if (present(icoupling_of_other_goveqns)) then
       cond%is_the_other_GE_coupled_via_int_auxvars(:) = icoupling_of_other_goveqns(:)
    endif

    call this%boundary_conditions%AddCondition(cond)

  end subroutine GoveqnBaseAddCouplingBC

  !------------------------------------------------------------------------
  subroutine ConditionMinimalSetup(cond, mesh, name, unit, cond_type, region_type, conn_set)
    !
    ! !DESCRIPTION:
    ! Does the minimal setup of a condition_type
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConditionType             , only : ConditionNew
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_type), pointer :: cond
    class(mesh_type), pointer        :: mesh
    character(len =*)                           :: name
    character(len =*)                           :: unit
    PetscInt                                    :: cond_type
    PetscInt, optional                          :: region_type
    type(connection_set_type),pointer, optional :: conn_set

    cond%name         = trim(name)
    cond%units        = trim(unit)
    cond%itype        = cond_type

    allocate(cond%conn_set)

    ! Check the optional arguments
    if (.not.present(conn_set) .and. .not.present(region_type)) then
       write(iulog,*) 'Neither region_type nor connection set is defined while adding a BC'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif
    if (present(conn_set) .and. present(region_type)) then
       write(iulog,*) 'Both region_type and connection set are defined while adding a BC.' // &
       ' Only one can be defined.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (present(region_type)) then
       cond%region_itype = region_type
       call MeshCreateConnectionSet(mesh, cond%region_itype, cond%conn_set, cond%ncells)
    else
       cond%conn_set => conn_set
       cond%ncells = conn_set%num_connections
       nullify(conn_set)
    endif

    allocate(cond%value(cond%ncells))
    cond%value(:) = 0.d0

  end subroutine ConditionMinimalSetup

  !------------------------------------------------------------------------
  subroutine GoveqnBaseUpdateConditionConnSet(this, icond, &
       ss_or_bc_type, nconn,  conn_id_up, conn_id_dn, &
       conn_dist_up, conn_dist_dn,  conn_unitvec, conn_area)
    !
    ! !DESCRIPTION:
    ! Updates connection of a source-sink/boundary condition
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use ConditionType             , only : ConditionNew
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use ConnectionSetType         , only : ConnectionSetNew
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)       :: this
    PetscInt                      :: icond
    PetscInt                      :: ss_or_bc_type
    PetscInt                      :: nconn
    PetscInt, pointer             :: conn_id_up(:)   !
    PetscInt, pointer             :: conn_id_dn(:)   !
    PetscReal, pointer            :: conn_dist_up(:) !
    PetscReal, pointer            :: conn_dist_dn(:) !
    PetscReal, pointer            :: conn_area(:)    !
    PetscReal, pointer            :: conn_type(:)    !
    PetscReal, pointer            :: conn_unitvec(:,:)
    !
    type(condition_type), pointer :: cond
    PetscInt                      :: iconn, ii
    PetscInt                      :: num_conditions


    num_conditions = 0
    select case (ss_or_bc_type)
    case (COND_BC)

       cond => this%boundary_conditions%first
       do
          if (.not.associated(cond)) exit
          num_conditions = num_conditions + 1
          cond => cond%next
       enddo

       if (icond > num_conditions) then
          write(iulog,*),'Could not find the BC whose ' // &
               'connection set needs to be updated.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       cond => this%boundary_conditions%first
       do ii = 1,icond-1
          cond => cond%next
       enddo

    case (COND_SS)
       cond => this%source_sinks%first
       do
          if (.not.associated(cond)) exit
          num_conditions = num_conditions + 1
          cond => cond%next
       enddo

       if (icond > num_conditions) then
          write(iulog,*),'Could not find the SS whose ' // &
               'connection set needs to be updated.'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       cond => this%source_sinks%first
       do ii = 1,icond-1
          cond => cond%next
       enddo

    case default
       write(iulog,*) 'Unknown condition type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

    call cond%conn_set%Destroy()

    allocate(cond%conn_set)
    cond%conn_set => ConnectionSetNew(nconn)

    do iconn = 1, nconn
       call cond%conn_set%conn(iconn)%SetIDUp  (conn_id_up(iconn))
       call cond%conn_set%conn(iconn)%SetIDDn  (conn_id_dn(iconn))
       call cond%conn_set%conn(iconn)%SetArea  (conn_area(iconn))
       call cond%conn_set%conn(iconn)%SetDistUp(conn_dist_up(iconn))
       call cond%conn_set%conn(iconn)%SetDistDn(conn_dist_dn(iconn))
       call cond%conn_set%conn(iconn)%SetDistUnitVec(conn_unitvec(iconn,1), conn_unitvec(iconn,2),conn_unitvec(iconn,3))
    enddo

  end subroutine GoveqnBaseUpdateConditionConnSet

  !------------------------------------------------------------------------
  function GoveqnBaseGetMeshIType(this)
    !
    ! !DESCRIPTION:
    ! Returns mesh itype
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    !
    PetscInt                :: GoveqnBaseGetMeshIType

    GoveqnBaseGetMeshIType = this%mesh%itype

  end function GoveqnBaseGetMeshIType

  !------------------------------------------------------------------------
  subroutine GoveqnBaseGetMeshGridCellIsActive(this, grid_cell_active)
    !
    ! !DESCRIPTION:
    ! Returns an array identifying if a grid cell is active
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type) :: this
    PetscBool, pointer      :: grid_cell_active(:)
    !
    PetscInt                :: icell

    allocate(grid_cell_active(this%mesh%ncells_all))

    do icell = 1, this%mesh%ncells_all
       if (this%mesh%is_active(icell)) then
          grid_cell_active(icell) = PETSC_TRUE
       else
          grid_cell_active(icell) = PETSC_FALSE
       end if
    end do

  end subroutine GoveqnBaseGetMeshGridCellIsActive

  !------------------------------------------------------------------------
  subroutine GoveqnBaseGetConnIDDnForCondsExcptCondItype(this, cond_type, &
              cond_itype_to_exclude, num_cells, cell_id_dn)
    !
    ! !DESCRIPTION:
    ! Returns an array with downwind cell ids associated with all conditions
    ! excluding conditions of itype = cond_type_to_exclude
    !
    ! !USES:
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)                     :: this
    PetscInt             , intent(in)           :: cond_type
    PetscInt             , intent(in)           :: cond_itype_to_exclude
    PetscInt             , intent(out)          :: num_cells
    PetscInt             , intent(out), pointer :: cell_id_dn(:)
    !
    ! !LOCAL VARIABLES:
    character(len=256)                          :: string

    ! Choose the condition type
    select case (cond_type)
    case (COND_BC)
       call this%boundary_conditions%GetConnIDDnForCondsExcptCondItype( &
            cond_itype_to_exclude, num_cells, cell_id_dn)

    case (COND_SS)
       call this%source_sinks%GetConnIDDnForCondsExcptCondItype( &
            cond_itype_to_exclude, num_cells, cell_id_dn)

    case default
       write(string,*) cond_type
       write(iulog,*) 'Unknown cond_type = ' // trim(string)
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine GoveqnBaseGetConnIDDnForCondsExcptCondItype

  !------------------------------------------------------------------------
  subroutine GoveqnIsCoupledToOtherEquation(this, rank_of_other_goveqn, is_coupled)
    !
    ! !DESCRIPTION:
    ! Returns if the given governing equation is coupled to other equation
    ! that the rank = rank_of_other_goveqn
    !
    use CouplingVariableType      , only : coupling_variable_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_base_type)               :: this
    PetscInt                              :: rank_of_other_goveqn
    !
    PetscBool                             :: is_coupled
    type(coupling_variable_type), pointer :: cpl_var

    is_coupled = PETSC_FALSE

    cpl_var => this%coupling_vars%first
    do
       if (.not.associated(cpl_var)) exit

       if (cpl_var%rank_of_coupling_goveqn == rank_of_other_goveqn) then
          is_coupled = PETSC_TRUE
          exit
       endif

       cpl_var => cpl_var%next
    enddo

  end subroutine GoveqnIsCoupledToOtherEquation

#endif

end module GoverningEquationBaseType
