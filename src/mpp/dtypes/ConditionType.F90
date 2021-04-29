module ConditionType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use ConnectionSetType  , only : connection_set_type
  use mpp_varctl         , only : iulog
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: condition_type

     character (len=256)                 :: name                                       ! name for condition
     character (len=256)                 :: units                                      ! units

     PetscInt                            :: id                                         ! identifier of condition within the list
     PetscInt                            :: itype                                      ! identifier for type of condition
     PetscInt                            :: region_itype                               ! identifier for region
     PetscInt                            :: ncells
     PetscReal                 , pointer :: value(:)                                   ! Magnitude of the condition

     PetscBool                           :: swap_order                                 ! FALSE(default): "upwind cell  " = BC; "downwind cell" = Internal cell
                                                                                       ! TRUE          : "downwind cell" = BC; "upwind cell  " = Internal cell

     PetscInt                            :: num_other_goveqs                           ! Number of other governing equations
     PetscInt                  , pointer :: rank_of_other_goveqs(:)                    ! Rank of other governing equations
     PetscInt                  , pointer :: itype_of_other_goveqs(:)                   ! Type of other governing equations
                                                                                       ! (e.g. GE_THERM_SSW_TBASED, GE_THERM_SNOW_TBASED, etc)
     PetscBool                 , pointer :: swap_order_of_other_goveqs(:)              !
     PetscBool                 , pointer :: is_the_other_GE_coupled_via_int_auxvars(:) ! TRUE : The i-th coupling governing equation is coupled via internal auxvars
                                                                                       ! FALSE: The i-th coupling governing equation is coupled via boundary auxvars

     type(connection_set_type) , pointer :: conn_set                                   ! Applicable to BC condition type
     type(condition_type)      , pointer :: next                                       ! Pointer to next condition

     contains
       procedure, public :: PrintInfo => ConditionPrintInfo
       procedure, public :: Destroy   => ConditionDestroy
  end type condition_type

  type, public :: condition_list_type
     PetscInt                         :: num_condition_list
     type(condition_type), pointer    :: first
     type(condition_type), pointer    :: last
   contains
     procedure, public :: Init                              => ConditionListInit
     procedure, public :: AddCondition                      => ConditionListAddCondition
     procedure, public :: Clean                             => ConditionListClean
     procedure, public :: PrintInfo                         => ConditionListPrintInfo
     procedure, public :: GetNumConds                       => CondListGetNumConds
     procedure, public :: GetNumCondsExcptCondItype         => CondListGetNumCondsExcptCondItype
     procedure, public :: GetNumCondsOfCondItype            => CondListGetNumCondsOfCondItype
     procedure, public :: GetNumCellsForConds               => CondListGetNumCellsForConds
     procedure, public :: GetNumCellsForCondsExcptCondItype => CondListGetNumCellsForCondsExcptCondItype
     procedure, public :: GetConnIDDnForCondsExcptCondItype => CondListGetConnIDDnForCondsExcptCondItype
     procedure, public :: GetNumCellsForCondsOfCondItype    => CondListGetNumCellsForCondsOfCondItype
     procedure, public :: GetCondNamesExcptCondItype        => CondListGetCondNamesExcptCondItype
  end type condition_list_type

  public :: ConditionNew
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  function ConditionNew()
    !
    ! !DESCRIPTION:
    ! Return a new condition
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_type),pointer :: ConditionNew
    type(condition_type),pointer :: cond

    allocate(cond)

    cond%name                   = ""
    cond%units                  = ""
    cond%id                     = -1
    cond%itype                  = -1
    cond%region_itype           = -1
    cond%ncells                 = 0

    cond%num_other_goveqs       = 0
    cond%swap_order             = PETSC_FALSE

    nullify(cond%value                                   )
    nullify(cond%rank_of_other_goveqs                    )
    nullify(cond%itype_of_other_goveqs                   )
    nullify(cond%conn_set                                )
    nullify(cond%next                                    )
    nullify(cond%swap_order_of_other_goveqs              )
    nullify(cond%is_the_other_GE_coupled_via_int_auxvars)

    ConditionNew => cond

  end function ConditionNew

  !------------------------------------------------------------------------
  subroutine ConditionListInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize an existing condition list
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type) :: this

    this%num_condition_list   = 0

    nullify(this%first)
    nullify(this%last )

  end subroutine ConditionListInit

  !------------------------------------------------------------------------
  subroutine ConditionPrintInfo(this)
    !
    ! !DESCRIPTION:
    ! Prints information about the condition
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_type) :: this

    write(iulog,*) '      Condition_name : ',trim(this%name)
    write(iulog,*) '      Condition_units: ',trim(this%units)

  end subroutine ConditionPrintInfo

  !------------------------------------------------------------------------
  subroutine ConditionListPrintInfo(this)
    !
    ! !DESCRIPTION:
    ! Prints information about all conditions present in the condition list
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in) :: this
    !
    type(condition_type), pointer         :: cur_cond

    cur_cond => this%first
    if (.not.associated(cur_cond)) return

    write(iulog,*) '    Condition-List_num : ',this%num_condition_list

    do
       if (.not.associated(cur_cond)) exit
       call cur_cond%PrintInfo()
       write(iulog,*) '      Number_of_Conn : ',cur_cond%conn_set%num_connections
       write(iulog,*)''
       cur_cond => cur_cond%next
    enddo
    write(iulog,*)''

  end subroutine ConditionListPrintInfo


  !------------------------------------------------------------------------
  subroutine ConditionListAddCondition(this, new_cond)
    !
    ! !DESCRIPTION:
    ! Add a condition to a condition list
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type)   :: this
    type(condition_type),pointer :: new_cond

    this%num_condition_list = this%num_condition_list + 1
    new_cond%id             = this%num_condition_list

    if (.not.associated(this%first)) then
       this%first => new_cond
    endif

    if (associated(this%last)) then
       this%last%next => new_cond
    endif

    this%last => new_cond

  end subroutine ConditionListAddCondition

  !------------------------------------------------------------------------
  subroutine CondListGetNumConds(this, num_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions
    !
    use MultiPhysicsProbConstants, only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in)  :: this
    PetscInt                  , intent(out) :: num_conds

    call this%GetNumCondsExcptCondItype(COND_NULL, num_conds)

  end subroutine CondListGetNumConds

  !------------------------------------------------------------------------
  subroutine CondListGetNumCondsExcptCondItype( &
       this, cond_itype, num_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions excluding conditions of
    ! itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in)  :: this
    PetscInt                  , intent(in)  :: cond_itype
    PetscInt                  , intent(out) :: num_conds
    !
    ! !LOCAL VARIABLES:
    type(condition_type), pointer            :: cur_cond

    num_conds = 0
    cur_cond => this%first
    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= cond_itype) then
          num_conds = num_conds + 1
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine CondListGetNumCondsExcptCondItype

  !------------------------------------------------------------------------
  subroutine CondListGetNumCondsOfCondItype( &
       this, cond_itype, num_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions of itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in)  :: this
    PetscInt                  , intent(in)  :: cond_itype
    PetscInt                  , intent(out) :: num_conds
    !
    ! !LOCAL VARIABLES:
    type(condition_type), pointer            :: cur_cond

    num_conds = 0
    cur_cond => this%first
    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype == cond_itype) then
          num_conds = num_conds + 1
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine CondListGetNumCondsOfCondItype

  !------------------------------------------------------------------------
  subroutine CondListGetNumCellsForConds(this, num_cells_for_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of cells associated with all conditions
    !
    use MultiPhysicsProbConstants, only : COND_NULL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in)           :: this
    PetscInt                  , intent(out), pointer :: num_cells_for_conds(:)

    call this%GetNumCellsForCondsExcptCondItype(COND_NULL, num_cells_for_conds)

  end subroutine CondListGetNumCellsForConds

  !------------------------------------------------------------------------
  subroutine CondListGetNumCellsForCondsExcptCondItype( &
       this, cond_itype, num_cells_for_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of cells associated with all conditions excluding
    ! conditions of itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in)           :: this
    PetscInt                  , intent(in)           :: cond_itype
    PetscInt                  , intent(out), pointer :: num_cells_for_conds(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                         :: num_conds
    type(condition_type)      , pointer              :: cur_cond

    call this%GetNumCondsExcptCondItype(cond_itype, num_conds)

    if (num_conds == 0) then

       nullify(num_cells_for_conds)

    else

       allocate(num_cells_for_conds(num_conds))

       num_conds = 0
       cur_cond => this%first
       do
          if (.not.associated(cur_cond)) exit
          if (cur_cond%itype /= cond_itype) then
             num_conds = num_conds + 1
             num_cells_for_conds(num_conds) = cur_cond%ncells
          endif
          cur_cond => cur_cond%next
       enddo

    end if

  end subroutine CondListGetNumCellsForCondsExcptCondItype

  !------------------------------------------------------------------------
  subroutine CondListGetNumCellsForCondsOfCondItype( &
       this, cond_itype, num_cells_for_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of cells associated with all conditions excluding
    ! conditions of itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in)           :: this
    PetscInt                  , intent(in)           :: cond_itype
    PetscInt                  , intent(out), pointer :: num_cells_for_conds(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                         :: num_conds
    type(condition_type)      , pointer              :: cur_cond

    call this%GetNumCondsOfCondItype(cond_itype, num_conds)

    if (num_conds == 0) then

       nullify(num_cells_for_conds)

    else

       allocate(num_cells_for_conds(num_conds))

       num_conds = 0
       cur_cond => this%first
       do
          if (.not.associated(cur_cond)) exit
          if (cur_cond%itype == cond_itype) then
             num_conds = num_conds + 1
             num_cells_for_conds(num_conds) = cur_cond%ncells
          endif
          cur_cond => cur_cond%next
       enddo

    end if

  end subroutine CondListGetNumCellsForCondsOfCondItype

  !------------------------------------------------------------------------
  subroutine CondListGetConnIDDnForCondsExcptCondItype( &
       this, cond_itype, num_cells, cell_id_dn)
    !
    ! !DESCRIPTION:
    ! Returns the downwind cell associated with all conditions excluding
    ! conditions of itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in)           :: this
    PetscInt                  , intent(in)           :: cond_itype
    PetscInt                  , intent(out)          :: num_cells
    PetscInt                  , intent(out), pointer :: cell_id_dn(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                         :: icond
    PetscInt                                         :: iconn
    PetscInt                                         :: count
    PetscInt                                         :: num_conds
    PetscInt                  , pointer              :: num_cells_for_conds(:)
    type(condition_type)      , pointer              :: cur_cond
    type(connection_set_type) , pointer              :: cur_conn_set

    call this%GetNumCondsExcptCondItype(cond_itype, num_conds)
    call this%GetNumCellsForCondsExcptCondItype(cond_itype, num_cells_for_conds)

    num_cells = 0
    do icond = 1, num_conds
       num_cells = num_cells + num_cells_for_conds(icond)
    end do

    if (num_cells == 0) then

       nullify(cell_id_dn)

    else

       allocate(cell_id_dn(num_cells))

       count = 0
       cur_cond => this%first
       do
          if (.not.associated(cur_cond)) exit

          if (cur_cond%itype /= cond_itype) then
             cur_conn_set => cur_cond%conn_set
             do iconn = 1, cur_conn_set%num_connections
                count = count + 1
                cell_id_dn(count) = cur_conn_set%conn(iconn)%GetIDDn()
             end do
          endif
          cur_cond => cur_cond%next
       enddo

    end if

  end subroutine CondListGetConnIDDnForCondsExcptCondItype

  !------------------------------------------------------------------------
  subroutine CondListGetCondNamesExcptCondItype( &
       this, cond_itype, cond_names)
    !
    ! !DESCRIPTION:
    ! Returns the name of all conditions except condition of
    ! itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(in) :: this
    PetscInt                  , intent(in) :: cond_itype
    character (len=256)       , pointer    :: cond_names(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                               :: num_conds
    type(condition_type)     , pointer     :: cur_cond

    call this%GetNumCondsExcptCondItype(cond_itype, num_conds)

    if (num_conds == 0) then

       nullify(cond_names)

    else

       allocate(cond_names(num_conds))

       num_conds = 0
       cur_cond => this%first
       do
          if (.not.associated(cur_cond)) exit
          if (cur_cond%itype /= cond_itype) then
             num_conds = num_conds + 1
             cond_names(num_conds) = trim(cur_cond%name)
          endif
          cur_cond => cur_cond%next
       enddo

    end if

  end subroutine CondListGetCondNamesExcptCondItype

  !------------------------------------------------------------------------
  subroutine ConditionListClean(this)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_list_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    type(condition_type), pointer            :: cur_cond, prev_cond
    
    cur_cond => this%first
    do
       if (.not.associated(cur_cond)) exit
       prev_cond => cur_cond
       cur_cond => cur_cond%next
       call prev_cond%Destroy()
    enddo

    call ConditionListInit(this)

  end subroutine ConditionListClean


  !------------------------------------------------------------------------
  subroutine ConditionDestroy(this)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(condition_type) :: this

    this%name           = ""
    this%units          = ""
    this%id             = -1
    this%itype          = -1
    this%region_itype   = -1
    this%ncells         = 0

    if (associated(this%value)) deallocate(this%value)
    nullify(this%value)

    call this%conn_set%Destroy()

    nullify(this%conn_set)

  end subroutine ConditionDestroy

#endif

end module ConditionType
