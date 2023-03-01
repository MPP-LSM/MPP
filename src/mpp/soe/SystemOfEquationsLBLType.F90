module SystemOfEquationsLblType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                  , only : iulog
  use mpp_abortutils              , only : endrun
  use mpp_shr_log_mod             , only : errMsg => shr_log_errMsg
  use SystemOfEquationsBaseType   , only : sysofeqns_base_type
  use SystemOfEquationsLblAuxType , only : sysofeqns_lbl_auxvar_type
  use GoverningEquationBaseType
  use SystemOfEquationsBaseType
  use petscsys
  use petscksp
  use petscvec
  use petscmat
  use petscdm
  use petscdmda
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(sysofeqns_base_type) :: sysofeqns_lbl_type

     type (sysofeqns_lbl_auxvar_type), pointer :: aux_vars_in(:) ! Internal state.

   contains

    procedure, public :: Init                  => LblSoeInit
    procedure, public :: AllocateAuxVars       => LblSoeAllocateAuxVars
    procedure, public :: PreSolve              => LblSoePreSolve
    !procedure, public :: PostSolve             => MlcSoePostSolve
    procedure, public :: ComputeRHS            => LblSoeComputeRhs
    procedure, public :: ComputeOperators      => LblSoeComputeOperators

 end type sysofeqns_lbl_type

contains

  !------------------------------------------------------------------------
  subroutine LblSoeInit (this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_lbl_type) :: this

    call SOEBaseInit(this)

    nullify(this%aux_vars_in)

  end subroutine LblSoeInit

  !------------------------------------------------------------------------
  subroutine LblSoeAllocateAuxVars (this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use GoverningEquationBaseType           , only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_lbl_type) :: this
    !
    class(goveqn_base_type)    , pointer :: cur_goveq

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       this%num_auxvars_in = this%num_auxvars_in + &
            cur_goveq%mesh%ncells_all

       cur_goveq => cur_goveq%next
    enddo

    ! Allocate memory
    allocate(this%aux_vars_in(this%num_auxvars_in))

  end subroutine LblSoeAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine LblSoePreSolve (this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_lbl_type) :: this
    !
    class(goveqn_base_type) , pointer :: cur_goveq

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       call cur_goveq%PreSolve()

       cur_goveq => cur_goveq%next
    enddo

  end subroutine LblSoePreSolve

  !------------------------------------------------------------------------
  subroutine LblSoeComputeRhs(this, ksp, B, ierr)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_TEMP
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_VAPOR
    use MultiPhysicsProbConstants, only : GE_CANOPY_LEAF_TEMP
    use GoverningEquationBaseType, only : goveqn_base_type
    use CanopyTurbulence         , only : ObukhovLength
    use CanopyTurbulence         , only : WindProfile
    use CanopyTurbulence         , only : AerodynamicConductances
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_lbl_type) :: this
    KSP                       :: ksp
    Vec                       :: B
    PetscErrorCode            :: ierr
    !
    class(goveqn_base_type) , pointer :: cur_goveq
    PetscInt                          :: offset, ii
    PetscInt                          :: nDM
    DM                      , pointer :: dms(:)
    Vec                     , pointer :: B_subvecs(:)
    PetscReal               , pointer :: v_p(:)
    class(goveqn_base_type), pointer :: cur_goveqn, cur_goveq_1, cur_goveq_2
    PetscInt :: row, col

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(B_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(this%solver%dm, &
         B, nDM, &
         PETSC_NULL_INTEGER, B_subvecs, ierr)

    ! 1) GE_1 <---> GE_2 exchange AuxVars()
    !do row = 1,nDM
    !   do col = row+1,nDM
    !      call this%SetPointerToIthGovEqn(row, cur_goveq_1)
    !      call this%SetPointerToIthGovEqn(col, cur_goveq_2)
    !      call GovEqnExchangeAuxVars(this, cur_goveq_1, cur_goveq_2)
    !      call GovEqnExchangeAuxVars(this, cur_goveq_2, cur_goveq_1)
    !   enddo
    !enddo

    ii = 0
    cur_goveqn => this%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       ii = ii + 1

       call VecZeroEntries(B_subvecs(ii), ierr); CHKERRQ(ierr)
       call cur_goveqn%ComputeRhs(B_subvecs(ii), ierr); CHKERRQ(ierr)
       cur_goveqn => cur_goveqn%next
    enddo

    call DMCompositeRestoreAccessArray(this%solver%dm, &
         this%solver%soln, nDM, &
         PETSC_NULL_INTEGER, B_subvecs, ierr)

    deallocate(dms)
    deallocate(B_subvecs)

  end subroutine LblSoeComputeRhs

  !------------------------------------------------------------------------
  subroutine LblSoeComputeOperators(this, ksp, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_TEMP
    use MultiPhysicsProbConstants, only : GE_CANOPY_AIR_VAPOR
    use MultiPhysicsProbConstants, only : GE_CANOPY_LEAF_TEMP
    use GoverningEquationBaseType, only : goveqn_base_type
    use CanopyTurbulence         , only : ObukhovLength
    use CanopyTurbulence         , only : WindProfile
    use CanopyTurbulence         , only : AerodynamicConductances
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_lbl_type) :: this
    KSP                       :: ksp
    Mat                       :: A, B
    PetscErrorCode            :: ierr

    !
    ! !LOCAL VARIABLES:
    PetscInt                          :: row, col
    PetscInt                          :: rank_1, rank_2
    PetscInt                          :: nDM

    IS                      , pointer :: is(:)
    DM                      , pointer :: dms(:)
    Mat                     , pointer :: B_submats(:,:)
    class(goveqn_base_type) , pointer :: cur_goveq_1
    class(goveqn_base_type) , pointer :: cur_goveq_2

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%solver%dm, dms, ierr); CHKERRQ(ierr)

    ! Get submatrices
    allocate(is(nDM))
    allocate(B_submats(nDM,nDM))
    call DMCompositeGetLocalISs(this%solver%dm, is, ierr); CHKERRQ(ierr)
    do row = 1,nDM
      do col = 1,nDM
        call MatGetLocalSubMatrix(B, is(row), is(col), B_submats(row,col), ierr); CHKERRQ(ierr)
      enddo
    enddo

    ! Operators and OperatorsOffDiag
    row = 0
    cur_goveq_1 => this%goveqns
    do
       if (.not.associated(cur_goveq_1)) exit

       row = row + 1

       call cur_goveq_1%ComputeOperatorsDiag(B_submats(row,row),    &
                                             B_submats(row,row),    &
                                             ierr); CHKERRQ(ierr)

       cur_goveq_2 => cur_goveq_1%next
       col = row
       do
          if (.not.associated(cur_goveq_2)) exit

          col = col + 1

          rank_1 = cur_goveq_1%rank_in_soe_list
          rank_2 = cur_goveq_2%rank_in_soe_list

          call cur_goveq_1%ComputeOperatorsOffDiag( &
               B_submats(row,col),                  &
               B_submats(row,col),                  &
               cur_goveq_2%itype,                   &
               rank_2,                              &
               ierr); CHKERRQ(ierr)

          call cur_goveq_2%ComputeOperatorsOffDiag( &
               B_submats(col,row),                  &
               B_submats(col,row),                  &
               cur_goveq_1%itype,                   &
               rank_1,                              &
               ierr); CHKERRQ(ierr)

          cur_goveq_2 => cur_goveq_2%next
       enddo

       cur_goveq_1 => cur_goveq_1%next
    enddo

    ! Restore submatrices
    do row = 1,nDM
      do col = 1,nDM
        call MatRestoreLocalSubMatrix(B, is(row), is(col), B_submats(row,col), ierr); CHKERRQ(ierr)
      enddo
    enddo

    ! Destroy IS
    do row = 1,nDM
      call ISDestroy(is(row), ierr); CHKERRQ(ierr)
    enddo

    ! Assemble matrix
    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    if ( A /= B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
       call MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    endif

    ! Free memory
    deallocate(dms         )
    deallocate(is          )
    deallocate(B_submats   )

  end subroutine LblSoeComputeOperators

#endif

end module SystemOfEquationsLblType
