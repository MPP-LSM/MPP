module SystemOfEquationsLongwaveRadiationType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                          , only : iulog
  use mpp_abortutils                      , only : endrun
  use mpp_shr_log_mod                     , only : errMsg => shr_log_errMsg
  use SystemOfEquationsBaseType           , only : sysofeqns_base_type
  use SystemOfEquationsLongwaveRadAuxType , only : sysofeqns_longwave_auxvar_type
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

  type, public, extends(sysofeqns_base_type) :: sysofeqns_longwave_rad_type

     type (sysofeqns_longwave_auxvar_type), pointer :: aux_vars_in(:) ! Internal state.

   contains

     procedure, public :: Init                  => LongwaveRadSoeInit
     procedure, public :: AllocateAuxVars       => LongwaveRadSoeAllocateAuxVars
     procedure, public :: PreSolve              => LongwaveRadSoePreSolve
     procedure, public :: ComputeRHS            => LongwaveRadSoeComputeRhs
     procedure, public :: ComputeOperators      => LongwaveRadSoeComputeOperators

  end type sysofeqns_longwave_rad_type

contains

  !------------------------------------------------------------------------
  subroutine LongwaveRadSoeInit (this)
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
    class(sysofeqns_longwave_rad_type) :: this

    call SOEBaseInit(this)

    nullify(this%aux_vars_in)

  end subroutine LongwaveRadSoeInit

  !------------------------------------------------------------------------
  subroutine LongwaveRadSoeAllocateAuxVars (this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use GoverningEquationBaseType           , only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_longwave_rad_type) :: this
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

  end subroutine LongwaveRadSoeAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine LongwaveRadSoePreSolve (this)
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
    class(sysofeqns_longwave_rad_type) :: this
    !
    class(goveqn_base_type) , pointer :: cur_goveq

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       call cur_goveq%PreSolve()

       cur_goveq => cur_goveq%next
    enddo

  end subroutine LongwaveRadSoePreSolve

  !------------------------------------------------------------------------
  subroutine LongwaveRadSoeComputeRhs(this, ksp, B, ierr)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_longwave_rad_type) :: this
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

  end subroutine LongwaveRadSoeComputeRhs

  !------------------------------------------------------------------------
  subroutine LongwaveRadSoeComputeOperators(this, ksp, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to system-of-equations
    !
    ! !USES:
    use SystemOfEquationsBaseType, only : SOEBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_longwave_rad_type) :: this
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

    ! Initialize the matrix
    call MatZeroEntries(B, ierr); CHKERRQ(ierr)

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

  end subroutine LongwaveRadSoeComputeOperators

#endif

end module SystemOfEquationsLongwaveRadiationType

