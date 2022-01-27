module SystemOfEquationsPhotosynthesisType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                                , only : iulog
  use mpp_abortutils                            , only : endrun
  use mpp_shr_log_mod                           , only : errMsg => shr_log_errMsg
  use SystemOfEquationsBaseType                 , only : sysofeqns_base_type
  use SystemOfEquationsPhotosynthesisAuxType , only : sysofeqns_photosynthesis_auxvar_type
  use GoverningEquationBaseType
  use SystemOfEquationsBaseType
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscdm
  use petscdmda
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  private

  type, public, extends(sysofeqns_base_type) :: sysofeqns_photosynthesis_type

     type (sysofeqns_photosynthesis_auxvar_type), pointer :: aux_vars_in(:) ! Internal state.

   contains

     procedure, public :: Init                  => PhotosynthesisSoeInit
     procedure, public :: AllocateAuxVars       => PhotosynthesisSoeAllocateAuxVars
     procedure, public :: PreSolve              => PhotosynthesisSoePreSolve
     procedure, public :: PostSolve             => PhotosynthesisSoePostSolve
     procedure, public :: Residual              => PhotosynthesisSoeResidual
     procedure, public :: Jacobian              => PhotosynthesisSoeJacobian

  end type sysofeqns_photosynthesis_type

contains

  !------------------------------------------------------------------------
  subroutine PhotosynthesisSoeInit (this)
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
    class(sysofeqns_photosynthesis_type) :: this

    call SOEBaseInit(this)

    nullify(this%aux_vars_in)

  end subroutine PhotosynthesisSoeInit


  !------------------------------------------------------------------------
  subroutine PhotosynthesisSoeAllocateAuxVars (this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_photosynthesis_type) :: this
    !
    class(goveqn_base_type)    , pointer :: cur_goveq

    call this%ComputeNumInternalAuxVars()

    allocate(this%aux_vars_in(this%num_auxvars_in))

  end subroutine PhotosynthesisSoeAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine PhotosynthesisSoePreSolve (this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : SOEBaseInit
    use MultiPhysicsProbConstants , only : SOE_PHOTOSYNTHESIS
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnPhotosynthesisType  , only : goveqn_photosynthesis_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_photosynthesis_type) :: this
    !
    class(goveqn_base_type),pointer :: cur_goveq
    PetscErrorCode :: ierr

  end subroutine PhotosynthesisSoePreSolve

  !------------------------------------------------------------------------
  subroutine PhotosynthesisSoePostSolve (this)
    !
    ! !DESCRIPTION:
    ! Initializes module variables and data structures
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : SOEBaseInit
    use MultiPhysicsProbConstants , only : SOE_PHOTOSYNTHESIS
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnPhotosynthesisType  , only : goveqn_photosynthesis_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_photosynthesis_type) :: this
    !
    class(goveqn_base_type),pointer :: cur_goveq
    PetscErrorCode :: ierr

    call VecCopy(this%solver%soln, this%solver%soln_prev, ierr); CHKERRQ(ierr)
    call this%SavePrimaryIndependentVar(this%solver%soln)

    select case (this%itype)
    case(SOE_PHOTOSYNTHESIS)

       cur_goveq => this%goveqns
       do
          if (.not.associated(cur_goveq)) exit
          select type(cur_goveq)
          class is (goveqn_photosynthesis_type)

             call cur_goveq%PostSolve()

          end select
          cur_goveq => cur_goveq%next
       enddo

    case default
       write(iulog,*) 'PhotosysnthesisSoE: Unknown soe_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine PhotosynthesisSoePostSolve

  !------------------------------------------------------------------------
  subroutine PhotosynthesisSoeResidual(this, snes, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Performs residual function evaluation for the VSFM
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnPhotosynthesisType , only : goveqn_photosynthesis_type
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_photosynthesis_type)      :: this
    SNES                            :: snes
    Vec                             :: X
    Vec                             :: F
    PetscErrorCode                  :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt                        :: dm_id
    PetscInt                        :: nDM
    DM, pointer                     :: dms(:)
    Vec, pointer                    :: X_subvecs(:)
    Vec, pointer                    :: F_subvecs(:)
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt                        :: row, col
    PetscViewer                     :: viewer
    character(len=256)              :: string

    ! 1) soln  ---> ge%aux_vars_in()
    call this%SavePrimaryIndependentVar(X)

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr); CHKERRQ(ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%solver%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(    nDM))
    allocate(F_subvecs(    nDM))

    ! Get vectors (X,F) for individual GEs
    call DMCompositeGetAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, X_subvecs, &
         ierr); CHKERRQ(ierr)
    call DMCompositeGetAccessArray(this%solver%dm, F, nDM, PETSC_NULL_INTEGER, F_subvecs, &
         ierr); CHKERRQ(ierr)

    ! 2) Update AuxVars
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       cur_goveq%nstep = this%nstep
       select type(cur_goveq)
       class is (goveqn_photosynthesis_type)

          call cur_goveq%UpdateAuxVars()

       end select

       cur_goveq => cur_goveq%next
    enddo

    ! Call Residual
    dm_id = 0
    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       dm_id = dm_id + 1

       call VecZeroEntries(F_subvecs(dm_id), ierr); CHKERRQ(ierr)

       call cur_goveq%ComputeResidual( &
            X_subvecs(dm_id),          &
            F_subvecs(dm_id),          &
            ierr); CHKERRQ(ierr)

       cur_goveq => cur_goveq%next
    enddo

    ! Restore vectors (u,udot,F) for individual GEs
    call DMCompositeRestoreAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)
    call DMCompositeRestoreAccessArray(this%solver%dm, F, nDM, PETSC_NULL_INTEGER, &
         F_subvecs, ierr); CHKERRQ(ierr)

    ! Free memory
    deallocate(dms)
    deallocate(X_subvecs)
    deallocate(F_subvecs)

  end subroutine PhotosynthesisSoeResidual

  !------------------------------------------------------------------------
  subroutine PhotosynthesisSoeJacobian(this, snes, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes jacobian for the VSFM
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    use MultiPhysicsProbConstants     , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_photosynthesis_type)      :: this
    SNES                            :: snes
    Vec                             :: X
    Mat                             :: A
    Mat                             :: B
    PetscErrorCode                  :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt                        :: row
    PetscInt                        :: col
    PetscInt                        :: nDM
    PetscInt                        :: icell
    PetscBool                       :: are_eqns_coupled

    IS,pointer                      :: is(:)
    DM, pointer                     :: dms(:)
    Vec, pointer                    :: X_subvecs(:)
    Mat, pointer                    :: B_submats(:,:)
    class(goveqn_base_type),pointer :: cur_goveq_1
    class(goveqn_base_type),pointer :: cur_goveq_2
    PetscViewer                     :: viewer
    character(len=256)              :: string

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%solver%dm, dms, ierr); CHKERRQ(ierr)

    ! Allocate vectors for individual GEs
    allocate(X_subvecs(    nDM))

    ! Get vectors (X) for individual GEs
    call DMCompositeGetAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Initialize the matrix
    call MatZeroEntries(B, ierr); CHKERRQ(ierr)

    ! Get submatrices
    allocate(is(nDM))
    allocate(B_submats(nDM,nDM))
    call DMCompositeGetLocalISs(this%solver%dm, is, ierr); CHKERRQ(ierr)
    do row = 1,nDM
       do col = 1,nDM
          call MatGetLocalSubMatrix(B, is(row), is(col), B_submats(row,col), &
               ierr); CHKERRQ(ierr)
       enddo
    enddo

    ! Jacobian and JacobianOffDiag
    row = 0
    cur_goveq_1 => this%goveqns
    do
       if (.not.associated(cur_goveq_1)) exit

       row = row + 1

       call cur_goveq_1%ComputeJacobian( &
            X_subvecs(row),              &
            B_submats(row,row),          &
            B_submats(row,row),          &
            ierr); CHKERRQ(ierr)

       cur_goveq_2 => cur_goveq_1%next
       col = row
       do
          if (.not.associated(cur_goveq_2)) exit

          col = col + 1

          call cur_goveq_1%IsCoupledToOtherEquation(cur_goveq_2%rank_in_soe_list, &
               are_eqns_coupled)

          if (are_eqns_coupled) then

             ! J = dF_1/dx_2
             call cur_goveq_1%ComputeOffDiagJacobian( &
                  X_subvecs(row),                     &
                  X_subvecs(col),                     &
                  B_submats(row,col),                 &
                  B_submats(row,col),                 &
                  cur_goveq_2%itype,                  &
                  cur_goveq_2%rank_in_soe_list,       &
                  ierr); CHKERRQ(ierr)

             ! J = dF_2/dx_1
             call cur_goveq_2%ComputeOffDiagJacobian( &
                  X_subvecs(col),                     &
                  X_subvecs(row),                     &
                  B_submats(col,row),                 &
                  B_submats(col,row),                 &
                  cur_goveq_1%itype,                  &
                  cur_goveq_1%rank_in_soe_list,       &
                  ierr); CHKERRQ(ierr)
          end if

          cur_goveq_2 => cur_goveq_2%next
       enddo

       cur_goveq_1 => cur_goveq_1%next
    enddo

    ! Restore vectors (X) for individual GEs
    call DMCompositeRestoreAccessArray(this%solver%dm, X, nDM, PETSC_NULL_INTEGER, &
         X_subvecs, ierr); CHKERRQ(ierr)

    ! Restore submatrices
    do row = 1,nDM
       do col = 1,nDM
          call MatRestoreLocalSubMatrix(B, is(row), is(col), B_submats(row,col), &
               ierr); CHKERRQ(ierr)
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
    deallocate(dms       )
    deallocate(X_subvecs )
    deallocate(is        )
    deallocate(B_submats )

  end subroutine PhotosynthesisSoeJacobian

#endif

end module SystemOfEquationsPhotosynthesisType
