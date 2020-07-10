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
     procedure, public :: Residual              => PhotosynthesisSoeResidual
     !procedure, public :: Jacobian              => PhotosynthesisSoeJacobian

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
    use SystemOfEquationsBaseType, only : SOEBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_photosynthesis_type) :: this

  end subroutine PhotosynthesisSoePreSolve

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

#endif

end module SystemOfEquationsPhotosynthesisType
