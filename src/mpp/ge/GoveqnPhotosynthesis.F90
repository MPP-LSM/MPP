module GoveqnPhotosynthesisType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use PhotosynthesisAuxType     , only : photosynthesis_auxvar_type
  use GoverningEquationBaseType , only : goveqn_base_type
  use petscvec
  use petscmat
  use petscsys
  !
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_photosynthesis_type

     type(photosynthesis_auxvar_type), pointer ::  aux_vars_in(:)

   contains

     procedure, public :: Setup                     => PhotosynthesisSetup
     procedure, public :: AllocateAuxVars           => PhotosynthesisAllocateAuxVars
     procedure, public :: UpdateAuxVars             => PhotosynthesisUpdateAuxVars
     procedure, public :: ComputeResidual           => PhotosynthesisComputeResidual
     !procedure, public :: ComputeOperatorsDiag     => PhotosynthesisComputeOperatorsDiag
     !procedure, public :: GetRValues               => PhotosynthesisGetRValues
     procedure, public :: SavePrimaryIndependentVar => PhotosynthesisSavePrmIndepVar

  end type goveqn_photosynthesis_type

contains

  !------------------------------------------------------------------------
  subroutine PhotosynthesisSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_PHOTOSYNTHESIS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_photosynthesis_type) :: this

    call this%Create()

    this%name  = ""
    this%itype = GE_PHOTOSYNTHESIS
    this%dof   = 1

    nullify(this%aux_vars_in)

  end subroutine PhotosynthesisSetup

  !---------------------------------------------------------------------
  subroutine PhotosynthesisAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates memory for storing auxiliary variables associated with:
    !   + Internal control volumes,
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_photosynthesis_type) :: this
    !
    PetscInt                          :: ghosted_id

    ! Allocate memory and initialize aux vars: For internal connections
    allocate(this%aux_vars_in(this%mesh%ncells_all))

    ! Update aux vars for internal cells
    do ghosted_id = 1, this%mesh%ncells_all
      if (this%mesh%is_active(ghosted_id)) then
         call this%aux_vars_in(ghosted_id)%Init()
      end if
   enddo

  end subroutine PhotosynthesisAllocateAuxVars

  !------------------------------------------------------------------------

  subroutine PhotosynthesisUpdateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Updates auxiliary variable associated with internal control volumes
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_photosynthesis_type)   :: this
    !
    ! !LOCAL VARIABLES
    PetscInt                            :: ghosted_id
    PetscInt                            :: sum_conn
    PetscInt                            :: iconn

    ! Update aux vars for internal cells
    do ghosted_id = 1, this%mesh%ncells_all
       if (this%mesh%is_active(ghosted_id)) then
          call this%aux_vars_in(ghosted_id)%AuxVarCompute()
       end if
    enddo

  end subroutine PhotosynthesisUpdateAuxVars

  !------------------------------------------------------------------------
  subroutine PhotosynthesisComputeResidual(this, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the residual equation for the discretized Richards equation
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_photosynthesis_type)         :: this
    Vec                                       :: X
    Vec                                       :: F
    PetscErrorCode                            :: ierr
    !
    ! !LOCAL VARIABLES
    type(photosynthesis_auxvar_type), pointer :: avars(:)
    PetscInt                                  :: icell
    PetscReal                                 :: gleaf
    PetscReal, pointer                        :: f_p(:)

    ! F(ci) = An(ci) - gleaf(ci) * (ca - ci)

    avars => this%aux_vars_in
    
    call VecGetArrayF90(F, f_p, ierr); CHKERRQ(ierr)

    do icell = 1, this%mesh%ncells_local
    
       gleaf = 1.d0 / (1.d0/avars(icell)%gbc + 1.6d0/avars(icell)%gs);

       f_p(icell) = avars(icell)%an - gleaf * (avars(icell)%cair - avars(icell)%ci)

    end do

    call VecRestoreArrayF90(F, f_p, ierr); CHKERRQ(ierr)

  end subroutine PhotosynthesisComputeResidual

  !------------------------------------------------------------------------
  subroutine PhotosynthesisSavePrmIndepVar (this, x)
   !
   ! !DESCRIPTION:
   !
   ! !USES:
   !
   implicit none
   !
   ! !ARGUMENTS
   class(goveqn_photosynthesis_type) :: this
   Vec :: x
   !
   PetscScalar, pointer :: x_p(:)
   PetscInt             :: ghosted_id, size
   PetscErrorCode       :: ierr
   
   call VecGetLocalSize(x, size, ierr); CHKERRQ(ierr)

   if (size /= this%mesh%ncells_local) then
      call endrun(msg="ERROR size of vector /= number of cells in the mesh "//errmsg(__FILE__, __LINE__))
   end if

   call VecGetArrayReadF90(x, x_p, ierr); CHKERRQ(ierr)

   do ghosted_id = 1, this%mesh%ncells_local
      this%aux_vars_in(ghosted_id)%ci = x_p(ghosted_id)
   end do

   call VecRestoreArrayReadF90(x, x_p, ierr); CHKERRQ(ierr)

 end subroutine PhotosynthesisSavePrmIndepVar

#endif

end module GoveqnPhotosynthesisType
