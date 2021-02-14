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
     procedure, public :: ComputeJacobian           => PhotosynthesisComputeJacobian
     !procedure, public :: ComputeOperatorsDiag     => PhotosynthesisComputeOperatorsDiag
     procedure, public :: GetRValues               => PhotosynthesisGetRValues
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
         call this%aux_vars_in(ghosted_id)%Init(this%dof)
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
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_MEDLYN
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_BBERRY
    use MultiPhysicsProbConstants , only : VAR_WUE
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
    PetscInt                                  :: icell, idof
    PetscReal                                 :: wl, term1, term2, term3
    PetscReal, pointer                        :: f_p(:)

    ! F(ci) = An(ci) - gleaf(ci) * (ca - ci)

    avars => this%aux_vars_in
    
    call VecGetArrayF90(F, f_p, ierr); CHKERRQ(ierr)

    do icell = 1, this%mesh%ncells_local
    
       select case (avars(icell)%gstype)
       case (VAR_STOMATAL_CONDUCTANCE_BBERRY)
          do idof = 1,this%dof
             f_p(icell) = avars(icell)%an(idof) - avars(icell)%gleaf_c(idof) * (avars(icell)%cair - avars(icell)%ci(idof))
          end do

       case (VAR_STOMATAL_CONDUCTANCE_MEDLYN)
          do idof = 1,this%dof
             f_p(icell) = avars(icell)%an(idof) - avars(icell)%gleaf_c(idof) * (avars(icell)%cair - avars(icell)%ci(idof))
          end do

       case (VAR_WUE)
          wl = (avars(icell)%esat - avars(icell)%eair)/avars(icell)%patm

          do idof = 1,this%dof
             term1 = (avars(icell)%cair - avars(icell)%ci(idof))/wl
             term2 = avars(icell)%dan_dci(idof) / (avars(icell)%dan_dci(idof) + avars(icell)%gleaf_c(idof))
             term3 = 1.6 * (avars(icell)%gleaf_c(idof)/avars(icell)%gleaf_w(idof))**2.d0

             f_p(icell) = avars(icell)%iota - term1 * term2 * term3
          end do

       end select
    end do

    call VecRestoreArrayF90(F, f_p, ierr); CHKERRQ(ierr)

  end subroutine PhotosynthesisComputeResidual

  !------------------------------------------------------------------------
  subroutine PhotosynthesisComputeJacobian(this, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Computes the jacobian matrix for the discretized Richards equation
    !
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_MEDLYN
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_BBERRY
    use MultiPhysicsProbConstants , only : VAR_WUE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_photosynthesis_type)        :: this
    Vec                                      :: X
    Mat                                      :: A
    Mat                                      :: B
    PetscErrorCode                           :: ierr
    !
    ! !LOCAL VARIABLES
    type(photosynthesis_auxvar_type), pointer :: avars(:)
    PetscInt                                  :: icell, idof
    PetscReal                                 :: value
    PetscReal                                 :: an_1, ci_1, gleaf_1
    PetscReal                                 :: an_2, ci_2, gleaf_2
    PetscReal                                 :: wl
    PetscReal                                 :: term1_1, term2_1, term3_1
    PetscReal                                 :: term1_2, term2_2, term3_2
    PetscReal                                 :: dterm1_dci, dterm2_dci, dterm3_dci
    PetscReal, pointer                        :: f_p(:)
    PetscReal                                 :: ci_perturb

    ! F(ci) = An(ci) - gleaf(ci) * (ca - ci)
    !
    ! dF/dci = dAn/dci - dgleaf/dci * (ca - ci) + gleaf

    
    avars => this%aux_vars_in

    do icell = 1, this%mesh%ncells_local
    
       wl = (avars(icell)%esat - avars(icell)%eair)/avars(icell)%patm

       do idof = 1,this%dof
          if (avars(icell)%gstype == VAR_WUE) then
             ci_perturb = 1.d-10
          else
             ci_perturb = 1.d-14
          end if

          an_1    = avars(icell)%an(idof)
          ci_1    = avars(icell)%ci(idof)
          gleaf_1 = avars(icell)%gleaf_c(idof)

          term1_1 = (avars(icell)%cair - avars(icell)%ci(idof))/wl
          term2_1 = avars(icell)%dan_dci(idof) / (avars(icell)%dan_dci(idof) + avars(icell)%gleaf_c(idof))
          term3_1 = 1.6 * (avars(icell)%gleaf_c(idof)/avars(icell)%gleaf_w(idof))**2.d0

          avars(icell)%ci(idof) = ci_1 - ci_perturb
          call avars(icell)%AuxVarCompute()

          an_2    = avars(icell)%an(idof)
          ci_2    = avars(icell)%ci(idof)
          gleaf_2 = avars(icell)%gleaf_c(idof)

          term1_2 = (avars(icell)%cair - avars(icell)%ci(idof))/wl
          term2_2 = avars(icell)%dan_dci(idof) / (avars(icell)%dan_dci(idof) + avars(icell)%gleaf_c(idof))
          term3_2 = 1.6 * (avars(icell)%gleaf_c(idof)/avars(icell)%gleaf_w(idof))**2.d0

          avars(icell)%ci(idof) = ci_1
          call avars(icell)%AuxVarCompute()

          select case (avars(icell)%gstype)
          case (VAR_STOMATAL_CONDUCTANCE_BBERRY)

             value = (an_1 - an_2)/ci_perturb &
                  - (gleaf_1 - gleaf_2)/ci_perturb * (avars(icell)%cair - ci_1) &
                  + gleaf_1

          case (VAR_STOMATAL_CONDUCTANCE_MEDLYN)
             value = (an_1 - an_2)/ci_perturb &
                  - (gleaf_1 - gleaf_2)/ci_perturb * (avars(icell)%cair - ci_1) &
                  + gleaf_1

          case (VAR_WUE)

             dterm1_dci = (term1_1 - term1_2)/ci_perturb
             dterm2_dci = (term2_1 - term2_2)/ci_perturb
             dterm3_dci = (term3_1 - term3_2)/ci_perturb

             !f_p = iota - term1 * term2 * term3
             value = &
                  - dterm1_dci * term2_1    * term3_1    &
                  - term1_1    * dterm2_dci * term3_1    &
                  - term1_1    * term2_1    * dterm3_dci

          end select


          call MatSetValuesLocal(B, 1, icell-1, 1, icell-1, value, ADD_VALUES, ierr); CHKERRQ(ierr)
       enddo
    end do

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    if ( A /= B) then
       call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
       call MatAssemblyEnd(  A, MAT_FINAL_ASSEMBLY, ierr); CHKERRQ(ierr)
    endif

  end subroutine PhotosynthesisComputeJacobian

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
   PetscInt             :: ghosted_id, size, idof
   PetscErrorCode       :: ierr
   
   call VecGetLocalSize(x, size, ierr); CHKERRQ(ierr)

   if (size /= this%mesh%ncells_local * this%dof) then
      call endrun(msg="ERROR size of vector /= number of cells in the mesh "//errmsg(__FILE__, __LINE__))
   end if

   call VecGetArrayReadF90(x, x_p, ierr); CHKERRQ(ierr)

   do ghosted_id = 1, this%mesh%ncells_local
      do idof = 1, this%dof
         this%aux_vars_in(ghosted_id)%ci(idof) = x_p((ghosted_id-1)*this%dof + idof)
      end do
   end do

   call VecRestoreArrayReadF90(x, x_p, ierr); CHKERRQ(ierr)

 end subroutine PhotosynthesisSavePrmIndepVar

  !------------------------------------------------------------------------
  subroutine PhotosynthesisGetRValues (this, auxvar_type, var_type, nauxvar, var_values)
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
    class(goveqn_photosynthesis_type) :: this
    PetscInt                     :: auxvar_type
    PetscInt                     :: var_type
    PetscInt                     :: nauxvar
    PetscReal, pointer           :: var_values(:)

    if (nauxvar > this%mesh%ncells_all * this%dof) then
      call endrun(msg="ERROR nauxvar exceeds the number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    endif

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call PhotosynthesisGetRValuesFromAuxVars(this%aux_vars_in, var_type, nauxvar, this%dof, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
    
  end subroutine PhotosynthesisGetRValues

  !------------------------------------------------------------------------
  subroutine PhotosynthesisGetRValuesFromAuxVars (aux_var, var_type, ncells, ndof, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_STOMATAL_CONDUCTANCE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(photosynthesis_auxvar_type) , pointer :: aux_var(:)
    PetscInt                                   :: var_type
    PetscInt                                   :: ncells
    PetscInt                                   :: ndof
    PetscReal                        , pointer :: var_values(:)
    !
    PetscInt                                   :: ghosted_id, idof

    select case(var_type)
    case(VAR_STOMATAL_CONDUCTANCE)

       do ghosted_id = 1, ncells
          do idof = 1, ndof
             var_values(ghosted_id) = aux_var(ghosted_id)%gs(idof)
          end do
       end do

    case default
       write(iulog,*) 'PhotosynthesisGetRValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine PhotosynthesisGetRValuesFromAuxVars

#endif

end module GoveqnPhotosynthesisType
