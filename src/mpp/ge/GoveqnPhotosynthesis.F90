module GoveqnPhotosynthesisType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use PhotosynthesisAuxType     , only : photosynthesis_auxvar_type
  use PhotosynthesisAuxType     , only : plant_auxvar_type
  use GoverningEquationBaseType , only : goveqn_base_type
  use petscvec
  use petscmat
  use petscsys
  use MultiPhysicsProbConstants , only : VAR_SCM_MEDLYN
  use MultiPhysicsProbConstants , only : VAR_SCM_BBERRY
  use MultiPhysicsProbConstants , only : VAR_SCM_WUE
  use MultiPhysicsProbConstants , only : VAR_SCM_MANZONI11
  use MultiPhysicsProbConstants , only : VAR_SCM_BONAN14
  use MultiPhysicsProbConstants , only : VAR_SCM_MODIFIED_BONAN14
  use MultiPhysicsProbConstants , only : VAR_SCM_OSMWANG
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
     procedure, public :: PreSolve                  => PhotosynthesisPreSolve
     procedure, public :: PostSolve                 => PhotosynthesisPostSolve

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
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_photosynthesis_type)         :: this
    Vec                                       :: X
    Vec                                       :: F
    PetscErrorCode                            :: ierr
    !
    ! !LOCAL VARIABLES
    type(photosynthesis_auxvar_type) , pointer :: avars(:)
    type(plant_auxvar_type)          , pointer :: plant
    PetscInt                                   :: icell, idof, idx, ileaf
    PetscReal                                  :: wl, term1, term2, term3, ci
    PetscReal, pointer                         :: f_p(:)

    ! F(ci) = An(ci) - gleaf(ci) * (ca - ci)

    avars => this%aux_vars_in
    
    call VecGetArrayF90(F, f_p, ierr); CHKERRQ(ierr)

    do icell = 1, this%mesh%ncells_local

       select case (avars(icell)%gstype)
       case (VAR_SCM_BBERRY)
          do idof = 1,this%dof
             idx = (icell-1)*this%dof + idof
             if ( (.not. this%mesh%is_active(icell)) .or. &
                  (.not. avars(icell)%soln_is_bounded(idof))) then
                f_p(idx) = 0.d0
             else
                if (avars(icell)%an(idof) > 0.d0) then
                   f_p(idx) = avars(icell)%an(idof) - avars(icell)%gleaf_c(idof) * (avars(icell)%cair - avars(icell)%ci(idof))
                else
                   f_p(idx) = 0.d0
                end if
             end if
          end do

       case (VAR_SCM_MEDLYN)
          do idof = 1,this%dof
             idx = (icell-1)*this%dof + idof

             if ( (.not. this%mesh%is_active(icell)) .or. &
                  (.not. avars(icell)%soln_is_bounded(idof))) then
                f_p(idx) = 0.d0
             else
                if (avars(icell)%an(idof) > 0.d0) then
                   f_p(idx) = avars(icell)%an(idof) - avars(icell)%gleaf_c(idof) * (avars(icell)%cair - avars(icell)%ci(idof))
                else
                   f_p(idx) = 0.d0
                endif
             end if
          end do

       case (VAR_SCM_WUE, VAR_SCM_MANZONI11)
          wl = (avars(icell)%esat - avars(icell)%eair)/avars(icell)%pref
          do idof = 1,this%dof
             idx = (icell-1)*this%dof + idof

             if ( (.not. this%mesh%is_active(icell)) .or. &
                  (.not. avars(icell)%soln_is_bounded(idof))) then
                f_p(idx) = 0.d0
             else
                f_p(idx) = avars(icell)%residual_wue(idof)
             end if
          end do

       case (VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14)
          wl = (avars(icell)%esat - avars(icell)%eair)/avars(icell)%pref

          do idof = 1,this%dof-1
             idx = (icell-1)*this%dof + idof

             if ( (.not. this%mesh%is_active(icell)) .or. &
                  (.not. avars(icell)%soln_is_bounded(idof))) then
                f_p(idx) = 0.d0
             else
                f_p(idx) = avars(icell)%residual_wue(idof)
             end if

          end do

          if (this%dof > 1) then
             plant => avars(icell)%plant

             idof = this%dof
             idx = (icell-1)*this%dof + idof
             ileaf = 1

             if ( (.not. this%mesh%is_active(icell)) .or. &
                  (.not. avars(icell)%soln_is_bounded(idof))) then
                f_p(idx) = 0.d0
             else
                ! psi_{t} + dpsi_{t+1} - leaf_minlwp
                f_p(idx) =  avars(icell)%residual_hyd(idof)
             end if
          end if

       case (VAR_SCM_OSMWANG)
          do idof = 1, this%dof

             idx = (icell-1)*this%dof + idof

             if ( (.not. this%mesh%is_active(icell)) .or. &
                  (.not. avars(icell)%soln_is_bounded(idof))) then
                f_p(idx) = 0.d0
             else
                f_p(idx) = avars(icell)%residual_wue(idof)
             end if

             if (this%aux_vars_in(idx)%gs(idof) < 0.d0) then
                f_p(idx) = 0.d0
                f_p(idx) = f_p(idx)/f_p(idx)
             end if
          end do
       case default
          write(*,*)'Unknown stomatal conductance model'
          call endrun(msg=errMsg(__FILE__, __LINE__))
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
    use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C3
    use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C4
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
    type(plant_auxvar_type)         , pointer :: plant
    PetscInt                                  :: icell, idof, idx, ileaf
    PetscReal                                 :: value
    PetscReal                                 :: an_1, ci_1, gleaf_1, psi_term_1
    PetscReal                                 :: an_2, ci_2, gleaf_2, psi_term_2
    PetscReal                                 :: wl
    PetscReal                                 :: term1_1, term2_1, term3_1
    PetscReal                                 :: term1_2, term2_2, term3_2
    PetscReal                                 :: dterm1_dci, dterm2_dci, dterm3_dci
    PetscReal, pointer                        :: f_p(:)
    PetscReal                                 :: ci_perturb, gs_perturb
    PetscReal                                 :: res_1, res_2, gs_1, gs_2

    ! F(ci) = An(ci) - gleaf(ci) * (ca - ci)
    !
    ! dF/dci = dAn/dci - dgleaf/dci * (ca - ci) + gleaf
    
    avars => this%aux_vars_in

    do icell = 1, this%mesh%ncells_local
    
       wl = (avars(icell)%esat - avars(icell)%eair)/avars(icell)%pref

       plant => avars(icell)%plant

       do idof = 1,this%dof
          if ( &
               (avars(icell)%gstype == VAR_SCM_BBERRY .and. (avars(icell)%c3psn == VAR_PHOTOSYNTHETIC_PATHWAY_C3)) .or. &
               (avars(icell)%gstype == VAR_SCM_MEDLYN .and. (avars(icell)%c3psn == VAR_PHOTOSYNTHETIC_PATHWAY_C3)) ) then
             ci_perturb = -1.e-7
             gs_perturb = -1.e-8
          elseif (avars(icell)%gstype == VAR_SCM_BONAN14 .or. avars(icell)%gstype == VAR_SCM_MODIFIED_BONAN14) then
             ci_perturb = -1.e-7
             gs_perturb = -1.e-8
          else
             ci_perturb = -1.e-7
             gs_perturb = -1.e-5
          endif

          if (.not. (avars(icell)%gstype == VAR_SCM_WUE .or. avars(icell)%gstype == VAR_SCM_MANZONI11)) then
             an_1    = avars(icell)%an(idof)
             ci_1    = avars(icell)%ci(idof)
             gleaf_1 = avars(icell)%gleaf_c(idof)

             term1_1 = (avars(icell)%cair - avars(icell)%ci(idof))/wl
             term2_1 = avars(icell)%dan_dci(idof) / (avars(icell)%dan_dci(idof) + avars(icell)%gleaf_c(idof))
             term3_1 = 1.6d0 * (avars(icell)%gleaf_c(idof)/avars(icell)%gleaf_w(idof))**2.d0

             ileaf = 1
             psi_term_1 = plant%leaf_psi(ileaf) + plant%dpsi_leaf(ileaf) - plant%leaf_minlwp(ileaf)

             avars(icell)%ci(idof) = ci_1 - ci_perturb
             call avars(icell)%AuxVarCompute()

             an_2    = avars(icell)%an(idof)
             ci_2    = avars(icell)%ci(idof)
             gleaf_2 = avars(icell)%gleaf_c(idof)

             term1_2 = (avars(icell)%cair - avars(icell)%ci(idof))/wl
             term2_2 = avars(icell)%dan_dci(idof) / (avars(icell)%dan_dci(idof) + avars(icell)%gleaf_c(idof))
             term3_2 = 1.6d0 * (avars(icell)%gleaf_c(idof)/avars(icell)%gleaf_w(idof))**2.d0

             psi_term_2 = plant%leaf_psi(ileaf) + plant%dpsi_leaf(ileaf) - plant%leaf_minlwp(ileaf)

             avars(icell)%ci(idof) = ci_1
             call avars(icell)%AuxVarCompute()
          endif

          select case (avars(icell)%gstype)
          case (VAR_SCM_BBERRY)

             if (avars(icell)%an(idof) > 0.d0) then
                value = (an_1 - an_2)/ci_perturb &
                     - (gleaf_1 - gleaf_2)/ci_perturb * (avars(icell)%cair - ci_1) &
                     + gleaf_1
             else
                value = 1.d0
             end if

          case (VAR_SCM_MEDLYN)
             if (avars(icell)%an(idof) > 0.d0) then
                value = (an_1 - an_2)/ci_perturb &
                     - (gleaf_1 - gleaf_2)/ci_perturb * (avars(icell)%cair - ci_1) &
                     + gleaf_1
             else
                value = 1.d0
             end if

          case (VAR_SCM_WUE, VAR_SCM_MANZONI11, VAR_SCM_OSMWANG)
             res_1 = avars(icell)%residual_wue(idof)
             gs_1  = avars(icell)%gs(idof)

             avars(icell)%gs = gs_1 - gs_perturb
             call avars(icell)%AuxVarCompute()
             res_2 = avars(icell)%residual_wue(idof)
             value = (res_1 - res_2)/gs_perturb

             avars(icell)%gs = gs_1
             call avars(icell)%AuxVarCompute()

          case (VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14)

             if (idof == 1) then
                res_1 = avars(icell)%residual_wue(idof)
                gs_1  = avars(icell)%gs(idof)

                avars(icell)%gs = gs_1 - gs_perturb
                call avars(icell)%AuxVarCompute()
                res_2 = avars(icell)%residual_wue(idof)
                value = (res_1 - res_2)/gs_perturb

                avars(icell)%gs = gs_1
                call avars(icell)%AuxVarCompute()
             else
                res_1 = avars(icell)%residual_hyd(idof)
                gs_1  = avars(icell)%gs(idof)

                avars(icell)%gs = gs_1 - gs_perturb
                call avars(icell)%AuxVarCompute()
                res_2 = avars(icell)%residual_hyd(idof)

                value = (res_1 - res_2)/gs_perturb

                avars(icell)%gs = gs_1
                call avars(icell)%AuxVarCompute()
                value = 1.d0

             endif

          end select

          idx = (icell-1)*this%dof + idof

          if ( (.not. this%mesh%is_active(icell)) .or. &
               (.not. avars(icell)%soln_is_bounded(idof))) then
             value = 1.d0
          end if

          call MatSetValuesLocal(B, 1, idx - 1, 1, idx - 1, value, ADD_VALUES, ierr); CHKERRQ(ierr)
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

      select case (this%aux_vars_in(ghosted_id)%gstype)
      case (VAR_SCM_BBERRY, VAR_SCM_MEDLYN)
         do idof = 1, this%dof
            if (this%aux_vars_in(ghosted_id)%soln_is_bounded(idof)) then
               this%aux_vars_in(ghosted_id)%ci(idof) = x_p((ghosted_id-1)*this%dof + idof)
            endif
         end do

      case (VAR_SCM_WUE, VAR_SCM_MANZONI11,VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14, VAR_SCM_OSMWANG)
         do idof = 1, this%dof
            if (this%aux_vars_in(ghosted_id)%soln_is_bounded(idof)) then
               this%aux_vars_in(ghosted_id)%gs(idof) = x_p((ghosted_id-1)*this%dof + idof)
            endif
         end do

      case default
          write(*,*)'Unknown stomatal conductance model'
          call endrun(msg=errMsg(__FILE__, __LINE__))
      end select
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
       call PhotosynthesisGetRValuesFromAuxVars(this%aux_vars_in, var_type, nauxvar, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
    
  end subroutine PhotosynthesisGetRValues

  !------------------------------------------------------------------------
  subroutine PhotosynthesisGetRValuesFromAuxVars (aux_var, var_type, ncells, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_STOMATAL_CONDUCTANCE
    use MultiPhysicsProbConstants, only : VAR_SCM_BONAN14
    use MultiPhysicsProbConstants, only : VAR_GROSS_PHOTOSYNTHESIS
    use MultiPhysicsProbConstants, only : VAR_NET_PHOTOSYNTHESIS
    !
    implicit none
    !
    ! !ARGUMENTS
    type(photosynthesis_auxvar_type) , pointer :: aux_var(:)
    PetscInt                                   :: var_type
    PetscInt                                   :: ncells
    PetscReal                        , pointer :: var_values(:)
    !
    PetscInt                                   :: ghosted_id, idof

    select case(var_type)
    case(VAR_STOMATAL_CONDUCTANCE)
       do ghosted_id = 1, ncells
          var_values(ghosted_id) = aux_var(ghosted_id)%gs_soln
       end do
    case(VAR_GROSS_PHOTOSYNTHESIS)
       do ghosted_id = 1, ncells
          var_values(ghosted_id) = aux_var(ghosted_id)%ag_soln
       end do
    case(VAR_NET_PHOTOSYNTHESIS)
       do ghosted_id = 1, ncells
          var_values(ghosted_id) = aux_var(ghosted_id)%an_soln
       end do
    case default
       write(iulog,*) 'PhotosynthesisGetRValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine PhotosynthesisGetRValuesFromAuxVars

  !------------------------------------------------------------------------
  subroutine PhotosynthesisPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_SCM_BONAN14
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_photosynthesis_type)            :: this

    ! !ARGUMENTS
    type(photosynthesis_auxvar_type) , pointer   :: avar(:)
    !
    PetscInt                                     :: icell, idof
    PetscInt                         , parameter :: idof_wue = 1
    PetscInt                         , parameter :: idof_hyd = 2

    avar => this%aux_vars_in

    do icell = 1, this%mesh%ncells_all
       call avar(icell)%PreSolve()
    end do

  end subroutine PhotosynthesisPreSolve

  !------------------------------------------------------------------------
  subroutine PhotosynthesisPostSolve(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_SCM_BONAN14
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_photosynthesis_type)            :: this

    ! !ARGUMENTS
    type(photosynthesis_auxvar_type) , pointer   :: avar(:)
    !
    PetscInt                                     :: icell, idof
    PetscInt                         , parameter :: idof_wue = 1
    PetscInt                         , parameter :: idof_hyd = 2

    avar => this%aux_vars_in

    do icell = 1, this%mesh%ncells_all

       call avar(icell)%PostSolve()
    end do

  end subroutine PhotosynthesisPostSolve

#endif

end module GoveqnPhotosynthesisType
