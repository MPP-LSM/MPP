
module SystemOfEquationsBaseType

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Base object for a system-of-equations
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  use petscsys
  use petscksp


  ! !USES:
  use mpp_abortutils            , only : endrun
  use mpp_varctl                , only : iulog
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MeshType                  , only : mesh_type
  use GoverningEquationBaseType , only : goveqn_base_type
  use SolverType                , only : solver_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: sysofeqns_base_type
     character(len =256)             :: name                         ! name for system-of-equations (SoE)
     PetscInt                        :: itype                        ! identifier for SoE

     class(goveqn_base_type),pointer :: goveqns                      ! pointer to governing equations within SoE
     PetscInt                        :: ngoveqns                     ! number of governing equations within SoE


     PetscInt                        :: mpi_rank                     ! [-]

     PetscReal                       :: time                         ! [sec]
     PetscReal                       :: dtime                        ! [sec]
     PetscReal                       :: nstep                        ! [-]

     PetscInt                        :: num_calls_to_ifunction
     PetscInt                        :: num_calls_to_ijacobian

     PetscInt, pointer               :: soe_auxvars_in_offset (:) ! Cummulative sum of number of control volumes associated with each internal cell.
     PetscInt, pointer               :: soe_auxvars_bc_offset (:) ! Cummulative sum of number of control volumes associated with each boundary condition.
     PetscInt, pointer               :: soe_auxvars_ss_offset (:) ! Cummulative sum of number of control volumes associated with each source-sink condition.
     PetscInt, pointer               :: soe_auxvars_in_ncells (:) ! Number of control volumes associated with each internal cell.
     PetscInt, pointer               :: soe_auxvars_bc_ncells (:) ! Number of control volumes associated with each boundary condition.
     PetscInt, pointer               :: soe_auxvars_ss_ncells (:) ! Number of control volumes associated with each source-sink condition.

     PetscInt                        :: num_auxvars_in            ! Number of auxvars associated with internal state.
     PetscInt                        :: num_auxvars_in_local      ! Number of auxvars associated with internal state.
     PetscInt                        :: num_auxvars_bc            ! Number of auxvars associated with boundary condition.
     PetscInt                        :: num_auxvars_ss            ! Number of auxvars associated with source-sink condition.

     type(solver_type)               :: solver

   contains
     procedure, public :: Init                         => SOEBaseInit
     procedure, public :: Clean                        => SOEBaseClean
     procedure, public :: Residual                     => SOEBaseResidual
     procedure, public :: Jacobian                     => SOEBaseJacobian
     procedure, public :: ComputeRHS                   => SOEComputeRHS
     procedure, public :: ComputeOperators             => SOEComputeOperators
     procedure, public :: StepDT                       => SOEBaseStepDT
     procedure, public :: SavePrimaryIndependentVar    => SOEBaseSavePrimaryIndependentVar
     procedure, public :: PreSolve                     => SOEBasePreSolve
     procedure, public :: PostSolve                    => SOEBasePostSolve
     procedure, public :: PreStepDT                    => SOEBasePreStepDT
     procedure, public :: PostStepDT                   => SOEBasePostStepDT
     procedure, public :: PrintInfo                    => SOEBasePrintInfo
     procedure, public :: SetPointerToIthGovEqn        => SOEBaseSetPointerToIthGovEqn
     procedure, public :: SetDtime                     => SOEBaseSetDtime
     procedure, public :: SetDataFromCLM               => SOEBaseSetDataFromCLM
     procedure, public :: GetDataForCLM                => SOEBaseGetDataForCLM
     procedure, public :: AddGovEqn                    => SOEBaseAddGovEqn
     procedure, public :: AddGovEqnWithMeshRank        => SOEBaseAddGovEqnWithMeshRank
     procedure, public :: SetMeshesOfGoveqns           => SOESetMeshesOfGoveqns
     procedure, public :: SetMeshesOfGoveqnsByMeshRank => SOESetMeshesOfGoveqnsByMeshRank
     procedure, public :: AddCouplingBCsInGovEqn       => SOEBaseAddCouplingBCsInGovEqn
     procedure, public :: AddConditionInGovEqn         => SOEBaseAddConditionInGovEqn
     procedure, public :: CreateVectorsForGovEqn       => SOBCreateVectorsForGovEqn
     procedure, public :: AllocateAuxVars              => SOEBaseAllocateAuxVars
     procedure, public :: ComputeNumInternalAuxVars    => SOEBaseComputeNumInternalAuxVars
  end type sysofeqns_base_type

  public :: SOEBaseInit
  public :: SOESetMeshesOfGoveqns

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine SOEBaseInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize a SoE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    this%name                         = ""
    this%itype                        = 0
    this%ngoveqns                     = 0

    this%mpi_rank                     = 0

    this%time                         = 0.d0
    this%dtime                        = 0.d0
    this%nstep                        = 0

    this%num_calls_to_ifunction = 0;
    this%num_calls_to_ijacobian = 0;

    this%num_auxvars_in         = 0
    this%num_auxvars_in_local   = 0
    this%num_auxvars_bc         = 0
    this%num_auxvars_ss         = 0

    call this%solver%Init()
    
    nullify(this%goveqns)
    nullify(this%soe_auxvars_bc_offset )
    nullify(this%soe_auxvars_ss_offset )
    nullify(this%soe_auxvars_bc_ncells )
    nullify(this%soe_auxvars_ss_ncells )


  end subroutine SOEBaseInit

  !------------------------------------------------------------------------
  subroutine SOEBaseResidual(this,snes, X, F, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc SNES.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    SNES                       :: snes
    PetscReal                  :: t
    Vec                        :: X
    Vec                        :: F
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseResidual must be extended')

  end subroutine SOEBaseResidual

  !------------------------------------------------------------------------
  subroutine SOEBaseJacobian(this, snes, X, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc SNES.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    SNES                       :: snes
    Vec                        :: X
    Mat                        :: A
    Mat                        :: B
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseJacobian must be extended')

  end subroutine SOEBaseJacobian

  !------------------------------------------------------------------------
  subroutine SOEComputeRHS(this, ksp, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc TS.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    KSP                        :: ksp
    Vec                        :: B
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SOEComputeRHS: '//&
         'SOEComputeRHS must be extended')

  end subroutine SOEComputeRHS

  !------------------------------------------------------------------------
  subroutine SOEComputeOperators(this, ksp, A, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine used when the SoE uses PETSc KSP.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    KSP                        :: ksp
    Mat                        :: A
    Mat                        :: B
    PetscErrorCode             :: ierr

    call endrun(msg='ERROR SOEComputeOperators: '//&
         'SOEComputeOperators must be extended')

  end subroutine SOEComputeOperators

  !------------------------------------------------------------------------
  subroutine SOEBasePreSolve(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that performs any required operations before calling
    ! PETSc solover.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBasePreSolve must be extended')

  end subroutine SOEBasePreSolve

  !------------------------------------------------------------------------
  subroutine SOEBaseSavePrimaryIndependentVar(this, X)
   !
   ! !DESCRIPTION:
   ! Save primary independent variable for each governing equation
   !
   implicit none
   !
   ! !ARGUMENTS
   class(sysofeqns_base_type) :: this
   Vec :: X
   !
   class(goveqn_base_type) , pointer :: cur_goveq
   PetscInt                          :: ii
   PetscInt                          :: nDM
   DM                      , pointer :: dms(:)
   Vec                     , pointer :: soln_subvecs(:)
   class(goveqn_base_type) , pointer :: cur_goveqn
   PetscErrorCode                    :: ierr

    ! Find number of GEs packed within the SoE
   call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr)

   ! Get DMs for each GE
   allocate (dms(nDM))
   call DMCompositeGetEntriesArray(this%solver%dm, dms, ierr)

   ! Allocate vectors for individual GEs
   allocate(soln_subvecs(nDM))

   ! Get solution vectors for individual GEs
   call DMCompositeGetAccessArray(this%solver%dm, &
        X, nDM, &
        PETSC_NULL_INTEGER, soln_subvecs, ierr)

   ii = 0
   cur_goveqn => this%goveqns
   do
      if (.not.associated(cur_goveqn)) exit
      ii = ii + 1
      call cur_goveqn%SavePrimaryIndependentVar(soln_subvecs(ii))
      cur_goveqn => cur_goveqn%next
   enddo

   call DMCompositeRestoreAccessArray(this%solver%dm, &
        X, nDM, &
        PETSC_NULL_INTEGER, soln_subvecs, ierr)

   ! Free memory
   deallocate(dms)
   deallocate(soln_subvecs)

      end subroutine SOEBaseSavePrimaryIndependentVar

 !------------------------------------------------------------------------
  subroutine SOEBasePreStepDT(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that performs any required operations before calling
    ! StepDT.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    call endrun(msg='ERROR SOEBasePreStepDT: '//&
         'SOEBasePreStepDT must be extended')

  end subroutine SOEBasePreStepDT

  !------------------------------------------------------------------------
  subroutine SOEBasePostStepDT(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that performs any required operations post StepDT.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    call endrun(msg='ERROR SOEBasePostStepDT: '//&
         'SOEBasePostStepDT must be extended')

  end subroutine SOEBasePostStepDT

  !------------------------------------------------------------------------
  subroutine SOEBaseStepDT(this, dt, nstep, converged, converged_reason, ierr)
    !
    ! !DESCRIPTION:
    ! Solves SoE by calling appropriate subroutine dependning on the choice
    ! of PETSc solver set for the SoE.
    !
    ! !USES
    use MultiPhysicsProbConstants, only : PETSC_TS
    use MultiPhysicsProbConstants, only : PETSC_SNES
    use MultiPhysicsProbConstants, only : PETSC_KSP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscReal                  :: dt
    PetscInt                   :: nstep
    PetscBool,intent(out)      :: converged
    PetscInt,intent(out)       :: converged_reason
    PetscErrorCode             :: ierr

    select case(this%solver%GetSolverType())
    case (PETSC_SNES)
       call SOEBaseStepDT_SNES(this, dt, nstep, converged, converged_reason, ierr)
    case (PETSC_KSP)
       call SOEBaseStepDT_KSP(this, dt, converged, ierr)
    case default
       write(iulog,*) 'VSFMMPPSetup: Unknown this%solver_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine SOEBaseStepDT

  !------------------------------------------------------------------------
  subroutine SOEBaseStepDT_SNES(soe, dt, nstep, converged, converged_reason, ierr)
    !
    ! !DESCRIPTION:
    ! Solves SoE via PETSc SNES
    !
    ! !USES
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: soe
    PetscErrorCode             :: ierr
    PetscInt                   :: nstep
    PetscBool,intent(out)      :: converged
    PetscInt,intent(out)       :: converged_reason
    PetscReal                  :: dt
    !
    ! !LOCAL VARIABLES:
    SNESConvergedReason        :: snes_reason
    PetscInt                   :: num_newton_iterations
    PetscInt                   :: num_linear_iterations
    PetscInt                   :: num_time_cuts
    PetscInt, parameter        :: max_num_time_cuts = 20
    PetscReal                  :: target_time
    PetscReal                  :: dt_iter
    PetscInt                   :: linesearch_iter
    SNESLineSearch             :: linesearch
    PetscInt                   :: max_linesearch_iter
    PetscInt, pointer          :: linesearch_iter_types(:)
    PetscInt, parameter        :: LS_BASIC = 1
    PetscInt, parameter        :: LS_BT    = 2
    PetscInt, parameter        :: LS_L2    = 3
    PetscInt, parameter        :: LS_CP    = 4
    PetscBool                  :: is_default_linesearch_basic
    PetscBool                  :: is_default_linesearch_bt
    PetscBool                  :: is_default_linesearch_l2
    PetscBool                  :: is_default_linesearch_cp
    character(len=32)          :: linesearch_name

    ! initialize
    linesearch_iter = 0
    num_time_cuts  = 0
    soe%time       = 0.d0
    target_time    = dt
    dt_iter        = dt
    soe%nstep      = nstep

    ! Determine the default linesearch option
    call SNESGetLineSearch(soe%solver%snes, linesearch, ierr); CHKERRQ(ierr)
    call PetscObjectTypeCompare(linesearch, SNESLINESEARCHBASIC, is_default_linesearch_basic, ierr); CHKERRQ(ierr)
    call PetscObjectTypeCompare(linesearch, SNESLINESEARCHBT, is_default_linesearch_bt, ierr); CHKERRQ(ierr)
    call PetscObjectTypeCompare(linesearch, SNESLINESEARCHL2, is_default_linesearch_l2, ierr); CHKERRQ(ierr)
    call PetscObjectTypeCompare(linesearch, SNESLINESEARCHCP, is_default_linesearch_cp, ierr); CHKERRQ(ierr)

    ! Check if the default linesearch option was a known option
    if ( (.not. is_default_linesearch_basic) .and. &
         (.not. is_default_linesearch_bt   ) .and. &
         (.not. is_default_linesearch_l2   ) .and. &
         (.not. is_default_linesearch_cp   ) ) then
       write(iulog,*) 'Unknown default linesearch option'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Create a list of linesearch options that will be tried before
    ! cutting the timestep
    if (is_default_linesearch_bt .or. is_default_linesearch_l2) then
       max_linesearch_iter = 2
       allocate(linesearch_iter_types(max_linesearch_iter))

       if (is_default_linesearch_bt) then
          linesearch_iter_types(1) = LS_BT
          linesearch_iter_types(2) = LS_L2
       else
          linesearch_iter_types(1) = LS_L2
          linesearch_iter_types(2) = LS_BT
       endif
    else
       max_linesearch_iter = 3
       allocate(linesearch_iter_types(max_linesearch_iter))

       if (is_default_linesearch_basic) linesearch_iter_types(1) = LS_BASIC
       if (is_default_linesearch_cp   ) linesearch_iter_types(1) = LS_CP
       linesearch_iter_types(2) = LS_L2
       linesearch_iter_types(3) = LS_BT
    endif

    do

       ! Set timestep
       call soe%SetDtime(dt_iter)

       ! Do any pre-solve operations
       call soe%PreSolve()

       select case (linesearch_iter_types(linesearch_iter+1))
       case (LS_BASIC)
          call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr); CHKERRQ(ierr)
          linesearch_name = 'linsearch_basic'
       case (LS_BT)
          call SNESLineSearchSetType(linesearch, SNESLINESEARCHBT, ierr); CHKERRQ(ierr)
          linesearch_name = 'linsearch_bt'
       case (LS_L2)
          call SNESLineSearchSetType(linesearch, SNESLINESEARCHL2, ierr); CHKERRQ(ierr)
          linesearch_name = 'linsearch_l2'
       case (LS_CP)
          call SNESLineSearchSetType(linesearch, SNESLINESEARCHCP, ierr); CHKERRQ(ierr)
          linesearch_name = 'linsearch_cp'
       end select

       ! Solve the nonlinear equation
       call SNESSolve(soe%solver%snes, PETSC_NULL_VEC, soe%solver%soln, ierr); CHKERRQ(ierr)

       ! Get reason why SNES iteration stopped
       call SNESGetConvergedReason(soe%solver%snes, snes_reason, ierr); CHKERRQ(ierr)

       converged_reason = snes_reason

       ! Did SNES converge?
       if (snes_reason < 0) then

          linesearch_iter = linesearch_iter + 1

          if (soe%solver%use_dynamic_linesearch .and. linesearch_iter < max_linesearch_iter) then
             ! Let's try another linesearch
             write(iulog,*),'On proc ', soe%mpi_rank, ' time_step = ', soe%nstep, &
                  linesearch_name // ' unsuccessful. Trying another one.'
             call VecCopy(soe%solver%soln_prev, soe%solver%soln, ierr); CHKERRQ(ierr)
          else
             ! All linesearch options failed, reset the linesearch iteration
             ! counter
             linesearch_iter = 0

             ! SNES diverged, so let's cut the timestep and try again.
             num_time_cuts = num_time_cuts + 1
             dt_iter = 0.5d0*dt_iter
             write(iulog,*),'On proc ', soe%mpi_rank, ' time_step = ', soe%nstep, &
                  'snes_reason = ',snes_reason,' cutting dt to ',dt_iter
          endif

          call VecCopy(soe%solver%soln_prev, soe%solver%soln, ierr); CHKERRQ(ierr)
       else
          ! SNES converged.
          converged = PETSC_TRUE
          soe%time = soe%time + dt_iter

          call SNESGetIterationNumber(soe%solver%snes,       &
               num_newton_iterations, ierr); CHKERRQ(ierr)
          call SNESGetLinearSolveIterations(soe%solver%snes, &
               num_linear_iterations, ierr); CHKERRQ(ierr)

          call soe%solver%IncrementNewtonIterCount(num_newton_iterations)
          call soe%solver%IncrementLinearIterCount(num_linear_iterations)
          !soe%cumulative_newton_iterations = soe%cumulative_newton_iterations + &
          !     num_newton_iterations
          !soe%cumulative_linear_iterations = soe%cumulative_linear_iterations + &
          !     num_linear_iterations

          ! Do any post-solve operations
          call soe%PostSolve()
       endif

       ! Do number of time cuts exceed maximum allowable number of cuts?
       if (num_time_cuts > max_num_time_cuts) then
          converged = PETSC_FALSE
          return
       endif

       if (soe%time >= target_time) exit
    enddo

    ! Set the linsearch type to be the default setting
    select case (linesearch_iter_types(1))
    case (LS_BASIC)
       call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr); CHKERRQ(ierr)
    case (LS_BT)
       call SNESLineSearchSetType(linesearch, SNESLINESEARCHBT, ierr); CHKERRQ(ierr)
    case (LS_L2)
       call SNESLineSearchSetType(linesearch, SNESLINESEARCHL2, ierr); CHKERRQ(ierr)
    case (LS_CP)
       call SNESLineSearchSetType(linesearch, SNESLINESEARCHCP, ierr); CHKERRQ(ierr)
    end select
    deallocate(linesearch_iter_types)


  end subroutine SOEBaseStepDT_SNES

  !------------------------------------------------------------------------
  subroutine SOEBaseStepDT_KSP(soe, dt, converged, ierr)
    !
    ! !DESCRIPTION:
    ! Solves SoE via PETSc KSP
    !
    ! !USES
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: soe
    PetscErrorCode             :: ierr
    PetscBool,intent(out)      :: converged
    PetscReal                  :: dt
    !
    ! !LOCAL VARIABLES:
    KSPConvergedReason         :: converged_reason
    PetscInt                   :: num_linear_iterations
    PetscReal                  :: target_time
    PetscReal                  :: dt_iter

    ! initialize
    soe%time       = 0.d0
    target_time    = dt
    dt_iter        = dt

    ! Set timestep
    call soe%SetDtime(dt_iter)

    ! Do any pre-solve operations
    call soe%PreSolve()

    ! Compute the 'b' vector
    call soe%ComputeRHS(soe%solver%ksp, soe%solver%rhs, ierr)

    ! Compute the 'A' matrix
    call soe%ComputeOperators(soe%solver%ksp, soe%solver%Amat, soe%solver%Amat, ierr);

    ! Set the A matrix
    call KSPSetOperators(soe%solver%ksp, soe%solver%Amat, soe%solver%Amat, ierr)

    ! Solve Ax = b
    call KSPSolve(soe%solver%ksp, soe%solver%rhs, soe%solver%soln, ierr)

    ! Did KSP converge?
    call KSPGetConvergedReason(soe%solver%ksp, converged_reason, ierr); CHKERRQ(ierr)

    ! Did KSP converge?
    if (converged_reason < 0) then
       write(iulog,*) 'KSP Diverged. Add new code'
       !call endrun(msg=errMsg(__FILE__, __LINE__))
       converged = PETSC_FALSE
    else
       converged = PETSC_TRUE
    endif

    ! Keep track of cumulative number of iterations
    call KSPGetIterationNumber(soe%solver%ksp, &
            num_linear_iterations, ierr); CHKERRQ(ierr)

    call soe%solver%IncrementLinearIterCount(num_linear_iterations)
    !soe%cumulative_linear_iterations = soe%cumulative_linear_iterations + &
    !     num_linear_iterations

    ! Do post-solve operations
    call soe%PostSolve()

  end subroutine SOEBaseStepDT_KSP

  !------------------------------------------------------------------------
  subroutine SOEBasePostSolve(this)
    !
    ! !DESCRIPTION:
    ! Subroutine that makes the copy of current solution after a successful
    ! PETSc solove.
    ! This subroutines may be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    !
    ! !LOCAL VARIABLES:
    PetscErrorCode             :: ierr

    call VecCopy(this%solver%soln, this%solver%soln_prev, ierr); CHKERRQ(ierr)

  end subroutine SOEBasePostSolve

  !------------------------------------------------------------------------
  subroutine SOEBasePrintInfo(this)
    !
    ! !DESCRIPTION:
    ! Displays information about SoE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer :: cur_goveqn

    write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    write(iulog,*)'  SystemOfEqns_name  : ',trim(this%name)
    write(iulog,*)'  SystemOfEqns_itype : ',this%itype
    write(iulog,*)''

    cur_goveqn => this%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       call cur_goveqn%PrintInfo()
       cur_goveqn => cur_goveqn%next
    enddo
    write(iulog,*)'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  end subroutine SOEBasePrintInfo

  !------------------------------------------------------------------------
  subroutine SOESetMeshesOfGoveqns(soe, meshes, nmesh)
    !
    ! !DESCRIPTION:
    ! Match the meshes in `meshes` with the governing equations in `soe`.
    !
    ! !USES
    use MeshType                     , only : mesh_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) , intent(inout)       :: soe
    class(mesh_type)           , pointer, intent(in) :: meshes(:)
    PetscInt                   , intent(in)          :: nmesh
    !
    ! !LOCAL VARIABLES:
    PetscInt                                         :: imesh
    PetscInt                                         :: mesh_itype
    PetscBool                                        :: mesh_found
    class(mesh_type)           , pointer             :: cur_mesh
    class(goveqn_base_type)    , pointer             :: cur_goveqn

    cur_goveqn => soe%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       mesh_itype = cur_goveqn%mesh_itype

       mesh_found = PETSC_FALSE
       do imesh = 1, nmesh
          cur_mesh => meshes(imesh)
          if (mesh_itype == cur_mesh%itype) then
             cur_goveqn%mesh => cur_mesh
             mesh_found = PETSC_TRUE
             exit
          endif
       enddo

       if (.not.mesh_found) then
          call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
               'Mesh not found for Governing equation within the list')
       endif

       cur_goveqn => cur_goveqn%next
    enddo

  end subroutine SOESetMeshesOfGoveqns

  !------------------------------------------------------------------------
  subroutine SOESetMeshesOfGoveqnsByMeshRank(soe, meshes, nmesh)
    !
    ! !DESCRIPTION:
    ! Match the meshes in `meshes` with the governing equations in `soe`.
    !
    ! !USES
    use MeshType                     , only : mesh_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) , intent(inout)       :: soe
    class(mesh_type)           , pointer, intent(in) :: meshes(:)
    PetscInt                   , intent(in)          :: nmesh
    !
    ! !LOCAL VARIABLES:
    PetscInt                                         :: imesh
    PetscInt                                         :: mesh_rank
    class(mesh_type)           , pointer             :: cur_mesh
    class(goveqn_base_type)    , pointer             :: cur_goveqn

    cur_goveqn => soe%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       mesh_rank = cur_goveqn%mesh_rank

       if (mesh_rank <= 0) then
          call endrun(msg='ERROR SystemOfEquationsBaseType: '// &
               'Invalid mesh rank associated with governing equation: ' // &
               trim(cur_goveqn%name))
       endif

       if (mesh_rank > nmesh) then
          call endrun(msg='ERROR SystemOfEquationsBaseType: '// &
               'Rank of mesh associated with governing equation ' // &
               'exceeds no. of meshes in the list')
       endif

       cur_mesh => meshes(mesh_rank)
       cur_goveqn%mesh => cur_mesh
       cur_goveqn%mesh_itype = cur_mesh%itype

       cur_goveqn => cur_goveqn%next
    enddo

  end subroutine SOESetMeshesOfGoveqnsByMeshRank

  !------------------------------------------------------------------------
  subroutine SOEBaseSetPointerToIthGovEqn(this, igoveqn, goveqn_ptr)
    !
    ! !DESCRIPTION:
    ! Returns pointer to the i-th governing equation present with SoE
    !
    ! !USES
    use MeshType , only : mesh_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type)      :: this
    PetscInt                        :: igoveqn
    class(goveqn_base_type),pointer :: goveqn_ptr
    !
    ! !LOCAL VARIABLES
    class(goveqn_base_type),pointer :: cur_goveq
    PetscInt                        :: sum_goveqn
    PetscBool                       :: found

    if (igoveqn > this%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    sum_goveqn = 0
    found      = PETSC_FALSE

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       sum_goveqn = sum_goveqn + 1
       if (sum_goveqn == igoveqn) then
          found = PETSC_TRUE
          goveqn_ptr => cur_goveq
          exit
       endif

       cur_goveq => cur_goveq%next
    enddo

    if (.not.found) then
       write(iulog,*) 'In SOEBaseSetPointerToIthGovEqn: ith-goveqn not found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

  end subroutine SOEBaseSetPointerToIthGovEqn

  !------------------------------------------------------------------------
  subroutine SOEBaseSetDtime(this, dtime)
    !
    ! !DESCRIPTION:
    ! Sets timestep for SoE and all governing equations present within SoE.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscReal                  :: dtime
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type),pointer   :: cur_goveqn

    this%dtime = dtime

    cur_goveqn => this%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       call cur_goveqn%SetDtime(dtime)
       cur_goveqn => cur_goveqn%next
    enddo

  end subroutine SOEBaseSetDtime

  !------------------------------------------------------------------------
  subroutine SOEBaseSetDataFromCLM(this, soe_auxvar_type, var_type, soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that sets data from CLM before proceeding to solve the SoE.
    ! E.g. Getting infilitration source from CLM for VSFM.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscInt, intent(in)       :: var_type
    PetscInt                   :: soe_auxvar_type
    PetscInt                   :: soe_auxvar_id
    PetscReal                  :: data_1d(:)
    !
    ! !LOCAL VARIABLES:

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseSetDataFromCLM must be extended')

  end subroutine SOEBaseSetDataFromCLM

  !------------------------------------------------------------------------
  subroutine SOEBaseGetDataForCLM(this, soe_auxvar_type, var_type, soe_auxvar_id, data_1d)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that extracts data from SoE's data structure for CLM after
    ! a successfull SoE solve. E.g. Getting updated soil moisture value from VSFM
    ! for CLM.
    ! This subroutines needs to be extended by a child class.
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscInt, intent(in)       :: var_type
    PetscInt                   :: soe_auxvar_type
    PetscInt                   :: soe_auxvar_id
    PetscReal                  :: data_1d(:)
    !
    ! !LOCAL VARIABLES:

    call endrun(msg='ERROR SystemOfEquationsBaseType: '//&
         'SOEBaseGetDataForCLM must be extended')

  end subroutine SOEBaseGetDataForCLM

  !------------------------------------------------------------------------
  subroutine SOEBaseAddGovEqn(this, geq_type, name, mesh_itype)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation
    !
   use MultiPhysicsProbConstants     , only : GE_THERM_SOIL_EBASED
   use MultiPhysicsProbConstants     , only : GE_RE
   use GoverningEquationBaseType     , only : goveqn_base_type
   use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
   use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
   !
   implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscInt                   :: geq_type
    character(len =*)          :: name
    PetscInt                   :: mesh_itype
   !
    class (goveqn_base_type)                  , pointer :: cur_goveqn
    class (goveqn_thermal_enthalpy_soil_type) , pointer :: goveq_soil
    class (goveqn_richards_ode_pressure_type) , pointer :: goveq_richards
    integer                                             :: igoveqn

    cur_goveqn => this%goveqns

    do igoveqn = 1, this%ngoveqns - 1
       cur_goveqn => cur_goveqn%next
    enddo

    this%ngoveqns = this%ngoveqns + 1

    select case(geq_type)
    case (GE_THERM_SOIL_EBASED)

       allocate(goveq_soil)
       call goveq_soil%Setup()

       goveq_soil%name                 = trim(name)
       goveq_soil%rank_in_soe_list     = this%ngoveqns
       goveq_soil%mesh_itype           = mesh_itype

       if (this%ngoveqns               == 1) then
          this%goveqns                 => goveq_soil
       else
          cur_goveqn%next              => goveq_soil
       endif

    case (GE_RE)

       allocate(goveq_richards)
       call goveq_richards%Setup()

       goveq_richards%name             = trim(name)
       goveq_richards%rank_in_soe_list = this%ngoveqns
       goveq_richards%mesh_itype       = mesh_itype

       if (this%ngoveqns == 1) then
          this%goveqns => goveq_richards
       else
          cur_goveqn%next => goveq_richards
       endif

    case default
       write(iulog,*) 'Unknown governing equation type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
    
    !write(iulog,*) 'SOEBaseAddGovEqn must be extended'
    !call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine SOEBaseAddGovEqn

  !------------------------------------------------------------------------
  subroutine SOEBaseAddGovEqnWithMeshRank(this, geq_type, name, mesh_rank)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation
    !
   use GoveqnCanopyAirTemperatureType  , only : goveqn_cair_temp_type
   use GoveqnCanopyAirVaporType        , only : goveqn_cair_vapor_type
   use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
   use GoveqnRichardsODEPressureType   , only : goveqn_richards_ode_pressure_type
   use GoveqnLeafBoundaryLayer         , only : goveqn_leaf_bnd_lyr_type
   use GoveqnPhotosynthesisType        , only : goveqn_photosynthesis_type
   use GoveqnLongwaveType              , only : goveqn_longwave_type
   use MultiPhysicsProbConstants       , only : GE_CANOPY_AIR_TEMP
   use MultiPhysicsProbConstants       , only : GE_CANOPY_AIR_VAPOR
   use MultiPhysicsProbConstants       , only : GE_CANOPY_LEAF_TEMP
   use MultiPhysicsProbConstants       , only : GE_RE
   use MultiPhysicsProbConstants       , only : GE_LEAF_BND_LAYER
   use MultiPhysicsProbConstants       , only : GE_PHOTOSYNTHESIS
   use MultiPhysicsProbConstants       , only : GE_LONGWAVE
   !
   implicit none
   !
   ! !ARGUMENTS
    class(sysofeqns_base_type) :: this
    PetscInt                   :: geq_type
    character(len =*)          :: name
    PetscInt                   :: mesh_rank
   !
    class (goveqn_base_type)                  , pointer :: cur_goveqn
    class (goveqn_cair_temp_type)             , pointer :: geq_air_temp
    class (goveqn_cair_vapor_type)            , pointer :: geq_air_vapor
    class (goveqn_cleaf_temp_type)            , pointer :: geq_leaf_temp
    class (goveqn_richards_ode_pressure_type) , pointer :: goveq_richards
    class (goveqn_leaf_bnd_lyr_type)          , pointer :: goveq_lbl
    class (goveqn_photosynthesis_type)        , pointer :: goveq_phtsyn
    class (goveqn_longwave_type)        , pointer :: goveq_longwave
    integer                                             :: igoveqn

    cur_goveqn => this%goveqns

    do igoveqn = 1, this%ngoveqns - 1
       cur_goveqn => cur_goveqn%next
    enddo

    this%ngoveqns = this%ngoveqns + 1

    select case(geq_type)
    case (GE_CANOPY_AIR_TEMP)
       allocate(geq_air_temp)

       call geq_air_temp%Setup()

       geq_air_temp%name              = trim(name)
       geq_air_temp%rank_in_soe_list  = this%ngoveqns
       geq_air_temp%mesh_rank         = mesh_rank

       if (this%ngoveqns == 1) then
          this%goveqns => geq_air_temp
       else
          cur_goveqn%next => geq_air_temp
       endif

    case (GE_CANOPY_AIR_VAPOR)
       allocate(geq_air_vapor)

       call geq_air_vapor%Setup()

       geq_air_vapor%name              = trim(name)
       geq_air_vapor%rank_in_soe_list  = this%ngoveqns
       geq_air_vapor%mesh_rank         = mesh_rank

       if (this%ngoveqns == 1) then
          this%goveqns => geq_air_vapor
       else
          cur_goveqn%next => geq_air_vapor
       endif

    case (GE_CANOPY_LEAF_TEMP)
       allocate(geq_leaf_temp)

       call geq_leaf_temp%Setup()

       geq_leaf_temp%name              = trim(name)
       geq_leaf_temp%rank_in_soe_list  = this%ngoveqns
       geq_leaf_temp%mesh_rank         = mesh_rank

       if (this%ngoveqns == 1) then
          this%goveqns => geq_leaf_temp
       else
          cur_goveqn%next => geq_leaf_temp
       endif

      case (GE_RE)

         allocate(goveq_richards)
         call goveq_richards%Setup()
  
         goveq_richards%name              = trim(name)
         goveq_richards%rank_in_soe_list  = this%ngoveqns
         goveq_richards%mesh_rank         = mesh_rank
  
         if (this%ngoveqns == 1) then
            this%goveqns => goveq_richards
         else
            cur_goveqn%next => goveq_richards
         endif
      case (GE_LEAF_BND_LAYER)
         allocate(goveq_lbl)
         call goveq_lbl%Setup()
         goveq_lbl%name              = trim(name)
         goveq_lbl%rank_in_soe_list  = this%ngoveqns
         goveq_lbl%mesh_rank         = mesh_rank
  
         if (this%ngoveqns == 1) then
            this%goveqns => goveq_lbl
         else
            cur_goveqn%next => goveq_lbl
         endif

      case (GE_PHOTOSYNTHESIS)
         allocate(goveq_phtsyn)
         call goveq_phtsyn%Setup()
         goveq_phtsyn%name              = trim(name)
         goveq_phtsyn%rank_in_soe_list  = this%ngoveqns
         goveq_phtsyn%mesh_rank         = mesh_rank
  
         if (this%ngoveqns == 1) then
            this%goveqns => goveq_phtsyn
         else
            cur_goveqn%next => goveq_phtsyn
         endif

      case (GE_LONGWAVE)
         allocate(goveq_longwave)
         call goveq_longwave%Setup()
         goveq_longwave%name              = trim(name)
         goveq_longwave%rank_in_soe_list  = this%ngoveqns
         goveq_longwave%mesh_rank         = mesh_rank
  
         if (this%ngoveqns == 1) then
            this%goveqns => goveq_longwave
         else
            cur_goveqn%next => goveq_longwave
         endif

      case default
       write(iulog,*) 'Unknown governing equation type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
    
    !write(iulog,*) 'SOEBaseAddGovEqnWithMeshRank must be extended'
    !call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine SOEBaseAddGovEqnWithMeshRank

  !------------------------------------------------------------------------
  subroutine SOEBaseAddCouplingBCsInGovEqn(this, igoveq, name, unit, &
       num_other_goveqs, id_of_other_goveqs, &
       icoupling_of_other_goveqns, region_type, conn_set)
    !
    ! !DESCRIPTION:
    ! Adds coupling boundary conditions to igoveq governing equation
    !
    use MultiPhysicsProbConstants, only : COND_BC
    use GoverningEquationBaseType, only : goveqn_base_type
    use ConnectionSetType        , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type)                      :: this
    PetscInt                  , intent(in)          :: igoveq
    character(len =*)         , intent(in)          :: name
    character(len =*)         , intent(in)          :: unit
    PetscInt                  , intent(in)          :: num_other_goveqs
    PetscInt                  , pointer, intent(in) :: id_of_other_goveqs(:)
    PetscBool                 , pointer, optional   :: icoupling_of_other_goveqns(:)
    PetscInt                  , intent(in), optional:: region_type
    type(connection_set_type) , pointer, optional   :: conn_set
    !
    ! !LOCAL VARIABLES
    class(goveqn_base_type)   , pointer             :: cur_goveq
    class(goveqn_base_type)   , pointer             :: other_goveq
    PetscInt                                        :: cond_type
    PetscInt                                        :: ss_or_bc_type
    PetscInt                                        :: jj
    PetscInt                  , pointer             :: itype_of_other_goveqs(:)

    if (igoveq > this%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    call this%SetPointerToIthGovEqn(igoveq, cur_goveq)

    allocate(itype_of_other_goveqs(num_other_goveqs))
    do jj = 1, num_other_goveqs
       call this%SetPointerToIthGovEqn(id_of_other_goveqs(jj), other_goveq)
       itype_of_other_goveqs(jj) = other_goveq%itype
    enddo

    call cur_goveq%AddCouplingBC(                   &
         name, unit, num_other_goveqs,              &
         id_of_other_goveqs, itype_of_other_goveqs, &
         icoupling_of_other_goveqns, region_type,   &
         conn_set )

    deallocate(itype_of_other_goveqs)

  end subroutine SOEBaseAddCouplingBCsInGovEqn

  !------------------------------------------------------------------------
  subroutine SOEBaseAddConditionInGovEqn(this, igoveqn, ss_or_bc_type, name, unit, &
       cond_type, region_type, conn_set)
    !
    ! !DESCRIPTION:
    ! Adds a boundary/source-sink condition to a governing equation
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    use ConnectionSetType        , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type)                  :: this
    PetscInt                                    :: igoveqn
    PetscInt                                    :: ss_or_bc_type
    character(len =*)                           :: name
    character(len =*)                           :: unit
    PetscInt                                    :: cond_type
    PetscInt, optional                          :: region_type
    type(connection_set_type),pointer, optional :: conn_set
    !
    class(goveqn_base_type),pointer             :: cur_goveq
    class(goveqn_base_type),pointer             :: other_goveq
    PetscInt                                    :: ii

    call this%SetPointerToIthGovEqn(igoveqn, cur_goveq)

    call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
         cond_type, region_type, conn_set)

  end subroutine SOEBaseAddConditionInGovEqn

  !------------------------------------------------------------------------
  subroutine SOEBaseClean(this)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

    !call this%goveqns%Clean

    nullify(this%goveqns)

  end subroutine SOEBaseClean

  !------------------------------------------------------------------------
  subroutine SOBCreateVectorsForGovEqn(this)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine that is extended by child SoE class
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_base_type) :: this

  end subroutine SOBCreateVectorsForGovEqn

  !------------------------------------------------------------------------
  subroutine SOEBaseAllocateAuxVars(this)
   !
   ! !DESCRIPTION:
   ! Dummy subroutine that is extended by child SoE class
   !
   implicit none
   !
   ! !ARGUMENTS
   class(sysofeqns_base_type) :: this

 end subroutine SOEBaseAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine SOEBaseComputeNumInternalAuxVars(this)
   !
   ! !DESCRIPTION:
   ! Computes the total number of grid cells in all governing equations
   !
   implicit none
   !
   ! !ARGUMENTS
   class(sysofeqns_base_type) :: this
   !
   class(goveqn_base_type)    , pointer :: cur_goveq

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       this%num_auxvars_in = this%num_auxvars_in + &
            cur_goveq%mesh%ncells_all

       cur_goveq => cur_goveq%next
    enddo

 end subroutine SOEBaseComputeNumInternalAuxVars

#endif

end module SystemOfEquationsBaseType
