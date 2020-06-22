
module SystemOfEquationsMlcType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                      , only : iulog
  use mpp_abortutils                  , only : endrun
  use mpp_shr_log_mod                 , only : errMsg => shr_log_errMsg
  use SystemOfEquationsBaseType       , only : sysofeqns_base_type
  use SystemOfEquationsMlcAuxType     , only : sysofeqns_mlc_auxvar_type
  use GoveqnCanopyAirTemperatureType  , only : goveqn_cair_temp_type
  use GoveqnCanopyAirVaporType        , only : goveqn_cair_vapor_type
  use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
  use CanopyTurbulenceAuxType         , only : canopy_turbulence_auxvar_type
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
  private

  type, public, extends(sysofeqns_base_type) :: sysofeqns_mlc_type

     type (sysofeqns_mlc_auxvar_type), pointer :: aux_vars_in(:) ! Internal state.
     type (sysofeqns_mlc_auxvar_type), pointer :: aux_vars_bc(:) ! Boundary conditions.
     type (sysofeqns_mlc_auxvar_type), pointer :: aux_vars_ss(:) ! Source-sink.

     type (canopy_turbulence_auxvar_type) :: cturb
     PetscInt :: ncair
     PetscInt, pointer :: leaftemp_goveqn_rank(:)

   contains

    procedure, public :: Init                  => MlcSoeInit
    procedure, public :: PreSolve              => MlcSoePreSolve
    procedure, public :: PostSolve             => MlcSoePostSolve
    procedure, public :: ComputeRHS            => MlcSoeComputeRhs
    procedure, public :: ComputeOperators      => MlcSoeComputeOperators

 end type sysofeqns_mlc_type

contains

  !------------------------------------------------------------------------
  subroutine MlcSoeInit (this)
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
    class(sysofeqns_mlc_type) :: this

    call SOEBaseInit(this)

    nullify(this%aux_vars_in)
    nullify(this%aux_vars_bc)
    nullify(this%aux_vars_ss)

    this%ncair = 0

  end subroutine MlcSoeInit

  !------------------------------------------------------------------------
  subroutine MlcSoePreSolve(this)
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
    class(sysofeqns_mlc_type)         :: this
    class(goveqn_base_type) , pointer :: cur_goveq
    PetscInt                          :: offset, ii
    PetscInt                          :: nDM
    DM                      , pointer :: dms(:)
    Vec                     , pointer :: soln_subvecs(:)
    PetscReal               , pointer :: v_p(:)
    class(goveqn_base_type) , pointer :: cur_goveqn
    PetscErrorCode                    :: ierr

    call ObukhovLength(this%cturb)
    call WindProfile(this%cturb)
    call AerodynamicConductances(this%cturb)

    cur_goveq => this%goveqns
    do
       if (.not.associated(cur_goveq)) exit
       select type(cur_goveq)

       class is (goveqn_cair_temp_type)
          call cur_goveq%GetFromSOEAuxVarsCturb(this%cturb)

       class is (goveqn_cair_vapor_type)
          call cur_goveq%GetFromSOEAuxVarsCturb(this%cturb)

       class is (goveqn_cleaf_temp_type)
          call cur_goveq%GetFromSOEAuxVarsCturb(this%cturb)
       end select

       cur_goveq => cur_goveq%next
    enddo

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(this%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(this%solver%dm, &
         this%solver%soln, nDM, &
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
         this%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

  end subroutine MlcSoePreSolve

  !------------------------------------------------------------------------
  subroutine MlcSoeComputeRhs(this, ksp, B, ierr)
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
    class(sysofeqns_mlc_type) :: this
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
    do row = 1,nDM
       do col = row+1,nDM
          call this%SetPointerToIthGovEqn(row, cur_goveq_1)
          call this%SetPointerToIthGovEqn(col, cur_goveq_2)
          call GovEqnExchangeAuxVars(this, cur_goveq_1, cur_goveq_2)
          call GovEqnExchangeAuxVars(this, cur_goveq_2, cur_goveq_1)
       enddo
    enddo

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

  end subroutine MlcSoeComputeRhs

  !------------------------------------------------------------------------

  subroutine GovEqnExchangeAuxVars(soe, cur_goveq_1, cur_goveq_2)
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType , only : goveqn_base_type
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use CouplingVariableType      , only : coupling_variable_type
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_mlc_type)         :: soe
    class(goveqn_base_type) , pointer :: cur_goveq_1
    class(goveqn_base_type) , pointer :: cur_goveq_2
    !
    ! !LOCAL VARIABLES:
    type(connection_set_type)     , pointer :: cur_conn_set_2
    type(condition_type)          , pointer :: cur_cond_2
    type (coupling_variable_type) , pointer :: cpl_var_1
    PetscInt                                :: idx
    PetscInt                                :: iauxvar
    PetscInt                                :: var_type
    PetscInt                                :: ghosted_id
    PetscInt                                :: bc_idx
    PetscInt                                :: bc_offset
    PetscInt                                :: bc_rank_in_cpl_eqn
    PetscReal                               :: var_value
    PetscBool                               :: bc_found
    PetscBool                               :: is_bc
    PetscInt                                :: geq_leaf_temp_rank
    PetscReal                     , pointer :: var_values(:)
    PetscInt                                :: nauxvar

    cpl_var_1 => cur_goveq_1%coupling_vars%first
    do
       if (.not.associated(cpl_var_1)) exit

       if (cpl_var_1%rank_of_coupling_goveqn == &
           cur_goveq_2%rank_in_soe_list) then

          var_type           = cpl_var_1%variable_type
          is_bc              = cpl_var_1%variable_is_bc_in_coupling_goveqn
          bc_offset          = cpl_var_1%offset_of_bc_in_current_goveqn
          bc_rank_in_cpl_eqn = cpl_var_1%rank_of_bc_in_coupling_goveqn

          if (.not.is_bc) then

             geq_leaf_temp_rank = 0
             

             select type(cur_goveq_2)
             class is (goveqn_cair_temp_type)
                nauxvar = size(cur_goveq_2%aux_vars_in)
                allocate(var_values(nauxvar))

                call cur_goveq_2%GetRValues(AUXVAR_INTERNAL, var_type, nauxvar, var_values)

             class is (goveqn_cair_vapor_type)
                nauxvar = size(cur_goveq_2%aux_vars_in)
                allocate(var_values(nauxvar))

                call cur_goveq_2%GetRValues(AUXVAR_INTERNAL, var_type, nauxvar, var_values)

             class is (goveqn_cleaf_temp_type)
                nauxvar = size(cur_goveq_2%aux_vars_in)
                allocate(var_values(nauxvar))
                call cur_goveq_2%GetRValues(AUXVAR_INTERNAL, var_type, nauxvar, var_values)
                geq_leaf_temp_rank = soe%leaftemp_goveqn_rank(cur_goveq_2%rank_in_soe_list)

             class default
                write(iulog,*)'EGovEqnExchangeAuxVars: Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select

             select type(cur_goveq_1)
             class is (goveqn_cair_temp_type)
                call cur_goveq_1%SetRValues(AUXVAR_INTERNAL, var_type, geq_leaf_temp_rank, nauxvar, var_values)

             class is (goveqn_cair_vapor_type)
                call cur_goveq_1%SetRValues(AUXVAR_INTERNAL, var_type, geq_leaf_temp_rank, nauxvar, var_values)

             class is (goveqn_cleaf_temp_type)
                call cur_goveq_1%SetRValues(AUXVAR_INTERNAL, var_type, nauxvar, var_values)

             class default
                write(iulog,*)'GovEqnExchangeAuxVars: Unknown class'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select


           else

             write(iulog,*) 'GovEqnExchangeAuxVars: Extend code to ' // &
                  'exchange boundary condition data'
             call endrun(msg=errMsg(__FILE__, __LINE__))

          endif ! if (.not.is_bc)

        endif

        cpl_var_1 => cpl_var_1%next
        
     enddo

   end subroutine GovEqnExchangeAuxVars

  !------------------------------------------------------------------------
  subroutine MlcSoeComputeOperators(this, ksp, A, B, ierr)
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
    class(sysofeqns_mlc_type) :: this
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

       cur_goveq_2 => cur_goveq_1%next
       col = row
       do
          if (.not.associated(cur_goveq_2)) exit

          col = col + 1

          rank_1 = cur_goveq_1%rank_in_soe_list
          rank_2 = cur_goveq_2%rank_in_soe_list

          select type(cur_goveq_1)
          class is (goveqn_cleaf_temp_type)
             rank_1 = this%leaftemp_goveqn_rank(cur_goveq_1%rank_in_soe_list)
          end select

          select type(cur_goveq_2)
          class is (goveqn_cleaf_temp_type)
             rank_2 = this%leaftemp_goveqn_rank(cur_goveq_2%rank_in_soe_list)
          end select

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

  end subroutine MlcSoeComputeOperators

  !------------------------------------------------------------------------
  subroutine MlcSoePostSolve(this)
    !
    ! !DESCRIPTION:
    ! Peform operations after a successful call to the PETSc solver.
    !
    ! !USES:
    use GoverningEquationBaseType     , only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_mlc_type) :: this
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type) , pointer :: cur_goveqn
    PetscInt                          :: ii, nDM
    Vec                     , pointer :: soln_subvecs(:)
    PetscErrorCode                    :: ierr

    call VecCopy(this%solver%soln, this%solver%soln_prev,ierr); CHKERRQ(ierr)

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(this%solver%dm, nDM, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(this%solver%dm, &
         this%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    ii = 0
    cur_goveqn => this%goveqns
    do
       if (.not.associated(cur_goveqn)) exit
       ii = ii + 1
       call cur_goveqn%SavePrimaryIndependentVar(soln_subvecs(ii))
       cur_goveqn => cur_goveqn%next
    enddo

    deallocate(soln_subvecs)

  end subroutine MlcSoePostSolve

#endif

end module SystemOfEquationsMlcType
