module GoveqnLongwaveType

#ifdef USE_PETSC_LIB
  !#define USE_BONAN_FORMULATION

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType , only : goveqn_base_type
  use LongwaveAuxType           , only : longwave_auxvar_type
  use petscvec
  use petscmat
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_longwave_type

     type(longwave_auxvar_type)      , pointer :: aux_vars_in(:)
     type(longwave_auxvar_type)      , pointer :: aux_vars_bc(:)

     ! For post-processing data to compute absorbed radiation as:
     ! Iabs = Mat * Iupdn + RHS
     Vec :: flux_up_dn
     Vec :: flux_absorbed_rhs
     Vec :: flux_absorbed
     Mat :: flux_absorbed_mat

     PetscInt                                   :: nLeaf, nleafGE

   contains

     procedure, public :: Setup                     => LongwaveSetup
     procedure, public :: AllocateAuxVars           => LongwaveAllocateAuxVars
     procedure, public :: PreSolve                  => LongwavePreSolve
     procedure, public :: PostSolve                 => LongwavePostSolve
     procedure, public :: UpdateAuxVars             => LongwaveUpdateAuxVars
     procedure, public :: SavePrimaryIndependentVar => LongwaveSavePrmIndepVar
     procedure, public :: ComputeRhs                => LongwaveComputeRhs
     procedure, public :: ComputeOperatorsDiag      => LongwaveComputeOperatorsDiag

  end type goveqn_longwave_type

contains
  !------------------------------------------------------------------------
  subroutine LongwaveSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_LONGWAVE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this

    call this%Create()

    this%name  = ""
    this%itype = GE_LONGWAVE
    this%dof   = 2

    nullify(this%aux_vars_in)
    nullify(this%aux_vars_bc)

    this%nLeaf = 1

  end subroutine LongwaveSetup

  !---------------------------------------------------------------------
  subroutine LongwaveAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates memory for storing auxiliary variables associated with:
    !   + Internal control volumes,
    !
    ! !USES:
    !
    use ConditionType             , only : condition_type
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type)    :: this
    !
    type(condition_type) , pointer :: cur_cond
    PetscInt                       :: ghosted_id, ncells_cond, icond
    PetscInt                       :: nrow, ncol
    PetscErrorCode                 :: ierr

    ! Allocate memory and initialize aux vars: For internal cells
    allocate(this%aux_vars_in(this%mesh%ncells_all))

    ! Update aux vars for internal cells
    do ghosted_id = 1, this%mesh%ncells_all
       if (this%mesh%is_active(ghosted_id)) then
          call this%aux_vars_in(ghosted_id)%Init(this%nLeaf)
       end if
    enddo

    ! Allocate memory and initialize aux vars: For boundary cells
    ncells_cond = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit
       ncells_cond = ncells_cond + cur_cond%ncells
       cur_cond => cur_cond%next
    enddo

    allocate(this%aux_vars_bc(ncells_cond))
    do icond = 1,ncells_cond
       call this%aux_vars_bc(icond)%Init(0)
    enddo

    nrow = this%mesh%ncells_all
    ncol = 2*nrow
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrow, ncol, &
         3, PETSC_NULL_INTEGER, this%flux_absorbed_mat, ierr); CHKERRQ(ierr)

    call VecCreateSeq(PETSC_COMM_SELF, nrow, this%flux_absorbed, ierr); CHKERRQ(ierr)
    call VecDuplicate(this%flux_absorbed, this%flux_absorbed_rhs, ierr); CHKERRQ(ierr);

    call VecCreateSeq(PETSC_COMM_SELF, ncol, this%flux_up_dn, ierr); CHKERRQ(ierr)

  end subroutine LongwaveAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine LongwavePreSolve(this)
    !
    ! !DESCRIPTION:
    ! Presolve
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this

    call this%UpdateAuxVars()
    call LongwaveSetupAbsorbedFluxMatAndVec(this)

  end subroutine LongwavePreSolve

  !------------------------------------------------------------------------
  subroutine LongwavePostSolve(this)
    !
    ! !DESCRIPTION:
    ! Presolve
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this
    !
    PetscInt             :: ghosted_id
    PetscScalar, pointer :: f_p(:)
    PetscErrorCode       :: ierr

    call VecGetArrayF90(this%flux_up_dn, f_p, ierr); CHKERRQ(ierr)

    do ghosted_id = 1, this%mesh%ncells_local
       f_p((ghosted_id-1)*2 + 1) = this%aux_vars_in(ghosted_id)%Iup
       f_p((ghosted_id-1)*2 + 2) = this%aux_vars_in(ghosted_id)%Idn
    end do

    call VecRestoreArrayF90(this%flux_up_dn, f_p, ierr); CHKERRQ(ierr)

    call MatMultAdd(this%flux_absorbed_mat, this%flux_up_dn, this%flux_absorbed_rhs, this%flux_absorbed, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(this%flux_absorbed, f_p, ierr); CHKERRQ(ierr)
    do ghosted_id = 1, this%mesh%ncells_local
       this%aux_vars_in(ghosted_id)%Iabs = f_p(ghosted_id)
    end do
    call VecRestoreArrayF90(this%flux_absorbed, f_p, ierr); CHKERRQ(ierr)

  end subroutine LongwavePostSolve

 !------------------------------------------------------------------------
  subroutine LongwaveSetupAbsorbedFluxMatAndVec(this)
    !
    ! !DESCRIPTION:
    !
    use ConditionType             , only : condition_type
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : STEFAN_BOLTZMAN_CONSTANT
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this
    Vec                                   :: B
    PetscErrorCode                        :: ierr
    !
    type(longwave_auxvar_type) , pointer  :: avars(:)
    class(connection_set_type) , pointer  :: cur_conn_set
    class(condition_type)      , pointer  :: cur_cond
    PetscScalar                , pointer  :: r_p(:)
    PetscInt                              :: icell, row, col
    PetscReal                             :: value

    avars => this%aux_vars_in

    call VecGetArrayF90(this%flux_absorbed_rhs, r_p, ierr)

    ! Set RHS for upwind and absorbed radiation
    do icell = 1, this%mesh%ncells_local

       if (avars(icell)%is_soil) then
          r_p((icell-1)*this%dof) = 0

          value = -1.d0; row = icell-1; col = 2*row;
          call MatSetValues(this%flux_absorbed_mat, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr);

          value = 1.d0; row = icell-1; col = 2*row + 1
          call MatSetValues(this%flux_absorbed_mat, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr);
       else
          ! Absorbed flux: -2 * (1-tau_{d,i}) * emiss * sigma * T^4
          r_p(icell) = -2.d0 * avars(icell)%rad_source

          value = avars(icell)%leaf_emiss * (1.d0 - avars(icell)%trans)
          row = icell-1; col = 2*row-2;
          call MatSetValues(this%flux_absorbed_mat, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr);

          row = icell-1; col = 2*row+1;
          call MatSetValues(this%flux_absorbed_mat, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr);
       end if

    end do

    call VecRestoreArrayF90(this%flux_absorbed_rhs, r_p, ierr)

    call MatAssemblyBegin(this%flux_absorbed_mat, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  this%flux_absorbed_mat, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  end subroutine LongwaveSetupAbsorbedFluxMatAndVec

  !------------------------------------------------------------------------
  subroutine LongwaveSavePrmIndepVar (this, x)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this
    Vec :: x
    !
    PetscScalar, pointer :: x_p(:)
    PetscInt             :: ghosted_id, size
    PetscErrorCode       :: ierr

    call VecGetLocalSize(x, size, ierr); CHKERRQ(ierr)

    if (size /= this%mesh%ncells_local * this%dof) then
       call endrun(msg="ERROR size of vector /= number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    end if

    call VecGetArrayF90(x, x_p, ierr); CHKERRQ(ierr)

    do ghosted_id = 1, this%mesh%ncells_local
       this%aux_vars_in(ghosted_id)%Iup  = x_p((ghosted_id-1)*this%dof + 1)
       this%aux_vars_in(ghosted_id)%Idn  = x_p((ghosted_id-1)*this%dof + 2)
       !this%aux_vars_in(ghosted_id)%Iabs = x_p((ghosted_id-1)*this%dof + 3)
    end do

    call VecRestoreArrayF90(x, x_p, ierr); CHKERRQ(ierr)

  end subroutine LongwaveSavePrmIndepVar

  !------------------------------------------------------------------------

  subroutine LongwaveUpdateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Updates auxiliary variable associated with internal control volumes
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this
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

  end subroutine LongwaveUpdateAuxVars

  !------------------------------------------------------------------------
  subroutine LongwaveComputeRhs(this, B, ierr)
    !
    ! !DESCRIPTION:
    !
    use ConditionType             , only : condition_type
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use MultiPhysicsProbConstants , only : STEFAN_BOLTZMAN_CONSTANT
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this
    Vec                                   :: B
    PetscErrorCode                        :: ierr
    !
    type(longwave_auxvar_type) , pointer  :: avars(:)
    class(connection_set_type) , pointer  :: cur_conn_set
    class(condition_type)      , pointer  :: cur_cond
    PetscScalar                , pointer  :: b_p(:)
    PetscInt                              :: icell    , ileaf, iconn
    PetscInt                              :: sum_conn
    PetscInt                              :: cell_id_up, cell_id_dn
    PetscInt                              :: cell_i, cell_i_plus_1
    PetscReal                             :: b_i, e_i

    avars => this%aux_vars_in

    call VecGetArrayF90(B, b_p, ierr)

    ! Set RHS for upwind and absorbed radiation
    do icell = 1, this%mesh%ncells_local

       if (avars(icell)%is_soil) then

          b_p((icell-1)*this%dof + 1) = avars(icell)%rad_source ! upward
          !b_p((icell-1)*this%dof + 3) = 0.d0                    ! absorbed

       else
          !
          ! Upward flux: c = (1-e_i) * (1-tau_{d,i}) * emiss * sigma * T^4
          b_p((icell-1)*this%dof + 1) = (1.d0 - avars(icell)%e ) * avars(icell)%rad_source ! upward

          ! Absorbed flux: -2 * (1-tau_{d,i}) * emiss * sigma * T^4
          !b_p((icell-1)*this%dof + 3) = -2.d0                    * avars(icell)%rad_source ! absorbed
       end if
    end do

    ! Set RHS for downward radiation for all cells except top of the canopy
    cur_conn_set => this%mesh%intrn_conn_set_list%first
    sum_conn   = 0
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(sum_conn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(sum_conn)%GetIDDn()

          call DetermineCellIDsTopAndBottom(this, cell_id_up, cell_id_dn, cell_i, cell_i_plus_1);

          ! Downward flux: d_i = (1 - b_i    ) * (...)
          !                    = (1 - e_{i+1}) * (...)
          icell = cell_i
          b_p((icell-1)*this%dof + 2) = (1.d0 - avars(cell_i_plus_1)%e) * avars(cell_i_plus_1)%rad_source

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    ! Set RHS for downward radiation at the top of the canopy
    sum_conn = 0
    cur_cond => this%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          icell = cur_conn_set%conn(iconn)%GetIDDn()

          if ( (.not. this%mesh%is_active(icell))) cycle

          b_p((icell-1)*this%dof + 2) = this%aux_vars_bc(sum_conn)%Idn ! downward

       enddo
       cur_cond => cur_cond%next
    enddo

    call VecRestoreArrayF90(B, b_p, ierr)

  end subroutine LongwaveComputeRhs


  !------------------------------------------------------------------------
  subroutine DetermineCellIDsTopAndBottom(this, icell_1, icell_2, icell_bot, icell_top)
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this
    PetscInt :: icell_1, icell_2
    PetscInt :: icell_top, icell_bot

    if (this%mesh%z(icell_1) > this%mesh%z(icell_2)) then
       icell_top = icell_1
       icell_bot = icell_2
    else
       icell_top = icell_2
       icell_bot = icell_1
    end if

  end subroutine DetermineCellIDsTopAndBottom

  !------------------------------------------------------------------------

  subroutine LongwaveComputeOperatorsDiag(this, A, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_longwave_type) :: this
    Mat                          :: A
    Mat                          :: B
    PetscErrorCode               :: ierr
    !
    ! !LOCAL VARIABLES
    class(connection_set_type) , pointer :: cur_conn_set
    type(condition_type)       , pointer :: cur_cond
    type(longwave_auxvar_type) , pointer :: avars(:)
    PetscInt                             :: icell, ileaf, iconn, idof, sum_conn
    PetscInt                             :: cell_id, cell_id_up, cell_id_dn
    PetscInt                              :: cell_i, cell_i_plus_1
    PetscInt                             :: row, col
    PetscReal                            :: value, ga, dist, gsw, gsa, gs0
    PetscReal                            :: qsat, si, lambda
    PetscBool                            :: found

    avars => this%aux_vars_in

    do icell = 1, this%mesh%ncells_local

       do idof = 1, this%dof
          row = (icell-1)*this%dof + idof - 1;
          col = row
          value = 1.d0
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr);
       end do

       if (avars(icell)%is_soil) then

          value = avars(icell)%f
          row = (icell-1)*this%dof
          col = row + 1
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr);

       end if
    end do

    cur_conn_set => this%mesh%intrn_conn_set_list%first
    sum_conn   = 0
    do
       if (.not.associated(cur_conn_set)) exit

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1

          cell_id_up = cur_conn_set%conn(sum_conn)%GetIDUp()
          cell_id_dn = cur_conn_set%conn(sum_conn)%GetIDDn()

          call DetermineCellIDsTopAndBottom(this, cell_id_up, cell_id_dn, cell_i, cell_i_plus_1);

          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! Insert a_i     for the i-th     downward flux
          ! Insert f_{i+1} for the {i+1}-th upward flux
          !
          ! Note: a_i = f_{i+1}
          value = -avars(cell_i_plus_1)%f

          row = (cell_i - 1)*this%dof + 1
          col = (cell_i - 1)*this%dof
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = (cell_i_plus_1 - 1)*this%dof
          col = (cell_i_plus_1 - 1)*this%dof + 1
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! Insert b_i     for the i-th     downward flux
          ! Insert e_{i+1} for the {i+1}-th upward flux
          !
          ! Note: b_i = e_{i+1}
          value = -avars(cell_i_plus_1)%e

          row = (cell_i        - 1)*this%dof + 1
          col = (cell_i_plus_1 - 1)*this%dof
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

          row = (cell_i_plus_1 - 1)*this%dof
          col = (cell_i        - 1)*this%dof + 1
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  end subroutine LongwaveComputeOperatorsDiag
#endif

end module GoveqnLongwaveType
