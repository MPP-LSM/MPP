module GoveqnShortwaveType

#ifdef USE_PETSC_LIB
  !#define USE_BONAN_FORMULATION

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType , only : goveqn_base_type
  use ShortwaveAuxType           , only : shortwave_auxvar_type
  use petscvec
  use petscmat
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_shortwave_type

     type(shortwave_auxvar_type)      , pointer :: aux_vars_in(:)
     type(shortwave_auxvar_type)      , pointer :: aux_vars_bc(:)

     PetscInt                                   :: nLeaf

   contains

     procedure, public :: Setup                     => ShortwaveSetup
     procedure, public :: AllocateAuxVars           => ShortwaveAllocateAuxVars
     procedure, public :: PreSolve                  => ShortwavePreSolve
     procedure, public :: UpdateAuxVars             => ShortwaveUpdateAuxVars
     procedure, public :: SavePrimaryIndependentVar => ShortwaveSavePrmIndepVar
     !procedure, public :: ComputeRhs                => ShortwaveComputeRhs
     !procedure, public :: ComputeOperatorsDiag      => ShortwaveComputeOperatorsDiag

  end type goveqn_shortwave_type

contains
  !------------------------------------------------------------------------
  subroutine ShortwaveSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_SHORTWAVE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_shortwave_type) :: this

    call this%Create()

    this%name  = ""
    this%itype = GE_SHORTWAVE
    this%dof   = 2

    nullify(this%aux_vars_in)
    nullify(this%aux_vars_bc)

    this%nLeaf = 1

  end subroutine ShortwaveSetup

  !---------------------------------------------------------------------
  subroutine ShortwaveAllocateAuxVars(this)
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
    class(goveqn_shortwave_type)   :: this
    !
    type(condition_type) , pointer :: cur_cond
    PetscInt                       :: ghosted_id, ncells_cond, icond

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
  end subroutine ShortwaveAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine ShortwavePreSolve(this)
    !
    ! !DESCRIPTION:
    ! Presolve
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_shortwave_type) :: this

    call this%UpdateAuxVars()

  end subroutine ShortwavePreSolve

  !------------------------------------------------------------------------
  subroutine ShortwaveSavePrmIndepVar (this, x)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_shortwave_type) :: this
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
    end do

    call VecRestoreArrayF90(x, x_p, ierr); CHKERRQ(ierr)

  end subroutine ShortwaveSavePrmIndepVar

  !------------------------------------------------------------------------

  subroutine ShortwaveUpdateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Updates auxiliary variable associated with internal control volumes
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_shortwave_type) :: this
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

  end subroutine ShortwaveUpdateAuxVars

  !------------------------------------------------------------------------
  subroutine ShortwaveComputeRhs(this, R, ierr)
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
    class(goveqn_shortwave_type) :: this
    Vec                                   :: R
    PetscErrorCode                        :: ierr
    !
    type(shortwave_auxvar_type) , pointer :: avars(:)
    class(connection_set_type) , pointer  :: cur_conn_set
    class(condition_type)      , pointer  :: cur_cond
    PetscScalar                , pointer  :: r_p(:)
    PetscInt                              :: icell, iband, iconn, idx
    PetscInt                              :: sum_conn
    PetscInt                              :: cell_id_up, cell_id_dn
    PetscInt                              :: cell_i, cell_i_plus_1
    PetscReal                             :: rho, tau, b, e

    avars => this%aux_vars_in

    call VecGetArrayF90(R, r_p, ierr)

    ! Set RHS for upward
    do icell = 1, this%mesh%ncells_local

       if (avars(icell)%is_soil) then

          do iband = 1, avars(icell)%nband
             idx = (icell-1)*this%dof + (iband-1)*avars(icell)%nband + 1
             !c = rho_{g,h} * I_{sky,b} * T_{b,0}
             r_p(idx) = avars(icell)%rad_source(iband)
          end do

       else
          do iband = 1, avars(icell)%nband
             idx = (icell-1)*this%dof + (iband-1)*avars(icell)%nband + 1

             !c_i = I_{sky,b} * (1 - tau_{b,i}) * (rho_{l,i} - tau_{leaf,i} * e_i)
             !    = rad_source_{b,i})           * (rho_{l,i} - tau_{leaf,i} * e_i)

             rho = avars(icell)%leaf_rho(iband)
             tau = avars(icell)%leaf_tau(iband)
             e   = avars(icell)%e(iband)

             r_p(idx) = avars(icell)%rad_source(iband) * (rho - tau * e)
          end do
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

          icell = cell_i
          do iband = 1, avars(cell_id_up)%nband
             idx = (icell-1)*this%dof + (iband-1)*avars(icell)%nband + 2

             ! Downward flux: d_i = I_{sky,b} * (1 - tau_{b,i+1}) * (tau_{leaf,i+1} - rho_{leaf,i+1}*b_i)
             !                    = rad_source(i+1)               * (tau_{leaf,i+1} - rho_{leaf,1+1}*b_i)

             tau = avars(cell_i_plus_1)%leaf_tau(iband)
             rho = avars(cell_i_plus_1)%leaf_rho(iband)
             b   = avars(cell_i_plus_1)%e(iband)
             r_p(idx) = avars(cell_i_plus_1)%rad_source(iband) * (tau - rho * b)
          end do

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

          do iband = 1, avars(icell)%nband
             idx = (icell-1)*this%dof + (iband-1)*avars(icell)%nband + 2
             r_p(idx) = this%aux_vars_bc(sum_conn)%Iskyd(iband)             
          end do

       enddo
       cur_cond => cur_cond%next
    enddo

    call VecRestoreArrayF90(R, r_p, ierr)

  end subroutine ShortwaveComputeRhs


  !------------------------------------------------------------------------
  subroutine DetermineCellIDsTopAndBottom(this, icell_1, icell_2, icell_bot, icell_top)
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_shortwave_type) :: this
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

  subroutine ShortwaveComputeOperatorsDiag(this, A, B, ierr)
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
    class(goveqn_shortwave_type) :: this
    Mat                          :: A
    Mat                          :: B
    PetscErrorCode               :: ierr
    !
    ! !LOCAL VARIABLES
    class(connection_set_type)  , pointer :: cur_conn_set
    type(condition_type)        , pointer :: cur_cond
    type(shortwave_auxvar_type) , pointer :: avars(:)
    PetscInt                              :: icell, iband, nband, iconn, idof, sum_conn
    PetscInt                              :: cell_id, cell_id_up, cell_id_dn
    PetscInt                              :: cell_i, cell_i_plus_1
    PetscInt                              :: row, col
    PetscReal                             :: value, ga, dist, gsw, gsa, gs0
    PetscReal                             :: qsat, si, lambda
    PetscBool                             :: found

    avars => this%aux_vars_in
    nband = avars(1)%nband

    do icell = 1, this%mesh%ncells_local

       do idof = 1, this%dof
          row = (icell-1)*this%dof + idof - 1;
          col = row
          value = 1.d0
          call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr);
       end do

       if (avars(icell)%is_soil) then

          do iband = 1, nband
             value = -avars(icell)%f(iband)

             row = (icell-1)*this%dof + (iband-1)*nband
             col = row + 1
             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr);
          end do

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

          do iband = 1, nband
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             ! Insert a_i     for the i-th     downward flux
             ! Insert f_{i+1} for the {i+1}-th upward flux
             !
             ! Note: a_i = f_{i+1}
             value = -avars(cell_i_plus_1)%f(iband)

             row = (cell_i - 1)*this%dof + (iband-1)*nband + 1
             col = row - 1
             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

             row = (cell_i_plus_1 - 1)*this%dof + (iband-1)*nband
             col = row + 1
             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             ! Insert b_i     for the i-th     downward flux
             ! Insert e_{i+1} for the {i+1}-th upward flux
             !
             ! Note: b_i = e_{i+1}
             value = -avars(cell_i_plus_1)%e(iband)

             row = (cell_i        - 1)*this%dof + (iband-1)*nband + 1
             col = (cell_i_plus_1 - 1)*this%dof + (iband-1)*nband
             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

             row = (cell_i_plus_1 - 1)*this%dof + (iband-1)*nband
             col = (cell_i        - 1)*this%dof + (iband-1)*nband + 1
             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          enddo

       enddo

       cur_conn_set => cur_conn_set%next
    enddo

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  end subroutine ShortwaveComputeOperatorsDiag

#endif

end module GoveqnShortwaveType
