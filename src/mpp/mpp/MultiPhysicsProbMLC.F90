module MultiPhysicsProbMLC

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for thermal model
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                       , only : iulog
  use mpp_abortutils                   , only : endrun
  use mpp_shr_log_mod                  , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbBaseType         , only : multiphysicsprob_base_type
  use SystemOfEquationsMLCType         , only : sysofeqns_mlc_type
  use SystemOfEquationsBasePointerType , only : sysofeqns_base_pointer_type
  use SystemOfEquationsBaseType        , only : sysofeqns_base_type
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda

  implicit none
  private

  type, public, extends(multiphysicsprob_base_type) :: mpp_mlc_type
   contains
     procedure, public :: Init            => MppMlcInit
     procedure, public :: AllocateAuxVars => MppMlcAllocateAuxVars
     procedure, public :: SetupProblem    => MppMlcSetupProblem

  end type mpp_mlc_type

contains

  !------------------------------------------------------------------------
  subroutine MppMlcInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize the MLC MPP
    !
    use MultiPhysicsProbBaseType , only : MPPBaseInit
    use MultiPhysicsProbConstants , only : SOE_MLC
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_mlc_type) :: this
    !
    class(sysofeqns_mlc_type), pointer :: sysofeqns

    call MPPBaseInit(this)

    allocate(sysofeqns)
    call sysofeqns%Init()

    this%soe => sysofeqns
    this%soe%itype = SOE_MLC

    allocate(this%soe_ptr)
    nullify(this%soe_ptr%ptr)

  end subroutine MppMlcInit

  !------------------------------------------------------------------------
  subroutine MppMlcAllocateAuxVars(this, ncair)
    !
    ! !DESCRIPTION:
    ! Allocates auxvars for governing equations and system-of-governing-eqns
    !
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_NULL
    use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_mlc_type) :: this
    PetscInt            :: ncair
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: igoveqn, rank

    base_soe => this%soe
    
    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    soe%ncair = ncair
    allocate(soe%leaftemp_goveqn_rank(soe%ngoveqns))
    soe%leaftemp_goveqn_rank(:) = 0
    
    call soe%cturb%Init(ncair)

    ! Allocate AuxVars for each governing equation
    igoveqn = 0
    rank = 0
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       igoveqn = igoveqn + 1
       call cur_goveq%AllocateAuxVars()

       select type(cur_goveq)
       class is (goveqn_cleaf_temp_type)
          rank = rank + 1
          soe%leaftemp_goveqn_rank(igoveqn) = rank
       end select
       cur_goveq => cur_goveq%next
    end do

  end subroutine MppMlcAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine MppMlcSetupProblem(this)
    !
    ! !DESCRIPTION:
    ! Initialize the MLC MPP
    !
    use MultiPhysicsProbBaseType, only : MPPSetupProblem
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_mlc_type) :: this
    
    call MPPSetupProblem(this)

  end subroutine MppMlcSetupProblem

#endif

end module MultiPhysicsProbMLC
