module MultiPhysicsProbShortwave

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

use mpp_varctl                       , only : iulog
use mpp_abortutils                   , only : endrun
use mpp_shr_log_mod                  , only : errMsg => shr_log_errMsg
use MultiPhysicsProbBaseType         , only : multiphysicsprob_base_type
use SystemOfEquationsShortwaveType   , only : sysofeqns_shortwave_type
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

type, public, extends(multiphysicsprob_base_type) :: mpp_shortwave_type
   
contains
  procedure, public :: Init            => MppShortwaveRadiationInit
  procedure, public :: AllocateAuxVars => MppShortwaveRadiationAllocateAuxVars
  procedure, public :: SetupProblem    => MppShortwaveRadiationSetupProblem
end type mpp_shortwave_type

contains

  !------------------------------------------------------------------------
subroutine MppShortwaveRadiationInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize the MLC MPP
    !
    use MultiPhysicsProbBaseType , only : MPPBaseInit
    use MultiPhysicsProbConstants , only : SOE_LONGWAVE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_shortwave_type) :: this
    !
    class(sysofeqns_shortwave_type), pointer :: sysofeqns

    call MPPBaseInit(this)

    allocate(sysofeqns)
    call sysofeqns%Init()

    this%soe => sysofeqns
    this%soe%itype = SOE_LONGWAVE

    allocate(this%soe_ptr)
    nullify(this%soe_ptr%ptr)

  end subroutine MppShortwaveRadiationInit

  !------------------------------------------------------------------------
  subroutine MPPShortwaveRadiationAllocateAuxVars(this)
   !
   ! !DESCRIPTION:
   ! Sets the PETSc KSP problem
   !
   use GoverningEquationBaseType , only : goveqn_base_type
   !
   implicit none
   !
   ! !ARGUMENTS
   class(mpp_shortwave_type) :: this
   !
   class(goveqn_base_type), pointer  :: cur_goveq

   call this%soe%AllocateAuxVars()

   cur_goveq => this%soe%goveqns
   do
      if (.not.associated(cur_goveq)) exit

      call cur_goveq%AllocateAuxVars()

      cur_goveq => cur_goveq%next
   end do

 end subroutine MPPShortwaveRadiationAllocateAuxVars

!------------------------------------------------------------------------
subroutine MppShortwaveRadiationSetupProblem(this)
    !
    ! !DESCRIPTION:
    ! Initialize the MLC MPP
    !
    use MultiPhysicsProbBaseType, only : MPPSetupProblem
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_shortwave_type) :: this
    
    call MPPSetupProblem(this)

  end subroutine MppShortwaveRadiationSetupProblem

#endif

end module MultiPhysicsProbShortwave
