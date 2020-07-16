module AuxVarType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>
  use petscvec
  use petscmat

  ! !USES:
  use mpp_varctl         , only : iulog
  use mpp_abortutils     , only : endrun
  use mpp_shr_log_mod    , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: auxvar_base_type
   contains
     procedure, public :: Init          => AuxVarBaseInit
     procedure, public :: SetValue      => AuxVarBaseSetValue
     procedure, public :: GetValue      => AuxVarBaseGetValue
     procedure, public :: AuxVarCompute => AuxVarBaseAuxVarCompute
  end type auxvar_base_type

contains

  !------------------------------------------------------------------------
  subroutine AuxVarBaseInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(auxvar_base_type) :: this

    write(iulog,*)'AuxVarBaseInit Must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine AuxVarBaseInit

  !------------------------------------------------------------------------
  subroutine AuxVarBaseSetValue(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(auxvar_base_type) , intent(inout) :: this
    PetscInt                , intent(in)    :: var_type
    PetscReal               , intent(in)    :: variable_value

    write(iulog,*)'AuxVarBaseSetValue Must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine AuxVarBaseSetValue

  !------------------------------------------------------------------------
  subroutine AuxVarBaseGetValue(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(auxvar_base_type) :: this
    PetscInt  , intent(in)  :: var_type
    PetscReal , intent(out) :: variable_value

    write(iulog,*)'AuxVarBaseGetValue Must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine AuxVarBaseGetValue

  !------------------------------------------------------------------------
  subroutine AuxVarBaseAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(auxvar_base_type) :: this

    write(iulog,*)'AuxVarBaseAuxVarCompute Must be extended by child class.'
    call endrun(msg=errMsg(__FILE__, __LINE__))

  end subroutine AuxVarBaseAuxVarCompute

#endif

end module AuxVarType
