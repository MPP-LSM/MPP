module ml_model_problem

  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbShortwave      , only : mpp_shortwave_type
  use MultiPhysicsProbLongwave       , only : mpp_longwave_type
  use MultiPhysicsProbLBL            , only : mpp_lbl_type
  use MultiPhysicsProbPhotosynthesis , only : mpp_photosynthesis_type
  use MultiPhysicsProbMLC            , only : mpp_mlc_type
  use ml_model_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none

  type(mpp_shortwave_type)      :: swv_mpp
  type(mpp_longwave_type)       :: lwv_mpp
  type(mpp_lbl_type)            :: lbl_mpp
  type(mpp_photosynthesis_type) :: psy_mpp
  type(mpp_mlc_type)            :: mlc_mpp

#include <petsc/finclude/petsc.h>

  public :: run_ml_model_problem

contains

  !------------------------------------------------------------------------
  subroutine read_namelist_file(namelist_filename)
    !
    implicit none
    !
    character(len=256), optional :: namelist_filename
    !
    character(len=256)           :: ioerror_msg
    character(len=2560)          :: namelist_buffer
    integer                      :: nml_unit, nml_error

    namelist / problem_options / ncair, ntree

    if (present(namelist_filename)) then
       nml_unit = 16
       open(unit=nml_unit, file=trim(namelist_filename), action='read', access='stream', &
            form='unformatted', iostat=nml_error)
       if (nml_error /= 0) then
          write(*,*)'ERROR: Unable to open namelist file: ',trim(namelist_filename)
          call exit(-1)
       endif

       read(unit=nml_unit, iostat=nml_error, iomsg=ioerror_msg) namelist_buffer
       if (.not. is_iostat_end(nml_error)) then
          write(*,*)"ERROR: Unable to read '",trim(namelist_filename),"' till EOF"
          call exit(-1)
       endif

       read(namelist_buffer, nml=problem_options, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          write(*,*)'ERROR: Unable to read "problem_options" in namelist file '
          call exit(-1)
       endif

       close(nml_unit)
    endif

  end subroutine read_namelist_file

  !------------------------------------------------------------------------
  subroutine read_command_options()
    !
    implicit none
    !
    PetscBool      :: flg
    PetscErrorCode :: ierr

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ncair',ncair,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ntree',ntree,flg,ierr)

  end subroutine read_command_options

  !------------------------------------------------------------------------
  subroutine init_mpps()
    !
    use ml_model_utils , only: compute_dpai_fssh
    use swv            , only : init_swv
    use lwv            , only : init_lwv
    !
    implicit none

    call compute_dpai_fssh()

    call init_swv(swv_mpp)
    call init_lwv(lwv_mpp)
    !call init_lbl(lbl_mpp)
    !call init_psy(psy_mpp)
    !call init_mlc(mlc_mpp)

  end subroutine init_mpps

  !------------------------------------------------------------------------
  subroutine run_ml_model_problem(namelist_filename)
    !
    implicit none
    !
    PetscReal                    :: dt
    PetscBool                    :: converged
    PetscInt                     :: istep
    PetscInt                     :: converged_reason
    character(len=256), optional :: namelist_filename

    ncair = 1;
    ntree = 1;

    call read_command_options()
    call read_namelist_file(namelist_filename)

    call init_mpps()

  end subroutine run_ml_model_problem

end module ml_model_problem
