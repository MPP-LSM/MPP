module ml_model_problem

  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbShortwave      , only : mpp_shortwave_type
  use MultiPhysicsProbLongwave       , only : mpp_longwave_type
  use MultiPhysicsProbLBL            , only : mpp_lbl_type
  use MultiPhysicsProbPhotosynthesis , only : mpp_photosynthesis_type
  use MultiPhysicsProbMLC            , only : mpp_mlc_type
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
    use ml_model_global_vars, only : ncair, ntree
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
    use ml_model_global_vars, only : ncair, ntree
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
    use swv            , only : init_swv
    use lwv            , only : init_lwv
    use lbl            , only : init_lbl
    use photosynthesis , only : init_photosynthesis
    use mlc            , only : init_mlc
    !
    implicit none

    call init_swv(swv_mpp)
    call init_lwv(lwv_mpp)
    call init_lbl(lbl_mpp)
    call init_photosynthesis(psy_mpp)
    call init_mlc(mlc_mpp)

  end subroutine init_mpps

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    use mlc            , only : mlc_set_initial_conditions
    !
    implicit none

    call mlc_set_initial_conditions(mlc_mpp)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine run_ml_model_problem(namelist_filename)
    !
    use ml_model_boundary_conditions , only : read_boundary_conditions, allocate_memory
    use ml_model_utils               , only : compute_dpai, compute_fssh_and_cumlai
    use mlc                          , only : extract_data_from_mlc
    use swv                          , only : extract_data_from_swv
    use lbl                          , only : extract_data_from_lbl
    use lwv                          , only : extract_data_from_lwv
    use photosynthesis               , only : extract_data_from_photosynthesis
    use swv                          , only : solve_swv
    use lwv                          , only : solve_lwv
    use lbl                          , only : solve_lbl
    use photosynthesis               , only : solve_photosynthesis
    use mlc                          , only : solve_mlc
    use ml_model_global_vars         , only : dpai, fssh, cumlai, sumpai, ncair, ntree, nz_cair, nbot, ntop
    !
    implicit none
    !
    character(len=256), optional :: namelist_filename
    !
    PetscReal                    :: dt
    PetscBool                    :: converged
    PetscInt                     :: istep
    PetscInt                     :: converged_reason

    ncair = 1;
    ntree = 1;

    call read_command_options()
    call read_namelist_file(namelist_filename)
    allocate(dpai  (nz_cair*ntree +1))
    allocate(fssh  (nz_cair*ntree +1))
    allocate(cumlai(nz_cair*ntree +1))
    allocate(sumpai(nz_cair*ntree +1))

    call compute_dpai(dpai, cumlai, sumpai)
    call compute_fssh_and_cumlai(nbot, ntop, dpai, fssh, cumlai, sumpai)

    call allocate_memory()
    call init_mpps()

    call set_initial_conditions()

    do istep = 1, 1
       call read_boundary_conditions(istep)

       dt = 30.d0 * 60.d0 ! [sec]
       write(*,*)'Solving shortwave radiation'
       call solve_swv(swv_mpp, istep, dt)
       call extract_data_from_swv(swv_mpp)

       call extract_data_from_mlc(mlc_mpp)
       write(*,*)'Solving longwave radiation'
       call solve_lwv(lwv_mpp, istep, dt)
       call extract_data_from_lwv(lwv_mpp)

       write(*,*)'Solving leaf boundary layer'
       call solve_lbl(lbl_mpp, istep, dt)
       call extract_data_from_lbl(lbl_mpp)

       write(*,*)'Solving photosynthesis'
       call solve_photosynthesis(psy_mpp, istep, dt)
       call extract_data_from_photosynthesis(psy_mpp)

       write(*,*)'Solving MLC'
       call solve_mlc(mlc_mpp, istep, dt)
    end do

  end subroutine run_ml_model_problem

end module ml_model_problem
