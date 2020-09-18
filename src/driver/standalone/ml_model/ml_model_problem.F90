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
    use ml_model_utils , only : compute_dpai_fssh
    use mlc            , only : mlc_set_initial_conditions
    !
    implicit none

    call mlc_set_initial_conditions(mlc_mpp)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine allocate_memory()
    !
    use ml_model_global_vars, only : condition_type
    use ml_model_utils      , only : allocate_memory_for_condition
    !
    implicit none

    call allocate_memory_for_condition(Iskyb_vis , ncair)
    call allocate_memory_for_condition(Iskyd_vis , ncair)
    call allocate_memory_for_condition(Iskyb_nir , ncair)
    call allocate_memory_for_condition(Iskyd_nir , ncair)
    call allocate_memory_for_condition(Irsky     , ncair)

    call allocate_memory_for_condition(Pref  , ncair)
    call allocate_memory_for_condition(Uref  , ncair)
    call allocate_memory_for_condition(Tref  , ncair)
    call allocate_memory_for_condition(Rhref , ncair)

    call allocate_memory_for_condition(gbh , ncair)
    call allocate_memory_for_condition(gbv , ncair)
    call allocate_memory_for_condition(gbc , ncair)

    call allocate_memory_for_condition(Tcan      , ncair)

    call allocate_memory_for_condition(Tair      , ncair*nz_cair)
    call allocate_memory_for_condition(eair      , ncair*nz_cair)

    call allocate_memory_for_condition(Tleaf_sun , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(Tleaf_shd , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(rn_sun    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(rn_shd    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(gs_sun    , ncair*(ntop-nbot+1) )
    call allocate_memory_for_condition(gs_shd    , ncair*(ntop-nbot+1) )
    
    call allocate_memory_for_condition(Tsoil  , ncair)
    call allocate_memory_for_condition(rn_soil, ncair)

  end subroutine allocate_memory

  !------------------------------------------------------------------------
  subroutine run_ml_model_problem(namelist_filename)
    !
    use ml_model_boundary_conditions , only : read_boundary_conditions
    use ml_model_utils               , only : compute_dpai_fssh
    use ml_model_utils               , only : extract_data_from_mlc
    use swv                          , only : solve_swv
    use lwv                          , only : solve_lwv
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

    call compute_dpai_fssh()

    call allocate_memory()
    call init_mpps()

    call set_initial_conditions()

    do istep = 1, 1
       call read_boundary_conditions(istep)

       dt = 30.d0 * 60.d0 ! [sec]
       call solve_swv(swv_mpp, istep, dt)

       call extract_data_from_mlc(mlc_mpp)
       call solve_lwv(lwv_mpp, istep, dt)
    end do

  end subroutine run_ml_model_problem

end module ml_model_problem
