module ml_model_problem

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petsclog.h>
  
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
    use ml_model_global_vars      , only : ncair, ntree, output_data, bc_file, mlc_ic_file, psy_ic_file, use_ic
    use ml_model_global_vars      , only : checkpoint_data, beg_step, end_step, nsubstep, gstype
    use MultiPhysicsProbConstants , only : VAR_SCM_MEDLYN
    use MultiPhysicsProbConstants , only : VAR_SCM_BBERRY
    use MultiPhysicsProbConstants , only : VAR_SCM_WUE
    use MultiPhysicsProbConstants , only : VAR_SCM_BONAN14
    use MultiPhysicsProbConstants , only : VAR_SCM_MODIFIED_BONAN14
    use MultiPhysicsProbConstants , only : VAR_SCM_MANZONI11
    use MultiPhysicsProbConstants , only : VAR_SCM_OSMWANG
    !
    implicit none
    !
    character(len=256) :: stomatal_conductance_model
    PetscBool          :: flg
    PetscErrorCode     :: ierr
    PetscBool          :: mlc_ic_file_specified, psy_ic_file_specified

    mlc_ic_file_specified = PETSC_FALSE
    psy_ic_file_specified = PETSC_FALSE

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ncair',ncair,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ntree',ntree,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-beg_step',beg_step,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-end_step',end_step,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nsubstep',nsubstep,flg,ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-output_data',output_data,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-checkpoint_data',checkpoint_data,flg,ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-mlc_ic_file',mlc_ic_file,flg,ierr)
    if (flg) mlc_ic_file_specified = PETSC_TRUE

    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-photosynthesis_ic_file',psy_ic_file,flg,ierr)
    if (flg) psy_ic_file_specified = PETSC_TRUE

    if (mlc_ic_file_specified .and. (.not. psy_ic_file_specified)) then
       write(*,*)'ERROR: IC file was specified for the MLC model, but not for the Photosynthesis model.'
       write(*,*)'       Need to specify IC for both models'
       call exit(-1)
    endif

    if ((.not. mlc_ic_file_specified) .and. psy_ic_file_specified) then
       write(*,*)'ERROR: IC file was specified for the Photosynthesis model, but not for the MLC model.'
       write(*,*)'       Need to specify IC for both models'
       call exit(-1)
    endif

    if (mlc_ic_file_specified .and. psy_ic_file_specified) then
       use_ic = PETSC_TRUE
    endif

    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-bc_file',bc_file,flg,ierr)
    if (.not.flg) then
       write(*,*)'ERROR: Need to specify the boundary condition file via -bc_file <filename>'
       call exit(-1)
    end if

    stomatal_conductance_model = 'medlyn'
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-stomatal_conductance_model', stomatal_conductance_model, flg, ierr)

    select case(trim(stomatal_conductance_model))
    case ('ball-berry')
       gstype = VAR_SCM_BBERRY
    case ('medlyn')
       gstype = VAR_SCM_MEDLYN
    case ('wue ')
       gstype = VAR_SCM_WUE
    case ('bonan14')
       gstype = VAR_SCM_BONAN14
    case ('modified_bonan14')
       gstype = VAR_SCM_MODIFIED_BONAN14
    case ('manzoni11')
       gstype = VAR_SCM_MANZONI11
    case ('osmwang')
       gstype = VAR_SCM_OSMWANG
    case default
       write(iulog,*) 'Invalid value for -stomatal_conductance_model and valid values are: '
       write(iulog,*) '  ball-berry'
       write(iulog,*) '  medlyn'
       write(iulog,*) '  wue'
       write(iulog,*) '  bonan14'
       write(iulog,*) '  modified_bonan14'
       write(iulog,*) '  manzoni'
       write(iulog,*) '  osmwang'
       call exit(-1)
    end select

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
    use mlc                       , only : mlc_set_initial_conditions
    use ml_model_global_vars      , only : nbot, ntop, ncair, ntree, nz_cair
    use ml_model_global_vars      , only : bnd_cond, int_cond
    use ml_model_utils            , only : get_value_from_condition
    use ml_model_utils            , only : set_value_in_condition
    use MultiPhysicsProbConstants , only : MM_H2O, MM_DRY_AIR
    !
    implicit none
    !
    PetscInt :: icair, itree, k, idx_leaf, idx_air
    PetscReal :: tleaf_value, tair_value, wind_value, qair_value, factor


    idx_leaf = 0
    idx_air  = 0
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, nz_cair+1
             if (k >= nbot .and. k <= ntop) then
                tleaf_value = get_value_from_condition(bnd_cond%tref, icair)

                idx_leaf = idx_leaf + 1
                call set_value_in_condition(int_cond%Tleaf_sun, idx_leaf, tleaf_value)
                call set_value_in_condition(int_cond%Tleaf_shd, idx_leaf, tleaf_value)

             end if

             if (k > 1) then
                tair_value = get_value_from_condition(bnd_cond%tref, icair)
                wind_value = get_value_from_condition(bnd_cond%uref, icair)
                qair_value = get_value_from_condition(bnd_cond%qref, icair)
                factor = 1.d0/(MM_H2O/MM_DRY_AIR + (1.d0 - MM_H2O/MM_DRY_AIR)*qair_value)

                idx_air = (icair-1)*ncair + (k-1)
                call set_value_in_condition(int_cond%Tair, idx_air, tair_value)
                call set_value_in_condition(int_cond%Wind, idx_air, wind_value)
                call set_value_in_condition(int_cond%qair, idx_air, qair_value*factor)
             end if

          end do
       end do
    end do

    call mlc_set_initial_conditions(mlc_mpp)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine initialize_from_checkpoint()
    !
    use mlc                       , only : mlc_initialize_from_checkpoint
    use photosynthesis            , only : photosynthesis_initialize_from_checkpoint
    use ml_model_global_vars      , only : nbot, ntop, ncair, ntree, nz_cair, mlc_ic_file, psy_ic_file
    use ml_model_global_vars      , only : bnd_cond, int_cond
    use ml_model_utils            , only : get_value_from_condition
    use ml_model_utils            , only : set_value_in_condition
    use MultiPhysicsProbConstants , only : MM_H2O, MM_DRY_AIR
    !
    implicit none
    !
    PetscInt           :: icair, itree, k, idx_leaf, idx_air
    PetscReal          :: tsun_value, tsha_value, tair_value, wind_value, qair_value, factor
    Vec                :: ic_data
    PetscInt           :: offset
    PetscViewer        :: viewer
    PetscErrorCode     :: ierr
    PetscReal, pointer :: ic_p(:)

    ! Load boundary condition data
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, mlc_ic_file, FILE_MODE_READ, viewer, ierr);CHKERRQ(ierr)
    call VecCreate(PETSC_COMM_WORLD, ic_data, ierr);CHKERRQ(ierr)
    call VecLoad(ic_data, viewer, ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)

    call VecGetArrayF90(ic_data, ic_p, ierr); CHKERRA(ierr)

    idx_leaf = 0
    idx_air  = 0
    icair = 1;
    do k = 1, nz_cair+1
       if (k > 1) then
          idx_air = (icair-1)*ncair + (k-1)
          offset = 0            ; wind_value = ic_p(k + offset);
          offset = (nz_cair+1)  ; tair_value = ic_p(k + offset);
          offset = (nz_cair+1)*2; qair_value = ic_p(k + offset);
          offset = (nz_cair+1)*3; tsun_value = ic_p(k + offset);
          offset = (nz_cair+1)*4; tsha_value = ic_p(k + offset);

          call set_value_in_condition(int_cond%Tair, idx_air, tair_value)
          call set_value_in_condition(int_cond%Wind, idx_air, wind_value)
          call set_value_in_condition(int_cond%qair, idx_air, qair_value)
          if (k >= nbot .and. k <= ntop) then
             idx_leaf = idx_leaf + 1
             call set_value_in_condition(int_cond%Tleaf_sun, idx_leaf, tsun_value)
             call set_value_in_condition(int_cond%Tleaf_shd, idx_leaf, tsha_value)
          end if
       end if
    end do

    call VecRestoreArrayF90(ic_data, ic_p, ierr); CHKERRA(ierr)
    call VecDestroy(ic_data, ierr); CHKERRA(ierr)

    call mlc_initialize_from_checkpoint(mlc_mpp)
    call photosynthesis_initialize_from_checkpoint(psy_mpp, psy_ic_file)

  end subroutine initialize_from_checkpoint

  !------------------------------------------------------------------------
  subroutine run_ml_model_problem(namelist_filename)
    !
    use ml_model_boundary_conditions , only : read_boundary_conditions, allocate_memory
    use ml_model_utils               , only : compute_vertical_veg_structure, compute_fssh
    use mlc                          , only : extract_data_from_mlc
    use swv                          , only : extract_data_from_swv
    use lbl                          , only : extract_data_from_lbl
    use lwv                          , only : extract_data_from_lwv
    use photosynthesis               , only : extract_data_from_photosynthesis
    use swv                          , only : solve_swv
    use lwv                          , only : solve_lwv
    use lbl                          , only : solve_lbl
    use photosynthesis               , only : solve_photosynthesis, checkpoint_photosynthesis
    use mlc                          , only : solve_mlc, checkpoint_mlc
    use ml_model_global_vars         , only : dsai, dlai, dpai, fssh, cumpai, sumpai, leaf_td, ncair, ntree, nz_cair, nbot, ntop
    use ml_model_global_vars         , only : nz_cair, ntree, output_data, bc_file, checkpoint_data, use_ic
    use ml_model_global_vars         , only : beg_step, end_step, nsubstep
    use petscvec
    !
    implicit none
    !
    character(len=256), optional :: namelist_filename
    !
    PetscReal      :: dt
    PetscBool      :: converged
    PetscInt       :: istep, isubstep
    PetscInt       :: converged_reason
    Vec            :: bc_data
    PetscViewer    :: viewer
    PetscErrorCode :: ierr
    PetscLogEvent :: event_swv
    PetscLogEvent :: event_lwv
    PetscLogEvent :: event_lbl
    PetscLogEvent :: event_phy
    PetscLogEvent :: event_mlc
    PetscClassId  :: classid
    PetscBool     :: is_first_substep

    classid = 0
    call PetscLogEventRegister('SWV', classid, event_swv, ierr)
    call PetscLogEventRegister('LWV', classid, event_lwv, ierr)
    call PetscLogEventRegister('LBL', classid, event_lbl, ierr)
    call PetscLogEventRegister('PHY', classid, event_phy, ierr)
    call PetscLogEventRegister('MLC', classid, event_mlc, ierr)

    ncair           = 1;
    ntree           = 1;
    beg_step        = 1
    end_step        = 1
    nsubstep        = 12
    output_data     = PETSC_FALSE
    checkpoint_data = PETSC_FALSE
    use_ic          = PETSC_FALSE

    call read_command_options()
    call read_namelist_file(namelist_filename)

    call compute_vertical_veg_structure(dlai, dsai, dpai, cumpai, sumpai, leaf_td)
    allocate(fssh(nz_cair*ntree + 1))

    call allocate_memory()
    call init_mpps()

    ! Load boundary condition data
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, bc_file, FILE_MODE_READ, viewer, ierr);CHKERRQ(ierr)
    call VecCreate(PETSC_COMM_WORLD, bc_data, ierr);CHKERRQ(ierr)
    call VecLoad(bc_data, viewer, ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)

    dt = 3600.d0/nsubstep
    istep = beg_step
    call read_boundary_conditions(1, bc_data)

    if (.not.use_ic) then
       call set_initial_conditions()
    else
       call initialize_from_checkpoint()
    end if

    do istep = beg_step, end_step
       call read_boundary_conditions(1, bc_data)
       write(*,*)'%istep: ',istep

       write(*,*)'%  Solving shortwave radiation'
       call PetscLogEventBegin(event_swv, ierr)
       call solve_swv(swv_mpp, istep, dt)
       call PetscLogEventEnd(event_swv, ierr)

       do isubstep = 1, nsubstep

          write(*,*)'%    isubstep:',isubstep
          dt = 300.d0 ! [sec]

          write(*,*)'%    Solving longwave radiation'
          call PetscLogEventBegin(event_lwv, ierr)
          call solve_lwv(lwv_mpp, istep, isubstep, dt)
          call PetscLogEventEnd(event_lwv, ierr)

          write(*,*)'%    Solving leaf boundary layer'
          call PetscLogEventBegin(event_lbl, ierr)
          call solve_lbl(lbl_mpp, istep, isubstep,  dt)
          call PetscLogEventEnd(event_lbl, ierr)

          write(*,*)'%    Solving photosynthesis'
          call PetscLogEventBegin(event_phy, ierr)
          if (istep == beg_step .and. isubstep == 1) then
            is_first_substep = PETSC_TRUE
          else
            is_first_substep = PETSC_FALSE
          endif
          call solve_photosynthesis(psy_mpp, istep, isubstep, is_first_substep, dt)
          call PetscLogEventEnd(event_phy, ierr)

          write(*,*)'%    Solving MLC'
          call PetscLogEventBegin(event_mlc, ierr)
          call solve_mlc(mlc_mpp, istep, isubstep, dt)
          call PetscLogEventEnd(event_mlc, ierr)
          write(*,*)''
       end do

       if (checkpoint_data) then
          call checkpoint_mlc(mlc_mpp, istep, isubstep-1)
          call checkpoint_photosynthesis(psy_mpp, istep, isubstep-1)
       end if
    end do

  end subroutine run_ml_model_problem

end module ml_model_problem
