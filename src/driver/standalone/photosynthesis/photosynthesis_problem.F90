module photosynthesis_problem

  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbPhotosynthesis , only : mpp_photosynthesis_type
  use photosynthesis_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none


#include <petsc/finclude/petsc.h>

  type(mpp_photosynthesis_type):: phtsyn_mpp

  public :: run_photosynthesis_problem
  public :: output_regression_photosynthesis_problem

contains

  !------------------------------------------------------------------------
  subroutine run_photosynthesis_problem(namelist_filename)
    !
    implicit none
    !
    character(len=256), optional :: namelist_filename
    !
    PetscReal          :: dt
    PetscBool          :: converged
    PetscInt           :: istep
    PetscInt           :: converged_reason
    PetscBool          :: flg
    PetscErrorCode     :: ierr

    ncair = 1;
    ntree = 1;

    call Init(phtsyn_mpp, namelist_filename)

    dt = 0.0
    istep = 1
    call phtsyn_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine run_photosynthesis_problem

  !------------------------------------------------------------------------
  subroutine Init(phtsyn_mpp, namelist_filename)
    !
    use photosynthesis_meshes     , only : setup_meshes
    use photosynthesis_parameters , only : set_parameters
    !
    implicit none
    !
    type(mpp_photosynthesis_type) :: phtsyn_mpp
    character(len=256), optional :: namelist_filename

    call read_command_line_options(namelist_filename)

    call initialize_mpp(phtsyn_mpp)

    call setup_meshes(phtsyn_mpp)

    call add_goveqn(phtsyn_mpp)

    call phtsyn_mpp%AllocateAuxVars()

    call phtsyn_mpp%SetupProblem()

    call set_parameters(phtsyn_mpp)

    call set_initial_condition(phtsyn_mpp)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine read_command_line_options(namelist_filename)
    !
    use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C3
    use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C4
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_MEDLYN
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_BBERRY
    use MultiPhysicsProbConstants , only : VAR_WUE
    !
    implicit none
    !
    character(len=256)           :: photosynthesis_pathway
    character(len=256)           :: stomatal_conductance_model
    PetscBool                    :: flag
    PetscErrorCode               :: ierr
    PetscOptions                 :: options
    character(len=256), optional :: namelist_filename
    character(len=256)           :: ioerror_msg
    character(len=2560)          :: namelist_buffer
    integer                      :: nml_unit, nml_error
    namelist / problem_options / photosynthesis_pathway, stomatal_conductance_model

    call PetscOptionsCreate(options, ierr); CHKERRQ(ierr);
    call PetscOptionsInsertString(options, '-photosynthesis_pathway     <c3|c4>', ierr); CHKERRQ(ierr)
    call PetscOptionsInsertString(options, '-stomatal_conductance_model <ball-berry|medlyn|wue>', ierr); CHKERRQ(ierr)
    call PetscOptionsView(options,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRQ(ierr)
    call PetscOptionsDestroy(options, ierr); CHKERRQ(ierr)

    photosynthesis_pathway     = 'c4'
    stomatal_conductance_model = 'medlyn'

    c3psn = VAR_PHOTOSYNTHETIC_PATHWAY_C4
    gstype = VAR_STOMATAL_CONDUCTANCE_MEDLYN

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-photosynthesis_pathway', photosynthesis_pathway, flag, ierr)

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-stomatal_conductance_model', stomatal_conductance_model, flag, ierr)

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

    select case(trim(photosynthesis_pathway))
    case ('c3')
       c3psn = VAR_PHOTOSYNTHETIC_PATHWAY_C3
    case ('c4')
       c3psn = VAR_PHOTOSYNTHETIC_PATHWAY_C4
    case default
       write(iulog,*) 'Invalid value for -photosynthesis_pathway and valid values are: c3, c4'
       call exit(0)
    end select

    select case(trim(stomatal_conductance_model))
    case ('ball-berry')
       gstype = VAR_STOMATAL_CONDUCTANCE_BBERRY
    case ('medlyn')
       gstype = VAR_STOMATAL_CONDUCTANCE_MEDLYN
    case ('wue')
       gstype = VAR_WUE
    case default
       write(iulog,*) 'Invalid value for -stomatal_conductance_model and valid values are: ball-berry, medlyn, wue'
       call exit(0)
    end select

  end subroutine read_command_line_options

 !------------------------------------------------------------------------
  subroutine initialize_mpp(phtsyn_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_PHOTOSYNTHESIS_SNES
    !
    implicit none

    type(mpp_photosynthesis_type) :: phtsyn_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call phtsyn_mpp%Init       ()
    call phtsyn_mpp%SetName    ('Photosynthesis model')
    call phtsyn_mpp%SetID      (MPP_PHOTOSYNTHESIS_SNES)
    call phtsyn_mpp%SetMPIRank (iam)

  end subroutine Initialize_Mpp

  !------------------------------------------------------------------------
  subroutine add_goveqn(phtsyn_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_PHOTOSYNTHESIS
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_photosynthesis_type) :: phtsyn_mpp

    PHTSYN_GE = 1

    call phtsyn_mpp%AddGovEqnWithMeshRank(GE_PHOTOSYNTHESIS, 'Photosynthesis model', PHTSYN_MESH)
    
    call phtsyn_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqn

  !------------------------------------------------------------------------
  subroutine set_initial_condition(phtsyn_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_PHOTOSYNTHESIS_SNES
    !
    implicit none

    type(mpp_photosynthesis_type) :: phtsyn_mpp
    PetscScalar, pointer          :: ci_p(:)
    PetscErrorCode                :: ierr

    call VecGetArrayF90(phtsyn_mpp%soe%solver%soln, ci_p, ierr); CHKERRQ(ierr)
    ci_p(:) = 152.d0
    call VecRestoreArrayF90(phtsyn_mpp%soe%solver%soln, ci_p, ierr); CHKERRQ(ierr)

  end subroutine set_initial_condition

 !------------------------------------------------------------------------
  subroutine output_regression_photosynthesis_problem(filename_base, num_cells)
    !
    use regression_mod            , only : regression_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnPhotosynthesisType  , only : goveqn_photosynthesis_type
    !
    implicit none
    !
    character(len=256)                :: filename_base
    PetscInt, intent(in)              :: num_cells
    !
    PetscInt                          :: output, ieqn, ncells, icell, k
    character(len=512)                :: filename
    character(len=64)                 :: category
    character(len=64)                 :: name
    PetscReal, pointer                :: ci(:)
    type(regression_type)             :: regression
    class(goveqn_base_type) , pointer :: goveq

    ncells = (nz_cair+1)*ncair

    allocate(ci(ncells))

    call regression%Init(filename_base, num_cells)
    call regression%OpenOutput()
    
    ieqn = 1

    call phtsyn_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

    select type(goveq)
    class is (goveqn_photosynthesis_type)

       do k = 1,nz_cair+1
          icell = k

          ci(icell) = goveq%aux_vars_in(icell)%ci

       end do

       category = 'general'
       name = 'ci'; call regression%WriteData(name, category, ci )
    end select

    call regression%CloseOutput()
    
    deallocate(ci)

  end subroutine output_regression_photosynthesis_problem

end module photosynthesis_problem
