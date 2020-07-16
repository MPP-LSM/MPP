module longwave_problem

  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbLongwave , only : mpp_longwave_type
  use longwave_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none


#include <petsc/finclude/petsc.h>

  type(mpp_longwave_type):: longwave_mpp

  public :: run_longwave_problem
  public :: output_regression_longwave_problem

contains

  !------------------------------------------------------------------------
  subroutine run_longwave_problem(namelist_filename)

    use MultiPhysicsProbMLC, only : mpp_mlc_type
    !
    implicit none
    !
    PetscReal                    :: dt
    PetscBool                    :: converged
    PetscInt                     :: istep
    PetscInt                     :: converged_reason
    PetscBool                    :: flg
    PetscErrorCode               :: ierr
    character(len=256), optional :: namelist_filename
    character(len=256)           :: ioerror_msg
    character(len=2560)          :: namelist_buffer
    integer                      :: nml_unit, nml_error
    namelist / problem_options / ncair, ntree

    ncair = 1;
    ntree = 1;

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ncair',ncair,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-ntree',ntree,flg,ierr)

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

    call Init(longwave_mpp)

    dt = 5.0 * 60.d0 ! sec
    istep = 1
    call longwave_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine run_longwave_problem

  !------------------------------------------------------------------------
  subroutine Init(longwave_mpp)
    !
    use longwave_conditions , only : add_conditions_to_goveqns
    use longwave_meshes     , only : setup_meshes
    use longwave_parameters , only : set_parameters
    !
    implicit none
    !
    type(mpp_longwave_type) :: longwave_mpp

    call initialize_mpp(longwave_mpp)

    call setup_meshes(longwave_mpp)

    call add_goveqn(longwave_mpp)

    call add_conditions_to_goveqns(longwave_mpp)

    call longwave_mpp%AllocateAuxVars()

    call longwave_mpp%SetupProblem()

    call set_parameters(longwave_mpp)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(longwave_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_LONGWAVE_KSP
    !
    implicit none

    type(mpp_longwave_type) :: longwave_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call longwave_mpp%Init       ()
    call longwave_mpp%SetName    ('Longwave radiation model')
    call longwave_mpp%SetID      (MPP_LONGWAVE_KSP)
    call longwave_mpp%SetMPIRank (iam)

  end subroutine Initialize_Mpp

  !------------------------------------------------------------------------
  subroutine add_goveqn(longwave_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_LONGWAVE
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_longwave_type) :: longwave_mpp

    LONGWAVE_GE = 1

    call longwave_mpp%AddGovEqnWithMeshRank(GE_LONGWAVE, 'leaf boundary layer', LONGWAVE_MESH)
    
    call longwave_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqn

 !------------------------------------------------------------------------
  subroutine output_regression_longwave_problem(filename_base, num_cells)
    !
    use regression_mod            , only : regression_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnLongwaveType        , only : goveqn_longwave_type
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
    PetscReal, pointer                :: iup(:),idn(:),iabs(:)
    type(regression_type)             :: regression
    class(goveqn_base_type) , pointer :: goveq

    ncells = (nz_cair+1)*ncair

    allocate(iup(ncells))
    allocate(idn(ncells))
    allocate(iabs(ncells))

    call regression%Init(filename_base, num_cells)
    call regression%OpenOutput()
    
    ieqn = 1

    call longwave_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

    select type(goveq)
    class is (goveqn_longwave_type)
       do k = 1,nz_cair+1
          icell = k
          iup(icell) = goveq%aux_vars_in(icell)%Iup
          idn(icell) = goveq%aux_vars_in(icell)%Idn
          iabs(icell) = goveq%aux_vars_in(icell)%Iabs
       end do       

       category = 'general'
       name = 'longwave_up'; call regression%WriteData(name, category, iup )
       name = 'longwave_dn'; call regression%WriteData(name, category, idn )
       name = 'longwave_abs';call regression%WriteData(name, category, iabs)
    end select

    call regression%CloseOutput()
    
    deallocate(iup)
    deallocate(idn)
    deallocate(iabs)

  end subroutine output_regression_longwave_problem

end module longwave_problem
