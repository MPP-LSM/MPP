module shortwave_problem

  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbShortwave , only : mpp_shortwave_type
  use shortwave_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none


#include <petsc/finclude/petsc.h>

  type(mpp_shortwave_type):: shortwave_mpp

  public :: run_shortwave_problem
  public :: output_regression_shortwave_problem

contains

  !------------------------------------------------------------------------
  subroutine run_shortwave_problem(namelist_filename)

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

    call Init(shortwave_mpp)

    dt = 5.0 * 60.d0 ! sec
    istep = 1
    call shortwave_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine run_shortwave_problem

  !------------------------------------------------------------------------
  subroutine Init(shortwave_mpp)
    !
    use shortwave_conditions , only : add_conditions_to_goveqns
    use shortwave_meshes     , only : setup_meshes
    use shortwave_parameters , only : set_parameters
    !
    implicit none
    !
    type(mpp_shortwave_type) :: shortwave_mpp

    call initialize_mpp(shortwave_mpp)

    call setup_meshes(shortwave_mpp)

    call add_goveqn(shortwave_mpp)

    call add_conditions_to_goveqns(shortwave_mpp)

    call shortwave_mpp%AllocateAuxVars()

    call shortwave_mpp%SetupProblem()

    call set_parameters(shortwave_mpp)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(shortwave_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_shortwave_KSP
    !
    implicit none

    type(mpp_shortwave_type) :: shortwave_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call shortwave_mpp%Init       ()
    call shortwave_mpp%SetName    ('Shortgwave radiation model')
    call shortwave_mpp%SetID      (MPP_shortwave_KSP)
    call shortwave_mpp%SetMPIRank (iam)

  end subroutine Initialize_Mpp

  !------------------------------------------------------------------------
  subroutine add_goveqn(shortwave_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_SHORTWAVE
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type) :: shortwave_mpp

    shortwave_GE = 1

    call shortwave_mpp%AddGovEqnWithMeshRank(GE_SHORTWAVE, 'leaf boundary layer', shortwave_MESH)
    
    call shortwave_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqn

 !------------------------------------------------------------------------
  subroutine output_regression_shortwave_problem(filename_base, num_cells)
    !
    use regression_mod            , only : regression_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnSHORTWAVEType        , only : goveqn_shortwave_type
    !
    implicit none
    !
    character(len=256)                :: filename_base
    PetscInt, intent(in)              :: num_cells
    !
    PetscInt                          :: output, ieqn, ncells, icell, k, iband
    character(len=512)                :: filename
    character(len=64)                 :: category
    character(len=64)                 :: name
    PetscReal, pointer                :: iupb_vis(:), iupb_nir(:)
    PetscReal, pointer                :: idnd_vis(:), idnd_nir(:)
    type(regression_type)             :: regression
    class(goveqn_base_type) , pointer :: goveq

    ncells = (nz_cair+1)*ncair

    allocate(iupb_vis(ncells))
    allocate(iupb_nir(ncells))
    allocate(idnd_vis(ncells))
    allocate(idnd_nir(ncells))

    call regression%Init(filename_base, num_cells)
    call regression%OpenOutput()
    
    ieqn = 1

    call shortwave_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

    select type(goveq)
    class is (goveqn_shortwave_type)
       do k = 1,nz_cair+1
          icell = k
          iupb_vis(icell) = goveq%aux_vars_in(icell)%Iup(1)
          iupb_nir(icell) = goveq%aux_vars_in(icell)%Iup(2)

          idnd_vis(icell) = goveq%aux_vars_in(icell)%Idn(1)
          idnd_nir(icell) = goveq%aux_vars_in(icell)%Idn(2)
       end do

       category = 'general'
       name = 'shortwave_up_beam_vis'   ; call regression%WriteData(name, category, iupb_vis )
       name = 'shortwave_up_beam_nir'   ; call regression%WriteData(name, category, iupb_nir )
       name = 'shortwave_dn_diffuse_vis'; call regression%WriteData(name, category, idnd_vis )
       name = 'shortwave_dn_diffuse_nir'; call regression%WriteData(name, category, idnd_nir )
    end select

    call regression%CloseOutput()
    
    deallocate(iupb_vis)
    deallocate(iupb_nir)
    deallocate(idnd_vis)
    deallocate(idnd_nir)

  end subroutine output_regression_shortwave_problem

end module shortwave_problem
