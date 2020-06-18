module mlc_problem

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use mlc_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none

  type(mpp_mlc_type) :: mlc_mpp

#include <petsc/finclude/petsc.h>

  public :: run_mlc_problem
  public :: output_regression_mlc_problem

contains

  !------------------------------------------------------------------------
  subroutine run_mlc_problem(namelist_filename)

    use MultiPhysicsProbMLC, only : mpp_mlc_type
    !
    implicit none
    !
    PetscReal          :: dt
    PetscBool          :: converged
    PetscInt           :: istep
    PetscInt           :: converged_reason
    PetscBool          :: flg
    PetscErrorCode     :: ierr
    character(len=256), optional :: namelist_filename
    character(len=256) :: ioerror_msg
    character(len=2560):: namelist_buffer
    integer            :: nml_unit, nml_error
    namelist / problem_options / ncair

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

    call Init(mlc_mpp)
    !call mlc_mpp%soe%PrintInfo()

    dt = 5.0 * 60.d0 ! sec
    istep = 1
    call mlc_mpp%soe%SetDtime(dt)
    call mlc_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine run_mlc_problem

  !------------------------------------------------------------------------
  subroutine Init(mlc_mpp)

    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use mlc_conditions            , only : add_conditions_to_goveqns
    use mlc_meshes                , only : setup_meshes
    use mlc_parameters            , only : set_parameters

    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp(mlc_mpp)

    ! 2. Add all meshes needed for the MPP
    call setup_meshes(mlc_mpp)
    !call add_meshes(mlc_mpp)
    !call add_connection_sets_within_meshes(mlc_mpp)

    ! 3. Add govering equations
    call add_multiple_goveqns(mlc_mpp)

    ! 5. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns(mlc_mpp)
    !call add_conditions_to_goveqns(mlc_mpp)

    ! 6. Add internal coupling variables
    !call add_internal_coupling_vars(mlc_mpp)

    ! 7. Allocate auxvars
    call mlc_mpp%AllocateAuxVars(ncair)

    ! 8. Setup Problem
    call mlc_mpp%SetupProblem()

    ! 9.
    call set_parameters(mlc_mpp)

    ! 10. 
    call set_initial_conditions(mlc_mpp)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(mlc_mpp)
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_MLC_KSP
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call mlc_mpp%Init       ()
    call mlc_mpp%SetName    ('Multi Layer Canopy Model')
    call mlc_mpp%SetID      (MPP_MLC_KSP)
    call mlc_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp


  !------------------------------------------------------------------------
  subroutine add_multiple_goveqns(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants      , only : GE_CANOPY_AIR_TEMP
    use MultiPhysicsProbConstants      , only : GE_CANOPY_AIR_VAPOR
    use MultiPhysicsProbConstants      , only : GE_CANOPY_LEAF_TEMP
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType , only : goveqn_cair_temp_type
    use GoveqnCanopyAirVaporType       , only : goveqn_cair_vapor_type
    use GoveqnCanopyLeafTemperatureType, only : goveqn_cleaf_temp_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type)                :: mlc_mpp
    class(goveqn_base_type) , pointer :: goveq
    PetscInt                , pointer :: leaf2cair(:)
    PetscInt                          :: nleaf2cair, icair, itree, k, count, ieqn

    CAIR_TEMP_GE = 1
    CAIR_VAPR_GE = 2 
    CLEF_TEMP_SUN_GE = 3
    CLEF_TEMP_SHD_GE = 4

    call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_AIR_TEMP, 'Canopy air temperature', CAIR_MESH)
    
    call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_AIR_VAPOR, 'Canopy air vapor', CAIR_MESH)

    call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_LEAF_TEMP, 'Sunlit canopy', CLEF_MESH)

    call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_LEAF_TEMP, 'Shaded canopy', CLEF_MESH)

    call mlc_mpp%SetMeshesOfGoveqnsByMeshRank()

    if (ntree == 1) then

       do ieqn = 1, 4
          call mlc_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

          select type(goveq)
          class is (goveqn_cair_temp_type)
             call goveq%SetDefaultLeaf2CAirMap()
          class is (goveqn_cair_vapor_type)
             call goveq%SetDefaultLeaf2CAirMap()
          class is (goveqn_cleaf_temp_type)
             call goveq%SetDefaultLeaf2CAirMap()
          end select
       end do
    else

       nleaf2cair = (nz_cair+1)*ncair*ntree

       allocate(leaf2cair(nleaf2cair))

       count = 0
       do icair = 1, ncair
          do itree = 1, ntree
             do k = 1, nz_cair + 1
                count = count + 1
                leaf2cair(count) = (nz_cair + 1)*(icair-1) + k
             end do
          end do
       end do

       do ieqn = 1, 4
          call mlc_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

          select type(goveq)
          class is (goveqn_cair_temp_type)
             call goveq%SetLeaf2CAirMap(leaf2cair, nleaf2cair)
          class is (goveqn_cair_vapor_type)
             call goveq%SetLeaf2CAirMap(leaf2cair, nleaf2cair)
          class is (goveqn_cleaf_temp_type)
             call goveq%SetLeaf2CAirMap(leaf2cair, nleaf2cair)
          end select
       end do

    end if
          
  end subroutine add_multiple_goveqns

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsMLCType  , only : sysofeqns_mlc_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscReal                            :: relhum, eref, esat, desatdt
    PetscInt                             :: p
    PetscInt                             :: nDM
    DM                         , pointer :: dms(:)
    Vec                        , pointer :: soln_subvecs(:)
    PetscReal                  , pointer :: v_p(:)
    PetscInt                             :: ii
    PetscViewer                          :: viewer
    PetscInt                             :: soe_auxvar_id
    PetscErrorCode                       :: ierr

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    p = 1
    soe%cturb%vcan(p) = soe%cturb%vref(p)
    soe%cturb%tcan(p) = soe%cturb%tref(p)

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(mlc_mpp%soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(mlc_mpp%soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(mlc_mpp%soe%solver%dm, &
         mlc_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)

       if (ii == CAIR_TEMP_GE .or. ii == CLEF_TEMP_SUN_GE .or. ii == CLEF_TEMP_SHD_GE) then
          v_p(:) = soe%cturb%tref(1)

       else if (ii == CAIR_VAPR_GE) then
          v_p(:) = soe%cturb%vref(1)          

       endif

       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(mlc_mpp%soe%solver%dm, &
         mlc_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)
    
    call VecCopy(mlc_mpp%soe%solver%soln, mlc_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(mlc_mpp%soe%solver%soln, mlc_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine output_regression_mlc_problem(filename_base, num_cells)
    !
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants       , only : VAR_PRESSURE
    use MultiPhysicsProbConstants       , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants       , only : VAR_WATER_VAPOR
    use MultiPhysicsProbConstants       , only : VAR_LEAF_TEMPERATURE
    use regression_mod                  , only : regression_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType  , only : goveqn_cair_temp_type
    use GoveqnCanopyAirVaporType        , only : goveqn_cair_vapor_type
    use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
    !
    implicit none
    !
    character(len=256)    :: filename_base
    PetscInt, intent(in)  :: num_cells
    !
    PetscInt              :: output, ieqn, ncells
    character(len=512)    :: filename
    character(len=64)     :: category
    character(len=64)     :: name
    PetscReal, pointer    :: data(:)
    type(regression_type) :: regression
    class(goveqn_base_type) , pointer :: goveq

    ncells = (nz_cair+1)*ncair

    allocate(data(ncells))

    call regression%Init(filename_base, num_cells)
    call regression%OpenOutput()
    
    do ieqn = 1,4

       call mlc_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

       select type(goveq)
       class is (goveqn_cair_temp_type)
          name = 'air_temperature'; category = 'temperature'
          call goveq%GetRValues(AUXVAR_INTERNAL, VAR_TEMPERATURE, ncells, data)
          
       class is (goveqn_cair_vapor_type)
          name = 'air_vapor'; category = 'general'
          call goveq%GetRValues(AUXVAR_INTERNAL, VAR_WATER_VAPOR, ncells, data)

       class is (goveqn_cleaf_temp_type)
          if (ieqn == 3) then
             name = 'sunlit_leaf_temperature'
          else
             name = 'shaded_leaf_temperature'
          end if
          category = 'temperature'
          call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_TEMPERATURE, ncells, data)
       end select

       call regression%WriteData(name, category, data)

    enddo

    call regression%CloseOutput()
    
    deallocate(data)

  end subroutine output_regression_mlc_problem

end module mlc_problem

