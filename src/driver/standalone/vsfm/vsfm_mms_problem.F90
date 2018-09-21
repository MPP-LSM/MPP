
module vsfm_mms_problem

#include <petsc/finclude/petsc.h>
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  use vsfm_mms_vars
  use vsfm_mms_steady_state_soil_only_1D
  use MultiPhysicsProbVSFM , only : mpp_vsfm_type
  use mpp_mesh_utils
  
  implicit none

  public :: run_vsfm_mms_problem

  
contains

  !------------------------------------------------------------------------
  
  subroutine run_vsfm_mms_problem()
    !
    implicit none
    !
    PetscBool           :: converged
    PetscInt            :: converged_reason
    PetscBool           :: flg
    PetscErrorCode      :: ierr
    character(len=256)  :: true_soln_fname
    character(len=256)  :: perm_fname
    character(len=256)  :: source_fname
    character(len=256)  :: liq_sat_fname
    type(mpp_vsfm_type) :: vsfm_mpp

    !
    ! Set default values
    !
    call set_default_prob_for_steady_soil_only_1D()
    true_soln_fname = ''
    perm_fname      = ''
    source_fname    = ''

    !
    ! Get command line options
    !
    call PetscOptionsGetInt( PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-problem_type', problem_type, flg, ierr)
    CHKERRQ(ierr)
    
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-nx', nx, flg, ierr)

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-ny', ny, flg, ierr)

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-nz', nz, flg, ierr)

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_true_solution', true_soln_fname, flg, ierr)
    
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_permeability', perm_fname, flg, ierr)
    
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_source', source_fname, flg, ierr)
    
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_liq_saturation', liq_sat_fname, flg, ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         "-fully_saturated", fully_saturated, flg, ierr)

    !
    ! Set up problem after command line options
    !
    select case(problem_type)
    case (STEADY_STATE_SOIL_ONLY_1D)
       call setup_prob_for_steady_soil_only_1D()
          
    case default
       write(iulog,*)'Unknown -problem_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    ! Initialize the problem
    call Init(vsfm_mpp)

    ! Preform Pre-StepDT operations
    call vsfm_mpp%soe%PreStepDT()

    call set_source_sink_conditions(vsfm_mpp)
    call set_boundary_conditions(vsfm_mpp)
    call vsfm_mpp%soe%StepDT(1.d0, 1, &
         converged, converged_reason, ierr); CHKERRQ(ierr)

    if (len(trim(adjustl(true_soln_fname)))>0) then
       call save_problem_variable(vsfm_mpp, true_soln_fname, DATA_PRESSURE)
    end if

    if (len(trim(adjustl(perm_fname)))>0) then
       call save_problem_variable(vsfm_mpp, perm_fname, DATA_PERMEABILITY)
    end if

    if (len(trim(adjustl(source_fname)))>0) then
       call save_problem_variable(vsfm_mpp, source_fname, DATA_MASS_SOURCE)
    end if

    if (len(trim(adjustl(liq_sat_fname)))>0) then
       call save_problem_variable(vsfm_mpp, liq_sat_fname, DATA_LIQUID_SATURATION)
    end if

  end subroutine run_vsfm_mms_problem

  !------------------------------------------------------------------------
  subroutine Init(vsfm_mpp)
    !
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp(vsfm_mpp)

    ! 2. Add all meshes needed for the MPP
    call add_meshes(vsfm_mpp)

    ! 3. Add all governing equations
    call add_goveqns(vsfm_mpp)

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns(vsfm_mpp)

    ! 5. Allocate memory to hold auxvars
    call vsfm_mpp%AllocateAuxVars()

    ! 6. Setup the MPP
    call vsfm_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties(vsfm_mpp)

    ! 8. Set initial conditions
    call set_initial_conditions(vsfm_mpp)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp(vsfm_mpp)
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
#include <petsc/finclude/petsc.h>
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscInt            :: iam
    PetscErrorCode      :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call vsfm_mpp%Init       ()
    call vsfm_mpp%SetName    ('VSFM for MMS')
    call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes(vsfm_mpp)
    !
    use MultiPhysicsProbConstants , only : MESH_ALONG_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_AGAINST_GRAVITY
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : VAR_XC
    use MultiPhysicsProbConstants , only : VAR_YC
    use MultiPhysicsProbConstants , only : VAR_ZC
    use MultiPhysicsProbConstants , only : VAR_DX
    use MultiPhysicsProbConstants , only : VAR_DY
    use MultiPhysicsProbConstants , only : VAR_DZ
    use MultiPhysicsProbConstants , only : VAR_AREA
    use MultiPhysicsProbConstants , only : CONN_SET_INTERNAL
    use MultiPhysicsProbConstants , only : CONN_IN_XYZ_DIR
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscInt            :: imesh
    PetscInt  , pointer :: soil_filter(:)
    PetscReal , pointer :: soil_xc_1d(:), soil_yc_1d(:), soil_zc_1d(:)
    PetscReal , pointer :: soil_dx_1d(:), soil_dy_1d(:), soil_dz_1d(:)
    PetscReal , pointer :: soil_area(:)
    PetscInt            :: nconn
    PetscInt  , pointer :: conn_id_up(:)
    PetscInt  , pointer :: conn_id_dn(:)
    PetscReal , pointer :: conn_dist_up(:)
    PetscReal , pointer :: conn_dist_dn(:)
    PetscReal , pointer :: conn_area(:)
    PetscInt  , pointer :: conn_type(:)
    PetscErrorCode      :: ierr

    call ComputeXYZCentroids(nx, ny, nz, dx, dy, dz, &
         x_min, y_min, z_min, &
         soil_xc_1d, soil_yc_1d, soil_zc_1d)

    allocate(soil_filter(ncells))
    allocate(soil_dx_1d (ncells))
    allocate(soil_dy_1d (ncells))
    allocate(soil_dz_1d (ncells))
    allocate(soil_area  (ncells))

    soil_filter (:) = 1
    soil_dx_1d  (:) = dx
    soil_dy_1d  (:) = dy
    soil_dz_1d  (:) = dz
    soil_area   (:) = dx*dy
    
    call vsfm_mpp%SetNumMeshes(1)

    imesh = 1
    call vsfm_mpp%MeshSetName        (imesh, 'Soil mesh')
    call vsfm_mpp%MeshSetOrientation (imesh, MESH_AGAINST_GRAVITY)
    call vsfm_mpp%MeshSetID          (imesh, MESH_CLM_SOIL_COL)
    call vsfm_mpp%MeshSetDimensions  (imesh, ncells, ncells_ghost, nz)

    call vsfm_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc_1d)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc_1d)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc_1d)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx_1d)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy_1d)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz_1d)
    call vsfm_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
    call vsfm_mpp%MeshComputeVolume          (imesh)

    call ComputeInternalConnections(nx, ny, nz, dx, dy, dz, CONN_IN_XYZ_DIR, &
         nconn, conn_id_up, conn_id_dn, conn_dist_up, conn_dist_dn, &
         conn_area, conn_type)

    call vsfm_mpp%CreateAndAddConnectionSet(imesh, CONN_SET_INTERNAL, &
         nconn,  conn_id_up, conn_id_dn,          &
         conn_dist_up, conn_dist_dn,  conn_area,  &
         conn_type)

    deallocate(soil_filter  )
    deallocate(soil_dx_1d   )
    deallocate(soil_dy_1d   )
    deallocate(soil_dz_1d   )
    deallocate(soil_area    )
    deallocate(conn_id_up   )
    deallocate(conn_id_dn   )
    deallocate(conn_dist_up )
    deallocate(conn_dist_dn )
    deallocate(conn_area    )
    deallocate(conn_type    )

  end subroutine add_meshes

  !------------------------------------------------------------------------

  subroutine add_goveqns(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : GE_RE
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp

    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE', MESH_CLM_SOIL_COL)

    call vsfm_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_SS
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_MASS_RATE
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscInt                            :: ieqn, ii, jj, kk
    PetscInt  :: count, nconn
    PetscInt                  , pointer :: id_up    (:)
    PetscInt                  , pointer :: id_dn    (:)
    PetscReal                 , pointer :: dist_up  (:)
    PetscReal                 , pointer :: dist_dn  (:)
    PetscReal                 , pointer :: area     (:)
    PetscInt                  , pointer :: itype    (:)
    PetscReal                 , pointer :: unit_vec (:,:)
    class(connection_set_type), pointer :: conn_set

    ieqn = 1
    
    call ComputeBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    call MeshCreateConnectionSet(vsfm_mpp%meshes(1), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call vsfm_mpp%soe%AddConditionInGovEqn( ieqn , &
         ss_or_bc_type = COND_BC                 , &
         name          = 'Pressure BC'           , &
         unit          = 'Pa'                    , &
         cond_type     = COND_DIRICHLET          , &
         conn_set      = conn_set)

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn  , &
         ss_or_bc_type = COND_SS                 , &
         name          = 'Source term for MMS'   , &
         unit          = 'kg/m^3'                , &
         cond_type     = COND_MASS_RATE          , &
         region_type   = ALL_CELLS)

    !call conn_set%Destroy()

    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

  end subroutine add_conditions_to_goveqns
  
  !------------------------------------------------------------------------
  subroutine set_material_properties(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM , only : VSFMMPPSetDensityType
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSoilPorosity
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSaturationFunction
    use MultiPhysicsProbVSFM , only : VSFMMPPSetSoilPermeability
    use MultiPhysicsProbVSFM , only : VSFMMPPSetRelativePermeability
    use SaturationFunction   , only : SAT_FUNC_VAN_GENUCHTEN
    use EOSWaterMod          , only : DENSITY_TGDPB01
    use vsfm_mms_vars
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscInt           :: eqn_id
    PetscReal, pointer :: por(:)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: lambda(:)
    PetscReal, pointer :: residual_sat(:)
    PetscReal, pointer :: perm(:)
    PetscInt , pointer :: satfunc_type(:)

    allocate(por             (ncells))
    allocate(perm            (ncells))
    allocate(alpha           (ncells))
    allocate(lambda          (ncells))
    allocate(residual_sat    (ncells))
    allocate(satfunc_type    (ncells))
    
    eqn_id = 1

    call VSFMMPPSetDensityType(vsfm_mpp, eqn_id, DENSITY_TGDPB01)
    
    satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN
    call set_variable_for_problem(problem_type , DATA_POROSITY      , por          )
    call set_variable_for_problem(problem_type , DATA_PERMEABILITY  , perm         )
    call set_variable_for_problem(problem_type , DATA_SATFUNC_ALPHA , alpha        )
    call set_variable_for_problem(problem_type , DATA_SATFUNC_LAMBDA, lambda       )
    call set_variable_for_problem(problem_type , DATA_RES_SAT       , residual_sat )

    call VSFMMPPSetSoilPorosity(vsfm_mpp, eqn_id, por)

    call VSFMMPPSetSoilPermeability(vsfm_mpp, eqn_id, perm, perm, perm)

    call VSFMMPPSetSaturationFunction(vsfm_mpp, eqn_id, satfunc_type, &
         alpha, lambda, residual_sat)

    deallocate(por             )
    deallocate(alpha           )
    deallocate(lambda          )
    deallocate(residual_sat    )
    deallocate(satfunc_type    )

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscReal, pointer :: pressure_ic(:)

    allocate(pressure_ic(ncells))

    call set_variable_for_problem(problem_type, DATA_INITIAL_PRESSURE, data_1D=pressure_ic)

    call vsfm_mpp%Restart(pressure_ic)
    
    deallocate(pressure_ic)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_boundary_conditions(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsVSFMType           , only : sysofeqns_vsfm_type
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use MultiPhysicsProbConstants           , only : AUXVAR_BC
    use MultiPhysicsProbConstants           , only : VAR_BC_SS_CONDITION
    use GoverningEquationBaseType           , only : goveqn_base_type
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscReal                  , pointer :: pressure_bc(:)
    class(sysofeqns_base_type) , pointer :: base_soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: nconn
    PetscInt                             :: soe_auxvar_id

    nconn = 0

    if (nx > 1) nconn = nconn + ny*nz*2
    if (ny > 1) nconn = nconn + nx*nz*2
    if (nz > 1) nconn = nconn + nx*ny*2

    allocate (pressure_bc   (nconn ))

    call set_variable_for_problem(prob_type=problem_type, &
         data_type=DATA_PRESSURE_BC, data_1D=pressure_bc)

    base_soe => vsfm_mpp%soe

    select type(base_soe)
    class is (sysofeqns_vsfm_type)

       soe_auxvar_id = 1
       
       call base_soe%SetDataFromCLM(AUXVAR_BC, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, pressure_bc)

    end select

    deallocate (pressure_bc  )

  end subroutine set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_source_sink_conditions(vsfm_mpp)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsVSFMType , only : sysofeqns_vsfm_type
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    !
    PetscInt                             :: soe_auxvar_id
    PetscReal                  , pointer :: source_sink(:)
    class(sysofeqns_base_type) , pointer :: base_soe
    PetscErrorCode                       :: ierr

    allocate(source_sink(ncells))

    call set_variable_for_problem(problem_type, DATA_MASS_SOURCE, data_1D=source_sink)

    base_soe => vsfm_mpp%soe

    select type(base_soe)
    class is (sysofeqns_vsfm_type)

       soe_auxvar_id = 1
       
       call base_soe%SetDataFromCLM(AUXVAR_SS, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, source_sink)
    end select

    deallocate(source_sink)
    
  end subroutine set_source_sink_conditions

  !------------------------------------------------------------------------
  subroutine set_variable_for_problem(prob_type, data_type, data_1D)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    PetscInt           :: prob_type
    PetscInt           :: data_type
    PetscReal, pointer :: data_1D(:)
    !
    PetscInt           :: ii, jj, kk
    PetscInt           :: count
    PetscReal          :: x, y, z
    !-----------------------------------------------------------------------

    select case(prob_type)
    case (STEADY_STATE_SOIL_ONLY_1D)
       call set_variable_for_steady_state_1D(data_type, data_1D)
       
    case default
       write(iulog,*)'Unsupported problem type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine set_variable_for_problem

  !------------------------------------------------------------------------
  subroutine save_problem_variable(vsfm_mpp, filename, variable_id)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants    , only : AUXVAR_SS
    use MultiPhysicsProbConstants    , only : VAR_BC_SS_CONDITION
    use SystemOfEquationsBaseType    , only : sysofeqns_base_type
    use SystemOfEquationsThermalType , only : sysofeqns_thermal_type
    use petscvec
    !
    implicit none
    !
    type(mpp_vsfm_type) :: vsfm_mpp
    character(len=256)                   :: filename
    integer                              :: variable_id
    !
    PetscInt                             :: soe_auxvar_id
    PetscReal                  , pointer :: val(:)
    Vec                                  :: variable
    PetscViewer                          :: viewer
    class(sysofeqns_base_type) , pointer :: base_soe
    character(len=256)                   :: string
    PetscErrorCode                       :: ierr

    base_soe => vsfm_mpp%soe

    string = trim(filename)

    call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)

    call VecDuplicate(vsfm_mpp%soe%solver%soln, variable, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(variable, val, ierr); CHKERRQ(ierr)

    call set_variable_for_problem(problem_type, variable_id, data_1D=val)

    call VecRestoreArrayF90(variable, val, ierr); CHKERRQ(ierr)

    call VecView(variable,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    call VecDestroy(variable, ierr); CHKERRQ(ierr)

  end subroutine save_problem_variable

end module vsfm_mms_problem
