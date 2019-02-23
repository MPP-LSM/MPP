module th_mms_problem

#include <petsc/finclude/petsc.h>

  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbConstants , only : PRESSURE_REF
  use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
  use EOSWaterMod               , only : DENSITY_TGDPB01,DENSITY_CONSTANT,DENSITY_IFC67
  use EOSWaterMod               , only : INT_ENERGY_ENTHALPY_CONSTANT
  use EOSWaterMod               , only : INT_ENERGY_ENTHALPY_IFC67
  use mpp_varcon                , only : denh2o
  use EOSWaterMod               , only : Density, Viscosity
  use MultiPhysicsProbTH        , only : th_mpp
  use mpp_mesh_utils
  use petscsys
  use petscvec
  use petscdm
  use petscdmda

  implicit none

  PetscReal , parameter :: PI = 4 * atan(1.0d0)

  PetscInt              :: problem_type
  
  PetscInt              :: nx,ny,nz
  PetscInt              :: ncells
  PetscInt              :: ncells_ghost

  PetscReal             :: dx, dy, dz
  PetscReal             :: xlim, ylim, zlim
  PetscReal             :: x_min , x_max
  PetscReal             :: y_min , y_max
  PetscReal             :: z_min , z_max
  
  PetscReal , pointer   :: soil_xc_3d(:,:,:)
  PetscReal , pointer   :: soil_yc_3d(:,:,:)
  PetscReal , pointer   :: soil_zc_3d(:,:,:)
  PetscInt  , pointer   :: soil_id_3d(:,:,:)
  
  
  PetscInt  , parameter :: STEADY_STATE_SOIL_ONLY_1D = 1

  PetscInt  , parameter :: DATA_XC                   = 1
  PetscInt  , parameter :: DATA_YC                   = 2
  PetscInt  , parameter :: DATA_ZC                   = 3
  PetscInt  , parameter :: DATA_PERMEABILITY         = 4
  PetscInt  , parameter :: DATA_PRESSURE_BC          = 5
  PetscInt  , parameter :: DATA_MASS_SOURCE          = 6
  PetscInt  , parameter :: DATA_PRESSURE             = 7
  PetscInt  , parameter :: DATA_POROSITY             = 8
  PetscInt  , parameter :: DATA_SATFUNC_TYPE         = 9
  PetscInt  , parameter :: DATA_SATFUNC_ALPHA        = 10
  PetscInt  , parameter :: DATA_SATFUNC_LAMBDA       = 11
  PetscInt  , parameter :: DATA_RES_SAT              = 12
  PetscInt  , parameter :: DATA_INITIAL_PRESSURE     = 13
  PetscInt  , parameter :: DATA_TEMPERATURE          = 14
  PetscInt  , parameter :: DATA_INITIAL_TEMPERATURE  = 15
  PetscInt  , parameter :: DATA_HEAT_SOURCE          = 16
  PetscInt  , parameter :: DATA_HEAT_CAP             = 17
  PetscInt  , parameter :: DATA_T_COND_DRY           = 18
  PetscInt  , parameter :: DATA_T_COND_WET           = 19
  PetscInt  , parameter :: DATA_T_ALPHA              = 20
  PetscInt  , parameter :: DATA_TEMPERATURE_BC       = 21

  PetscInt :: SOIL_MESH
  PetscInt :: SOIL_MASS_GE
  PetscInt :: SOIL_ENTH_GE

  PetscInt :: num_icpl
  PetscInt :: max_Iauxvar_cpl
  PetscInt :: num_bcpl
  PetscInt :: max_Bauxvar_cpl

  PetscInt :: density_type
  PetscInt :: int_energy_enthalpy_type

  PetscInt, pointer :: i_cpl_data(:,:)
  PetscInt, pointer :: b_cpl_data(:,:)

  public :: run_th_mms_problem

contains

  !------------------------------------------------------------------------
  subroutine run_th_mms_problem()

    implicit none

    PetscBool          :: converged
    PetscInt           :: converged_reason
    PetscBool           :: flg
    PetscErrorCode     :: ierr
    PetscReal          :: dtime
    character(len=256) :: true_soln_fname
    character(len=256) :: source_fname
    character(len=256) :: perm_fname

    call set_default_problem()
    true_soln_fname = ''
    source_fname = ''

    !
    ! Get command line options
    !
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-nx', nx, flg, ierr)
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_true_solution', true_soln_fname, flg, ierr)
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_source', source_fname, flg, ierr)
    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-view_permeability', perm_fname, flg, ierr)

    call set_problem()
    
    call Init()

    call set_source_sink_conditions()
    call set_boundary_conditions()
    
    ! Run the model
    call th_mpp%soe%StepDT(1.d0, 1, &
         converged, converged_reason, ierr); CHKERRQ(ierr)

    if (len(trim(adjustl(true_soln_fname)))>0) then
       call save_true_solution(true_soln_fname)
    end if

    if (len(trim(adjustl(source_fname)))>0) then
       call save_source(source_fname)
    end if

    if (len(trim(adjustl(perm_fname)))>0) then
       call save_perm(perm_fname)
    end if

  end subroutine run_th_mms_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    !
    implicit none

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp()

    ! 2. Add all meshes needed for the MPP
    call add_meshes()

    ! 3. Add all governing equations
    call add_goveqns()
    !call add_connection_sets_to_meshes()

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns()
    !call setup_goveqn_connectivity()

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars()

    ! 6. Setup the MPP
    call th_mpp%SetupProblem()

    ! 7. Add material properities associated with all governing equations
    call set_material_properties()

    ! 8. Set initial conditions
    call set_initial_conditions()

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine initialize_mpp()
    !
    ! !DESCRIPTION:
    ! Initialization VSFM
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_TH_SNES_CLM
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call th_mpp%Init       ()
    call th_mpp%SetName    ('Thermal-Hydrology For SPAC')
    call th_mpp%SetID      (MPP_TH_SNES_CLM)
    call th_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes()
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
    
    call th_mpp%SetNumMeshes(1)

    SOIL_MESH = 1
    imesh = SOIL_MESH
    call th_mpp%MeshSetName        (imesh, 'Soil mesh')
    call th_mpp%MeshSetOrientation (imesh, MESH_AGAINST_GRAVITY)
    call th_mpp%MeshSetID          (imesh, MESH_CLM_SOIL_COL)
    call th_mpp%MeshSetDimensions  (imesh, ncells, ncells_ghost, nz)

    call th_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
    call th_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc_1d)
    call th_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc_1d)
    call th_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc_1d)
    call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx_1d)
    call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy_1d)
    call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz_1d)
    call th_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
    call th_mpp%MeshComputeVolume          (imesh)

    call ComputeInternalConnections(nx, ny, nz, dx, dy, dz, CONN_IN_XYZ_DIR, &
         nconn, conn_id_up, conn_id_dn, conn_dist_up, conn_dist_dn, &
         conn_area, conn_type)

    call th_mpp%CreateAndAddConnectionSet(imesh, CONN_SET_INTERNAL, &
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
  subroutine add_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_RE
    use MultiPhysicsProbConstants, only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants, only : GE_THERM_SOIL_EBASED
    use MultiPhysicsProbConstants, only : MESH_CLM_THERMAL_SOIL_COL
    !
    !
    ! !ARGUMENTS
    implicit none

    SOIL_MASS_GE = 1
    call th_mpp%AddGovEqn(GE_RE, 'Mass Equation ODE for Soil', &
         MESH_CLM_SOIL_COL)

    SOIL_ENTH_GE = 2
    call th_mpp%AddGovEqn(GE_THERM_SOIL_EBASED, &
         'Enthalpy-based ODE for heat transport', MESH_CLM_SOIL_COL)

    call th_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
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
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn, imesh
    PetscInt  :: count, nconn
    PetscInt                  , pointer :: id_up    (:)
    PetscInt                  , pointer :: id_dn    (:)
    PetscReal                 , pointer :: dist_up  (:)
    PetscReal                 , pointer :: dist_dn  (:)
    PetscReal                 , pointer :: area     (:)
    PetscInt                  , pointer :: itype    (:)
    PetscReal                 , pointer :: unit_vec (:,:)
    class(connection_set_type), pointer :: conn_set

    call ComputeBoundaryDomainConnection(nx, ny, nz, dx, dy, dz, &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec)

    ieqn = SOIL_MASS_GE
    imesh = SOIL_MESH

    call MeshCreateConnectionSet(th_mpp%meshes(imesh), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call th_mpp%soe%AddConditionInGovEqn( ieqn , &
         ss_or_bc_type = COND_BC                 , &
         name          = 'Pressure BC'           , &
         unit          = 'Pa'                    , &
         cond_type     = COND_DIRICHLET          , &
         conn_set      = conn_set)

    call th_mpp%soe%AddConditionInGovEqn(ieqn  , &
         ss_or_bc_type = COND_SS                 , &
         name          = 'Source term for MMS'   , &
         unit          = 'kg/m^3'                , &
         cond_type     = COND_MASS_RATE          , &
         region_type   = ALL_CELLS)

    nullify(conn_set)

    ieqn  = SOIL_ENTH_GE
    imesh = SOIL_MESH
    
    call MeshCreateConnectionSet(th_mpp%meshes(imesh), &
         nconn, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call th_mpp%soe%AddConditionInGovEqn( ieqn , &
         ss_or_bc_type = COND_BC                 , &
         name          = 'Temperature BC'        , &
         unit          = 'K'                     , &
         cond_type     = COND_DIRICHLET          , &
         conn_set      = conn_set)

    call th_mpp%soe%AddConditionInGovEqn(ieqn  , &
         ss_or_bc_type = COND_SS                 , &
         name          = 'Source term for MMS'   , &
         unit          = 'W/m^3'                , &
         cond_type     = COND_HEAT_RATE          , &
         region_type   = ALL_CELLS)

    nullify(conn_set)

    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine setup_goveqn_connectivity()
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    PetscInt :: ii

    num_icpl        = 1
    max_Iauxvar_cpl = 2
    num_bcpl        = 0
    max_Bauxvar_cpl = 0

    allocate(i_cpl_data(num_icpl, max_Iauxvar_cpl))
    nullify (b_cpl_data)

    ! 'Soil mass PDE' is coupled with 'Soil energy PDE' via internal auxvars
    ii = 1; i_cpl_data(ii,1) = SOIL_MASS_GE; i_cpl_data(ii,2) = SOIL_ENTH_GE

  end subroutine setup_goveqn_connectivity

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    integer           :: ieqn
    integer           :: nvars_for_coupling
    integer, pointer  :: var_ids_for_coupling(:)
    integer, pointer  :: goveqn_ids_for_coupling(:)
    integer, pointer  :: is_bc(:)

    !
    ! Set variables to be exchanged among the coupled equations
    !

    ! mass equation (id=1) needs temperature from energy equation (id=2)
    ieqn = 1
    nvars_for_coupling = 1
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))
    allocate (is_bc                   (nvars_for_coupling))

    var_ids_for_coupling    (1) = VAR_TEMPERATURE
    goveqn_ids_for_coupling (1) = 2
    is_bc                   (1) = 0

    call th_mpp%GovEqnSetInternalCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

    ! energy equation (id=2) needs pressure from energy equation (id=1)
    ieqn = 2
    nvars_for_coupling = 1
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))

    var_ids_for_coupling    (1) = VAR_PRESSURE
    goveqn_ids_for_coupling (1) = 1
    is_bc                   (1) = 0

    call th_mpp%GovEqnSetInternalCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling)

    deallocate(var_ids_for_coupling   )
    deallocate(goveqn_ids_for_coupling)

    !
    ! Allocate auxvars
    !
    call th_mpp%AllocateAuxVars()

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine ComputeCouplingGEInfo(ieqn, num_icpl, max_Iauxvar_cpl, &
       num_bcpl, max_Bauxvar_cpl, i_cpl_data, b_cpl_data, nvars_for_coupling, &
       var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)
    !
    use MultiPhysicsProbTH            , only : th_mpp
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use MultiPhysicsProbConstants     , only : VAR_PRESSURE
    use MultiPhysicsProbConstants     , only : VAR_TEMPERATURE
    use GoverningEquationBaseType
    use GoveqnRichardsODEPressureType
    use RichardsODEPressureAuxType
    !
    implicit none
    !
    integer                         :: ieqn
    PetscInt                        :: num_icpl
    PetscInt                        :: max_Iauxvar_cpl
    PetscInt                        :: num_bcpl
    PetscInt                        :: max_Bauxvar_cpl
    PetscInt , pointer              :: i_cpl_data(:,:)
    PetscInt , pointer              :: b_cpl_data(:,:)
    integer                         :: nvars_for_coupling
    integer  , pointer              :: var_ids_for_coupling(:)
    integer  , pointer              :: goveqn_ids_for_coupling(:)
    integer  , pointer              :: is_bc(:)
    !
    PetscInt                        :: ii, jj
    PetscBool                       :: ieqn_found
    class(goveqn_base_type),pointer :: cur_goveq

    nvars_for_coupling = 0

    do ii = 1, num_icpl
       ieqn_found = PETSC_FALSE
       do jj = 1, max_Iauxvar_cpl
          if (i_cpl_data(ii,jj) == ieqn) ieqn_found = PETSC_TRUE
       end do
       if (ieqn_found) then
          do jj = 1, max_Iauxvar_cpl
             if (i_cpl_data(ii,jj) /= ieqn) nvars_for_coupling = nvars_for_coupling + 1
          end do
       end if
    end do

    do ii = 1, num_bcpl
       if (b_cpl_data(ii,1) == ieqn) &
            nvars_for_coupling = nvars_for_coupling + b_cpl_data(ii,4)
    end do

    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))
    allocate (is_bc                   (nvars_for_coupling))

    nvars_for_coupling = 0

    do ii = 1, num_icpl
       ieqn_found = PETSC_FALSE
       do jj = 1, max_Iauxvar_cpl
          if (i_cpl_data(ii,jj) == ieqn) ieqn_found = PETSC_TRUE
       end do
       if (ieqn_found) then
          do jj = 1, max_Iauxvar_cpl
             if (i_cpl_data(ii,jj) /= ieqn) then
                nvars_for_coupling = nvars_for_coupling + 1
                goveqn_ids_for_coupling(nvars_for_coupling) = i_cpl_data(ii,jj)
                is_bc                  (nvars_for_coupling) = 0
             end if
          end do
       end if
    end do

    do ii = 1, num_bcpl
       if (b_cpl_data(ii,1) == ieqn) then
          do jj = 1, b_cpl_data(ii,4)
             nvars_for_coupling = nvars_for_coupling + 1
             goveqn_ids_for_coupling(nvars_for_coupling) = b_cpl_data(ii,4+jj)
             is_bc                  (nvars_for_coupling) = 1
          end do
       end if
    end do

    do ii = 1, nvars_for_coupling

       cur_goveq => th_mpp%soe%goveqns
       do jj = 1, goveqn_ids_for_coupling(ii) - 1
          cur_goveq => cur_goveq%next
       end do

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)          
          var_ids_for_coupling(ii) = VAR_PRESSURE

       class is (goveqn_thermal_enthalpy_soil_type)
          var_ids_for_coupling(ii) = VAR_TEMPERATURE

       class default
          write(*,*)'SPACMPPSetupPetscSNESSetup: Unknown goveq type.'
          stop
       end select

    enddo
    
  end subroutine ComputeCouplingGEInfo

  !------------------------------------------------------------------------
  subroutine set_material_properties()
    !
    ! !DESCRIPTION:
    !
    use GoverningEquationBaseType    , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType, only : goveqn_thermal_enthalpy_soil_type
    use GoveqnRichardsODEPressureType, only : goveqn_richards_ode_pressure_type
    !
    implicit none
    !
    class(goveqn_base_type), pointer :: cur_goveq
    !-----------------------------------------------------------------------

    cur_goveq => th_mpp%soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_richards_ode_pressure_type)
          call set_material_properties_for_mass_equation(cur_goveq)

       class is (goveqn_thermal_enthalpy_soil_type)
          call set_material_properties_for_energy_equation(cur_goveq)

       class default
          write(*,*)'SPACMPPSetupPetscSNESSetup: Unknown goveq type.'
          stop
       end select

       cur_goveq => cur_goveq%next
    enddo

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_mass_equation(goveq_richards)
    !
    ! !DESCRIPTION:
    !
    use GoveqnRichardsODEPressureType
    use SaturationFunction   , only : SAT_FUNC_VAN_GENUCHTEN
    !
    implicit none

    class(goveqn_richards_ode_pressure_type),pointer    :: goveq_richards

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
    
    call goveq_richards%SetDensityType(density_type)

    satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN
    call set_variable_for_problem(DATA_POROSITY      , por          )
    call set_variable_for_problem(DATA_PERMEABILITY  , perm         )
    call set_variable_for_problem(DATA_SATFUNC_ALPHA , alpha        )
    call set_variable_for_problem(DATA_SATFUNC_LAMBDA, lambda       )
    call set_variable_for_problem(DATA_RES_SAT       , residual_sat )

    call goveq_richards%SetSoilPermeability(perm, perm, perm)
    call goveq_richards%SetSoilPorosity(por)
    call goveq_richards%SetSaturationFunction( &
           satfunc_type, alpha, lambda, residual_sat)

    deallocate(por             )
    deallocate(alpha           )
    deallocate(lambda          )
    deallocate(residual_sat    )
    deallocate(satfunc_type    )

  end subroutine set_material_properties_for_mass_equation

  !------------------------------------------------------------------------
  subroutine set_material_properties_for_energy_equation(goveq_enthalpy)
    !
    ! !DESCRIPTION:
    !
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use SaturationFunction            , only : SAT_FUNC_VAN_GENUCHTEN
    !
    implicit none
    !
    class(goveqn_thermal_enthalpy_soil_type) , pointer   :: goveq_enthalpy

    PetscReal, pointer :: por(:)
    PetscReal, pointer :: alpha(:)
    PetscReal, pointer :: lambda(:)
    PetscReal, pointer :: residual_sat(:)
    PetscReal, pointer :: perm(:)
    PetscInt , pointer :: satfunc_type(:)
    PetscReal, pointer :: heat_cap(:)
    PetscReal, pointer :: t_cond_dry(:)
    PetscReal, pointer :: t_cond_wet(:)
    PetscReal, pointer :: t_alpha(:)

    allocate(por          (ncells))
    allocate(perm         (ncells))
    allocate(alpha        (ncells))
    allocate(lambda       (ncells))
    allocate(residual_sat (ncells))
    allocate(satfunc_type (ncells))
    allocate(heat_cap     (ncells))
    allocate(t_cond_dry   (ncells))
    allocate(t_cond_wet   (ncells))
    allocate(t_alpha      (ncells))
    
    call goveq_enthalpy%SetDensityType(density_type)
    call goveq_enthalpy%SetIntEnergyEnthalpyType(int_energy_enthalpy_type)

    satfunc_type(:) = SAT_FUNC_VAN_GENUCHTEN
    call set_variable_for_problem(DATA_POROSITY      , por          )
    call set_variable_for_problem(DATA_PERMEABILITY  , perm         )
    call set_variable_for_problem(DATA_SATFUNC_ALPHA , alpha        )
    call set_variable_for_problem(DATA_SATFUNC_LAMBDA, lambda       )
    call set_variable_for_problem(DATA_RES_SAT       , residual_sat )
    call set_variable_for_problem(DATA_HEAT_CAP      , heat_cap     )
    call set_variable_for_problem(DATA_T_COND_DRY    , t_cond_dry )
    call set_variable_for_problem(DATA_T_COND_WET    , t_cond_wet )
    call set_variable_for_problem(DATA_T_ALPHA       , t_alpha )

    call goveq_enthalpy%SetSoilPermeability(perm, perm, perm)
    call goveq_enthalpy%SetSoilPorosity(por)
    call goveq_enthalpy%SetSaturationFunction( &
           satfunc_type, alpha, lambda, residual_sat)

    call goveq_enthalpy%SetHeatCapacity(heat_cap)
    call goveq_enthalpy%SetThermalCondDry(t_cond_dry)
    call goveq_enthalpy%SetThermalCondWet(t_cond_wet)
    call goveq_enthalpy%SetThermalAlpha(t_alpha)

    deallocate(por          )
    deallocate(perm         )
    deallocate(alpha        )
    deallocate(lambda       )
    deallocate(residual_sat )
    deallocate(satfunc_type )
    deallocate(heat_cap     )
    deallocate(t_cond_dry   )
    deallocate(t_cond_wet   )
    deallocate(t_alpha      )

  end subroutine set_material_properties_for_energy_equation

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbTH, only : th_mpp
    !
    implicit none
    !
    PetscErrorCode      :: ierr
    PetscInt            :: nDM
    DM        , pointer :: dms(:)
    Vec       , pointer :: soln_subvecs(:)
    PetscReal , pointer :: v_p(:)
    PetscInt            :: ii
    PetscViewer         :: viewer
    PetscInt            :: soe_auxvar_id

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(th_mpp%soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(th_mpp%soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(th_mpp%soe%solver%dm, &
         th_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)
       if (ii == 1) then
          call set_variable_for_problem(DATA_INITIAL_PRESSURE, v_p)
       else
          call set_variable_for_problem(DATA_INITIAL_TEMPERATURE, v_p)
       end if
       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(th_mpp%soe%solver%dm, &
         th_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)
    
    call VecCopy(th_mpp%soe%solver%soln, th_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(th_mpp%soe%solver%soln, th_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_boundary_conditions()
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    use SystemOfEquationsTHType       , only : sysofeqns_th_type
    use SystemOfEquationsBaseType     , only : sysofeqns_base_type
    use MultiPhysicsProbConstants     , only : AUXVAR_BC
    use MultiPhysicsProbConstants     , only : VAR_BC_SS_CONDITION
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnRichardsODEPressureType , only : goveqn_richards_ode_pressure_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    implicit none
    !
    PetscInt                             :: ieqn
    PetscInt                             :: iconn
    PetscInt                             :: sum_conn
    PetscReal                  , pointer :: pres_bc(:)
    PetscReal                  , pointer :: temp_bc(:)
    type(condition_type)       , pointer :: cur_cond
    class(connection_set_type) , pointer :: cur_conn_set
    class(goveqn_base_type)    , pointer :: cur_goveq

    allocate (pres_bc(2))
    allocate (temp_bc(2))

    call set_variable_for_problem(DATA_PRESSURE_BC, pres_bc)
    call set_variable_for_problem(DATA_TEMPERATURE_BC, temp_bc)

    ieqn = SOIL_MASS_GE
    call th_mpp%soe%SetDataFromCLM(AUXVAR_BC, VAR_BC_SS_CONDITION, ieqn, pres_bc)

    ieqn = SOIL_ENTH_GE
    call th_mpp%soe%SetDataFromCLM(AUXVAR_BC, VAR_BC_SS_CONDITION, ieqn, temp_bc)

    cur_goveq => th_mpp%soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_enthalpy_soil_type)
          sum_conn = 0

          cur_cond => cur_goveq%boundary_conditions%first
          do
             if (.not.associated(cur_cond)) exit
             cur_conn_set => cur_cond%conn_set

             do iconn = 1, cur_conn_set%num_connections
                sum_conn = sum_conn + 1
                cur_goveq%aux_vars_bc(sum_conn)%pressure = pres_bc(iconn)

             enddo
             cur_cond => cur_cond%next
          enddo

       class is (goveqn_richards_ode_pressure_type)
          sum_conn = 0

          cur_cond => cur_goveq%boundary_conditions%first
          do
             if (.not.associated(cur_cond)) exit
             cur_conn_set => cur_cond%conn_set

             do iconn = 1, cur_conn_set%num_connections
                sum_conn = sum_conn + 1
                cur_goveq%aux_vars_bc(sum_conn)%temperature = temp_bc(iconn)

             enddo
             cur_cond => cur_cond%next
          enddo

       class default

          write(*,*)'Unknown goveqn type'
          stop

       end select
       cur_goveq => cur_goveq%next
    enddo
    
    deallocate (pres_bc)
    deallocate (temp_bc)

  end subroutine set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_source_sink_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsTHType   , only : sysofeqns_th_type
    !
    PetscInt                             :: ieqn
    PetscReal                  , pointer :: source_sink(:)
    class(sysofeqns_base_type) , pointer :: base_soe
    PetscErrorCode                       :: ierr

    allocate(source_sink(ncells))

    base_soe => th_mpp%soe

    select type(base_soe)
    class is (sysofeqns_th_type)

       ieqn = SOIL_MASS_GE
       call set_variable_for_problem(DATA_MASS_SOURCE, source_sink)
       call base_soe%SetDataFromCLM(AUXVAR_SS, VAR_BC_SS_CONDITION, ieqn, source_sink)

       ieqn = SOIL_ENTH_GE
       call set_variable_for_problem(DATA_HEAT_SOURCE, source_sink)
       call base_soe%SetDataFromCLM(AUXVAR_SS, VAR_BC_SS_CONDITION, ieqn, source_sink)

    end select

    deallocate(source_sink)

  end subroutine set_source_sink_conditions

  !------------------------------------------------------------------------
  subroutine set_default_problem()

    implicit none

    problem_type = STEADY_STATE_SOIL_ONLY_1D

    nx = 20
    ny = 1
    nz = 1

    x_min = 0.d0; x_max = 10.d0
    y_min = 0.d0; y_max = 1.d0
    z_min = 0.d0; z_max = 1.d0

    xlim = x_max - x_min
    ylim = y_max - y_min
    zlim = z_max - z_min

  end subroutine set_default_problem

  !------------------------------------------------------------------------
  subroutine set_problem()

    implicit none

    problem_type = STEADY_STATE_SOIL_ONLY_1D
    density_type = DENSITY_IFC67
    density_type = DENSITY_CONSTANT
    int_energy_enthalpy_type = INT_ENERGY_ENTHALPY_IFC67
    !int_energy_enthalpy_type = INT_ENERGY_ENTHALPY_CONSTANT

    xlim = x_max - x_min
    ylim = y_max - y_min
    zlim = z_max - z_min

    ncells = nx*ny*nz
    ncells_ghost = 0

    dx = (x_max - x_min)/nx
    dy = (y_max - y_min)/ny
    dz = (z_max - z_min)/nz

    allocate(soil_xc_3d(nx,ny,nz))
    allocate(soil_yc_3d(nx,ny,nz))
    allocate(soil_zc_3d(nx,ny,nz))
    allocate(soil_id_3d(nx,ny,nz))
    
    call ComputeXYZCentroids(nx, ny, nz, dx, dy, dz, &
         x_min, y_min, z_min, &
         soil_xc_3d, soil_yc_3d, soil_zc_3d)

  end subroutine set_problem

  !------------------------------------------------------------------------
   subroutine compute_pressure_or_deriv(x, val, dval_dx, d2val_dx2)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     PetscReal, intent(out), optional :: d2val_dx2
     !
     PetscReal, parameter             :: g   = GRAVITY_CONSTANT
     PetscReal, parameter             :: a0 =  15000.d0
     PetscReal, parameter             :: a1 = -20000.d0

     if (present(val      )) val       =  a0                *sin((x-x_min)/xlim*pi) + a1 + PRESSURE_REF
     if (present(dval_dx  )) dval_dx   =  a0*PI/xlim        *cos((x-x_min)/xlim*pi)
     if (present(d2val_dx2)) d2val_dx2 = -a0*PI*PI/xlim/xlim*sin((x-x_min)/xlim*pi)

     !if (present(val      )) val       =  a0                *cos((x-x_min)/xlim*pi) + a1 + PRESSURE_REF
     !if (present(dval_dx  )) dval_dx   = -a0*PI/xlim        *sin((x-x_min)/xlim*pi)
     !if (present(d2val_dx2)) d2val_dx2 = -a0*PI*PI/xlim/xlim*cos((x-x_min)/xlim*pi)

     !if (present(val      )) val       =  a0                          *cos((x-x_min)/xlim*pi*4.d0) + a1 + PRESSURE_REF
     !if (present(dval_dx  )) dval_dx   = -a0*PI*4.d0/xlim             *sin((x-x_min)/xlim*pi*4.d0)
     !if (present(d2val_dx2)) d2val_dx2 = -a0*PI*PI*4.d0*4.d0/xlim/xlim*cos((x-x_min)/xlim*pi*4.d0)

     !if (present(val      )) val       =  a0*((x-xlim)/xlim) + a1 + PRESSURE_REF
     !if (present(dval_dx  )) dval_dx   =  a0*(1.d0/xlim)
     !if (present(d2val_dx2)) d2val_dx2 =  0.d0
     
   end subroutine compute_pressure_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_temperature_or_deriv(x, val, dval_dx, d2val_dx2)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     PetscReal, intent(out), optional :: d2val_dx2
     !
     PetscReal, parameter             :: a0 = 5.d0
     PetscReal, parameter             :: a1 = 290.d0

     if (present(val      )) val       =  a0                *cos((x-x_min)/xlim*pi) + a1
     if (present(dval_dx  )) dval_dx   = -a0*PI/xlim        *sin((x-x_min)/xlim*pi)
     if (present(d2val_dx2)) d2val_dx2 = -a0*PI*PI/xlim/xlim*cos((x-x_min)/xlim*pi)
     
     if (present(val      )) val       =  a0                *sin((x-x_min)/xlim*pi)*1.d0 + a1
     if (present(dval_dx  )) dval_dx   =  a0*PI/xlim        *cos((x-x_min)/xlim*pi)*1.d0
     if (present(d2val_dx2)) d2val_dx2 = -a0*PI*PI/xlim/xlim*sin((x-x_min)/xlim*pi)*1.d0
     
   end subroutine compute_temperature_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_permeability_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: p0 = 1.d-11 ! [m^2]

     if (present(val    )) val     = p0        *(2.d0 - cos((x-x_min)/xlim * pi))
     if (present(dval_dx)) dval_dx = p0*PI/xlim*(     + sin((x-x_min)/xlim * pi))

   end subroutine compute_permeability_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_alpha_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 1.d0/4000.d0 ! [Pa^{-1}]

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_alpha_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_lambda_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 0.5d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_lambda_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_residualsat_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 0.0d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_residualsat_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_dryconductivity_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 0.25d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_dryconductivity_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_wetconductivity_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 1.3d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_wetconductivity_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_thermalalpha_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 0.45d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_thermalalpha_or_deriv

  !------------------------------------------------------------------------
  subroutine set_variable_for_problem(data_type, data_1D)
    !
    ! !DESCRIPTION:
    !
    use EOSWaterMod              , only : Density, Viscosity, InternalEnergyAndEnthalpy
    use MultiPhysicsProbConstants, only : FMWH2O
    use SaturationFunction       , only : saturation_params_type, SatFunc_Set_VG, SatFunc_PressToSat, SatFunc_PressToRelPerm
    !
    implicit none
    !
    PetscInt              :: data_type
    PetscReal, pointer    :: data_1D(:)
    !
    PetscInt              :: ii,jj,kk
    PetscInt              :: count
    PetscReal             :: x
    PetscReal             :: P     , dP_dx, d2P_dx2
    PetscReal             :: T     , dT_dx, d2T_dx2
    PetscReal             :: p0    , dp0_dx
    PetscReal             :: k     , dk_dx
    PetscReal             :: ds0   , ds0_dx
    PetscReal             :: dsr   , dsr_dx
    PetscReal             :: kr    , dkr_dx, dkr_dse
    PetscReal             :: rho   , drho_dP, d2rho_dP2, drho_dT, d2rho_dT2, drho_dx, d2rho_dx2
    PetscReal             :: mu    , dmu_dP
    PetscReal             :: dse_dP, dse_dp0
    PetscReal             :: m
    PetscReal             :: dmu_dT
    PetscReal             :: sat_res, se
    PetscReal             :: dkr_dP
    PetscReal             :: rhoq, drhoq_dx
    PetscReal             :: kappa, dkappa_dx
    PetscReal             :: alpha
    PetscReal             :: dH_dx
    PetscReal             :: Ke, dKe_dx, kdry, dkdry_dx, kwet, dkwet_dx
    PetscReal             :: U, H, dU_dT, dH_dT, dU_dP, dH_dP
    PetscReal             :: x_ppert, P_ppert, T_ppert
    PetscReal             :: x_npert, P_npert, T_npert
    PetscReal             :: rho_ppert, drho_ppert_dP, drho_ppert_dT
    PetscReal             :: rho_npert, drho_npert_dP, drho_npert_dT
    PetscReal             :: U_ppert, H_ppert, dU_ppert_dT, dH_ppert_dT, dU_ppert_dP, dH_ppert_dP, se_ppert
    PetscReal             :: U_npert, H_npert, dU_npert_dT, dH_npert_dT, dU_npert_dP, dH_npert_dP, se_npert
    PetscReal , parameter :: g = GRAVITY_CONSTANT
    PetscReal , parameter :: pert = 1.d-6
    type(saturation_params_type) :: satParam
    !-----------------------------------------------------------------------

    jj = 1; kk = 1;

    select case(data_type)
    case (DATA_PRESSURE)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_pressure_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_INITIAL_PRESSURE)
       P = 0.d0
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_pressure_or_deriv(x, val=data_1D(ii))
          P = P + 1.d0/nx*data_1D(ii)
       end do
       data_1D(:) = P

    case (DATA_POROSITY)
       do ii = 1, nx
          data_1D(ii) = 0.d0
       end do

    case (DATA_PERMEABILITY)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_permeability_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_SATFUNC_ALPHA)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_alpha_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_SATFUNC_LAMBDA)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_lambda_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_RES_SAT)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_residualsat_or_deriv(x, val=data_1D(ii))
       end do

     case (DATA_PRESSURE_BC)
        count = 0

        ii = 1; jj = 1; kk = 1;
        x     = soil_xc_3d(ii,jj,kk)-dx/2.d0
        count = count + 1;
        call compute_pressure_or_deriv(x,val=data_1D(count))

        ii = nx; jj = 1; kk = 1;
        x     = soil_xc_3d(ii,jj,kk)+dx/2.d0
        count = count + 1;
        call compute_pressure_or_deriv(x,val=data_1D(count))

    case (DATA_MASS_SOURCE)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          x_ppert= x + pert
          x_npert= x - pert

          call compute_permeability_or_deriv(x, val=k , dval_dx=dk_dx)
          call compute_alpha_or_deriv       (x, val=p0, dval_dx=dp0_dx)
          call compute_lambda_or_deriv      (x, val=m)
          call compute_pressure_or_deriv    (x, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2)
          call compute_temperature_or_deriv (x, val=T, dval_dx=dT_dx)
          call compute_residualsat_or_deriv (x, val=sat_res)

          call Viscosity(P, T, mu, dmu_dP, dmu_dT)
          call Density(  P, T, density_type,  rho, drho_dP, drho_dT)

          rho     = rho     * FMWH2O
          drho_dP = drho_dP * FMWH2O
          d2rho_dP2 = 0.d0
          d2rho_dT2 = 0.d0

          call SatFunc_Set_VG(satParam, sat_res, p0, m)
          call SatFunc_PressToSat(satParam, P, se, dse_dP)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)
          
          dkr_dse   = &
               0.5d0 * se**(-0.5d0) *( 1.d0 - (1.d0 - se **(1.d0/m))**m)**2.d0 + &
               se**(0.5d0) * 2.d0   *( 1.d0 - (1.d0 - se **(1.d0/m))**m) * (1.d0 - se**(1.d0/m)) * se**(1.d0/m - 1.d0)
          dse_dp0   = 0.d0
          
          dkr_dx    = dkr_dP * dP_dx + dkr_dse * dse_dp0 * dp0_dx

          call compute_pressure_or_deriv    (x_ppert, val=P_ppert)
          call compute_temperature_or_deriv (x_ppert, val=T_ppert)
          call Density(  P_ppert, T_ppert, density_type,  rho_ppert, drho_ppert_dP, drho_ppert_dT)
          call compute_pressure_or_deriv    (x_npert, val=P_npert)
          call compute_temperature_or_deriv (x_npert, val=T_npert)
          call Density(  P_npert, T_npert, density_type,  rho_npert, drho_npert_dP, drho_npert_dT)

          rho_ppert     = rho_ppert * FMWH2O
          rho_npert     = rho_npert * FMWH2O

          drho_dx   = (rho_ppert - rho_npert)/pert/2.d0
          d2rho_dx2 = (rho_ppert + rho_npert - 2.d0*rho)/pert/pert

          data_1D(ii) = &
               -((k*kr/mu)*drho_dx + (rho*kr/mu)*dk_dx + (rho*k/mu)*dkr_dx)*(dP_dx) &
               -(rho*k*kr/mu)*(d2P_dx2)
          data_1D(ii) = data_1D(ii)*dx
          
       end do

    case (DATA_HEAT_CAP)
       do ii = 1, nx
          data_1D(ii) = 0.d0
       end do

    case (DATA_T_COND_DRY)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_dryconductivity_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_T_COND_WET)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_wetconductivity_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_T_ALPHA)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_thermalalpha_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_TEMPERATURE)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_temperature_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_INITIAL_TEMPERATURE)
       T = 0.d0
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_temperature_or_deriv(x, val=data_1D(ii))
          T = T + 1.d0/nx*data_1D(ii)
       end do
       data_1D(:) = T

     case (DATA_TEMPERATURE_BC)
        count = 0

        ii = 1; jj = 1; kk = 1;
        x     = soil_xc_3d(ii,jj,kk)-dx/2.d0
        count = count + 1;
        call compute_temperature_or_deriv(x,val=data_1D(count))

        ii = nx; jj = 1; kk = 1;
        x     = soil_xc_3d(ii,jj,kk)+dx/2.d0
        count = count + 1;
        call compute_temperature_or_deriv(x,val=data_1D(count))

    case (DATA_HEAT_SOURCE)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          x_ppert= x + pert
          x_npert= x - pert

          call compute_permeability_or_deriv(x, val=k , dval_dx=dk_dx)
          call compute_alpha_or_deriv       (x, val=p0, dval_dx=dp0_dx)
          call compute_lambda_or_deriv      (x, val=m)
          call compute_pressure_or_deriv    (x, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2)
          call compute_temperature_or_deriv (x, val=T, dval_dx=dT_dx, d2val_dx2=d2T_dx2)
          call compute_residualsat_or_deriv (x, val=sat_res)
          call compute_dryconductivity_or_deriv(x, val=kdry, dval_dx=dkdry_dx)
          call compute_wetconductivity_or_deriv(x, val=kwet, dval_dx=dkwet_dx)
          call compute_thermalalpha_or_deriv(x, val=alpha)

          call Viscosity(P, T, mu, dmu_dP, dmu_dT)
          call Density(  P, T, density_type,  rho, drho_dP, drho_dT)

          rho     = rho     * FMWH2O
          drho_dP = drho_dP * FMWH2O
          drho_dT = drho_dT * FMWH2O
          d2rho_dP2 = 0.d0
          d2rho_dT2 = 0.d0

          call InternalEnergyAndEnthalpy(P, T, int_energy_enthalpy_type, rho, drho_dT, drho_dP, &
               U, H, dU_dT, dH_dT, dU_dP, dH_dP)

          call SatFunc_Set_VG(satParam, sat_res, p0, m)
          call SatFunc_PressToSat(satParam, P, se, dse_dP)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)

          dkr_dse   = &
               0.5d0 * se**(-0.5d0) *( 1.d0 - (1.d0 - se **(1.d0/m))**m)**2.d0 + &
               se**(0.5d0) * 2.d0   *( 1.d0 - (1.d0 - se **(1.d0/m))**m) * (1.d0 - se**(1.d0/m)) * se**(1.d0/m - 1.d0)
          dse_dp0   = 0.d0
          
          dkr_dx    = dkr_dP * dP_dx !+ dkr_dse * dse_dp0 * dp0_dx

          call compute_pressure_or_deriv    (x_ppert, val=P_ppert)
          call compute_temperature_or_deriv (x_ppert, val=T_ppert)
          call Density(  P_ppert, T_ppert, density_type,  rho_ppert, drho_ppert_dP, drho_ppert_dT)
          call compute_pressure_or_deriv    (x_npert, val=P_npert)
          call compute_temperature_or_deriv (x_npert, val=T_npert)
          call Density(  P_npert, T_npert, density_type,  rho_npert, drho_npert_dP, drho_npert_dT)

          rho_ppert     = rho_ppert * FMWH2O
          rho_npert     = rho_npert * FMWH2O
          drho_ppert_dP = drho_ppert_dP * FMWH2O
          drho_ppert_dT = drho_ppert_dT * FMWH2O
          drho_npert_dP = drho_npert_dP * FMWH2O
          drho_npert_dT = drho_npert_dT * FMWH2O

          drho_dx   = (rho_ppert - rho_npert)/pert/2.d0
          d2rho_dx2 = (rho_ppert + rho_npert - 2.d0*rho)/pert/pert

          rhoq      = -rho*(k*kr/mu*dP_dx)
          drhoq_dx  = &
               -((k*kr/mu)*drho_dx + (rho*kr/mu)*dk_dx + (rho*k/mu)*dkr_dx)*(dP_dx) &
               -(rho*k*kr/mu)*(d2P_dx2)


          Ke        = (se + 1.d-6)**(alpha)
          dKe_dx    = alpha*((se+1.d-6)**(alpha-1.d0))*dse_dP*dP_dx
          call SatFunc_PressToSat(satParam, P_ppert, se_ppert, dse_dP)
          call SatFunc_PressToSat(satParam, P_npert, se_npert, dse_dP)
          dKe_dx    = ( (se_ppert + 1.d-6)**(alpha) - (se_npert + 1.d-6)**(alpha))/pert/2.d0

          kappa     = kwet*Ke + kdry*(1.d0-Ke)
          dkappa_dx = dkwet_dx*Ke + dkdry_dx*(1.d0-Ke) + (kwet-kdry)*dKe_dx

          call InternalEnergyAndEnthalpy(P_ppert, T_ppert, int_energy_enthalpy_type, rho_ppert, drho_ppert_dT, drho_ppert_dP, &
               U_ppert, H_ppert, dU_ppert_dT, dH_ppert_dT, dU_ppert_dP, dH_ppert_dP)
          call InternalEnergyAndEnthalpy(P_npert, T_npert, int_energy_enthalpy_type, rho_npert, drho_npert_dT, drho_npert_dP, &
               U_npert, H_npert, dU_npert_dT, dH_npert_dT, dU_npert_dP, dH_npert_dP)

          dH_dx     = dH_dP*dP_dx + dH_dT*dT_dx
          dH_dx     = (H_ppert - H)/pert
          dH_dx     = (H_ppert - H_npert)/pert/2.d0

          data_1D(ii) = -(drhoq_dx*H/FMWH2O + rhoq*dH_dx/FMWH2O - dkappa_dx*dT_dx - kappa*d2T_dx2)*dx

       end do

    case default
       write(*,*)data_type
       write(iulog,*)'Unsupported data type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine set_variable_for_problem

  !------------------------------------------------------------------------
  subroutine save_true_solution(true_soln_filename)
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
    character(len=256)                   :: true_soln_filename
    PetscReal , pointer :: v_p(:)
    Vec                                  :: true_soln
    PetscViewer                          :: viewer
    character(len=256)                   :: string
    PetscErrorCode                       :: ierr
    PetscInt            :: nDM
    DM        , pointer :: dms(:)
    Vec       , pointer :: soln_subvecs(:)
    PetscInt            :: ii

    string = trim(true_soln_filename)

    call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)

    call VecDuplicate(th_mpp%soe%solver%soln, true_soln, ierr); CHKERRQ(ierr)

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(th_mpp%soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(th_mpp%soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(th_mpp%soe%solver%dm, &
         true_soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)
       if (ii == 1) then
          call set_variable_for_problem(DATA_PRESSURE, v_p)
       else
          call set_variable_for_problem(DATA_TEMPERATURE, v_p)
       end if
       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(th_mpp%soe%solver%dm, &
         true_soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)

    call VecView(true_soln,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  end subroutine save_true_solution

  !------------------------------------------------------------------------
  subroutine save_source(source_filename)
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
    character(len=256)                   :: source_filename
    PetscReal , pointer :: v_p(:)
    Vec                                  :: source
    PetscViewer                          :: viewer
    character(len=256)                   :: string
    PetscErrorCode                       :: ierr
    PetscInt            :: nDM
    DM        , pointer :: dms(:)
    Vec       , pointer :: soln_subvecs(:)
    PetscInt            :: ii

    string = trim(source_filename)

    call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)

    call VecDuplicate(th_mpp%soe%solver%soln, source, ierr); CHKERRQ(ierr)

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(th_mpp%soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(th_mpp%soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(th_mpp%soe%solver%dm, &
         source, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)
       if (ii == 1) then
          call set_variable_for_problem(DATA_MASS_SOURCE, v_p)
       else
          call set_variable_for_problem(DATA_HEAT_SOURCE, v_p)
       end if
       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(th_mpp%soe%solver%dm, &
         source, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)

    call VecView(source,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  end subroutine save_source

  !------------------------------------------------------------------------
  subroutine save_perm(perm_filename)
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
    character(len=256)                   :: perm_filename
    PetscReal , pointer :: v_p(:)
    Vec                                  :: source
    PetscViewer                          :: viewer
    character(len=256)                   :: string
    PetscErrorCode                       :: ierr
    PetscInt            :: nDM
    DM        , pointer :: dms(:)
    Vec       , pointer :: soln_subvecs(:)
    PetscInt            :: ii

    string = trim(perm_filename)

    call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)

    call VecCreate(PETSC_COMM_SELF, source, ierr); CHKERRQ(ierr)
    call VecSetSizes(source, ncells, PETSC_DECIDE, ierr); CHKERRQ(ierr)
    call VecSetFromOptions(source, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(source, v_p, ierr)
    call set_variable_for_problem(DATA_PERMEABILITY, v_p)
    call VecRestoreArrayF90(source, v_p, ierr)

    call VecView(source,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  end subroutine save_perm

end module th_mms_problem
