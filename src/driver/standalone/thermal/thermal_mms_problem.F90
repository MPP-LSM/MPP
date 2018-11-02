
module thermal_mms_problem

  use mpp_varctl              , only : iulog
  use MultiPhysicsProbThermal , only : thermal_mpp
  use mpp_abortutils          , only : endrun
  use mpp_shr_log_mod         , only : errMsg => shr_log_errMsg
  use thermal_mms_vars
  use thermal_mms_steady_state_problem_1D
  use thermal_mms_steady_state_problem_2D
  use thermal_mms_steady_state_problem_3D

  implicit none

#include <petsc/finclude/petsc.h>

  !
  PetscInt              :: problem_type
  
  public :: run_thermal_mms_problem
  public :: output_regression_th_mms_problem

contains
!------------------------------------------------------------------------

  subroutine run_thermal_mms_problem(namelist_filename)
    !
#include <petsc/finclude/petsc.h>
    !
    use petscsys
    use petscvec
    use petscmat
    use petscts
    use petscsnes
    use petscdm
    use petscdmda
    use mpp_varcon                , only : cnfac
    !
    implicit none
    !
    character(len=256), optional :: namelist_filename
    !
    !
    PetscBool          :: converged
    PetscInt           :: converged_reason
    PetscErrorCode     :: ierr
    PetscReal          :: dtime
    PetscInt           :: istep, nstep
    PetscBool          :: flg
    character(len=256) :: string
    character(len=256) :: true_soln_fname
    character(len=256) :: source_fname
    character(len=256) :: thermal_cond_fname
    PetscViewer        :: viewer
    character(len=256) :: ioerror_msg
    character(len=2560):: namelist_buffer
    integer            :: nml_unit, nml_error

    namelist / problem_options / nx, ny, nz, problem_type

    ! Set default settings
    problem_type      = STEADY_STATE_1D
    nx                = 20
    ny                = 1
    nz                = 1
    dtime             = 3600.d0
    nstep             = 1
    true_soln_fname   = ''
    source_fname      = ''
    thermal_cond_fname= ''

    cnfac = 0.d0 ! For steady state equation need to set this to 0.d0
    
    ! Get some command line options

    call PetscOptionsGetInt    (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-problem_type      ',problem_type      ,flg,ierr)
    select case(problem_type)
    case (STEADY_STATE_1D)
       nx = 20; ny = 1; nz = 1;
    case (STEADY_STATE_2D)
       nx = 20; ny = 20; nz = 1;       
    case (STEADY_STATE_3D)
       nx = 20; ny = 20; nz = 20;
    case default
       write(iulog,*)'Unknown -problem_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
    
    call PetscOptionsGetInt    (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-nx                ',nx                ,flg,ierr)
    call PetscOptionsGetInt    (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-ny                ',ny                ,flg,ierr)
    call PetscOptionsGetInt    (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-nz                ',nz                ,flg,ierr)
    call PetscOptionsGetReal   (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-dt                ',dtime             ,flg,ierr)
    call PetscOptionsGetInt    (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-nstep             ',nstep             ,flg,ierr)
    call PetscOptionsGetString (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-view_true_solution',true_soln_fname   ,flg,ierr)
    call PetscOptionsGetString (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-view_source'       ,source_fname      ,flg,ierr)
    call PetscOptionsGetString (PETSC_NULL_OPTIONS ,PETSC_NULL_CHARACTER ,'-view_thermal_cond' ,thermal_cond_fname,flg,ierr)

    ! Read namelist, if a namelist file is provided
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

    select case(problem_type)
    case (STEADY_STATE_1D)
       if (nx <= 1) then
          write(iulog,*)'For STEADY_STATE_1D problem, need to set nx such that  nx > 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    
       if (ny /= 1) then
          write(iulog,*)'For STEADY_STATE_1D problem, need to set ny = 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    
       if (nz /= 1) then
          write(iulog,*)'For STEADY_STATE_1D problem, need to set nz = 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if

    case (STEADY_STATE_2D)
       if (nx <= 1) then
          write(iulog,*)'For STEADY_STATE_2D problem, need to set nx such that  nx > 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    
       if (ny < 1) then
          write(iulog,*)'For STEADY_STATE_2D problem, need to set ny such that  ny > 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    
       if (nz /= 1) then
          write(iulog,*)'For STEADY_STATE_2D problem, need to set nz = 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
       
    case (STEADY_STATE_3D)
       if (nx < 1) then
          write(iulog,*)'For STEADY_STATE_3D problem, need to set nx such that  nx > 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    
       if (ny < 1) then
          write(iulog,*)'For STEADY_STATE_3D problem, need to set ny such that  ny > 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    
       if (nz < 1) then
          write(iulog,*)'For STEADY_STATE_3D problem, need to set nz such that  nz > 1'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end if
    end select
       
    
    ! Initialize the problem
    call Init()

    ! Preform Pre-StepDT operations
    call thermal_mpp%soe%PreStepDT()

    do istep = 1, nstep

       if (istep == nstep) then
          if (len(trim(adjustl(true_soln_fname)))>0) then
             call save_problem_variable(true_soln_fname, DATA_TEMPERATURE)
          end if
       end if

       if (istep == nstep) then
          if (len(trim(adjustl(source_fname)))>0) then
             call save_problem_variable(source_fname, DATA_HEAT_SOURCE)
          end if
       end if

       if (istep == nstep) then
          if (len(trim(adjustl(thermal_cond_fname)))>0) then
             call save_problem_variable(thermal_cond_fname, DATA_THERMAL_CONDUCTIVITY_1D_VEC)
          end if
       end if

       call set_source_sink_conditions()
       call set_boundary_conditions()

       call thermal_mpp%soe%StepDT(1.d0, 1, &
           converged, converged_reason, ierr); CHKERRQ(ierr)
    end do

  end subroutine run_thermal_mms_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    implicit none
    !

    ! 1. Initialize the multi-physics-problem (MPP)
    call initialize_mpp()

    ! 2. Add all meshes needed for the MPP
    call add_meshes()

    ! 3. Add all governing equations
    call add_goveqns()

    ! 4. Add boundary and source-sink conditions to all governing equations
    call add_conditions_to_goveqns()

    ! 5. Allocate memory to hold auxvars
    call allocate_auxvars()

    ! 6. Setup the MPP
    call thermal_mpp%SetupProblem()

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
    !
#include <petsc/finclude/petsc.h>
    ! !USES:
    use MultiPhysicsProbConstants , only : MPP_THERMAL_TBASED_KSP_CLM
    use petscsys
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
    call thermal_mpp%Init       ()
    call thermal_mpp%SetName    ('Thermal model for MMS')
    call thermal_mpp%SetID      (MPP_THERMAL_TBASED_KSP_CLM)
    call thermal_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes()
    !
    use MultiPhysicsProbConstants , only : MESH_CLM_THERMAL_SOIL_COL
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    use MeshType                  , only : mesh_type, MeshCreate
    use MultiPhysicsProbConstants , only : CONN_IN_XYZ_DIR
    use mpp_mesh_utils            , only : ComputeXYZCentroids, ComputeCellID
    !
    implicit none
    !
    class (mesh_type), pointer :: mesh
    PetscInt                   :: imesh
    PetscInt :: ii,jj,kk

    call mpp_varpar_set_nlevsoi(nz)
    call mpp_varpar_set_nlevgrnd(nz)

    imesh        = 1
    ncells_local = nx*ny*nz
    ncells_ghost = 0

    call thermal_mpp%SetNumMeshes(1)

    allocate(mesh)

    call MeshCreate(mesh, 'Soil mesh', x_column, y_column, z_column, &
         nx, ny, nz, CONN_IN_XYZ_DIR)
    call mesh%SetID(MESH_CLM_THERMAL_SOIL_COL)

    call thermal_mpp%AddMesh(imesh, mesh)

    call mesh%Clean()

    deallocate(mesh)

    dx = x_column/nx
    dy = y_column/ny
    dz = z_column/nz

    call ComputeXYZCentroids(nx, ny, nz, dx, dy, dz, &
         soil_xc_3d, soil_yc_3d, soil_zc_3d)

    call ComputeCellID(nx, ny, nz, soil_id)

  end subroutine add_meshes
  !------------------------------------------------------------------------

  subroutine add_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    use MultiPhysicsProbConstants , only : GE_THERM_SOIL_TBASED
    use MultiPhysicsProbConstants , only : MESH_CLM_THERMAL_SOIL_COL
    !
    ! !ARGUMENTS
    implicit none

    call thermal_mpp%AddGovEqn(GE_THERM_SOIL_TBASED, &
         'Thermal equation using temprature formulation in soil (KSP formulation)', &
         MESH_CLM_THERMAL_SOIL_COL)

    call thermal_mpp%SetMeshesOfGoveqns()

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
    use MultiPhysicsProbConstants , only : COND_HEAT_RATE
    use MultiPhysicsProbConstants , only : CONN_HORIZONTAL
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt                            :: ieqn, ii, jj, kk
    PetscInt  :: count, nconn
    PetscInt                  , pointer :: id_up    (:)
    PetscInt                  , pointer :: id_dn    (:)
    PetscReal                 , pointer :: dist_up  (:)
    PetscReal                 , pointer :: dist_dn  (:)
    PetscReal                 , pointer :: area     (:)
    PetscInt                  , pointer :: itype    (:)
    PetscInt                  , pointer :: tmp_3d      (:,:,:)
    PetscReal                 , pointer :: unit_vec (:,:)
    class(connection_set_type) , pointer :: conn_set

    ieqn = 1
    
    nconn = 0

    if (nx > 1) nconn = nconn + ny*nz*2
    if (ny > 1) nconn = nconn + nx*nz*2
    if (nz > 1) nconn = nconn + nx*ny*2
        
    allocate(id_up    (nconn ))
    allocate(id_dn    (nconn ))
    allocate(dist_up  (nconn ))
    allocate(dist_dn  (nconn ))
    allocate(area     (nconn ))
    allocate(itype    (nconn ))
    allocate(unit_vec (nconn,3 ))

    id_up(:)      = 0
    unit_vec(:,:) = 0.d0
    count = 0
    
    ! X-direction
    if (nx > 1) then
       do kk = 1,nz
          do jj = 1,ny

             ! At the beginning
             ii = 1;
             count             = count + 1
             id_dn(count)      = soil_id(ii,jj,kk)
             dist_up(count)    = dx/2.d0*0.d0
             dist_dn(count)    = dx/2.d0
             area(count)       = dy*dz
             unit_vec(count,1) = 1.d0
             itype(count)      = CONN_HORIZONTAL

             ! At the end
             ii = nx
             count             = count + 1
             id_dn(count)      = soil_id(ii,jj,kk)
             dist_up(count)    = dx/2.d0*0.d0
             dist_dn(count)    = dx/2.d0
             area(count)       = dy*dz
             unit_vec(count,1) = -1.d0
             itype(count)      = CONN_HORIZONTAL
          end do
       end do
    end if

    ! Y-direction
    if (ny > 1) then
       do kk = 1,nz
          do ii = 1,nx

             ! At the beginning
             jj = 1;
             count             = count + 1
             id_dn(count)      = soil_id(ii,jj,kk)
             dist_up(count)    = 0.d0
             dist_dn(count)    = dy/2.d0
             area(count)       = dx*dz
             unit_vec(count,1) = 1.d0
             itype(count)      = CONN_HORIZONTAL

             ! At the end
             jj = ny
             count             = count + 1
             id_dn(count)      = soil_id(ii,jj,kk)
             dist_up(count)    = 0.d0
             dist_dn(count)    = dy/2.d0
             area(count)       = dx*dz
             unit_vec(count,1) = -1.d0
             itype(count)      = CONN_HORIZONTAL
          end do
       end do
    end if

    ! Z-direction
    if (nz > 1) then
       do jj = 1,ny
          do ii = 1,nx

             ! At the beginning
             kk = 1;
             count             = count + 1
             id_dn(count)      = soil_id(ii,jj,kk)
             dist_up(count)    = 0.d0
             dist_dn(count)    = dz/2.d0
             area(count)       = dx*dy
             unit_vec(count,1) = 1.d0
             itype(count)      = CONN_VERTICAL

             ! At the end
             kk = nz
             count             = count + 1
             id_dn(count)      = soil_id(ii,jj,kk)
             dist_up(count)    = 0.d0
             dist_dn(count)    = dz/2.d0
             area(count)       = dx*dy
             unit_vec(count,1) = -1.d0
             itype(count)      = CONN_VERTICAL
          end do
       end do
    end if

    call MeshCreateConnectionSet(thermal_mpp%meshes(1), &
         count, id_up, id_dn, &
         dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call thermal_mpp%soe%AddConditionInGovEqn( ieqn , &
         ss_or_bc_type      = COND_BC               , &
         name               = 'Temp BC'             , &
         unit               = 'T'                   , &
         cond_type          = COND_DIRICHLET        , &
         conn_set           = conn_set)

    call thermal_mpp%soe%AddConditionInGovEqn(ieqn, COND_SS, &
         'Source term for MMS', 'W/m^2', COND_HEAT_RATE,     &
         ALL_CELLS)

    deallocate(id_up    )
    deallocate(id_dn    )
    deallocate(dist_up  )
    deallocate(dist_dn  )
    deallocate(area     )
    deallocate(itype    )
    deallocate(unit_vec )

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    !
    implicit none

    !
    ! Allocate auxvars
    !
    call thermal_mpp%AllocateAuxVars()

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbThermal   , only : MPPThermalSetSoils
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use mpp_varcon                , only : denh2o
    use mpp_varcon                , only : grav
    use mpp_varpar                , only : nlevgrnd, nlevsno
    use mpp_varcon                , only : istsoil
    !
    implicit none
    !
    integer   , pointer :: thermal_filter(:)
    integer   , pointer :: thermal_lun_type(:)
    PetscReal , pointer :: thermal_watsat(:,:)
    PetscReal , pointer :: thermal_csol(:,:)
    PetscReal , pointer :: thermal_tkmg(:,:)
    PetscReal , pointer :: thermal_tkdry(:,:)
    PetscInt            :: begc , endc
    !-----------------------------------------------------------------------

    begc = 1
    endc = nx*ny

    nlevgrnd = nz

    ! Soil properties
    allocate (thermal_filter   (begc:endc      ))
    allocate (thermal_lun_type (begc:endc      ))
    allocate (thermal_watsat   (begc:endc,nlevgrnd ))
    allocate (thermal_csol     (begc:endc,nlevgrnd ))
    allocate (thermal_tkmg     (begc:endc,nlevgrnd ))
    allocate (thermal_tkdry    (begc:endc,nlevgrnd ))

    thermal_filter(:)   = 1
    thermal_lun_type(:) = istsoil
    thermal_watsat(:,:) = 0.1d0
    thermal_csol(:,:)   = 0.d0
    thermal_tkmg(:,:)   = 0.d0
    thermal_tkdry(:,:)  = 0.d0

    call set_variable_for_problem(problem_type, DATA_THERMAL_CONDUCTIVITY, data_2D=thermal_tkdry)
    
    call MPPThermalSetSoils(thermal_mpp, begc, &
                           endc,               &
                           thermal_filter,     &
                           thermal_lun_type,   &
                           thermal_watsat,     &
                           thermal_csol,       &
                           thermal_tkmg,       &
                           thermal_tkdry       &
                           )

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use SystemOfEquationsThermalType        , only : sysofeqns_thermal_type
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use MultiPhysicsProbConstants           , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants           , only : AUXVAR_BC
    use MultiPhysicsProbConstants           , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants           , only : VAR_TUNING_FACTOR
    use MultiPhysicsProbConstants           , only : VAR_LIQ_AREAL_DEN
    use MultiPhysicsProbConstants           , only : VAR_ACTIVE
    use MultiPhysicsProbConstants           , only : VAR_FRAC
    use ThermalKSPTemperatureSoilAuxMod     , only : ThermKSPTempSoilAuxVarSetRValues
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnThermalKSPTemperatureSoilType , only : goveqn_thermal_ksp_temp_soil_type
    !
    implicit none
    !
    PetscReal                  , pointer :: temperature_1d(:)
    PetscReal                  , pointer :: tuning_factor_1d(:)
    PetscReal                  , pointer :: liq_areal_den_1d(:)
    PetscInt                   , pointer :: ids(:)
    PetscReal                  , pointer :: is_active_1d(:)
    PetscReal                  , pointer :: fraction_1d(:)
    PetscReal                  , pointer :: temperature_bc(:)
    PetscReal                  , pointer :: heat_source(:)
    class(sysofeqns_base_type) , pointer :: base_soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: ii,jj,kk, count, nconn
    PetscInt                             :: soe_auxvar_id, gov_eqn_id
    PetscReal                            :: x, y

    allocate (temperature_1d   (nx*ny*nz ))
    allocate (tuning_factor_1d (nx*ny*nz ))
    allocate (liq_areal_den_1d (nx*ny*nz ))
    allocate (heat_source      (nx*ny*nz ))

    nconn = 0

    if (nx > 1) nconn = nconn + ny*nz*2
    if (ny > 1) nconn = nconn + nx*nz*2
    if (nz > 1) nconn = nconn + nx*ny*2

    allocate (temperature_bc   (nconn ))
    allocate (is_active_1d     (nconn ))
    allocate (fraction_1d      (nconn ))
    allocate (ids              (nconn ))

    temperature_1d(:)   = 290.d0
    tuning_factor_1d(:) = 1.d0
    liq_areal_den_1d(:) = 0.d0
    is_active_1d(:)     = 1.d0
    fraction_1d(:)      = 1.d0

    do ii = 1, nconn
       ids(ii) = ii
    end do

    call set_variable_for_problem(prob_type=problem_type, &
         data_type=DATA_TEMPERATURE_BC, data_1D=temperature_bc)

    base_soe => thermal_mpp%soe

    select type(base_soe)
    class is (sysofeqns_thermal_type)

       soe_auxvar_id = 1
       
       ! Set temperature
       call base_soe%SetSolnPrevCLM(temperature_1d)

       call base_soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_TUNING_FACTOR, soe_auxvar_id, tuning_factor_1d)

       call base_soe%SetRDataFromCLM(AUXVAR_INTERNAL, &
            VAR_LIQ_AREAL_DEN, soe_auxvar_id, liq_areal_den_1d)

       call base_soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, temperature_bc)


       gov_eqn_id = 1
       call base_soe%SetPointerToIthGovEqn(gov_eqn_id, cur_goveq)
       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_soil_type)
          call ThermKSPTempSoilAuxVarSetRValues(cur_goveq%aux_vars_bc, &
               VAR_ACTIVE, nconn, ids, is_active_1d)

          call ThermKSPTempSoilAuxVarSetRValues(cur_goveq%aux_vars_bc, &
               VAR_FRAC, nconn, ids, fraction_1d)
       end select

    end select

    deallocate (temperature_1d  )
    deallocate (tuning_factor_1d)
    deallocate (liq_areal_den_1d)    
    deallocate (temperature_bc  )
    deallocate (is_active_1d    )
    deallocate (fraction_1d     )
    deallocate (ids             )

  end subroutine set_initial_conditions

    !------------------------------------------------------------------------
  subroutine set_boundary_conditions()
    !
    ! !DESCRIPTION:
    !
    use SystemOfEquationsThermalType        , only : sysofeqns_thermal_type
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use MultiPhysicsProbConstants           , only : AUXVAR_BC
    use MultiPhysicsProbConstants           , only : VAR_BC_SS_CONDITION
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnThermalKSPTemperatureSoilType , only : goveqn_thermal_ksp_temp_soil_type
    !
    implicit none
    !
    PetscReal                  , pointer :: temperature_bc(:)
    class(sysofeqns_base_type) , pointer :: base_soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: nconn
    PetscInt                             :: soe_auxvar_id

    nconn = 0

    if (nx > 1) nconn = nconn + ny*nz*2
    if (ny > 1) nconn = nconn + nx*nz*2
    if (nz > 1) nconn = nconn + nx*ny*2

    allocate (temperature_bc   (nconn ))

    call set_variable_for_problem(prob_type=problem_type, &
         data_type=DATA_TEMPERATURE_BC, data_1D=temperature_bc)

    base_soe => thermal_mpp%soe

    select type(base_soe)
    class is (sysofeqns_thermal_type)

       soe_auxvar_id = 1
       
       call base_soe%SetRDataFromCLM(AUXVAR_BC, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, temperature_bc)

    end select

    deallocate (temperature_bc  )

  end subroutine set_boundary_conditions
  
  !------------------------------------------------------------------------
  subroutine set_source_sink_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants    , only : AUXVAR_SS
    use MultiPhysicsProbConstants    , only : VAR_BC_SS_CONDITION
    use SystemOfEquationsBaseType    , only : sysofeqns_base_type
    use SystemOfEquationsThermalType , only : sysofeqns_thermal_type
    !
    PetscInt                             :: soe_auxvar_id
    PetscReal                  , pointer :: source_sink(:)
    class(sysofeqns_base_type) , pointer :: base_soe
    PetscErrorCode                       :: ierr

    allocate(source_sink(nx*ny*nz))

    call set_variable_for_problem(problem_type, DATA_HEAT_SOURCE, data_1D=source_sink)

    base_soe => thermal_mpp%soe

    select type(base_soe)
    class is (sysofeqns_thermal_type)

       soe_auxvar_id = 1
       
       call base_soe%SetRDataFromCLM(AUXVAR_SS, &
            VAR_BC_SS_CONDITION, soe_auxvar_id, source_sink)
    end select

    deallocate(source_sink)
    
  end subroutine set_source_sink_conditions

  !------------------------------------------------------------------------
  subroutine save_problem_variable(true_soln_filename, variable_id)
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
    PetscInt                             :: variable_id
    !
    PetscInt                             :: soe_auxvar_id
    PetscReal                  , pointer :: val(:)
    Vec                                  :: true_soln
    PetscViewer                          :: viewer
    class(sysofeqns_base_type) , pointer :: base_soe
    character(len=256)                   :: string
    PetscErrorCode                       :: ierr

    base_soe => thermal_mpp%soe

    string = trim(true_soln_filename)

    call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)

    call VecDuplicate(thermal_mpp%soe%solver%soln, true_soln, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(true_soln, val, ierr); CHKERRQ(ierr)

    call set_variable_for_problem(problem_type, variable_id, data_1D=val)

    call VecRestoreArrayF90(true_soln, val, ierr); CHKERRQ(ierr)

    call VecView(true_soln,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  end subroutine save_problem_variable

  !------------------------------------------------------------------------
  subroutine set_variable_for_problem(prob_type, data_type, data_1D, data_2D)
    !
    ! !DESCRIPTION:
    !
    !
    implicit none
    !
    PetscInt                     :: prob_type
    PetscInt                     :: data_type
    PetscReal, pointer, optional :: data_1D(:)
    PetscReal, pointer, optional :: data_2D(:,:)
    !
    PetscInt                     :: ii, jj, kk
    PetscInt                     :: count
    PetscReal                    :: x, y, z
    !-----------------------------------------------------------------------

    select case(prob_type)
    case (STEADY_STATE_1D)
       call set_variable_for_steady_state_1D(data_type, data_1D, data_2D)
       
    case (STEADY_STATE_2D)
       call set_variable_for_steady_state_2D(data_type, data_1D, data_2D)
       
    case (STEADY_STATE_3D)
       call set_variable_for_steady_state_3D(data_type, data_1D, data_2D)
       
    case default
       write(iulog,*)'Unsupported problem type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine set_variable_for_problem


  !------------------------------------------------------------------------
  subroutine output_regression_th_mms_problem(filename_base, num_cells)
    !
    use MultiPhysicsProbConstants    , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants    , only : VAR_TEMPERATURE
    use regression_mod               , only : regression_type
    use SystemOfEquationsBaseType    , only : sysofeqns_base_type
    use SystemOfEquationsThermalType , only : sysofeqns_thermal_type
    !
    implicit none
    !
    class(sysofeqns_base_type) , pointer :: base_soe

    character(len=256)    :: filename_base
    character(len=512)    :: filename
    PetscInt, intent(in)  :: num_cells
    !
    PetscInt              :: output
    character(len=64)     :: category
    character(len=64)     :: name
    PetscReal, pointer    :: data(:)
    type(regression_type) :: regression

    allocate(data(ncells_local))

    call regression%Init(filename_base, num_cells)
    call regression%OpenOutput()

    name = 'temperature'
    category = 'general'
    base_soe => thermal_mpp%soe

    select type(base_soe)
    class is (sysofeqns_thermal_type)
       call base_soe%GetSoln(data)
    end select
    !call thermal_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,  &
    !     VAR_TEMPERATURE, -1, data)
    call regression%WriteData(name, category, data)

    call regression%CloseOutput()
    
    deallocate(data)

  end subroutine output_regression_th_mms_problem
  
end module thermal_mms_problem
  
