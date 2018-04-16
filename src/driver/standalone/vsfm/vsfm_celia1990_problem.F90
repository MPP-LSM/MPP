
module vsfm_celia1990_problem

  implicit none

#include <petsc/finclude/petsc.h>
  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscInt              :: nz
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost

  public :: run_vsfm_celia1990_problem
  public :: output_regression_vsfm_celia1990_problem

contains
!------------------------------------------------------------------------

  subroutine run_vsfm_celia1990_problem()
    !
#include <petsc/finclude/petsc.h>
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use mpp_varpar           , only : mpp_varpar_init
    use petscsys
    use petscvec
    use petscmat
    use petscts
    use petscsnes
    use petscdm
    use petscdmda
    !
    implicit none
    !
    !
    !
    PetscBool          :: converged
    PetscInt           :: converged_reason
    PetscErrorCode     :: ierr
    PetscReal          :: dtime
    PetscInt           :: istep, nstep
    PetscBool           :: flg
    PetscBool          :: save_initial_soln, save_final_soln
    character(len=256) :: string
    character(len=256) :: output_suffix
    PetscViewer        :: viewer

    ! Set default settings
    nz                = 100
    dtime             = 3600.d0
    nstep             = 24
    save_initial_soln = PETSC_FALSE
    save_final_soln   = PETSC_FALSE
    output_suffix     = ''

    ! Get some command line options

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nz',nz,flg,ierr)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-dt',dtime,flg,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_initial_soln',save_initial_soln,flg,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-save_final_soln',save_final_soln,flg,ierr)
    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-output_suffix',output_suffix,flg,ierr)

    ! Initialize the problem
    call Init()  

    if (save_initial_soln) then
       if (len(trim(adjustl(output_suffix))) ==  0) then
          string = trim(output_suffix) // 'initial_soln.bin'
       else
          string = 'initial_soln_' // trim(output_suffix) // '.bin'
       endif

       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(vsfm_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

    do istep = 1, nstep
       ! Update BC
       call set_bondary_conditions()

       ! Run the model
       call vsfm_mpp%soe%StepDT(dtime, istep, &
            converged, converged_reason, ierr); CHKERRQ(ierr)
    enddo

    if (save_final_soln) then
       if (len(trim(adjustl(output_suffix))) ==  0) then
          string = trim(output_suffix) // 'final_soln.bin'
       else
          string = 'final_soln_' // trim(output_suffix) // '.bin'
       endif
       call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
       call VecView(vsfm_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
       call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif

  end subroutine run_vsfm_celia1990_problem

  !------------------------------------------------------------------------
  subroutine Init()
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
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
    call vsfm_mpp%SetupProblem()

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
    use MultiPhysicsProbConstants , only : MPP_VSFM_SNES_CLM
    use MultiPhysicsProbVSFM , only : vsfm_mpp
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
    call vsfm_mpp%Init       ()
    call vsfm_mpp%SetName    ('Variably-Saturated-Flow-Model')
    call vsfm_mpp%SetID      (MPP_VSFM_SNES_CLM)
    call vsfm_mpp%SetMPIRank (iam)

  end subroutine initialize_mpp

  !------------------------------------------------------------------------
  subroutine add_meshes()
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : CONN_IN_Z_DIR
    use mpp_varpar                , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    use MeshType                  , only : mesh_type, MeshCreate
    !
    implicit none
    !
    class(mesh_type), pointer :: mesh
    PetscInt                  :: imesh
    PetscInt                  :: nlev

    call mpp_varpar_set_nlevsoi(nz)
    call mpp_varpar_set_nlevgrnd(nz)

    imesh        = 1
    ncells_local = nx*ny*nz
    
    !
    ! Set up the meshes
    !    
    call vsfm_mpp%SetNumMeshes(1)

    allocate(mesh)
    
    ! Configuration-1
    call MeshCreate(mesh, 'Soil mesh', x_column, y_column, z_column, &
         nx, ny, nz, &
         CONN_IN_Z_DIR)

    call vsfm_mpp%AddMesh(imesh, mesh)
    
    call mesh%Clean()

    deallocate(mesh)

  end subroutine add_meshes

  !------------------------------------------------------------------------

  subroutine add_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : GE_RE
    use MultiPhysicsProbConstants , only : MESH_CLM_SOIL_COL
    !
    ! !ARGUMENTS
    implicit none

    call vsfm_mpp%AddGovEqn(GE_RE, 'Richards Equation ODE', MESH_CLM_SOIL_COL)

    call vsfm_mpp%SetMeshesOfGoveqns()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns()
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : SOIL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    !
    ! !ARGUMENTS
    implicit none
    !
    PetscInt :: ieqn

    ieqn = 1

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Constant head condition at top', 'Pa', COND_DIRICHLET, &
         SOIL_TOP_CELLS)

    call vsfm_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
         'Constant head condition at bottom', 'Pa', COND_DIRICHLET, &
         SOIL_BOTTOM_CELLS)

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine allocate_auxvars()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM , only : vsfm_mpp
    !
    implicit none

    !
    ! Allocate auxvars
    !
    call vsfm_mpp%AllocateAuxVars()

  end subroutine allocate_auxvars

  !------------------------------------------------------------------------
  subroutine set_material_properties()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbVSFM      , only : VSFMMPPSetSoils
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
    use EOSWaterMod               , only : DENSITY_TGDPB01
    use mpp_varcon                , only : denh2o
    use mpp_varcon              , only : grav
    !
    implicit none
    !
    PetscReal        , parameter :: porosity            = 0.368d0
    PetscReal        , parameter :: lambda              = 0.5d0
    PetscReal        , parameter :: alpha               = 3.4257d-4
    PetscReal        , parameter :: perm                = 8.3913d-12
    !
    PetscReal        , pointer   :: vsfm_watsat(:,:)
    PetscReal        , pointer   :: vsfm_hksat(:,:)
    PetscReal        , pointer   :: vsfm_bsw(:,:)
    PetscReal        , pointer   :: vsfm_sucsat(:,:)
    PetscReal        , pointer   :: vsfm_eff_porosity(:,:)
    PetscReal        , pointer   :: vsfm_residual_sat(:,:)
    PetscReal        , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
    PetscInt                     :: begc , endc
    integer          , pointer   :: vsfm_filter(:)
    character(len=32)            :: satfunc_type
    !-----------------------------------------------------------------------

    begc = 1
    endc = 1

    satfunc_type = 'van_genuchten'

    allocate(vsfm_watsat       (1,nz))
    allocate(vsfm_hksat        (1,nz))
    allocate(vsfm_bsw          (1,nz))
    allocate(vsfm_sucsat       (1,nz))
    allocate(vsfm_eff_porosity (1,nz))
    allocate(vsfm_residual_sat (1,nz))
    allocate(vsfm_filter       (nz))

    ! Soil properties
    vsfm_filter(:)         = 1
    vsfm_watsat(:,:)       = porosity
    vsfm_hksat(:,:)        = perm /vish2o * (denh2o * grav) / 0.001d0
    vsfm_bsw(:,:)          = 1.d0/lambda
    vsfm_sucsat(:,:)       = 1.d0/(alpha*GRAVITY_CONSTANT)
    vsfm_eff_porosity(:,:) = 0.d0
    vsfm_residual_sat(:,:) = 0.2772d0

    call VSFMMPPSetSoils(vsfm_mpp, begc, endc, &
         ncells_ghost, vsfm_filter, &
         vsfm_watsat, vsfm_hksat, vsfm_bsw, &
         vsfm_sucsat, vsfm_residual_sat, &
         satfunc_type, DENSITY_TGDPB01)

    deallocate(vsfm_watsat)
    deallocate(vsfm_hksat   )
    deallocate(vsfm_bsw   )
    deallocate(vsfm_sucsat)
    deallocate(vsfm_eff_porosity)
    deallocate(vsfm_filter      )
    deallocate(vsfm_residual_sat)

  end subroutine set_material_properties

  !------------------------------------------------------------------------
  subroutine set_initial_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    !
    implicit none
    !
    PetscReal, pointer :: press_ic(:)

    allocate(press_ic(nz))

    press_ic(:) = 3.5355d3

    call vsfm_mpp%Restart(press_ic)

    deallocate(press_ic)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine set_bondary_conditions()
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_BC, VAR_BC_SS_CONDITION
    !
    implicit none
    !
    PetscReal, pointer :: top_pressure_bc(:)
    PetscReal, pointer :: bot_pressure_bc(:)
    PetscInt           :: soe_auxvar_id

    allocate(top_pressure_bc(1))
    allocate(bot_pressure_bc(1))

    top_pressure_bc(:) = 9.3991d4
    bot_pressure_bc(:) = 3.5355d3

    soe_auxvar_id = 1
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, top_pressure_bc)

    soe_auxvar_id = 2
    call vsfm_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
         VAR_BC_SS_CONDITION, soe_auxvar_id, bot_pressure_bc)

    deallocate(top_pressure_bc)
    deallocate(bot_pressure_bc)

  end subroutine set_bondary_conditions

  !------------------------------------------------------------------------
  subroutine output_regression_vsfm_celia1990_problem(filename_base, num_cells)
    !
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use regression_mod            , only : regression_type
    !
    implicit none
    !
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

    name = 'liquid_pressure'
    category = 'pressure'
    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,  &
         VAR_PRESSURE, -1, data)
    call regression%WriteData(name, category, data)

    name = 'liquid_saturation'
    category = 'general'
    call vsfm_mpp%soe%GetDataForCLM(AUXVAR_INTERNAL,  &
         VAR_LIQ_SAT, -1, data)
    call regression%WriteData(name, category, data)

    call regression%CloseOutput()
    
    deallocate(data)

  end subroutine output_regression_vsfm_celia1990_problem
  
end module vsfm_celia1990_problem
  
