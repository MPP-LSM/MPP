
module mesh_info
#include <petsc/finclude/petsc.h>
  PetscInt  , parameter :: nx       = 1
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscInt              :: nz
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost
end module mesh_info

!------------------------------------------------------------------------
program heat_transport_1D

#include <petsc/finclude/petsc.h>

  use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
  use mpp_varpar                      , only : mpp_varpar_init
  use mesh_info                       , only : nz
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscdm
  use petscdmda
  !
  implicit none
  !
  !
#include <petsc/finclude/petsc.h>
  !
  !
  PetscBool          :: converged
  PetscInt           :: converged_reason
  PetscErrorCode     :: ierr
  PetscReal          :: dtime
  PetscInt           :: istep, nstep
  PetscBool          :: flg
  PetscBool          :: save_initial_soln, save_final_soln
  character(len=256) :: string
  character(len=256) :: output_suffix
  PetscViewer        :: viewer

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  ! Set default settings
  nz                = 100
  dtime             = 3600.d0
  nstep             = 2
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
  call init()

  if (save_initial_soln) then
     if (len(trim(adjustl(output_suffix))) ==  0) then
        string = trim(output_suffix) // 'initial_soln.bin'
     else
        string = 'initial_soln_' // trim(output_suffix) // '.bin'
     endif

     call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
     call VecView(thermal_enthalpy_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  do istep = 1, nstep
     ! Update BC
     call set_bondary_conditions()

     ! Run the model
     call thermal_enthalpy_mpp%soe%StepDT(dtime, istep, &
          converged, converged_reason, ierr); CHKERRQ(ierr)
  enddo

  if (save_final_soln) then
     if (len(trim(adjustl(output_suffix))) ==  0) then
        string = trim(output_suffix) // 'final_soln.bin'
     else
        string = 'final_soln_' // trim(output_suffix) // '.bin'
     endif
     call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
     call VecView(thermal_enthalpy_mpp%soe%solver%soln,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end program heat_transport_1D

!------------------------------------------------------------------------
subroutine Init()
  !
  use MultiPhysicsProbThermalEnthalpy, only : thermal_enthalpy_mpp
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
  call thermal_enthalpy_mpp%SetupProblem()

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
#include <petsc/finclude/petsc.h>
  !
  ! !USES:
  use MultiPhysicsProbConstants , only : MPP_THERMAL_EBASED_SNES_CLM
  use MultiPhysicsProbThermalEnthalpy, only : thermal_enthalpy_mpp
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
  call thermal_enthalpy_mpp%Init       ()
  call thermal_enthalpy_mpp%SetName    ('1D heat conduction')
  call thermal_enthalpy_mpp%SetID      (MPP_THERMAL_EBASED_SNES_CLM)
  call thermal_enthalpy_mpp%SetMPIRank (iam)

end subroutine initialize_mpp

!------------------------------------------------------------------------
subroutine add_meshes()
  !
#include <petsc/finclude/petsc.h>
  !
  use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
    use MultiPhysicsProbConstants     , only : CONN_IN_Z_DIR, MESH_CLM_THERMAL_SOIL_COL
    use mpp_varpar                    , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    use MeshType                      , only : mesh_type, MeshCreate
    use mpp_varpar                    , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
    use mesh_info
    !
    implicit none
    !
    class(mesh_type), pointer :: mesh
    PetscInt                  :: imesh
    PetscInt                  :: nlev

    call mpp_varpar_set_nlevsoi(nz)
    call mpp_varpar_set_nlevgrnd(nz)

    imesh        = 1
    nlev         = nz
    ncells_local = nx*ny*nz

    !
    ! Set up the meshes
    !
    call thermal_enthalpy_mpp%SetNumMeshes(1)

    allocate(mesh)

    call MeshCreate(mesh, 'Soil mesh', x_column, y_column, z_column, &
         nx, ny, nz, CONN_IN_Z_DIR)
    call mesh%SetID(MESH_CLM_THERMAL_SOIL_COL)

    call thermal_enthalpy_mpp%AddMesh(imesh, mesh)

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
  use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
  use MultiPhysicsProbConstants , only : GE_THERM_SOIL_EBASED
  use MultiPhysicsProbConstants       , only : MESH_CLM_THERMAL_SOIL_COL
  !
  ! !ARGUMENTS
  implicit none

  call thermal_enthalpy_mpp%AddGovEqn(GE_THERM_SOIL_EBASED, 'Heat transport based on enthalpy ODE', &
       MESH_CLM_THERMAL_SOIL_COL)

  call thermal_enthalpy_mpp%SetMeshesOfGoveqns()
    
end subroutine add_goveqns

!------------------------------------------------------------------------
subroutine add_conditions_to_goveqns()
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
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

  call thermal_enthalpy_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
       'Constant temperature condition at top', 'K', COND_DIRICHLET, &
       SOIL_TOP_CELLS)

  call thermal_enthalpy_mpp%soe%AddConditionInGovEqn(ieqn, COND_BC,   &
       'Constant temperature condition at bottom', 'K', COND_DIRICHLET, &
       SOIL_BOTTOM_CELLS)

end subroutine add_conditions_to_goveqns

!------------------------------------------------------------------------
subroutine allocate_auxvars()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
  !
  implicit none
  
  !
  ! Allocate auxvars
  !
  call thermal_enthalpy_mpp%AllocateAuxVars()

end subroutine allocate_auxvars

!------------------------------------------------------------------------
subroutine set_material_properties()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
  use MultiPhysicsProbThermalEnthalpy , only : MPPThermalSetSoils
  use MultiPhysicsProbConstants       , only : GRAVITY_CONSTANT
  use EOSWaterMod                     , only : DENSITY_TGDPB01, DENSITY_CONSTANT
  use EOSWaterMod                     , only : INT_ENERGY_ENTHALPY_CONSTANT
  use mesh_info                       , only : nz, ncells_local, ncells_ghost
  use mpp_varcon                      , only : denh2o
  use mpp_varcon                      , only : grav
  !
  implicit none
  !
  PetscReal        , parameter :: porosity            = 0.368d0
  PetscReal        , parameter :: lambda              = 0.5d0
  PetscReal        , parameter :: alpha               = 3.4257d-4
  PetscReal        , parameter :: perm                = 8.3913d-12
  !
  PetscReal        , pointer   :: watsat(:,:)
  PetscReal        , pointer   :: hksat(:,:)
  PetscReal        , pointer   :: bsw(:,:)
  PetscReal        , pointer   :: sucsat(:,:)
  PetscReal        , pointer   :: eff_porosity(:,:)
  PetscReal        , pointer   :: residual_sat(:,:)
  PetscInt         , pointer   :: lun_type(:)
  PetscReal        , pointer   :: csol(:,:)
  PetscReal        , pointer   :: tkmg(:,:)
  PetscReal        , pointer   :: tkdry(:,:)
  PetscReal        , parameter :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC
  PetscInt                     :: begc , endc
  integer          , pointer   :: filter(:)
  character(len=32)            :: satfunc_type
  !-----------------------------------------------------------------------

  begc = 1
  endc = 1

  satfunc_type = 'van_genuchten'
  
  allocate(watsat       (1,nz))
  allocate(hksat        (1,nz))
  allocate(bsw          (1,nz))
  allocate(sucsat       (1,nz))
  allocate(eff_porosity (1,nz))
  allocate(residual_sat (1,nz))
  allocate(filter       (nz))
  allocate(lun_type     (nz))
  allocate(csol         (1,nz))
  allocate(tkmg         (1,nz))
  allocate(tkdry        (1,nz))

  ! Soil properties
  filter(:)         = 1
  watsat(:,:)       = porosity
  hksat(:,:)        = perm /vish2o * (denh2o * grav) / 0.001d0
  bsw(:,:)          = 1.d0/lambda
  sucsat(:,:)       = 1.d0/(alpha*GRAVITY_CONSTANT)
  eff_porosity(:,:) = 0.d0
  residual_sat(:,:) = 0.2772d0

  lun_type(:)       = 1
  csol(:,:)         = 837.d0
  tkmg(:,:)         = 0.5d0
  tkdry(:,:)        = 0.25d0
  
  call MPPThermalSetSoils(thermal_enthalpy_mpp, begc, endc, filter, &
       watsat, csol, tkdry,                           &
       hksat, bsw, sucsat, residual_sat,                &
       satfunc_type, DENSITY_CONSTANT, INT_ENERGY_ENTHALPY_CONSTANT)

  deallocate(watsat       )
  deallocate(hksat        )
  deallocate(bsw          )
  deallocate(sucsat       )
  deallocate(eff_porosity )
  deallocate(filter       )
  deallocate(residual_sat )
  deallocate(lun_type     )
  deallocate(csol         )
  deallocate(tkmg         )
  deallocate(tkdry        )

end subroutine set_material_properties

!------------------------------------------------------------------------
subroutine set_initial_conditions()
  !
  ! !DESCRIPTION:
  !
#include <petsc/finclude/petsc.h>
  !
  use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
  use mesh_info                       , only : nz, ncells_local, ncells_ghost
  use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL, VAR_PRESSURE
  use petscsys
  use petscvec
  use petscdm
  !
  implicit none
  !
  PetscErrorCode      :: ierr
  PetscInt            :: nDM
  DM        , pointer :: dms(:)
  Vec       , pointer :: soln_subvecs(:)
  PetscReal , pointer :: temp_p(:)
  PetscInt            :: ii
  PetscViewer         :: viewer
  PetscReal , pointer :: pressure(:)
  PetscInt            :: soe_auxvar_id

  ! Find number of GEs packed within the SoE
  call DMCompositeGetNumberDM(thermal_enthalpy_mpp%soe%solver%dm, nDM, ierr)

  ! Get DMs for each GE
  allocate (dms(nDM))
  call DMCompositeGetEntriesArray(thermal_enthalpy_mpp%soe%solver%dm, dms, ierr)

  ! Allocate vectors for individual GEs
  allocate(soln_subvecs(nDM))

  ! Get solution vectors for individual GEs
  call DMCompositeGetAccessArray(thermal_enthalpy_mpp%soe%solver%dm, thermal_enthalpy_mpp%soe%solver%soln, nDM, &
       PETSC_NULL_INTEGER, soln_subvecs, ierr)

  do ii = 1, nDM
     call VecGetArrayF90(soln_subvecs(ii), temp_p, ierr)
     temp_p(:) = 283.15d0
     call VecRestoreArrayF90(soln_subvecs(ii), temp_p, ierr)
  enddo

  ! Restore solution vectors for individual GEs
  call DMCompositeRestoreAccessArray(thermal_enthalpy_mpp%soe%solver%dm, thermal_enthalpy_mpp%soe%solver%soln, nDM, &
       PETSC_NULL_INTEGER, soln_subvecs, ierr)

  call VecCopy(thermal_enthalpy_mpp%soe%solver%soln, thermal_enthalpy_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
  call VecCopy(thermal_enthalpy_mpp%soe%solver%soln, thermal_enthalpy_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

  allocate(pressure(ncells_local))
  pressure(:) = 091325.d0
  soe_auxvar_id = -1
  call thermal_enthalpy_mpp%soe%SetDataFromCLM(AUXVAR_INTERNAL,  &
       VAR_PRESSURE, soe_auxvar_id, pressure)
  deallocate(pressure)

end subroutine set_initial_conditions

!------------------------------------------------------------------------
subroutine set_bondary_conditions()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbThermalEnthalpy , only : thermal_enthalpy_mpp
  use mesh_info                       , only : nz, ncells_local, ncells_ghost
  use MultiPhysicsProbConstants       , only : AUXVAR_BC, VAR_BC_SS_CONDITION
  use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL, VAR_PRESSURE
  !
  implicit none
  !
  PetscReal, pointer :: top_bc(:)
  PetscReal, pointer :: bot_bc(:)
  PetscReal, pointer :: pressure(:)
  PetscInt           :: soe_auxvar_id

  allocate(top_bc(1))
  allocate(bot_bc(1))
  allocate(pressure(ncells_local))

  top_bc(:) = 303.15d0
  bot_bc(:) = 293.15d0
  pressure(:) = 091325.d0

  soe_auxvar_id = 1
  call thermal_enthalpy_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
       VAR_BC_SS_CONDITION, soe_auxvar_id, top_bc)
  
  soe_auxvar_id = 2
  call thermal_enthalpy_mpp%soe%SetDataFromCLM(AUXVAR_BC,  &
       VAR_BC_SS_CONDITION, soe_auxvar_id, bot_bc)
  
  soe_auxvar_id = -1
  call thermal_enthalpy_mpp%soe%SetDataFromCLM(AUXVAR_INTERNAL,  &
       VAR_PRESSURE, soe_auxvar_id, pressure)

  deallocate(top_bc)
  deallocate(bot_bc)
  deallocate(pressure)

end subroutine set_bondary_conditions
