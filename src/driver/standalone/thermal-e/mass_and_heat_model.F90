
module mesh_info
#include "finclude/petscsys.h"
  PetscInt  , parameter :: nx       = 100
  PetscInt  , parameter :: ny       = 1
  PetscReal , parameter :: x_column = 1.d0
  PetscReal , parameter :: y_column = 1.d0
  PetscReal , parameter :: z_column = 1.d0
  PetscInt              :: nz
  PetscInt              :: ncells_local
  PetscInt              :: ncells_ghost
end module mesh_info

!------------------------------------------------------------------------
program heat_transport_1D_with_advection

  use MultiPhysicsProbTH , only : th_mpp
  use mpp_varpar                      , only : mpp_varpar_init
  use mesh_info                       , only : nz
  !
  implicit none
  !
  !
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscts.h"
#include "finclude/petscts.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscsnes.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"
#include "finclude/petscviewer.h"
  !
  !
  PetscBool          :: converged
  PetscInt           :: converged_reason
  PetscErrorCode     :: ierr
  PetscReal          :: dtime
  PetscInt           :: istep, nstep
  PetscInt           :: flg
  PetscBool          :: save_initial_soln, save_final_soln
  character(len=256) :: string
  character(len=256) :: output_suffix
  PetscViewer        :: viewer

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  ! Set default settings
  nz                = 1
  dtime             = 3600.d0
  nstep             = 1
  save_initial_soln = PETSC_FALSE
  save_final_soln   = PETSC_FALSE
  output_suffix     = ''

  ! Get some command line options

  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-nz',nz,flg,ierr)
  call PetscOptionsGetReal(PETSC_NULL_CHARACTER,'-dt',dtime,flg,ierr)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-nstep',nstep,flg,ierr)
  call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-save_initial_soln',save_initial_soln,flg,ierr)
  call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-save_final_soln',save_final_soln,flg,ierr)
  call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-output_suffix',output_suffix,flg,ierr)

  ! Initialize the problem
  call init()

  if (save_initial_soln) then
     if (len(trim(adjustl(output_suffix))) ==  0) then
        string = trim(output_suffix) // 'initial_soln.bin'
     else
        string = 'initial_soln_' // trim(output_suffix) // '.bin'
     endif

     call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
     call VecView(th_mpp%sysofeqns%soln,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  do istep = 1, nstep
     ! Update BC
     call set_bondary_conditions()

     ! Run the model
     call th_mpp%sysofeqns%StepDT(dtime, istep, &
          converged, converged_reason, ierr); CHKERRQ(ierr)
  enddo

  if (save_final_soln) then
     if (len(trim(adjustl(output_suffix))) ==  0) then
        string = trim(output_suffix) // 'final_soln.bin'
     else
        string = 'final_soln_' // trim(output_suffix) // '.bin'
     endif
     call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(string),FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
     call VecView(th_mpp%sysofeqns%soln,viewer,ierr);CHKERRQ(ierr)
     call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end program heat_transport_1D_with_advection

!------------------------------------------------------------------------
subroutine Init()
  !
  use MultiPhysicsProbTH, only : th_mpp
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
  use MultiPhysicsProbTH        , only : th_mpp
  !
  ! !ARGUMENTS
  implicit none
  !
#include "finclude/petscsys.h"
  !
  PetscInt       :: iam
  PetscErrorCode :: ierr

  call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

  !
  ! Set up the multi-physics problem
  !
  call th_mpp%Init       ()
  call th_mpp%SetName    ('1D heat conduction')
  call th_mpp%SetID      (MPP_TH_SNES_CLM)
  call th_mpp%SetMPIRank (iam)

end subroutine initialize_mpp

!------------------------------------------------------------------------
subroutine add_meshes()
  !
  use MultiPhysicsProbTH , only : th_mpp
  use MultiPhysicsProbConstants       , only : MESH_ALONG_GRAVITY
  use MultiPhysicsProbConstants       , only : MESH_AGAINST_GRAVITY
  use MultiPhysicsProbConstants       , only : MESH_CLM_THERMAL_SOIL_COL
  use MultiPhysicsProbConstants       , only : VAR_XC
  use MultiPhysicsProbConstants       , only : VAR_YC
  use MultiPhysicsProbConstants       , only : VAR_ZC
  use MultiPhysicsProbConstants       , only : VAR_DX
  use MultiPhysicsProbConstants       , only : VAR_DY
  use MultiPhysicsProbConstants       , only : VAR_DZ
  use MultiPhysicsProbConstants       , only : VAR_AREA
  use MultiPhysicsProbConstants       , only : CONN_SET_INTERNAL
  use MultiPhysicsProbConstants       , only : CONN_SET_LATERAL
  use MultiPhysicsProbConstants       , only : CONN_VERTICAL
  use mesh_info                       , only : nx, ny, x_column, y_column, z_column, nz
  use mesh_info                       , only : ncells_local, ncells_ghost
  use mpp_varpar                      , only : mpp_varpar_set_nlevsoi, mpp_varpar_set_nlevgrnd
  !
  implicit none
  !
#include "finclude/petscsys.h"
  !
  PetscReal          :: dx, dy, dz
  PetscInt           :: imesh, kk
  PetscInt           :: nlev
  PetscInt           :: iconn, vert_nconn
  PetscReal, pointer :: soil_xc(:)           ! x-position of grid cell [m]
  PetscReal, pointer :: soil_yc(:)           ! y-position of grid cell [m]
  PetscReal, pointer :: soil_zc(:)           ! z-position of grid cell [m]
  PetscReal, pointer :: soil_dx(:)           ! layer thickness of grid cell [m]
  PetscReal, pointer :: soil_dy(:)           ! layer thickness of grid cell [m]
  PetscReal, pointer :: soil_dz(:)           ! layer thickness of grid cell [m]
  PetscReal, pointer :: soil_area(:)         ! area of grid cell [m^2]
  PetscInt , pointer :: soil_filter(:)       ! 
  
  PetscInt, pointer  :: vert_conn_id_up(:)   !
  PetscInt, pointer  :: vert_conn_id_dn(:)   !
  PetscReal, pointer :: vert_conn_dist_up(:) !
  PetscReal, pointer :: vert_conn_dist_dn(:) !
  PetscReal, pointer :: vert_conn_area(:)    !
  PetscReal, pointer :: vert_conn_type(:)    !

  PetscErrorCode     :: ierr

  call mpp_varpar_set_nlevsoi(nz)
  call mpp_varpar_set_nlevgrnd(nz)

  dx = x_column/nx
  dy = y_column/ny
  dz = z_column/nz

  imesh        = 1
  nlev         = nz
  ncells_local = nx*ny*nz
  ncells_ghost = 0

  allocate(soil_xc(ncells_local))
  allocate(soil_yc(ncells_local))
  allocate(soil_zc(ncells_local))
  allocate(soil_dx(ncells_local))
  allocate(soil_dy(ncells_local))
  allocate(soil_dz(ncells_local))
  allocate(soil_area(ncells_local))
  allocate(soil_filter(ncells_local))
  
  soil_filter (:) = 1
  soil_area   (:) = dx*dy
  soil_dx     (:) = dx
  soil_dy     (:) = dy
  soil_dz     (:) = dz
  soil_yc     (:) = dy/2.d0
  soil_zc     (:) = dz/2.d0

  do kk = 1,nx
     soil_xc(kk) = dx/2.d0 + dx * (kk - 1)
  enddo
  
  !
  ! Set up the meshes
  !    
  call th_mpp%SetNumMeshes(1)

  call th_mpp%MeshSetName        (imesh, 'Soil mesh')
  call th_mpp%MeshSetOrientation (imesh, MESH_AGAINST_GRAVITY)
  call th_mpp%MeshSetID          (imesh, MESH_CLM_THERMAL_SOIL_COL)
  call th_mpp%MeshSetDimensions  (imesh, ncells_local, ncells_ghost, nlev)

  call th_mpp%MeshSetGridCellFilter      (imesh, soil_filter)
  call th_mpp%MeshSetGeometricAttributes (imesh, VAR_XC   , soil_xc)
  call th_mpp%MeshSetGeometricAttributes (imesh, VAR_YC   , soil_yc)
  call th_mpp%MeshSetGeometricAttributes (imesh, VAR_ZC   , soil_zc)
  call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DX   , soil_dx)
  call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DY   , soil_dy)
  call th_mpp%MeshSetGeometricAttributes (imesh, VAR_DZ   , soil_dz)
  call th_mpp%MeshSetGeometricAttributes (imesh, VAR_AREA , soil_area)
  call th_mpp%MeshComputeVolume          (imesh)

  allocate (vert_conn_id_up   (nx-1))
  allocate (vert_conn_id_dn   (nx-1))
  allocate (vert_conn_dist_up (nx-1))
  allocate (vert_conn_dist_dn (nx-1))
  allocate (vert_conn_area    (nx-1))
  allocate (vert_conn_type    (nx-1))

  iconn = 0
  do kk = 1, nx-1

     iconn = iconn + 1
     vert_conn_id_up(iconn)   = kk
     vert_conn_id_dn(iconn)   = vert_conn_id_up(iconn) + 1
     vert_conn_dist_up(iconn) = 0.5d0*dx
     vert_conn_dist_dn(iconn) = 0.5d0*dx
     vert_conn_area(iconn)    = dy*dz
     vert_conn_type(iconn)    = CONN_VERTICAL

  end do
  vert_nconn = iconn

  call th_mpp%MeshSetConnectionSet(imesh, CONN_SET_INTERNAL, &
       vert_nconn,  vert_conn_id_up, vert_conn_id_dn,                      &
       vert_conn_dist_up, vert_conn_dist_dn,  vert_conn_area)
  
  deallocate(soil_xc)
  deallocate(soil_yc)
  deallocate(soil_zc)
  deallocate(soil_dx)
  deallocate(soil_dy)
  deallocate(soil_dz)
  deallocate(soil_area)
  deallocate(soil_filter)  

end subroutine add_meshes

  !------------------------------------------------------------------------

subroutine add_goveqns()
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  use MultiPhysicsProbTH , only : th_mpp
  use MultiPhysicsProbConstants       , only : GE_THERM_SOIL_EBASED
  use MultiPhysicsProbConstants       , only : GE_RE
  !
  ! !ARGUMENTS
  implicit none

  call th_mpp%AddGovEqn(GE_RE, 'Mass equation')
  call th_mpp%AddGovEqn(GE_THERM_SOIL_EBASED, 'Heat transport based on enthalpy')

  call th_mpp%SetMeshesOfGoveqns()
    
end subroutine add_goveqns

!------------------------------------------------------------------------
subroutine add_conditions_to_goveqns()
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
  use MultiPhysicsProbTH , only : th_mpp
  use MultiPhysicsProbConstants       , only : ALL_CELLS
  use MultiPhysicsProbConstants       , only : SOIL_CELLS
  use MultiPhysicsProbConstants       , only : SOIL_TOP_CELLS
  use MultiPhysicsProbConstants       , only : SOIL_BOTTOM_CELLS
  use MultiPhysicsProbConstants       , only : COND_BC
  use MultiPhysicsProbConstants       , only : COND_DIRICHLET
  use mesh_info                       , only : nx, ny, x_column, y_column, z_column, nz
  use ConnectionSetType               , only : connection_set_type
  use MeshType                        , only : MeshCreateConnectionSet
  !
  ! !ARGUMENTS
  implicit none
  !
  PetscInt                            :: ieqn
  PetscInt                            :: ieqn_1, ieqn_2
  PetscInt                            :: nconn
  PetscReal                           :: dx, dy, dz
  PetscInt                  , pointer :: id_up(:)
  PetscInt                  , pointer :: id_dn(:)
  PetscReal                 , pointer :: dist_up(:)
  PetscReal                 , pointer :: dist_dn(:)
  PetscReal                 , pointer :: area(:)
  PetscReal                 , pointer :: unit_vec(:,:)
  type(connection_set_type) , pointer :: conn_set
  
!  ieqn_1 = 1
!  ieqn_2 = 2
!  call th_mpp%GovEqnAddCouplingCondition(ieqn_1, ieqn_2, &
!       ALL_CELLS, ALL_CELLS)
  
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  nconn         = 1

  allocate(id_up    (nconn   ))
  allocate(id_dn    (nconn   ))
  allocate(dist_up  (nconn   ))
  allocate(dist_dn  (nconn   ))
  allocate(area     (nconn   ))
  allocate(unit_vec (nconn,3 ))

  dx = x_column/nx
  dy = y_column/ny
  dz = z_column/nz

  id_up(1)      = 0
  id_dn(1)      = 1
  dist_up(1)    = 0.d0
  dist_dn(1)    = 0.5d0*dx
  area(1)       = dy*dz
  unit_vec(1,1) = 1.d0
  unit_vec(1,2) = 0.d0
  unit_vec(1,3) = 0.d0

  allocate(conn_set)
  call MeshCreateConnectionSet(th_mpp%meshes(1), &
       -1, nconn, id_up, id_dn, &
       dist_up, dist_dn, area, unit_vec, conn_set)

  ieqn = 2

  call th_mpp%GovEqnAddCondition(ieqn, COND_BC,   &
       'Constant temperature condition at top', 'K', COND_DIRICHLET, &
       SOIL_TOP_CELLS, conn_set=conn_set)

  id_up(1)      = 0
  id_dn(1)      = nx
  dist_up(1)    = 0.d0
  dist_dn(1)    = 0.5d0*dx
  area(1)       = dy*dz
  unit_vec(1,1) = -1.d0
  unit_vec(1,2) = 0.d0
  unit_vec(1,3) = 0.d0

  allocate(conn_set)
  call MeshCreateConnectionSet(th_mpp%meshes(1), &
       -1, nconn, id_up, id_dn, &
       dist_up, dist_dn, area, unit_vec, conn_set)

  call th_mpp%GovEqnAddCondition(ieqn, COND_BC,   &
       'Constant temperature condition at bottom', 'K', COND_DIRICHLET, &
       SOIL_BOTTOM_CELLS, conn_set=conn_set)

  deallocate(id_up   )
  deallocate(id_dn   )
  deallocate(dist_up )
  deallocate(dist_dn )
  deallocate(area    )
  deallocate(unit_vec)

end subroutine add_conditions_to_goveqns

!------------------------------------------------------------------------
subroutine allocate_auxvars()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbTH , only : th_mpp
  use MultiPhysicsProbConstants , only : VAR_PRESSURE
  use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
  !
  implicit none
  !
  integer           :: ieqn
  integer           :: nvars_for_coupling
  integer, pointer  :: var_ids_for_coupling(:)
  integer, pointer  :: goveqn_ids_for_coupling(:)
  integer, pointer  :: is_bc(:)
  
  !
  ! Allocate auxvars
  !
  call th_mpp%AllocateAuxVars()

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

  call th_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
       var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)
    
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

  call th_mpp%GovEqnSetCouplingVars(ieqn, nvars_for_coupling, &
       var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)
    
  deallocate(var_ids_for_coupling   )
  deallocate(goveqn_ids_for_coupling)

end subroutine allocate_auxvars

!------------------------------------------------------------------------
subroutine set_material_properties()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbTH        , only : th_mpp
  use MultiPhysicsProbTH        , only : MPPTHSetSoils
  use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
  use EOSWaterMod               , only : DENSITY_TGDPB01, DENSITY_CONSTANT, DENSITY_IFC67
  use EOSWaterMod               , only : INT_ENERGY_ENTHALPY_IFC67
  use mesh_info                 , only : nx,ny, nz, ncells_local, ncells_ghost
  use mpp_varcon                , only : denh2o
  use mpp_varcon                , only : grav
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
  endc = nx*2

  satfunc_type = 'van_genuchten'
  
  allocate(watsat       (nx*ny*2,nz))
  allocate(hksat        (nx*ny*2,nz))
  allocate(bsw          (nx*ny*2,nz))
  allocate(sucsat       (nx*ny*2,nz))
  allocate(eff_porosity (nx*ny*2,nz))
  allocate(residual_sat (nx*ny*2,nz))
  allocate(filter       (ncells_local*2))
  allocate(lun_type     (ncells_local*2))
  allocate(csol         (nx*ny*2,nz))
  allocate(tkmg         (nx*ny*2,nz))
  allocate(tkdry        (nx*ny*2,nz))

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
  
  call MPPTHSetSoils(th_mpp, begc, endc, filter, &
       lun_type, watsat, csol, tkmg, tkdry,                           &
       hksat, bsw, sucsat, eff_porosity, residual_sat,                &
       satfunc_type, DENSITY_IFC67, INT_ENERGY_ENTHALPY_IFC67)

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
  use MultiPhysicsProbTH , only : th_mpp
  use mesh_info                       , only : nz, ncells_local, ncells_ghost
  use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL, VAR_PRESSURE
  !
  implicit none
  !
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"
  !
  PetscErrorCode      :: ierr
  PetscInt            :: nDM
  DM        , pointer :: dms(:)
  Vec       , pointer :: soln_subvecs(:)
  PetscReal , pointer :: v_p(:)
  PetscInt            :: ii
  PetscViewer         :: viewer
  Vec                 :: pressure_ic
  PetscInt            :: soe_auxvar_id

  ! Find number of GEs packed within the SoE
  call DMCompositeGetNumberDM(th_mpp%sysofeqns%dm, nDM, ierr)

  ! Get DMs for each GE
  allocate (dms(nDM))
  call DMCompositeGetEntriesArray(th_mpp%sysofeqns%dm, dms, ierr)

  ! Allocate vectors for individual GEs
  allocate(soln_subvecs(nDM))

  ! Get solution vectors for individual GEs
  call DMCompositeGetAccessArray(th_mpp%sysofeqns%dm, &
       th_mpp%sysofeqns%soln, nDM, &
       PETSC_NULL_INTEGER, soln_subvecs, ierr)

  do ii = 1, nDM
     call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)
     if (ii == 1) then
        v_p(:) = 091325.d0
     else
        v_p(:) = 283.15d0
     endif
     call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
  enddo

  ! Restore solution vectors for individual GEs
  call DMCompositeRestoreAccessArray(th_mpp%sysofeqns%dm, &
       th_mpp%sysofeqns%soln, nDM, &
       PETSC_NULL_INTEGER, soln_subvecs, ierr)

  call VecCopy(th_mpp%sysofeqns%soln, &
       th_mpp%sysofeqns%soln_prev, ierr); CHKERRQ(ierr)
  call VecCopy(th_mpp%sysofeqns%soln, &
       th_mpp%sysofeqns%soln_prev_clm, ierr); CHKERRQ(ierr)

end subroutine set_initial_conditions

!------------------------------------------------------------------------
subroutine set_bondary_conditions()
  !
  ! !DESCRIPTION:
  !
  use MultiPhysicsProbTH , only : th_mpp
  use mesh_info                       , only : nz, ncells_local, ncells_ghost
  use MultiPhysicsProbConstants       , only : AUXVAR_BC, VAR_BC_SS_CONDITION
  !
  implicit none
  !
  PetscReal, pointer :: top_bc(:)
  PetscReal, pointer :: bot_bc(:)
  PetscInt           :: soe_auxvar_id

  allocate(top_bc(1))
  allocate(bot_bc(1))

  top_bc(:) = 303.15d0
  bot_bc(:) = 293.15d0

  soe_auxvar_id = 1
  call th_mpp%sysofeqns%SetDataFromCLM(AUXVAR_BC,  &
       VAR_BC_SS_CONDITION, soe_auxvar_id, top_bc)
  
  soe_auxvar_id = 2
  call th_mpp%sysofeqns%SetDataFromCLM(AUXVAR_BC,  &
       VAR_BC_SS_CONDITION, soe_auxvar_id, bot_bc)
  
  deallocate(top_bc)
  deallocate(bot_bc)

end subroutine set_bondary_conditions