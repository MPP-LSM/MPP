module mlc_mms_conditions

#include <petsc/finclude/petsc.h>

  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use mlc_mms_global_vars
  use petscsys

  implicit none

  public :: add_conditions_to_goveqns

contains

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS
    !
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp

    call add_non_coupling_conditions_to_goveqns(mlc_mpp)
    call add_internal_coupling_vars(mlc_mpp)

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine add_non_coupling_conditions_to_goveqns(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use ConnectionSetType               , only : connection_set_type
    use MeshType                        , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants       , only : CONN_VERTICAL
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    PetscInt            :: kk , nconn
    PetscInt  , pointer :: id_up(:)
    PetscInt  , pointer :: id_dn(:)
    PetscReal , pointer :: dist_up(:)
    PetscReal , pointer :: dist_dn(:)
    PetscReal , pointer :: area(:)
    PetscInt  , pointer :: itype(:)
    PetscReal , pointer :: unit_vec(:,:)
    PetscReal                  :: dz_cair
    class(connection_set_type) , pointer :: conn_set


    nconn  = ncair

    allocate(id_up    (nconn   ))
    allocate(id_dn    (nconn   ))
    allocate(dist_up  (nconn   ))
    allocate(dist_dn  (nconn   ))
    allocate(area     (nconn   ))
    allocate(itype    (nconn   ))
    allocate(unit_vec (nconn,3 ))

    dz_cair = z_cair/nz_cair;

    do kk = 1, ncair
       id_up(kk)      = 0
       id_dn(kk)      = (nz_cair + 1)*kk
       dist_up(kk)    =  0.d0
       dist_dn(kk)    = dz_cair
       area(kk)       = 1.d0
       unit_vec(kk,1) = -1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_VERTICAL
    end do

    allocate(conn_set)
    call MeshCreateConnectionSet(mlc_mpp%meshes(CAIR_MESH), &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call mlc_mpp%soe%AddConditionInGovEqn(CAIR_TEMP_GE, &
         ss_or_bc_type = COND_BC,   &
         name = 'Atmospheric forcing', &
         unit = 'K', &
         cond_type = COND_DIRICHLET, &
         conn_set = conn_set)

    do kk = 1, ncair
       id_up(kk)      = 0
       id_dn(kk)      = (nz_cair + 1)*kk
       dist_up(kk)    =  0.d0
       dist_dn(kk)    = dz_cair
       area(kk)       = 1.d0
       unit_vec(kk,1) = -1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_VERTICAL
    end do

    allocate(conn_set)
    call MeshCreateConnectionSet(mlc_mpp%meshes(CAIR_MESH), &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call mlc_mpp%soe%AddConditionInGovEqn(CAIR_VAPR_GE, &
         ss_or_bc_type = COND_BC,   &
         name = 'Atmospheric forcing', &
         unit = 'K', &
         cond_type = COND_DIRICHLET, &
         conn_set = conn_set)

  end subroutine add_non_coupling_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine add_internal_coupling_vars(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : ALL_CELLS
    use MultiPhysicsProbConstants , only : SOIL_TOP_CELLS
    use MultiPhysicsProbConstants , only : SOIL_BOTTOM_CELLS
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants     , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants     , only : VAR_WATER_VAPOR
!    use MultiPhysicsProbConstants     , only : VAR_SUNLIT_TEMPERATURE
!    use MultiPhysicsProbConstants     , only : VAR_SHADED_TEMPERATURE
    use MultiPhysicsProbConstants     , only : VAR_LEAF_TEMPERATURE
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    PetscInt                             :: ieqn, num_ieqn_others
    PetscInt                   , pointer :: ieqn_others(:)
    PetscBool                  , pointer :: icoupling_others(:)
    PetscInt                             :: ibc
    PetscInt                             :: jj
    PetscInt                             :: mesh_id
    PetscInt                             :: region_id
    class(connection_set_type) , pointer :: conn_set
    character(len=256)                   :: bc_name
    character(len=256)                   :: eqn_1_string
    character(len=256)                   :: eqn_2_string
    integer                         :: nvars_for_coupling
    integer  , pointer              :: var_ids_for_coupling(:)
    integer  , pointer              :: goveqn_ids_for_coupling(:)
    integer  , pointer              :: is_bc(:)

    nvars_for_coupling = 3
    allocate (var_ids_for_coupling    (nvars_for_coupling))
    allocate (goveqn_ids_for_coupling (nvars_for_coupling))
    allocate (is_bc                   (nvars_for_coupling))

    is_bc(:) = 0 ! Coupling among all GEs is internal

    ! Air Temp <--- Tsun
    !          <--- Tshd
    !          <--- Air vapor
    ieqn = CAIR_TEMP_GE
    goveqn_ids_for_coupling(1) = CLEF_TEMP_SUN_GE; var_ids_for_coupling(1) = VAR_LEAF_TEMPERATURE
    goveqn_ids_for_coupling(2) = CLEF_TEMP_SHD_GE; var_ids_for_coupling(2) = VAR_LEAF_TEMPERATURE
    goveqn_ids_for_coupling(3) = CAIR_VAPR_GE    ; var_ids_for_coupling(3) = VAR_WATER_VAPOR

    call mlc_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

    ! Air Vapor <--- Tsun
    !           <--- Tshd
    ieqn = CAIR_VAPR_GE
    goveqn_ids_for_coupling(1) = CLEF_TEMP_SUN_GE; var_ids_for_coupling(1) = VAR_LEAF_TEMPERATURE
    goveqn_ids_for_coupling(2) = CLEF_TEMP_SHD_GE; var_ids_for_coupling(2) = VAR_LEAF_TEMPERATURE
    goveqn_ids_for_coupling(3) = CAIR_TEMP_GE    ; var_ids_for_coupling(3) = VAR_TEMPERATURE

    call mlc_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

    ! Tsun <--- Air Temp
    !      <--- Air Vapor
    ieqn = CLEF_TEMP_SUN_GE
    nvars_for_coupling = 2
    goveqn_ids_for_coupling(1) = CAIR_TEMP_GE; var_ids_for_coupling(1) = VAR_TEMPERATURE
    goveqn_ids_for_coupling(2) = CAIR_VAPR_GE; var_ids_for_coupling(2) = VAR_WATER_VAPOR

    call mlc_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

    ! Tsun <--- Air Temp
    !      <--- Air Vapor
    ieqn = CLEF_TEMP_SHD_GE
    nvars_for_coupling = 2
    goveqn_ids_for_coupling(1) = CAIR_TEMP_GE; var_ids_for_coupling(1) = VAR_TEMPERATURE
    goveqn_ids_for_coupling(2) = CAIR_VAPR_GE; var_ids_for_coupling(2) = VAR_WATER_VAPOR

    call mlc_mpp%GovEqnSetBothCouplingVars(ieqn, nvars_for_coupling, &
         var_ids_for_coupling, goveqn_ids_for_coupling, is_bc)

  end subroutine add_internal_coupling_vars


end module mlc_mms_conditions
