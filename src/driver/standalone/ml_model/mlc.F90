module mlc

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use ml_model_utils      , only : get_value_from_condition, set_value_in_condition
  use ml_model_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none

#include <petsc/finclude/petsc.h>

  public :: init_mlc
  public :: mlc_set_initial_conditions
  public :: mlc_initialize_from_checkpoint
  public :: extract_data_from_mlc
  public :: checkpoint_mlc

contains

  !------------------------------------------------------------------------
  subroutine initialize(mlc_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_MLC_KSP
    !
    implicit none

    class(mpp_mlc_type) :: mlc_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call mlc_mpp%Init       ()
    call mlc_mpp%SetName    ('MLC radiation model')
    call mlc_mpp%SetID      (MPP_MLC_KSP)
    call mlc_mpp%SetMPIRank (iam)

  end subroutine initialize

  !------------------------------------------------------------------------
  subroutine setup_meshes(mlc_mpp)
    !
    use MeshType             , only : mesh_type
    use ml_model_meshes      , only : create_canopy_airspace_mesh
    use ml_model_meshes      , only : create_canopy_mesh
    use ml_model_global_vars , only : CAIR_MESH, CLEF_MESH
    !
    implicit none

    class(mpp_mlc_type) :: mlc_mpp

    class(mesh_type), pointer :: mesh

    call mlc_mpp%SetNumMeshes(2)

    CAIR_MESH = 1; call create_canopy_airspace_mesh(mesh); call mlc_mpp%AddMesh(CAIR_MESH, mesh)
    CLEF_MESH = 2; call create_canopy_airspace_mesh(mesh); call mlc_mpp%AddMesh(CLEF_MESH, mesh)

    deallocate(mesh)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine setup_leaf2cair_map(mlc_mpp)
    !
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType , only : goveqn_cair_temp_type
    use GoveqnCanopyAirVaporType       , only : goveqn_cair_vapor_type
    use GoveqnCanopyLeafTemperatureType, only : goveqn_cleaf_temp_type
    !
    implicit none
    !
    class(mpp_mlc_type)          :: mlc_mpp
    !
    class(goveqn_base_type) , pointer :: goveq
    PetscInt                , pointer :: leaf2cair(:)
    PetscInt                          :: nleaf2cair, icair, itree, k, count, ieqn

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

  end subroutine setup_leaf2cair_map

  !------------------------------------------------------------------------
  subroutine add_goveqns(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_CANOPY_AIR_TEMP
    use MultiPhysicsProbConstants , only : GE_CANOPY_AIR_VAPOR
    use MultiPhysicsProbConstants , only : GE_CANOPY_LEAF_TEMP
    use ml_model_global_vars      , only : CAIR_TEMP_GE, CAIR_VAPR_GE
    use ml_model_global_vars      , only : CLEF_TEMP_SUN_GE, CLEF_TEMP_SUN_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp

    CAIR_TEMP_GE     = 1; call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_AIR_TEMP , 'Canopy air temperature' , CAIR_MESH)
    CAIR_VAPR_GE     = 2; call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_AIR_VAPOR, 'Canopy air vapor'       , CAIR_MESH)
    CLEF_TEMP_SUN_GE = 3; call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_LEAF_TEMP, 'Sunlit canopy'          , CLEF_MESH)
    CLEF_TEMP_SHD_GE = 4; call mlc_mpp%AddGovEqnWithMeshRank(GE_CANOPY_LEAF_TEMP, 'Shaded canopy'          , CLEF_MESH)

    call mlc_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine init_mlc(mlc_mpp)
    !
    use mlc_conditions, only : mlc_add_conditions_to_goveqns
    use mlc_parameters, only : mlc_set_parameters
    !
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp

    call initialize(mlc_mpp)

    call setup_meshes(mlc_mpp)

    call add_goveqns(mlc_mpp)

    call setup_leaf2cair_map(mlc_mpp)

    call mlc_add_conditions_to_goveqns(mlc_mpp)

    call mlc_mpp%AllocateAuxVars(ncair)

    call mlc_mpp%SetupProblem()

    call mlc_set_parameters(mlc_mpp)

  end subroutine init_mlc

  !------------------------------------------------------------------------
  subroutine mlc_set_initial_conditions(mlc_mpp)
    !
    use MultiPhysicsProbConstants , only : MM_H2O, MM_DRY_AIR
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsMLCType  , only : sysofeqns_mlc_type
    use ml_model_utils            , only : get_value_from_condition
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscReal                            :: relhum, eref, esat, desatdt
    PetscInt                             :: icair
    PetscInt                             :: nDM
    DM                         , pointer :: dms(:)
    Vec                        , pointer :: soln_subvecs(:)
    PetscReal                  , pointer :: v_p(:)
    PetscInt                             :: ii
    PetscViewer                          :: viewer
    PetscInt                             :: soe_auxvar_id
    PetscErrorCode                       :: ierr
    PetscReal :: qref_value, factor

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do icair = 1, ncair
       soe%cturb%pref(icair) = get_value_from_condition(bnd_cond%pref, icair)
       soe%cturb%uref(icair) = get_value_from_condition(bnd_cond%uref, icair)
       soe%cturb%tref(icair) = get_value_from_condition(bnd_cond%tref, icair)
       soe%cturb%rhref(icair)= 80.d0    !get_value_from_condition(Rhref, icair)

       soe%cturb%wind(icair,:) = get_value_from_condition(bnd_cond%uref, icair)

       qref_value = get_value_from_condition(bnd_cond%qref, icair)
       factor     = 1.d0/ (MM_H2O / MM_DRY_AIR + (1.d0 - MM_H2O / MM_DRY_AIR)*qref_value)
       soe%cturb%qref(icair) = qref_value! * factor
       soe%cturb%qcan(icair) = qref_value! * factor
       
       call soe%cturb%ComputeDerivedAtmInputs(icair)

       soe%cturb%tcan(icair) = soe%cturb%tref(icair)
       !soe%cturb%qref(icair) = qref_value * factor
    end do

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

    icair = 1;
    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)

       if (ii == CAIR_TEMP_GE .or. ii == CLEF_TEMP_SUN_GE .or. ii == CLEF_TEMP_SHD_GE) then
          v_p(:) = get_value_from_condition(bnd_cond%tref, icair)

       else if (ii == CAIR_VAPR_GE) then
          qref_value = get_value_from_condition(bnd_cond%qref, icair)
          factor     = 1.d0/ (MM_H2O / MM_DRY_AIR + (1.d0 - MM_H2O / MM_DRY_AIR)*qref_value)
          v_p(:) = qref_value * factor
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

    call mlc_mpp%soe%PreSolve()

  end subroutine mlc_set_initial_conditions

  !------------------------------------------------------------------------
  subroutine mlc_initialize_from_checkpoint(mlc_mpp)
    !
    use MultiPhysicsProbConstants , only : MM_H2O, MM_DRY_AIR
    use ml_model_global_vars      , only : ic_file
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsMLCType       , only : sysofeqns_mlc_type
    use ml_model_utils                 , only : get_value_from_condition
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnCanopyAirVaporType       , only : goveqn_cair_vapor_type
    use GoveqnCanopyAirTemperatureType , only : goveqn_cair_temp_type
    use ml_model_global_vars           , only : nz_cair
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscReal                            :: relhum, eref, esat, desatdt
    PetscInt                             :: icair
    PetscInt                             :: nDM
    DM                         , pointer :: dms(:)
    Vec                        , pointer :: soln_subvecs(:)
    PetscReal                  , pointer :: v_p(:), ic_p(:)
    Vec                                  :: ic_data
    PetscInt                             :: ii, k, offset_ic
    PetscViewer                          :: viewer
    PetscInt                             :: soe_auxvar_id
    PetscErrorCode                       :: ierr
    PetscReal                            :: qref_value, factor

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, ic_file, FILE_MODE_READ, viewer, ierr);CHKERRQ(ierr)
    call VecCreate(PETSC_COMM_WORLD, ic_data, ierr);CHKERRQ(ierr)
    call VecLoad(ic_data, viewer, ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)

    call VecGetArrayF90(ic_data, ic_p, ierr); CHKERRA(ierr)

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    call base_soe%SetPointerToIthGovEqn(CAIR_TEMP_GE, cur_goveq)
    select type(cur_goveq)
    class is (goveqn_cair_temp_type)
       cur_goveq%aux_vars_in(ntop)%temperature = get_value_from_condition(int_cond%Tair, ntop)
    end select

    call base_soe%SetPointerToIthGovEqn(CAIR_VAPR_GE, cur_goveq)
    select type(cur_goveq)
    class is (goveqn_cair_vapor_type)
       cur_goveq%aux_vars_in(ntop)%qair = get_value_from_condition(int_cond%qair, ntop)
    end select

    do icair = 1, ncair
       soe%cturb%pref(icair) = get_value_from_condition(bnd_cond%pref, icair)
       soe%cturb%uref(icair) = get_value_from_condition(bnd_cond%uref, icair)
       soe%cturb%tref(icair) = get_value_from_condition(bnd_cond%tref, icair)
       soe%cturb%rhref(icair)= 80.d0

       offset_ic = 0
       do k = 1, nz_cair + 1
          soe%cturb%wind(icair,k) = ic_p(offset_ic + k)
       end do

       qref_value = get_value_from_condition(bnd_cond%qref, icair)
       soe%cturb%qref(icair) = qref_value
       soe%cturb%qcan(icair) = get_value_from_condition(int_cond%qair, ntop)

       call soe%cturb%ComputeDerivedAtmInputs(icair)

       soe%cturb%tcan(icair) = get_value_from_condition(int_cond%tair, ntop)
    end do

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

    icair = 1;
    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)

       if (ii == CAIR_TEMP_GE) then
          offset_ic = (nz_cair + 1)*1
          do k = 1, nz_cair + 1
             v_p(k) = ic_p(offset_ic + k)
          end do
          call set_value_in_condition(bnd_cond%tg, icair, v_p(1))
       else if (ii == CAIR_VAPR_GE) then
          offset_ic = (nz_cair + 1)*2
          do k = 1, nz_cair + 1
             v_p(k) = ic_p(offset_ic + k)
          end do
       else if (ii == CLEF_TEMP_SUN_GE) then
          offset_ic = (nz_cair + 1)*3
          do k = 1, nz_cair + 1
             v_p(k) = ic_p(offset_ic + k)
          end do
       else if (ii == CLEF_TEMP_SHD_GE) then
          offset_ic = (nz_cair + 1)*4
          do k = 1, nz_cair + 1
             v_p(k) = ic_p(offset_ic + k)
          end do
       endif
       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(mlc_mpp%soe%solver%dm, &
         mlc_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)

    call VecRestoreArrayF90(ic_data, ic_p, ierr); CHKERRA(ierr)

    call VecCopy(mlc_mpp%soe%solver%soln, mlc_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(mlc_mpp%soe%solver%soln, mlc_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    call mlc_mpp%soe%PreSolve()

  end subroutine mlc_initialize_from_checkpoint

 !------------------------------------------------------------------------
  subroutine set_boundary_conditions(mlc_mpp, istep, isubstep)
    !
    use MultiPhysicsProbConstants      , only : MM_H2O, MM_DRY_AIR, RGAS
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsMLCType       , only : sysofeqns_mlc_type
    use ml_model_utils                 , only : get_value_from_condition
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnCanopyAirVaporType       , only : goveqn_cair_vapor_type
    use GoveqnCanopyAirTemperatureType , only : goveqn_cair_temp_type
    use ml_model_utils                 , only : accumulate_data
    use ml_model_global_vars           , only : canp_lev_vars
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    PetscInt :: istep, isubstep
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: icair
    PetscReal                            :: factor, qref_value, tcan_value, qcan_value, eair
    PetscReal                            :: rhomol

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    call base_soe%SetPointerToIthGovEqn(CAIR_TEMP_GE, cur_goveq)
    select type(cur_goveq)
    class is (goveqn_cair_temp_type)
       tcan_value = cur_goveq%aux_vars_in(ntop)%temperature
    end select

    call base_soe%SetPointerToIthGovEqn(CAIR_VAPR_GE, cur_goveq)
    select type(cur_goveq)
    class is (goveqn_cair_vapor_type)
       qcan_value = cur_goveq%aux_vars_in(ntop)%qair
    end select


    do icair = 1, ncair
       soe%cturb%pref(icair) = get_value_from_condition(bnd_cond%pref, icair)
       soe%cturb%uref(icair) = get_value_from_condition(bnd_cond%uref, icair)
       soe%cturb%tref(icair) = get_value_from_condition(bnd_cond%tref, icair)

       qref_value = get_value_from_condition(bnd_cond%qref, icair)
       factor     = 1.d0/ (MM_H2O / MM_DRY_AIR + (1.d0 - MM_H2O / MM_DRY_AIR)*qref_value)
       soe%cturb%qref(icair) = qref_value !* factor

       call soe%cturb%ComputeDerivedAtmInputs(icair)

       if (istep == 1 .and. isubstep == 1) then
          soe%cturb%qcan(icair) = qref_value! * factor
          soe%cturb%tcan(icair) = soe%cturb%tref(icair)
       else
          factor = 1.d0/ (MM_H2O / MM_DRY_AIR + (1.d0 - MM_H2O / MM_DRY_AIR)*qcan_value)

          eair = qcan_value * soe%cturb%pref(icair)
          factor = (MM_H2O / MM_DRY_AIR) / ( soe%cturb%pref(icair) - (1.d0 - MM_H2O / MM_DRY_AIR)*eair )
          soe%cturb%qcan(icair) = eair * factor
          soe%cturb%tcan(icair) = tcan_value
       end if
       soe%cturb%soil_temperature(icair) = get_value_from_condition(bnd_cond%soil_t,icair)

       soe%cturb%soil_rn(icair) = &
            get_value_from_condition(int_cond%Isoil_vis, icair) + &
            get_value_from_condition(int_cond%Isoil_nir, icair) + &
            get_value_from_condition(int_cond%Labs_soil, icair)

       soe%cturb%soil_rhg(icair) = get_value_from_condition(bnd_cond%rhg, icair)
       soe%cturb%soil_res(icair) = get_value_from_condition(bnd_cond%soilres, icair)
       soe%cturb%soil_tk(icair) = get_value_from_condition(bnd_cond%soil_tk, icair)

       call accumulate_data(canp_lev_vars%rnabs_soi, soe%cturb%soil_rn(icair), icair, isubstep)

       rhomol = soe%cturb%pref(icair) / (RGAS * soe%cturb%tref(icair))
       call accumulate_data(canp_lev_vars%gac0_soi, rhomol/soe%cturb%soil_res(icair), icair, isubstep)

    end do

    call set_air_temp_ge_parameters(mlc_mpp)
    call set_air_vapor_ge_parameters(mlc_mpp)
    call set_canopy_leaf_parameters(mlc_mpp, CLEF_TEMP_SUN_GE, istep, isubstep)
    call set_canopy_leaf_parameters(mlc_mpp, CLEF_TEMP_SHD_GE, istep, isubstep)

  end subroutine set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_air_temp_ge_parameters(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsMLCType       , only : sysofeqns_mlc_type
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType , only : goveqn_cair_temp_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, ileaf, icair, icell, gb_count, offset

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CAIR_TEMP_GE, cur_goveq)

    offset = ncair * ntree * (ntop - nbot + 1)

    select type(cur_goveq)
    class is (goveqn_cair_temp_type)

       if (ntree > 1) then
          write(*,*)'set_air_temp_ge_parameters: Need to extend the code to support ntree > 1'
       end if

       gb_count = 0
       do icair = 1, ncair
         do k = 1, nz_cair + 1
            if (k >= nbot .and. k<= ntop) then
               gb_count = gb_count + 1
               icell = (icair-1)*(nz_cair+1) + k

               cur_goveq%aux_vars_in(icell)%gbh(1)  = get_value_from_condition(int_cond%gbh, gb_count         )
               cur_goveq%aux_vars_in(icell)%gbh(2)  = get_value_from_condition(int_cond%gbh, gb_count + offset)

               do ileaf = 1, ntree
                  cur_goveq%aux_vars_in(icell)%leaf_gs(        ileaf) = get_value_from_condition(int_cond%gs_sun, gb_count)
                  cur_goveq%aux_vars_in(icell)%leaf_gs(ntree + ileaf) = get_value_from_condition(int_cond%gs_shd, gb_count)
                  cur_goveq%aux_vars_in(icell)%leaf_fssh(        ileaf) = fssh(k)
                  cur_goveq%aux_vars_in(icell)%leaf_fssh(ntree + ileaf) = 1.d0 - fssh(k)
               enddo
            end if
         end do ! k-loop

      enddo

    end select

  end subroutine set_air_temp_ge_parameters

  !------------------------------------------------------------------------
  subroutine set_air_vapor_ge_parameters(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyAirVaporType  , only : goveqn_cair_vapor_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, ileaf, icair, icell, gb_count, offset

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CAIR_VAPR_GE, cur_goveq)

    offset = ncair * ntree * (ntop - nbot + 1)

    select type(cur_goveq)
    class is (goveqn_cair_vapor_type)

       if (ntree > 1) then
          write(*,*)'set_air_vapor_ge_parameters: Need to extend the code to support ntree > 1'
       end if

       gb_count = 0
       do icair = 1, ncair
         do k = 1, nz_cair + 1
            if (k >= nbot .and. k <= ntop) then

               gb_count = gb_count + 1
               icell = (icair-1)*(nz_cair+1) + k

               cur_goveq%aux_vars_in(icell)%gbv(1)  = get_value_from_condition(int_cond%gbv, gb_count         )
               cur_goveq%aux_vars_in(icell)%gbv(2)  = get_value_from_condition(int_cond%gbv, gb_count + offset)

               do ileaf = 1, ntree
                  cur_goveq%aux_vars_in(icell)%leaf_gs(        ileaf) = get_value_from_condition(int_cond%gs_sun, gb_count)
                  cur_goveq%aux_vars_in(icell)%leaf_gs(ntree + ileaf) = get_value_from_condition(int_cond%gs_shd, gb_count)
                  cur_goveq%aux_vars_in(icell)%leaf_fssh(        ileaf) = fssh(k)
                  cur_goveq%aux_vars_in(icell)%leaf_fssh(ntree + ileaf) = 1.d0 - fssh(k)
               end do
            end if
         end do

      end do
   end select

  end subroutine set_air_vapor_ge_parameters

  !------------------------------------------------------------------------
  subroutine set_canopy_leaf_parameters(mlc_mpp, ge_rank, istep, isubstep)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
    use ml_model_utils                  , only : accumulate_data
    use ml_model_global_vars            , only : vert_lev_vars
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp
    PetscInt           :: ge_rank
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, i, num_int, icell, icair, itree, offset, count
    character(len=20)                    :: step_string, substep_string
    PetscInt :: istep, isubstep

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(ge_rank, cur_goveq)

    if (output_data .and. ge_rank == CLEF_TEMP_SUN_GE) then
       write(step_string,*)istep
       write(substep_string,*)isubstep
       write(*,*)'mpp.rn{' // trim(adjustl(step_string)) // ',' //trim(adjustl(substep_string)) // '} = ['
    end if
    select type(cur_goveq)
    class is (goveqn_cleaf_temp_type)

       offset = ncair * ntree * (ntop - nbot + 1)
       
       count = 0
       icell = 0
       do icair = 1, ncair
          do itree = 1, ntree
             do k = 1, nz_cair + 1

                icell = icell + 1
                if (k >= nbot .and. k <= ntop) then

                   count = count + 1
                   if (cur_goveq%rank_in_soe_list == CLEF_TEMP_SUN_GE) then
                      cur_goveq%aux_vars_in(icell)%gbh  = get_value_from_condition(int_cond%gbh, count)
                      cur_goveq%aux_vars_in(icell)%gbv  = get_value_from_condition(int_cond%gbv, count)
                      cur_goveq%aux_vars_in(icell)%gs   = get_value_from_condition(int_cond%gs_sun  , count)
                      cur_goveq%aux_vars_in(icell)%fssh = fssh(k)

                      cur_goveq%aux_vars_in(icell)%rn   = &
                           get_value_from_condition(int_cond%Ileaf_sun_vis, count) + &
                           get_value_from_condition(int_cond%Ileaf_sun_nir, count) + &
                           get_value_from_condition(int_cond%Labs_leaf_sun       , count)

                      call accumulate_data(vert_lev_vars%rn_leaf_sun, &
                           cur_goveq%aux_vars_in(icell)%rn, count, isubstep)

                      if (output_data) write(*,*)icell,cur_goveq%aux_vars_in(icell)%rn, &
                           get_value_from_condition(int_cond%Ileaf_shd_vis, count), &
                           get_value_from_condition(int_cond%Ileaf_shd_nir, count), &
                           get_value_from_condition(int_cond%Labs_leaf_shd       , count)
                   else
                      cur_goveq%aux_vars_in(icell)%gbh  = get_value_from_condition(int_cond%gbh, count + offset)
                      cur_goveq%aux_vars_in(icell)%gbv  = get_value_from_condition(int_cond%gbv, count + offset)
                      cur_goveq%aux_vars_in(icell)%gs   = get_value_from_condition(int_cond%gs_shd, count)
                      cur_goveq%aux_vars_in(icell)%fssh = 1.d0 - fssh(k)

                      cur_goveq%aux_vars_in(icell)%rn   = &
                           get_value_from_condition(int_cond%Ileaf_shd_vis, count) + &
                           get_value_from_condition(int_cond%Ileaf_shd_nir, count) + &
                           get_value_from_condition(int_cond%Labs_leaf_shd       , count)

                      call accumulate_data(vert_lev_vars%rn_leaf_shd, &
                           cur_goveq%aux_vars_in(icell)%rn, count, isubstep)

                      if (output_data) write(*,*)icell,cur_goveq%aux_vars_in(icell)%rn, &
                           get_value_from_condition(int_cond%Ileaf_shd_vis, count), &
                           get_value_from_condition(int_cond%Ileaf_shd_nir, count), &
                           get_value_from_condition(int_cond%Labs_leaf_shd       , count)
                   end if
                end if

             end do
          end do
       end do
       if (output_data .and. ge_rank == CLEF_TEMP_SHD_GE) write(*,*)'];'

    end select

  end subroutine set_canopy_leaf_parameters

 !------------------------------------------------------------------------
  subroutine get_data_from_mlc_eqn(mlc_mpp, ieqn, var_id, ncells, data)
   !
   ! !DESCRIPTION:
   !
   ! !USES:
   use SystemOfEquationsBaseType       , only : sysofeqns_base_type
   use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
   use MultiPhysicsProbMLC             , only : mpp_mlc_type
   use GoverningEquationBaseType       , only : goveqn_base_type
   use GoveqnCanopyAirTemperatureType  , only : goveqn_cair_temp_type
   use GoveqnCanopyAirVaporType        , only : goveqn_cair_vapor_type
   use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
   use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    !
    ! !ARGUMENTS
   implicit none
   !
   class(mpp_mlc_type)  :: mlc_mpp
   PetscInt             :: ieqn
   PetscInt             :: var_id
   PetscInt             :: ncells
   PetscReal,  pointer  :: data(:)
   !
   class(goveqn_base_type) , pointer :: goveq
   PetscErrorCode       :: ierr

   call mlc_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

   select type(goveq)
   class is (goveqn_cair_temp_type)
      call goveq%GetRValues(AUXVAR_INTERNAL, var_id, ncells, data)
   class is (goveqn_cair_vapor_type)
      call goveq%GetRValues(AUXVAR_INTERNAL, var_id, ncells, data)
   class is (goveqn_cleaf_temp_type)
      call goveq%GetRValues(AUXVAR_INTERNAL, var_id, ncells, data)
   end select

   end subroutine get_data_from_mlc_eqn

   !------------------------------------------------------------------------
  subroutine extract_data_from_mlc(mlc_mpp, istep, isubstep, dt)
    !
    ! !DESCRIPTION:
    !   Extracts following variables from the MLC model:
    !     - Tleaf_sun
    !     - Tleaf_shd
    !     - Tsoil
    !     - Tair
    !     - eair
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use MultiPhysicsProbMLC             , only : mpp_mlc_type
    use ml_model_global_vars            , only : nbot, ntop, ncair, ntree, nz_cair, output_data
    use ml_model_global_vars            , only : CLEF_TEMP_SUN_GE, CLEF_TEMP_SHD_GE, CAIR_TEMP_GE, CAIR_VAPR_GE
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType  , only : goveqn_cair_temp_type
    use GoveqnCanopyAirVaporType        , only : goveqn_cair_vapor_type
    use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants       , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants       , only : VAR_LEAF_TEMPERATURE
    use MultiPhysicsProbConstants       , only : VAR_LEAF_HEAT_STORAGE
    use MultiPhysicsProbConstants       , only : MM_H2O, MM_DRY_AIR
    use MultiPhysicsProbConstants       , only : VAR_SENSIBLE_HEAT_FLUX
    use MultiPhysicsProbConstants       , only : VAR_LATENT_HEAT_FLUX
    use MultiPhysicsProbConstants       , only : VAR_LEAF_TEMPERATURE
    use MultiPhysicsProbConstants       , only : VAR_LEAF_TRANSPIRATION
    use MultiPhysicsProbConstants       , only : VAR_WATER_VAPOR
    use ml_model_utils                  , only : accumulate_data
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type)                  :: mlc_mpp
    PetscInt                             :: istep, isubstep
    PetscReal                            :: dt
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscScalar                , pointer :: soln_p(:)
    PetscInt                             :: idx_leaf, idx_data, idx_soil, idx_air
    PetscInt                             :: ileaf, icair, itree, k, ieqn
    PetscInt                             :: ncells
    PetscReal,  pointer                  :: tleaf_data(:), tair_data(:), qair_data(:)
    PetscReal,  pointer                  :: sh_data(:), lh_data(:), st_data(:), tr_data(:)
    PetscErrorCode                       :: ierr
    PetscReal                            :: factor
    character(len=20)                    :: step_string, substep_string

    ncells = ncair*ntree*(nz_cair+1)
    allocate(tleaf_data (ncells))

    ncells = ncair*ntree*(nz_cair+1)
    allocate(tair_data  (ncells))
    allocate(qair_data  (ncells))

    ncells = ncair*ntree*(nz_cair+1)*2
    allocate(sh_data    (ncells))
    allocate(lh_data    (ncells))
    allocate(st_data    (ncells))
    allocate(tr_data    (ncells))

    if (output_data) then
       write(step_string,*)istep
       write(substep_string,*)isubstep
       write(*,*)'mpp.tleaf{' // trim(adjustl(step_string)) // ',' //trim(adjustl(substep_string)) // '} = ['
    end if
    ncells = ncair*ntree*(nz_cair+1)
    do ileaf = 1, 2
       if (ileaf == 1) then
          ieqn = CLEF_TEMP_SUN_GE
       else
          ieqn = CLEF_TEMP_SHD_GE
       end if

       call get_data_from_mlc_eqn(mlc_mpp, ieqn, VAR_LEAF_TEMPERATURE , ncells, tleaf_data)
       call get_data_from_mlc_eqn(mlc_mpp, ieqn, VAR_LEAF_HEAT_STORAGE, ncells, st_data)

       idx_leaf = 0
       idx_data = 0
       do icair = 1, ncair
          do itree = 1, ntree
             do k = 1, nz_cair+1
                if (k>=nbot .and. k<=ntop) then
                   idx_leaf = idx_leaf + 1
                   idx_data = (icair-1)*ntree*(nz_cair+1) + (ntree-1)*(nz_cair+1) + k
                   if (output_data) then
                      write(*,*)idx_data, tleaf_data(idx_data)
                   endif
                   if (ileaf == 1) then
                      call set_value_in_condition(int_cond%Tleaf_sun, idx_leaf, tleaf_data(idx_data))

                      call accumulate_data(vert_lev_vars%sh_leaf_sun, sh_data(idx_data)/dt, idx_leaf, isubstep)
                   else
                      call set_value_in_condition(int_cond%Tleaf_shd, idx_leaf, tleaf_data(idx_data))

                      call accumulate_data(vert_lev_vars%sh_leaf_shd, sh_data(idx_data)/dt, idx_leaf, isubstep)
                   endif
                end if
             end do
          end do
       end do
    end do
    if (output_data) then
       write(*,*)'];'
    end if

    ncells = ncair*(nz_cair+1)
    call get_data_from_mlc_eqn(mlc_mpp, CAIR_TEMP_GE, VAR_TEMPERATURE       , ncells, tair_data)
    call get_data_from_mlc_eqn(mlc_mpp, CAIR_TEMP_GE, VAR_SENSIBLE_HEAT_FLUX, ncells, sh_data)

    call get_data_from_mlc_eqn(mlc_mpp, CAIR_VAPR_GE, VAR_WATER_VAPOR       , ncells, qair_data)
    call get_data_from_mlc_eqn(mlc_mpp, CAIR_VAPR_GE, VAR_LATENT_HEAT_FLUX  , ncells, lh_data)
    call get_data_from_mlc_eqn(mlc_mpp, CAIR_VAPR_GE, VAR_LEAF_TRANSPIRATION, ncells, tr_data)

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    idx_soil = 0
    idx_data = 0
    idx_air  = 0
    idx_leaf = 0

    if (output_data) then
       write(step_string,*)istep
       write(substep_string,*)isubstep
       write(*,*)'mpp.air{' // trim(adjustl(step_string)) // ',' //trim(adjustl(substep_string)) // '} = ['
    end if
    do icair = 1, ncair
       do k = 1, nz_cair+1
          idx_data = idx_data + 1
          if (k == 1) then
             idx_soil = idx_soil + 1
             call set_value_in_condition(bnd_cond%tg, idx_soil, tair_data(idx_data))
          else
             idx_air = idx_air + 1
             call set_value_in_condition(int_cond%Tair, idx_air, tair_data(idx_data))
             call set_value_in_condition(int_cond%qair, idx_air, qair_data(idx_data))
             if (output_data) then
                write(*,*)idx_data,tair_data(idx_data),qair_data(idx_data)
             end if
             call set_value_in_condition(int_cond%wind, idx_air, soe%cturb%wind(icair,k))
          end if
          if (k >= nbot .and. k <= ntop) then
             idx_leaf = idx_leaf + 1

             call accumulate_data(vert_lev_vars%sh_leaf_sun, sh_data((idx_data-1)*2 + 1)   , idx_leaf, isubstep)
             call accumulate_data(vert_lev_vars%sh_leaf_shd, sh_data((idx_data-1)*2 + 2)   , idx_leaf, isubstep)

             call accumulate_data(vert_lev_vars%lh_leaf_sun, lh_data((idx_data-1)*2 + 1)   , idx_leaf, isubstep)
             call accumulate_data(vert_lev_vars%lh_leaf_shd, lh_data((idx_data-1)*2 + 2)   , idx_leaf, isubstep)

             call accumulate_data(vert_lev_vars%tr_leaf_sun, tr_data((idx_data-1)*2 + 1)   , idx_leaf, isubstep)
             call accumulate_data(vert_lev_vars%tr_leaf_shd, tr_data((idx_data-1)*2 + 2)   , idx_leaf, isubstep)
          end if
       end do
    end do
    if (output_data) then
       write(*,*)'];'
    endif

    deallocate(tleaf_data)
    deallocate(tair_data)
    deallocate(qair_data)
    deallocate(sh_data)
    deallocate(lh_data)
    deallocate(tr_data)

  end subroutine extract_data_from_mlc

  !------------------------------------------------------------------------
  subroutine solve_mlc(mlc_mpp, istep, isubstep, dt)
    !
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp
    PetscInt            :: isubstep, istep
    PetscReal           :: dt
    !
    PetscBool           :: converged
    PetscInt            :: converged_reason
    PetscBool           :: flg
    PetscErrorCode      :: ierr

    call set_boundary_conditions(mlc_mpp, istep, isubstep)

    call mlc_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

    call extract_data_from_mlc(mlc_mpp, istep, isubstep, dt)

  end subroutine solve_mlc

  !------------------------------------------------------------------------
  subroutine checkpoint_mlc(mlc_mpp, istep, isubstep)
    !
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use ml_model_global_vars            , only : nbot, ntop, ncair, ntree, nz_cair, output_data
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp
    PetscInt            :: isubstep, istep
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscViewer                          :: viewer
    PetscErrorCode                       :: ierr
    PetscInt                             :: nDM
    DM                         , pointer :: dms(:)
    Vec                        , pointer :: soln_subvecs(:)
    PetscReal                  , pointer :: c_p(:), s_p(:)
    character(len=240)                   :: step_string , substep_string, filename
    Vec                                  :: checkpoint_vec
    PetscInt                             :: ii, kk, icair, offset

    write(step_string,'(I0.3)')istep
    write(substep_string,*)isubstep
    write(filename,*)'mlc_checkpoint.' // trim(adjustl(step_string)) // '.' //trim(adjustl(substep_string)) // '.bin'
    write(*,*)'%filename: ',trim(adjustl(filename))
    
    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    call VecCreate(PETSC_COMM_SELF, checkpoint_vec, ierr); CHKERRA(ierr)
    call VecSetSizes(checkpoint_vec, 93*5, PETSC_DECIDE, ierr); CHKERRA(ierr)
    call VecSetFromOptions(checkpoint_vec, ierr); CHKERRA(ierr)

    call VecGetArrayF90(checkpoint_vec, c_p, ierr)
    icair = 1;
    offset = 0
    do kk = 1, nz_cair + 1
       c_p(kk + offset) = soe%cturb%wind(icair,kk)
    end do

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

    offset = nz_cair + 1
    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), s_p, ierr)
       do kk = 1, nz_cair + 1
          c_p(kk + offset) = s_p(kk)
       end do
       offset = offset + nz_cair + 1
       call VecRestoreArrayF90(soln_subvecs(ii), s_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(mlc_mpp%soe%solver%dm, &
         mlc_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)

    call VecRestoreArrayF90(checkpoint_vec, c_p, ierr)

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(adjustl(filename)), FILE_MODE_WRITE, viewer, ierr); CHKERRA(ierr)
    call VecView(checkpoint_vec, viewer, ierr); CHKERRA(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRA(ierr)

  end subroutine checkpoint_mlc

end module mlc
