module photosynthesis

  use mpp_varctl                          , only : iulog
  use mpp_abortutils                      , only : endrun
  use mpp_shr_log_mod                     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbPhotosynthesis      , only : mpp_photosynthesis_type
  use SystemOfEquationsPhotosynthesisType , only : sysofeqns_photosynthesis_type
  use SystemOfEquationsBaseType           , only : sysofeqns_base_type
  use GoverningEquationBaseType           , only : goveqn_base_type
  use GoveqnPhotosynthesisType            , only : goveqn_photosynthesis_type
  use ml_model_global_vars
  use petscsys
  use petscdm
  use petscdmda

  implicit none

#include <petsc/finclude/petsc.h>

  public :: init_photosynthesis
  public :: photosynthesis_set_boundary_conditions
  public :: solve_photosynthesis

contains

  !------------------------------------------------------------------------
  subroutine initialize(psy_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_PHOTOSYNTHESIS_SNES
    !
    implicit none

    class(mpp_photosynthesis_type) :: psy_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call psy_mpp%Init       ()
    call psy_mpp%SetName    ('Photosynthesis model')
    call psy_mpp%SetID      (MPP_PHOTOSYNTHESIS_SNES)
    call psy_mpp%SetMPIRank (iam)

  end subroutine initialize

  !------------------------------------------------------------------------
  subroutine setup_meshes(psy_mpp)
    !
    use MeshType             , only : mesh_type
    use ml_model_meshes      , only : create_canopy_mesh_for_leaf
    use ml_model_global_vars , only : PHOTOSYNTHESIS_MESH
    !
    implicit none

    class(mpp_photosynthesis_type) :: psy_mpp

    class(mesh_type), pointer :: mesh

    PHOTOSYNTHESIS_MESH = 1
    call create_canopy_mesh_for_leaf(mesh)
    call psy_mpp%SetNumMeshes(1)
    call psy_mpp%AddMesh(PHOTOSYNTHESIS_MESH, mesh)

    deallocate(mesh)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine add_goveqns(psy_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_PHOTOSYNTHESIS
    use ml_model_global_vars      , only : PHOTOSYNTHESIS_MESH, PHOTOSYNTHESIS_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_photosynthesis_type) :: psy_mpp

    PHOTOSYNTHESIS_GE = 1
    call psy_mpp%AddGovEqnWithMeshRank(GE_PHOTOSYNTHESIS, 'Photosynthesis model', PHOTOSYNTHESIS_MESH)
    call psy_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine set_parameters(psy_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use ConditionType        , only : condition_type
    use ConnectionSetType    , only : connection_set_type
    use ml_model_global_vars , only : PHOTOSYNTHESIS_MESH, PHOTOSYNTHESIS_GE, c3psn, gstype
    use MultiPhysicsProbConstants , only : VAR_PHOTOSYNTHETIC_PATHWAY_C4
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE_MEDLYN
    use WaterVaporMod        , only : SatVap
    use ml_model_meshes      , only : nleaf
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_photosynthesis_type)        :: psy_mpp
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    class(connection_set_type) , pointer :: cur_conn_set
    class(condition_type)      , pointer :: cur_cond
    PetscInt                             :: k, icol, icell, ileaf, iconn, sum_conn, nz, ncol

    PetscReal                , parameter :: tau = 0.1d0
    PetscReal                , parameter :: rho = 0.1d0

    c3psn  = VAR_PHOTOSYNTHETIC_PATHWAY_C4
    gstype = VAR_STOMATAL_CONDUCTANCE_MEDLYN

    call psy_mpp%soe%SetPointerToIthGovEqn(PHOTOSYNTHESIS_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_photosynthesis_type)

       ncol = ncair * ntree
       nz   = (ntop-nbot+1)

       icell = 0
       do icol = 1, ncol
          do ileaf = 1, nleaf
             do k = 1, nz
                icell = icell + 1

                cur_goveq%aux_vars_in(icell)%colim  = 1
                cur_goveq%aux_vars_in(icell)%c3psn  = c3psn
                cur_goveq%aux_vars_in(icell)%gstype = gstype

                cur_goveq%aux_vars_in(icell)%cair  = 380.d0                        ! [mmol/mol]
                cur_goveq%aux_vars_in(icell)%o2ref = 0.209d0 * 1000.d0             ! [mmol/mol]
                cur_goveq%aux_vars_in(icell)%apar  = 2000.d0 * ( 1.d0 - rho - tau) ! [mmol photon/m2/s]

                cur_goveq%aux_vars_in(icell)%btran = 1.d0
                cur_goveq%aux_vars_in(icell)%dpai  = dpai(nbot-1 + k)
                cur_goveq%aux_vars_in(icell)%dpai  = 1.d0
             end do
          end do
       end do

    end select

  end subroutine set_parameters

  !------------------------------------------------------------------------
  subroutine photosynthesis_set_boundary_conditions(psy_mpp)

    ! !DESCRIPTION
    !
    ! !USES:
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use SystemOfEquationsPhotosynthesisType , only : sysofeqns_photosynthesis_type
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnPhotosynthesisType            , only : goveqn_photosynthesis_type
    use ConditionType                       , only : condition_type
    use ConnectionSetType                   , only : connection_set_type
    use MultiPhysicsProbConstants           , only : VAR_PHOTOSYNTHETIC_PATHWAY_C3, VAR_PHOTOSYNTHETIC_PATHWAY_C4
    use MultiPhysicsProbConstants           , only : VAR_STOMATAL_CONDUCTANCE_MEDLYN, VAR_STOMATAL_CONDUCTANCE_BBERRY
    use MultiPhysicsProbConstants           , only : TFRZ
    use ml_model_global_vars                , only : Tleaf_sun, Tleaf_shd, Tair, gbh, gbc, Ileaf_sun_vis, Ileaf_shd_vis
    use WaterVaporMod                       , only : SatVap
    use ml_model_utils                      , only : get_value_from_condition
    use ml_model_meshes                     , only : nleaf
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_photosynthesis_type)        :: psy_mpp
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    class(connection_set_type) , pointer :: cur_conn_set
    class(condition_type)      , pointer :: cur_cond
    PetscInt                             :: k, icol, icell, iconn, sum_conn, nz, ncol, leaf_count, ileaf, idx_data, icair
    PetscReal                            :: eair, vpd_tleaf
    PetscReal                            :: esat_tair, desat_tair
    PetscReal                            :: esat_tleaf, desat_tleaf
    PetscReal                  , pointer :: tair_local(:), tleaf_local(:)
    PetscReal                            :: Isun_vis, Ishd_vis
    PetscReal                , parameter :: relhum = 80.d0
    PetscReal                , parameter :: unit_conversion = 4.6 ! w/m2 to mmol_photons/m2/s

    call psy_mpp%soe%SetPointerToIthGovEqn(PHOTOSYNTHESIS_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_photosynthesis_type)

       ncol = ncair
       nz   = (ntop-nbot+1)

       allocate(tair_local (ncol * nz * 2))
       allocate(tleaf_local(ncol * nz * 2))

       icell = 0
       ileaf = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1
             ileaf = ileaf + 1

             tleaf_local(ileaf           ) = get_value_from_condition(Tleaf_sun, ileaf) ! [K]
             tleaf_local(ileaf + ncol*nz ) = get_value_from_condition(Tleaf_shd, ileaf) ! [K]
             !tleaf_local(icell)            = TFRZ + 11.d0
             !tleaf_local(icell + ncol*nz ) = TFRZ + 11.d0

             Isun_vis = get_value_from_condition(Ileaf_sun_vis, icell) ! w/m2
             Ishd_vis = get_value_from_condition(Ileaf_shd_vis, icell) ! w/m2

             cur_goveq%aux_vars_in(ileaf          )%apar = Isun_vis * unit_conversion ! [mmol_photon/m2/s]
             cur_goveq%aux_vars_in(ileaf + ncol*nz)%apar = Ishd_vis * unit_conversion ! [mmol_photon/m2/s]
             cur_goveq%aux_vars_in(ileaf          )%apar = 2000.d0 *(1.d0 - 0.1d0 - 0.1d0)
             cur_goveq%aux_vars_in(ileaf + ncol*nz)%apar = 2000.d0 *(1.d0 - 0.1d0 - 0.1d0)
          end do
       end do

       icell = 0
       idx_data = 0
       do icair = 1, ncair
          do k = 1, nz_cair
              idx_data = idx_data + 1
              if (k >= nbot .and. k<=ntop) then
                 icell = icell + 1
                 tair_local(icell           )  = get_value_from_condition(Tair, idx_data)  ! [K]
                 tair_local(icell + ncol*nz )  = get_value_from_condition(Tair, idx_data)  ! [K]
                 tair_local(icell)            = TFRZ + 25.d0
                 tair_local(icell + ncol*nz ) = TFRZ + 25.d0
              endif
            enddo
       enddo

       do icell = 1, ncol*nz*nleaf

          call SatVap (tair_local( icell) , esat_tair, desat_tair )
          call SatVap (tleaf_local(icell), esat_tleaf, desat_tleaf)
          eair = esat_tair * relhum /100.d0
          vpd_tleaf = esat_tair - eair

          cur_goveq%aux_vars_in(icell)%tleaf = tleaf_local(icell)
           cur_goveq%aux_vars_in(icell)%gbv   = get_value_from_condition(gbh, icell)
           cur_goveq%aux_vars_in(icell)%gbc   = get_value_from_condition(gbc, icell)
           cur_goveq%aux_vars_in(icell)%gbv   = 2.224407920268566d0
           cur_goveq%aux_vars_in(icell)%gbc   = 1.637448199187622d0

             if (cur_goveq%aux_vars_in(icell)%gstype == VAR_STOMATAL_CONDUCTANCE_MEDLYN) then
                cur_goveq%aux_vars_in(icell)%eair = esat_tleaf - vpd_tleaf
             else
                cur_goveq%aux_vars_in(icell)%eair = esat_tleaf * relhum/100.d0
             end if

       end do

    end select

  end subroutine photosynthesis_set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine init_photosynthesis(psy_mpp)
    !
    implicit none
    !
    class(mpp_photosynthesis_type) :: psy_mpp
    !
    PetscScalar, pointer          :: ci_p(:)
    PetscBool                :: flg
    PetscErrorCode           :: ierr

    call initialize(psy_mpp)

    call setup_meshes(psy_mpp)

    call add_goveqns(psy_mpp)

    call psy_mpp%AllocateAuxVars()

    call psy_mpp%SetupProblem()

    call set_parameters(psy_mpp)

    call set_initial_conditions(psy_mpp)

  end subroutine init_photosynthesis

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(psy_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use ml_model_utils            , only : get_value_from_condition
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_photosynthesis_type) :: psy_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_photosynthesis_type)  , pointer :: soe
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

    base_soe => psy_mpp%soe

    select type(base_soe)
    class is (sysofeqns_photosynthesis_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    ! Find number of GEs packed within the SoE
    call DMCompositeGetNumberDM(psy_mpp%soe%solver%dm, nDM, ierr)

    ! Get DMs for each GE
    allocate (dms(nDM))
    call DMCompositeGetEntriesArray(psy_mpp%soe%solver%dm, dms, ierr)

    ! Allocate vectors for individual GEs
    allocate(soln_subvecs(nDM))

    ! Get solution vectors for individual GEs
    call DMCompositeGetAccessArray(psy_mpp%soe%solver%dm, &
         psy_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)
       v_p(:) = 152.d0
       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(psy_mpp%soe%solver%dm, &
         psy_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)

    call VecCopy(psy_mpp%soe%solver%soln, psy_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(psy_mpp%soe%solver%soln, psy_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

    call psy_mpp%soe%PreSolve()

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine solve_photosynthesis(psy_mpp, istep, dt)
    !
    implicit none
    !
    class(mpp_photosynthesis_type) :: psy_mpp
    PetscInt                 :: istep
    PetscReal                :: dt
    !
    PetscScalar, pointer          :: ci_p(:)
    PetscBool                :: converged
    PetscInt                 :: converged_reason
    PetscBool                :: flg
    PetscErrorCode           :: ierr

    call photosynthesis_set_boundary_conditions(psy_mpp)

    !call VecGetArrayF90(psy_mpp%soe%solver%soln, ci_p, ierr); CHKERRQ(ierr)
    !ci_p(:) = 152.d0
    !call VecRestoreArrayF90(psy_mpp%soe%solver%soln, ci_p, ierr); CHKERRQ(ierr)

    call psy_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

  end subroutine solve_photosynthesis

end module photosynthesis
