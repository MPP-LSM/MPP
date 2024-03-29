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
  use MultiPhysicsProbConstants           , only : VAR_PHOTOSYNTHETIC_PATHWAY_C3, VAR_PHOTOSYNTHETIC_PATHWAY_C4
  use MultiPhysicsProbConstants           , only : VAR_SCM_MEDLYN, VAR_SCM_BBERRY, VAR_SCM_WUE
  use MultiPhysicsProbConstants           , only : VAR_SCM_BONAN14
  use MultiPhysicsProbConstants           , only : VAR_SCM_MODIFIED_BONAN14
  use MultiPhysicsProbConstants           , only : VAR_SCM_MANZONI11
  use MultiPhysicsProbConstants           , only : VAR_SCM_OSMWANG

  implicit none

#include <petsc/finclude/petsc.h>

  public :: init_photosynthesis
  public :: photosynthesis_set_boundary_conditions
  public :: solve_photosynthesis
  public :: extract_data_from_photosynthesis
  public :: checkpoint_photosynthesis
  public :: photosynthesis_initialize_from_checkpoint

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
    c3psn  = VAR_PHOTOSYNTHETIC_PATHWAY_C3

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
                
                cur_goveq%aux_vars_in(icell)%fwet  = 0.d0
                cur_goveq%aux_vars_in(icell)%fdry  = 0.8218390792391702d0

                ! Root parameters
                cur_goveq%aux_vars_in(icell)%root%biomass = 500.d0
                cur_goveq%aux_vars_in(icell)%root%radius  = 0.29d-3
                cur_goveq%aux_vars_in(icell)%root%density = 0.31d6
                cur_goveq%aux_vars_in(icell)%root%resist  = 25.d0

                ! Soil paramters
                call set_soil_parameters(cur_goveq%aux_vars_in(icell)%soil)

                ! Plant parameters
                cur_goveq%aux_vars_in(icell)%plant%nleaf = 1
                call cur_goveq%aux_vars_in(icell)%plant%AllocateMemory()

                cur_goveq%aux_vars_in(icell)%plant%leaf_psi(:)    = -2.4d0
                cur_goveq%aux_vars_in(icell)%plant%leaf_height(:) = (k-1)*0.5d0 + 2.75d0
                cur_goveq%aux_vars_in(icell)%plant%leaf_capc(:)   = 2500.d0
                cur_goveq%aux_vars_in(icell)%plant%leaf_minlwp(:) = -1.2d0
                cur_goveq%aux_vars_in(icell)%plant%leaf_lai(:)    = 4.1516127586364746d0
                cur_goveq%aux_vars_in(icell)%plant%k_stem2leaf(:) = 4.d0

                call cur_goveq%aux_vars_in(icell)%SetDefaultParameters()

                select case(gstype)
                case (VAR_SCM_MEDLYN)
                   cur_goveq%aux_vars_in(icell)%g0opt = 0.0001d0
                   cur_goveq%aux_vars_in(icell)%g1opt = 4.0d0

                case (VAR_SCM_BBERRY)
                   cur_goveq%aux_vars_in(icell)%g0opt = 0.027d0
                   cur_goveq%aux_vars_in(icell)%g1opt = 9.0d0

                case (VAR_SCM_WUE)
                   cur_goveq%aux_vars_in(icell)%iota = 820.0d0
                   cur_goveq%aux_vars_in(icell)%plant%leaf_minlwp(:) = -2.5d0

                case (VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14)
                   cur_goveq%aux_vars_in(icell)%plant%leaf_minlwp(:) = -2.5d0
                   cur_goveq%aux_vars_in(icell)%iota = 820.0d0

                case (VAR_SCM_MANZONI11)
                  cur_goveq%aux_vars_in(icell)%plant%leaf_minlwp(:) = -2.5d0
                  cur_goveq%aux_vars_in(icell)%manzoni11_beta = -0.001d0
                  cur_goveq%aux_vars_in(icell)%iota = 820.0d0

                end select

             end do
          end do
       end do

    end select

  end subroutine set_parameters

  !------------------------------------------------------------------------
  subroutine set_soil_parameters(soil)
   !
   use PhotosynthesisAuxType, only : soil_auxvar_type
   !
   implicit none
   !
   type(soil_auxvar_type), pointer :: soil
   !
   PetscInt                            :: j, texture
   PetscReal                           :: beta_param, z1, z2
   PetscReal, parameter, dimension(11) :: theta_sat = [0.395d0, 0.410d0, 0.435d0, 0.485d0, 0.451d0, 0.420d0, 0.477d0, 0.476d0, 0.426d0, 0.492d0, 0.482d0]
   PetscReal, parameter, dimension(11) :: psi_sat   = [-121.d0, -90.d0, -218d0, -786.d0, -478.d0, -299.d0, -356.d0, -630.d0, -153.d0, -490.d0, -405.d0]
   PetscReal, parameter, dimension(11) :: b         = [4.05d0, 4.38d0, 4.90d0, 5.30d0, 5.39d0, 7.12d0, 7.75d0, 8.52d0, 10.40d0, 10.40d0, 11.40d0];
   PetscReal, parameter, dimension(11) :: k_sat     = [1.056d0, 0.938d0, 0.208d0, 0.0432d0, 0.0417d0, 0.0378d0, 0.0102d0, 0.0147d0, 0.0130d0, 0.0062d0, 0.0077d0];
   PetscReal                           :: zi(0:10)
   PetscReal, parameter                :: m_to_cm = 1.d2

   texture = 1;

   soil%nlevsoi = 11
   soil%nlevsoi = 10

   call soil%AllocateMemory()

   soil%dz( 1 ) =  1.7512817916255204d-002
   soil%dz( 2 ) =  2.7578969259676251d-002
   soil%dz( 3 ) =  4.5470033242413201d-002
   soil%dz( 4 ) =  7.4967410986208557d-002
   soil%dz( 5 ) =  0.12360036510228053d0
   soil%dz( 6 ) =  0.20378255101043175d0
   soil%dz( 7 ) =  0.33598062644843263d0
   soil%dz( 8 ) =  0.55393840536868488d0
   soil%dz( 9 ) =  0.91329003158906108d0
   soil%dz(10 ) =  1.5057607013992766d0

   zi(           0 ) =    0.0000000000000000d0
   zi(           1 ) =    1.7512817916255204d-002
   zi(           2 ) =    4.5091787175931458d-002
   zi(           3 ) =    9.0561820418344652d-002
   zi(           4 ) =   0.16552923140455322d0
   zi(           5 ) =   0.28912959650683373d0
   zi(           6 ) =   0.49291214751726548d0
   zi(           7 ) =   0.82889277396569816d0
   zi(           8 ) =    1.3828311793343830d0
   zi(           9 ) =    2.2961212109234443d0
   zi(          10 ) =    3.8018819123227208d0

   beta_param = 0.90d0;  ! root profile parameter: shallow profile
   !beta_param = 0.97d0;  ! root profile parameter: deep profile
   beta_param  = 0.966d0

   do j = 1, soil%nlevsoi

      soil%rootfr(j) = ( beta_param ** (zi(j-1)*m_to_cm) - beta_param ** (zi(j)*m_to_cm) )

      soil%watsat(j)     = theta_sat(texture)
      soil%hksat(j)      = k_sat(texture) * 10.d0/60.d0
      soil%bsw(j)        = b(texture)

      soil%h2osoi_vol(j) = 0.5d0*soil%watsat(j)
      soil%psi_sat(j)    = psi_sat(texture)
      soil%psi(j)        = psi_sat(texture) * (soil%h2osoi_vol(j)/soil%watsat(j))**(-soil%bsw(j))
   end do

 end subroutine set_soil_parameters

 !------------------------------------------------------------------------
  subroutine photosynthesis_set_boundary_conditions(psy_mpp, istep, isubstep, is_first_substep)

    ! !DESCRIPTION
    !
    ! !USES:
    use MultiPhysicsProbConstants           , only : MM_H2O, MM_DRY_AIR
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use SystemOfEquationsPhotosynthesisType , only : sysofeqns_photosynthesis_type
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnPhotosynthesisType            , only : goveqn_photosynthesis_type
    use ConditionType                       , only : condition_type
    use ConnectionSetType                   , only : connection_set_type
    use MultiPhysicsProbConstants           , only : TFRZ
    use WaterVaporMod                       , only : SatVap
    use ml_model_utils                      , only : get_value_from_condition
    use ml_model_meshes                     , only : nleaf
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_photosynthesis_type)        :: psy_mpp
    PetscInt :: istep, isubstep
    PetscBool :: is_first_substep
    !
    class(goveqn_base_type)    , pointer   :: cur_goveq
    class(connection_set_type) , pointer   :: cur_conn_set
    class(condition_type)      , pointer   :: cur_cond
    PetscInt                               :: k, icol, icell, iconn, sum_conn, nz, ncol, leaf_count, ileaf, idx_data, icair, j
    PetscReal                              :: Isun_vis, Ishd_vis
    PetscReal                              :: qair_value , pref_value, esat, desat
    PetscReal                  , pointer   :: tleaf_local(:)
    PetscReal                  , parameter :: unit_conversion = 4.6d0 ! w/m2 to mmol_photons/m2/s
    PetscInt                   , parameter :: nlev = 10
    character(len=20)                      :: step_string, substep_string

    call psy_mpp%soe%SetPointerToIthGovEqn(PHOTOSYNTHESIS_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_photosynthesis_type)

       ncol = ncair
       nz   = (ntop-nbot+1)

       allocate(tleaf_local(ncol * nz * 2))

       icell = 0
       ileaf = 0
       do icol = 1, ncol
          pref_value = get_value_from_condition(bnd_cond%pref, icol)
          do k = 1, nz
             icell = icell + 1
             ileaf = ileaf + 1

             tleaf_local(ileaf           ) = get_value_from_condition(int_cond%Tleaf_sun, ileaf) ! [K]
             tleaf_local(ileaf + ncol*nz ) = get_value_from_condition(int_cond%Tleaf_shd, ileaf) ! [K]

             Isun_vis = get_value_from_condition(int_cond%Ileaf_sun_vis, icell) ! w/m2
             Ishd_vis = get_value_from_condition(int_cond%Ileaf_shd_vis, icell) ! w/m2

             cur_goveq%aux_vars_in(ileaf          )%apar = Isun_vis * unit_conversion ! [mmol_photon/m2/s]
             cur_goveq%aux_vars_in(ileaf + ncol*nz)%apar = Ishd_vis * unit_conversion ! [mmol_photon/m2/s]
             cur_goveq%aux_vars_in(icell)%eair = get_value_from_condition(int_cond%qair, (icol-1)*nz + nbot + k - 2) * pref_value
          end do
       end do

       icell = 0

       write(step_string,*)istep
       write(substep_string,*)isubstep
       do icair = 1, ncair
          do ileaf = 1, nleaf
             do k = 1, nz
                icell = icell + 1

                if (.not.(istep == 1 .and. isubstep == 1)) then
                   cur_goveq%aux_vars_in(icell)%tleaf_prev = cur_goveq%aux_vars_in(icell)%tleaf
                endif
                cur_goveq%aux_vars_in(icell)%tleaf = tleaf_local(icell)

                cur_goveq%aux_vars_in(icell)%gbv   = get_value_from_condition(int_cond%gbv, icell)
                cur_goveq%aux_vars_in(icell)%gbc   = get_value_from_condition(int_cond%gbc, icell)

                cur_goveq%aux_vars_in(icell)%cair  = get_value_from_condition(bnd_cond%co2ref, icair)
                cur_goveq%aux_vars_in(icell)%o2ref = get_value_from_condition(bnd_cond%o2ref, icair)

                pref_value = get_value_from_condition(bnd_cond%pref, icair)
                if (isubstep == 1) then
                   pref_value = get_value_from_condition(bnd_cond%pref_prev, icair)
                end if
                cur_goveq%aux_vars_in(icell)%pref = pref_value

                qair_value = get_value_from_condition(int_cond%qair,  (icair-1)*nz + nbot + k - 2)
                cur_goveq%aux_vars_in(icell)%eair = qair_value * pref_value

                cur_goveq%aux_vars_in(icell)%tleaf     = tleaf_local(icell)

                cur_goveq%aux_vars_in(icell)%ci = get_value_from_condition(bnd_cond%co2ref, icair)

                ! - For run that uses IC, call PreSolve for all substeps
                ! - For run that does not uses IC, skip calling PreSolve only on the first substep
                if (use_ic .or. (.not. is_first_substep)) then
                   call cur_goveq%aux_vars_in(icell)%PreSolve()
                endif

                do j = 1, nlev
                   cur_goveq%aux_vars_in(icell)%soil%h2osoi_vol(j) = get_value_from_condition(bnd_cond%h2osoi_vol, j)
                end do

                call cur_goveq%aux_vars_in(icell)%DetermineIfSolutionIsBounded()

             end do
          end do
       end do

    end select

  end subroutine photosynthesis_set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine init_photosynthesis(psy_mpp)
    !
    use ml_model_global_vars      , only : gstype
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

    select case(gstype)
    case (VAR_SCM_BONAN14, VAR_SCM_MODIFIED_BONAN14)
       call psy_mpp%soe%SetDofsForGovEqn(PHOTOSYNTHESIS_GE,2)

    case (VAR_SCM_BBERRY, VAR_SCM_MEDLYN,VAR_SCM_WUE, VAR_SCM_MANZONI11)
        call psy_mpp%soe%SetDofsForGovEqn(PHOTOSYNTHESIS_GE,1)
    end select
    call psy_mpp%AllocateAuxVars()

    call psy_mpp%SetupProblem()

    call set_parameters(psy_mpp)

  end subroutine init_photosynthesis

  !------------------------------------------------------------------------
  subroutine set_initial_conditions(psy_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use ml_model_utils            , only : get_value_from_condition
    use ml_model_meshes           , only : nleaf
    use ml_model_global_vars      , only : gstype
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
    PetscInt                             :: ii, icair, k, ileaf, idx_data, icell, nz
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

    nz   = (ntop-nbot+1)

    do ii = 1, nDM
       call VecGetArrayF90(soln_subvecs(ii), v_p, ierr)
       call VecGetLocalSize(soln_subvecs(ii), icair, ierr)

       icell = 0
       do icair = 1, ncair
          do k = 1, nz
             do ileaf = 1, nleaf
                icell = icell + 1
                select case(gstype)
                   case (VAR_SCM_WUE, VAR_SCM_MANZONI11, VAR_SCM_OSMWANG)
                      v_p(icell) = 0.005d0
                   case (VAR_SCM_BONAN14,VAR_SCM_MODIFIED_BONAN14)
                      v_p((icell-1)*2 + 1) = 0.002d0
                      v_p((icell-1)*2 + 2) = 0.002d0
                   case (VAR_SCM_MEDLYN, VAR_SCM_BBERRY)
                      v_p(icell) = 0.9d0 * get_value_from_condition(bnd_cond%co2ref, icair)
                   case default
                      write(*,*)'Unknown gstype to set initial condition'
                   end select
             end do
          end do
       end do
       call VecRestoreArrayF90(soln_subvecs(ii), v_p, ierr)
    enddo

    ! Restore solution vectors for individual GEs
    call DMCompositeRestoreAccessArray(psy_mpp%soe%solver%dm, &
         psy_mpp%soe%solver%soln, nDM, &
         PETSC_NULL_INTEGER, soln_subvecs, ierr)

    deallocate(dms)

    call VecCopy(psy_mpp%soe%solver%soln, psy_mpp%soe%solver%soln_prev, ierr); CHKERRQ(ierr)
    call VecCopy(psy_mpp%soe%solver%soln, psy_mpp%soe%solver%soln_prev_clm, ierr); CHKERRQ(ierr)

  end subroutine set_initial_conditions

  !------------------------------------------------------------------------
  subroutine extract_data_from_photosynthesis(psy_mpp, istep, isubstep)
    !
    ! !DESCRIPTION:
    !   Extracts following variables from the photsynthesis model
    !     - gs
    !
    ! !USES:
    use ml_model_global_vars      , only : nbot, ntop, ncair, ntree, nz_cair, output_data
    use ml_model_meshes           , only : nleaf
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnPhotosynthesisType  , only : goveqn_photosynthesis_type
    use MultiPhysicsProbPhotosynthesis  , only : mpp_photosynthesis_type
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE
    use MultiPhysicsProbConstants , only : VAR_GROSS_PHOTOSYNTHESIS
    use MultiPhysicsProbConstants , only : VAR_NET_PHOTOSYNTHESIS
    use ml_model_utils            , only : set_value_in_condition, accumulate_data
    use ml_model_global_vars      , only : vert_lev_vars
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_photosynthesis_type)      :: psy_mpp
    PetscInt :: istep, isubstep
    !
    PetscInt                            :: icell, offset
    PetscInt                            :: ncells
    PetscReal               , pointer   :: gs_data(:), ag_data(:), an_data(:)
    class(goveqn_base_type) , pointer   :: goveq
    PetscErrorCode                      :: ierr
    character(len=20)                   :: step_string, substep_string

    ncells = ncair * ntree * (ntop - nbot + 1) * nleaf

    allocate(gs_data(ncells))
    allocate(ag_data(ncells))
    allocate(an_data(ncells))

    call psy_mpp%soe%SetPointerToIthGovEqn(1, goveq)

    select type(goveq)
    class is (goveqn_photosynthesis_type)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_STOMATAL_CONDUCTANCE, ncells, gs_data)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_GROSS_PHOTOSYNTHESIS, ncells, ag_data)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_NET_PHOTOSYNTHESIS  , ncells, an_data)
    end select

    if (output_data .and. isubstep == 12) then
       write(step_string,*)istep
       write(substep_string,*)isubstep
       write(*,*)'mpp.gs{' // trim(adjustl(step_string)) // ',' //trim(adjustl(substep_string)) // '} = ['
    end if
    offset = ncair * ntree * (ntop - nbot + 1)
    do icell = 1, ncair * ntree * (ntop - nbot + 1)
       call set_value_in_condition(int_cond%gs_sun, icell, gs_data(icell         ))
       call set_value_in_condition(int_cond%gs_shd, icell, gs_data(icell + offset))

       call accumulate_data(vert_lev_vars%gs_leaf_sun, gs_data(icell         ), icell, isubstep)
       call accumulate_data(vert_lev_vars%gs_leaf_shd, gs_data(icell + offset), icell, isubstep)

       call accumulate_data(vert_lev_vars%agross_leaf_sun, ag_data(icell         ), icell, isubstep)
       call accumulate_data(vert_lev_vars%agross_leaf_shd, ag_data(icell + offset), icell, isubstep)

       call accumulate_data(vert_lev_vars%anet_leaf_sun, an_data(icell         ), icell, isubstep)
       call accumulate_data(vert_lev_vars%anet_leaf_shd, an_data(icell + offset), icell, isubstep)

       if (output_data .and. isubstep == 12) then
          write(*,*)icell, gs_data(icell), gs_data(icell + offset)
       end if
    end do
    if (output_data .and. isubstep == 12) then
       write(*,*)'];'
    endif

    deallocate(gs_data)
    deallocate(ag_data)
    deallocate(an_data)

  end subroutine extract_data_from_photosynthesis

  !------------------------------------------------------------------------
  subroutine checkpoint_photosynthesis(psy_mpp, istep, isubstep)
    !
    use ml_model_global_vars , only : nbot, ntop, ncair, ntree, nz_cair, output_data
    use ml_model_meshes      , only : nleaf
    !
    implicit none
    !
    class(mpp_photosynthesis_type) :: psy_mpp
    PetscInt                       :: istep, isubstep
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscViewer                          :: viewer
    PetscErrorCode                       :: ierr
    PetscReal                  , pointer :: c_p(:)
    character(len=240)                   :: step_string , substep_string, filename
    Vec                                  :: checkpoint_vec
    PetscInt                             :: ncells, count, icell, ileaf, nvars

    nvars = &
      1 + & ! leaf_psi
      1 + & ! tleaf
      1     ! gleaf_w_soln

    ncells = ncair * ntree * (ntop - nbot + 1) * nleaf * nvars

    write(step_string,'(I0.3)')istep
    write(substep_string,*)isubstep
    write(filename,*)'photosynthesis_checkpoint.' // trim(adjustl(step_string)) // '.' //trim(adjustl(substep_string)) // '.bin'

    call VecCreate(PETSC_COMM_SELF, checkpoint_vec, ierr); CHKERRA(ierr)
    call VecSetSizes(checkpoint_vec, ncells, PETSC_DECIDE, ierr); CHKERRA(ierr)
    call VecSetFromOptions(checkpoint_vec, ierr); CHKERRA(ierr)

    call VecGetArrayF90(checkpoint_vec, c_p, ierr)

    call psy_mpp%soe%SetPointerToIthGovEqn(PHOTOSYNTHESIS_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_photosynthesis_type)

       count = 0
       do icell = 1,  ncair * ntree * (ntop - nbot + 1) * nleaf
          do ileaf = 1, 1!nleaf
             count = count + 1
             c_p(count) = cur_goveq%aux_vars_in(icell)%plant%leaf_psi(ileaf)

             count = count + 1
             c_p(count) = cur_goveq%aux_vars_in(icell)%tleaf

             count = count + 1
             c_p(count) = cur_goveq%aux_vars_in(icell)%gleaf_w_soln

          end do
       end do
    end select

    call VecRestoreArrayF90(checkpoint_vec, c_p, ierr)

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(adjustl(filename)), FILE_MODE_WRITE, viewer, ierr); CHKERRA(ierr)
    call VecView(checkpoint_vec, viewer, ierr); CHKERRA(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRA(ierr)

  end subroutine checkpoint_photosynthesis

  !------------------------------------------------------------------------
  subroutine photosynthesis_initialize_from_checkpoint(psy_mpp, filename)
    !
    use ml_model_global_vars , only : nbot, ntop, ncair, ntree, nz_cair, output_data
    use ml_model_meshes      , only : nleaf
    !
    implicit none
    !
    class(mpp_photosynthesis_type) :: psy_mpp
    character(len=*)                   :: filename
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscViewer                          :: viewer
    PetscErrorCode                       :: ierr
    PetscReal                  , pointer :: c_p(:)
    Vec                                  :: checkpoint_vec
    PetscInt                             :: ncells, count, icell, ileaf, nvars

    nvars = &
      1 + & ! leaf_psi
      1 + & ! tleaf
      1     ! gleaf_w_soln

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, viewer, ierr);CHKERRQ(ierr)
    call VecCreate(PETSC_COMM_WORLD, checkpoint_vec, ierr);CHKERRQ(ierr)
    call VecLoad(checkpoint_vec, viewer, ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)

    call VecGetArrayF90(checkpoint_vec, c_p, ierr)

    call psy_mpp%soe%SetPointerToIthGovEqn(PHOTOSYNTHESIS_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_photosynthesis_type)

       count = 0
       do icell = 1,  ncair * ntree * (ntop - nbot + 1) * nleaf
          do ileaf = 1, 1
             count = count + 1
             cur_goveq%aux_vars_in(icell)%plant%leaf_psi(ileaf) = c_p(count)

             count = count + 1
             cur_goveq%aux_vars_in(icell)%tleaf_prev = c_p(count)

             count = count + 1
             cur_goveq%aux_vars_in(icell)%gleaf_w_soln = c_p(count)

          end do
       end do
    end select

    call VecRestoreArrayF90(checkpoint_vec, c_p, ierr); CHKERRA(ierr)
    call VecDestroy(checkpoint_vec, ierr); CHKERRA(ierr)

  end subroutine photosynthesis_initialize_from_checkpoint

  !------------------------------------------------------------------------
  subroutine solve_photosynthesis(psy_mpp, istep, isubstep, is_first_substep, dt)
    !
    implicit none
    !
    class(mpp_photosynthesis_type) :: psy_mpp
    PetscInt                 :: istep, isubstep
    PetscBool                :: is_first_substep
    PetscReal                :: dt
    !
    PetscScalar, pointer          :: ci_p(:)
    PetscBool                :: converged
    PetscInt                 :: converged_reason
    PetscBool                :: flg
    PetscErrorCode           :: ierr

    call set_initial_conditions(psy_mpp)

    call photosynthesis_set_boundary_conditions(psy_mpp, istep, isubstep, is_first_substep)

    call psy_mpp%soe%StepDT(dt, (istep-1)*12 + isubstep, converged, converged_reason, ierr)

    if (.not.converged) then
       write(*,*)'Photosynthesis model did not converge'
       call exit(0)
    endif

    call extract_data_from_photosynthesis(psy_mpp, istep, isubstep)

  end subroutine solve_photosynthesis

end module photosynthesis
