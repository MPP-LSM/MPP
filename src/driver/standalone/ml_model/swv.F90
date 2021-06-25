module swv
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbShortwave , only : mpp_shortwave_type
  use mpp_abortutils            , only : endrun
  use mpp_varctl                , only : iulog
  use ml_model_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: init_swv
  public :: swv_set_boundary_conditions
  public :: solve_swv
  public :: extract_data_from_swv

contains

  !------------------------------------------------------------------------
  subroutine initialize(swv_mpp)
    !
    use MultiPhysicsProbConstants , only : MPP_shortwave_KSP
    !
    implicit none

    class(mpp_shortwave_type) :: swv_mpp

    PetscInt       :: iam
    PetscErrorCode :: ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, iam, ierr)

    !
    ! Set up the multi-physics problem
    !
    call swv_mpp%Init       ()
    call swv_mpp%SetName    ('Shortgwave radiation model')
    call swv_mpp%SetID      (MPP_shortwave_KSP)
    call swv_mpp%SetMPIRank (iam)

  end subroutine initialize

  !------------------------------------------------------------------------
  subroutine setup_meshes(swv_mpp)
    !
    use MeshType             , only : mesh_type
    use ml_model_meshes      , only : create_canopy_and_soil_mesh
    use ml_model_global_vars , only : SHORTWAVE_MESH
    !
    implicit none

    class(mpp_shortwave_type) :: swv_mpp

    class(mesh_type), pointer :: mesh

    SHORTWAVE_MESH = 1
    call create_canopy_and_soil_mesh(mesh)
    call swv_mpp%SetNumMeshes(1)
    call swv_mpp%AddMesh(SHORTWAVE_MESH, mesh)

    deallocate(mesh)

  end subroutine setup_meshes

  !------------------------------------------------------------------------
  subroutine add_goveqns(swv_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_SHORTWAVE
    use ml_model_global_vars      , only : SHORTWAVE_MESH, SHORTWAVE_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_shortwave_type) :: swv_mpp

    SHORTWAVE_GE = 1
    call swv_mpp%AddGovEqnWithMeshRank(GE_SHORTWAVE, 'Shortwave radiation model', SHORTWAVE_MESH)
    call swv_mpp%SetMeshesOfGoveqnsByMeshRank()

  end subroutine add_goveqns

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(swv_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : GE_SHORTWAVE
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use ml_model_global_vars      , only : SHORTWAVE_MESH, SHORTWAVE_GE
    use ConnectionSetType         , only : connection_set_type
    use ml_model_meshes           , only : create_connection_set_to_canopy_and_soil_mesh
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_shortwave_type)            :: swv_mpp
    class(connection_set_type) , pointer :: conn_set    

    call create_connection_set_to_canopy_and_soil_mesh(swv_mpp%meshes(SHORTWAVE_MESH), conn_set)

    call swv_mpp%soe%AddConditionInGovEqn( &
         SHORTWAVE_GE                 , &
         ss_or_bc_type = COND_BC      , &
         name = 'Atmospheric forcing' , &
         unit = 'K'                   , &
         cond_type = COND_DIRICHLET   , &
         conn_set = conn_set)
    
  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine set_parameters(swv_mpp, dlai, dsai, dpai)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsShortwaveType , only : sysofeqns_shortwave_type
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnShortwaveType            , only : goveqn_shortwave_type
    use ConditionType                  , only : condition_type
    use ConnectionSetType              , only : connection_set_type
    use ml_model_global_vars           , only : SHORTWAVE_MESH, SHORTWAVE_GE
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type)               :: swv_mpp
    PetscReal, pointer, intent(in)         :: dlai(:), dsai(:), dpai(:)
    !
    class(goveqn_base_type)    , pointer   :: cur_goveq
    class(connection_set_type) , pointer   :: cur_conn_set
    class(condition_type)      , pointer   :: cur_cond
    PetscInt                               :: k, icol, icell, ileaf, iconn, sum_conn, nz, ncol, idx
    PetscReal                              :: wl, ws, inc

    !PetscReal                  , parameter :: clumpfac  = 1.d0
    !PetscReal                  , parameter :: Kb        = 1.7628174450198393d0
    !PetscReal                  , parameter :: td        = 0.913235689378651d0
    PetscReal                  , parameter :: rho_l_vis = 0.1d0
    PetscReal                  , parameter :: rho_s_vis = 0.16d0
    PetscReal                  , parameter :: tau_l_vis = 5.d-2
    PetscReal                  , parameter :: tau_s_vis = 1.d-3

    PetscReal                  , parameter :: rho_l_nir = 0.45d0
    PetscReal                  , parameter :: rho_s_nir = 0.39d0
    PetscReal                  , parameter :: tau_l_nir = 0.25d0
    PetscReal                  , parameter :: tau_s_nir = 1.d-3

    PetscReal                  , parameter :: rho_min = 1.d-6
    PetscReal                  , parameter :: tau_min = 1.d-6
    PetscReal                  , parameter :: albedo_sat_vis = 0.07d0
    PetscReal                  , parameter :: albedo_sat_nir = 0.14d0
    PetscReal                  , parameter :: albedo_dry_vis = 0.14d0
    PetscReal                  , parameter :: albedo_dry_nir = 0.25d0
    !
    PetscReal                  , parameter :: h2osoi_vol = 0.10914649814367294

    call swv_mpp%soe%SetPointerToIthGovEqn(SHORTWAVE_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_shortwave_type)

       ncol = ncair * ntree
       nz   = (ntop-nbot+1) + 1

       icell = 0
       do icol = 1, ncol
          do k = 1, nz
             icell = icell + 1

             if (k == 1) then
                cur_goveq%aux_vars_in(icell)%is_soil = PETSC_TRUE

                inc = max(0.11d0 - 0.4d0*h2osoi_vol, 0.d0)
                cur_goveq%aux_vars_in(icell)%soil_albedo_b(1) = min(albedo_sat_vis + inc, albedo_dry_vis)
                cur_goveq%aux_vars_in(icell)%soil_albedo_b(2) = min(albedo_sat_nir + inc, albedo_dry_nir)

                cur_goveq%aux_vars_in(icell)%soil_albedo_d(1) = min(albedo_sat_vis + inc, albedo_dry_vis)
                cur_goveq%aux_vars_in(icell)%soil_albedo_d(2) = min(albedo_sat_nir + inc, albedo_dry_nir)

             else
                idx = nbot + k - 2
                wl = dlai(idx)/dpai(idx)
                ws = dsai(idx)/dpai(idx)

                cur_goveq%aux_vars_in(icell)%leaf_rho(1)   = max(rho_l_vis*wl + rho_s_vis*ws, rho_min)
                cur_goveq%aux_vars_in(icell)%leaf_rho(2)   = max(rho_l_nir*wl + rho_s_nir*ws, rho_min)

                cur_goveq%aux_vars_in(icell)%leaf_tau(1)   = max(tau_l_vis*wl + tau_s_vis*ws, tau_min)
                cur_goveq%aux_vars_in(icell)%leaf_tau(2)   = max(tau_l_nir*wl + tau_s_nir*ws, tau_min)

                cur_goveq%aux_vars_in(icell)%leaf_omega(1) = cur_goveq%aux_vars_in(icell)%leaf_rho(1) + cur_goveq%aux_vars_in(icell)%leaf_tau(1)
                cur_goveq%aux_vars_in(icell)%leaf_omega(2) = cur_goveq%aux_vars_in(icell)%leaf_rho(2) + cur_goveq%aux_vars_in(icell)%leaf_tau(2)

                cur_goveq%aux_vars_in(icell)%leaf_dpai        = dpai(k + nbot - 2)
             end if
          end do
       end do

    end select

  end subroutine set_parameters

  !------------------------------------------------------------------------
  subroutine compute_kb(xl, sza, phi1, phi2, kb)
    !
    implicit none
    !
    PetscReal, intent(in) :: xl, sza
    PetscReal, intent(out) :: phi1, phi2, kb
    !
    PetscReal :: chil, gdir
    PetscReal, parameter :: kb_max = 40.d0
    PetscReal, parameter :: chil_min = -0.4d0
    PetscReal, parameter :: chil_max =  0.6d0

    chil = min(max(xl, chil_min),chil_max)

    if (abs(chil) <= 0.01d0) chil = 0.01d0

    phi1 = 0.5d0 - 0.633d0 * chil - 0.330d0 * chil * chil
    phi2 = 0.877d0 * (1.d0 - 2.d0 * phi1)

    gdir = phi1 + phi2 * cos(sza)

    ! Direct beam extinction coefficient
    kb = gdir / cos(sza)
    kb = min(kb, kb_max)

  end subroutine compute_kb

  !------------------------------------------------------------------------
  subroutine compute_transmittance_coefficents(xl, sza, dpai, clump_fac, tb, td)
    !
    implicit none
    !
    PetscReal, intent(in) :: xl, sza, dpai, clump_fac
    PetscReal, intent(out) :: tb, td
    !
    PetscReal :: chil, phi1, phi2, gdir, gdirj, kb, angle
    PetscInt  :: j
    PetscReal, parameter :: pi = 4.d0 * atan(1.0d0)

    ! Direct beam extinction coefficient
    call compute_kb(xl, sza, phi1, phi2, kb)

    ! Canopy layer transmittance of direct beam radiation
    tb = exp(-kb * dpai * clump_fac)

    ! Canopy layer transmittance of diffuse radiation
    td = 0.d0
    do j = 1, 9
       angle = (5.d0 + (j - 1) * 10.d0) * pi / 180.d0
       gdirj = phi1 + phi2 * cos(angle)
       td = td + exp(-gdirj / cos(angle) * dpai * clump_fac) * sin(angle) * cos(angle)
    end do
    td = td * 2.d0 * (10.d0 * pi / 180.d0)

  end subroutine compute_transmittance_coefficents

  !------------------------------------------------------------------------
  subroutine swv_set_boundary_conditions(swv_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType      , only : sysofeqns_base_type
    use SystemOfEquationsShortwaveType , only : sysofeqns_shortwave_type
    use GoverningEquationBaseType      , only : goveqn_base_type
    use GoveqnShortwaveType            , only : goveqn_shortwave_type
    use ConditionType                  , only : condition_type
    use ConnectionSetType              , only : connection_set_type
    use ml_model_utils                 , only : get_value_from_condition
    use ml_model_global_vars           , only : dpai, sumpai, cumpai
    use ml_model_utils                 , only : compute_fssh
    use ShortwaveAuxType               , only : shortwave_auxvar_type
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type)             :: swv_mpp
    !
    class(goveqn_base_type)    , pointer :: cur_goveq
    class(connection_set_type) , pointer :: cur_conn_set
    class(condition_type)      , pointer :: cur_cond
    PetscInt                             :: icair, itree, k
    PetscInt                             :: icell, iconn, sum_conn, nz
    PetscReal                            :: sza_value, sumpai_value, dpai_value, cumpai_value
    PetscReal                            :: Kb, tb, tbcum, td
    PetscReal                            :: phi1, phi2
    type(shortwave_auxvar_type), pointer :: avars(:)
    !
    PetscReal, parameter                 :: xl        = 0.25d0
    PetscReal, parameter                 :: clump_fac = 1.d0

    call swv_mpp%soe%SetPointerToIthGovEqn(SHORTWAVE_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_shortwave_type)

       avars => cur_goveq%aux_vars_in

       nz   = (ntop-nbot+1) + 1

       icell = 0
       do icair = 1, ncair
          sza_value = get_value_from_condition(bnd_cond%sza, icair)

          do itree = 1, ntree
             do k = 1, nz
                icell = icell + 1

                avars(icell)%Iskyb(1) = get_value_from_condition(bnd_cond%Iskyb_vis, icair)
                avars(icell)%Iskyb(2) = get_value_from_condition(bnd_cond%Iskyb_nir, icair)
                avars(icell)%Iskyd(1) = get_value_from_condition(bnd_cond%Iskyd_vis, icair)
                avars(icell)%Iskyd(2) = get_value_from_condition(bnd_cond%Iskyd_nir, icair)

                call compute_kb(xl, sza_value, phi1, phi2, kb) ! xl = xl(icell)
                call compute_fssh(nbot, ntop, sumpai, kb, fssh)

                ! Direct beam transmittance (tb) and diffuse transmittance (td)
                ! through a single layer
                !
                ! Transmittance (tbcum) of unscattered direct beam onto layer i

                if ( k == 1) then ! soil layer
                   dpai_value   = 0.d0
                   avars(icell)%leaf_tb = 0.d0
                   avars(icell)%leaf_td = 0.d0

                   cumpai_value = cumpai(nbot) ! i-th
                   avars(icell)%leaf_tbcum = exp(-kb * cumpai_value * clump_fac)
                else
                   dpai_value   = dpai  (nbot + k - 2) ! i-th
                   call compute_transmittance_coefficents(xl, sza_value, dpai_value, clump_fac, tb, td)

                   avars(icell)%leaf_tb = tb
                   avars(icell)%leaf_td = td
                   leaf_td(nbot + k - 2) = td

                   if (k == ntop) then
                      cumpai = 0.d0
                   else
                      cumpai_value   = cumpai(nbot + k - 1) ! (i+1)-th
                   end if
                   avars(icell)%leaf_tbcum = exp(-kb * cumpai_value * clump_fac)

                endif

                if (k > 1) then
                   !avars(icell)%leaf_fssh(1) = avars(icell)%leaf_tbcum/(kb * dpai_value) * (1.d0 - exp(-kb * clump_fac * dpai_value))
                   avars(icell)%leaf_fssh(1) = fssh(nbot + k - 2)
                   avars(icell)%leaf_fssh(2) = 1.d0 - avars(icell)%leaf_fssh(1)
                end if

             end do ! k-loop
          end do ! itree-loop
       end do ! icair-loop

       sum_conn = 0
       cur_cond => cur_goveq%boundary_conditions%first
       do
          if (.not.associated(cur_cond)) exit

          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             if (sum_conn > ncair) then
                write(iulog,*)'Num. of cells in BC is greater than ncair'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end if

             cur_goveq%aux_vars_bc(sum_conn)%Iskyb(1) = get_value_from_condition(bnd_cond%Iskyb_vis, sum_conn)
             cur_goveq%aux_vars_bc(sum_conn)%Iskyb(2) = get_value_from_condition(bnd_cond%Iskyb_nir, sum_conn)
             cur_goveq%aux_vars_bc(sum_conn)%Iskyd(1) = get_value_from_condition(bnd_cond%Iskyd_vis, sum_conn)
             cur_goveq%aux_vars_bc(sum_conn)%Iskyd(2) = get_value_from_condition(bnd_cond%Iskyd_nir, sum_conn)

          enddo
          cur_cond => cur_cond%next
       enddo

    end select

  end subroutine swv_set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_number_of_leaves(swv_mpp, nleaves)
    !
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnShortwaveType       , only : goveqn_shortwave_type
    !
    implicit none
    !
    class(mpp_shortwave_type)        :: swv_mpp
    PetscInt                         :: nleaves
    !
    class(goveqn_base_type), pointer :: cur_goveq

    cur_goveq => swv_mpp%soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_shortwave_type)
          cur_goveq%nleaf = nleaves
       end select

       cur_goveq => cur_goveq%next
    end do

  end subroutine set_number_of_leaves

  !------------------------------------------------------------------------
  subroutine init_swv(swv_mpp)
    !
    use ml_model_global_vars, only : dlai, dsai, dpai
    !
    implicit none
    !
    class(mpp_shortwave_type) :: swv_mpp

    call initialize(swv_mpp)

    call setup_meshes(swv_mpp)

    call add_goveqns(swv_mpp)

    call add_conditions_to_goveqns(swv_mpp)

    call set_number_of_leaves(swv_mpp, 2)

    call swv_mpp%AllocateAuxVars()

    call swv_mpp%SetupProblem()

    call set_parameters(swv_mpp, dlai, dsai, dpai)

  end subroutine init_swv

  !------------------------------------------------------------------------
  subroutine extract_data_from_swv(swv_mpp, istep)
    !
    ! !DESCRIPTION:
    !   Extracts following variables from the LBL model:
    !     - Iabs_leaf:
    !         - sunlit + VIS
    !         - sunlit + NIR
    !         - shaded + VIS
    !         - shaded + NIR
    !     - Iabs_soil
    !         - VIS
    !         - NIR
    !
    ! !USES:
    use ml_model_global_vars      , only : nbot, ntop, ncair, ntree, nz_cair
    use ml_model_meshes           , only : nleaf
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnShortwaveType       , only : goveqn_shortwave_type
    use MultiPhysicsProbShortwave , only : mpp_shortwave_type
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_LEAF_ABSORBED_SHORTWAVE_RAD_PER_LAI
    use MultiPhysicsProbConstants , only : VAR_SOIL_ABSORBED_SHORTWAVE_RAD_PER_GROUND
    use ml_model_utils            , only : set_value_in_condition
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_shortwave_type)           :: swv_mpp
    PetscInt                            :: istep
    !
    PetscInt                            :: idx_leaf, idx_data, idx_soil, idx_air
    PetscInt                            :: ileaf, icair, itree, k, ieqn, icell, iband
    PetscInt                            :: ncells, nctz, count
    PetscReal               , pointer   :: Iabs_leaf(:), Iabs_soil(:)
    class(goveqn_base_type) , pointer   :: goveq
    PetscInt                , parameter :: nband = 2      ! Visible + NIR
    PetscErrorCode                      :: ierr
    character(len=20)                   :: step_string

    nctz   = ncair * ntree * (ntop - nbot + 1 + 1) ! num of canopy airspace x num. of tree x num. of levels
    ncells = nctz  * nleaf * nband
    allocate(Iabs_leaf(ncells))

    ncells = nctz * nband
    allocate(Iabs_soil(ncells))

    call swv_mpp%soe%SetPointerToIthGovEqn(1, goveq)

    select type(goveq)
    class is (goveqn_shortwave_type)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_ABSORBED_SHORTWAVE_RAD_PER_LAI, nctz, Iabs_leaf)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_SOIL_ABSORBED_SHORTWAVE_RAD_PER_GROUND, nctz, Iabs_soil)
    end select

    if (output_data) then
       write(step_string,*)istep
       write(*,*)'mpp.iabs{' // trim(adjustl(step_string)) // '} = ['
    end if
    count = 0
    icell = 0
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, ntop - nbot + 1 + 1
             if (k > 1) then
                icell = icell + 1
                count = count + 1; call set_value_in_condition(int_cond%Ileaf_sun_vis, icell, Iabs_leaf(count))
                count = count + 1; call set_value_in_condition(int_cond%Ileaf_shd_vis, icell, Iabs_leaf(count))
                count = count + 1; call set_value_in_condition(int_cond%Ileaf_sun_nir, icell, Iabs_leaf(count))
                count = count + 1; call set_value_in_condition(int_cond%Ileaf_shd_nir, icell, Iabs_leaf(count))
                if (output_data) then
                   write(*,*)icell, Iabs_leaf(count-3), Iabs_leaf(count-2), Iabs_leaf(count-1), Iabs_leaf(count)
                end if
             else
                count = count + 4
             end if
          end do
       end do
    end do
    if (output_data) then
       write(*,*)'];'
    end if

    count = 0
    icell = 0
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, ntop - nbot + 1 + 1
                if (k == 1) then
                   icell = icell + 1;

                   count = count + 1; call set_value_in_condition(int_cond%Isoil_vis, icell, Iabs_soil(count))
                   count = count + 1; call set_value_in_condition(int_cond%Isoil_nir, icell, Iabs_soil(count))
                else
                   count = count + 2
                end if
          end do
       end do
    end do

    deallocate(Iabs_leaf)
    deallocate(Iabs_soil)

  end subroutine extract_data_from_swv
 
  !------------------------------------------------------------------------
  subroutine solve_swv(swv_mpp, istep, dt)
    !
    implicit none
    !
    class(mpp_shortwave_type) :: swv_mpp
    PetscInt                  :: istep
    PetscReal                 :: dt
    !
    PetscBool                 :: converged
    PetscInt                  :: converged_reason
    PetscBool                 :: flg
    PetscErrorCode            :: ierr

    call swv_set_boundary_conditions(swv_mpp)

    call swv_mpp%soe%StepDT(dt, istep, converged, converged_reason, ierr)

    call extract_data_from_swv(swv_mpp, istep)

  end subroutine solve_swv

end module swv
