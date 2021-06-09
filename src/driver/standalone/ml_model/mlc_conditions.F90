module mlc_conditions

#include <petsc/finclude/petsc.h>

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use ml_model_global_vars
  use petscsys

  implicit none

  public :: mlc_add_conditions_to_goveqns
  public :: mlc_set_boundary_conditions

contains

  !------------------------------------------------------------------------
  subroutine mlc_add_conditions_to_goveqns(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS
    !
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp

    call add_non_coupling_conditions_to_goveqns(mlc_mpp)
    call add_internal_coupling_vars(mlc_mpp)

  end subroutine mlc_add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine add_non_coupling_conditions_to_goveqns(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use ConnectionSetType         , only : connection_set_type
    use MeshType                  , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants , only : CONN_VERTICAL
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type)                  :: mlc_mpp
    !
    PetscInt                             :: kk , nconn
    PetscInt                   , pointer :: id_up(:)
    PetscInt                   , pointer :: id_dn(:)
    PetscReal                  , pointer :: dist_up(:)
    PetscReal                  , pointer :: dist_dn(:)
    PetscReal                  , pointer :: area(:)
    PetscInt                   , pointer :: itype(:)
    PetscReal                  , pointer :: unit_vec(:,:)
    PetscReal                            :: dz_cair
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
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_WATER_VAPOR
    use MultiPhysicsProbConstants , only : VAR_LEAF_TEMPERATURE
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp
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


   !------------------------------------------------------------------------
  subroutine mlc_set_boundary_conditions(mlc_mpp , &
       gs_sun, rn_sun                                 , &
       gs_shd, rn_shd                                 , &
       gbh, gbv                                       , &
       Pref, Uref, Rhref, Tref, Tcan                  , &
       rn_soil, tg)

    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type)  :: mlc_mpp
    PetscReal , pointer :: gs_sun(:), gs_shd(:)
    PetscReal , pointer :: rn_sun(:), rn_shd(:)
    PetscReal , pointer :: gbh(:), gbv(:)
    PetscReal , pointer :: Pref(:), Uref(:), Rhref(:), Tref(:), Tcan(:)
    PetscReal , pointer :: rn_soil(:), tg(:)

    call set_air_temp_boundary_conditions( mlc_mpp, gbh, gs_sun, gs_shd)
    call set_air_vapor_boundary_conditions(mlc_mpp, gbv, gs_sun, gs_shd)

    call set_canopy_leaf_boundary_conditions(mlc_mpp, CLEF_TEMP_SUN_GE, gbh, gbv, gs_sun, rn_sun)
    call set_canopy_leaf_boundary_conditions(mlc_mpp, CLEF_TEMP_SHD_GE, gbh, gbv, gs_shd, rn_shd)

    call set_atmospheric_boundary_conditions(mlc_mpp, Pref, Uref, Rhref, Tref, Tcan)
    call set_soil_boundary_conditions(mlc_mpp, tg, rn_soil)
    
  end subroutine mlc_set_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_air_temp_boundary_conditions(mlc_mpp, gbh, gs_sun, gs_shd)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType , only : goveqn_cair_temp_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp
    PetscReal, pointer :: gbh(:) 
    PetscReal, pointer :: gs_sun(:), gs_shd(:)
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, ileaf, icair, icell

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CAIR_TEMP_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_cair_temp_type)


       do icair = 1, ncair
          do k = 1, nz_cair
             icell = (icair-1)*(nz_cair+1) + k
             cur_goveq%aux_vars_in(icell)%gbh  = gbh(icell)

             do ileaf = 1, ntree
                cur_goveq%aux_vars_in(icell)%leaf_gs(        ileaf) = gs_sun(k)
                cur_goveq%aux_vars_in(icell)%leaf_gs(ntree + ileaf) = gs_shd(k)
             enddo
          end do
       enddo

    end select

  end subroutine set_air_temp_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_air_vapor_boundary_conditions(mlc_mpp, gbv, gs_sun, gs_shd)
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
    class(mpp_mlc_type) :: mlc_mpp
    PetscReal                  , pointer :: gbv(:)
    PetscReal                  , pointer :: gs_sun(:), gs_shd(:)
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, ileaf, icair, icell

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CAIR_VAPR_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_cair_vapor_type)

      do icair = 1, ncair
         do k = 1, nz_cair
            icell = (icair-1)*(nz_cair+1) + k
            cur_goveq%aux_vars_in(icell)%gbv  = gbv(icell)

            do ileaf = 1, ntree
               cur_goveq%aux_vars_in(icell)%leaf_gs(        ileaf) = gs_sun(k)
               cur_goveq%aux_vars_in(icell)%leaf_gs(ntree + ileaf) = gs_shd(k)

            end do
         end do
      end do
    end select

  end subroutine set_air_vapor_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_canopy_leaf_boundary_conditions(mlc_mpp, ge_rank, gbh, gbv, gs, rn)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type) :: mlc_mpp
    PetscInt           :: ge_rank
    PetscReal, pointer :: gbh(:), gbv(:)
    PetscReal, pointer :: gs(:), rn(:)
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, i, num_int, icell, icair, itree

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(ge_rank, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_cleaf_temp_type)


       do icair = 1, ncair
          do itree = 1, ntree
             do k = 1, nz_cair + 1
                icell = (icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cair+1) + k
                cur_goveq%aux_vars_in(icell)%gs = gs(icell)
                cur_goveq%aux_vars_in(icell)%rn = rn(icell)

                icell = (icair-1)*(nz_cair+1) + k
                cur_goveq%aux_vars_in(icell)%gbh = gbh(icell)
                cur_goveq%aux_vars_in(icell)%gbv = gbv(icell)
             end do
          end do
       end do

    end select

  end subroutine set_canopy_leaf_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_atmospheric_boundary_conditions(mlc_mpp, Pref, Uref, Rhref, Tref, Tcan)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsMLCType  , only : sysofeqns_mlc_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type)                  :: mlc_mpp
    PetscReal                  , pointer :: Pref(:), Uref(:), Rhref(:), Tref(:), Tcan(:)
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscInt                             :: icair

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do icair = 1, ncair

       soe%cturb%pref(icair)  = Pref(icair)
       soe%cturb%uref(icair)  = Uref(icair)
       soe%cturb%tref(icair)  = Tref(icair)
       soe%cturb%rhref(icair) = Rhref(icair)
       soe%cturb%tcan(icair)  = Tcan(icair)

       call soe%cturb%ComputeDerivedAtmInputs(icair)

       soe%cturb%qcan(icair) = soe%cturb%qref(icair)
    end do

  end subroutine set_atmospheric_boundary_conditions

  !------------------------------------------------------------------------
  subroutine set_soil_boundary_conditions(mlc_mpp, tg, rn_soil)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use SystemOfEquationsMLCType  , only : sysofeqns_mlc_type
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type)                   :: mlc_mpp
    PetscReal                  , pointer :: rn_soil(:), tg(:)
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscInt                             :: icair

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do icair = 1, ncair

       soe%cturb%tsoi(icair)   = tg(icair)
       soe%cturb%rnsoi(icair)  = rn_soil(icair)
    end do

  end subroutine set_soil_boundary_conditions

end module
