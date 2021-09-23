module photosynthesis_parameters

  use mpp_varctl                     , only : iulog
  use mpp_abortutils                 , only : endrun
  use mpp_shr_log_mod                , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbPhotosynthesis , only : mpp_photosynthesis_type
  use WaterVaporMod                  , only : SatVap
  use photosynthesis_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: set_parameters

contains

  !------------------------------------------------------------------------
  subroutine set_parameters(phtsyn_mpp)

    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use SystemOfEquationsPhotosynthesisType , only : sysofeqns_photosynthesis_type
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnPhotosynthesisType            , only : goveqn_photosynthesis_type
    use PhotosynthesisAuxType               , only : photosynthesis_auxvar_type
    use MultiPhysicsProbConstants           , only : VAR_PHOTOSYNTHETIC_PATHWAY_C3, VAR_PHOTOSYNTHETIC_PATHWAY_C4
    use MultiPhysicsProbConstants           , only : VAR_SCM_MEDLYN, VAR_SCM_BBERRY
    use MultiPhysicsProbConstants           , only : VAR_SCM_BONAN14
    use MultiPhysicsProbConstants           , only : VAR_SCM_WUE
    use MultiPhysicsProbConstants           , only : TFRZ
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_photosynthesis_type)                         :: phtsyn_mpp

    class(sysofeqns_base_type)                , pointer   :: base_soe
    class(sysofeqns_photosynthesis_type)      , pointer   :: soe
    class(goveqn_base_type)                   , pointer   :: cur_goveq
    type(photosynthesis_auxvar_type)          , pointer   :: auxvar
    PetscReal                                 , parameter :: tau = 0.1d0
    PetscReal                                 , parameter :: rho = 0.1d0
    PetscReal                                             :: tair, esat_tair, desat_tair, relhum, vpd_tleaf
    PetscReal :: esat_25C, desat_25C, eair
    PetscReal :: esat_current, desat_current
    PetscInt                                              :: k , icair, icell
    PetscBool                              :: bounded

    call phtsyn_mpp%soe%SetPointerToIthGovEqn(PHTSYN_GE, cur_goveq)

    call SatVap (273.15d0 + 25.d0, esat_25C, desat_25C)

    tair = TFRZ + 25.d0
    relhum = 80.d0
    call SatVap (tair, esat_tair, desat_tair)
    eair = esat_tair * relhum /100.d0
    vpd_tleaf = esat_25C - eair
   
    select type(cur_goveq)
    class is (goveqn_photosynthesis_type)

       do icair = 1, ncair
          do k = 1, nz_cair+1
             icell = (icair-1)*(nz_cair+1) + k

             auxvar => cur_goveq%aux_vars_in(icell)

             cur_goveq%aux_vars_in(icell)%colim  = 1
             cur_goveq%aux_vars_in(icell)%c3psn  = c3psn
             cur_goveq%aux_vars_in(icell)%gstype = gstype

             cur_goveq%aux_vars_in(icell)%cair  = 380.d0                       ! [mmol/mol]
             cur_goveq%aux_vars_in(icell)%o2ref = 0.209d0 * 1000.d0            ! [mmol/mol]
             cur_goveq%aux_vars_in(icell)%apar  = 2000.d0 * ( 1.d0 - rho - tau) ! [mmol photon/m2/s]
             
             cur_goveq%aux_vars_in(icell)%tleaf = TFRZ + 11.d0 + 0.25d0 * (k - 1)

             ! Values at tleaf = 11C and tair = 25C
             cur_goveq%aux_vars_in(icell)%gbv   = 2.224407920268566d0
             cur_goveq%aux_vars_in(icell)%gbc   = 1.637448199187622d0

             ! Values at tleaf = 25C and tair = 25C
             !cur_goveq%aux_vars_in(icell)%gbv   = 2.304789280798361d0
             !cur_goveq%aux_vars_in(icell)%gbc   = 1.694493818233533d0

             call SatVap (cur_goveq%aux_vars_in(icell)%tleaf, esat_current, desat_current)

             if (cur_goveq%aux_vars_in(icell)%gstype == VAR_SCM_MEDLYN .or. &
                 cur_goveq%aux_vars_in(icell)%gstype == VAR_SCM_WUE) then
                cur_goveq%aux_vars_in(icell)%eair = esat_current - vpd_tleaf
             else
                cur_goveq%aux_vars_in(icell)%eair = esat_current * relhum/100.d0
             end if

             cur_goveq%aux_vars_in(icell)%btran = 1.d0
             cur_goveq%aux_vars_in(icell)%dpai  = 1.d0

             if (cur_goveq%aux_vars_in(icell)%gstype == VAR_SCM_WUE) then
                ! Set cell active/inactive if the solution is bounded/unbounded
                call cur_goveq%aux_vars_in(icell)%IsWUESolutionBounded(bounded)
                cur_goveq%mesh%is_active(icell) = bounded
             elseif (cur_goveq%aux_vars_in(icell)%gstype == VAR_SCM_BONAN14) then
                call cur_goveq%aux_vars_in(icell)%DetermineIfSolutionIsBounded()
             end if

             ! Root parameters
             auxvar%root%biomass = 500.d0
             auxvar%root%radius  = 0.29d-3
             auxvar%root%density = 0.31d6
             auxvar%root%resist  = 25.d0

             ! Soil paramters
             call set_soil_parameters(auxvar%soil)

             ! Plant parameters
             auxvar%plant%nleaf = 1
             call auxvar%plant%AllocateMemory()

             auxvar%plant%leaf_psi(:)    = -1.5d0
             auxvar%plant%leaf_height(:) = 15.d0
             auxvar%plant%leaf_capc(:)   = 2500.d0
             auxvar%plant%leaf_minlwp(:) = -2.d0
             auxvar%plant%leaf_lai(:)    = 500.d0
             auxvar%plant%k_stem2leaf(:) = 4.d0

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
    PetscInt  :: j, texture
    PetscReal :: beta_param, z1, z2
    PetscReal, parameter, dimension(11) :: theta_sat = [0.395d0, 0.410d0, 0.435d0, 0.485d0, 0.451d0, 0.420d0, 0.477d0, 0.476d0, 0.426d0, 0.492d0, 0.482d0]
    PetscReal, parameter, dimension(11) :: psi_sat   = [-121.d0, -90.d0, -218d0, -786.d0, -478.d0, -299.d0, -356.d0, -630.d0, -153.d0, -490.d0, -405.d0]
    PetscReal, parameter, dimension(11) :: b         = [4.05d0, 4.38d0, 4.90d0, 5.30d0, 5.39d0, 7.12d0, 7.75d0, 8.52d0, 10.40d0, 10.40d0, 11.40d0];
    PetscReal, parameter, dimension(11) :: k_sat     = [1.056d0, 0.938d0, 0.208d0, 0.0432d0, 0.0417d0, 0.0378d0, 0.0102d0, 0.0147d0, 0.0130d0, 0.0062d0, 0.0077d0];

    texture = 5;

    soil%nlevsoi = 11

    call soil%AllocateMemory()

    soil%dz(01) = 0.050
    soil%dz(02) = 0.050
    soil%dz(03) = 0.100
    soil%dz(04) = 0.100
    soil%dz(05) = 0.200
    soil%dz(06) = 0.200
    soil%dz(07) = 0.200
    soil%dz(08) = 0.300
    soil%dz(09) = 0.400
    soil%dz(10) = 0.400
    soil%dz(11) = 0.500

    beta_param = 0.90;  ! root profile parameter: shallow profile
    !beta_param = 0.97;  ! root profile parameter: deep profile

    do j = 1, soil%nlevsoi
       if (j == 1) then
          z2 = soil%dz(j) * 100;
          soil%rootfr(j) = 1 - beta_param**z2;
       else
          z1 = z2;
          z2 = z1 + soil%dz(j) * 100;
          soil%rootfr(j) = beta_param**z1 - beta_param**z2;
       end if

       soil%watsat(j)     = theta_sat(texture)
       soil%hksat(j)      = k_sat(texture) * 10.d0/60.d0
       soil%bsw(j)        = b(texture)

       soil%h2osoi_vol(j) = 0.5*soil%watsat(j)
       soil%psi_sat(j)    = psi_sat(texture)
       soil%psi(j)        = psi_sat(texture) * (soil%h2osoi_vol(j)/soil%watsat(j))**(-soil%bsw(j))
    end do

  end subroutine set_soil_parameters

end module photosynthesis_parameters
