module ml_model_utils
  !
  use mpp_varctl           , only : iulog
  use mpp_abortutils       , only : endrun
  use mpp_shr_log_mod      , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbLBL  , only : mpp_lbl_type
  !
  implicit none
  !
#include <petsc/finclude/petsc.h>
  PetscInt, parameter :: CANOPY_AIR_MESH      = 1
  PetscInt, parameter :: CANOPY_MESH          = 2
  PetscInt, parameter :: CANOPY_AND_SOIL_MESH = 3

  public :: compute_dpai_fssh
  !public :: extract_data_from_mlc
  public :: extract_data_from_swv
  public :: extract_data_from_lbl
  public :: extract_data_from_lwv
  public :: extract_data_from_photosynthesis
  public :: set_value_in_condition
  public :: get_value_from_condition

contains

  !------------------------------------------------------------------------
  subroutine compute_dpai_fssh()
    !
    use ml_model_global_vars , only : hc, nveg, nbot, ntop, nz_cair, ncair, dz_cair, dpai, dpai, fssh, cumlai
    !
    implicit none
    !
    PetscInt  :: k, i, num_int, ic_bot
    PetscReal :: Kb, sumpai
    PetscReal :: dz_leaf, qbeta, pbeta, pai
    PetscReal :: zl, zu, z_int, dz_int, zrel, beta_pdf, pad
    PetscReal :: pai_sum, pai_miss, pai_old, pai_new

    allocate(dpai(nz_cair+1))
    allocate(fssh(nz_cair+1))
    allocate(cumlai(nz_cair+1))

    dpai(:) = 0.d0
    fssh(:) = 0.d0
    cumlai(:) = 0.d0

    ! Determine plant area index increment for each layer by numerically
    ! integrating the plant area density (beta distribution) between
    ! the bottom and top heights for that layer

    dz_leaf = dz_cair;
    pbeta   = 3.5d0;              ! Parameter for beta distribution
    qbeta   = 2.d0;               ! Parameter for beta distribution
    pai     = 5.051612734794617d0

    nbot = 6 + 1     ! 1st layer is soil layer
    ntop = nveg + 1  ! 

    do k = nbot, ntop
       zl = dz_leaf * (k-2);
       zu = dz_leaf * (k-1);
       dpai(k) = 0.d0

       ! Numerical integration between zl and zu using 100 sublayers
       num_int = 100
       dz_int = (zu-zl)/num_int
       do i = 1, num_int
          if (i == 1) then
             z_int = zl + 0.5d0 * dz_int
          else
             z_int = z_int + dz_int
          end if

          zrel = min(z_int/hc,1.d0)
          beta_pdf = (zrel**(pbeta-1) * (1.d0 - zrel)**(qbeta-1))/exp(lgamma(pbeta) + lgamma(qbeta) - lgamma(pbeta+qbeta));
          pad = (pai / hc) * beta_pdf

          dpai(k) = dpai(k) + pad * dz_int

       end do
    end do

    pai_sum = 0;
    do k = nbot, ntop
       pai_sum = pai_sum + dpai(k)
    end do

    ! Set layers with small plant area index to zero
    pai_miss = 0;
    do k = nbot, ntop
       if (dpai(k) < 0.01d0) then
          pai_miss = pai_miss + dpai(k)
          dpai(k) = 0.d0;
       end if
    end do

    ! Distribute the missing plant area across vegetation layers
    ! in proportion to the plant area profile

    if (pai_miss > 0.d0) then
       pai_old = pai_sum;
       pai_new = pai_old - pai_miss;
       do k = nbot, ntop
          dpai(k) = dpai(k) + pai_miss * (dpai(k) / pai_new);
       end do
    end if

    ! Find the lowest vegetation layer
    do k = ntop, nbot,-1
       if (dpai(k) > 0.d0) then
          ic_bot = k
       end if
    end do
    nbot = ic_bot

    Kb = 1.762817445019839d0;
    do k = ntop, nbot, -1
       if (k == ntop) then
          sumpai = 0.5d0 * dpai(k);
       else
          sumpai = sumpai + &
               0.5d0 * (dpai(k+1) + dpai(k));
       end if
       fssh(k) = exp(-Kb * sumpai);
       cumlai(k) = sumpai
    end do

  end subroutine compute_dpai_fssh

  !------------------------------------------------------------------------
  subroutine allocate_memory_for_condition(cond, ndata)
    !
    use ml_model_global_vars, only : condition_type
    !
    implicit none
    !
    type(condition_type) :: cond
    PetscInt             :: ndata

    cond%ndata = ndata
    allocate(cond%data(ndata))
    cond%data(:) = -999.d0

  end subroutine allocate_memory_for_condition  

  !------------------------------------------------------------------------
  subroutine set_value_in_condition(cond, idata, value)
    !
    use ml_model_global_vars, only : condition_type
    !
    implicit none
    !
    type(condition_type) :: cond
    PetscInt             :: idata
    PetscReal            :: value

    if (idata > cond%ndata) then
       write(*,*)'idata = ',idata,'cond%ndata = ',cond%ndata
       call endrun(msg=' ERROR: Attempting to set data for a index that exceeds size'//&
            errMsg(__FILE__, __LINE__))
    end if

    cond%data(idata) = value

  end subroutine set_value_in_condition

  !------------------------------------------------------------------------
  function get_value_from_condition(cond, idata) result (value)
    !
    use ml_model_global_vars, only : condition_type
    !
    implicit none
    !
    type(condition_type) :: cond
    PetscInt             :: idata
    !
    PetscReal            :: value

    if (idata > cond%ndata) then
       call endrun(msg=' ERROR: Attempting to get data for a index that exceeds size'//&
            errMsg(__FILE__, __LINE__))
    end if

    value = cond%data(idata)

  end function get_value_from_condition

  !------------------------------------------------------------------------
  subroutine extract_data_from_lbl(lbl_mpp)
    !
    ! !DESCRIPTION:
    !   Extracts following variables from the LBL model:
    !     - gbh
    !     - gbv
    !     - gbc
    !
    ! !USES:
    use ml_model_global_vars            , only : nbot, ntop, ncair, ntree, nz_cair
    use ml_model_global_vars            , only : gbh, gbv, gbc
    use ml_model_meshes                 , only : nleaf
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnLeafBoundaryLayer         , only : goveqn_leaf_bnd_lyr_type
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_HEAT
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_H2O
    use MultiPhysicsProbConstants       , only : VAR_LEAF_BDN_LYR_COND_CO2
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_lbl_type)               :: lbl_mpp
    !
    PetscScalar             , pointer :: soln_p(:)
    PetscInt                          :: idx_leaf, idx_data, idx_soil, idx_air
    PetscInt                          :: ileaf, icair, itree, k, ieqn, icell
    PetscInt                          :: ncells
    PetscReal               , pointer :: gbh_data(:), gbv_data(:), gbc_data(:)
    class(goveqn_base_type) , pointer :: goveq
    PetscErrorCode                    :: ierr

    ncells = ncair * ntree * (ntop - nbot + 1) * nleaf

    allocate(gbh_data(ncells))
    allocate(gbv_data(ncells))
    allocate(gbc_data(ncells))

    call lbl_mpp%soe%SetPointerToIthGovEqn(1, goveq)

    select type(goveq)
    class is (goveqn_leaf_bnd_lyr_type)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_HEAT, ncells, gbh_data)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_H2O , ncells, gbv_data)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_BDN_LYR_COND_CO2 , ncells, gbc_data)
    end select

    do icell = 1, ncells
       call set_value_in_condition(gbh, icell, gbh_data(icell))
       call set_value_in_condition(gbv, icell, gbv_data(icell))
       call set_value_in_condition(gbc, icell, gbc_data(icell))
    end do

    deallocate(gbh_data)
    deallocate(gbv_data)
    deallocate(gbc_data)

  end subroutine extract_data_from_lbl
 
  !------------------------------------------------------------------------
  subroutine extract_data_from_swv(swv_mpp)
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
    use ml_model_global_vars      , only : Ileaf_sun_vis, Ileaf_shd_vis
    use ml_model_global_vars      , only : Ileaf_sun_nir, Ileaf_shd_nir
    use ml_model_global_vars      , only : Isoil_vis, Isoil_nir
    use ml_model_meshes           , only : nleaf
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnShortwaveType       , only : goveqn_shortwave_type
    use MultiPhysicsProbShortwave , only : mpp_shortwave_type
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_LEAF_ABSORBED_SHORTWAVE_RAD_PER_LAI
    use MultiPhysicsProbConstants , only : VAR_SOIL_ABSORBED_SHORTWAVE_RAD_PER_GROUND
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_shortwave_type)           :: swv_mpp
    !
    PetscInt                            :: idx_leaf, idx_data, idx_soil, idx_air
    PetscInt                            :: ileaf, icair, itree, k, ieqn, icell, iband
    PetscInt                            :: ncells, nctz, count
    PetscReal               , pointer   :: Iabs_leaf(:), Iabs_soil(:)
    class(goveqn_base_type) , pointer   :: goveq
    PetscInt                , parameter :: nband = 2      ! Visible + NIR
    PetscErrorCode                      :: ierr

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

    count = 0
    icell = 0
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, ntop - nbot + 1 + 1
             if (k > 1) then
                icell = icell + 1
                count = count + 1; call set_value_in_condition(Ileaf_sun_vis, icell, Iabs_leaf(count))
                count = count + 1; call set_value_in_condition(Ileaf_shd_vis, icell, Iabs_leaf(count))
                count = count + 1; call set_value_in_condition(Ileaf_sun_nir, icell, Iabs_leaf(count))
                count = count + 1; call set_value_in_condition(Ileaf_shd_nir, icell, Iabs_leaf(count))
             else
                count = count + 4
             end if
          end do
       end do
    end do

    count = 0
    icell = 0
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, ntop - nbot + 1 + 1
                if (k == 1) then
                   icell = icell + 1;

                   count = count + 1; call set_value_in_condition(Isoil_vis, icell, Iabs_soil(count))
                   count = count + 1; call set_value_in_condition(Isoil_vis, icell, Iabs_soil(count))
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
  subroutine extract_data_from_lwv(lwv_mpp)
    !
    ! !DESCRIPTION:
    !   Extracts following variables from the longwave radiation model:
    !     - Iabs_leaf:
    !     - Iabs_soil
    !
    ! !USES:
    use ml_model_global_vars      , only : nbot, ntop, ncair, ntree, nz_cair
    use ml_model_global_vars      , only : Labs_leaf_sun, Labs_leaf_shd, Labs_soil
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnLongwaveType        , only : goveqn_longwave_type
    use MultiPhysicsProbLongwave  , only : mpp_longwave_type
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_LEAF_ABSORBED_LONGWAVE_RAD_PER_LAI
    use MultiPhysicsProbConstants , only : VAR_SOIL_ABSORBED_LONGWAVE_RAD_PER_GROUND
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_longwave_type)           :: lwv_mpp
    !
    PetscInt                            :: idx_leaf, idx_data, idx_soil, idx_air
    PetscInt                            :: icair, itree, k, leaf_icell, soil_icell
    PetscInt                            :: ncells, count
    PetscReal               , pointer   :: Labs_leaf_data(:), Labs_soil_data(:)
    class(goveqn_base_type) , pointer   :: goveq
    PetscErrorCode                      :: ierr

    ncells = ncair * ntree * (ntop - nbot + 1 + 1) ! num of canopy airspace x num. of tree x num. of levels
    allocate(Labs_leaf_data(ncells))
    allocate(Labs_soil_data(ncells))

    call lwv_mpp%soe%SetPointerToIthGovEqn(1, goveq)

    select type(goveq)
    class is (goveqn_longwave_type)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_ABSORBED_LONGWAVE_RAD_PER_LAI   , ncells, Labs_leaf_data)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_SOIL_ABSORBED_LONGWAVE_RAD_PER_GROUND, ncells, Labs_soil_data)
    end select

    count = 0
    leaf_icell = 0
    soil_icell = 0
    do icair = 1, ncair
       do itree = 1, ntree
          do k = 1, ntop - nbot + 1 + 1
             count = count + 1;
             if (k == 1) then
                soil_icell = soil_icell + 1
                call set_value_in_condition(Labs_soil, soil_icell, Labs_soil_data(count))
             else
                leaf_icell = leaf_icell + 1
                call set_value_in_condition(Labs_leaf_sun, leaf_icell, Labs_leaf_data(count))
                call set_value_in_condition(Labs_leaf_shd, leaf_icell, Labs_leaf_data(count))
            end if
          end do
       end do
    end do

    deallocate(Labs_leaf_data)
    deallocate(Labs_soil_data)

  end subroutine extract_data_from_lwv

  !------------------------------------------------------------------------
  subroutine extract_data_from_photosynthesis(psy_mpp)
    !
    ! !DESCRIPTION:
    !   Extracts following variables from the photsynthesis model
    !     - gs
    !
    ! !USES:
    use ml_model_global_vars      , only : nbot, ntop, ncair, ntree, nz_cair
    use ml_model_meshes           , only : nleaf
    use ml_model_global_vars      , only : gs_sun, gs_shd
    use GoverningEquationBaseType , only : goveqn_base_type
    use GoveqnPhotosynthesisType  , only : goveqn_photosynthesis_type
    use MultiPhysicsProbPhotosynthesis  , only : mpp_photosynthesis_type
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_STOMATAL_CONDUCTANCE
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_photosynthesis_type)      :: psy_mpp
    !
    PetscInt                            :: icell, offset
    PetscInt                            :: ncells
    PetscReal               , pointer   :: gs_data(:)
    class(goveqn_base_type) , pointer   :: goveq
    PetscErrorCode                      :: ierr

    ncells = ncair * ntree * (ntop - nbot + 1) * nleaf

    allocate(gs_data(ncells))

    call psy_mpp%soe%SetPointerToIthGovEqn(1, goveq)

    select type(goveq)
    class is (goveqn_photosynthesis_type)
       call goveq%GetRValues(AUXVAR_INTERNAL, VAR_STOMATAL_CONDUCTANCE, ncells, gs_data)
    end select

    offset = ncair * ntree * (ntop - nbot + 1)
    do icell = 1, ncair * ntree * (ntop - nbot + 1)
       call set_value_in_condition(gs_sun, icell, gs_data(icell         ))
       call set_value_in_condition(gs_shd, icell, gs_data(icell + offset))
    end do

    deallocate(gs_data)

  end subroutine extract_data_from_photosynthesis

end module ml_model_utils
