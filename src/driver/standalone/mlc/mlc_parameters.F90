module mlc_parameters

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use mlc_global_vars
  use petscsys

  implicit none

#include <petsc/finclude/petsc.h>

  public :: set_parameters

contains

  !------------------------------------------------------------------------
  subroutine set_parameters(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_mlc_type) :: mlc_mpp

    call set_air_temp_ge_parameters(mlc_mpp)
    call set_air_vapor_ge_parameters(mlc_mpp)
    call set_common_canopy_leaf_parameters(mlc_mpp, CLEF_TEMP_SUN_GE)
    call set_common_canopy_leaf_parameters(mlc_mpp, CLEF_TEMP_SHD_GE)
    call set_sunlit_canopy_parameters(mlc_mpp)
    call set_shaded_canopy_parameters(mlc_mpp)
    call set_turbulence_parameters(mlc_mpp)
    call set_soil_parameters(mlc_mpp)

  end subroutine set_parameters

  !------------------------------------------------------------------------
  subroutine set_air_temp_ge_parameters(mlc_mpp)
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
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, ileaf, icair, icell
    PetscReal                  , pointer :: dpai(:), fssh(:), gs_sun(:), gs_shd(:)

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CAIR_TEMP_GE, cur_goveq)

    allocate(dpai(nz_cleaf+1))
    allocate(fssh(nz_cleaf+1))
    allocate(gs_sun(nz_cleaf+1))
    allocate(gs_shd(nz_cleaf+1))

    select type(cur_goveq)
    class is (goveqn_cair_temp_type)

      call get_dpai_fssh(dpai, fssh)
      call get_sunlit_canopy_gs(gs_sun)
      call get_shaded_canopy_gs(gs_shd)

      do icair = 1, ncair
         do k = 1, nz_cair
            icell = (icair-1)*(nz_cair+1) + k
            cur_goveq%aux_vars_in(icell)%gbh  = 2.268731551029694d0

            cur_goveq%aux_vars_in(icell)%leaf_dpai(:) = dpai(k)/ntree
            cur_goveq%aux_vars_in(icell)%leaf_fwet(:) = 0.d0
            cur_goveq%aux_vars_in(icell)%leaf_fdry(:) = 0.8218390792391702d0

            do ileaf = 1, ntree
               cur_goveq%aux_vars_in(icell)%leaf_gs(        ileaf) = gs_sun(k)
               cur_goveq%aux_vars_in(icell)%leaf_gs(ntree + ileaf) = gs_shd(k)
               cur_goveq%aux_vars_in(icell)%leaf_fssh(        ileaf) = fssh(k)
               cur_goveq%aux_vars_in(icell)%leaf_fssh(ntree + ileaf) = 1.d0 - fssh(k)
            enddo
         end do
         cur_goveq%aux_vars_in((nz_cair+1)*(icair-1)+1)%is_soil = PETSC_TRUE
      enddo

    end select

    deallocate(dpai)
    deallocate(fssh)
    deallocate(gs_sun)
    deallocate(gs_shd)

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
    PetscInt                             :: k, ileaf, icair, icell
    PetscReal                  , pointer :: dpai(:), fssh(:), gs_sun(:), gs_shd(:)

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CAIR_VAPR_GE, cur_goveq)

    allocate(dpai(nz_cleaf+1))
    allocate(fssh(nz_cleaf+1))
    allocate(gs_sun(nz_cleaf+1))
    allocate(gs_shd(nz_cleaf+1))

    select type(cur_goveq)
    class is (goveqn_cair_vapor_type)

      call get_dpai_fssh(dpai, fssh)
      call get_sunlit_canopy_gs(gs_sun)
      call get_shaded_canopy_gs(gs_shd)

      do icair = 1, ncair
         do k = 1, nz_cair
            icell = (icair-1)*(nz_cair+1) + k
            cur_goveq%aux_vars_in(icell)%gbv  = 2.496430918408511d0

            cur_goveq%aux_vars_in(icell)%leaf_dpai(:) = dpai(k)/ntree
            cur_goveq%aux_vars_in(icell)%leaf_fwet(:) = 0.d0
            cur_goveq%aux_vars_in(icell)%leaf_fdry(:) = 0.8218390792391702d0

            do ileaf = 1, ntree
               cur_goveq%aux_vars_in(icell)%leaf_gs(        ileaf) = gs_sun(k)
               cur_goveq%aux_vars_in(icell)%leaf_gs(ntree + ileaf) = gs_shd(k)

               cur_goveq%aux_vars_in(icell)%leaf_fssh(        ileaf) = fssh(k)
               cur_goveq%aux_vars_in(icell)%leaf_fssh(ntree + ileaf) = 1.d0 - fssh(k)
            end do
         end do
         cur_goveq%aux_vars_in((nz_cair+1)*(icair-1)+1)%is_soil = PETSC_TRUE
      end do
    end select

    deallocate(dpai)
    deallocate(fssh)
    deallocate(gs_sun)
    deallocate(gs_shd)

  end subroutine set_air_vapor_ge_parameters

  !------------------------------------------------------------------------
  subroutine set_common_canopy_leaf_parameters(mlc_mpp, ge_rank)
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
    type(mpp_mlc_type) :: mlc_mpp
    PetscInt           :: ge_rank
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt                             :: k, i, num_int, icell, icair, itree
    PetscReal                  , pointer :: dpai(:),fssh(:)

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(ge_rank, cur_goveq)

    allocate(dpai(nz_cleaf+1))
    allocate(fssh(nz_cleaf+1))

    select type(cur_goveq)
    class is (goveqn_cleaf_temp_type)

       call get_dpai_fssh(dpai, fssh)

       do icair = 1, ncair
          do itree = 1, ntree
          do k = 1, nz_cleaf + 1
             icell = (icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k
             cur_goveq%aux_vars_in(icell)%gbh  = 2.268731551029694d0
             cur_goveq%aux_vars_in(icell)%gbv  = 2.496430918408511d0
             cur_goveq%aux_vars_in(icell)%cp   = 744.5333333333334d0
             cur_goveq%aux_vars_in(icell)%fwet = 0.d0
             cur_goveq%aux_vars_in(icell)%fdry = 0.8218390792391702d0
             cur_goveq%aux_vars_in(icell)%dpai = dpai(k)/ntree
             if (cur_goveq%rank_in_soe_list == CLEF_TEMP_SUN_GE) then
                cur_goveq%aux_vars_in(icell)%fssh = fssh(k)
             else
                cur_goveq%aux_vars_in(icell)%fssh = 1.d0 - fssh(k)
             end if
          end do
          end do
       end do

    end select

    deallocate(dpai)
    deallocate(fssh)

  end subroutine set_common_canopy_leaf_parameters

  !------------------------------------------------------------------------
  subroutine get_dpai_fssh(dpai, fssh)

    implicit none
    PetscReal, pointer :: dpai(:), fssh(:)


    PetscInt                             :: k, i, num_int
    PetscReal                            :: Kb, sumpai
    PetscReal                            :: dz_leaf, qbeta, pbeta, pai
    PetscReal                            :: zl, zu, z_int, dz_int, zrel, beta_pdf, pad
    PetscReal                            :: pai_sum, pai_miss, pai_old, pai_new

    dpai(:) = 0.d0
    fssh(:) = 0.d0

    ! Determine plant area index increment for each layer by numerically
    ! integrating the plant area density (beta distribution) between
    ! the bottom and top heights for that layer

    dz_leaf = z_cleaf/nz_cleaf;
    pbeta   = 3.5d0;              ! Parameter for beta distribution
    qbeta   = 2.d0;               ! Parameter for beta distribution
    pai     = 5.051612734794617d0

    do k = 2, 42 + 1
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
          beta_pdf = (zrel**(pbeta-1) * (1.d0 - zrel)**(qbeta-1))/exp(log_gamma(pbeta) + log_gamma(qbeta) - log_gamma(pbeta+qbeta));
          pad = (pai / hc) * beta_pdf

          dpai(k) = dpai(k) + pad * dz_int

       end do
    end do

    pai_sum = 0;
    do k = 2,42+1
       pai_sum = pai_sum + dpai(k)
    end do

    ! Set layers with small plant area index to zero
    pai_miss = 0;
    do k = 2, 42+1
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
       do k = 2, 42+1
          dpai(k) = dpai(k) + pai_miss * (dpai(k) / pai_new);
       end do
    end if

    Kb = 1.762817445019839d0;
    do k = 43,7,-1
       if (k == 43) then
          sumpai = 0.5d0 * dpai(k);
       else
          sumpai = sumpai + &
               0.5d0 * (dpai(k+1) + dpai(k));
       end if
       fssh(k) = exp(-Kb * sumpai);
    end do

  end subroutine get_dpai_fssh

  !------------------------------------------------------------------------
  subroutine set_sunlit_canopy_parameters(mlc_mpp)
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
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt :: k, i, num_int, icair, icell, itree
    PetscReal :: Kb, sumpai
    PetscReal :: dz_leaf, qbeta, pbeta, pai
    PetscReal :: zl, zu, z_int, dz_int, zrel, beta_pdf, pad
    PetscReal :: pai_sum, pai_miss, pai_old, pai_new
    PetscReal, pointer :: gs(:)

    base_soe => mlc_mpp%soe

    call base_soe%SetPointerToIthGovEqn(CLEF_TEMP_SUN_GE, cur_goveq)

    allocate(gs(93))
    select type(cur_goveq)
    class is (goveqn_cleaf_temp_type)

       call get_sunlit_canopy_gs(gs)
       do icair = 1, ncair
          do itree = 1, ntree
          do k = 7,43
             icell = (icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k
             cur_goveq%aux_vars_in(icell)%gs = gs(k)
          end do

          k =  7;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.9869857739781d0
          k =  8;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.8100113537029d0
          k =  9;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.7147998761629d0
          k = 10;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.6645467566822d0
          k = 11;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.6422035725484d0
          k = 12;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.6392966303582d0
          k = 13;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.6514847604817d0
          k = 14;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.6766021357984d0
          k = 15;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.7137254254163d0
          k = 16;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.7627019640728d0
          k = 17;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.8238999626867d0
          k = 18;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.8980702313243d0
          k = 19;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.9862631887909d0
          k = 20;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 140.0897684653183d0
          k = 21;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 140.2100538053315d0
          k = 22;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 140.3486818847138d0
          k = 23;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 140.5071806149416d0
          k = 24;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 140.6868352048059d0
          k = 25;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 140.8883584829672d0
          k = 26;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 141.1113792315132d0
          k = 27;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 141.3536664423189d0
          k = 28;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 141.6099822011559d0
          k = 29;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 141.8704336551236d0
          k = 30;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 142.1181913093900d0
          k = 31;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 142.3264909566734d0
          k = 32;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 142.4550034158019d0
          k = 33;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 142.4460421185886d0
          k = 34;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 142.2218178601452d0
          k = 35;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 141.6851596824207d0
          k = 36;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 140.7277716843982d0
          k = 37;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 139.2518034108234d0
          k = 38;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 137.2114197261891d0
          k = 39;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 134.6805463548995d0
          k = 40;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 131.9550915266485d0
          k = 41;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 129.7361873094630d0
          k = 42;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 129.7993862020948d0
          k = 43;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = 143.7045065806239d0

       end do
       end do

    end select

    deallocate(gs)

  end subroutine set_sunlit_canopy_parameters

  !------------------------------------------------------------------------
  subroutine get_sunlit_canopy_gs(gs)
    !
    implicit none
    !
    PetscReal, pointer :: gs(:)
    !
    PetscInt :: k

    gs(:) = 0

    k =  7;  gs(k) = 0.1056193510550169d0
    k =  8;  gs(k) = 0.1058669704208841d0
    k =  9;  gs(k) = 0.1062166035088956d0
    k = 10;  gs(k) = 0.1066846074875817d0
    k = 11;  gs(k) = 0.1072854387286280d0
    k = 12;  gs(k) = 0.1080315168674592d0
    k = 13;  gs(k) = 0.1089335362366439d0
    k = 14;  gs(k) = 0.1100012607812562d0
    k = 15;  gs(k) = 0.1112447128077408d0
    k = 16;  gs(k) = 0.1126755044648808d0
    k = 17;  gs(k) = 0.1138467165585616d0
    k = 18;  gs(k) = 0.1170524695200598d0
    k = 19;  gs(k) = 0.1186451281076514d0
    k = 20;  gs(k) = 0.1206859738130298d0
    k = 21;  gs(k) = 0.1228219389652392d0
    k = 22;  gs(k) = 0.1263235652964973d0
    k = 23;  gs(k) = 0.1300019677357508d0
    k = 24;  gs(k) = 0.1322680545506565d0
    k = 25;  gs(k) = 0.1367071935229807d0
    k = 26;  gs(k) = 0.1408216759258680d0
    k = 27;  gs(k) = 0.1452273039039047d0
    k = 28;  gs(k) = 0.1499262843535941d0
    k = 29;  gs(k) = 0.1549264640058029d0
    k = 30;  gs(k) = 0.1611234013632947d0
    k = 31;  gs(k) = 0.1668845999057947d0
    k = 32;  gs(k) = 0.1727971327085968d0
    k = 33;  gs(k) = 0.1788628079180081d0
    k = 34;  gs(k) = 0.1850771375553107d0
    k = 35;  gs(k) = 0.1934140277837149d0
    k = 36;  gs(k) = 0.1998116684650200d0
    k = 37;  gs(k) = 0.2061626747018590d0
    k = 38;  gs(k) = 0.2124795008223110d0
    k = 39;  gs(k) = 0.2173241738995193d0
    k = 40;  gs(k) = 0.2228796106202699d0
    k = 41;  gs(k) = 0.2272584280787935d0
    k = 42;  gs(k) = 0.2303662043528620d0
    k = 43;  gs(k) = 0.2315636153119537d0

  end subroutine

  !------------------------------------------------------------------------
  subroutine set_shaded_canopy_parameters(mlc_mpp)
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
    type(mpp_mlc_type) :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    class(goveqn_base_type)    , pointer :: cur_goveq
    PetscInt :: k, icair, icell, itree
    PetscReal :: Kb, sumpai
    PetscReal, pointer :: gs(:)

    base_soe => mlc_mpp%soe

    allocate(gs(93))

    call base_soe%SetPointerToIthGovEqn(CLEF_TEMP_SHD_GE, cur_goveq)

    select type(cur_goveq)
    class is (goveqn_cleaf_temp_type)

       call get_shaded_canopy_gs(gs)

       do icair = 1, ncair
          do itree = 1, ntree
          do k = 7,43
             icell = (icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k
             cur_goveq%aux_vars_in(icell)%gs = gs(k)
          end do

          k =  7;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.411488333307743d0
          k =  8;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.234513913032590d0
          k =  9;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.139302435492522d0
          k = 10;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.089049316011852d0
          k = 11;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.066706131878055d0
          k = 12;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.063799189687813d0
          k = 13;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.075987319811279d0
          k = 14;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.101104695127997d0
          k = 15;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.138227984745972d0
          k = 16;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.187204523402388d0
          k = 17;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.248402522016342d0
          k = 18;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.322572790653995d0
          k = 19;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.410765748120540d0
          k = 20;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.514271024647946d0
          k = 21;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.634556364661143d0
          k = 22;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.773184444043382d0
          k = 23;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  1.931683174271291d0
          k = 24;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  2.111337764135555d0
          k = 25;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  2.312861042296822d0
          k = 26;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  2.535881790842842d0
          k = 27;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  2.778169001648555d0
          k = 28;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  3.034484760485563d0
          k = 29;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  3.294936214453243d0
          k = 30;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  3.542693868719610d0
          k = 31;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  3.750993516003037d0
          k = 32;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  3.879505975131512d0
          k = 33;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  3.870544677918208d0
          k = 34;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  3.646320419474797d0
          k = 35;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  3.109662241750371d0
          k = 36;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  2.152274243727867d0
          k = 37;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  0.676305970153017d0
          k = 38;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = -1.364077714481233d0
          k = 39;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = -3.894951085770870d0
          k = 40;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = -6.620405914021802d0
          k = 41;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = -8.839310131207370d0
          k = 42;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn = -8.776111238575538d0
          k = 43;  cur_goveq%aux_vars_in((icair-1)*(nz_cair+1)*ntree + (itree-1)*(nz_cleaf+1) + k)%rn =  5.129009139953610d0
       end do
       end do

    end select

  end subroutine set_shaded_canopy_parameters

  !------------------------------------------------------------------------
  subroutine get_shaded_canopy_gs(gs)
    !
    implicit none
    !
    PetscReal, pointer :: gs(:)
    !
    PetscInt :: k

    gs(:) = 0

    k =  7;  gs(k) = 2.0000000000000000d-003
    k =  8;  gs(k) = 2.0000000000000000d-003
    k =  9;  gs(k) = 2.0000000000000000d-003
    k = 10;  gs(k) = 2.0000000000000000d-003
    k = 11;  gs(k) = 2.0000000000000000d-003
    k = 12;  gs(k) = 2.0000000000000000d-003
    k = 13;  gs(k) = 2.0000000000000000d-003
    k = 14;  gs(k) = 2.0000000000000000d-003
    k = 15;  gs(k) = 2.0000000000000000d-003
    k = 16;  gs(k) = 2.0000000000000000d-003
    k = 17;  gs(k) = 2.0000000000000000d-003
    k = 18;  gs(k) = 2.0000000000000000d-003
    k = 19;  gs(k) = 2.0000000000000000d-003
    k = 20;  gs(k) = 2.0000000000000000d-003
    k = 21;  gs(k) = 2.0000000000000000d-003
    k = 22;  gs(k) = 2.0000000000000000d-003
    k = 23;  gs(k) = 2.0000000000000000d-003
    k = 24;  gs(k) = 5.2146013029975334d-003
    k = 25;  gs(k) = 5.5227688387169205d-003
    k = 26;  gs(k) = 9.2945439124555301d-003
    k = 27;  gs(k) = 9.4101275089266447d-003
    k = 28;  gs(k) = 1.2582674218550544d-002
    k = 29;  gs(k) = 1.6999874421743270d-002
    k = 30;  gs(k) = 2.3036435105984941d-002
    k = 31;  gs(k) = 2.7903866816023401d-002
    k = 32;  gs(k) = 3.7385308959493969d-002
    k = 33;  gs(k) = 4.6808450662473224d-002
    k = 34;  gs(k) = 5.9036977283335762d-002
    k = 35;  gs(k) = 7.1890808689035093d-002
    k = 36;  gs(k) = 8.7547774703355410d-002
    k = 37;  gs(k) = 0.1059444058487105d0
    k = 38;  gs(k) = 0.1228398700721039d0
    k = 39;  gs(k) = 0.1416660859387607d0
    k = 40;  gs(k) = 0.1584170776550386d0
    k = 41;  gs(k) = 0.1712280540285039d0
    k = 42;  gs(k) = 0.1801048624092090d0
    k = 43;  gs(k) = 0.1844507421254655d0

  end subroutine

  !------------------------------------------------------------------------
  subroutine set_turbulence_parameters(mlc_mpp)
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
    type(mpp_mlc_type)                   :: mlc_mpp
    !
    class(sysofeqns_base_type) , pointer :: base_soe
    class(sysofeqns_mlc_type)  , pointer :: soe
    PetscInt                             :: p

    base_soe => mlc_mpp%soe

    select type(base_soe)
    class is (sysofeqns_mlc_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do p = 1, ncair

       soe%cturb%pai(p)   = 5.051612734794617d0

       soe%cturb%hc(p)    = 21.d0
       soe%cturb%zref(p)  = 46.d0

       soe%cturb%pref(p)  = 98620.d0
       soe%cturb%uref(p)  = 5.169d0
       soe%cturb%tref(p)  = 295.9349938964844d0
       soe%cturb%rhref(p) = 53.871d0
       soe%cturb%tcan(p)  = soe%cturb%tref(p)

       call soe%cturb%ComputeDerivedAtmInputs(p)

       soe%cturb%qcan(p) = soe%cturb%qref(p)
    end do

  end subroutine set_turbulence_parameters

  !------------------------------------------------------------------------
  subroutine set_soil_parameters(mlc_mpp)
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
    type(mpp_mlc_type)                   :: mlc_mpp
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

       soe%cturb%soil_tk(icair)  = 1.261326601469150d0
       soe%cturb%soil_dz(icair)  = 7.1006354171935350d-003
       soe%cturb%soil_temperature(icair)   = 294.8492736816406d0
       soe%cturb%soil_res(icair) = 3361.509423807650d0
       soe%cturb%soil_rhg(icair) = 0.9984057411945876d0
       soe%cturb%soil_rn(icair)  = 1.896127799819662d0
    end do

  end subroutine set_soil_parameters


end module mlc_parameters
