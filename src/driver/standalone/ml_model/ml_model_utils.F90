module ml_model_utils
  !
  use mpp_varctl           , only : iulog
  use mpp_abortutils       , only : endrun
  use mpp_shr_log_mod      , only : errMsg => shr_log_errMsg
  !
  implicit none
  !
#include <petsc/finclude/petsc.h>
  PetscInt, parameter :: CANOPY_AIR_MESH      = 1
  PetscInt, parameter :: CANOPY_MESH          = 2
  PetscInt, parameter :: CANOPY_AND_SOIL_MESH = 3

  public :: compute_dpai_fssh
  public :: extract_data_from_mlc
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

    nbot = 1 + 1     ! 1st layer is soil layer
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
  subroutine extract_data_from_mlc(mlc_mpp)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use SystemOfEquationsBaseType       , only : sysofeqns_base_type
    use SystemOfEquationsMLCType        , only : sysofeqns_mlc_type
    use MultiPhysicsProbMLC             , only : mpp_mlc_type
    use ml_model_global_vars            , only : nbot, ntop, ncair, ntree, nz_cair
    use ml_model_global_vars            , only : Tleaf_sun, Tleaf_shd, Tsoil, Tair, eair
    use ml_model_global_vars            , only : CLEF_TEMP_SUN_GE, CLEF_TEMP_SHD_GE, CAIR_TEMP_GE, CAIR_VAPR_GE
    use GoverningEquationBaseType       , only : goveqn_base_type
    use GoveqnCanopyAirTemperatureType  , only : goveqn_cair_temp_type
    use GoveqnCanopyAirVaporType        , only : goveqn_cair_vapor_type
    use GoveqnCanopyLeafTemperatureType , only : goveqn_cleaf_temp_type
    use MultiPhysicsProbConstants       , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants       , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants       , only : VAR_LEAF_TEMPERATURE
    use petscvec
    !
    ! !ARGUMENTS
    implicit none
    !
    class(mpp_mlc_type)  :: mlc_mpp
    !
    PetscScalar, pointer :: soln_p(:)
    PetscInt             :: idx_leaf, idx_data, idx_soil, idx_air
    PetscInt             :: ileaf, icair, itree, k, ieqn
    PetscInt             :: ncells
    PetscReal,  pointer  :: tleaf_data(:), tair_data(:), eair_data(:)
    PetscErrorCode       :: ierr

    ncells = ncair*ntree*(nz_cair+1)
    allocate(tleaf_data(ncells))
    allocate(tair_data(ncells))
    allocate(eair_data(ncells))

    do ileaf = 1, 2
       ncells = ncair*ntree*(ntop-nbot+1)
       if (ileaf == 1) then
          ieqn = CLEF_TEMP_SUN_GE
       else
          ieqn = CLEF_TEMP_SHD_GE
       end if

       call get_data_from_mlc_eqn(mlc_mpp, ieqn, ncells, tleaf_data)

       idx_leaf = 0
       idx_data = 0
       do icair = 1, ncair
          do itree = 1, ntree
             do k = 1, nz_cair+1
                if (k>=nbot .and. k<=ntop) then
                   idx_leaf = idx_leaf + 1
                   idx_data = idx_data + 1
                   if (ileaf == 1) then
                      call set_value_in_condition(Tleaf_sun, idx_leaf, tleaf_data(idx_data))
                   else
                      call set_value_in_condition(Tleaf_shd, idx_leaf, tleaf_data(idx_data))
                   endif
                end if
             end do
          end do
       end do
    end do

    call get_data_from_mlc_eqn(mlc_mpp, CAIR_TEMP_GE, ncells, tair_data)
    call get_data_from_mlc_eqn(mlc_mpp, CAIR_VAPR_GE, ncells, eair_data)

    idx_soil = 0
    idx_data = 0
    idx_air = 0
    do icair = 1, ncair
       do k = 1, nz_cair+1
          idx_data = idx_data + 1
          if (k == 1) then
             idx_soil = idx_soil + 1
             call set_value_in_condition(Tsoil, idx_soil, tair_data(idx_data))
          else
             idx_air = idx_air + 1
             call set_value_in_condition(Tair, idx_air, tair_data(idx_data))
             call set_value_in_condition(eair, idx_air, eair_data(idx_data))
          end if
       end do
    end do

  end subroutine extract_data_from_mlc


  !------------------------------------------------------------------------
  subroutine get_data_from_mlc_eqn(mlc_mpp, ieqn, ncells, data)
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
   use MultiPhysicsProbConstants       , only : VAR_TEMPERATURE
   use MultiPhysicsProbConstants       , only : VAR_LEAF_TEMPERATURE
   use MultiPhysicsProbConstants       , only : VAR_WATER_VAPOR
    !
    ! !ARGUMENTS
   implicit none
   !
   class(mpp_mlc_type)  :: mlc_mpp
   PetscInt             :: ieqn
   PetscInt             :: ncells
   PetscReal,  pointer  :: data(:)
   !
   class(goveqn_base_type) , pointer :: goveq
   PetscErrorCode       :: ierr

   call mlc_mpp%soe%SetPointerToIthGovEqn(ieqn, goveq)

   select type(goveq)
   class is (goveqn_cair_temp_type)
      call goveq%GetRValues(AUXVAR_INTERNAL, VAR_TEMPERATURE, ncells, data)
   class is (goveqn_cair_vapor_type)
      call goveq%GetRValues(AUXVAR_INTERNAL, VAR_WATER_VAPOR, ncells, data)
   class is (goveqn_cleaf_temp_type)
      call goveq%GetRValues(AUXVAR_INTERNAL, VAR_LEAF_TEMPERATURE, ncells, data)
   end select

   end subroutine get_data_from_mlc_eqn

end module ml_model_utils
