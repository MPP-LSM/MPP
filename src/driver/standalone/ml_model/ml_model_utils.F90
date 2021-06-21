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

  public :: compute_vertical_veg_structure
  public :: compute_fssh
  public :: set_value_in_condition
  public :: get_value_from_condition

contains

  !------------------------------------------------------------------------
  function cummulative_area_index(pbeta, qbeta, z_u, z_l, hc) result(cumm_area_index)
    !
    implicit none
    !
    PetscReal, intent(in) :: pbeta, qbeta, z_u, z_l, hc
    PetscReal :: cumm_area_index
    !
    PetscReal :: dz_int, z_int, beta_pdf, zrel
    PetscInt  :: i
    !
    PetscInt, parameter :: num_int = 100

    ! Numerical integration between zl and zu using 100 sublayers

    cumm_area_index = 0.d0

    dz_int = (z_u - z_l)/num_int

    do i = 1, num_int
       if (i == 1) then
          z_int = z_l + 0.5d0 * dz_int
       else
          z_int = z_int + dz_int
       end if
       zrel = min(z_int/hc, 1.d0)

       beta_pdf = (zrel**(pbeta-1) * (1.d0 - zrel)**(qbeta-1))/exp(lgamma(pbeta) + lgamma(qbeta) - lgamma(pbeta+qbeta));

       cumm_area_index = cumm_area_index + beta_pdf * dz_int
    end do

  end function cummulative_area_index

  !------------------------------------------------------------------------
  subroutine compute_vertical_veg_structure(dlai, dsai, dpai, cumpai, sumpai, leaf_td)
    !
    use ml_model_global_vars , only : hc, nveg, nbot, ntop, nz_cair, ncair, dz_cair, ntree
    !
    implicit none
    !
    PetscReal, pointer, intent(inout) :: dlai(:), dsai(:), dpai(:), cumpai(:), sumpai(:), leaf_td(:)
    !
    PetscInt  :: k, i, num_int, ic_bot
    PetscReal :: Kb
    PetscReal :: dz_leaf
    PetscReal :: zl, zu, z_int, dz_int, zrel, beta_pdf, pad
    PetscReal :: lai_sum, lai_miss, lai_old, lai_new
    PetscReal :: sai_sum, sai_miss, sai_old, sai_new
    !
    PetscReal, parameter :: lai_pbeta = 3.5d0
    PetscReal, parameter :: lai_qbeta = 2.0d0
    PetscReal, parameter :: sai_pbeta = 3.5d0
    PetscReal, parameter :: sai_qbeta = 2.0d0
    PetscReal, parameter :: lai      = 4.1516127586364746d0
    PetscReal, parameter :: sai      = 0.89999997615814209d0

    allocate(dlai  (nz_cair*ntree +1))
    allocate(dsai  (nz_cair*ntree +1))
    allocate(dpai  (nz_cair*ntree +1))
    allocate(cumpai(nz_cair*ntree +1))
    allocate(sumpai(nz_cair*ntree +1))
    allocate(leaf_td(nz_cair*ntree +1))

    dlai(:) = 0.d0
    dsai(:) = 0.d0
    dpai(:) = 0.d0
    cumpai(:) = 0.d0
    sumpai(:) = 0.d0
    leaf_td(:)= 0.d0

    ! Determine plant area index increment for each layer by numerically
    ! integrating the plant area density (beta distribution) between
    ! the bottom and top heights for that layer

    dz_leaf = dz_cair;

    nbot = 1 + 1     ! 1st layer is soil layer
    ntop = nveg + 1  ! 

    do k = nbot, ntop
       zl = dz_leaf * (k-2);
       zu = dz_leaf * (k-1);

       dlai(k) = cummulative_area_index(lai_pbeta, lai_qbeta, zu, zl, hc) * (lai/hc)
       dsai(k) = cummulative_area_index(sai_pbeta, sai_qbeta, zu, zl, hc) * (sai/hc)
    end do

    lai_sum = 0;
    sai_sum = 0;
    do k = nbot, ntop
       lai_sum = lai_sum + dlai(k)
       sai_sum = sai_sum + dsai(k)
    end do

    ! Set layers with small plant area index to zero
    lai_miss = 0.d0
    sai_miss = 0.d0
    do k = nbot, ntop
       if (dlai(k) + dsai(k) < 0.01d0) then

          lai_miss = lai_miss + dlai(k)
          sai_miss = sai_miss + dsai(k)

          dlai(k) = 0.d0;
          dsai(k) = 0.d0;
       end if
    end do

    ! Distribute the missing plant area across vegetation layers
    ! in proportion to the plant area profile

    if (lai_miss > 0.d0) then
       lai_old = lai_sum;
       lai_new = lai_old - lai_miss;
       do k = nbot, ntop
          dlai(k) = dlai(k) + lai_miss * (dlai(k) / lai_new);
       end do
    end if

    if (sai_miss > 0.d0) then
       sai_old = sai_sum;
       sai_new = sai_old - sai_miss;
       do k = nbot, ntop
          dsai(k) = dsai(k) + sai_miss * (dsai(k) / sai_new);
       end do
    end if

    ! Find the lowest vegetation layer
    ic_bot = 0
    do k = ntop, nbot,-1
       if (dlai(k) + dsai(k) > 0.d0) then
          ic_bot = k
       end if
    end do
    nbot = ic_bot
    if (nbot == 0) then
       call endrun (msg=' ERROR: compute_vertical_veg_structure: nbot not defined')
    endif

    do k = ntop, nbot, -1
       dpai(k) = dlai(k) + dsai(k)
       if (k == ntop) then
          sumpai(k) = 0.5d0 * dpai(k);
          cumpai(k) = dpai(k);
       else
          sumpai(k) = sumpai(k+1) + 0.5d0 * (dpai(k+1) + dpai(k));
          cumpai(k) = cumpai(k+1) + dpai(k)
       end if
    end do

  end subroutine compute_vertical_veg_structure

  !------------------------------------------------------------------------
  subroutine compute_fssh(nbot, ntop, sumpai, fssh)
    !
    use ml_model_global_vars, only : nz_cair, ntree
    !
    implicit none
    !
    PetscInt  , intent(in)          :: nbot, ntop
    PetscReal , pointer, intent(in) :: sumpai(:)
    PetscReal , pointer, intent(inout)       :: fssh(:)
    !
    PetscInt :: k
    PetscReal, parameter :: Kb = 1.7628168203899233d0;

    allocate(fssh  (nz_cair*ntree +1))

    fssh(:) = 0.d0

    do k = ntop, nbot, -1
       fssh(k) = exp(-Kb * sumpai(k));
    end do

  end subroutine compute_fssh

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

end module ml_model_utils
