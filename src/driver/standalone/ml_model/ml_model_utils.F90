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

  public :: compute_dpai
  public :: compute_fssh_and_cumlai
  public :: set_value_in_condition
  public :: get_value_from_condition

contains

  !------------------------------------------------------------------------
  subroutine compute_dpai(dpai, cumlai, sumpai)
    !
    use ml_model_global_vars , only : hc, nveg, nbot, ntop, nz_cair, ncair, dz_cair
    !
    implicit none
    !
    PetscReal, pointer :: dpai(:), fssh(:), cumlai(:), sumpai(:)
    !
    PetscInt  :: k, i, num_int, ic_bot
    PetscReal :: Kb
    PetscReal :: dz_leaf, qbeta, pbeta, pai
    PetscReal :: zl, zu, z_int, dz_int, zrel, beta_pdf, pad
    PetscReal :: pai_sum, pai_miss, pai_old, pai_new

    dpai(:) = 0.d0
    cumlai(:) = 0.d0
    sumpai(:) = 0.d0

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

  end subroutine compute_dpai

  !------------------------------------------------------------------------
  subroutine compute_fssh_and_cumlai(nbot, ntop, dpai, fssh, cumlai, sumpai)
    !
    implicit none
    !
    PetscInt  , intent(in)          :: nbot, ntop
    PetscReal , pointer, intent(in) :: dpai(:)
    PetscReal , intent(inout)       :: fssh(:), cumlai(:), sumpai(:)
    !
    PetscInt :: k
    PetscReal :: Kb

    fssh(:) = 0.d0
    cumlai(:) = 0.d0
    sumpai(:) = 0.d0

    Kb = 1.7628174450198393d0;
    do k = ntop, nbot, -1
       if (k == ntop) then
          sumpai(k) = 0.5d0 * dpai(k);
          cumlai(k) = dpai(k);
       else
          sumpai(k) = sumpai(k+1) + &
               0.5d0 * (dpai(k+1) + dpai(k));
          cumlai(k) = cumlai(k+1) + dpai(k)
       end if
       fssh(k) = exp(-Kb * sumpai(k));
    end do

  end subroutine compute_fssh_and_cumlai

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
