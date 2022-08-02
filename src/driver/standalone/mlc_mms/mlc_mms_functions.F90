module mlc_mms_functions

#include <petsc/finclude/petsc.h>
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbMLC , only : mpp_mlc_type
  use mlc_mms_global_vars
  use petscsys
  use mlc_mms_global_vars, only : zmax_l, zmin_l, z_l

  implicit none

  public :: mms_gbv
  public :: mms_gbh
  public :: mms_gs_sun
  public :: mms_gs_shd
  public :: mms_tair
  public :: mms_qair
  public :: mms_tl_sun
  public :: mms_tl_shd
  public :: mms_rn_sun
  public :: mms_rn_shd
contains
  
  !------------------------------------------------------------------------
  subroutine mms_conductance(z,coeff,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactor solution for conductances
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z, coeff
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal :: arg, a, b
    PetscErrorCode :: ierr

    if (z > zmax_l .or. z < zmin_l) then
       f     = 0.d0
       df_dz = 0.d0
       df_dt = 0.d0
    else       
       arg = (z - zmin_l) * PETSC_PI/z_l
       a = sin(arg)
       b = cos(arg)

       f     = coeff*(a**2.d0)
       df_dz = coeff*(2.d0*a )*PETSC_PI/z_l*b
       df_dt = 0.d0
    end if
    
  end subroutine mms_conductance

  !------------------------------------------------------------------------
  subroutine mms_gbv(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for boundary layer conductance for water vapor
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal, parameter :: coeff = 2.0d0

    call mms_conductance(z,coeff,f,df_dz,df_dt)

  end subroutine mms_gbv

  !------------------------------------------------------------------------
  subroutine mms_gbh(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for boundary layer conductance for heat
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal, parameter :: coeff = 2.5d0

    call mms_conductance(z,coeff,f,df_dz,df_dt)

  end subroutine mms_gbh

  !------------------------------------------------------------------------
  subroutine mms_gs_sun(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for stomatal conductance for sunlit leaf
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal, parameter :: coeff = 0.3d0

    call mms_conductance(z,coeff,f,df_dz,df_dt)

  end subroutine mms_gs_sun

  !------------------------------------------------------------------------
  subroutine mms_gs_shd(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for stomatal conductance for sunlit leaf
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal, parameter :: coeff = 0.1d0

    call mms_conductance(z,coeff,f,df_dz,df_dt)

  end subroutine mms_gs_shd

  !------------------------------------------------------------------------
  subroutine mms_tair(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for air temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal :: arg, a, b, da_dz
    PetscErrorCode :: ierr

    arg = z*PETSC_PI*2.d0
    a = sin(arg)
    b = cos(arg)

    da_dz = (PETSC_PI*2.d0)*b

    f     = 300.d0 + 10.d0 * a
    df_dz =          10.d0 * da_dz
    df_dt = 0.d0

  end subroutine mms_tair

  !------------------------------------------------------------------------
  subroutine mms_qair(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for water vapor
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal :: arg, a, b, da_dz
    PetscErrorCode :: ierr

    arg = z*PETSC_PI*2.d0
    a = cos(arg)
    b = sin(arg)

    da_dz = (PETSC_PI*2.d0)*(-b)
    
    f     = 0.02d0 + 0.01d0 * a
    df_dz =          0.01d0 * da_dz
    df_dt = 0.d0

  end subroutine mms_qair

  !------------------------------------------------------------------------
  subroutine mms_tl_sun(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for sunlit leaf temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal :: arg, a, b, da_dz
    PetscErrorCode :: ierr

    if (z > zmax_l .or. z < zmin_l) then
       f     = 0.d0
       df_dz = 0.d0
       df_dt = 0.d0
    else       
       arg = (z - zmin_l) * PETSC_PI/z_l
       a = cos(arg)
       b = sin(arg)
       da_dz = (PETSC_PI/z_l)*(-b)
       
       f     = 280.d0 + 5.d0 * a
       df_dz =          5.d0 * da_dz
       df_dt = 0.d0
    end if

  end subroutine mms_tl_sun

  !------------------------------------------------------------------------
  subroutine mms_tl_shd(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for shaded leaf temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt
    !
    PetscReal :: arg, a, b, da_dz
    PetscErrorCode :: ierr

    if (z > zmax_l .or. z < zmin_l) then
       f     = 0.d0
       df_dz = 0.d0
       df_dt = 0.d0
    else       
       arg = 2.d0*(z - zmin_l) * PETSC_PI/z_l
       a = cos(arg)
       b = sin(arg)
       da_dz = (2.d0*PETSC_PI/z_l)*(-b)
       
       f     = 275.d0 - 5.d0 * a
       df_dz =        - 5.d0 * da_dz
       df_dt = 0.d0
    end if

  end subroutine mms_tl_shd

  !------------------------------------------------------------------------
  subroutine mms_rn_sun(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for shaded leaf temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt

    if (z > zmax_l .or. z < zmin_l) then
       f     = 0.d0
       df_dz = 0.d0
       df_dt = 0.d0
    else       
       f     = 200.d0 - 100.d0 * exp(-z)
       df_dz =        + 100.d0 * exp(-z)
       df_dt = 0.d0
    end if

  end subroutine mms_rn_sun

  !------------------------------------------------------------------------
  subroutine mms_rn_shd(z,f,df_dz,df_dt)
    !
    ! !DESCRIPTION:
    ! Manufactured solution for shaded leaf temperature
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in) :: z
    PetscReal, intent(out) :: f, df_dz, df_dt

    if (z > zmax_l .or. z < zmin_l) then
       f     = 0.d0
       df_dz = 0.d0
       df_dt = 0.d0
    else       
       f     = 100.d0 - 70.d0 * exp(-z)
       df_dz =        + 70.d0 * exp(-z)
       df_dt = 0.d0
    end if

  end subroutine mms_rn_shd

end module mlc_mms_functions
