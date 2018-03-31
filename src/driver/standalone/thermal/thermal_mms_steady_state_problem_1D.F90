module thermal_mms_steady_state_problem_1D

   use thermal_mms_vars
   use mpp_varctl              , only : iulog
   use mpp_abortutils          , only : endrun
   use mpp_shr_log_mod         , only : errMsg => shr_log_errMsg
   
   implicit none

#include <petsc/finclude/petsc.h>

   public :: set_variable_for_steady_state_1D

   private :: compute_temperature_or_deriv
   private :: compute_thermal_conductivity_or_deriv

 contains
   
   !------------------------------------------------------------------------
   subroutine compute_temperature_or_deriv(x,y,z,T,dT_dx,d2T_dx2)
     !
     ! !DESCRIPTION:
     !
     implicit none
     !
     PetscReal, intent(in)            :: x,y,z
     PetscReal, intent(out), optional :: T
     PetscReal, intent(out), optional :: dT_dx
     PetscReal, intent(out), optional :: d2T_dx2

     if (present(T       )) T       = 10 * sin (PI*x) + 270.d0
     if (present(dT_dx   )) dT_dx   = 10.d0*PI*cos(PI*x)
     if (present(d2T_dx2 )) d2T_dx2 = -10.d0*PI*PI*sin(PI*x)

   end subroutine compute_temperature_or_deriv
   
   !------------------------------------------------------------------------
   subroutine compute_thermal_conductivity_or_deriv(x,y,z,l,dl_dx)
     !
     ! !DESCRIPTION:
     !
     implicit none
     !
     PetscReal, intent(in)            :: x,y,z
     PetscReal, intent(out), optional :: l
     PetscReal, intent(out), optional :: dl_dx

     if (present(l     )) l     = exp(x)
     if (present(dl_dx )) dl_dx = exp(x)
     
   end subroutine compute_thermal_conductivity_or_deriv
   
   !------------------------------------------------------------------------
   subroutine set_variable_for_steady_state_1D(data_type, data_1D, data_2D)
     !
     ! !DESCRIPTION:
     !
     !
     implicit none
     !
     PetscInt                     :: data_type
     PetscReal, pointer, optional :: data_1D(:)
     PetscReal, pointer, optional :: data_2D(:,:)
     !
     PetscInt                     :: ii, jj, kk
     PetscInt                     :: count
     PetscReal                    :: x, y, z
     PetscReal                    :: l, dl_dx
     PetscReal                    :: dT_dx, d2T_dx2
     !-----------------------------------------------------------------------

     ! Soil properties

     select case(data_type)
     case (DATA_XC)
        count = 0
        do kk = 1,nz
           do jj = 1,ny
              do ii = 1, nx
                 count = count + 1
                 data_1D(count) = dx/2.d0 + dx * (ii - 1)
              end do
           end do
        end do

     case (DATA_YC)
        count = 0
        do kk = 1,nz
           do jj = 1,ny
              do ii = 1, nx
                 count = count + 1
                 data_1D(count) = dy/2.d0 + dy * (jj - 1)
              end do
           end do
        end do

     case (DATA_ZC)
        count = 0
        do kk = 1,nz
           do jj = 1,ny
              do ii = 1, nx
                 count = count + 1
                 data_1D(count) = dz/2.d0 + dz * (kk - 1)
              end do
           end do
        end do

     case(DATA_THERMAL_CONDUCTIVITY)
        do kk = 1, nz
           count = 0
           do jj = 1, ny
              do ii = 1, nx
                 count = count + 1
                 x     = soil_xc_3d(ii,jj,kk)
                 y     = soil_yc_3d(ii,jj,kk)
                 z     = soil_zc_3d(ii,jj,kk)
                 call compute_thermal_conductivity_or_deriv(x,y,z,l=data_2D(count,kk))
              end do
           end do
        end do

     case(DATA_TEMPERATURE)
        jj = 1; kk = 1
        count = 0

        do ii = 1, nx
           x     = soil_xc_3d(ii,jj,kk)
           y     = soil_yc_3d(ii,jj,kk)
           z     = soil_zc_3d(ii,jj,kk)

           call compute_temperature_or_deriv(x,y,z,T=data_1D(ii))
        end do

     case (DATA_TEMPERATURE_BC)
        jj = 1; kk = 1
        count = 0

        ii    = 1
        x     = soil_xc_3d(ii,jj,kk)-dx/2.d0
        y     = soil_yc_3d(ii,jj,kk)
        z     = soil_zc_3d(ii,jj,kk)
        count = count + 1;          
        call compute_temperature_or_deriv(x,y,z,T=data_1D(count))

        ii    = nx
        x     = soil_xc_3d(ii,jj,kk)+dx/2.d0
        y     = soil_yc_3d(ii,jj,kk)
        z     = soil_zc_3d(ii,jj,kk)
        count = count + 1;          
        call compute_temperature_or_deriv(x,y,z,T=data_1D(count))

     case(DATA_HEAT_SOURCE)
        jj = 1; kk = 1;
        do ii = 1, nx
           x = soil_xc_3d(ii,1,1)
           y = soil_yc_3d(ii,jj,kk)
           z = soil_zc_3d(ii,jj,kk)

           call compute_thermal_conductivity_or_deriv(x,y,z,l=l,dl_dx=dl_dx)
           call compute_temperature_or_deriv(x,y,z,dT_dx=dT_dx,d2T_dx2=d2T_dx2)

           data_1D(ii) = (-dl_dx * dT_dx - l * d2T_dx2)*dx*dy*dz
        end do

     case default
        write(iulog,*)'Unsupported problem in material'
        call endrun(msg=errMsg(__FILE__, __LINE__))
     end select

   end subroutine set_variable_for_steady_state_1D

end module thermal_mms_steady_state_problem_1D
