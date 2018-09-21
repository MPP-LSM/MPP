module thermal_mms_steady_state_problem_3D

   use thermal_mms_vars
   use mpp_varctl              , only : iulog
   use mpp_abortutils          , only : endrun
   use mpp_shr_log_mod         , only : errMsg => shr_log_errMsg
   
   implicit none

#include <petsc/finclude/petsc.h>

   public :: set_variable_for_steady_state_3D

   private :: compute_temperature_or_deriv
   private :: compute_thermal_conductivity_or_deriv

 contains

   !------------------------------------------------------------------------
   subroutine compute_temperature_or_deriv(x, y, z, T, dT_dx, dT_dy, dT_dz, &
        d2T_dx2, d2T_dy2, d2T_dz2)
     !
     ! !DESCRIPTION:
     !
     implicit none
     !
     PetscReal, intent(in)            :: x, y, z
     PetscReal, intent(out), optional :: T
     PetscReal, intent(out), optional :: dT_dx, dT_dy, dT_dz
     PetscReal, intent(out), optional :: d2T_dx2, d2T_dy2, d2T_dz2

     if (present(T       )) T       =  10.d0      *sin(x*PI)*cos(2.d0*PI*y)*sin(3.d0*PI*z) + 270.d0
     if (present(dT_dx   )) dT_dx   =  10.d0*PI   *cos(x*PI)*cos(2.d0*PI*y)*sin(3.d0*PI*z)
     if (present(dT_dy   )) dT_dy   = -20.d0*PI   *sin(x*PI)*sin(2.d0*PI*y)*sin(3.d0*PI*z)
     if (present(dT_dz   )) dT_dz   =  30.d0*PI   *sin(x*PI)*cos(2.d0*PI*y)*cos(3.d0*PI*z)
     if (present(d2T_dx2 )) d2T_dx2 = -10.d0*PI*PI*sin(x*PI)*cos(2.d0*PI*y)*sin(3.d0*PI*z)
     if (present(d2T_dy2 )) d2T_dy2 = -40.d0*PI*PI*sin(x*PI)*cos(2.d0*PI*y)*sin(3.d0*PI*z)
     if (present(d2T_dz2 )) d2T_dz2 = -90.d0*PI*PI*sin(x*PI)*cos(2.d0*PI*y)*sin(3.d0*PI*z)
     
   end subroutine compute_temperature_or_deriv
   
   !------------------------------------------------------------------------
   subroutine compute_thermal_conductivity_or_deriv(x, y, z, l, &
        dl_dx, dl_dy, dl_dz)
     !
     ! !DESCRIPTION:
     !
     implicit none
     !
     PetscReal, intent(in)            :: x, y, z
     PetscReal, intent(out), optional :: l
     PetscReal, intent(out), optional :: dl_dx, dl_dy, dl_dz
     !
     PetscReal :: alpha

     alpha   = exp(x + y + z - 1.d0)
    
     !if (present(l     )) l       =  (y + 0.5d0)         * ((z + 0.5d0)**2.d0                   ) * alpha
     !if (present(dl_dx )) dl_dx   = ((y + 0.5d0)       ) * ((z + 0.5d0)**2.d0                   ) * alpha
     !if (present(dl_dy )) dl_dy   = ((y + 0.5d0) + 1.d0) * ((z + 0.5d0)**2.d0                   ) * alpha
     !if (present(dl_dz )) dl_dz   = ((y + 0.5d0)       ) * ((z + 0.5d0)**2.d0 + 2.d0*(z + 0.5d0)) * alpha
     
     !if (present(l     )) l       = (x*y + 0.5d0)*alpha
     !if (present(dl_dx )) dl_dx   = (x*y + y    )*alpha
     !if (present(dl_dy )) dl_dy   = (x*y + x    )*alpha
     !if (present(dl_dz )) dl_dz   = (x*y        )*alpha
     
     if (present(l     )) l       = alpha
     if (present(dl_dx )) dl_dx   = alpha
     if (present(dl_dy )) dl_dy   = alpha
     if (present(dl_dz )) dl_dz   = alpha
     
   end subroutine compute_thermal_conductivity_or_deriv
   
  !------------------------------------------------------------------------
  subroutine set_variable_for_steady_state_3D(data_type, data_1D, data_2D)
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
    PetscReal                    :: alpha
    PetscReal                    :: x, y, z
    PetscReal                    :: l, dl_dx, dl_dy, dl_dz
    PetscReal                    :: dT_dx, dT_dy, dT_dz
    PetscReal                    :: d2T_dx2, d2T_dy2, d2T_dz2
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
                count   = count + 1
                x       = soil_xc_3d(ii,jj,kk)
                y       = soil_yc_3d(ii,jj,kk)
                z       = soil_zc_3d(ii,jj,kk)

                call compute_thermal_conductivity_or_deriv(x, y, z, l=data_2D(count,kk))
             end do
          end do
       end do
       
    case(DATA_THERMAL_CONDUCTIVITY_1D_VEC)
       count = 0
       do kk = 1, nz
          do jj = 1, ny
             do ii = 1, nx
                count   = count + 1
                x       = soil_xc_3d(ii,jj,kk)
                y       = soil_yc_3d(ii,jj,kk)
                z       = soil_zc_3d(ii,jj,kk)

                call compute_thermal_conductivity_or_deriv(x, y, z, l=data_1D(count))
             end do
          end do
       end do
       
    case (DATA_TEMPERATURE)
       count = 0

       do kk = 1, nz
          do jj = 1, ny
             do ii = 1, nx
                count = count + 1
                x     = soil_xc_3d(ii,jj,kk)
                y     = soil_yc_3d(ii,jj,kk)
                z     = soil_zc_3d(ii,jj,kk)
             
                call compute_temperature_or_deriv(x, y, z, T=data_1D(count))
             end do
          end do
       end do

    case(DATA_TEMPERATURE_BC)

       count = 0;

       do kk = 1,nz
          do jj = 1,ny
             ! At the beginning
             ii    = 1;
             count = count + 1
             x     = soil_xc_3d(ii,jj,kk)-dx/2.d0 
             y     = soil_yc_3d(ii,jj,kk)
             z     = soil_zc_3d(ii,jj,kk)

             call compute_temperature_or_deriv(x, y, z, T=data_1D(count))

             ! At the end
             ii    = nx
             count = count + 1
             x     = soil_xc_3d(ii,jj,kk)+dx/2.d0
             y     = soil_yc_3d(ii,jj,kk)
             z     = soil_zc_3d(ii,jj,kk)

             call compute_temperature_or_deriv(x, y, z, T=data_1D(count))
          end do
       end do

       do kk = 1,nz
          do ii = 1,nx
             ! At the beginning
             jj    = 1;
             count = count + 1
             x     = soil_xc_3d(ii,jj,kk)
             y     = soil_yc_3d(ii,jj,kk)-dy/2.d0
             z     = soil_zc_3d(ii,jj,kk)

             call compute_temperature_or_deriv(x, y, z, T=data_1D(count))

             ! At the end
             jj    = ny
             count = count + 1
             x     = soil_xc_3d(ii,jj,kk)
             y     = soil_yc_3d(ii,jj,kk)+dy/2.d0
             z     = soil_zc_3d(ii,jj,kk)

             call compute_temperature_or_deriv(x, y, z, T=data_1D(count))
          end do
       end do

       do jj = 1,ny
          do ii = 1,nx
             ! At the beginning
             kk    = 1;
             count = count + 1
             x     = soil_xc_3d(ii,jj,kk)
             y     = soil_yc_3d(ii,jj,kk)
             z     = soil_zc_3d(ii,jj,kk)-dz/2.d0

             call compute_temperature_or_deriv(x, y, z, T=data_1D(count))

             ! At the end
             kk    = nz
             count = count + 1
             x     = soil_xc_3d(ii,jj,kk)
             y     = soil_yc_3d(ii,jj,kk)
             z     = soil_zc_3d(ii,jj,kk)+dz/2.d0

             call compute_temperature_or_deriv(x, y, z, T=data_1D(count))
          end do
       end do

       case (DATA_HEAT_SOURCE)
          count = 0
          do kk = 1, nz
             do jj = 1, ny
                do ii = 1, nx
                   count   = count + 1
                   x       = soil_xc_3d(ii,jj,kk)
                   y       = soil_yc_3d(ii,jj,kk)
                   z       = soil_zc_3d(ii,jj,kk)

                   call compute_thermal_conductivity_or_deriv(x, y, z, l, &
                        dl_dx, dl_dy, dl_dz)
                   
                   call compute_temperature_or_deriv(x, y, z, &
                        dT_dx=dT_dx, dT_dy=dT_dy, dT_dz=dT_dz, &
                        d2T_dx2=d2T_dx2,d2T_dy2=d2T_dy2,d2T_dz2=d2T_dz2)


                   data_1D(count) = (&
                        -dl_dx*dT_dx - l * d2T_dx2 &
                        -dl_dy*dT_dy - l * d2T_dy2 &
                        -dl_dz*dT_dz - l * d2T_dz2 &
                        )*dx*dy*dz
                end do
             end do
          end do

       case default
          write(iulog,*)'Unsupported problem in material'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       end select

     end subroutine set_variable_for_steady_state_3D

end module thermal_mms_steady_state_problem_3D
