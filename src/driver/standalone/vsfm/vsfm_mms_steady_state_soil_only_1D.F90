module vsfm_mms_steady_state_soil_only_1D

#include <petsc/finclude/petsc.h>

  use vsfm_mms_vars
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbConstants , only : PRESSURE_REF
  use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT
  use EOSWaterMod              , only : DENSITY_TGDPB01,DENSITY_CONSTANT
  use mpp_varcon                , only : denh2o
  use EOSWaterMod              , only : Density, Viscosity
  
  implicit none

  public :: set_default_prob_for_steady_soil_only_1D
  public :: setup_prob_for_steady_soil_only_1D
  public :: set_variable_for_steady_state_1D

contains

  !------------------------------------------------------------------------
  subroutine set_default_prob_for_steady_soil_only_1D()

    implicit none

    problem_type = STEADY_STATE_SOIL_ONLY_1D

    nx = 20
    ny = 1
    nz = 1

    x_min = 0.d0; x_max = 10.d0
    y_min = 0.d0; y_max = 1.d0
    z_min = 0.d0; z_max = 1.d0

    xlim = x_max - x_min
    ylim = y_max - y_min
    zlim = z_max - z_min

  end subroutine set_default_prob_for_steady_soil_only_1D
    
  !------------------------------------------------------------------------
   subroutine compute_pressure_or_deriv(x, val, dval_dx, d2val_dx2)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     PetscReal, intent(out), optional :: d2val_dx2
     !
     PetscReal, parameter             :: g   = GRAVITY_CONSTANT
     PetscReal, parameter             :: a0 =  15000.d0
     PetscReal, parameter             :: a1 = -20000.d0

     if (fully_saturated) then
        if (present(val      )) val       =  a0                *sin((x-x_min)/xlim*pi) - a1 + PRESSURE_REF
     else
        if (present(val      )) val       =  a0                *sin((x-x_min)/xlim*pi) + a1 + PRESSURE_REF
     end if
     if (present(dval_dx  )) dval_dx   =  a0*PI/xlim        *cos((x-x_min)/xlim*pi)
     if (present(d2val_dx2)) d2val_dx2 = -a0*PI*PI/xlim/xlim*sin((x-x_min)/xlim*pi)
     
   end subroutine compute_pressure_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_permeability_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: p0 = 1.d-11 ! [m^2]

     if (present(val    )) val     = p0        *(2.d0 + cos((x-x_min)/xlim * pi))
     if (present(dval_dx)) dval_dx = p0*PI/xlim*(      -sin((x-x_min)/xlim * pi))

   end subroutine compute_permeability_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_alpha_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 1.d0/4000.d0 ! [Pa^{-1}]

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_alpha_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_lambda_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 0.5d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_lambda_or_deriv

  !------------------------------------------------------------------------
   subroutine compute_residualsat_or_deriv(x, val, dval_dx)
     !
     implicit none
     !
     PetscReal, intent(in)            :: x
     PetscReal, intent(out), optional :: val
     PetscReal, intent(out), optional :: dval_dx
     !
     PetscReal, parameter :: a0 = 0.5d0

     if (present(val    )) val     =  a0
     if (present(dval_dx)) dval_dx =  0.d0
     
   end subroutine compute_residualsat_or_deriv

  !------------------------------------------------------------------------
  subroutine setup_prob_for_steady_soil_only_1D()
    !
    !
    use mpp_mesh_utils, only : ComputeXYZCentroids
    !
    implicit none

    if (nx <= 1) then
       write(iulog,*)'For STEADY_STATE_SOIL_ONLY_1D, nx should be greater than 1'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (ny /= 1) then
       write(iulog,*)'For STEADY_STATE_SOIL_ONLY_1D, ny should be set to 1'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    if (nz /= 1) then
       write(iulog,*)'For STEADY_STATE_SOIL_ONLY_1D, nz should be set to 1'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    dx = (x_max - x_min)/nx
    dy = (y_max - y_min)/ny
    dz = (z_max - z_min)/nz

    ncells = nx*ny*nz
    ncells_ghost = 0
    
    allocate(soil_xc_3d(nx,ny,nz))
    allocate(soil_yc_3d(nx,ny,nz))
    allocate(soil_zc_3d(nx,ny,nz))
    allocate(soil_id_3d(nx,ny,nz))
    
    call ComputeXYZCentroids(nx, ny, nz, dx, dy, dz, &
         x_min, y_min, z_min, &
         soil_xc_3d, soil_yc_3d, soil_zc_3d)

  end subroutine setup_prob_for_steady_soil_only_1D

  !------------------------------------------------------------------------
  subroutine set_variable_for_steady_state_1D(data_type, data_1D)
    !
    ! !DESCRIPTION:
    !
    use EOSWaterMod              , only : DENSITY_TGDPB01,DENSITY_CONSTANT
    use EOSWaterMod              , only : Density, Viscosity
    use MultiPhysicsProbConstants, only : FMWH2O
    use SaturationFunction       , only : saturation_params_type, SatFunc_Set_VG, SatFunc_PressToSat, SatFunc_PressToRelPerm
    !
    implicit none
    !
    PetscInt              :: data_type
    PetscReal, pointer    :: data_1D(:)
    !
    PetscInt              :: ii,jj,kk
    PetscInt              :: count
    PetscReal             :: x
    PetscReal             :: P     , dP_dx, d2P_dx2
    PetscReal             :: p0    , dp0_dx
    PetscReal             :: k     , dk_dx
    PetscReal             :: ds0   , ds0_dx
    PetscReal             :: dsr   , dsr_dx
    PetscReal             :: kr    , dkr_dx, dkr_dse
    PetscReal             :: rho   , drho_dP, d2rho_dP2, drho_dx, d2rho_dx2
    PetscReal             :: mu    , dmu_dP
    PetscReal             :: dse_dP, dse_dp0
    PetscReal             :: m
    PetscReal , parameter :: g = GRAVITY_CONSTANT
    PetscReal             :: dden_dT, dmu_dT
    PetscReal             :: sat_res, se
    PetscReal             :: dkr_dP
    type(saturation_params_type) :: satParam
    !-----------------------------------------------------------------------

    jj = 1; kk = 1;

    select case(data_type)
    case (DATA_PRESSURE)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_pressure_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_INITIAL_PRESSURE)
       P = 0.d0
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_pressure_or_deriv(x, val=data_1D(ii))
          P = P + 1.d0/nx*data_1D(ii)
       end do
       data_1D(:) = P
    case (DATA_POROSITY)
       do ii = 1, nx
          data_1D(ii) = 0.d0
       end do

    case (DATA_PERMEABILITY)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_permeability_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_SATFUNC_ALPHA)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_alpha_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_SATFUNC_LAMBDA)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_lambda_or_deriv(x, val=data_1D(ii))
       end do

    case (DATA_RES_SAT)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)
          call compute_residualsat_or_deriv(x, val=data_1D(ii))
       end do

     case (DATA_PRESSURE_BC)
        count = 0

        ii = 1; jj = 1; kk = 1;
        x     = soil_xc_3d(ii,jj,kk)-dx/2.d0
        count = count + 1;
        call compute_pressure_or_deriv(x,val=data_1D(count))

        ii = nx; jj = 1; kk = 1;
        x     = soil_xc_3d(ii,jj,kk)+dx/2.d0
        count = count + 1;
        call compute_pressure_or_deriv(x,val=data_1D(count))

    case (DATA_LIQUID_SATURATION)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)

          call compute_alpha_or_deriv       (x, val=p0, dval_dx=dp0_dx)
          call compute_lambda_or_deriv      (x, val=m)
          call compute_pressure_or_deriv    (x, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2)
          call compute_residualsat_or_deriv (x, val=sat_res)

          call SatFunc_Set_VG(satParam, sat_res, p0, m)
          call SatFunc_PressToSat(satParam, P, se, dse_dP)
          data_1D(ii) = se

       end do

    case (DATA_MASS_SOURCE)
       do ii = 1, nx
          x     = soil_xc_3d(ii,jj,kk)

          call compute_permeability_or_deriv(x, val=k , dval_dx=dk_dx)
          call compute_alpha_or_deriv       (x, val=p0, dval_dx=dp0_dx)
          call compute_lambda_or_deriv      (x, val=m)
          call compute_pressure_or_deriv    (x, val=P, dval_dx=dP_dx, d2val_dx2=d2P_dx2)
          call compute_residualsat_or_deriv (x, val=sat_res)
          
          call Viscosity(P, 298.15d0, mu, dmu_dP, dmu_dT)
          call Density(P, 298.15d0, DENSITY_TGDPB01,  rho, drho_dP, dden_dT)

          rho     = rho     * FMWH2O
          drho_dP = drho_dP * FMWH2O
          d2rho_dP2 = 0.d0

          call SatFunc_Set_VG(satParam, sat_res, p0, m)
          call SatFunc_PressToSat(satParam, P, se, dse_dP)
          call SatFunc_PressToRelPerm(satParam, P, 1.d0, kr, dkr_dP)
          
          dkr_dse   = &
               0.5d0 * se**(-0.5d0) *( 1.d0 - (1.d0 - se **(1.d0/m))**m)**2.d0 + &
               se**(0.5d0) * 2.d0   *( 1.d0 - (1.d0 - se **(1.d0/m))**m) * (1.d0 - se**(1.d0/m)) * se**(1.d0/m - 1.d0)
          dse_dp0   = 0.d0
          
          dkr_dx    = dkr_dP * dP_dx + dkr_dse * dse_dp0 * dp0_dx
          drho_dx   = drho_dP * dP_dx
          d2rho_dx2 = d2rho_dP2 * dP_dx + drho_dP * d2P_dx2

          data_1D(ii) = &
               -((k*kr/mu)*drho_dx + (rho*kr/mu)*dk_dx + (rho*k/mu)*dkr_dx)*(dP_dx) &
               -(rho*k*kr/mu)*(d2P_dx2)
          data_1D(ii) = data_1D(ii)*dx
          
       end do

    case default
       write(*,*)data_type
       write(iulog,*)'Unsupported data type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine set_variable_for_steady_state_1D

 end module vsfm_mms_steady_state_soil_only_1D
