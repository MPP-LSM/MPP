module CanopyAirTemperatureAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: cair_temp_auxvar_type

     ! primary unknown independent variable
     PetscReal          :: temperature         ! Air temperature profile for previous timestep (K)

     PetscReal          :: gbh                 ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
     PetscReal          :: cpair               ! Specific heat of air at constant pressure (J/mol/K)
     PetscReal          :: rhomol              ! Molar density (mol/m3)
     PetscReal          :: pref                ! Atmospheric pressure (Pa)

     PetscReal          :: water_vapor         ! Water vapor (mol/mol)

     PetscInt           :: nleaf               ! Number of canopy leaves
     PetscReal, pointer :: leaf_temperature(:) ! Vegetation from previous timestep (K)
     PetscReal, pointer :: leaf_gs(:)          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
     PetscReal, pointer :: leaf_fwet(:)        ! Fraction of plant area index that is wet
     PetscReal, pointer :: leaf_fdry(:)        ! Fraction of plant area index that is green and dry
     PetscReal, pointer :: leaf_fssh(:)        ! Sunlit or shaded fraction of canopy layer
     PetscReal, pointer :: leaf_dpai(:)        ! Layer plant area index (m2/m2)

     PetscBool          :: is_soil             ! PETSC_TRUE if the grid cell is a soil grid cell
     PetscReal          :: soil_rhg            ! Relative humidity of airspace at soil surface (fraction)
     PetscReal          :: soil_rn             ! Net radiation at ground (W/m2)
     PetscReal          :: soil_tk             ! Soil thermal conductivity (W/m/K)
     PetscReal          :: soil_dz             ! Soil layer depth (m)
     PetscReal          :: soil_resis          ! Soil evaporative resistance (s/m)
     PetscReal          :: soil_temperature    ! Soil temperature (K)
   contains

     procedure, public :: Init => CAirTempAuxVarInit

  end type cair_temp_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine CAirTempAuxVarInit(this, nleaf)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(cair_temp_auxvar_type) :: this
    PetscInt                     :: nleaf

    this%temperature = 0.d0
    this%gbh         = 0.d0
    this%cpair       = 0.d0
    this%rhomol      = 0.d0
    this%pref        = 0.d0

    this%water_vapor = 0.d0

    this%nleaf = nleaf
    allocate(this%leaf_temperature(nleaf))
    allocate(this%leaf_gs(nleaf))
    allocate(this%leaf_fwet(nleaf))
    allocate(this%leaf_fdry(nleaf))
    allocate(this%leaf_fssh(nleaf))
    allocate(this%leaf_dpai(nleaf))

    this%leaf_temperature(:) = 0.d0
    this%leaf_gs(:)          = 0.d0
    this%leaf_fwet(:)        = 0.d0
    this%leaf_fdry(:)        = 0.d0
    this%leaf_fssh(:)        = 0.d0
    this%leaf_dpai(:)        = 0.d0

    this%is_soil             = PETSC_FALSE
    this%soil_rhg            = 0.d0
    this%soil_rn             = 0.d0
    this%soil_tk             = 0.d0
    this%soil_dz             = 0.d0
    this%soil_resis          = 0.d0

  end subroutine CAirTempAuxVarInit

#endif

end module CanopyAirTemperatureAuxType
