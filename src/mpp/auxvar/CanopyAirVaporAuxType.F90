module CanopyAirVaporAuxType

#ifdef USE_PETSC_LIB

  ! !USES:
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include <petsc/finclude/petsc.h>

  type, public :: cair_vapor_auxvar_type

     ! primary unknown independent variable
     PetscReal          :: qair                ! Water vapor (mol/mol)

     PetscReal          :: qair_prev           ! Water vapor for previous timestep (mol/mol)
     PetscReal          :: temperature         ! Air temperature (K)
     PetscReal          :: cpair               ! Specific heat of air at constant pressure (J/mol/K)
     PetscReal          :: rhomol              ! Molar density (mol/m3)
     PetscReal          :: pref                ! Atmospheric pressure (Pa)

     PetscBool          :: is_soil             ! PETSC_TRUE if the grid cell is a soil grid cell
     PetscReal          :: soil_rhg            ! Relative humidity of airspace at soil surface (fraction)
     PetscReal          :: soil_rn             ! Net radiation at ground (W/m2)
     PetscReal          :: soil_tk             ! Soil thermal conductivity (W/m/K)
     PetscReal          :: soil_dz             ! Soil layer depth (m)
     PetscReal          :: soil_resis          ! Soil evaporative resistance (s/m)
     PetscReal          :: soil_temperature    ! Soil temperature (K)

     PetscInt           :: num_leaves
     PetscReal, pointer :: gbv(:)              ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
     PetscReal, pointer :: leaf_temperature(:) ! Vegetation temperature from previous timestep (K)
     PetscReal, pointer :: leaf_gs(:)          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
     PetscReal, pointer :: leaf_fwet(:)        ! Fraction of plant area index that is wet
     PetscReal, pointer :: leaf_fdry(:)        ! Fraction of plant area index that is green and dry
     PetscReal, pointer :: leaf_fssh(:)        ! Sunlit or shaded fraction of canopy layer
     PetscReal, pointer :: leaf_dpai(:)        ! Layer plant area index (m2/m2)
   contains
     procedure, public :: Init     => CAirVaporAuxVarInit
     procedure, public :: PreSolve => CAirVaporAuxVarPreSolve
  end type cair_vapor_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine CAirVaporAuxVarInit(this, num_leaves)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(cair_vapor_auxvar_type)   :: this
    PetscInt :: num_leaves

    this%num_leaves = num_leaves
    allocate(this%gbv              (num_leaves))
    allocate(this%leaf_temperature (num_leaves))
    allocate(this%leaf_gs          (num_leaves))
    allocate(this%leaf_fwet        (num_leaves))
    allocate(this%leaf_fdry        (num_leaves))
    allocate(this%leaf_fssh        (num_leaves))
    allocate(this%leaf_dpai        (num_leaves))

    this%qair                = 0.d0

    this%temperature         = 0.d0
    this%cpair               = 0.d0
    this%rhomol              = 0.d0
    this%pref                = 0.d0

    this%is_soil             = PETSC_FALSE
    this%soil_rhg            = 0.d0
    this%soil_rn             = 0.d0
    this%soil_tk             = 0.d0
    this%soil_dz             = 0.d0
    this%soil_resis          = 0.d0

    this%gbv(:)              = 0.d0
    this%leaf_temperature(:) = 0.d0
    this%leaf_gs(:)          = 0.d0
    this%leaf_fwet(:)        = 0.d0
    this%leaf_fdry(:)        = 0.d0
    this%leaf_fssh(:)        = 0.d0
    this%leaf_dpai(:)        = 0.d0

  end subroutine CAirVaporAuxVarInit

  !------------------------------------------------------------------------
  subroutine CAirVaporAuxVarPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Make a copy of solution from previous timestep
    !
    implicit none
    !
    ! !ARGUMENTS
    class(cair_vapor_auxvar_type)   :: this

    this%qair_prev = this%qair

  end subroutine CAirVaporAuxVarPreSolve
#endif

end module CanopyAirVaporAuxType
