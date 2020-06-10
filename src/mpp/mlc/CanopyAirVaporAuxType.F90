module CanopyAirVaporAuxType

#ifdef USE_PETSC_LIB

  ! !USES:
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include <petsc/finclude/petsc.h>

  type, public :: cair_vapor_auxvar_type

     ! primary unknown independent variable
     PetscReal          :: water_vapor         ! Water vapor for previous timestep (mol/mol)

     PetscReal          :: temperature         !
     PetscReal          :: gbv                 ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
     PetscReal          :: cpair               ! Specific heat of air at constant pressure (J/mol/K)
     PetscReal          :: rhomol              ! Molar density (mol/m3)
     PetscReal          :: pref                ! Atmospheric pressure (Pa)

     PetscReal          :: soil_rhg            ! Relative humidity of airspace at soil surface (fraction)
     PetscReal          :: soil_resis          ! Soil evaporative resistance (s/m)

     PetscInt           :: nleaf
     PetscReal, pointer :: leaf_temperature(:) ! Vegetation temperature from previous timestep (K)
     PetscReal, pointer :: leaf_gs(:)          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
     PetscReal, pointer :: leaf_fwet(:)        ! Fraction of plant area index that is wet
     PetscReal, pointer :: leaf_fdry(:)        ! Fraction of plant area index that is green and dry
     PetscReal, pointer :: leaf_fssh(:)        ! Sunlit or shaded fraction of canopy layer
     PetscReal, pointer :: leaf_dpai(:)        ! Layer plant area index (m2/m2)
   contains
     procedure, public :: Init => CAirVaporAuxVarInit
  end type cair_vapor_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine CAirVaporAuxVarInit(this, nleaf)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(cair_vapor_auxvar_type)   :: this
    PetscInt :: nleaf

    this%nleaf = nleaf
    allocate(this%leaf_temperature(nleaf))
    allocate(this%leaf_gs(nleaf))
    allocate(this%leaf_fwet(nleaf))
    allocate(this%leaf_fdry(nleaf))
    allocate(this%leaf_fssh(nleaf))
    allocate(this%leaf_dpai(nleaf))

    this%water_vapor         = 0.d0

    this%temperature         = 0.d0
    this%gbv                 = 0.d0
    this%cpair               = 0.d0
    this%rhomol              = 0.d0
    this%pref                = 0.d0
    this%soil_rhg            = 0.d0
    this%soil_resis          = 0.d0

    this%leaf_temperature(:) = 0.d0
    this%leaf_gs(:)          = 0.d0
    this%leaf_fwet(:)        = 0.d0
    this%leaf_fdry(:)        = 0.d0
    this%leaf_fssh(:)        = 0.d0
    this%leaf_dpai(:)        = 0.d0

  end subroutine CAirVaporAuxVarInit

#endif

end module CanopyAirVaporAuxType
