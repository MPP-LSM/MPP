module CanopyLeafTemperatureAuxType

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

  type, public :: cleaf_temp_auxvar_type

     ! primary unknown independent variable
     PetscReal :: temperature     ! Vegetation temperature from previous timestep (K)

     PetscReal :: air_temperature ! Air temperature for previous timestep (K)
     PetscReal :: pref            ! Atmospheric pressure (Pa)
     PetscReal :: vcan            ! Water vapor for previous timestep (mol/mol)
     PetscReal :: cpair           ! Specific heat of air at constant pressure (J/mol/K)
     PetscReal :: gbh             ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
     PetscReal :: gbv             ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
     PetscReal :: gs              ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
     PetscReal :: rn              ! Leaf net radiation (W/m2 leaf)
     PetscReal :: cp              ! Leaf heat capacity (J/m2 leaf/K)
     PetscReal :: fssh            ! Sunlit or shaded fraction of canopy layer
     PetscReal :: dpai            ! Layer plant area index (m2/m2)
     PetscReal :: fwet            ! Fraction of plant area index that is wet
     PetscReal :: fdry            ! Fraction of plant area index that is green and dry

   contains
     procedure, public :: Init => CLeafTempAuxVarInit
  end type cleaf_temp_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine CLeafTempAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(cleaf_temp_auxvar_type)   :: this

    this%temperature = 0.d0

    this%air_temperature = 0.d0
    this%pref            = 0.d0
    this%vcan            = 0.d0
    this%gbh             = 0.d0
    this%gbv             = 0.d0
    this%gs              = 0.d0
    this%rn              = 0.d0
    this%cp              = 0.d0
    this%fssh            = 0.d0
    this%dpai            = 0.d0
    this%fwet            = 0.d0
    this%fdry            = 0.d0

  end subroutine CLeafTempAuxVarInit

#endif

end module CanopyLeafTemperatureAuxType
