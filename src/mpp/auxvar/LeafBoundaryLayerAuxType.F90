module LeafBoundaryLayerAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl      , only : iulog
  use mpp_abortutils  , only : endrun
  use mpp_shr_log_mod , only : errMsg => shr_log_errMsg

  implicit none
  private

  type, public :: leaf_bnd_lyr_auxvar_type

     ! Primary indepedent variables
     PetscReal :: gbh    ! boundary layer conductance for heat (mol/m2 leaf/s)      
     PetscReal :: gbv    ! boundary layer conductance for water (mol H2O/m2 leaf/s)
     PetscReal :: gbc    ! boundary layer conductance for CO2 (mol CO2/m2 leaf/s)

     PetscReal :: patm   ! atmospheric pressure (Pa)
     PetscReal :: rhomol ! Molar density (mol/m3)
     PetscReal :: wind   ! wind speed (m/s)
     PetscReal :: tair   ! air temperature (K)
     PetscReal :: tleaf  ! leaf temperature (K)
     PetscReal :: dleaf  ! leaf dimension (m)
   contains
     procedure, public :: Init => LeafBndLyrInit

  end type leaf_bnd_lyr_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine LeafBndLyrInit(this)
    !
    implicit none
    !
    class(leaf_bnd_lyr_auxvar_type) :: this

    this%gbh    = 0.d0
    this%gbv    = 0.d0
    this%gbc    = 0.d0

    this%patm   = 0.d0
    this%rhomol = 0.d0
    this%wind   = 0.d0
    this%tair   = 0.d0
    this%tleaf  = 0.d0
    this%dleaf  = 0.d0

  end subroutine LeafBndLyrInit

#endif

end module LeafBoundaryLayerAuxType
