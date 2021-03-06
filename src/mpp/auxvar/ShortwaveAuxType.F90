module ShortwaveAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use petscsys
  !
  implicit none

  private

  type, public :: shortwave_auxvar_type

     ! Primary indepedent variables
     PetscReal, pointer :: Iup(:)               ! upward longwave flux above layer (size = nband)
     PetscReal, pointer :: Idn(:)               ! downward longwave flux onto layer (size = nband)

     PetscInt           :: nband

     PetscReal, pointer :: Iskyb(:)             ! atmospheric direct beam solar radiation (W/m2)
     PetscReal, pointer :: Iskyd(:)             ! atmospheric diffuse solar radiation (W/m2)

     PetscInt           :: nleaf                ! number of leaves (= 1 for a single leaf; = 2 for sunlit and shaded leaf)
     PetscReal, pointer :: leaf_rho(:)          ! leaf reflectance (size = nband)
     PetscReal, pointer :: leaf_tau(:)          ! leaf transmittance (size = nband)
     PetscReal, pointer :: leaf_omega(:)        ! leaf scattering coefficient (size = nband)
     PetscReal          :: leaf_td              ! exponential transmittances of diffuse radiation through a single leaf
     PetscReal          :: leaf_tb              ! exponential transmittances of direct radiation through a single leaf
     PetscReal          :: leaf_tbcum           ! cumulative exponential transmittance of direct beam onto a canopy layer
     PetscReal          :: leaf_dlai            ! layer leaf area index (m^2/m^2) (size = nleaf)
     PetscReal, pointer :: leaf_fraction(:)     ! fraction of leaf(size = nleaf)
     

     PetscBool          :: is_soil              ! TRUE if the grid cell is a soil grid cell
     PetscReal, pointer :: soil_albedo(:)       ! ground albedo (size = nleaf*nband)

     ! entries of the matrix to setup the matrix and rhs of the linear system
     PetscReal, pointer :: f(:)                 ! (size = nband)
     PetscReal, pointer :: e(:)                 ! (size = nband)
     PetscReal, pointer :: rad_source(:)        ! (size = nband)

     ! Post-processing
     PetscReal, pointer :: Iabs_leaf(:)         ! radiation absorbed by leaf [W/m^2_leaf] (size: nleaf * nband)
     PetscReal, pointer :: Iabs_veg(:)          ! radiation absorbed by vegetation [W/m^2_soil] (size: nband)
     PetscReal, pointer :: Iabs_veg_leaftype(:) ! radiation absorbed by vegetation for different leaf types [W/m^2_soil] (size: nleaf * nband)
     PetscReal, pointer :: Iabs_soil(:)         ! radiation absorbed by soil [W/m^2_soil] (size: nband)

   contains

     procedure, public :: Init          => ShortwaveAuxVarInit
     procedure, public :: AuxVarCompute => ShortwaveAuxVarCompute

  end type shortwave_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine ShortwaveAuxVarInit(this, nleaf)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(shortwave_auxvar_type) :: this
    PetscInt, intent (in)        :: nleaf
    !
    PetscInt                     :: nband

    nband = 2
    allocate(this%Iup               (nband       ))
    allocate(this%Idn               (nband       ))

    allocate(this%Iskyb             (nband       ))
    allocate(this%Iskyd             (nband       ))

    allocate(this%leaf_rho          (nband       ))
    allocate(this%leaf_tau          (nband       ))
    allocate(this%leaf_omega        (nband       ))
    allocate(this%leaf_fraction     (nleaf       ))

    allocate(this%soil_albedo       (nleaf*nband ))

    allocate(this%f                 (nband       ))
    allocate(this%e                 (nband       ))
    allocate(this%rad_source        (nband       ))

    allocate(this%Iabs_leaf         (nleaf*nband ))
    allocate(this%Iabs_veg          (nband       ))
    allocate(this%Iabs_veg_leaftype (nleaf*nband ))
    allocate(this%Iabs_soil         (nband       ))


    this%Iup(:)               = 0.d0
    this%Idn(:)               = 0.d0

    this%nband                = nband
    this%nleaf                = nleaf

    this%leaf_rho(:)          = 0.d0
    this%leaf_tau(:)          = 0.d0
    this%leaf_omega(:)        = 0.d0
    this%leaf_tb              = 0.d0
    this%leaf_td              = 0.d0
    this%leaf_tbcum           = 0.d0
    this%leaf_dlai            = 0.d0
    this%leaf_fraction(:)     = 0.d0

    this%is_soil              = PETSC_FALSE
    this%soil_albedo(:)       = 0.d0

    this%f(:)                 = 0.d0
    this%e(:)                 = 0.d0
    this%rad_source(:)        = 0.d0

    this%Iabs_leaf(:)         = 0.d0
    this%Iabs_veg(:)          = 0.d0
    this%Iabs_veg_leaftype(:) = 0.d0
    this%Iabs_soil(:)         = 0.d0

  end subroutine ShortwaveAuxVarInit

  !------------------------------------------------------------------------
  subroutine ShortwaveAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants , only : STEFAN_BOLTZMAN_CONSTANT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(shortwave_auxvar_type) :: this
    !
    PetscInt  :: iband
    PetscReal :: aa, bb, cc

    if (this%is_soil) then

       do iband = 1, this%nband
          this%e(iband) = 0.d0

          this%f(iband)   = this%soil_albedo(iband)
          this%rad_source = this%Iskyb(iband) * this%leaf_tbcum * this%soil_albedo(iband)
       end do

    else

       do iband = 1, this%nband

          aa = (1.d0 - this%leaf_td) * this%leaf_rho(iband)
          bb = (1.d0 - this%leaf_td) * this%leaf_tau(iband) + this%leaf_td

          this%f(iband) = aa - bb**2.d0/aa;
          this%e(iband) = bb/aa
          this%rad_source(iband) = this%Iskyb(iband) * this%leaf_tbcum * (1.d0 - this%leaf_tb)

       end do

    end if

  end subroutine ShortwaveAuxVarCompute

#endif

end module ShortwaveAuxType
