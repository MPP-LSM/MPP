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
     PetscReal          :: leaf_dpai            ! layer plant area index (m^2/m^2) (size = nleaf)
     PetscReal, pointer :: leaf_fssh(:)         ! fraction of leaf(size = nleaf)
     

     PetscBool          :: is_soil              ! TRUE if the grid cell is a soil grid cell
     PetscReal, pointer :: soil_albedo_b(:)     ! beam ground albedo (size = nleaf*nband)
     PetscReal, pointer :: soil_albedo_d(:)     ! diffuse ground albedo (size = nleaf*nband)

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
     procedure, public :: ComputePostSolve => ShortwaveComputePostSolve

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
    !PetscInt, intent (in)        :: nleaf
    PetscInt        :: nleaf
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
    allocate(this%leaf_fssh     (nleaf       ))

    allocate(this%soil_albedo_b     (nband       ))
    allocate(this%soil_albedo_d     (nband       ))

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
    this%leaf_dpai            = 0.d0
    this%leaf_fssh(:)     = 0.d0

    this%is_soil              = PETSC_FALSE
    this%soil_albedo_b(:)     = 0.d0
    this%soil_albedo_d(:)     = 0.d0

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

          this%f(iband)   = this%soil_albedo_b(iband)
          this%rad_source(iband) = this%Iskyb(iband) * this%leaf_tbcum * this%soil_albedo_d(iband)
       end do

    else

       do iband = 1, this%nband

          aa = (1.d0 - this%leaf_td) * this%leaf_rho(iband)
          bb = (1.d0 - this%leaf_td) * this%leaf_tau(iband) + this%leaf_td

          this%f(iband) = aa - bb*bb/aa;
          this%e(iband) = bb/aa
          this%rad_source(iband) = this%Iskyb(iband) * this%leaf_tbcum * (1.d0 - this%leaf_tb)

       end do

    end if

  end subroutine ShortwaveAuxVarCompute

  !------------------------------------------------------------------------
  subroutine ShortwaveComputePostSolve(this)
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
          !this%Ig(iband) = (1.d0 - this%soil_albedo_d(iband)) * this%Idn(iband) + (1.d0 - this%albedo_b(iband)) * this%leaf_tbcum
       end do

    else

       do iband = 1, this%nband
          !!this%Icd(iband) = (this%Idn(iband) - this%Iup() ) * (1.0d - this%tau(iband) * (1.d0 - this%leaf_omega(iband))
          !this%Icb(iband) = this%Iskyb(iband) * this%leaf_tbcum * (1.0d - this%tau(iband) * (1.d0 - this%leaf_omega(iband))
       end do

    end if

  end subroutine ShortwaveComputePostSolve

#endif

end module ShortwaveAuxType
