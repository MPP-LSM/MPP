module LongwaveAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use AuxVarType          , only : auxvar_base_type
  use petscsys
  !
  implicit none

  private

  !type, public, extends(auxvar_base_type) :: longwave_auxvar_type
  type, public :: longwave_auxvar_type

     ! Primary indepedent variables
     PetscReal          :: Iup                 ! upward longwave flux above layer
     PetscReal          :: Idn                 ! downward longwave flux onto layer
     PetscReal          :: Iabs                ! absorbed radiation by layer

     PetscReal          :: trans               ! transmittance of diffuse radiation through a single leaf layer

     PetscInt           :: nleaf               ! number of leaves
     PetscReal          :: leaf_rho            ! leaf reflectance
     PetscReal          :: leaf_tau            ! leaf transmittance
     PetscReal          :: leaf_emiss          ! leaf emissivity
     PetscReal, pointer :: leaf_temperature(:) ! leaf temperature
     PetscReal, pointer :: leaf_fraction(:)    ! fraction of leaf
     PetscReal, pointer :: leaf_lai(:)         ! leaf area index
     

     PetscBool          :: is_soil             ! PETSC_TRUE if the grid cell is a soil grid cell
     PetscReal          :: soil_temperature    ! temperature of soil
     PetscReal          :: soil_emiss          ! soil_emissivity

     ! entries of the matrix to setup the matrix and rhs of the linear system
     PetscReal          :: f
     PetscReal          :: e
     PetscReal          :: rad_source
   contains

     procedure, public :: Init          => LongwaveAuxVarInit
     procedure, public :: AuxVarCompute => LongwaveAuxVarCompute

  end type longwave_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine LongwaveAuxVarInit(this, nleaf)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(longwave_auxvar_type) :: this
    PetscInt                    :: nleaf

    allocate(this%leaf_temperature(nleaf));
    allocate(this%leaf_fraction(   nleaf));
    allocate(this%leaf_lai(        nleaf));

    this%Idn                 = 0.d0
    this%Idn                 = 0.d0
    this%Idn                 = 0.d0

    this%trans               = 0.d0

    this%nleaf               = nleaf
    this%leaf_rho            = 0.d0
    this%leaf_tau            = 0.d0
    this%leaf_emiss          = 0.d0
    this%leaf_temperature(:) = 0.d0
    this%leaf_fraction(:)    = 0.d0
    this%leaf_lai(:)         = 0.d0

    this%is_soil             = PETSC_FALSE
    this%soil_temperature    = 0.d0
    this%soil_emiss          = 0.d0

    this%f                   = 0.d0
    this%e                   = 0.d0
    this%rad_source          = 0.d0

  end subroutine LongwaveAuxVarInit

  !------------------------------------------------------------------------
  subroutine LongwaveAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants , only : STEFAN_BOLTZMAN_CONSTANT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(longwave_auxvar_type) :: this
    !
    PetscInt  :: ileaf
    PetscReal :: aa, bb, cc

    if (this%is_soil) then

       aa = this%trans + (1 - this%trans) * this%leaf_tau
       bb = (1.d0 - this%trans) * this%leaf_rho
       this%e = aa/bb

       this%f = 1.d0 - this%soil_emiss
       this%rad_source = STEFAN_BOLTZMAN_CONSTANT * this%soil_emiss * (this%soil_temperature ** 4.d0)

    else

       aa = this%trans + (1 - this%trans) * this%leaf_tau
       bb = (1.d0 - this%trans) * this%leaf_rho

       this%f = (1.d0 - this%trans) * this%leaf_rho - &
            aa**2.d0/bb

       this%e = aa/bb

       this%rad_source = 0.d0
       do ileaf = 1, this%nleaf

          ! frac * (1 - tau) * emiss * sigma * T^4
          cc = this%leaf_fraction(ileaf)    * &
               (1.d0 - this%trans)          * &
               STEFAN_BOLTZMAN_CONSTANT     * &
               this%leaf_emiss              * &
               this%leaf_temperature(ileaf)**4.d0

          this%rad_source = this%rad_source + cc
       end do

    end if

  end subroutine LongwaveAuxVarCompute

#endif

end module LongwaveAuxType
