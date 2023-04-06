module LongwaveAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use petscsys
  !
  implicit none

  private

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
     PetscReal, pointer :: leaf_fssh(:)        ! fraction of leaf
     PetscReal, pointer :: leaf_dpai(:)        ! plant area index
     

     PetscBool          :: is_soil             ! PETSC_TRUE if the grid cell is a soil grid cell
     PetscReal          :: ground_temperature    ! temperature of soil
     PetscReal          :: ground_emiss          ! ground_emissivity

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
    allocate(this%leaf_fssh(   nleaf));
    allocate(this%leaf_dpai(       nleaf));

    this%Idn                 = 0.d0
    this%Idn                 = 0.d0
    this%Idn                 = 0.d0

    this%trans               = 0.d0

    this%nleaf               = nleaf
    this%leaf_rho            = 0.d0
    this%leaf_tau            = 0.d0
    this%leaf_emiss          = 0.d0
    this%leaf_temperature(:) = 0.d0
    this%leaf_fssh(:)    = 0.d0
    this%leaf_dpai(:)        = 0.d0

    this%is_soil             = PETSC_FALSE
    this%ground_temperature    = 0.d0
    this%ground_emiss          = 0.d0

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

       aa = (1.d0 - this%trans) * this%leaf_tau + this%trans
       bb = (1.d0 - this%trans) * this%leaf_rho
       this%e = aa/bb

       this%f = 1.d0 - this%ground_emiss
       this%rad_source = STEFAN_BOLTZMAN_CONSTANT * this%ground_emiss * (this%ground_temperature ** 4.d0)

    else

       aa = (1.d0 - this%trans) * this%leaf_tau + this%trans
       bb = (1.d0 - this%trans) * this%leaf_rho

       this%f = bb - aa*aa/bb

       this%e = aa/bb

       this%rad_source = 0.d0
       do ileaf = 1, this%nleaf

          ! frac * (1 - tau) * emiss * sigma * T^4
          cc = &
               this%leaf_emiss              * &
               STEFAN_BOLTZMAN_CONSTANT     * &
               this%leaf_temperature(ileaf)**4.d0 * &
               this%leaf_fssh(ileaf)

          this%rad_source = this%rad_source + cc
       end do
       this%rad_source = this%rad_source * (1.d0 - this%trans)

    end if

  end subroutine LongwaveAuxVarCompute

#endif

end module LongwaveAuxType
