module GoveqnLeafBoundaryLayer

#ifdef USE_PETSC_LIB
  
#include <petsc/finclude/petsc.h>
  
  ! !USES:
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use LeafBoundaryLayerAuxType  , only : leaf_bnd_lyr_auxvar_type
  use GoverningEquationBaseType , only : goveqn_base_type
  use petscvec
  use petscmat
  use petscsys
  !
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_leaf_bnd_lyr_type

     type(leaf_bnd_lyr_auxvar_type), pointer ::  aux_vars_in(:)

   contains

     procedure, public :: Setup                     => LeafBndLyrSetup
     procedure, public :: AllocateAuxVars           => LeafBndLyrAllocateAuxVars
     procedure, public :: PreSolve                  => LeafBndLyrPreSolve
     procedure, public :: ComputeRHS                => LeafBndLyrComputeRHS
     procedure, public :: ComputeOperatorsDiag      => LeafBndLyrComputeOperatorsDiag
     procedure, public :: GetRValues                => LeafBndLyrGetRValues
     procedure, public :: SavePrimaryIndependentVar => LeafBndLyrSavePrmIndepVar

  end type goveqn_leaf_bnd_lyr_type

contains

  !------------------------------------------------------------------------
  subroutine LeafBndLyrSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_LEAF_BND_LAYER
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_leaf_bnd_lyr_type) :: this

    call this%Create()

    this%name  = ""
    this%itype = GE_LEAF_BND_LAYER
    this%dof   = 3

    nullify(this%aux_vars_in)

  end subroutine LeafBndLyrSetup

  !------------------------------------------------------------------------
  subroutine LeafBndLyrAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates memory for storing auxiliary variables associated with:
    !   + Internal control volumes,
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_leaf_bnd_lyr_type) :: this
    !

    ! Allocate memory and initialize aux vars: For internal connections
    allocate(this%aux_vars_in(this%mesh%ncells_all))

  end subroutine LeafBndLyrAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine LeafBndLyrPreSolve(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : RGAS
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_leaf_bnd_lyr_type) :: this
    !
    PetscInt :: icell
    type(leaf_bnd_lyr_auxvar_type), pointer :: avars(:)

    avars => this%aux_vars_in

    do icell = 1, this%mesh%ncells_all
        avars(icell)%rhomol = avars(icell)%patm / (RGAS * avars(icell)%tair)
    enddo
  end subroutine LeafBndLyrPreSolve

  !------------------------------------------------------------------------
  subroutine LeafBndLyrSavePrmIndepVar (this, x)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use CanopyTurbulenceAuxType, only : canopy_turbulence_auxvar_type 
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_leaf_bnd_lyr_type) :: this
    Vec :: x
    !
    PetscScalar, pointer :: x_p(:)
    PetscInt             :: ghosted_id, size
    PetscErrorCode       :: ierr
    
    call VecGetLocalSize(x, size, ierr); CHKERRQ(ierr)

    if (size /= this%mesh%ncells_local*this%dof) then
       call endrun(msg="ERROR size of vector /= number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    end if

    call VecGetArrayF90(x, x_p, ierr); CHKERRQ(ierr)

    do ghosted_id = 1, this%mesh%ncells_local
       this%aux_vars_in(ghosted_id)%gbh = x_p((ghosted_id-1)*3 + 1)
       this%aux_vars_in(ghosted_id)%gbv = x_p((ghosted_id-1)*3 + 2)
       this%aux_vars_in(ghosted_id)%gbc = x_p((ghosted_id-1)*3 + 3)
   end do

    call VecRestoreArrayF90(x, x_p, ierr); CHKERRQ(ierr)

  end subroutine LeafBndLyrSavePrmIndepVar

  !------------------------------------------------------------------------
  subroutine LeafBndLyrComputeRHS(this, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    !
    use WaterVaporMod, only : SatVap
    use MultiPhysicsProbConstants , only : GRAVITY_CONSTANT, TFRZ
    use MultiPhysicsProbConstants , only : VISC_0C, MOD_DIFF_CO2_OC, MOD_DIFF_H2O_OC, MOD_DIFF_HEAT_OC
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_leaf_bnd_lyr_type) :: this
    Vec                           :: B
    PetscErrorCode                :: ierr
    !
    PetscInt                                :: icell, ileaf
    PetscScalar                   , pointer :: b_p(:)
    PetscReal                               :: qsat, si, gleaf, gleaf_et, lambda
    PetscReal :: b1, factor
    PetscReal :: visc, Dh, Dv, Dc
    PetscReal :: Re, Pr, Scv, Scc, Gr, dT
    PetscReal :: Nu_lam, Shv_lam, Shc_lam
    PetscReal :: Nu_turb, Shv_turb, Shc_turb
    PetscReal :: Nu_free, Shv_free, Shc_free
    PetscReal :: Nu_forced, Shv_forced, Shc_forced
    PetscReal :: Nu, Shv, Shc
    type(leaf_bnd_lyr_auxvar_type), pointer :: avars(:)

    b1 = 1.5d0 ! empirical correction factor for Nu and Sh
    
    avars => this%aux_vars_in

    call VecGetLocalSize(B, icell, ierr);
    call VecGetArrayF90(B, b_p, ierr)

    do icell = 1, this%mesh%ncells_local
        ! adjust diffusivity for temperature and pressure
        factor = 101325.d0 / avars(icell)%patm * (avars(icell)%tair/TFRZ)**1.81d0

        visc = VISC_0C * factor
        Dh = MOD_DIFF_HEAT_OC * factor
        Dv = MOD_DIFF_H2O_OC * factor
        Dc = MOD_DIFF_CO2_OC * factor

        ! Dimensionless numbers
        Re = avars(icell)%wind * avars(icell)%dleaf / visc
        Pr = visc/Dh
        Scv = visc/Dv
        Scc = visc/Dc

        ! Grashof number
        dT = max(avars(icell)%tleaf - avars(icell)%tair, 0.d0)
        Gr = GRAVITY_CONSTANT * (avars(icell)%dleaf**3.d0) * dT /( avars(icell)%tair * visc**2.d0 )

        ! Choose correct flow regime for forced convection
        Nu_lam  = b1 * 0.66d0 *  Pr**0.33d0 * Re**0.5d0;     ! Nusselt number
        Shv_lam = b1 * 0.66d0 * Scv**0.33d0 * Re**0.5d0;     ! Sherwood number, H2O
        Shc_lam = b1 * 0.66d0 * Scc**0.33d0 * Re**0.5d0;     ! Sherwood number, CO2
        
        ! Forced convection - turbulent flow
        Nu_turb  = b1 * 0.036d0 *  Pr**0.33 * Re**0.8d0;   ! Nusselt number
        Shv_turb = b1 * 0.036d0 * Scv**0.33 * Re**0.8d0;   ! Sherwood number, H2O
        Shc_turb = b1 * 0.036d0 * Scc**0.33 * Re**0.8d0;   ! Sherwood number, CO2
        
        ! Choose correct flow regime for forced convection        
        Nu_forced = max(Nu_lam, Nu_turb);
        Shv_forced = max(Shv_lam, Shv_turb);
        Shc_forced = max(Shc_lam, Shc_turb);
        
        ! Free convection
        Nu_free  = 0.54d0 *  Pr**0.25d0 * Gr**0.25d0;        ! Nusselt number
        Shv_free = 0.54d0 * Scv**0.25d0 * Gr**0.25d0;        ! Sherwood number, H2O
        Shc_free = 0.54d0 * Scc**0.25d0 * Gr**0.25d0;        ! Sherwood number, CO2
        
        ! Both forced and free convection regimes occur together
        Nu = Nu_forced + Nu_free;
        Shv = Shv_forced + Shv_free;
        Shc = Shc_forced + Shc_free;
        
        ! Boundary layer conductances (m/s)
        avars(icell)%gbh = Dh *  Nu / avars(icell)%dleaf;
        avars(icell)%gbv = Dv * Shv / avars(icell)%dleaf;
        avars(icell)%gbc = Dc * Shc / avars(icell)%dleaf;
        
        ! Convert conductance (m/s) to (mol/m2/s)
        avars(icell)%gbh = avars(icell)%gbh * avars(icell)%rhomol;
        avars(icell)%gbv = avars(icell)%gbv * avars(icell)%rhomol;
        avars(icell)%gbc = avars(icell)%gbc * avars(icell)%rhomol;

        b_p((icell-1)*3 + 1) = avars(icell)%gbh
        b_p((icell-1)*3 + 2) = avars(icell)%gbv
        b_p((icell-1)*3 + 3) = avars(icell)%gbc

    end do

    call VecRestoreArrayF90(B, b_p, ierr)

  end subroutine LeafBndLyrComputeRHS
  
  !------------------------------------------------------------------------

  subroutine LeafBndLyrComputeOperatorsDiag(this, A, B, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use mpp_varcon                , only : cnfac
    use MultiPhysicsProbConstants , only : COND_HEAT_FLUX
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use MultiPhysicsProbConstants , only : GE_THERM_SSW_TBASED
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    !
    use WaterVaporMod, only : SatVap
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_leaf_bnd_lyr_type) :: this
    Mat                           :: A
    Mat                           :: B
    PetscErrorCode                :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                             :: icell, row, col
    PetscReal                            :: value

    do icell = 1, this%mesh%ncells_local
       value = 1.d0

       row = (icell-1)*3; col = row;
       call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

       row = (icell-1)*3+1; col = row;
       call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

       row = (icell-1)*3+2; col = row;
       call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)

    end do

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  end subroutine LeafBndLyrComputeOperatorsDiag

   !------------------------------------------------------------------------
  subroutine LeafBndLyrGetRValues (this, auxvar_type, var_type, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    use MultiPhysicsProbConstants, only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants, only : AUXVAR_BC
    use MultiPhysicsProbConstants, only : AUXVAR_SS
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_leaf_bnd_lyr_type) :: this
    PetscInt                        :: auxvar_type
    PetscInt                        :: var_type
    PetscInt                        :: nauxvar
    PetscReal, pointer              :: var_values(:)

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call LeafBndLyrGetRValuesFromAuxVars(this%aux_vars_in, var_type, nauxvar, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine LeafBndLyrGetRValues

  !------------------------------------------------------------------------
  subroutine LeafBndLyrGetRValuesFromAuxVars (aux_var, var_type, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : VAR_LEAF_BDN_LYR_COND_HEAT
    use MultiPhysicsProbConstants , only : VAR_LEAF_BDN_LYR_COND_H2O
    use MultiPhysicsProbConstants , only : VAR_LEAF_BDN_LYR_COND_CO2
    !
    implicit none
    !
    ! !ARGUMENTS
    type(leaf_bnd_lyr_auxvar_type), pointer ::  aux_var(:)
    PetscInt                                :: var_type
    PetscInt                                :: nauxvar
    PetscReal                     , pointer :: var_values(:)
    !
    PetscInt :: iauxvar

    select case(var_type)
    case(VAR_LEAF_BDN_LYR_COND_HEAT)
       do iauxvar = 1,nauxvar
          var_values(iauxvar) = aux_var(iauxvar)%gbh
       end do

    case(VAR_LEAF_BDN_LYR_COND_H2O)
        do iauxvar = 1,nauxvar
           var_values(iauxvar) = aux_var(iauxvar)%gbv
        end do

    case(VAR_LEAF_BDN_LYR_COND_CO2)
        do iauxvar = 1,nauxvar
           var_values(iauxvar) = aux_var(iauxvar)%gbc
        end do
   case default
       write(iulog,*) 'CLeafTempGetRValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine LeafBndLyrGetRValuesFromAuxVars

#endif

end module GoveqnLeafBoundaryLayer
