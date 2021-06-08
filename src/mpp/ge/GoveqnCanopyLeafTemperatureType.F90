module GoveqnCanopyLeafTemperatureType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                   , only : iulog
  use mpp_abortutils               , only : endrun
  use mpp_shr_log_mod              , only : errMsg => shr_log_errMsg
  use GoverningEquationBaseType    , only : goveqn_base_type
  use CanopyLeafTemperatureAuxType , only : cleaf_temp_auxvar_type
  use petscvec
  use petscmat
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public, extends(goveqn_base_type) :: goveqn_cleaf_temp_type

     type(cleaf_temp_auxvar_type), pointer :: aux_vars_in(:)
     PetscInt                    , pointer :: Leaf2CAir(:)

   contains

     procedure, public :: Setup                     => CLeafTempSetup
     procedure, public :: AllocateAuxVars           => CLeafTempAllocateAuxVars
     procedure, public :: GetFromSoeAuxVarsCturb    => CLeafTempGetFromSoeAuxVarsCturb
     procedure, public :: SavePrimaryIndependentVar => CLeafTempSavePrmIndepVar
     procedure, public :: GetRValues                => CLeafTempGetRValues
     procedure, public :: SetRValues                => CLeafTempSetRValues
     procedure, public :: ComputeRHS                => CLeafTempComputeRHS
     procedure, public :: ComputeOperatorsDiag      => CLeafTempComputeOperatorsDiag
     procedure, public :: ComputeOperatorsOffDiag   => CLeafTempComputeOperatorsOffDiag
     procedure, public :: SetLeaf2CAirMap
     procedure, public :: SetDefaultLeaf2CAirMap

  end type goveqn_cleaf_temp_type

contains

  !------------------------------------------------------------------------
  subroutine CLeafTempSetup(this)
    !
    ! !DESCRIPTION:
    ! Default setup of governing equation
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : GE_CANOPY_LEAF_TEMP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cleaf_temp_type) :: this

    call this%Create()

    this%name  = ""
    this%itype = GE_CANOPY_LEAF_TEMP

    nullify(this%aux_vars_in)
    nullify(this%Leaf2CAir)

  end subroutine CLeafTempSetup

  !------------------------------------------------------------------------
  subroutine SetLeaf2CAirMap(this, CAirIdForLeaf, nCAirIdForLeaf)

    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cleaf_temp_type) :: this
    PetscInt, pointer :: CAirIdForLeaf(:)
    PetscInt          :: nCAirIdForLeaf
    !
    PetscInt :: ileaf

    if (this%mesh%ncells_all /= nCAirIdForLeaf) then
       call endrun(msg="size of Leaf2CAir /= number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    endif

    allocate(this%Leaf2CAir(nCAirIdForLeaf))

    do ileaf = 1, nCAirIdForLeaf
       this%Leaf2CAir(ileaf) = CAirIdForLeaf(ileaf) ! Canopy airspace id for the ileaf-th leaf
    enddo

  end subroutine SetLeaf2CAirMap

  !------------------------------------------------------------------------
  subroutine SetDefaultLeaf2CAirMap(this)

    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cleaf_temp_type) :: this
    !
    PetscInt :: ileaf

    allocate(this%Leaf2CAir(this%mesh%ncells_all))

    do ileaf = 1, this%mesh%ncells_all
       this%Leaf2CAir(ileaf) = ileaf  ! Canopy airspace id for the ileaf-th leaf
    enddo

  end subroutine SetDefaultLeaf2CAirMap

  !------------------------------------------------------------------------
  subroutine CLeafTempAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates memory for storing auxiliary variables associated with:
    !   + Internal control volumes,
    !   + Boundary condtions,
    !   + Source-sink condition.
    !
    ! !USES:
    use ConditionType, only : condition_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cleaf_temp_type) :: this
    !
    type(condition_type),pointer             :: cur_cond
    PetscInt                                 :: ncells_cond
    PetscInt                                 :: ncond
    PetscInt                                 :: icond

    ! Allocate memory and initialize aux vars: For internal connections
    allocate(this%aux_vars_in(this%mesh%ncells_all))
    do icond = 1,this%mesh%ncells_all
       call this%aux_vars_in(icond)%Init()
    enddo

  end subroutine CLeafTempAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine CLeafTempGetFromSoeAuxVarsCturb(this, cturb)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use CanopyTurbulenceAuxType, only : canopy_turbulence_auxvar_type 
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cleaf_temp_type)        :: this
    class(canopy_turbulence_auxvar_type) :: cturb
    !
    PetscInt :: icair, icell, level, iconn

    ! all layers
    do icair = 1, cturb%ncair
       do icell = 1, this%mesh%ncells_local
          level = icell
          this%aux_vars_in(icell)%cpair = cturb%cpair(icair)
          this%aux_vars_in(icell)%pref  = cturb%pref(icair)
       end do
    end do
  end subroutine CLeafTempGetFromSoeAuxVarsCturb

  !------------------------------------------------------------------------
  subroutine CLeafTempSavePrmIndepVar (this, x)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use CanopyTurbulenceAuxType, only : canopy_turbulence_auxvar_type 
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cleaf_temp_type) :: this
    Vec                           :: x
    !
    PetscScalar, pointer :: x_p(:)
    PetscInt             :: ghosted_id, size
    PetscErrorCode       :: ierr
    
    call VecGetLocalSize(x, size, ierr); CHKERRQ(ierr)

    if (size /= this%mesh%ncells_local) then
       call endrun(msg="ERROR size of vector /= number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    end if

    call VecGetArrayF90(x, x_p, ierr); CHKERRQ(ierr)

    do ghosted_id = 1, this%mesh%ncells_local
       this%aux_vars_in(ghosted_id)%temperature = x_p(ghosted_id)
    end do

    call VecRestoreArrayF90(x, x_p, ierr); CHKERRQ(ierr)

  end subroutine CLeafTempSavePrmIndepVar

  !------------------------------------------------------------------------
  subroutine CLeafTempGetRValues (this, auxvar_type, var_type, nauxvar, var_values)
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
    class(goveqn_cleaf_temp_type) :: this
    PetscInt                     :: auxvar_type
    PetscInt                     :: var_type
    PetscInt                     :: nauxvar
    PetscReal, pointer           :: var_values(:)

    if (nauxvar > this%mesh%ncells_all) then
      call endrun(msg="ERROR nauxvar exceeds the number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    endif

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call CLeafTempGetRValuesFromAuxVars(this%aux_vars_in, var_type, nauxvar, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CLeafTempGetRValues

  !------------------------------------------------------------------------
  subroutine CLeafTempSetRValues (this, auxvar_type, var_type, nauxvar, var_values)
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
    class(goveqn_cleaf_temp_type) :: this
    PetscInt                     :: auxvar_type
    PetscInt                     :: var_type
    PetscInt                     :: nauxvar
    PetscReal, pointer           :: var_values(:)

    if (nauxvar > this%mesh%ncells_all) then
      call endrun(msg="ERROR nauxvar exceeds the number of cells in the mesh "//errmsg(__FILE__, __LINE__))
    endif

    select case(auxvar_type)
    case(AUXVAR_INTERNAL)
       call CLeafTempSetRValuesFromAuxVars(this%aux_vars_in, var_type, nauxvar, var_values)
    case default
       write(*,*)'Unknown auxvar_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CLeafTempSetRValues

  !------------------------------------------------------------------------
  subroutine CLeafTempGetRValuesFromAuxVars (aux_var, var_type, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_LEAF_TEMPERATURE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(cleaf_temp_auxvar_type) , pointer :: aux_var(:)
    PetscInt                              :: var_type
    PetscInt                              :: nauxvar
    PetscReal                   , pointer :: var_values(:)
    !
    PetscInt :: iauxvar

    select case(var_type)
    case(VAR_LEAF_TEMPERATURE)
       do iauxvar = 1,nauxvar
          var_values(iauxvar) = aux_var(iauxvar)%temperature
       end do
    case default
       write(iulog,*) 'CLeafTempGetRValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine CLeafTempGetRValuesFromAuxVars
       
  !------------------------------------------------------------------------
  subroutine CLeafTempSetRValuesFromAuxVars (aux_var, var_type, nauxvar, var_values)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_WATER_VAPOR
    !
    implicit none
    !
    ! !ARGUMENTS
    type(cleaf_temp_auxvar_type) , pointer :: aux_var(:)
    PetscInt                              :: var_type
    PetscInt                              :: nauxvar
    PetscReal                   , pointer :: var_values(:)
    !
    PetscInt :: iauxvar

    select case(var_type)
    case(VAR_TEMPERATURE)
       do iauxvar = 1,nauxvar
          aux_var(iauxvar)%air_temperature = var_values(iauxvar)
       end do
    case(VAR_WATER_VAPOR)
       do iauxvar = 1,nauxvar
          aux_var(iauxvar)%water_vapor_canopy = var_values(iauxvar)
       end do
    case default
       write(iulog,*) 'CLeafTempGetSValuesFromAuxVars: Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select
       
  end subroutine CLeafTempSetRValuesFromAuxVars

  !------------------------------------------------------------------------
  subroutine CLeafTempComputeRHS(this, B, ierr)
    !
    ! !DESCRIPTION:
    ! Dummy subroutine for PETSc TS RSHFunction
    !
    use WaterVaporMod, only : SatVap
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cleaf_temp_type) :: this
    Vec                           :: B
    PetscErrorCode                :: ierr
    !
    PetscInt                                :: icell, ileaf
    PetscScalar                   , pointer :: b_p(:)
    class(cleaf_temp_auxvar_type) , pointer :: auxvar(:)
    PetscReal                               :: qsat, dqsat, esat, desat, gleaf, gleaf_et, lambda

    lambda = HVAP * MM_H2O

    auxvar => this%aux_vars_in

    call VecGetArrayF90(B, b_p, ierr)

    do icell = 1, this%mesh%ncells_local

       if (auxvar(icell)%dpai > 0.d0) then

          call SatVap(auxvar(icell)%temperature, esat, desat)
          qsat  = esat/this%aux_vars_in(icell)%pref
          dqsat = desat/this%aux_vars_in(icell)%pref

          gleaf = &
                auxvar(icell)%gs * auxvar(icell)%gbv / &
               (auxvar(icell)%gs + auxvar(icell)%gbv )

          gleaf_et = &
               gleaf             * auxvar(icell)%fdry + &
               auxvar(icell)%gbv * auxvar(icell)%fwet

          b_p(icell) = auxvar(icell)%rn &
               + auxvar(icell)%cp/this%dtime * auxvar(icell)%temperature &
               - lambda * (qsat - dqsat * auxvar(icell)%temperature) * gleaf_et

       end if
    end do
    call VecRestoreArrayF90(B, b_p, ierr)

  end subroutine CLeafTempComputeRHS

  !------------------------------------------------------------------------

  subroutine CLeafTempComputeOperatorsDiag(this, A, B, ierr)
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
    class(goveqn_cleaf_temp_type) :: this
    Mat                           :: A
    Mat                           :: B
    PetscErrorCode                :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt                             :: icell, ileaf, iconn, sum_conn
    PetscInt                             :: cell_id, cell_id_up, cell_id_dn
    PetscInt                             :: row, col
    PetscReal                            :: value, ga, dist
    PetscReal                            :: qsat, dqsat, esat, desat, gleaf, gleaf_et, lambda
    class(connection_set_type) , pointer :: cur_conn_set
    type(condition_type)       , pointer :: cur_cond

    lambda = HVAP * MM_H2O
    ! For interior and top cell
    do icell = 1, this%mesh%ncells_local
       row = icell-1; col = icell-1
       if (this%aux_vars_in(icell)%dpai > 0.d0) then

          call SatVap(this%aux_vars_in(icell)%temperature, esat, desat)
          qsat  = esat /this%aux_vars_in(icell)%pref
          dqsat = desat/this%aux_vars_in(icell)%pref

          gleaf = &
               this%aux_vars_in(icell)%gs * this%aux_vars_in(icell)%gbv / &
               (this%aux_vars_in(icell)%gs + this%aux_vars_in(icell)%gbv )

          gleaf_et = &
               gleaf                       * this%aux_vars_in(icell)%fdry + &
               this%aux_vars_in(icell)%gbv * this%aux_vars_in(icell)%fwet

          value = this%aux_vars_in(icell)%cp/this%dtime + &
               2.d0 * this%aux_vars_in(icell)%cpair * this%aux_vars_in(icell)%gbh + &
               lambda * dqsat * gleaf_et
       else
          value = 1.d0
       end if

       call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
    end do

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  end subroutine CLeafTempComputeOperatorsDiag

  !------------------------------------------------------------------------

  subroutine CLeafTempComputeOperatorsOffDiag(this, A, B, &
       itype_of_other_goveq, rank_of_other_goveq, ierr)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ConnectionSetType         , only : connection_set_type
    use ConditionType             , only : condition_type
    use MultiPhysicsProbConstants , only : HVAP, MM_H2O
    use MultiPhysicsProbConstants , only : GE_CANOPY_AIR_TEMP
    use MultiPhysicsProbConstants , only : GE_CANOPY_AIR_VAPOR
    use WaterVaporMod             , only : SatVap
    !
    implicit none
    !
    ! !ARGUMENTS
    class(goveqn_cleaf_temp_type) :: this
    Mat                           :: A
    Mat                           :: B
    PetscInt                      :: itype_of_other_goveq
    PetscInt                      :: rank_of_other_goveq
    PetscErrorCode                :: ierr
    !
    ! !LOCAL VARIABLES
    PetscInt :: icell, ileaf, row, col
    PetscReal :: value, gleaf, gleaf_et, qsat, dqsat, esat, desat, lambda

    lambda = HVAP * MM_H2O
    select case (itype_of_other_goveq)

    case (GE_CANOPY_AIR_TEMP)
       do icell = 1, this%mesh%ncells_local
          row = icell-1; col = icell-1
          col = this%Leaf2CAir(icell) - 1
          if (this%aux_vars_in(icell)%dpai > 0.d0) then
             value = -2.d0 * this%aux_vars_in(icell)%cpair * this%aux_vars_in(icell)%gbh
             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
          end if
       end do

    case (GE_CANOPY_AIR_VAPOR)
       do icell = 1, this%mesh%ncells_local
          row = icell-1; col = icell-1
          col = this%Leaf2CAir(icell) - 1
          if (this%aux_vars_in(icell)%dpai > 0.d0) then
             call SatVap(this%aux_vars_in(icell)%temperature, esat, desat)
             qsat  = esat /this%aux_vars_in(icell)%pref
             dqsat = desat/this%aux_vars_in(icell)%pref

             gleaf = &
                  this%aux_vars_in(icell)%gs * this%aux_vars_in(icell)%gbv / &
                  (this%aux_vars_in(icell)%gs + this%aux_vars_in(icell)%gbv )

             gleaf_et = &
                  gleaf                       * this%aux_vars_in(icell)%fdry + &
                  this%aux_vars_in(icell)%gbv * this%aux_vars_in(icell)%fwet

             value = -lambda * gleaf_et

             call MatSetValuesLocal(B, 1, row, 1, col, value, ADD_VALUES, ierr); CHKERRQ(ierr)
          end if
       end do

    end select

    call MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(  B, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

  end subroutine CLeafTempComputeOperatorsOffDiag

#endif

end module GoveqnCanopyLeafTemperatureType
