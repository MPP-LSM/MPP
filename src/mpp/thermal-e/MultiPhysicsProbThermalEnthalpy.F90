
module MultiPhysicsProbThermalEnthalpy

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for thermal model
  !-----------------------------------------------------------------------

  ! !USES:
  use mpp_varctl                           , only : iulog
  use mpp_abortutils                       , only : endrun
  use mpp_shr_log_mod                      , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbBaseType             , only : multiphysicsprob_base_type
  use SystemOfEquationsThermalEnthalpyType , only : sysofeqns_thermal_enthalpy_type
  use SystemOfEquationsBasePointerType     , only : sysofeqns_base_pointer_type

  implicit none
  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscts.h"
#include "finclude/petscts.h90"
#include "finclude/petscsnes.h"
#include "finclude/petscsnes.h90"
#include "finclude/petscdm.h"
#include "finclude/petscdm.h90"
#include "finclude/petscdmda.h"
#include "finclude/petscdmda.h90"
#include "finclude/petscviewer.h"

  type, public, extends(multiphysicsprob_base_type) :: mpp_thermal_type
     class(sysofeqns_thermal_enthalpy_type),pointer :: sysofeqns
     type(sysofeqns_base_pointer_type), pointer     :: sysofeqns_ptr
   contains
     procedure, public :: Init                        => ThermalEnthalpyMPPInit
     procedure, public :: AddGovEqn                   => ThermalEnthalpyMPPAddGovEqn
     procedure, public :: GovEqnAddCondition          => ThermalEnthalpyMPPGovEqnAddCondition
     procedure, public :: SetMeshesOfGoveqns          => ThermalEnthalpyMPPSetMeshesOfGoveqns
     procedure, public :: GovEqnAddCouplingCondition  => ThermalEnthalpyMPPGovEqnAddCouplingCondition
     procedure, public :: AllocateAuxVars             => ThermalEnthalpyMPPAllocateAuxVars
     procedure, public :: GovEqnSetCouplingVars       => ThermalEnthalpyMPPGovEqnSetCouplingVars
     procedure, public :: SetupProblem                => ThermalEnthalpyMPPSetupProblem
     procedure, public :: GovEqnUpdateBCConnectionSet => ThermalEnthalpyMPPGovEqnUpdateBCConnectionSet
     procedure, public :: SetMPIRank                  => ThermalEnthalpyMPPSetMPIRank

  end type mpp_thermal_type

  public :: MPPThermalSetSoils

  type(mpp_thermal_type), public, target :: thermal_enthalpy_mpp

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize the thermal MPP
    !
    use MultiPhysicsProbBaseType , only : MPPBaseInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this

    call MPPBaseInit(this)

    allocate(this%sysofeqns)
    call this%sysofeqns%Init()

    allocate(this%sysofeqns_ptr)
    nullify(this%sysofeqns_ptr%ptr)

  end subroutine ThermalEnthalpyMPPInit

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPSetMPIRank(this, rank)
    !
    ! !DESCRIPTION:
    ! Sets MPI rank
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt                :: rank

    if (associated(this%sysofeqns)) then
       this%sysofeqns%mpi_rank = rank
    endif

  end subroutine ThermalEnthalpyMPPSetMPIRank

  !------------------------------------------------------------------------
  subroutine MPPThermalSetSoils(therm_enth_mpp, begc, endc, filter_thermal, &
       lun_type, watsat, csol, tkmg, tkdry,                                 &
       hksat, bsw, sucsat, eff_porosity, residual_sat,                      &
       vsfm_satfunc_type, density_type)
    !
    ! !DESCRIPTION:
    ! Sets soil thermal properties
    !
    use GoverningEquationBaseType     , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType , only : goveqn_thermal_enthalpy_soil_type
    use ThermalEnthalpySoilAuxType    , only : therm_enthalpy_soil_auxvar_type
    use mpp_varpar                    , only : nlevgrnd, nlevsoi, nlevsno
    use SaturationFunction            , only : SatFunc_Set_BC
    use SaturationFunction            , only : SatFunc_Set_SBC_bz2
    use SaturationFunction            , only : SatFunc_Set_SBC_bz3
    use SaturationFunction            , only : SatFunc_Set_VG
    use PorosityFunctionMod           , only : PorosityFunctionSetConstantModel
    use mpp_varcon                    , only : grav, denh2o, tkwat
    use ConditionType                 , only : condition_type
    use ConnectionSetType             , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                             :: therm_enth_mpp
    integer              , intent(in)                   :: begc,endc
    integer, intent(in)                                 :: filter_thermal(:)
    PetscInt, pointer, intent(in)                       :: lun_type(:)
    PetscReal, pointer, intent(in)                      :: watsat(:,:)
    PetscReal, pointer, intent(in)                      :: csol(:,:)
    PetscReal, pointer, intent(in)                      :: tkmg(:,:)
    PetscReal, pointer, intent(in)                      :: tkdry(:,:)
    PetscReal, intent(in), pointer                      :: hksat(:,:)
    PetscReal, intent(in), pointer                      :: bsw(:,:)
    PetscReal, intent(in), pointer                      :: sucsat(:,:)
    PetscReal, intent(in), pointer                      :: eff_porosity(:,:)
    PetscReal, intent(in), pointer                      :: residual_sat(:,:)
    character(len=32), intent(in)                       :: vsfm_satfunc_type
    PetscInt                                            :: density_type

    !
    ! !LOCAL VARIABLES:
    class (goveqn_thermal_enthalpy_soil_type) , pointer :: goveq_soil
    class(sysofeqns_thermal_enthalpy_type)    , pointer :: therm_soe
    class(goveqn_base_type)                   , pointer :: cur_goveq
    type (therm_enthalpy_soil_auxvar_type)    , pointer :: aux_vars_in(:)
    type (therm_enthalpy_soil_auxvar_type)    , pointer :: aux_vars_bc(:)
    type (therm_enthalpy_soil_auxvar_type)    , pointer :: aux_vars_ss(:)
    type(condition_type),pointer                        :: cur_cond
    type(connection_set_type), pointer                  :: cur_conn_set
    PetscInt                                            :: j,c,g,l
    PetscInt                                            :: icell
    PetscInt                                            :: ghosted_id
    PetscInt                                            :: iconn
    PetscInt                                            :: sum_conn
    PetscInt                                            :: col_id
    PetscInt                                            :: first_active_col_id
    PetscBool                                           :: found
    PetscReal                                           :: perm
    PetscReal                                           :: por
    PetscReal                                           :: sat_res
    PetscReal                                           :: alpha
    PetscReal                                           :: lambda
    PetscReal, parameter                                :: vish2o = 0.001002d0    ! [N s/m^2] @ 20 degC

    therm_soe => therm_enth_mpp%sysofeqns

    found = PETSC_FALSE
    cur_goveq => therm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_enthalpy_soil_type)
          goveq_soil => cur_goveq
          found = PETSC_TRUE
       end select

       cur_goveq => cur_goveq%next
    enddo

    if (.not.found) then
       write(iulog,*)'Thermal governing equation for soil NOT found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    aux_vars_in => goveq_soil%aux_vars_in

    first_active_col_id = -1
    do c = begc, endc

       if (filter_thermal(c) == 1) then
          if (first_active_col_id == -1) then
             first_active_col_id = c
             exit
          endif
       endif
    enddo

    if (first_active_col_id == -1) then
       write(iulog,*)'No active soil column found'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Set thermal properties

    call goveq_soil%SetDensityType(density_type)

    icell = 0
    do c = begc, endc

       ! Soil layers
      do j = 1, nlevgrnd

         icell = icell + 1

         if (filter_thermal(c) == 1 ) then
            ! Columns on which thermal model is active
            col_id = c
         else
            col_id = first_active_col_id
         endif

         ! perm = hydraulic-conductivity * viscosity / ( density * gravity )
         ![m^2]  [mm/sec]   * [Ns/m^2] /([kg/m^3] * [m/s^2]) * [m/mm]
         perm  = hksat(col_id,j) * vish2o   /(denh2o   * grav   ) * 0.001d0

         ! alpha [1/Pa]; while sucsat [mm of H20]
         ! [Pa] = [mm of H20] * 0.001 [m/mm] * 1000 [kg/m^3] * 9.81 [m/sec^2]
         alpha = 1.d0/(sucsat(col_id,j)*grav)

         ! lambda = 1/bsw
         lambda = 1.d0/bsw(col_id,j)

         sat_res = residual_sat(col_id,j)
         por     = watsat(col_id,j)

         !aux_vars_in(icell)%perm(1:3) = perm
         aux_vars_in(icell)%por       = por

         call PorosityFunctionSetConstantModel(aux_vars_in(icell)%porParams, &
              aux_vars_in(icell)%por)

         if (vsfm_satfunc_type == 'brooks_corey') then
            call SatFunc_Set_BC(aux_vars_in(icell)%satParams,      &
                 sat_res,                               &
                 alpha,                                 &
                 lambda)
         elseif (vsfm_satfunc_type == 'smooth_brooks_corey_bz2') then
            call SatFunc_Set_SBC_bz2(aux_vars_in(icell)%satParams, &
                 sat_res,                          &
                 alpha,                            &
                 lambda,                           &
                 -0.9d0/alpha)
         elseif (vsfm_satfunc_type == 'smooth_brooks_corey_bz3') then
            call SatFunc_Set_SBC_bz3(aux_vars_in(icell)%satParams, &
                 sat_res,                          &
                 alpha,                            &
                 lambda,                           &
                 -0.9d0/alpha)
         elseif (vsfm_satfunc_type == 'van_genuchten') then
            call SatFunc_Set_VG(aux_vars_in(icell)%satParams,      &
                 sat_res,                               &
                 alpha,                                 &
                 lambda)
         else
            call endrun(msg='ERROR:: Unknown vsfm_satfunc_type = '//vsfm_satfunc_type//&
                 errMsg(__FILE__, __LINE__))
         endif

         !aux_vars_in(icell)%por                   = watsat(col_id,j)
         !aux_vars_in(icell)%therm_cond_minerals   = tkmg(col_id,j)
         !aux_vars_in(icell)%therm_cond_dry        = tkdry(col_id,j)
         !aux_vars_in(icell)%heat_cap_minerals_puv = csol(col_id,j)

         aux_vars_in(icell)%therm_alpha     = 0.45d0
         aux_vars_in(icell)%therm_cond_wet  = 1.3d0 !tkwat
         aux_vars_in(icell)%therm_cond_dry  = tkdry(col_id,j)
         aux_vars_in(icell)%heat_cap_soil   = csol(col_id,j)
         aux_vars_in(icell)%den_soil        = 2700.d0

      enddo
   enddo

   ! Set soil properties for boundary-condition auxvars
   sum_conn = 0
   aux_vars_bc => goveq_soil%aux_vars_bc
   cur_cond    => goveq_soil%boundary_conditions%first
   
   do
      if (.not.associated(cur_cond)) exit
      cur_conn_set => cur_cond%conn_set
      
      do iconn = 1, cur_conn_set%num_connections
         sum_conn = sum_conn + 1
         ghosted_id = cur_conn_set%id_dn(iconn)

         !aux_vars_bc(sum_conn)%perm(:)        = aux_vars_in(ghosted_id)%perm(:)
         aux_vars_bc(sum_conn)%por            = aux_vars_in(ghosted_id)%por
         aux_vars_bc(sum_conn)%satParams      = aux_vars_in(ghosted_id)%satParams
         aux_vars_bc(sum_conn)%porParams      = aux_vars_in(ghosted_id)%porParams

         aux_vars_bc(sum_conn)%therm_alpha    = aux_vars_in(ghosted_id)%therm_alpha
         aux_vars_bc(sum_conn)%therm_cond_wet = aux_vars_in(ghosted_id)%therm_cond_wet
         aux_vars_bc(sum_conn)%therm_cond_dry = aux_vars_in(ghosted_id)%therm_cond_dry
         aux_vars_bc(sum_conn)%heat_cap_soil  = aux_vars_in(ghosted_id)%heat_cap_soil
         aux_vars_bc(sum_conn)%den_soil       = aux_vars_in(ghosted_id)%den_soil
         
         call aux_vars_bc(sum_conn)%satParams%Copy(aux_vars_in(ghosted_id)%satParams)
         
      enddo
      cur_cond => cur_cond%next
   enddo
   
   ! Set soil properties for source-sink auxvars
   sum_conn = 0
   
   aux_vars_ss => goveq_soil%aux_vars_ss
   cur_cond    => goveq_soil%source_sinks%first
   
   do
      if (.not.associated(cur_cond)) exit
      cur_conn_set => cur_cond%conn_set
      
      do iconn = 1, cur_conn_set%num_connections
         sum_conn = sum_conn + 1
         ghosted_id = cur_conn_set%id_dn(iconn)
         
         !aux_vars_ss(sum_conn)%perm(:)        = aux_vars_in(ghosted_id)%perm(:)
         aux_vars_ss(sum_conn)%por            = aux_vars_in(ghosted_id)%por
         aux_vars_ss(sum_conn)%satParams      = aux_vars_in(ghosted_id)%satParams
         aux_vars_ss(sum_conn)%porParams      = aux_vars_in(ghosted_id)%porParams

         aux_vars_ss(sum_conn)%therm_alpha    = aux_vars_in(ghosted_id)%therm_alpha
         aux_vars_ss(sum_conn)%therm_cond_wet = aux_vars_in(ghosted_id)%therm_cond_wet
         aux_vars_ss(sum_conn)%therm_cond_dry = aux_vars_in(ghosted_id)%therm_cond_dry
         aux_vars_ss(sum_conn)%heat_cap_soil  = aux_vars_in(ghosted_id)%heat_cap_soil
         aux_vars_ss(sum_conn)%den_soil       = aux_vars_in(ghosted_id)%den_soil
         
         !aux_vars_ss(sum_conn)%pressure_prev = 3.5355d3
         
      enddo
      cur_cond => cur_cond%next
   enddo

  end subroutine MPPThermalSetSoils

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPAddGovEqn(this, geq_type, name)
    !
    ! !DESCRIPTION:
    ! Adds a governing equation to the MPP
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt                :: geq_type
    character(len =*)       :: name

    call this%sysofeqns%AddGovEqn(geq_type, name)

  end subroutine ThermalEnthalpyMPPAddGovEqn

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPGovEqnAddCondition(this, igoveqn, ss_or_bc_type, name, unit, &
       cond_type, region_type, id_of_other_goveq)
    !
    ! !DESCRIPTION:
    ! Adds a boundary/source-sink condition to a governing equation
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt                          :: igoveqn
    PetscInt                          :: ss_or_bc_type
    character(len =*)                 :: name
    character(len =*)                 :: unit
    PetscInt                          :: cond_type
    PetscInt                          :: region_type
    PetscInt, optional                :: id_of_other_goveq
    !
    class(goveqn_base_type),pointer   :: cur_goveq
    class(goveqn_base_type),pointer   :: other_goveq
    PetscInt                          :: ii

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to add condition for governing equation ' // &
            'that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    if (.not.present(id_of_other_goveq)) then
       call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
            cond_type, region_type)
    else

       other_goveq => this%sysofeqns%goveqns
       do ii = 1,id_of_other_goveq-1
          other_goveq => other_goveq%next
       enddo

       call cur_goveq%AddCondition(ss_or_bc_type, name, unit, &
            cond_type, region_type, id_of_other_goveq, other_goveq%id )
    endif

  end subroutine ThermalEnthalpyMPPGovEqnAddCondition

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPSetMeshesOfGoveqns(this)
    !
    ! !DESCRIPTION:
    ! Set association of governing equations and meshes
    !
    use GoverningEquationBaseType, only : goveqn_base_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this

    call this%sysofeqns%SetMeshesOfGoveqns(this%meshes, this%nmesh)

  end subroutine ThermalEnthalpyMPPSetMeshesOfGoveqns

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPGovEqnAddCouplingCondition(this, ieqn_1, ieqn_2, &
       iregion_1, iregion_2)
    !
    ! !DESCRIPTION:
    ! Adds a boundary condition to couple ieqn_1 and ieqn_2
    !
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt                :: ieqn_1
    PetscInt                :: ieqn_2
    PetscInt                :: iregion_1
    PetscInt                :: iregion_2
    !
    character(len=256)        :: name

    write(name,*) ieqn_2
    name = 'BC_for_coupling_with_equation_' // trim(adjustl(name))
    call this%GovEqnAddCondition(ieqn_1, COND_BC, &
         name, '[K]', COND_DIRICHLET_FRM_OTR_GOVEQ, iregion_1, ieqn_2)

    write(name,*) ieqn_1
    name = 'BC_for_coupling_with_equation_' // trim(adjustl(name))
    call this%GovEqnAddCondition(ieqn_2, COND_BC, &
         name, '[K]', COND_DIRICHLET_FRM_OTR_GOVEQ, iregion_2, ieqn_1)

  end subroutine ThermalEnthalpyMPPGovEqnAddCouplingCondition

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates auxvars for governing equations and system-of-governing-eqns
    !
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnThermalEnthalpySoilType       , only : goveqn_thermal_enthalpy_soil_type
    use MultiPhysicsProbConstants           , only : COND_BC
    use MultiPhysicsProbConstants           , only : COND_SS
    use MultiPhysicsProbConstants           , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                :: this
    !
    class(sysofeqns_base_type), pointer    :: soe_base
    class(sysofeqns_thermal_enthalpy_type), pointer :: soe
    class(goveqn_base_type), pointer       :: cur_goveq
    PetscInt                               :: igoveqn
    PetscInt                               :: num_bc
    PetscInt                               :: num_ss
    PetscInt                               :: icond
    PetscInt                               :: iauxvar
    PetscInt                               :: iauxvar_beg, iauxvar_end
    PetscInt                               :: iauxvar_beg_bc, iauxvar_end_bc
    PetscInt                               :: iauxvar_beg_ss, iauxvar_end_ss
    PetscInt                               :: count_bc, count_ss
    PetscInt                               :: offset_bc, offset_ss
    PetscInt, pointer                      :: ncells_for_bc(:)
    PetscInt, pointer                      :: ncells_for_ss(:)
    PetscInt, pointer                      :: offsets_bc(:)
    PetscInt, pointer                      :: offsets_ss(:)

    soe_base => this%sysofeqns

    select type(soe_base)
    class is(sysofeqns_thermal_enthalpy_type)
       soe => this%sysofeqns
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    !
    ! Pass-1: Determine total number of BCs (excluding BCs
    !         needed for coupling various governing equations)
    !         and SSs for all governing equations
    !
    igoveqn = 0
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_enthalpy_soil_type)
          call cur_goveq%AllocateAuxVars()
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       end select

       igoveqn = igoveqn + 1

       soe%num_auxvars_in = soe%num_auxvars_in + &
            cur_goveq%mesh%ncells_all

       do icond = 1, num_bc
          soe%num_auxvars_bc = soe%num_auxvars_bc + &
               ncells_for_bc(icond)
       enddo

       do icond = 1, num_ss
          soe%num_auxvars_ss = soe%num_auxvars_ss + &
               ncells_for_ss(icond)
       enddo

       if (num_bc > 0) deallocate(ncells_for_bc)
       if (num_ss > 0) deallocate(ncells_for_ss)

       cur_goveq => cur_goveq%next
    enddo

    ! Allocate memory
    allocate(soe%aux_vars_in           (soe%num_auxvars_in))
    allocate(soe%aux_vars_bc           (soe%num_auxvars_bc))
    allocate(soe%aux_vars_ss           (soe%num_auxvars_ss))

    allocate(soe%soe_auxvars_bc_offset (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_ss_offset (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_bc_ncells (soe%num_auxvars_bc))
    allocate(soe%soe_auxvars_ss_ncells (soe%num_auxvars_ss))

    igoveqn        = 0
    iauxvar_beg    = 0
    iauxvar_end    = 0
    iauxvar_beg_bc = 0
    iauxvar_end_bc = 0
    iauxvar_beg_ss = 0
    iauxvar_end_ss = 0
    count_bc       = 0
    count_ss       = 0
    offset_bc      = 0
    offset_ss      = 0

    !
    ! Pass-2: Set values for auxvars of SoE
    !
    cur_goveq => soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_enthalpy_soil_type)
          call cur_goveq%GetNumCellsInConditions(COND_BC, &
               COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)
          call cur_goveq%GetNumCellsInConditions(COND_SS, -9999, &
               num_ss, ncells_for_ss)

       end select

       igoveqn = igoveqn + 1

       iauxvar_beg = iauxvar_end + 1
       iauxvar_end = iauxvar_end + cur_goveq%mesh%ncells_all

       do iauxvar = iauxvar_beg, iauxvar_end
          call soe%aux_vars_in(iauxvar)%Init()

          soe%aux_vars_in(iauxvar)%is_in     = PETSC_TRUE
          soe%aux_vars_in(iauxvar)%goveqn_id = igoveqn
       enddo

       allocate(offsets_bc(num_bc))
       allocate(offsets_ss(num_ss))

       do icond = 1, num_bc
          count_bc = count_bc + 1

          soe%soe_auxvars_bc_offset(count_bc) = offset_bc
          offsets_bc(icond)                   = offset_bc
          soe%soe_auxvars_bc_ncells(count_bc) = ncells_for_bc(icond)
          offset_bc                           = offset_bc + ncells_for_bc(icond)

          iauxvar_beg_bc = iauxvar_end_bc + 1
          iauxvar_end_bc = iauxvar_end_bc + ncells_for_bc(icond)

          do iauxvar = iauxvar_beg_bc, iauxvar_end_bc
             call soe%aux_vars_bc(iauxvar)%Init()

             soe%aux_vars_bc(iauxvar)%is_bc        = PETSC_TRUE
             soe%aux_vars_bc(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_bc(iauxvar)%condition_id = icond
          enddo
       enddo

       do icond = 1, num_ss
          count_ss = count_ss + 1
          soe%soe_auxvars_ss_offset(count_ss) = offset_ss
          offsets_ss(icond)                   = offset_ss
          soe%soe_auxvars_ss_ncells(count_ss) = ncells_for_ss(icond)
          offset_ss                           = offset_ss + ncells_for_ss(icond)

          iauxvar_beg_ss = iauxvar_end_ss + 1
          iauxvar_end_ss = iauxvar_end_ss + ncells_for_ss(icond)

          do iauxvar = iauxvar_beg_ss, iauxvar_end_ss
             call soe%aux_vars_ss(iauxvar)%Init()

             soe%aux_vars_ss(iauxvar)%is_ss        = PETSC_TRUE
             soe%aux_vars_ss(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_ss(iauxvar)%condition_id = icond
          enddo
       enddo

       select type(cur_goveq)
       class is (goveqn_thermal_enthalpy_soil_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)

       end select

       if (num_bc > 0) deallocate(ncells_for_bc)
       if (num_ss > 0) deallocate(ncells_for_ss)

       deallocate(offsets_bc)
       deallocate(offsets_ss)

       cur_goveq => cur_goveq%next
    enddo

  end subroutine ThermalEnthalpyMPPAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPGovEqnSetCouplingVars(this, igoveqn, nvars, &
       var_ids, goveqn_ids)
    !
    ! !DESCRIPTION:
    ! In order to couple the given governing equation, add:
    ! - ids of variables needed, and
    ! - ids of governing equations from which variables are needed.
    ! needed for coupling
    ! 
    ! !USES:
    use ConditionType             , only : condition_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)           :: this
    PetscInt                          :: igoveqn
    PetscInt                          :: nvars
    PetscInt, pointer                 :: var_ids(:)
    PetscInt, pointer                 :: goveqn_ids(:)
    !
    class(goveqn_base_type) , pointer :: cur_goveq_1
    class(goveqn_base_type) , pointer :: cur_goveq_2
    type(condition_type)    , pointer :: cur_cond_1
    type(condition_type)    , pointer :: cur_cond_2
    PetscInt                          :: ii
    PetscInt                          :: ivar
    PetscInt                          :: bc_idx_1
    PetscInt                          :: bc_idx_2
    PetscInt                          :: bc_offset_1
    PetscBool                         :: bc_found

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to set coupling vars for governing ' // &
            'equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq_1 => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq_1 => cur_goveq_1%next
    end do

    call cur_goveq_1%AllocVarsFromOtherGEs(nvars)

    do ivar = 1, nvars

       if (goveqn_ids(ivar) > this%sysofeqns%ngoveqns) then
          write(iulog,*) 'Attempting to set coupling vars to a governing ' // &
               'equation that is not in the list'
          call endrun(msg=errMsg(__FILE__, __LINE__))
       endif

       bc_found    = PETSC_FALSE
       bc_idx_1    = 1
       bc_offset_1 = 0

       cur_cond_1 => cur_goveq_1%boundary_conditions%first
       do
          if (.not.associated(cur_cond_1)) exit

          ! Is this the appropriate BC?
          if (cur_cond_1%itype == COND_DIRICHLET_FRM_OTR_GOVEQ .and. &
               cur_cond_1%list_id_of_other_goveq == goveqn_ids(ivar) ) then
             bc_found = PETSC_TRUE
             exit
          endif

          bc_idx_1    = bc_idx_1    + 1
          bc_offset_1 = bc_offset_1 + cur_cond_1%conn_set%num_connections

          cur_cond_1 => cur_cond_1%next
       enddo

       if (.not.bc_found) then
          write(iulog,*)'For goveqn%name = ',trim(cur_goveq_1%name) // &
               ', no coupling boundary condition found to copule it with ' // &
               'equation_number = ', goveqn_ids(ivar)
       endif

       cur_goveq_2 => this%sysofeqns%goveqns
       do ii = 1, goveqn_ids(ivar)-1
          cur_goveq_2 => cur_goveq_2%next
       enddo

       bc_found    = PETSC_FALSE
       bc_idx_2    = 1

       cur_cond_2 => cur_goveq_2%boundary_conditions%first
       do
          if (.not.associated(cur_cond_2)) exit

          ! Is this the appropriate BC?
          if (cur_cond_2%itype == COND_DIRICHLET_FRM_OTR_GOVEQ .and. &
               cur_cond_2%list_id_of_other_goveq == igoveqn ) then
             bc_found = PETSC_TRUE
             exit
          endif

          bc_idx_2    = bc_idx_2    + 1

          cur_cond_2 => cur_cond_2%next
       enddo

       if (.not.bc_found) then
          write(iulog,*)'For goveqn%name = ',trim(cur_goveq_2%name) // &
               ', no coupling boundary condition found to copule it with ' // &
               'equation_number = ', bc_idx_2
       endif

       cur_goveq_1%var_ids_needed_from_other_goveqns (ivar) = var_ids(ivar)
       cur_goveq_1%ids_of_other_goveqns              (ivar) = goveqn_ids(ivar)
       cur_goveq_1%is_bc_auxvar_type                 (ivar) = PETSC_TRUE
       cur_goveq_1%bc_auxvar_offset                  (ivar) = bc_offset_1
       cur_goveq_1%bc_auxvar_ncells                  (ivar) = cur_cond_1%conn_set%num_connections
       cur_goveq_1%bc_auxvar_idx                     (ivar) = bc_idx_1
       cur_goveq_1%bc_auxvar_idx_of_other_goveqn     (ivar) = bc_idx_2

    enddo

  end subroutine ThermalEnthalpyMPPGovEqnSetCouplingVars

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPSetupProblem(this)
    !
    ! !DESCRIPTION:
    ! Sets up the thermal MPP
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : SOE_THERMAL_EBASED
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)           :: this

    call ThermalEnthalpyMPPKSPSetup(this)

    this%sysofeqns%solver_type = this%solver_type
    this%sysofeqns%itype       = SOE_THERMAL_EBASED

  end subroutine ThermalEnthalpyMPPSetupProblem

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPKSPSetup(therm_enth_mpp)
    !
    ! !DESCRIPTION:
    ! Sets the PETSc KSP for the thermal mpp
    !
    use GoverningEquationBaseType        , only : goveqn_base_type
    use SystemOfEquationsBasePointerType , only : SOEResidual
    use SystemOfEquationsBasePointerType , only : SOEJacobian
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                            :: therm_enth_mpp
    !
    ! !LOCAL VARIABLES:
    class(goveqn_base_type)                , pointer   :: cur_goveq
    class(sysofeqns_thermal_enthalpy_type) , pointer   :: therm_soe
    PetscInt                                           :: size
    PetscInt                                           :: igoveq
    PetscErrorCode                                     :: ierr
    DM                                     , pointer   :: dms(:)
    PetscReal                              , parameter :: atol    = PETSC_DEFAULT_REAL
    PetscReal                              , parameter :: rtol    = PETSC_DEFAULT_REAL
    PetscReal                              , parameter :: stol    = 1.d-10
    PetscInt                               , parameter :: max_it  = PETSC_DEFAULT_INTEGER
    PetscInt                               , parameter :: max_f   = PETSC_DEFAULT_INTEGER
    character(len=256)                                 :: name

    therm_soe => therm_enth_mpp%sysofeqns
    therm_enth_mpp%sysofeqns_ptr%ptr => therm_enth_mpp%sysofeqns

    ! Create PETSc DM for each governing equation

    allocate(dms(therm_soe%ngoveqns))

    igoveq = 0
    cur_goveq => therm_soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       igoveq = igoveq + 1
       size   = cur_goveq%mesh%ncells_local

       call DMDACreate1d(PETSC_COMM_SELF, &
            DM_BOUNDARY_NONE, size, 1, 0, &
            PETSC_NULL_INTEGER, dms(igoveq), ierr);
       CHKERRQ(ierr)

       write(name,*) igoveq
       name = 'fgoveq_' // trim(adjustl(name))
       call DMSetOptionsPrefix(dms(igoveq), name, ierr); CHKERRQ(ierr)

       write(name,*) igoveq
       name = 'goveq_' // trim(adjustl(name))
       call DMDASetFieldName(dms(igoveq), 0, name, ierr); CHKERRQ(ierr)

       call DMSetFromOptions(dms(igoveq), ierr); CHKERRQ(ierr)

       cur_goveq => cur_goveq%next
    enddo

    ! DM-Composite approach

    ! Create DMComposite: temperature
    call DMCompositeCreate(PETSC_COMM_SELF, therm_soe%dm, ierr); CHKERRQ(ierr)
    call DMSetOptionsPrefix(therm_soe%dm, "temperature_", ierr); CHKERRQ(ierr)

    ! Add DMs to DMComposite
    do igoveq = 1, therm_soe%ngoveqns
       call DMCompositeAddDM(therm_soe%dm, dms(igoveq), ierr); CHKERRQ(ierr)
    enddo

    ! Setup DM
    call DMSetUp(therm_soe%dm, ierr); CHKERRQ(ierr)

    ! Create matrix
    call DMCreateMatrix    (therm_soe%dm   , therm_soe%Amat, ierr); CHKERRQ(ierr)

    call MatSetOption      (therm_soe%Amat , MAT_NEW_NONZERO_LOCATION_ERR , &
         PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption      (therm_soe%Amat , MAT_NEW_NONZERO_ALLOCATION_ERR, &
         PETSC_FALSE, ierr); CHKERRQ(ierr)

    call MatSetFromOptions (therm_soe%Amat , ierr); CHKERRQ(ierr)

    call DMCreateMatrix     (therm_soe%dm , therm_soe%jac, ierr ); CHKERRQ(ierr)
    call MatSetOption       (therm_soe%jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)
    call MatSetOption       (therm_soe%jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHKERRQ(ierr)

    ! Create vectors
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln         , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%rhs          , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln_prev    , ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(therm_soe%dm, therm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! Initialize vectors
    call VecZeroEntries(therm_soe%soln          , ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%soln_prev     , ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%soln_prev_clm , ierr); CHKERRQ(ierr)

    call VecZeroEntries(therm_soe%soln         , ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%soln_prev    , ierr); CHKERRQ(ierr)
    call VecZeroEntries(therm_soe%soln_prev_clm, ierr); CHKERRQ(ierr)

    ! Create SNES
    call SNESCreate             (PETSC_COMM_SELF , therm_soe%snes, ierr); CHKERRQ(ierr)
    call SNESSetOptionsPrefix   (therm_soe%snes  , "temperature_", ierr); CHKERRQ(ierr)

    call SNESSetTolerances(therm_soe%snes, atol, rtol, stol, &
                           max_it, max_f, ierr); CHKERRQ(ierr)

    call SNESSetFunction(therm_soe%snes, PETSC_NULL_OBJECT, SOEResidual, &
         thermal_enthalpy_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)

    call SNESSetJacobian(therm_soe%snes, therm_soe%jac, therm_soe%jac,     &
         SOEJacobian, thermal_enthalpy_mpp%sysofeqns_ptr, ierr); CHKERRQ(ierr)

    call SNESSetFromOptions(therm_soe%snes, ierr); CHKERRQ(ierr)

    ! Get pointers to governing-equations
    call therm_soe%CreateVectorsForGovEqn()

    ! Cleanup
    do igoveq = 1, therm_soe%ngoveqns
       call DMDestroy(dms(igoveq), ierr); CHKERRQ(ierr)
    enddo
    deallocate(dms)

  end subroutine ThermalEnthalpyMPPKSPSetup

  !------------------------------------------------------------------------
  subroutine ThermalEnthalpyMPPGovEqnUpdateBCConnectionSet(this, igoveqn, icond, &
       var_type, nval, values)
    !
    ! !DESCRIPTION:
    ! For a boundary condition of a given governing equation, update distance
    ! for a downstream cell.
    !
    use ConditionType             , only : condition_type
    use ConnectionSetType         , only : connection_set_type
    use GoverningEquationBaseType, only : goveqn_base_type
    use MultiPhysicsProbConstants, only : VAR_DIST_DN
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    PetscInt :: igoveqn
    PetscInt :: icond
    PetscInt :: nval
    PetscInt :: var_type
    PetscReal, pointer :: values (:)
    !
    class(goveqn_base_type),pointer   :: cur_goveq
    type(condition_type)    , pointer :: cur_cond
    type(connection_set_type)     , pointer :: cur_conn_set
    PetscInt :: ii
    PetscInt :: iconn
    PetscInt :: bc_idx
    PetscBool :: bc_found

    if (igoveqn > this%sysofeqns%ngoveqns) then
       write(iulog,*) 'Attempting to access governing equation that is not in the list'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    cur_goveq => this%sysofeqns%goveqns
    do ii = 1, igoveqn-1
       cur_goveq => cur_goveq%next
    enddo

    bc_found = PETSC_FALSE
    cur_cond => cur_goveq%boundary_conditions%first
    do
       if (.not.associated(cur_cond)) exit

       bc_idx = bc_idx + 1
       if (bc_idx == icond) then
          bc_found = PETSC_TRUE

          cur_conn_set => cur_cond%conn_set
          if (nval /= cur_conn_set%num_connections) then
             write(iulog,*) 'Number of values to update connections ' // &
                  'do not match number of connections.'
             call endrun(msg=errMsg(__FILE__, __LINE__))
          endif

          do iconn = 1, cur_conn_set%num_connections

             select case(var_type)
             case (VAR_DIST_DN)
                cur_conn_set%dist_dn(iconn) = values(iconn)
             case default
                write(iulog,*) 'Unknown variable type'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             end select
          enddo

          exit

       end if

       cur_cond => cur_cond%next
    enddo

    if (.not.bc_found) then
       write(iulog,*) 'Failed to find icond = ',icond,' in the boundary condition list.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif


  end subroutine ThermalEnthalpyMPPGovEqnUpdateBCConnectionSet

#endif

end module MultiPhysicsProbThermalEnthalpy
