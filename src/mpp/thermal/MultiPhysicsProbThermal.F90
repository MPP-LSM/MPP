
module MultiPhysicsProbThermal

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Object for thermal model
  !-----------------------------------------------------------------------

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                         , only : iulog
  use mpp_abortutils                     , only : endrun
  use mpp_shr_log_mod                    , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbBaseType           , only : multiphysicsprob_base_type
  use SystemOfEquationsThermalType       , only : sysofeqns_thermal_type
  use SystemOfEquationsBasePointerType   , only : sysofeqns_base_pointer_type
  use SystemOfEquationsBaseType          , only : sysofeqns_base_type
  use petscsys
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda

  implicit none
  private

  type, public, extends(multiphysicsprob_base_type) :: mpp_thermal_type
   contains
     procedure, public :: Init            => ThermalMPPInit
     procedure, public :: AllocateAuxVars => ThermalMPPAllocateAuxVars
     procedure, public :: SetupProblem    => ThermalMPPSetupProblem

  end type mpp_thermal_type

  public :: MPPThermalSetSoils

  type(mpp_thermal_type), public, target :: thermal_mpp

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermalMPPInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize the thermal MPP
    !
    use MultiPhysicsProbBaseType , only : MPPBaseInit
    use MultiPhysicsProbConstants, only : MPP_THERMAL_TBASED_KSP_CLM
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this
    !
    class(sysofeqns_thermal_type), pointer :: sysofeqns

    call MPPBaseInit(this)

    allocate(sysofeqns)
    call sysofeqns%Init()

    this%soe => sysofeqns

    allocate(this%soe_ptr)
    nullify(this%soe_ptr%ptr)

  end subroutine ThermalMPPInit

  !------------------------------------------------------------------------
  subroutine MPPThermalSetSoils(therm_mpp, begc, endc, filter_thermal, &
       lun_type, watsat, csol, tkmg, tkdry)
    !
    ! !DESCRIPTION:
    ! Sets soil thermal properties
    !
    use GoverningEquationBaseType
    use GoveqnThermalKSPTemperatureSoilType
    use ThermalKSPTemperatureSoilAuxType
    use mpp_varpar       , only : nlevgrnd, nlevsoi, nlevsno
    use ConditionType    , only : condition_type
    use ConnectionSetType, only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                             :: therm_mpp
    integer              , intent(in)                   :: begc,endc
    integer, pointer, intent(in)                        :: filter_thermal(:)
    PetscInt, pointer, intent(in)                       :: lun_type(:)
    PetscReal, pointer, intent(in)                      :: watsat(:,:)
    PetscReal, pointer, intent(in)                      :: csol(:,:)
    PetscReal, pointer, intent(in)                      :: tkmg(:,:)
    PetscReal, pointer, intent(in)                      :: tkdry(:,:)

    !
    ! !LOCAL VARIABLES:
    class (goveqn_thermal_ksp_temp_soil_type) , pointer :: goveq_soil
    class(goveqn_base_type)                   , pointer :: cur_goveq
    type (therm_ksp_temp_soil_auxvar_type)    , pointer :: aux_vars_in(:)
    type (therm_ksp_temp_soil_auxvar_type)    , pointer :: aux_vars_bc(:)
    type(condition_type)                      , pointer :: cur_cond
    type(connection_set_type)                 , pointer :: cur_conn_set
    PetscInt                                            :: j,c,g,l
    PetscInt                                            :: sum_conn
    PetscInt                                            :: icell
    PetscInt                                            :: iconn
    PetscInt                                            :: col_id
    PetscInt                                            :: first_active_col_id
    PetscBool                                           :: found

    found = PETSC_FALSE
    cur_goveq => therm_mpp%soe%goveqns
    do
       if (.not.associated(cur_goveq)) exit

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_soil_type)
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

    icell = 0
    do c = begc, endc

       ! Soil layers
      do j = 1, nlevgrnd

         icell = icell + 1

         if (filter_thermal(c) == 1) then
            col_id                       = c
            aux_vars_in(icell)%is_active = PETSC_TRUE
         else
            col_id                       = first_active_col_id 
            aux_vars_in(icell)%is_active = PETSC_FALSE
         endif

         if (j > nlevsoi) then
            aux_vars_in(icell)%is_soil_shallow = PETSC_FALSE
         else
            aux_vars_in(icell)%is_soil_shallow = PETSC_TRUE
         endif         

         aux_vars_in(icell)%itype                 = lun_type(col_id)
         aux_vars_in(icell)%por                   = watsat(col_id,j)
         aux_vars_in(icell)%therm_cond_minerals   = tkmg(col_id,j)
         aux_vars_in(icell)%therm_cond_dry        = tkdry(col_id,j)
         aux_vars_in(icell)%heat_cap_minerals_puv = csol(col_id,j)

      enddo
   enddo

   sum_conn = 0
   aux_vars_bc => goveq_soil%aux_vars_bc
   cur_cond    => goveq_soil%boundary_conditions%first
   do
      if (.not.associated(cur_cond)) exit
      cur_conn_set => cur_cond%conn_set

      do iconn = 1, cur_conn_set%num_connections
         sum_conn = sum_conn + 1
         icell    = cur_conn_set%conn(iconn)%GetIDDn()

         aux_vars_bc(sum_conn)%itype                 = aux_vars_in(icell)%itype
         aux_vars_bc(sum_conn)%por                   = aux_vars_in(icell)%por
         aux_vars_bc(sum_conn)%therm_cond_minerals   = aux_vars_in(icell)%therm_cond_minerals
         aux_vars_bc(sum_conn)%therm_cond_dry        = aux_vars_in(icell)%therm_cond_dry
         aux_vars_bc(sum_conn)%heat_cap_minerals_puv = aux_vars_in(icell)%heat_cap_minerals_puv
      end do

      cur_cond => cur_cond%next
   end do

  end subroutine MPPThermalSetSoils

  !------------------------------------------------------------------------
  subroutine ThermalMPPAllocateAuxVars(this)
    !
    ! !DESCRIPTION:
    ! Allocates auxvars for governing equations and system-of-governing-eqns
    !
    use SystemOfEquationsBaseType           , only : sysofeqns_base_type
    use GoverningEquationBaseType           , only : goveqn_base_type
    use GoveqnThermalKSPTemperatureSnowType , only : goveqn_thermal_ksp_temp_snow_type
    use GoveqnThermalKSPTemperatureSSWType  , only : goveqn_thermal_ksp_temp_ssw_type
    use GoveqnThermalKSPTemperatureSoilType , only : goveqn_thermal_ksp_temp_soil_type
    use MultiPhysicsProbConstants           , only : COND_BC
    use MultiPhysicsProbConstants           , only : COND_SS
    use MultiPhysicsProbConstants           , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                :: this
    !
    class(sysofeqns_base_type), pointer    :: base_soe
    class(sysofeqns_thermal_type), pointer :: soe
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

    base_soe => this%soe
    
    select type(base_soe)
    class is(sysofeqns_thermal_type)
       soe => base_soe
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
       class is (goveqn_thermal_ksp_temp_snow_type)
          call cur_goveq%AllocateAuxVars()

       class is (goveqn_thermal_ksp_temp_ssw_type)
          call cur_goveq%AllocateAuxVars()

       class is (goveqn_thermal_ksp_temp_soil_type)
          call cur_goveq%AllocateAuxVars()

       end select

       call cur_goveq%GetNCellsInCondsExcptCondItype(COND_BC, &
            COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)

       call cur_goveq%GetNCellsInCondsExcptCondItype(COND_SS, -9999, &
            num_ss, ncells_for_ss)

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

       call cur_goveq%GetNCellsInCondsExcptCondItype(COND_BC, &
            COND_DIRICHLET_FRM_OTR_GOVEQ, num_bc, ncells_for_bc)

       call cur_goveq%GetNCellsInCondsExcptCondItype(COND_SS, -9999, &
            num_ss, ncells_for_ss)

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

             soe%aux_vars_bc(iauxvar)%is_in        = PETSC_TRUE
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

             soe%aux_vars_ss(iauxvar)%is_in        = PETSC_TRUE
             soe%aux_vars_ss(iauxvar)%goveqn_id    = igoveqn
             soe%aux_vars_ss(iauxvar)%condition_id = icond
          enddo
       enddo

       select type(cur_goveq)
       class is (goveqn_thermal_ksp_temp_snow_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)

       class is (goveqn_thermal_ksp_temp_ssw_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)

       class is (goveqn_thermal_ksp_temp_soil_type)
          call cur_goveq%SetSOEAuxVarOffsets(num_bc, offsets_bc, num_ss, offsets_ss)

       end select

       if (num_bc > 0) deallocate(ncells_for_bc)
       if (num_ss > 0) deallocate(ncells_for_ss)

       deallocate(offsets_bc)
       deallocate(offsets_ss)

       cur_goveq => cur_goveq%next
    enddo

  end subroutine ThermalMPPAllocateAuxVars

  !------------------------------------------------------------------------
  subroutine ThermalMPPSetupProblem(this)
    !
    ! !DESCRIPTION:
    ! Sets up the thermal MPP
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : SOE_THERMAL_TBASED
    use MultiPhysicsProbBaseType  , only : MPPSetupProblem
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type) :: this

    call ThermalMPPUpdatCouplingBCConnections(this)

    call MPPSetupProblem(this)

    this%soe%itype = SOE_THERMAL_TBASED

  end subroutine ThermalMPPSetupProblem

  !------------------------------------------------------------------------
  subroutine ThermalMPPUpdatCouplingBCConnections(this)
    !
    ! !DESCRIPTION:
    ! For BCs used to coupling two governing equations, updates
    ! ID and distance of upwind grid cell.
    !
    use ConditionType             , only : condition_type
    use GoverningEquationBaseType , only : goveqn_base_type
    use MultiPhysicsProbConstants , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    use SystemOfEquationsBaseType , only : sysofeqns_base_type
    use ConnectionSetType         , only : connection_set_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(mpp_thermal_type)                 :: this
    !
    class(sysofeqns_base_type)    , pointer :: base_soe
    class(sysofeqns_thermal_type) , pointer :: soe
    class(goveqn_base_type)       , pointer :: cur_goveq_1
    class(goveqn_base_type)       , pointer :: cur_goveq_2
    type(condition_type)          , pointer :: cur_cond_1
    type(condition_type)          , pointer :: cur_cond_2
    type(connection_set_type)     , pointer :: cur_conn_set_1
    type(connection_set_type)     , pointer :: cur_conn_set_2
    PetscInt                                :: igoveqn
    PetscInt                                :: ieqn
    PetscInt                                :: iconn
    PetscInt                                :: ii, jj
    PetscInt                                :: bc_idx_1
    PetscInt                                :: bc_idx_2
    PetscInt                                :: bc_offset_1
    PetscBool                               :: bc_found

    base_soe => this%soe

    select type(base_soe)
    class is(sysofeqns_thermal_type)
       soe => base_soe
    class default
       write(iulog,*) 'Unsupported class type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

    do ii = 1, soe%ngoveqns

       cur_goveq_1 => soe%goveqns
       do igoveqn = 1, ii-1
          cur_goveq_1 => cur_goveq_1%next
       end do


       do jj = ii+1, soe%ngoveqns

          cur_goveq_2 => soe%goveqns
          do igoveqn = 1, jj-1
             cur_goveq_2 => cur_goveq_2%next
          end do

          bc_found    = PETSC_FALSE
          cur_cond_1 => cur_goveq_1%boundary_conditions%first
          do
             if (.not.associated(cur_cond_1)) exit

             ! Is this the appropriate BC?
             if (cur_cond_1%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
                do ieqn = 1, cur_cond_1%num_other_goveqs
                   if (cur_cond_1%rank_of_other_goveqs(ieqn) == jj ) then
                      bc_found = PETSC_TRUE
                      exit
                   endif
                enddo
             endif

             if (bc_found) exit
             cur_cond_1 => cur_cond_1%next
          enddo

          if (bc_found) then

             bc_found    = PETSC_FALSE
             cur_cond_2 => cur_goveq_2%boundary_conditions%first
             do
                if (.not.associated(cur_cond_2)) exit

                ! Is this the appropriate BC?
                if (cur_cond_2%itype == COND_DIRICHLET_FRM_OTR_GOVEQ) then
                   do ieqn = 1, cur_cond_2%num_other_goveqs
                      if (cur_cond_2%rank_of_other_goveqs(ieqn) == ii ) then
                         bc_found = PETSC_TRUE
                         exit
                      endif
                   enddo
                endif

                if (bc_found) exit
                cur_cond_2 => cur_cond_2%next
             enddo

             if (.not.bc_found) then
                write(iulog,*)'cur_goveq_1: ',trim(cur_goveq_1%name)
                write(iulog,*)'cur_goveq_2: ',trim(cur_goveq_2%name)
                write(iulog,*)'cur_goveq_1 has a coupling BC for cur_eqn_2, but not vice-versa.'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             cur_conn_set_1 => cur_cond_1%conn_set
             cur_conn_set_2 => cur_cond_2%conn_set

             if (cur_conn_set_1%num_connections /= cur_conn_set_2%num_connections) then
                write(iulog,*) 'Number of connections in two equations are not same.'
                call endrun(msg=errMsg(__FILE__, __LINE__))
             endif

             do iconn = 1, cur_conn_set_1%num_connections
                call cur_conn_set_2%conn(iconn)%SetIDUp(cur_conn_set_1%conn(iconn)%GetIDDn())
                call cur_conn_set_2%conn(iconn)%SetDistUp(cur_conn_set_1%conn(iconn)%GetDistDn())

                call cur_conn_set_1%conn(iconn)%SetIDUp(cur_conn_set_2%conn(iconn)%GetIDDn())
                call cur_conn_set_1%conn(iconn)%SetDistUp(cur_conn_set_2%conn(iconn)%GetDistDn())
             enddo

          endif

       enddo
    enddo

  end subroutine ThermalMPPUpdatCouplingBCConnections

  !------------------------------------------------------------------------
#endif

end module MultiPhysicsProbThermal
