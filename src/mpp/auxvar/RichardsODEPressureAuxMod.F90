module RichardsODEPressureAuxMod

#ifdef USE_PETSC_LIB

  use mpp_varctl                 , only : iulog
  use mpp_abortutils             , only : endrun
  use mpp_shr_log_mod            , only : errMsg => shr_log_errMsg
  use RichardsODEPressureAuxType , only : rich_ode_pres_auxvar_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include <petsc/finclude/petsc.h>

  public :: RichardsODEPressureAuxVarSetRValues
  public :: RichardsODEPressureAuxVarGetRValues
  public :: RichardsODEPressureAuxVarSetAbsPerm
  public :: RichardsODEPressureAuxVarSetRelPerm
  public :: RichardsODEPressureAuxVarSetPorosity
  public :: RichardsODEPressureAuxVarSetSatFunc

contains
  
  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAuxVarSetRValues(auxvars, var_type, nauxvar, &
       auxvar_ids, data_1d)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(rich_ode_pres_auxvar_type) , pointer    :: auxvars(:)
    PetscInt                        , intent(in) :: var_type
    PetscInt                        , intent(in) :: nauxvar
    PetscInt                        , pointer    :: auxvar_ids(:)
    PetscReal                       , pointer    :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                     :: iauxvar

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'ERROR: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_type)
    case (VAR_TEMPERATURE)
       do iauxvar = 1, size(data_1d)
          auxvars(auxvar_ids(iauxvar))%temperature = data_1d(iauxvar)
       enddo


    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichardsODEPressureAuxVarSetRValues

  
  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAuxVarGetRValues(auxvars, var_type, nauxvar, &
       auxvar_ids, data_1d)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(rich_ode_pres_auxvar_type) , pointer    :: auxvars(:)
    PetscInt                        , intent(in) :: var_type
    PetscInt                        , intent(in) :: nauxvar
    PetscInt                        , pointer    :: auxvar_ids(:)
    PetscReal                       , pointer    :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                     :: iauxvar

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'ERROR: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_type)
    case (VAR_PRESSURE)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(auxvar_ids(iauxvar))%pressure
       enddo

    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichardsODEPressureAuxVarGetRValues

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAuxVarSetAbsPerm(auxvars_in, auxvars_bc, auxvars_ss, &
              boundary_cond, source_sinks, perm_x, perm_y, perm_z)
    !
    ! !DESCRIPTION:
    ! Set soil permeability values
    !
    ! !USES:
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type, condition_list_type
    use ConnectionSetType             , only : connection_set_type
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_in(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_bc(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_ss(:)
    type(condition_list_type)                                       :: boundary_cond, source_sinks
    PetscReal                                 , pointer, intent(in) :: perm_x(:)
    PetscReal                                 , pointer, intent(in) :: perm_y(:)
    PetscReal                                 , pointer, intent(in) :: perm_z(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                                        :: icell
    PetscInt                                                        :: ii
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    ! Set soil properties for internal auxvars
    do icell = 1, size(perm_x)
       auxvars_in(icell)%perm(1) = perm_x(icell)
       auxvars_in(icell)%perm(2) = perm_y(icell)
       auxvars_in(icell)%perm(3) = perm_z(icell)
    enddo

    ! Set soil properties for boundary-condition auxvars
    sum_conn = 0
    cur_cond => boundary_cond%first

    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()
             auxvars_bc(sum_conn)%perm(:) = auxvars_in(ghosted_id)%perm(:)
          enddo

       else
          cur_conn_set => cur_cond%conn_set

          ghosted_id = 1
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             auxvars_bc(sum_conn)%perm(:) = auxvars_in(ghosted_id)%perm(:)
          enddo

       endif
       cur_cond => cur_cond%next
    enddo

    ! Set soil properties for source-sink auxvars
    sum_conn = 0
    cur_cond => source_sinks%first

    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          auxvars_ss(sum_conn)%perm (:) = auxvars_in(ghosted_id)%perm(:)

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureAuxVarSetAbsPerm

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAuxVarSetRelPerm(auxvars_in, auxvars_bc, auxvars_ss, &
              boundary_cond, source_sinks, relperm_type, param_1, param_2)
    !
    ! !DESCRIPTION:
    ! Set soil permeability values
    !
    ! !USES:
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type, condition_list_type
    use ConnectionSetType             , only : connection_set_type
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_in(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_bc(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_ss(:)
    type(condition_list_type)                                       :: boundary_cond, source_sinks
    PetscInt                                  , pointer, intent(in) :: relperm_type(:)
    PetscReal                                 , pointer, intent(in) :: param_1(:)
    PetscReal                                 , pointer, intent(in) :: param_2(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                                        :: icell
    PetscInt                                                        :: igov
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    ! Set soil properties for internal auxvars
    do icell = 1, size(param_1)

       call auxvars_in(icell)%satParams%SetRelPerm( &
            relperm_type(icell), param_1(icell), param_2(icell))
    end do

    ! Set soil properties for boundary-condition auxvars
    sum_conn = 0
    cur_cond => boundary_cond%first

    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             call auxvars_bc(sum_conn)%satParams%Copy(auxvars_in(ghosted_id)%satParams)

          enddo
       else
          cur_conn_set => cur_cond%conn_set

          ghosted_id = 1
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             call auxvars_bc(sum_conn)%satParams%Copy(auxvars_in(ghosted_id)%satParams)
          enddo

       endif
       cur_cond => cur_cond%next
    enddo

    ! Set satParams for source-sink auxvars
    sum_conn = 0
    cur_cond => source_sinks%first

    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()
          call auxvars_ss(sum_conn)%satParams%Copy(auxvars_in(ghosted_id)%satParams)

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureAuxVarSetRelPerm

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAuxVarSetPorosity(auxvars_in, auxvars_bc, auxvars_ss, &
              boundary_cond, source_sinks, por)
    !
    ! !DESCRIPTION:
    ! Set soil porosity value
    !
    ! !USES:
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type, condition_list_type
    use ConnectionSetType             , only : connection_set_type
    use PorosityFunctionMod           , only : PorosityFunctionSetConstantModel
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_in(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_bc(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_ss(:)
    type(condition_list_type)                                       :: boundary_cond, source_sinks
    PetscReal                                 , pointer, intent(in) :: por(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                                        :: icell
    PetscInt                                                        :: ii
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    ! Set soil properties for internal auxvars
    do icell = 1, size(por)
       auxvars_in(icell)%por = por(icell)
       call PorosityFunctionSetConstantModel(auxvars_in(icell)%porParams, &
            auxvars_in(icell)%por)
    enddo

    ! Set soil properties for boundary-condition auxvars
    sum_conn = 0
    cur_cond => boundary_cond%first

    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then

          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             auxvars_bc(sum_conn)%por       = auxvars_in(ghosted_id)%por
             auxvars_bc(sum_conn)%porParams = auxvars_in(ghosted_id)%porParams

          enddo
       else
          cur_conn_set => cur_cond%conn_set

          ghosted_id = 1
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             auxvars_bc(sum_conn)%por       = auxvars_in(ghosted_id)%por
             auxvars_bc(sum_conn)%porParams = auxvars_in(ghosted_id)%porParams

          enddo
       end if
       cur_cond => cur_cond%next
    enddo

    ! Set soil properties for source-sink auxvars
    sum_conn = 0
    cur_cond => source_sinks%first

    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          auxvars_ss(sum_conn)%por       = auxvars_in(ghosted_id)%por
          auxvars_ss(sum_conn)%porParams = auxvars_in(ghosted_id)%porParams

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureAuxVarSetPorosity

  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAuxVarSetSatFunc(auxvars_in, auxvars_bc, auxvars_ss, &
              boundary_cond, source_sinks, satfunc_type, alpha, lambda, sat_res)
    !
    ! !DESCRIPTION:
    ! Set saturation function
    !
    ! !USES:
    use RichardsODEPressureAuxType    , only : rich_ode_pres_auxvar_type
    use ConditionType                 , only : condition_type, condition_list_type
    use ConnectionSetType             , only : connection_set_type
    use MultiPhysicsProbConstants     , only : COND_DIRICHLET_FRM_OTR_GOVEQ
    !
    implicit none
    !
    !
    ! !ARGUMENTS
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_in(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_bc(:)
    type (rich_ode_pres_auxvar_type)          , pointer             :: auxvars_ss(:)
    type(condition_list_type)                                       :: boundary_cond, source_sinks
    PetscInt                                  , pointer, intent(in) :: satfunc_type(:)
    PetscReal                                 , pointer, intent(in) :: alpha(:)
    PetscReal                                 , pointer, intent(in) :: lambda(:)
    PetscReal                                 , pointer, intent(in) :: sat_res(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                                        :: icell
    PetscInt                                                        :: ii
    type(condition_type)                      , pointer             :: cur_cond
    type(connection_set_type)                 , pointer             :: cur_conn_set
    PetscInt                                                        :: ghosted_id
    PetscInt                                                        :: sum_conn
    PetscInt                                                        :: iconn

    ! Set soil properties for internal auxvars
    do icell = 1, size(alpha)

       call auxvars_in(icell)%satParams%SetSatFunc( &
            satfunc_type(icell), alpha(icell), lambda(icell), sat_res(icell))
    end do

    ! Set soil properties for boundary-condition auxvars
    sum_conn = 0
    cur_cond => boundary_cond%first

    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= COND_DIRICHLET_FRM_OTR_GOVEQ) then
          cur_conn_set => cur_cond%conn_set

          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1
             ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

             call auxvars_bc(sum_conn)%satParams%Copy(auxvars_in(ghosted_id)%satParams)

          enddo
       else
          cur_conn_set => cur_cond%conn_set

          ghosted_id = 1
          do iconn = 1, cur_conn_set%num_connections
             sum_conn = sum_conn + 1

             call auxvars_bc(sum_conn)%satParams%Copy(auxvars_in(ghosted_id)%satParams)

          enddo
       endif
       cur_cond => cur_cond%next
    enddo

    ! Set soil properties for source-sink auxvars
    sum_conn = 0
    cur_cond => source_sinks%first

    do
       if (.not.associated(cur_cond)) exit
       cur_conn_set => cur_cond%conn_set

       do iconn = 1, cur_conn_set%num_connections
          sum_conn = sum_conn + 1
          ghosted_id = cur_conn_set%conn(iconn)%GetIDDn()

          call auxvars_ss(sum_conn)%satParams%Copy(auxvars_in(ghosted_id)%satParams)

       enddo
       cur_cond => cur_cond%next
    enddo

  end subroutine RichardsODEPressureAuxVarSetSatFunc

#endif

end module RichardsODEPressureAuxMod
