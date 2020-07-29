module shortwave_conditions

#include <petsc/finclude/petsc.h>

  use MultiPhysicsProbShortwave , only : mpp_shortwave_type
  use shortwave_global_vars
  use petscsys

  implicit none

  public :: add_conditions_to_goveqns

contains

  !------------------------------------------------------------------------
  subroutine add_conditions_to_goveqns(shortwave_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS
    !
    implicit none
    !
    type(mpp_shortwave_type) :: shortwave_mpp

    call add_non_coupling_conditions_to_goveqns(shortwave_mpp)

  end subroutine add_conditions_to_goveqns

  !------------------------------------------------------------------------
  subroutine add_non_coupling_conditions_to_goveqns(shortwave_mpp)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : COND_BC
    use MultiPhysicsProbConstants , only : COND_DIRICHLET
    use ConnectionSetType               , only : connection_set_type
    use MeshType                        , only : MeshCreateConnectionSet
    use MultiPhysicsProbConstants       , only : CONN_VERTICAL
    !
    ! !ARGUMENTS
    implicit none
    !
    type(mpp_shortwave_type) :: shortwave_mpp
    !
    PetscInt            :: kk , nconn
    PetscInt  , pointer :: id_up(:)
    PetscInt  , pointer :: id_dn(:)
    PetscReal , pointer :: dist_up(:)
    PetscReal , pointer :: dist_dn(:)
    PetscReal , pointer :: area(:)
    PetscInt  , pointer :: itype(:)
    PetscReal , pointer :: unit_vec(:,:)
    PetscReal                  :: dz_cair
    class(connection_set_type) , pointer :: conn_set


    nconn  = ncair

    allocate(id_up    (nconn   ))
    allocate(id_dn    (nconn   ))
    allocate(dist_up  (nconn   ))
    allocate(dist_dn  (nconn   ))
    allocate(area     (nconn   ))
    allocate(itype    (nconn   ))
    allocate(unit_vec (nconn,3 ))

    dz_cair = z_cair/nz_cair;

    do kk = 1, ncair
       id_up(kk)      = 0
       id_dn(kk)      = (nz_cair + 1)*kk
       dist_up(kk)    =  0.d0
       dist_dn(kk)    = dz_cair
       area(kk)       = 1.d0
       unit_vec(kk,1) = -1.d0
       unit_vec(kk,2) = 0.d0
       unit_vec(kk,3) = 0.d0
       itype(kk)      = CONN_VERTICAL
    end do

    allocate(conn_set)
    call MeshCreateConnectionSet(shortwave_mpp%meshes(shortwave_MESH), &
         nconn, id_up, id_dn, dist_up, dist_dn, area, itype, unit_vec, conn_set)

    call shortwave_mpp%soe%AddConditionInGovEqn(shortwave_GE, &
         ss_or_bc_type = COND_BC,   &
         name = 'Atmospheric forcing', &
         unit = 'K', &
         cond_type = COND_DIRICHLET, &
         conn_set = conn_set)

  end subroutine add_non_coupling_conditions_to_goveqns

end module shortwave_conditions
