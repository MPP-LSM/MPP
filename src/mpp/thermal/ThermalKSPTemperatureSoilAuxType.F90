module ThermalKSPTemperatureSoilAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !
  ! !USES:
  use ThermalKSPTemperatureBaseAuxType
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none

  private

  type, public, extends(therm_ksp_temp_base_auxvar_type)  :: therm_ksp_temp_soil_auxvar_type

     private

     PetscReal :: liq_areal_den         ! [kg/m^2] = h2osoi_liq
     PetscReal :: ice_areal_den         ! [kg/m^2] = h2osoi_ice
     PetscReal :: snow_water            ! snow water (mm H2O)
     PetscInt  :: num_snow_layer        ! number of snow layer
     PetscReal :: tuning_factor         !
     PetscReal :: dz                    !

     ! parameters
     PetscReal :: por                   ! [m^3/m^3]
     PetscReal :: therm_cond_minerals   ! [W/(m K)]
     PetscReal :: therm_cond_dry        ! [W/(m K)]
     PetscReal :: heat_cap_minerals_puv ! [J/(m3 K)]
     PetscBool :: is_soil_shallow
     PetscInt  :: itype

   contains
     procedure, public :: Init                              => ThermKSPTempSoilAuxVarInit
     procedure, public :: AuxVarCompute                     => ThermKSPTempSoilAuxVarCompute

     procedure, public :: SetLiqArealDensity                => ThermalKSPTemperatureSoilAuxSetLiqArealDensity
     procedure, public :: SetIceArealDensity                => ThermalKSPTemperatureSoilAuxSetIceArealDensity
     procedure, public :: SetSnowWater                      => ThermalKSPTemperatureSoilAuxSetSnowWater
     procedure, public :: SetNumberOfSnowLayers             => ThermalKSPTemperatureSoilAuxSetNumberOfSnowLayers
     procedure, public :: SetTuningFactor                   => ThermalKSPTemperatureSoilAuxSetTuningFactor
     procedure, public :: SetDepth                          => ThermalKSPTemperatureSoilAuxSetDepth
     procedure, public :: SetPorosity                       => ThermalKSPTemperatureSoilAuxSetPorosity
     procedure, public :: SetThermalCondMinerals            => ThermalKSPTemperatureSoilAuxSetThermalCondMinerals
     procedure, public :: SetThermalCondDrySoil             => ThermalKSPTemperatureSoilAuxSetThermalCondDrySoil
     procedure, public :: SetVolumetricHeatCapacityMinerals => ThermalKSPTemperatureSoilAuxSetVolumetricHeatCapacityMinerals
     procedure, public :: SetSoilShallow                    => ThermalKSPTemperatureSoilAuxSetSoilShallow
     procedure, public :: SetSoilDeep                       => ThermalKSPTemperatureSoilAuxSetSoilDeep
     procedure, public :: SetColumnType                     => ThermalKSPTemperatureSoilAuxSetColumnType

     procedure, public :: GetLiqArealDensity                => ThermalKSPTemperatureSoilAuxGetLiqArealDensity
     procedure, public :: GetIceArealDensity                => ThermalKSPTemperatureSoilAuxGetIceArealDensity
     procedure, public :: GetSnowWater                      => ThermalKSPTemperatureSoilAuxGetSnowWater
     procedure, public :: GetNumberOfSnowLayers             => ThermalKSPTemperatureSoilAuxGetNumberOfSnowLayers
     procedure, public :: GetTuningFactor                   => ThermalKSPTemperatureSoilAuxGetTuningFactor
     procedure, public :: GetDepth                          => ThermalKSPTemperatureSoilAuxGetDepth
     procedure, public :: GetPorosity                       => ThermalKSPTemperatureSoilAuxGetPorosity
     procedure, public :: GetThermalCondMinerals            => ThermalKSPTemperatureSoilAuxGetThermalCondMinerals
     procedure, public :: GetThermalCondDrySoil             => ThermalKSPTemperatureSoilAuxGetThermalCondDrySoil
     procedure, public :: GetVolumetricHeatCapacityMinerals => ThermalKSPTemperatureSoilAuxGetVolumetricHeatCapacityMinerals
     procedure, public :: GetColumnType                     => ThermalKSPTemperatureSoilAuxGetColumnType

  end type therm_ksp_temp_soil_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_soil_auxvar_type)   :: this

    call this%BaseInit()

    this%liq_areal_den         = 0.d0
    this%ice_areal_den         = 0.d0
    this%snow_water            = 0.d0
    this%num_snow_layer        = 0
    this%tuning_factor         = 1.d0
    this%dz                    = 0.d0

    this%por                   = 0.d0
    this%therm_cond_minerals   = 0.d0
    this%therm_cond_dry        = 0.d0
    this%heat_cap_minerals_puv = 0.d0

    this%is_soil_shallow       = PETSC_FALSE
    this%itype                 = -1

  end subroutine ThermKSPTempSoilAuxVarInit

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSoilAuxVarCompute(this, dz, vol)
    !
    ! !DESCRIPTION:
    !
    use mpp_varcon , only : denh2o, denice, tfrz, tkwat, tkice, cpice,  cpliq, thk_bedrock
    use mpp_varcon , only : istice, istice_mec, istwet, istsoil, istcrop
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_soil_auxvar_type)   :: this
    PetscReal                                :: dz
    PetscReal                                :: vol
    !
    ! LOCAL VARIABLES
    PetscReal :: satw
    PetscReal :: fl
    PetscReal :: dke
    PetscReal :: dksat
    PetscReal :: heat_cap_pva

    !select case(this%itype)
    !case (istsoil, istcrop)
    if (this%itype == istsoil .or. this%itype == istcrop) then
       if (this%is_soil_shallow) then

          satw = (this%liq_areal_den/denh2o + this%ice_areal_den/denice)/(dz*this%por)
          satw = min(1.d0, satw)

          if (satw > .1e-6) then
             if (this%GetTemperature() >= tfrz) then       ! Unfrozen soil
                dke = max(0.d0, log10(satw) + 1.0d0)
             else                                      ! Frozen soil
                dke = satw
             end if

             fl = (this%liq_areal_den/(denh2o*dz)) / (this%liq_areal_den/(denh2o*dz) + &
                                                      this%ice_areal_den/(denice*dz))

             dksat = this%therm_cond_minerals*    &
                  tkwat**(fl*this%por)*        &
                  tkice**((1.d0-fl)*this%por)

             call this%SetThermalConductivity( dke*dksat + (1.d0-dke)*this%therm_cond_dry)

          else
             call this%SetThermalConductivity(this%therm_cond_dry)
          endif

          heat_cap_pva   = this%heat_cap_minerals_puv*(1.d0 - this%por)*dz + &
                           this%ice_areal_den*cpice                         + &
                           this%liq_areal_den*cpliq

          if (this%num_snow_layer == 0) then
             heat_cap_pva = heat_cap_pva + this%snow_water*cpice
          endif

       else

          call this%SetThermalConductivity(thk_bedrock)
          heat_cap_pva   = this%heat_cap_minerals_puv*(1.d0 - this%por)*dz + &
                           this%ice_areal_den*cpice                         + &
                           this%liq_areal_den*cpliq
       endif

       call this%SetVolumetricHeatCapacity(heat_cap_pva/dz)

    else if (this%itype == istwet) then

       if (this%is_soil_shallow) then
          if (this%GetTemperature() < tfrz) then
             call this%SetThermalConductivity(tkice)
          else
             call this%SetThermalConductivity(tkwat)
          endif
          heat_cap_pva = this%ice_areal_den*cpice + this%liq_areal_den*cpliq
          if (this%num_snow_layer == 0) then
             heat_cap_pva = heat_cap_pva + this%snow_water*cpice
          endif

          call this%SetVolumetricHeatCapacity(heat_cap_pva/dz)
       else
          call this%SetThermalConductivity(thk_bedrock)
          call this%SetVolumetricHeatCapacity(this%heat_cap_minerals_puv)
       endif

    else if (this%itype == istice .or. this%itype == istice_mec) then
          if (this%GetTemperature() < tfrz) then
             call this%SetThermalConductivity(tkice)
          else
             call this%SetThermalConductivity(tkwat)
          endif
          heat_cap_pva = this%ice_areal_den*cpice + this%liq_areal_den*cpliq
          if (this%num_snow_layer == 0) then
             heat_cap_pva = heat_cap_pva + this%snow_water*cpice
          endif

          call this%SetVolumetricHeatCapacity(heat_cap_pva/dz)
    end if
       

  end subroutine ThermKSPTempSoilAuxVarCompute

  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetLiqArealDensity(this, val)
    !
    ! !DESCRIPTION:
    ! Set liquid areal density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%liq_areal_den = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetLiqArealDensity
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetIceArealDensity(this, val)
    !
    ! !DESCRIPTION:
    ! Set ice areal density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%ice_areal_den = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetIceArealDensity
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetSnowWater(this, val)
    !
    ! !DESCRIPTION:
    ! Set snow water
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%snow_water = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetSnowWater
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetNumberOfSnowLayers(this, val)
    !
    ! !DESCRIPTION:
    ! Set of number snow layers
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscInt                        :: val
 
    this%num_snow_layer = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetNumberOfSnowLayers
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetTuningFactor(this, val)
    !
    ! !DESCRIPTION:
    ! Set tuning factor
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%tuning_factor = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetTuningFactor
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetDepth(this, val)
    !
    ! !DESCRIPTION:
    ! Set depth
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%dz = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetDepth
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetPorosity(this, val)
    !
    ! !DESCRIPTION:
    ! Set porosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%por = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetPorosity
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetThermalCondMinerals(this, val)
    !
    ! !DESCRIPTION:
    ! Set thermal conductivity of minerals
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%therm_cond_minerals = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetThermalCondMinerals
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetThermalCondDrySoil(this, val)
    !
    ! !DESCRIPTION:
    ! Set thermal conductivity of dry soil
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%therm_cond_dry = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetThermalCondDrySoil
 
  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetVolumetricHeatCapacityMinerals(this, val)
    !
    ! !DESCRIPTION:
    ! Set volumetric heat capacity of minerals
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscReal                        :: val
 
    this%heat_cap_minerals_puv = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetVolumetricHeatCapacityMinerals

  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetColumnType(this, val)
    !
    ! !DESCRIPTION:
    ! Set of number snow layers
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    PetscInt                        :: val
 
    this%itype = val
 
  end subroutine ThermalKSPTemperatureSoilAuxSetColumnType

  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetSoilShallow(this)
    !
    ! !DESCRIPTION:
    ! Set soil is shallow
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
 
    this%is_soil_shallow = PETSC_TRUE
 
  end subroutine ThermalKSPTemperatureSoilAuxSetSoilShallow

  !------------------------------------------------------------------------
  subroutine ThermalKSPTemperatureSoilAuxSetSoilDeep(this)
    !
    ! !DESCRIPTION:
    ! Set soil is deep
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
 
    this%is_soil_shallow = PETSC_FALSE
 
  end subroutine ThermalKSPTemperatureSoilAuxSetSoilDeep

  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetLiqArealDensity(this)
    !
    ! !DESCRIPTION:
    ! Get liquid areal density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetLiqArealDensity
 
    ThermalKSPTemperatureSoilAuxGetLiqArealDensity = this%liq_areal_den
 
  end function ThermalKSPTemperatureSoilAuxGetLiqArealDensity
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetIceArealDensity(this)
    !
    ! !DESCRIPTION:
    ! Get ice areal density
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetIceArealDensity
 
    ThermalKSPTemperatureSoilAuxGetIceArealDensity = this%ice_areal_den
 
  end function ThermalKSPTemperatureSoilAuxGetIceArealDensity
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetSnowWater(this)
    !
    ! !DESCRIPTION:
    ! Get snow water
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetSnowWater
 
    ThermalKSPTemperatureSoilAuxGetSnowWater = this%snow_water
 
  end function ThermalKSPTemperatureSoilAuxGetSnowWater
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetNumberOfSnowLayers(this)
    !
    ! !DESCRIPTION:
    ! Get snow water
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscInt                         :: ThermalKSPTemperatureSoilAuxGetNumberOfSnowLayers
 
    ThermalKSPTemperatureSoilAuxGetNumberOfSnowLayers = this%num_snow_layer
 
  end function ThermalKSPTemperatureSoilAuxGetNumberOfSnowLayers
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetTuningFactor(this)
    !
    ! !DESCRIPTION:
    ! Get tuning factor
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetTuningFactor
 
    ThermalKSPTemperatureSoilAuxGetTuningFactor = this%tuning_factor
 
  end function ThermalKSPTemperatureSoilAuxGetTuningFactor
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetDepth(this)
    !
    ! !DESCRIPTION:
    ! Get depth
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetDepth
 
    ThermalKSPTemperatureSoilAuxGetDepth = this%dz
 
  end function ThermalKSPTemperatureSoilAuxGetDepth
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetPorosity(this)
    !
    ! !DESCRIPTION:
    ! Get porosity
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetPorosity
 
    ThermalKSPTemperatureSoilAuxGetPorosity = this%por
 
  end function ThermalKSPTemperatureSoilAuxGetPorosity
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetThermalCondMinerals(this)
    !
    ! !DESCRIPTION:
    ! Get thermal conductivity of minerals
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetThermalCondMinerals
 
    ThermalKSPTemperatureSoilAuxGetThermalCondMinerals = this%therm_cond_minerals
 
  end function ThermalKSPTemperatureSoilAuxGetThermalCondMinerals
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetThermalCondDrySoil(this)
    !
    ! !DESCRIPTION:
    ! Get thermal conductivity of dry soil
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetThermalCondDrySoil
 
    ThermalKSPTemperatureSoilAuxGetThermalCondDrySoil = this%therm_cond_dry
 
  end function ThermalKSPTemperatureSoilAuxGetThermalCondDrySoil
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetVolumetricHeatCapacityMinerals(this)
    !
    ! !DESCRIPTION:
    ! Get volumetric heat capacity of minerals
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscReal                        :: ThermalKSPTemperatureSoilAuxGetVolumetricHeatCapacityMinerals
 
    ThermalKSPTemperatureSoilAuxGetVolumetricHeatCapacityMinerals = this%heat_cap_minerals_puv
 
  end function ThermalKSPTemperatureSoilAuxGetVolumetricHeatCapacityMinerals

  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetColumnType(this)
    !
    ! !DESCRIPTION:
    ! Get snow water
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscInt                        :: ThermalKSPTemperatureSoilAuxGetColumnType
 
    ThermalKSPTemperatureSoilAuxGetColumnType = this%itype
 
  end function ThermalKSPTemperatureSoilAuxGetColumnType
 
  !------------------------------------------------------------------------
  function ThermalKSPTemperatureSoilAuxGetIsSoilShallow(this)
    !
    ! !DESCRIPTION:
    ! Returns true if soil is shallow
    !
    implicit none
    !
    ! !ARGUMENTS
    !
    class(therm_ksp_temp_soil_auxvar_type) :: this
    !
    PetscBool                       :: ThermalKSPTemperatureSoilAuxGetIsSoilShallow
 
    ThermalKSPTemperatureSoilAuxGetIsSoilShallow = this%is_soil_shallow
 
  end function ThermalKSPTemperatureSoilAuxGetIsSoilShallow

#endif
  
end module ThermalKSPTemperatureSoilAuxType
