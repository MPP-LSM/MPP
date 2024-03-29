
module MultiPhysicsProbConstants

#ifdef USE_PETSC_LIB

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Constants for multi-physcis problems
  !-----------------------------------------------------------------------
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include <petsc/finclude/petsc.h>

  PetscInt, parameter, public :: DISCRETIZATION_VERTICAL_ONLY     = 1
  PetscInt, parameter, public :: DISCRETIZATION_HORIZONTAL_ONLY   = 2
  PetscInt, parameter, public :: DISCRETIZATION_THREE_DIM         = 3
  PetscInt, parameter, public :: DISCRETIZATION_VERTICAL_WITH_SS  = 4

  ! mpp_itype
  PetscInt, parameter, public :: MPP_VSFM_SNES_CLM                 = 11
  PetscInt, parameter, public :: MPP_THERMAL_TBASED_KSP_CLM        = 12
  PetscInt, parameter, public :: MPP_THERMAL_EBASED_SNES_CLM       = 13
  PetscInt, parameter, public :: MPP_TH_SNES_CLM                   = 14
  PetscInt, parameter, public :: MPP_MLC_KSP                       = 15
  PetscInt, parameter, public :: MPP_LBL_KSP                       = 16
  PetscInt, parameter, public :: MPP_PHOTOSYNTHESIS_SNES           = 17
  PetscInt, parameter, public :: MPP_LONGWAVE_KSP                  = 18
  PetscInt, parameter, public :: MPP_SHORTWAVE_KSP                 = 19

  ! soe_itype
  PetscInt, parameter, public :: SOE_RE_ODE                        = 101
  PetscInt, parameter, public :: SOE_THERMAL_TBASED                = 102
  PetscInt, parameter, public :: SOE_THERMAL_EBASED                = 103
  PetscInt, parameter, public :: SOE_TH                            = 104
  PetscInt, parameter, public :: SOE_MLC                           = 105
  PetscInt, parameter, public :: SOE_LBL                           = 106
  PetscInt, parameter, public :: SOE_PHOTOSYNTHESIS                = 107
  PetscInt, parameter, public :: SOE_LONGWAVE                      = 108
  PetscInt, parameter, public :: SOE_SHORTWAVE                     = 109

  ! ge_itype
  PetscInt, parameter, public :: GE_RE                             = 201
  PetscInt, parameter, public :: GE_THERM_SOIL_TBASED              = 202
  PetscInt, parameter, public :: GE_THERM_SNOW_TBASED              = 203
  PetscInt, parameter, public :: GE_THERM_SSW_TBASED               = 204
  PetscInt, parameter, public :: GE_THERM_SOIL_EBASED              = 205
  PetscInt, parameter, public :: GE_CANOPY_AIR_TEMP                = 206
  PetscInt, parameter, public :: GE_CANOPY_AIR_VAPOR               = 207
  PetscInt, parameter, public :: GE_CANOPY_LEAF_TEMP               = 208
  PetscInt, parameter, public :: GE_LEAF_BND_LAYER                 = 209
  PetscInt, parameter, public :: GE_PHOTOSYNTHESIS                 = 210
  PetscInt, parameter, public :: GE_LONGWAVE                       = 211
  PetscInt, parameter, public :: GE_SHORTWAVE                      = 212

  ! mesh_itype
  PetscInt, parameter, public :: MESH_CLM_SOIL_COL                 = 301
  PetscInt, parameter, public :: MESH_CLM_THERMAL_SOIL_COL         = 302
  PetscInt, parameter, public :: MESH_CLM_SNOW_COL                 = 303
  PetscInt, parameter, public :: MESH_CLM_SSW_COL                  = 304
  PetscInt, parameter, public :: MESH_SPAC_ROOT_COL                = 305
  PetscInt, parameter, public :: MESH_SPAC_XYLEM_COL               = 306
  PetscInt, parameter, public :: MESH_ALONG_GRAVITY                = 311
  PetscInt, parameter, public :: MESH_AGAINST_GRAVITY              = 312

  ! region_itype
  PetscInt, parameter, public :: SOIL_TOP_CELLS                    = 401
  PetscInt, parameter, public :: SOIL_BOTTOM_CELLS                 = 402
  PetscInt, parameter, public :: SOIL_CELLS                        = 403
  PetscInt, parameter, public :: SOIL_CELLS_OF_NEIGH_MESH          = 404
  PetscInt, parameter, public :: SNOW_TOP_CELLS                    = 405
  PetscInt, parameter, public :: SNOW_BOTTOM_CELLS                 = 406
  PetscInt, parameter, public :: SSW_TOP_CELLS                     = 407
  PetscInt, parameter, public :: ALL_CELLS                         = 408
  PetscInt, parameter, public :: DEFINED_BY_CELL_ID                = 409
  PetscInt, parameter, public :: FACE_TOP                          = 410
  PetscInt, parameter, public :: FACE_BOTTOM                       = 411

  ! condition%itype
  PetscInt, parameter, public :: COND_NULL                         = 500
  PetscInt, parameter, public :: COND_BC                           = 501
  PetscInt, parameter, public :: COND_SS                           = 502
  PetscInt, parameter, public :: COND_MASS_RATE                    = 503
  PetscInt, parameter, public :: COND_MASS_FLUX                    = 504
  PetscInt, parameter, public :: COND_DIRICHLET                    = 505
  PetscInt, parameter, public :: COND_DIRICHLET_FRM_OTR_GOVEQ      = 506
  PetscInt, parameter, public :: COND_HEAT_FLUX                    = 507
  PetscInt, parameter, public :: COND_DARCY_RATE                   = 508
  PetscInt, parameter, public :: COND_SEEPAGE_BC                   = 509
  PetscInt, parameter, public :: COND_HEAT_RATE                    = 511
  PetscInt, parameter, public :: COND_DOWNREG_MASS_RATE_CAMPBELL   = 512
  PetscInt, parameter, public :: COND_DOWNREG_MASS_RATE_FETCH2     = 513

  !
  PetscInt, parameter, public :: VAR_XI                                     = 601
  PetscInt, parameter, public :: VAR_DXI_DP                                 = 602
  PetscInt, parameter, public :: VAR_DENxSATxPOR                            = 611
  PetscInt, parameter, public :: VAR_DDENxSATxPOR_DP                        = 612
  PetscInt, parameter, public :: VAR_DXI_DTIME                              = 603
  PetscInt, parameter, public :: VAR_PRESSURE                               = 604
  PetscInt, parameter, public :: VAR_TEMPERATURE                            = 605
  PetscInt, parameter, public :: VAR_PRESSURE_PREV                          = 606
  PetscInt, parameter, public :: VAR_BC_SS_CONDITION                        = 607
  PetscInt, parameter, public :: VAR_LIQ_SAT                                = 608
  PetscInt, parameter, public :: VAR_DENSITY_TYPE                           = 609
  PetscInt, parameter, public :: VAR_MASS                                   = 610
  PetscInt, parameter, public :: VAR_SOIL_MATRIX_POT                        = 611
  PetscInt, parameter, public :: VAR_FRAC_LIQ_SAT                           = 612
  PetscInt, parameter, public :: VAR_LATERAL_MASS_EXCHANGED                 = 613
  PetscInt, parameter, public :: VAR_BC_MASS_EXCHANGED                      = 614
  PetscInt, parameter, public :: VAR_LIQ_AREAL_DEN                          = 615
  PetscInt, parameter, public :: VAR_ICE_AREAL_DEN                          = 617
  PetscInt, parameter, public :: VAR_FRAC                                   = 618
  PetscInt, parameter, public :: VAR_SNOW_WATER                             = 619
  PetscInt, parameter, public :: VAR_NUM_SNOW_LYR                           = 620
  PetscInt, parameter, public :: VAR_DHS_DT                                 = 621
  PetscInt, parameter, public :: VAR_THERMAL_COND                           = 622
  PetscInt, parameter, public :: VAR_HEAT_CAP                               = 623
  PetscInt, parameter, public :: VAR_ACTIVE                                 = 624
  PetscInt, parameter, public :: VAR_DX                                     = 625
  PetscInt, parameter, public :: VAR_DY                                     = 626
  PetscInt, parameter, public :: VAR_DZ                                     = 627
  PetscInt, parameter, public :: VAR_DIST_UP                                = 628
  PetscInt, parameter, public :: VAR_DIST_DN                                = 629
  PetscInt, parameter, public :: VAR_TUNING_FACTOR                          = 630
  PetscInt, parameter, public :: VAR_XC                                     = 631
  PetscInt, parameter, public :: VAR_YC                                     = 632
  PetscInt, parameter, public :: VAR_ZC                                     = 633
  PetscInt, parameter, public :: VAR_AREA                                   = 634
  PetscInt, parameter, public :: VAR_VOLUME                                 = 635
  PetscInt, parameter, public :: VAR_CONDUCTANCE                            = 636
  PetscInt, parameter, public :: VAR_FLUX_TYPE                              = 637
  PetscInt, parameter, public :: VAR_POT_MASS_SINK_PRESSURE                 = 638
  PetscInt, parameter, public :: VAR_POT_MASS_SINK_EXPONENT                 = 639
  PetscInt, parameter, public :: VAR_PRESSURE_UP                            = 640
  PetscInt, parameter, public :: VAR_PRESSURE_DN                            = 641
  PetscInt, parameter, public :: VAR_CAMPBELL_HE                            = 642
  PetscInt, parameter, public :: VAR_CAMPBELL_N                             = 643
  PetscInt, parameter, public :: VAR_MASS_FLUX                              = 644
  PetscInt, parameter, public :: VAR_CONDUCTANCE_TYPE                       = 645
  PetscInt, parameter, public :: VAR_CONDUCTANCE_UP                         = 646
  PetscInt, parameter, public :: VAR_CONDUCTANCE_DN                         = 647
  PetscInt, parameter, public :: VAR_WATER_VAPOR                            = 648
  PetscInt, parameter, public :: VAR_LEAF_TEMPERATURE                       = 649
  PetscInt, parameter, public :: VAR_LEAF_BDN_LYR_COND_HEAT                 = 650
  PetscInt, parameter, public :: VAR_LEAF_BDN_LYR_COND_H2O                  = 651
  PetscInt, parameter, public :: VAR_LEAF_BDN_LYR_COND_CO2                  = 652
  PetscInt, parameter, public :: VAR_SCM_MEDLYN                             = 653
  PetscInt, parameter, public :: VAR_SCM_BBERRY                             = 654
  PetscInt, parameter, public :: VAR_SCM_WUE                                = 655
  PetscInt, parameter, public :: VAR_SCM_BONAN14                            = 656
  PetscInt, parameter, public :: VAR_SCM_MANZONI11                          = 657
  PetscInt, parameter, public :: VAR_SCM_MODIFIED_BONAN14                   = 658
  PetscInt, parameter, public :: VAR_PHOTOSYNTHETIC_PATHWAY_C4              = 659
  PetscInt, parameter, public :: VAR_PHOTOSYNTHETIC_PATHWAY_C3              = 660
  PetscInt, parameter, public :: VAR_STOMATAL_CONDUCTANCE                   = 661
  PetscInt, parameter, public :: VAR_LEAF_ABSORBED_SHORTWAVE_RAD_PER_LAI    = 662
  PetscInt, parameter, public :: VAR_SOIL_ABSORBED_SHORTWAVE_RAD_PER_GROUND = 663
  PetscInt, parameter, public :: VAR_LEAF_ABSORBED_LONGWAVE_RAD_PER_LAI     = 664
  PetscInt, parameter, public :: VAR_SOIL_ABSORBED_LONGWAVE_RAD_PER_GROUND  = 665
  PetscInt, parameter, public :: VAR_GROSS_PHOTOSYNTHESIS                   = 666
  PetscInt, parameter, public :: VAR_NET_PHOTOSYNTHESIS                     = 667
  PetscInt, parameter, public :: VAR_LEAF_HEAT_STORAGE                      = 668
  PetscInt, parameter, public :: VAR_LATENT_HEAT_FLUX                       = 669
  PetscInt, parameter, public :: VAR_SENSIBLE_HEAT_FLUX                     = 670
  PetscInt, parameter, public :: VAR_LEAF_TRANSPIRATION                     = 671
  PetscInt, parameter, public :: VAR_POT_SINK_DOWNREG_FACTOR                = 672
  PetscInt, parameter, public :: VAR_SCM_OSMWANG                            = 673
  !
  PetscInt, parameter, public :: AUXVAR_INTERNAL                   = 701
  PetscInt, parameter, public :: AUXVAR_BC                         = 702
  PetscInt, parameter, public :: AUXVAR_SS                         = 703
  PetscInt, parameter, public :: AUXVAR_CONN_INTERNAL              = 704
  PetscInt, parameter, public :: AUXVAR_CONN_BC                    = 705
  PetscInt, parameter, public :: AUXVAR_BC_OTR_GOVEQ               = 706

  PetscInt, parameter, public :: PETSC_TS                          = 801
  PetscInt, parameter, public :: PETSC_SNES                        = 802
  PetscInt, parameter, public :: PETSC_KSP                         = 803

  PetscInt, parameter, public :: CONN_VERTICAL                     = 901
  PetscInt, parameter, public :: CONN_HORIZONTAL                   = 902
  PetscInt, parameter, public :: CONN_SET_INTERNAL                 = 903
  PetscInt, parameter, public :: CONN_SET_LATERAL                  = 904
  PetscInt, parameter, public :: CONN_SET_CONDITIONS               = 905
  PetscInt, parameter, public :: CONN_IN_X_DIR                     = 906
  PetscInt, parameter, public :: CONN_IN_Y_DIR                     = 907
  PetscInt, parameter, public :: CONN_IN_Z_DIR                     = 908
  PetscInt, parameter, public :: CONN_IN_XYZ_DIR                   = 909

  PetscInt, parameter, public :: DARCY_FLUX_TYPE                   = 1001
  PetscInt, parameter, public :: CONDUCTANCE_FLUX_TYPE             = 1002
  PetscInt, parameter, public :: CONDUCTANCE_CAMPBELL_TYPE         = 1003
  PetscInt, parameter, public :: CONDUCTANCE_MANOLI_TYPE           = 1004

  !
  PetscReal, parameter, public :: PRESSURE_REF                     = 101325.d0     ! [Pa]
  !PetscReal, parameter, public :: GRAVITY_CONSTANT                 = 9.8068d0      ! [m s^{-2}]
  PetscReal, parameter, public :: GRAVITY_CONSTANT                 = 9.80665d0      ! [m s^{-2}]
  PetscReal, parameter, public :: FMWH2O                           = 18.01534d0    ! [kg kmol^{-1}]
  PetscReal, parameter, public :: STEFAN_BOLTZMAN_CONSTANT         = 5.67d-08

  PetscReal, parameter, public :: VKC                              = 0.4d0         ! von Karman constant [-]
  PetscReal, parameter, public :: TFRZ                             = 273.15d0      ! [K]
  PetscReal, parameter, public :: MM_H2O                           = 18.02d-3      ! [kgmol{-1}]
  PetscReal, parameter, public :: MM_DRY_AIR                       = 28.97d-3      ! [kg mol^{-1}]
  PetscReal, parameter, public :: HVAP                             = 2.501e6       ! [J/kg]
  PetscReal, parameter, public :: HSUB                             = 2.8347d6      ! latent heat of evaporation (J/kg)
  PetscReal, parameter, public :: CPD                              = 1005.d0       ! [J/kg/K]
  PetscReal, parameter, public :: CPW                              = 1846.d0       ! [J/kg/K]
  PetscReal, parameter, public :: RGAS                             = 8.31446d0     ! [J/K/mol]
  PetscReal, parameter, public :: VISC_0C                          = 13.3d-06      ! [m^2/s]
  PetscReal, parameter, public :: MOD_DIFF_HEAT_OC                 = 18.9d-6       ! [m^2/s]
  PetscReal, parameter, public :: MOD_DIFF_H2O_OC                  = 21.8d-6       ! [m^2/s]
  PetscReal, parameter, public :: MOD_DIFF_CO2_OC                  = 13.8d-6       ! [m^2/s]
#endif


end module MultiPhysicsProbConstants
