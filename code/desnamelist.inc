        NAMELIST / DES_INPUT_DATA_LOCKED /&

! USE discretelement
        DES_ONEWAY_COUPLED, DES_REPORT_MASS_INTERP, DES_INTG_METHOD, &
        PARTICLES, NFACTOR, USE_VDH_DEM_MODEL, &
        GENER_PART_CONFIG, DES_NEIGHBOR_SEARCH, NEIGHBOR_SEARCH_N,&
        FACTOR_RLM, NEIGHBOR_SEARCH_RAD_RATIO, DES_USR_VAR_SIZE,&

        DES_INTERP_SCHEME, DES_INTERP_WIDTH, DES_INTERP_ON, &
        DES_INTERP_MEAN_FIELDS, DES_EXPLICITLY_COUPLED, &
        DES_DIFFUSE_WIDTH, &

! Used for grid/cell based search methods
        DESGRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, DESGRIDSEARCH_KMAX, &

! Particle properties
        KN, KN_W, KT_FAC, KT_W_FAC, MEW, MEW_W, MEW_R, MEW_RW, &
        DES_EN_INPUT, DES_ET_INPUT, DES_ETAT_FAC, &
        DES_EN_WALL_INPUT, DES_ET_WALL_INPUT, DES_ETAT_W_FAC, &
!
! Additional properties for Hertzian model:
        DES_COLL_MODEL, EW_YOUNG, VW_POISSON, E_YOUNG, V_POISSON,&

! Actual material properties for particles/wall.  Used for computing
! conduction correction terms to account for artificial softening
    EW_YOUNG_ACTUAL, E_YOUNG_ACTUAL, V_POISSON_ACTUAL, &
  VW_POISSON_ACTUAL, &

! Lees Edwards BC (for DEM only)
        DES_LE_BC, DES_LE_REL_VEL, DES_LE_SHEAR_DIR,&


! DES wall BC: wall translational velocity
! These need to be improved to work for arbitrary geoms
        DES_BC_Uw_s, DES_BC_Vw_s, DES_BC_Ww_s,&

! DES flow boundary condition information: J.Musser
        FORCE_ORD_BC, &

! DES output control information
        DEBUG_DES, FOCUS_PARTICLE, VTP_DIR, &
        PRINT_DES_DATA, DES_OUTPUT_TYPE,  &

! DES Cohesion Information: M.Weber
        USE_COHESION,VAN_DER_WAALS, HAMAKER_CONSTANT, &
        WALL_HAMAKER_CONSTANT, &
        VDW_INNER_CUTOFF, VDW_OUTER_CUTOFF, ASPERITIES, &
        WALL_VDW_INNER_CUTOFF, WALL_VDW_OUTER_CUTOFF,&

! JKR adhesive DEM model (spherical particles only)
        USE_JKR, DES_SURFACE_ENERGY, DES_SURFACE_ENERGY_WALL, &
        JKR_THETAR_CRIT, WALL_JKR_THETAR_CRIT, &
        JKR_ALPHA_SAFETY, &


! DES Energy Equations: J.Musser
        DES_CONV_CORR, FLPC, DES_Em, DES_MIN_COND_DIST, &

! MP-PIC related flags
         PSFAC_FRIC_PIC, FRIC_EXP_PIC, FRIC_NON_SING_FAC,&
         MPPIC_COEFF_EN1, MPPIC_VELFAC_COEFF, &
         MPPIC_COEFF_EN_WALL, MPPIC_COEFF_ET_WALL, &
         PIC_CFL, PIC_CFL_CONTROL, PIC_CFL_PARCEL_FRACTION, &
         PIC_COLLISION_DAMPING, PIC_CD_E, &

! Coarse Grain DEM
         CGP_STAT_WT, CGP_D_P0, &

! SuperDEM
      SQP_a,  SQP_b,  SQP_c,  SQP_m, SQP_n, &
      SQP_q1, SQP_q2, SQP_q3, SQP_q4, &
      SQP_contact_evolution, SQP_init_contact_evolution, &

!GluedSphereDEM
      GSP_q1, GSP_q2, GSP_q3, GSP_q4, NGluedParticles, init_theta, &
      GSP_MODEL, flexible_gsp,&

!GluedSphereDEM with virtual bonds
      E_YOUNG_BOND, C_ratio_bond, V_POISSON_BOND,&
      TENSILE_STRENGTH_BOND, SHEAR_STRENGTH_BOND,&
      FLEXIBLE_GSP_BOND_BREAK,&

! Flag to turn on/off optimizing the list of facets at each des grid cell
         MINIMIZE_DES_FACET_LIST, &

! Time interval at which particle dynamic load balance is called, Eulerian Grid Weight
         DLB_DT,DLB_EGW, &

! Advances settings to control DEM timestep
         DTSOLID_FACTOR, DTSOLID_UPDATE_DT, &

! Flag to prevent particles from falling back in a top exit

         FREEZE_DES_EXIT_VEL, &

! Flag to remove rogue particles
           REMOVE_ROGUE_PARTICLES, &

! Growth factor when resizing send/recv buffers
         DES_BUFF_RESIZE_FACTOR, &

! Void fraction clipping
         DES_EPG_CLIP, &

! Specified particle velocity
         DES_RIGID_MOTION, &
         DES_RIGID_MOTION_VEL, DES_RIGID_MOTION_OMEGA, DES_RIGID_MOTION_ROT_CENTER
