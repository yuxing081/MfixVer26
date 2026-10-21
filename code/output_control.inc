!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: OUTPUT_CONTROL_LOCKED                                     !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   NAMELIST /OUTPUT_CONTROL_LOCKED/ &
      SPX_DT, nRR, PHIP_OUT_JJ, bDist_IO, bStart_with_one_RES, &
      DBGPRN_LAYOUT, ENABLE_DMP_LOG, RES_BACKUPS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: OUTPUT_CONTROL_UNLOCKED                                   !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   NAMELIST /OUTPUT_CONTROL_UNLOCKED/ &
      RES_DT, OUT_DT, NLOG, FULL_LOG, RESID_STRING, GROUP_RESID, &
      TSS_HDR, &
      REPORT_MASS_BALANCE_DT, REPORT_NEG_DENSITY, RES_BACKUP_DT, &
      REPORT_NEG_SPECIFICHEAT, REPORT_SOLID_INVENTORY, &
      REPORT_SOLID_INVENTORY_DT, BREAKDOWN_SOLID_INVENTORY_BY_PHASE, &
      VTK_DT, VTK_T_START, VTK_T_END, VTK_DBG_FILE, VTK_GSP_OUTPUT_TYPE,&
      VTK_VAR, FRAME, VTU_DIR, &
      VTK_X_w, VTK_X_e, VTK_Y_s, VTK_Y_n, VTK_Z_b, VTK_Z_t, VTK_VARLIST, &
      VTK_EP_g, VTK_RO_g, VTK_Mu_g, VTK_k_g, &
      VTK_P_g, VTK_P_star, VTK_P_s, VTK_VEL_g, VTK_VEL_S, &
      VTK_U_G, VTK_V_G, VTK_W_G, VTK_U_S, VTK_V_S, VTK_W_S, &
      VTK_ROP_s, VTK_RO_S, VTK_EP_S, &
      VTK_T_g, VTK_T_s, VTK_C_Pg, VTK_C_Ps, VTK_X_g, VTK_X_s, &
      VTK_MW_MIX_g, VTK_Y_g, &
      VTK_Theta_m, VTK_Scalar, VTK_RRate, VTK_RRate_label, &
      VTK_FLUID_RRATE, VTK_DES_RRATE, &
      VTK_K_Turb_G, VTK_E_Turb_G, &
      VTK_VORTICITY, VTK_LAMBDA_2, VTK_PARTITION, VTK_BC_ID, VTK_DWALL, &
      VTK_FACET_COUNT_DES, VTK_NB_FACET_DES, VTK_IJK, &
      VTK_SCALAR_CC, VTK_NORMAL, VTK_DEBUG, &
      VTK_PART_DIAMETER, VTK_PART_VEL, VTK_PART_ANGULAR_VEL, VTK_PART_ORIENTATION, &
      VTK_PART_PHYSICAL_DIAMETER, VTK_PART_CGP_STAT_WT, &
      VTK_PART_USR_VAR, VTK_PART_TEMP, VTK_PART_C_Ps, &
      VTK_PART_X_s, VTK_PART_COHESION, &
      VTK_PART_RRATE, VTK_PART_DENSITY, &
      VTK_PART_MASS, VTK_PART_VOLUME, VTK_PART_DRAG_COEFF, &
      VTK_PART_RANK, VTK_PART_ID, VTK_PART_PHASE, VTK_PART_PHASE_ID, &
      VTK_DOMAIN_DECOMPOSITION, VTK_PART_RESIDENCE_TIME, &
      VTK_FCMAX_ONLY, VTK_FC_JKR_STATE, VTK_FC_JKR_RADIUS, &
      VTK_FC_JKR_GAMMA, VTK_FC_TANG_DISP, VTK_FC_ROLL_DISP, &
      VTK_FC_TWIST_ANGLE, &
      VTK_PART_COLL_FORCE, &
      VTK_PART_GLUEID, VTK_PART_GLUEDIAMETER, VTK_PART_GLUEBOUNDING, VTK_PART_GLUEMASS, &
      VTK_PART_GLUEVOLUME, VTK_PART_GLUEQUAT, VTK_FILEBASE, VTK_DATA, VTK_GEO, VTK_CUTCELL_ONLY, &
      VTK_PART_COORDINATION, VTK_PART_OVERLAP_STATS, &
      VTK_NXS, VTK_NYS, VTK_NZS, VTK_SLICE_TOL, VTK_SELECT_MODE, &
      WRITE_PART_OUT, PART_OUT_ZERO_VEL, &
      PART_OUT_X_MIN, PART_OUT_X_MAX,PART_OUT_X_EXCLUDE, &
      PART_OUT_Y_MIN, PART_OUT_Y_MAX,PART_OUT_Y_EXCLUDE, &
      PART_OUT_Z_MIN, PART_OUT_Z_MAX,PART_OUT_Z_EXCLUDE, &
      PART_OUT_PHASE, &
      PART_OUT_DIAMETER_MIN, PART_OUT_DIAMETER_MAX, PART_OUT_DIAMETER_EXCLUDE, &
      PART_OUT_DENSITY_MIN, PART_OUT_DENSITY_MAX, PART_OUT_DENSITY_EXCLUDE, &
      PART_OUT_U_MIN, PART_OUT_U_MAX, PART_OUT_U_EXCLUDE, &
      PART_OUT_V_MIN, PART_OUT_V_MAX, PART_OUT_V_EXCLUDE, &
      PART_OUT_W_MIN, PART_OUT_W_MAX, PART_OUT_W_EXCLUDE, &
      PART_OUT_TEMP_MIN, PART_OUT_TEMP_MAX, PART_OUT_TEMP_EXCLUDE, &
      PART_OUT_X_S_MIN, PART_OUT_X_S_MAX, PART_OUT_X_S_EXCLUDE, &
      PART_OUT_USR_VAR_MIN, PART_OUT_USR_VAR_MAX, PART_OUT_USR_VAR_EXCLUDE, &
      STL_TORQUE_PIVOT_BC, STL_TORQUE_PIVOT_IS, STL_TORQUE_PIVOT_DEFAULT_BC, &
      VTK_STL_GROUP_ID, VTK_STL_COLL_COUNT, VTK_STL_Q, VTK_STL_AREA, &
      VTK_STL_WEAR, VTK_STL_NORMAL_FORCE, VTK_STL_TANG_FORCE, &
      VTK_STL_FORCE, VTK_STL_TORQUE, &
      VTK_STL_NORMAL_VEL, VTK_STL_TANG_VEL, VTK_STL_VEL
