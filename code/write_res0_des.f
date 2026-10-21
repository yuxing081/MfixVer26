MODULE WRITE_RES0_DES_MOD

   use des_bc, only: dem_mi, dem_bcmi, dem_mi_time
   use des_rxns, only: des_x_s
   use des_thermo, only: des_t_s
   use discretelement
   use error_manager
   use geometry, only: no_k
   use mfix_pic, only: DES_STAT_WT, MPPIC
   use param, only: dimension_n_s
   use run, only: energy_eq, run_name
   use write_res1_des, only: init_write_res_des, finl_write_res_des
   use write_res1_des, only: write_res_carray, write_res_des, write_res_parray, write_res_parray_gsp
   use usr
CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
!
!  module name: des_write_restart
!  purpose: writing des data for restart
!
!  Author : Pradeep G
!  Purpose : Reads either single restart file or multiple restart files
!            (based on bdist_io) flag
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^c
   SUBROUTINE WRITE_RES0_DES
      use functions, only: is_nonexistent
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      INTEGER :: LC1
      INTEGER :: lNEXT_REC
      INTEGER :: lDIMN

      DOUBLE PRECISION :: VERSION

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FN_MAGTMP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MN_MAGTMP

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FTbondXTMP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FTbondYTMP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FTbondZTMP

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MTbondXTMP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MTbondYTMP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MTbondZTMP

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: n_histxTMP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: n_histyTMP
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: n_histzTMP

      INTEGER :: I, LL, cc, CC_START, CC_END, nbid, size1
      INTEGER, ALLOCATABLE :: PAIR_OWNER_GID(:),PAIR_NEIGH_GID(:)
! Set the version of the DES RES file.
      VERSION = 1.2
      IF(USE_JKR) VERSION = 2.2D0

! Set the output dimension.
      lDIMN = merge(2,3,NO_K)

      CALL INIT_WRITE_RES_DES(trim(RUN_NAME), VERSION, lNEXT_REC)

      CALL WRITE_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, TECPLOT_FINDEX)
      CALL WRITE_RES_DES(lNEXT_REC, DTSOLID)

      DO LC1 = 1, lDIMN
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_POS_NEW(:,LC1))
      ENDDO

      CALL WRITE_RES_pARRAY(lNEXT_REC, iGLOBAL_ID)

      CALL WRITE_RES_pARRAY(lNEXT_REC, particle_state)

      DO LC1 = 1, lDIMN
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_VEL_NEW(:,LC1))
      ENDDO

      DO LC1 = 1, merge(1,3,NO_K)
         CALL WRITE_RES_pARRAY(lNEXT_REC, OMEGA_NEW(:,LC1))
      ENDDO

      CALL WRITE_RES_pARRAY(lNEXT_REC, DES_RADIUS)
      CALL WRITE_RES_pARRAY(lNEXT_REC, RO_SOL)

      IF(MPPIC) &
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_STAT_WT)

      IF(CGDEM) THEN
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_CGP_STW)
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_CGP_RPR)
      ENDIF

      IF(ENERGY_EQ) &
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_T_s)

      CALL WRITE_RES_pARRAY(lNEXT_REC, PIJK(:,5))

      IF(allocated(DES_X_s)) THEN
         DO LC1=1, DIMENSION_N_S
            CALL WRITE_RES_pARRAY(lNEXT_REC, DES_X_s(:,LC1))
         ENDDO
      ENDIF

! DES User defined variable :: added for VERSION>= 1.1
      CALL WRITE_RES_DES(lNEXT_REC, DES_USR_VAR_SIZE)
      DO LC1=1,DES_USR_VAR_SIZE
         CALL WRITE_RES_pARRAY(lNEXT_REC, DES_USR_VAR(LC1,:))
      ENDDO
! Residence time
     CALL WRITE_RES_pARRAY(lNEXT_REC, RESIDENCE_TIME)
!SuperDEM
      IF(SuperDEM) THEN
         DO LC1=1,3
            CALL WRITE_RES_pARRAY(lNEXT_REC, super_r(:,LC1))
         ENDDO
         DO LC1=1,2
            CALL WRITE_RES_pARRAY(lNEXT_REC, super_mn(:,LC1))
         ENDDO
         DO LC1=1,4
            CALL WRITE_RES_pARRAY(lNEXT_REC, super_q(:,LC1))
         ENDDO
      ENDIF
!GluedSphereDEM
      IF(gsp_explicit) THEN
      ! @renjieke, in gsp implicit, glueSurface does not need a save, it will recompute
         CALL WRITE_RES_pARRAY(lNEXT_REC, gid_list)
         CALL WRITE_RES_pARRAY(lNEXT_REC, gp_kps)
         CALL WRITE_RES_pARRAY(lNEXT_REC, gp_sa)

         DO LC1=1, 6
           CALL WRITE_RES_pARRAY(lNEXT_REC, gp_neighsid(:, LC1))
         ENDDO
         DO LC1=1, 6
            CALL WRITE_RES_pARRAY(lNEXT_REC, gp_neighsa(:, LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY(lNEXT_REC, sc2gpc_vec(:,LC1))
         ENDDO

         if(flexible_gsp) then
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bfn(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bmn(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bftx(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bfty(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bftz(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bmtx(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bmty(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bmtz(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bnx(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bny(:, LC1))
            ENDDO
            DO LC1=1, 6
               CALL WRITE_RES_pARRAY(lNEXT_REC, bnz(:, LC1))
            ENDDO
         endif

         CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, particle_state_gsp)
         CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueMass)
	      CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueVolume)
         CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueDiameter)
         CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueBounding)
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueForce(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueTorque(:,LC1))
         ENDDO
         ! DO LC1=1, lDIMN
         !    CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueForcePP(:,LC1))
         ! ENDDO
         ! DO LC1=1, lDIMN
         !    CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueTorquePP(:,LC1))
         ! ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueAngMom(:,LC1))
         ENDDO
         DO LC1=1, 4
           CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueQuat(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueInertia(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueOmg(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
           CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, glueVel(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY_GSP(lNEXT_REC, gluePos(:,LC1))
         ENDDO
      ELSEIF(GSP_IMPLICIT) THEN
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY(lNEXT_REC, glueEX(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY(lNEXT_REC, glueEY(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY(lNEXT_REC, glueEZ(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL WRITE_RES_pARRAY(lNEXT_REC, glueAngMom(:,LC1))
         ENDDO
         DO LC1=1, 4
           CALL WRITE_RES_pARRAY(lNEXT_REC, glueQuat(:,LC1))
         ENDDO
      ENDIF

      IF(.NOT.MPPIC) THEN
         ! neighbor_index will be rebuilt after read, no need to save at this time.
         ! neighbors will be used to build the restart map, it is needed.
         CALL WRITE_RES_cARRAY(lNEXT_REC, NEIGHBORS(:), pLOC2GLB=.TRUE.)
      ENDIF

      DO LC1=1, lDIMN
         CALL WRITE_RES_cARRAY(lNEXT_REC,PFT_NEIGHBOR(LC1,:))

         IF (SuperDEM) THEN
              CALL WRITE_RES_cARRAY(lNEXT_REC,CONTACT_POINT_A(LC1,:))
              CALL WRITE_RES_cARRAY(lNEXT_REC,CONTACT_POINT_B(LC1,:))
              CALL WRITE_RES_cARRAY(lNEXT_REC,CONTACT_POINT_A_old(LC1,:))
              CALL WRITE_RES_cARRAY(lNEXT_REC,CONTACT_POINT_B_old(LC1,:))
         ENDIF
      ENDDO

      IF (SuperDEM) CALL WRITE_RES_cARRAY(lNEXT_REC,CONTACT_LAMBDA_A(:))
      IF (SuperDEM) CALL WRITE_RES_cARRAY(lNEXT_REC,CONTACT_LAMBDA_B(:))
      IF (SuperDEM) CALL WRITE_RES_cARRAY(lNEXT_REC,CONTACT_LAMBDA_A_OLD(:))
      IF (SuperDEM) CALL WRITE_RES_cARRAY(lNEXT_REC,CONTACT_LAMBDA_B_OLD(:))

! Version 2.1 appends spherical JKR histories after the official 1.2
! collision layout. Neighbor slots retain the official global-ID map.
      IF(USE_JKR) THEN
         ALLOCATE(PAIR_OWNER_GID(SIZE(NEIGHBORS)))
         ALLOCATE(PAIR_NEIGH_GID(SIZE(NEIGHBORS)))
         PAIR_OWNER_GID = 0
         PAIR_NEIGH_GID = 0
         DO LL=1,MAX_PIP
            IF(IS_NONEXISTENT(LL)) CYCLE
            CC_START = 1
            IF(LL > 1) CC_START = NEIGHBOR_INDEX(LL-1)
            CC_END = NEIGHBOR_INDEX(LL)
            DO CC=CC_START,CC_END-1
               IF(IS_NONEXISTENT(NEIGHBORS(CC))) CYCLE
               PAIR_OWNER_GID(CC) = IGLOBAL_ID(LL)
               PAIR_NEIGH_GID(CC) = IGLOBAL_ID(NEIGHBORS(CC))
            ENDDO
         ENDDO
         CALL WRITE_RES_cARRAY(lNEXT_REC,PAIR_OWNER_GID)
         CALL WRITE_RES_cARRAY(lNEXT_REC,PAIR_NEIGH_GID)
         DO LC1=1,lDIMN
            CALL WRITE_RES_cARRAY(lNEXT_REC,PFR_NEIGHBOR(LC1,:))
         ENDDO
         CALL WRITE_RES_cARRAY(lNEXT_REC,JKR_CONTACT_STATE(:))
         CALL WRITE_RES_DES(lNEXT_REC,COLLISION_ARRAY_MAX)
         DO CC=1,COLLISION_ARRAY_MAX
            CALL WRITE_RES_pARRAY(lNEXT_REC,WALL_COLLISION_FACET_ID(CC,:))
            CALL WRITE_RES_pARRAY(lNEXT_REC,JKR_WALL_CONTACT_STATE(CC,:))
            DO LC1=1,lDIMN
               CALL WRITE_RES_pARRAY(lNEXT_REC,WALL_COLLISION_PFT(LC1,CC,:))
               CALL WRITE_RES_pARRAY(lNEXT_REC,WALL_COLLISION_PFR(LC1,CC,:))
            ENDDO
         ENDDO
         CALL WRITE_RES_cARRAY(lNEXT_REC,PFTW_NEIGHBOR(:))
         DO CC=1,COLLISION_ARRAY_MAX
            CALL WRITE_RES_pARRAY(lNEXT_REC,WALL_COLLISION_PFTW(CC,:))
         ENDDO
! Contact geometry belongs to the last force evaluation, before the
! position/omega update represented by this checkpoint's particle state.
         DO LC1=1,4
            CALL WRITE_RES_cARRAY(lNEXT_REC,JKR_FRAME_NEIGHBOR(LC1,:))
            DO CC=1,COLLISION_ARRAY_MAX
               CALL WRITE_RES_pARRAY(lNEXT_REC,WALL_COLLISION_FRAME(LC1,CC,:))
            ENDDO
         ENDDO
         DEALLOCATE(PAIR_OWNER_GID)
         DEALLOCATE(PAIR_NEIGH_GID)
      ENDIF

      CALL WRITE_RES_DES(lNEXT_REC, DEM_BCMI)
      DO LC1=1, DEM_BCMI
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI_TIME(LC1))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%VACANCY)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%OCCUPANTS)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%WINDOW)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%OFFSET)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%L)
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%W(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%H(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%P(:))
         CALL WRITE_RES_DES(lNEXT_REC, DEM_MI(LC1)%Q(:))
      ENDDO

      !CALL WRITE_RES_DES(lNEXT_REC, small_particles_removed_10s_counter) ! Number of small particles removed last 10s based on threshold  remove_small_particle_dp
      CALL WRITE_RES_DES(lNEXT_REC, small_particles_removed) ! Total number of small particles removed during the simulation based on threshold remove_small_particle_dp
      !CALL WRITE_RES_DES(lNEXT_REC, report_small_particle_time) !Reporting time to report number of small particles removed last 10s

      CALL FINL_WRITE_RES_DES

      RETURN

   END SUBROUTINE WRITE_RES0_DES

END MODULE WRITE_RES0_DES_MOD
