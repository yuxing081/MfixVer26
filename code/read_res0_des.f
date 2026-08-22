MODULE READ_RES0_DES_MOD

   use des_allocate, only: allocate_dem_mi, particle_grow_GluedSphereDEM
   use des_bc, only: dem_mi, dem_bcmi, dem_mi_time
   use des_rxns, only: des_x_s
   use des_thermo, only: des_t_s
   use discretelement
   use error_manager
   use geometry, only: no_k
   use mfix_pic, only: MPPIC, DES_STAT_WT
   use param, only: dimension_n_s
   use read_res1_des, only: init_read_res_des, finl_read_res_des, read_par_col, read_par_pos
   use read_res1_des, only: read_res_carray, read_res_des, read_res_parray, read_res_parray_gsp
   use run, only: energy_eq, run_name, run_type, time, any_solids_species_eq

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DES_READ_RESTART                                        !
!  Purpose : Reads either single restart file or multiple restart      !
!  fles (based on bdist_io) flag.                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_RES0_DES
      use mpi_utility, only: global_all_max, bcast
      use functions, only: is_nonexistent
      implicit none

      INTEGER :: LC1, LC2, LL, II
      INTEGER :: lDIMN, lNEXT_REC

      INTEGER :: lVAR_SIZE
      DOUBLE PRECISION :: VERSION
      DOUBLE PRECISION :: DTSOLID_DUMMY

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

      INTEGER :: size1, nbid

      lDIMN = merge(2,3,NO_K)

      CALL INIT_READ_RES_DES(trim(RUN_NAME), VERSION, lNEXT_REC)

      CALL READ_RES_DES(lNEXT_REC, VTP_FINDEX)
      CALL READ_RES_DES(lNEXT_REC, TECPLOT_FINDEX)
! DTSOLID is not used anymore upon restart so that changes in particle
! properties (say spring stiffness) generate a new solids time step
      CALL READ_RES_DES(lNEXT_REC, DTSOLID_DUMMY)

! Position data is read and used to setup pARRAY reads.
      CALL READ_PAR_POS(lNEXT_REC)

      CALL READ_RES_pARRAY(lNEXT_REC, iGLOBAL_ID)

      CALL READ_RES_pARRAY(lNEXT_REC, particle_state)

      DO LC1 = 1, lDIMN
         CALL READ_RES_pARRAY(lNEXT_REC, DES_VEL_NEW(:,LC1))
      ENDDO

      DO LC1 = 1, merge(1,3,NO_K)
         CALL READ_RES_pARRAY(lNEXT_REC, OMEGA_NEW(:,LC1))
      ENDDO

      CALL READ_RES_pARRAY(lNEXT_REC, DES_RADIUS)
      CALL READ_RES_pARRAY(lNEXT_REC, RO_SOL)

      IF(MPPIC) CALL READ_RES_pARRAY(lNEXT_REC, DES_STAT_WT)
      IF(CGDEM) CALL READ_RES_pARRAY(lNEXT_REC, DES_CGP_STW)
      IF(CGDEM) CALL READ_RES_pARRAY(lNEXT_REC, DES_CGP_RPR)
      IF(ENERGY_EQ) CALL READ_RES_pARRAY(lNEXT_REC, DES_T_s)

      IF(VERSION >= 1.2) THEN
         CALL READ_RES_pARRAY(lNEXT_REC, PIJK(:,5))
         IF(allocated(DES_X_s)) THEN
            DO LC1=1, DIMENSION_N_S
               CALL READ_RES_pARRAY(lNEXT_REC, DES_X_s(:,LC1))
            ENDDO
         ENDIF
      ELSE
         IF(ANY_SOLIDS_SPECIES_EQ) THEN
            CALL READ_RES_pARRAY(lNEXT_REC, PIJK(:,5))
            DO LC1=1, DIMENSION_N_S
               CALL READ_RES_pARRAY(lNEXT_REC, DES_X_s(:,LC1))
            ENDDO
         ENDIF
      ENDIF

      IF(VERSION >= 1.1) THEN
         CALL READ_RES_DES(lNEXT_REC, lVAR_SIZE)
         DO LC1=1, lVAR_SIZE
            if(lVAR_SIZE <= DES_USR_VAR_SIZE) &
            CALL READ_RES_pARRAY(lNEXT_REC, DES_USR_VAR(LC1,:))
         ENDDO
      ENDIF

! Residence time
      CALL READ_RES_pARRAY(lNEXT_REC, RESIDENCE_TIME)

!SuperDEM
      IF(SuperDEM) THEN
          DO LC1 = 1, 3
              CALL READ_RES_pARRAY(lNEXT_REC, super_r(:,LC1))
          ENDDO
          DO LC1 = 1, 2
              CALL READ_RES_pARRAY(lNEXT_REC, super_mn(:,LC1))
          ENDDO
          DO LC1 = 1, 4
              CALL READ_RES_pARRAY(lNEXT_REC, super_q(:,LC1))
          ENDDO
      ENDIF

!GluedSphereDEM
      IF(gsp_explicit) THEN
         CALL READ_RES_pARRAY(lNEXT_REC, gid_list)
         CALL READ_RES_pARRAY(lNEXT_REC, gp_kps)
         CALL READ_RES_pARRAY(lNEXT_REC, gp_sa)
         DO LC1 = 1, 6
            CALL READ_RES_pARRAY(lNEXT_REC, gp_neighsid(:, LC1))
         ENDDO
         DO LC1 = 1, 6
            CALL READ_RES_pARRAY(lNEXT_REC, gp_neighsa(:, LC1))
         ENDDO
         DO LC1 = 1, lDIMN
            CALL READ_RES_pARRAY(lNEXT_REC, sc2gpc_vec(:,LC1))
         ENDDO

         IF(flexible_gsp) THEN
            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bfn(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bmn(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bftx(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bfty(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bftz(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bmtx(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bmty(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bmtz(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bnx(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bny(:, LC1))
            ENDDO

            DO LC1=1, 6
               CALL READ_RES_pARRAY(lNEXT_REC, bnz(:, LC1))
            ENDDO
         ENDIF

         ! @renjieke 5-13-2025 GSP explicit model relies on NGluedparticles in many places
         ! Careful when changing the initialized value of NGluedparticles, it will cause bugs in many places
         ! Only PE_IO reads NGluedparticles from .RES file, so it should broadcast to other ranks
         CALL BCAST(NGluedparticles)
         ! Other ranks should also increase their array sizes, respectively
         CALL PARTICLE_GROW_GluedSphereDEM(NGluedParticles)

         CALL READ_RES_pARRAY_GSP(lNEXT_REC, particle_state_gsp)
         CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueMass)
	      CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueVolume)
         CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueDiameter)
         CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueBounding)
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueForce(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueTorque(:,LC1))
         ENDDO
         ! DO LC1=1, lDIMN
         !    CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueForcePP(:,LC1))
         ! ENDDO
         ! DO LC1=1, lDIMN
         !    CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueTorquePP(:,LC1))
         ! ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueAngMom(:,LC1))
         ENDDO
         DO LC1=1, 4
            CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueQuat(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueInertia(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueOmg(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY_GSP(lNEXT_REC, glueVel(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY_GSP(lNEXT_REC, gluePos(:,LC1))
         ENDDO
      ELSEIF(gsp_implicit) THEN
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY(lNEXT_REC, glueEX(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY(lNEXT_REC, glueEY(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY(lNEXT_REC, glueEZ(:,LC1))
         ENDDO
         DO LC1=1, lDIMN
            CALL READ_RES_pARRAY(lNEXT_REC, glueAngMom(:,LC1))
         ENDDO
         DO LC1=1, 4
            CALL READ_RES_pARRAY(lNEXT_REC, glueQuat(:,LC1))
         ENDDO
      ENDIF

! RES2 does not need the collision of BC information.
      IF(RUN_TYPE == 'RESTART_2') RETURN

! Collision/neighbor data is read and used to setup cARRAY reads.
      IF(.NOT.MPPIC) THEN
! Read the NEIGHBORS array from the restart file.
         CALL READ_PAR_COL(lNEXT_REC)
         DO LC1=1, lDIMN
            CALL READ_RES_cARRAY(lNEXT_REC, PFT_NEIGHBOR(LC1,:))
              IF (SuperDEM) THEN
               CALL READ_RES_cARRAY(lNEXT_REC, CONTACT_POINT_A(LC1,:))
               CALL READ_RES_cARRAY(lNEXT_REC, CONTACT_POINT_B(LC1,:))
               CALL READ_RES_cARRAY(lNEXT_REC, CONTACT_POINT_A_old(LC1,:))
               CALL READ_RES_cARRAY(lNEXT_REC, CONTACT_POINT_B_old(LC1,:))
              ENDIF
         ENDDO
         IF(SuperDEM) THEN
             CALL READ_RES_cARRAY(lNEXT_REC, CONTACT_LAMBDA_A(:))
             CALL READ_RES_cARRAY(lNEXT_REC, CONTACT_LAMBDA_B(:))
             CALL READ_RES_cARRAY(lNEXT_REC, CONTACT_LAMBDA_A_old(:))
             CALL READ_RES_cARRAY(lNEXT_REC, CONTACT_LAMBDA_B_old(:))
         ENDIF
      ENDIF

! Save the number of BCMI's read from input file, then read the
! value from the restart file.
      CALL READ_RES_DES(lNEXT_REC, DEM_BCMI)

! Allocation is done here to ignore keyword changes during RES1.
      IF(DEM_BCMI > 0) CALL ALLOCATE_DEM_MI

! Only save the number of mass inflows for RESTART_1. This allows
! for mass inflows to be added/removed with RESTART_2.
! Todo: Prune entering/exiting flagged particles for RESTART_2.
      DO LC1=1, DEM_BCMI
         CALL READ_RES_DES(lNEXT_REC, DEM_MI_TIME(LC1))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%VACANCY)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%OCCUPANTS)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%WINDOW)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%OFFSET)
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%L)

         LC2 = DEM_MI(LC1)%OCCUPANTS

         allocate(DEM_MI(LC1)%W(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%W(:))
         allocate(DEM_MI(LC1)%H(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%H(:))
         allocate(DEM_MI(LC1)%P(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%P(:))
         allocate(DEM_MI(LC1)%Q(LC2))
         CALL READ_RES_DES(lNEXT_REC, DEM_MI(LC1)%Q(:))
      ENDDO

      !CALL READ_RES_DES(lNEXT_REC, small_particles_removed_10s_counter) ! Number of small particles removed last 10s based on threshold  remove_small_particle_dp
      CALL READ_RES_DES(lNEXT_REC, small_particles_removed) ! Total number of small particles removed during the simulation based on threshold remove_small_particle_dp
      !CALL READ_RES_DES(lNEXT_REC, report_small_particle_time) !Reporting time to report number of small particles removed last 10s


      CALL FINL_READ_RES_DES


      WRITE(ERR_MSG,"('DES restart file read at Time = ',g12.5)") TIME
      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

      RETURN

   END SUBROUTINE READ_RES0_DES

END MODULE READ_RES0_DES_MOD
