#include "error.inc"

MODULE vtp

   USE cdist, only: bdist_io
   USE compar, only: myPE, pe_io, numpes, mpierr
   USE des_rxns, only: des_x_s
   USE des_thermo, only: DES_T_s, des_c_ps
   USE desmpi, only: dprocbuf, drootbuf, iprocbuf, irootbuf, igath_sendcnt, igathercnts, idispls
   use des_rxns, only: NO_OF_DES_RXNS
   USE discretelement, only: DES_VEL_NEW, IGHOST_CNT, ORIENTATION, DES_USR_VAR_SIZE, DES_USR_VAR
   USE discretelement, only: MAX_PIP, PIP, DES_POS_NEW, PIJK, VTP_DIR, USE_COHESION, PARTICLE_ORIENTATION
   USE discretelement, only: POSTCOHESIVE, S_TIME, IGLOBAL_ID, DES_RADIUS, VTP_FINDEX, OMEGA_NEW, RO_SOL
   USE discretelement, only: PMASS, PVOL, F_GP
   USE discretelement, only: CGDEM, DES_CGP_RPR, DES_CGP_STW
   USE discretelement, only: super_r, super_mn, super_q, SuperDEM, SQP_POLY
   USE discretelement, only: pijk

   USE discretelement, only: GluedSphereDEM, gid_list, glueDiameter, glueBounding, gluePos, glueMass, glueVolume, glueQuat, NGluedParticles
   USE discretelement, only: gsp_explicit, gsp_implicit
   USE discretelement, only: glueEX, glueEY, glueEZ
   USE discretelement, only: FCHAINC,FCHAIN_MIDPOINT,FCHAIN_CONTACT_POINT
   USE discretelement, only: FCHAIN_ORIENT, FCHAIN_GLOBAL_ID1
   USE discretelement, only: FCHAIN_GLOBAL_ID1, FCHAIN_GLOBAL_ID2
   USE discretelement, only: FCHAIN_LOCAL_ID1, FCHAIN_LOCAL_ID2
   USE discretelement, only: FCHAIN_OVERLAP
   USE discretelement, only: FCHAIN_FN,FCHAIN_FT,FCHAIN_FC
   USE discretelement, only: FCHAIN_LENGTH,FCHAIN_FN_MAG, FCHAIN_FCMAX
   USE discretelement, only: VTP_LUB,VTP_P_DATA, VTP_FC_DATA, VTP_GLUE_DATA
   USE discretelement, only: FCHAIN_JKR_STATE,FCHAIN_JKR_RADIUS,FCHAIN_JKR_GAMMA
   USE discretelement, only: FCHAIN_TANGENTIAL_DISPLACEMENT
   USE discretelement, only: FCHAIN_ROLLING_DISPLACEMENT
   USE discretelement, only: FCHAIN_TWIST_ANGLE
   USE discretelement, only: RESIDENCE_TIME
   USE discretelement, only: COORDINATION, MIN_OVERLAP, MAX_OVERLAP, MEAN_OVERLAP
   USE discretelement, only: DES_COLL_FORCE
   USE discretelement, only: SuperDEM, sqp_a, sqp_b, sqp_c, sqp_m, sqp_n
   USE discretelement, only: sphere_radius_out, sphere_vel_out, sphere_rot_out, sphere_orient_out
   USE discretelement, only: sphere_usr_var_out, sphere_temp_out, sphere_cps_out, sphere_xs_out
   USE discretelement, only: sphere_prates_out, sphere_rho_out, sphere_mass_out, sphere_vol_out
   USE discretelement, only: sphere_rank_out, sphere_pid_out, sphere_id_out, sphere_res_time_out
   USE discretelement, only: sphere_coordination_out, sphere_min_overlap_out, sphere_max_overlap_out
   USE discretelement, only: sphere_mean_overlap_out
   USE discretelement, only: sphere_coll_force_out, sphere_drag_out
   USE discretelement, only: sphere_total_mass_out, sphere_bounding_out, sphere_eqv_dia_out
   USE error_manager, only: err_msg, loglevel_error, log_message, loglevel_status, ival, ivar
   USE fs_util, only: create_dir
   USE functions, only: is_normal, is_rigid_motion, is_nonexistent
   USE functions, only: is_any_ghost, is_entering_ghost, is_exiting_ghost, is_ghost
   USE funits, only: unit_vtp, unit_pvd, unit_pvtp, dmp_log, unit_log
   USE mfix_pic, only: MPPIC, des_stat_wt
   USE mpi_comm_des, only: desmpi_gatherv, des_gather
   USE mpi_utility, only: allgather_1i, global_sum, global_all_sum, global_all_max, bcast
   USE output, only: FULL_LOG
   USE param, only: dim_i, dim_j, dim_k, dimension_n_s, dim_m
   USE param1, only: UNDEFINED, ZERO
   USE parse, only: des_rxn_name
   USE physprop, only: d_p0
   USE pvd_mod, only: read_pvd_frames, open_pvd_file, update_and_close_pvd_file, write_pvd_frames
   USE run, only: RUN_TYPE, RUN_NAME, TIME, ENERGY_EQ, ANY_SOLIDS_SPECIES_EQ, SOLIDS_MODEL
   USE rxns, only: SPECIES_ALIAS_s
   USE vtk, only: BELONGS_TO_VTK_SUBDOMAIN
   USE vtk, only: BUFFER, VTU_OFFSET, VTU_FILENAME,VTP_FILENAME, VTK_HEADER_SIZE
   USE vtk, only: DIMENSION_VTK, FRAME, VTK_REGION, VTK_DEFINED
   USE vtk, only: END_REC, BELONGS_TO_VTK_SUBDOMAIN
   USE vtk, only: NUMBER_OF_POINTS, BUFFER, END_REC, VTU_OFFSET, BELONGS_TO_VTK_SUBDOMAIN
   USE vtk, only: RESET_FRAME_AT_TIME_ZERO, PVTU_FILENAME, BUFFER, END_REC
   USE vtk, only: TIME_DEPENDENT_FILENAME, STRIP, ADD_VTU_DIRECTORY, FIND_VTK_BASENAME
   USE vtk, only: VTK_DBG_FILE, VTU_DIR
   USE vtk, only: VTK_T_START, VTK_T_END
   USE vtk, only: VTK_DEFINED, VTK_DATA, FRAME
   USE vtk, only: VTK_FILEBASE
   USE vtk, only: VTK_NXS, VTK_NYS, VTK_NZS
   USE vtk, only: VTK_PART_ANGULAR_VEL, VTK_PART_ORIENTATION
   USE vtk, only: VTK_PART_COHESION, VTK_PART_RANK, VTK_PART_ID
   USE vtk, only: VTK_PART_VEL, VTK_PART_USR_VAR, VTK_PART_TEMP
   USE vtk, only: VTK_PART_PHYSICAL_DIAMETER, VTK_PART_CGP_STAT_WT
   USE vtk, only: VTK_PART_PHASE_ID, VTK_PART_X_S, VTK_PART_C_PS
   USE vtk, only: VTK_PHASE_FOR_DES_X_S, VTK_PART_DENSITY
   USE vtk, only: VTK_PART_MASS, VTK_PART_VOLUME, VTK_PART_DRAG_COEFF
   USE vtk, only: VTK_SLICE_TOL, VTK_SELECT_MODE, VTK_PART_PHASE
   USE vtk, only: VTK_PART_RESIDENCE_TIME
   USE vtk, only: VTK_PART_COORDINATION, VTK_PART_OVERLAP_STATS
   USE vtk, only: VTK_FC_JKR_STATE, VTK_FC_JKR_RADIUS
   USE vtk, only: VTK_FC_JKR_GAMMA, VTK_FC_TANG_DISP, VTK_FC_ROLL_DISP
   USE vtk, only: VTK_FC_TWIST_ANGLE
   USE vtk, only: VTK_X_E, VTK_X_W, VTK_Y_S, VTK_Y_N, VTK_Z_B, VTK_Z_T
   USE VTK, ONLY: SAVE_PART_RRATES, PART_RRATES_OUT, VTK_PART_RRATE
   USE VTK, ONLY: VTK_FCMAX_ONLY, VTK_PART_COLL_FORCE
   USE VTK, ONLY: VTK_PART_GLUEID, VTK_PART_GLUEDIAMETER, VTK_PART_GLUEVOLUME, VTK_PART_GLUEBOUNDING, VTK_PART_GLUEMASS, VTK_PART_GLUEQUAT
   USE VTK, ONLY: VTK_GSP_OUTPUT_TYPE
   USE, intrinsic :: iso_c_binding
   use SQ_ROTATION_MOD
   use usr

#ifdef MPI
   USE mpi, only: mpi_comm_world
#endif

   IMPLICIT NONE

   INTEGER, PRIVATE :: GLOBAL_CNT
   INTEGER, PRIVATE :: LOCAL_CNT

   INTEGER :: DES_UNIT = 2000

! formatted file name
   CHARACTER(LEN=511) :: FNAME_VTP

   INTERFACE VTP_WRITE_DATA
      MODULE PROCEDURE VTP_WRITE_DP1
      MODULE PROCEDURE VTP_WRITE_DP2
      MODULE PROCEDURE VTP_WRITE_I1
   END INTERFACE VTP_WRITE_DATA

CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_DP1                                            !
!                                                                      !
! Purpose: Collect and write 1D double precision arrays to the VTP     !
! file. This routine is designed to collect the data for parallel and  !
! serial runs. This routine also manages the distribted IO case.       !
!``````````````````````````````````````````````````````````````````````!
   SUBROUTINE VTP_WRITE_DP1(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      DOUBLE PRECISION, INTENT(in) :: DATA(:)

      INTEGER :: LC, PC

      IF(bDist_IO) THEN

         WRITE(DES_UNIT,1000) NAME

         PC = 1
         DO LC = 1, MAX_PIP
            IF(PC > PIP) EXIT
            IF(IS_NONEXISTENT(LC)) CYCLE
            PC = PC+1
            IF(IS_ANY_GHOST(LC)) CYCLE
            WRITE(DES_UNIT, 1001,ADVANCE="NO") real(DATA(LC))
         ENDDO
         WRITE(DES_UNIT,1002)

      ELSE

         allocate (dProcBuf(LOCAL_CNT) )
         allocate (dRootBuf(GLOBAL_CNT))

         CALL DES_GATHER(DATA)

         IF(myPE == PE_IO) THEN
            WRITE(DES_UNIT,1000) NAME
            DO LC=1, GLOBAL_CNT
               WRITE(DES_UNIT,1001,ADVANCE="NO") real(drootbuf(LC))
            ENDDO
            WRITE(DES_UNIT,1002)
         ENDIF

         deallocate(dProcBuf, dRootBuf)

      ENDIF

 1000 FORMAT('<DataArray type="Float32" Name="',A,'" format="ascii">')
 1001 FORMAT(ES14.6,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_DP1

!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_DP2                                            !
!                                                                      !
! Purpose: Collect and write 2D double precision arrays to the VTP     !
! file. This routine is designed to collect the data for parallel and  !
! serial runs. This routine also manages the distribted IO case.       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_DP2(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      DOUBLE PRECISION, INTENT(in) :: DATA(:,:)

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)

      CHARACTER(len=16) :: NOC
      INTEGER :: LB, UB, UBB
      INTEGER :: PC, LC1, LC2
      DOUBLE PRECISION :: qc(4), qrot(3,3),qrott(9),qrottt(9)

      UBB = 9
      if (NAME .EQ. 'Q_ROTATION')  THEN
          LB = LBOUND(DATA,2) ! LBOUND(DATA,1)
          UB = UBOUND(DATA,2) ! UBOUND(DATA,1) !
          NOC=''; WRITE(NOC,*) UBB
      ELSE
          LB = LBOUND(DATA,2)
          UB = UBOUND(DATA,2)
          NOC=''; WRITE(NOC,*) (UB-LB)+1
      ENDIF


      IF(bDist_IO) THEN


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! for superquradric rotation matrix output,UB=4
          IF(UB>=4) then
              WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))
              PC = 1
              DO LC1 = 1, MAX_PIP
                  IF(PC > PIP) EXIT
                  IF(IS_NONEXISTENT(LC1)) CYCLE
                  PC = PC+1
                  IF(IS_ANY_GHOST(LC1)) CYCLE
                  Qc(1)=DATA(LC1,1)  !DATA(1,LC1)
                  Qc(2)=DATA(LC1,2)  !DATA(2,LC1)
                  Qc(3)=DATA(LC1,3)  !DATA(3,LC1)
                  Qc(4)=DATA(LC1,4)  !DATA(4,LC1)
                  call  QROTATIONMATRIX(Qc,Qrot)

                  Qrott(1)=2.0D0*Qrot(1,1)
                  Qrott(2)=2.0D0*Qrot(1,2)
                  Qrott(3)=2.0D0*Qrot(1,3)
                  Qrott(4)=2.0D0*Qrot(2,1)
                  Qrott(5)=2.0D0*Qrot(2,2)
                  Qrott(6)=2.0D0*Qrot(2,3)
                  Qrott(7)=2.0D0*Qrot(3,1)
                  Qrott(8)=2.0D0*Qrot(3,2)
                  Qrott(9)=2.0D0*Qrot(3,3)

                  DO LC2=LB, UBB
                      WRITE(DES_UNIT,1001,ADVANCE="NO") real(QROTT(LC2))
                  ENDDO
              ENDDO
              WRITE(DES_UNIT,1002)
          ELSE
              WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))
              PC = 1
              DO LC1 = 1, MAX_PIP
                  IF(PC > PIP) EXIT
                  IF(IS_NONEXISTENT(LC1)) CYCLE
                  PC = PC+1
                  IF(IS_ANY_GHOST(LC1)) CYCLE
                  DO LC2=LB, UB
                      WRITE(DES_UNIT,1001,ADVANCE="NO") real(DATA(LC1,LC2))
                  ENDDO
              ENDDO
              WRITE(DES_UNIT,1002)
          ENDIF

      ELSE

          allocate (dProcBuf(LOCAL_CNT) )
          allocate (dRootBuf(GLOBAL_CNT))
          allocate (ltemp_array((UB-LB)+1,GLOBAL_CNT))

          DO LC1 = LB, UB
              CALL DES_GATHER(DATA(:,LC1))
              ltemp_array(LC1,:) = drootbuf(:)
          ENDDO

          IF(myPE == PE_IO) THEN

              IF(NAME .EQ. 'Q_ROTATION') then

                  WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))
                  DO LC1=1, GLOBAL_CNT

                      Qc(1)=ltemp_array(1,LC1)
                      Qc(2)=ltemp_array(2,LC1)
                      Qc(3)=ltemp_array(3,LC1)
                      Qc(4)=ltemp_array(4,LC1)

                      call  QROTATIONMATRIX(Qc,Qrot)

                      Qrott(1)=2.0D0*Qrot(1,1)
                      Qrott(2)=2.0D0*Qrot(1,2)
                      Qrott(3)=2.0D0*Qrot(1,3)
                      Qrott(4)=2.0D0*Qrot(2,1)
                      Qrott(5)=2.0D0*Qrot(2,2)
                      Qrott(6)=2.0D0*Qrot(2,3)
                      Qrott(7)=2.0D0*Qrot(3,1)
                      Qrott(8)=2.0D0*Qrot(3,2)
                      Qrott(9)=2.0D0*Qrot(3,3)

                      DO LC2=LB, UBB
                          WRITE(DES_UNIT,1001,ADVANCE="NO") real(QROTt(LC2))
                      ENDDO
                  ENDDO
                  WRITE(DES_UNIT,1002)
              ELSE
                  WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))
                  DO LC1=1, GLOBAL_CNT
                      DO LC2=LB, UB
                          WRITE(DES_UNIT,1001,ADVANCE="NO") &
                              real(ltemp_array(LC2,LC1))
                      ENDDO
                  ENDDO
                  WRITE(DES_UNIT,1002)
              ENDIF
          ENDIF
          deallocate (dProcBuf, dRootBuf, ltemp_array)

      ENDIF


 1000 FORMAT('<DataArray type="Float32" Name="',A,'" NumberOf',        &
         'Components="',A,'" format="ascii">')
 1001 FORMAT(ES14.6,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_DP2



!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_I1                                             !
!                                                                      !
! Purpose: Collect and write 1D integer arrays to the VTP file. This   !
! routine is designed to collect the data for parallel and serial      !
! runs. This routine also manages the distribted IO case.              !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_I1(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      INTEGER, INTENT(in) :: DATA(:)

      INTEGER :: LC, PC

      IF(bDist_IO) THEN

         WRITE(DES_UNIT,1000) NAME

         PC = 1
         DO LC = 1, MAX_PIP
            IF(PC > PIP) EXIT
            IF(IS_NONEXISTENT(LC)) CYCLE
            PC = PC+1
            IF(IS_ANY_GHOST(LC)) CYCLE
            WRITE(DES_UNIT, 1001,ADVANCE="NO") DATA(LC)
         ENDDO
         WRITE(DES_UNIT,1002)

      ELSE

         allocate (iProcBuf(LOCAL_CNT) )
         allocate (iRootBuf(GLOBAL_CNT))

         CALL DES_GATHER(DATA)

         IF(myPE == PE_IO) THEN
            WRITE(DES_UNIT,1000) NAME
            DO LC=1, GLOBAL_CNT
               WRITE(DES_UNIT,1001,ADVANCE="NO") irootbuf(LC)
            ENDDO
            WRITE(DES_UNIT,1002)
         ENDIF

         deallocate(iProcBuf, iRootBuf)

      ENDIF

 1000 FORMAT('<DataArray type="Float32" Name="',A,'" format="ascii">')
 1001 FORMAT(I10,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_I1


!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_ELEMENT                                        !
!                                                                      !
! Purpose: Write a string to the VTP file. It masks the need to check  !
! the logical before flushing.                                         !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_ELEMENT(ELEMENT)

      CHARACTER(len=*), INTENT(in) :: ELEMENT

      IF(bDist_IO .OR. myPE == PE_IO) &
         WRITE(DES_UNIT,"(A)") ELEMENT

      RETURN
      END SUBROUTINE VTP_WRITE_ELEMENT


!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_OPEN_FILE                                            !
!                                                                      !
! Purpose: This routine opens the VTP file and calculates the offsets  !
! for dmp data collection.                                             !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_OPEN_FILE(NoPc)

      IMPLICIT NONE

      CHARACTER(len=*) :: NoPc

      INTEGER :: NumberOfPoints

! Variables related to gather
      integer lgathercnts(0:numpes-1), lproc

! check whether an error occurs in opening a file
      INTEGER :: IOS
! Integer error flag.
      INTEGER :: IER

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS_VTP
! status of the vtp file to be written
      CHARACTER(LEN=8) :: STATUS_VTP

      IF(TRIM(VTP_DIR)/='.'.AND.TRIM(VTP_DIR)/='' .AND. myPE == PE_IO) THEN
         CALL CREATE_DIR(trim(VTP_DIR))
      END IF

! Initial the global count.
      GLOBAL_CNT = 10
! Calculate the number of 'real' particles on the local process.
      LOCAL_CNT = PIP - iGHOST_CNT

! Distributed IO
      IF(bDIST_IO) THEN
         NumberOfPoints = LOCAL_CNT
         WRITE(NoPc,"(I10.10)") NumberOfPoints

         IF(TRIM(VTP_DIR)/='.'.AND.TRIM(VTP_DIR)/='') THEN
            WRITE(fname_vtp,'(A,"/",A,"_DES",I4.4,"_",I5.5,".vtp")') &
               trim(VTP_DIR), trim(run_name), vtp_findex, mype
         ELSE
            WRITE(fname_vtp,'(A,"_DES",I4.4,"_",I5.5,".vtp")') &
               trim(run_name), vtp_findex, mype
         ENDIF

! Serial IO
      ELSE

! Calculate the total number of particles system-wide.
         call global_sum(LOCAL_CNT, GLOBAL_CNT)
         NumberOfPoints = GLOBAL_CNT
         WRITE(NoPc,"(I10.10)") NumberOfPoints

! Set the send count from the local process.
         igath_sendcnt = LOCAL_CNT

! Collect the number of particles on each rank.all ranks.
         lgathercnts = 0
         lgathercnts(myPE) = LOCAL_CNT
         call global_sum(lgathercnts,igathercnts)

! Calculate the rank displacements.
         idispls(0) = 0
         DO lPROC = 1,NUMPEs-1
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
         ENDDO

! set the file name and unit number and open file
         IF(TRIM(VTP_DIR)/='.'.AND.TRIM(VTP_DIR)/='') THEN
            WRITE(fname_vtp,'(A,"/",A,"_DES_",I5.5,".vtp")') &
               trim(VTP_DIR),trim(run_name), vtp_findex
         ELSE
            WRITE(fname_vtp,'(A,"_DES_",I5.5,".vtp")') &
               trim(run_name), vtp_findex
         ENDIF
      ENDIF

      IER = 0
      IF(bDIST_IO .OR. myPE == PE_IO) THEN

! The file should be new but could exist due to restarting.
         STATUS_VTP = 'NEW'
! Check to see if the file already exists.
         INQUIRE(FILE=FNAME_VTP,EXIST=EXISTS_VTP)
! The given file should not exist if the run type is NEW.
         IF(EXISTS_VTP)THEN
! The VTP should never exist for a NEW run.
            IF(RUN_TYPE == 'NEW')THEN
               IER = 1
! The file may exist during a RESTART.
            ELSE
               STATUS_VTP = 'REPLACE'
            ENDIF
         ENDIF

! Open the file and record any errors.
         IF(IER == 0) THEN
            OPEN(UNIT=DES_UNIT, FILE=FNAME_VTP,   &
               STATUS=STATUS_VTP, IOSTAT=IOS)
            IF(IOS /= 0) IER = 2
         ENDIF
      ENDIF

      CALL GLOBAL_ALL_MAX(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG, 1100) IER
         CALL LOG_ERROR()
      ENDIF

 1100 FORMAT('Error 1100: Unable to open VTP file. This could be ',    &
         'caused by a VTP',/'file with the same file name already ',   &
         'existing, or an error code',/' returned by the OPEN ',       &
         'function.'/'Error code: ',I2,4x,'Aborting.')

      END SUBROUTINE VTP_OPEN_FILE

!......................................................................!
! SUBROUTINE: VTP_CLOSE_FILE                                           !
!                                                                      !
! Purpose: This routine closes the vtp file.                           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_CLOSE_FILE

      VTP_FINDEX=VTP_FINDEX+1

      IF(bDist_io .OR. (myPE .eq.pe_IO)) CLOSE(des_unit)

      END SUBROUTINE VTP_CLOSE_FILE


!......................................................................!
! SUBROUTINE: ADD_VTP_TO_PVD                                           !
!                                                                      !
! Purpose: This routine opens the pvd file.                            !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE ADD_VTP_TO_PVD

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Index position of desired character
      INTEGER IDX_f, IDX_b
! logical used for testing is the data file already exists
      LOGICAL :: EXISTS_PVD
! Generic input limited to 256 characters
      CHARACTER(LEN=256) INPUT

! formatted file name
      CHARACTER(LEN=64) :: FNAME_PVD = ''
! formatted time
      CHARACTER(LEN=64) :: cTIME = ''

      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

! IO Status flag
      INTEGER :: IOS

! Variables related to gather
      integer :: IER

!-----------------------------------------------

! Initialize the error flag.
      IER = 0

! Obtain the file name and open the pvd file
      FNAME_PVD = TRIM(RUN_NAME)//'_DES.pvd'

! The PVD file is only written by PE_IO with serial IO.
      IF(myPE == PE_IO .AND. .NOT.bDist_IO) THEN

! Check to see if the file already exists.
         INQUIRE(FILE=FNAME_PVD,EXIST=EXISTS_PVD)

         IF(FIRST_PASS) THEN

! Open the "NEW" file and write the necessary header information.
            IF(RUN_TYPE /= 'RESTART_1')THEN

! The file exists but first_pass is also true so most likely an existing
! file from an earlier/other run is present in the directory. Exit to
! prevent accidentally overwriting the existing file.
               IF(EXISTS_PVD) THEN
                  IER = 1
               ELSE
                  OPEN(UNIT=UNIT_PVD,FILE=FNAME_PVD,STATUS='NEW')
                  WRITE(UNIT_PVD,"(A)")'<?xml version="1.0"?>'
                  WRITE(UNIT_PVD,"(A)")'<VTKFile type="Collection" &
                     &version="0.1" byte_order="BigEndian">'
                  WRITE(UNIT_PVD,"(3X,'<Collection>')")
               ENDIF

! This is the first pass of a restart run. Extra care is needed to make
! sure that the pvd file is ready to accept new data.
            ELSE ! a restart run
               IF(EXISTS_PVD) THEN
! Open the file at the beginning.
                  OPEN(UNIT=UNIT_PVD,FILE=FNAME_PVD,&
                     POSITION="REWIND",STATUS='OLD',IOSTAT=IOS)
                  IF(IOS /= 0) IER = 2
               ELSE ! a pvd file does not exist
                  IER = 3
               ENDIF

               IF(IER == 0) THEN
! Loop over the entries in the PVD file, looking for a match to the
! file that is being written. If no match is found, the data will be
! appended to the end of the pvd file, otherwise, the old data will
! be over-written.
                  DO
! Read in the entries of the PVD file.
                     READ(UNIT_PVD,"(A)",IOSTAT=IOS)INPUT
                     IF(IOS > 0) THEN
                        IER = 4
                        EXIT
                     ELSEIF(IOS<0)THEN
! The end of the pvd file has been reached without finding an entry
! matching the current record. Exit the loop.
                        BACKSPACE(UNIT_PVD)
                        BACKSPACE(UNIT_PVD)
                        BACKSPACE(UNIT_PVD)
                        EXIT
                     ENDIF
! Find the first instances of file=" and "/> in the read data.
                     IDX_f = INDEX(INPUT,'file="')
                     IDX_b = INDEX(INPUT,'"/>')
! Skip rows that do not contain file data
                     IF(IDX_f == 0 .AND. IDX_b == 0) CYCLE
! Truncate the file name from the read data
                     WRITE (INPUT,"(A)") INPUT(IDX_f+6:IDX_b-1)
! If the file name matches the current VTP record, break the loop to
! over-write this record.
                     IF(TRIM(FNAME_VTP) == TRIM(INPUT)) THEN
                        BACKSPACE(UNIT_PVD)
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF ! No errors
            ENDIF ! run_type new or restart

         ELSE ! not FIRST_PASS
            OPEN(UNIT=UNIT_PVD,FILE=FNAME_PVD,&
               POSITION="APPEND",STATUS='OLD',IOSTAT=IOS)
            IF (IOS /= 0) IER = 2
         ENDIF

      ENDIF ! if myPE == PE_IO and not distributed IO


      CAlL GLOBAL_ALL_SUM(IER)
      IF(IER /= 0) THEN
         SELECT CASE(IER)
         CASE(1); WRITE(ERR_MSG,1101) trim(FNAME_PVD)
         CASE(2); WRITE(ERR_MSG,1102) trim(FNAME_PVD)
         CASE(3); WRITE(ERR_MSG,1103) trim(FNAME_PVD)
         CASE(4); WRITE(ERR_MSG,1104) trim(FNAME_PVD)
         CASE DEFAULT; WRITE(ERR_MSG,1105) trim(FNAME_PVD)
         END SELECT
         CALL LOG_ERROR()
      ENDIF

 1101 FORMAT('Error 1101: A PVD file was detected in the run ',        &
         'directory which should',/'not exist for a NEW run.',/        &
         'File: ',A)

 1102 FORMAT('Error 1102: Fatal error status returned while OPENING ', &
         'PVD file.',/'File: ', A)

 1103 FORMAT('Error 1103: PVD file MISSING from run directory.',/      &
         'File: ',A)

 1104 FORMAT('Error 1104: Fatal error status returned while READING ', &
         'PVD file.',/'File: ', A)

 1105 FORMAT('Error 1105:: Fatal unclassified error when processing ', &
         'PVD file.',/'File: ', A)


! If there were no errors, updated the file.
      IF(myPE == PE_IO .AND. .NOT.bDist_IO) THEN

! Remove the last two lines written so that additional data can be added
         IF(.NOT.FIRST_PASS) THEN
            BACKSPACE(UNIT_PVD)
            BACKSPACE(UNIT_PVD)
         ENDIF

         WRITE(cTIME,"(F12.6)") S_TIME
! Write the data to the file
         WRITE(UNIT_PVD,"(6X,A,A,A,A,A,A,A)")&
         '<DataSet timestep="',trim(adjustl(cTIME)),'" ',&
         'group="" part="0" ',& ! necessary file data
         'file="',TRIM(FNAME_VTP),'"/>' ! file name of vtp

! Write the closing tags
         WRITE(UNIT_PVD,"(3X,A)")'</Collection>'
         WRITE(UNIT_PVD,"(A)")'</VTKFile>'

         CLOSE(UNIT_PVD)
      ENDIF
! Identify that the files has been created and opened for next pass
      FIRST_PASS = .FALSE.

! Return to the calling routine
      RETURN

      END SUBROUTINE ADD_VTP_TO_PVD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VTP_FILE                                         C
!  Purpose: Writes particles data in VTK format (Polydata VTP)         C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_VTP_FILE(LCV,MODE,TIMESTAMP)
      use discretelement, only: sphereNumber
      use constant, only: pi
      use time_mod, only: start_io, end_io
      IMPLICIT NONE
      INTEGER :: M,N,LCV
      INTEGER(c_int64_t) :: L

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)

      DOUBLE PRECISION :: TIMESTAMP

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PARTICLE_RANK
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PARTICLE_PHASE_ID
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FCHAIN_FN_MAG_ND
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SQP_POLY_SCALE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OVERLAP_ND
      DOUBLE PRECISION :: SUM_FN,GLOBAL_SUM_FN,AVG_FN
      INTEGER :: GLOBAL_FCHAINC
      INTEGER :: cur_n_sphere, total_n_sphere, lcurpar, MM
      INTEGER :: TMP_LOCAL_CNT, TMP_GLOBAL_CNT
      INTEGER :: lGatherCnts(0:NUMPEs-1)
      INTEGER :: lproc

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: des_bounding
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gsp_equivalent_diameter
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: des_total_mass
      INTEGER :: LL


      DOUBLE PRECISION, PARAMETER  :: THIRD = (1.0d0/3.0d0)

      VTK_REGION = LCV
! There is nothing to write if we are not in a defined vtk region
      IF(.NOT.VTK_DEFINED(VTK_REGION)) RETURN
      IF(TIMESTAMP<VTK_T_START(VTK_REGION)) RETURN
      IF(TIMESTAMP>VTK_T_END(VTK_REGION)) RETURN

      CALL START_IO

      VTP_P_DATA = .FALSE.
      VTP_FC_DATA = .FALSE.
      VTP_GLUE_DATA = .FALSE.

      SELECT CASE (VTK_DATA(VTK_REGION))

      CASE('P') ! Regular Particle data
         VTP_LUB = MAX_PIP
         VTP_P_DATA = .TRUE.
      CASE('F') ! Force chain data
         VTP_LUB = FCHAINC
         VTP_FC_DATA =.TRUE.
      CASE('N') ! Glued particle data
         VTP_LUB = NGluedParticles
	      VTP_GLUE_DATA = .TRUE.
      CASE DEFAULT
         CALL END_IO
!        print*,'should not be here (vtk_data)',VTK_REGION,VTK_DATA(VTK_REGION)
         RETURN

      END SELECT


! Regular particle data
      IF(VTK_DATA(VTK_REGION)/='P'.AND.VTK_DATA(VTK_REGION)/='F'.AND.VTK_DATA(VTK_REGION)/='N') THEN
          CALL END_IO
          RETURN
      END IF

      IF(MODE==0.AND.(VTK_DBG_FILE(VTK_REGION))) THEN
          CALL END_IO
          RETURN
      END IF

      IF(MODE==1.AND.(.NOT.VTK_DBG_FILE(VTK_REGION))) THEN
          CALL END_IO
          RETURN
      END IF

      CALL SETUP_VTK_REGION_PARTICLES

      CALL OPEN_VTP_FILE_BIN(MODE, TIMESTAMP)

! Only open pvd file when there are particles in vtk region
      IF(GLOBAL_CNT>0.AND.MODE==0) CALL OPEN_PVD_FILE

! First thing to check if NSphereGSP changed in current PE
      cur_n_sphere = 0
      IF(GSP_IMPLICIT .and. VTK_GSP_OUTPUT_TYPE(VTK_REGION) == 1) THEN
         DO lcurpar = 1, MAX_PIP
            !@renjieke, mass inflow particles are first set to entering particles
            ! they can not be written into vtp so use is_normal to cycle them
            IF(.NOT.(IS_NORMAL(lcurpar).OR.IS_RIGID_MOTION(lcurpar))) CYCLE
            MM = PIJK(lcurpar,5)
            cur_n_sphere = cur_n_sphere + sphereNumber(MM)
         ENDDO
         total_n_sphere = 0
         call global_all_sum(cur_n_sphere, total_n_sphere)

! pe 94 is the place where we have problem, it has a send cnt as 9 but it stuck somewhere
! need to find out why it stuck

! Set the send count from the local process.
         igath_sendcnt = cur_n_sphere
! Collect the number of particles on each rank.all ranks.
         lgathercnts = 0
         lgathercnts(myPE) = cur_n_sphere
         call global_sum(lgathercnts,igathercnts)
! Calculate the rank displacements.
         idispls(0) = 0
         DO lproc = 1,NUMPEs-1
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
         ENDDO
      ENDIF

! First pass write the data header.
! Second pass writes the data (appended binary format).

      DO PASS=WRITE_HEADER,WRITE_DATA

! @renjieke, write the x,y,z coordinates through this subroutine
         IF(GSP_IMPLICIT .and. VTK_GSP_OUTPUT_TYPE(VTK_REGION) == 1) THEN
            CALL WRITE_GEOMETRY_IN_VTP_BIN(PASS,cur_n_sphere,total_n_sphere)
         ELSE
            CALL WRITE_GEOMETRY_IN_VTP_BIN(PASS)
         ENDIF

         IF(VTP_P_DATA) THEN ! Regular particle data

! For project2, there are no components sphere save explicitly, so I have to save
! a component sphere list here, but I don't want them to be dynamic allocated when writing
! vtp, each time allocate and de allocate will cause a lot of trouble right
! @renjieke，make the program run first then go back to consider how to process the write vtp
            IF(GSP_IMPLICIT .and. VTK_GSP_OUTPUT_TYPE(VTK_REGION) == 1) THEN
               CALL WRITE_COMP_SPHERE_IMPLICIT(PASS, cur_n_sphere, total_n_sphere)
! Always save particle diameter, the GUI needs it to scale particles in
! the visualization viewport
            ELSE
               IF(ALLOCATED(DES_RADIUS)) THEN
                  IF(MPPIC.AND.ALLOCATED(DES_STAT_WT)) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN('Diameter',2.0D0*DES_RADIUS*DES_STAT_WT**(THIRD),PASS)
                  ELSE
                     IF(GSP_IMPLICIT .and. VTK_GSP_OUTPUT_TYPE(VTK_REGION) == 2) THEN
                     ! @Renjieke show particles in the bounding diameter in vtp file when using the implicit method with output type = 2
                        CALL WRITE_SCALAR_IN_VTP_BIN('Diameter',2.0D0*DES_RADIUS,PASS)
                     ELSE
                        CALL WRITE_SCALAR_IN_VTP_BIN('Diameter',2.0D0*DES_RADIUS,PASS)
                     ENDIF
                  ENDIF

! Dump whole gsp bounding diameter to each component spheres
                  IF(gsp_explicit .and. VTK_PART_GLUEBOUNDING(VTK_REGION)) THEN
                     ALLOCATE(des_bounding(size(des_radius)))
                     des_bounding = 0.0
                     DO LL = 1, size(des_radius)
                        des_bounding(LL) = glueBounding(gid_list(LL))
                     ENDDO
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Bounding_Diameter',des_bounding,PASS)
                     deallocate(des_bounding)
                  ENDIF

                  IF(gsp_explicit .and. VTK_PART_GLUEDIAMETER(VTK_REGION)) THEN
                     ALLOCATE(gsp_equivalent_diameter(size(des_radius)))
                     gsp_equivalent_diameter = 0.0
                     DO LL = 1, size(des_radius)
                        gsp_equivalent_diameter(LL) = glueDiameter(gid_list(LL))
                     ENDDO
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Equivalent_Diameter',gsp_equivalent_diameter,PASS)
                     deallocate(gsp_equivalent_diameter)
                  ENDIF

                  IF(gsp_implicit .and. VTK_PART_GLUEBOUNDING(VTK_REGION) .and. VTK_GSP_OUTPUT_TYPE(VTK_REGION) == 2) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Bounding_Diameter',glueBounding,PASS)
                  ENDIF

                  IF(gsp_implicit .and. VTK_PART_GLUEDIAMETER(VTK_REGION) .and. VTK_GSP_OUTPUT_TYPE(VTK_REGION) == 2) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Equivalent_Diameter',glueDiameter,PASS)
                  ENDIF

               ENDIF

               IF(VTK_PART_VEL(VTK_REGION).AND.ALLOCATED(DES_VEL_NEW)) &
                  CALL WRITE_VECTOR_IN_VTP_BIN('Velocity',DES_VEL_NEW,PASS)

               IF(CGDEM) THEN
                  IF(VTK_PART_PHYSICAL_DIAMETER(VTK_REGION).AND.ALLOCATED(DES_CGP_RPR)) &
                     CALL WRITE_SCALAR_IN_VTP_BIN('Physical_Diameter',2.0D0*DES_CGP_RPR,PASS)
                  IF(VTK_PART_CGP_STAT_WT(VTK_REGION).AND.ALLOCATED(DES_CGP_STW)) &
                     CALL WRITE_SCALAR_IN_VTP_BIN('Statistical_Weight',DES_CGP_STW,PASS)
               ENDIF

               IF(VTK_PART_ANGULAR_VEL(VTK_REGION).AND.ALLOCATED(OMEGA_NEW)) &
                  CALL WRITE_VECTOR_IN_VTP_BIN('Angular_velocity', OMEGA_NEW,PASS)

               IF(PARTICLE_ORIENTATION) THEN
                  IF(VTK_PART_ORIENTATION(VTK_REGION).AND.ALLOCATED(ORIENTATION)) &
                     CALL WRITE_VECTOR_IN_VTP_BIN('Orientation', ORIENTATION,PASS)
               ENDIF

               DO N=1, DES_USR_VAR_SIZE
                  IF(VTK_PART_USR_VAR(VTK_REGION,N).AND.ALLOCATED(DES_USR_VAR)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('User Defined Var '//trim(iVal(N)),DES_USR_VAR(N,:),PASS)
               ENDDO

               IF(ENERGY_EQ.AND.VTK_PART_TEMP(VTK_REGION).AND.ALLOCATED(DES_T_s)) &
               CALL WRITE_SCALAR_IN_VTP_BIN('Temperature', DES_T_s,PASS)

               IF(ENERGY_EQ.AND.VTK_PART_C_Ps(VTK_REGION).AND.ALLOCATED(DES_C_Ps)) &
               CALL WRITE_SCALAR_IN_VTP_BIN('Specific_heat', DES_C_Ps,PASS)

               IF(ANY_SOLIDS_SPECIES_EQ) THEN
                  DO N=1, DIMENSION_N_S
                     IF(VTK_PART_X_s(VTK_REGION,N).AND.ALLOCATED(DES_X_s)) THEN
                        M= VTK_PHASE_FOR_DES_X_S(VTK_REGION)
                        IF(M/=0) THEN
                           CALL WRITE_SCALAR_IN_VTP_BIN(trim(SPECIES_ALIAS_s(M,N)), DES_X_s(:,N),PASS)
                        ELSE
                           CALL WRITE_SCALAR_IN_VTP_BIN(trim(iVar('X_s',N)), DES_X_s(:,N),PASS)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF

               IF(SAVE_PART_RRATES) THEN
                  DO N=1, NO_OF_DES_RXNS
                     IF(VTK_PART_RRATE(VTK_REGION,N).AND.ALLOCATED(PART_RRATES_OUT)) THEN
                        CALL WRITE_SCALAR_IN_VTP_BIN('RRate_'//trim(DES_RXN_NAME(N)), PART_RRATES_OUT(:,N),PASS)
                     ENDIF
                  ENDDO
               ENDIF

               IF(VTK_PART_DENSITY(VTK_REGION).AND.ALLOCATED(RO_SOL)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('Density', RO_SOL,PASS)

               IF(VTK_PART_MASS(VTK_REGION).AND.ALLOCATED(PMASS)) THEN
                  CALL WRITE_SCALAR_IN_VTP_BIN('Mass', PMASS,PASS)
! Dump total gsp mass to each component spheres
                  IF(ALLOCATED(glueMass) .and. gsp_explicit .and. VTK_PART_GLUEMASS(VTK_REGION)) THEN
                     ALLOCATE(des_total_mass(size(PMASS)))
                     des_total_mass = 0.0
                     DO LL = 1, size(PMASS)
                        des_total_mass(LL) = glueMass(gid_list(LL))
                     ENDDO
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Mass',des_total_mass,PASS)
                     deallocate(des_total_mass)
                  ENDIF
               ENDIF

               IF(VTK_PART_VOLUME(VTK_REGION).AND.ALLOCATED(PVOL)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('Volume', PVOL,PASS)

               IF(USE_COHESION.AND.VTK_PART_COHESION(VTK_REGION).AND.ALLOCATED(PostCohesive)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('CohesiveForce', PostCohesive,PASS)

               IF(VTK_PART_COLL_FORCE(VTK_REGION).AND.ALLOCATED(DES_COLL_FORCE)) &
                  CALL WRITE_VECTOR_IN_VTP_BIN('Collision_force', DES_COLL_FORCE,PASS)

               IF(VTK_PART_DRAG_COEFF(VTK_REGION).AND.ALLOCATED(F_GP)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('Drag_coefficient', F_GP,PASS)

               IF(VTK_PART_RANK(VTK_REGION)) THEN
                  IF(PASS==WRITE_DATA) THEN
                     ALLOCATE(PARTICLE_RANK(MAX_PIP))
                     DO L = 1, MAX_PIP
                        PARTICLE_RANK(L) = DBLE(MyPE)
                     ENDDO
                  ENDIF

                  CALL WRITE_SCALAR_IN_VTP_BIN('Particle_Rank', PARTICLE_RANK,PASS)
                  IF(PASS==WRITE_DATA) DEALLOCATE(PARTICLE_RANK)
               ENDIF

               IF(VTK_PART_PHASE_ID(VTK_REGION)) THEN
                  IF(PASS==WRITE_DATA) THEN
                     ALLOCATE(PARTICLE_PHASE_ID(MAX_PIP))
                     DO L = 1, MAX_PIP
                        PARTICLE_PHASE_ID(L) = DBLE(PIJK(L,5))
                     ENDDO
                  ENDIF

                  CALL WRITE_SCALAR_IN_VTP_BIN('Particle_Phase_ID', PARTICLE_PHASE_ID,PASS)
                  IF(PASS==WRITE_DATA) DEALLOCATE(PARTICLE_PHASE_ID)
               ENDIF

               IF(VTK_PART_ID(VTK_REGION).AND.ALLOCATED(iGLOBAL_ID)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('Particle_ID', DBLE(iGLOBAL_ID),PASS)

               IF(VTK_PART_RESIDENCE_TIME(VTK_REGION).AND.ALLOCATED(RESIDENCE_TIME)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('Residence_Time', RESIDENCE_TIME,PASS)

               IF(VTK_PART_COORDINATION(VTK_REGION).AND.ALLOCATED(COORDINATION)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('Coordination', COORDINATION,PASS)

               IF(VTK_PART_OVERLAP_STATS(VTK_REGION).AND.ALLOCATED(MIN_OVERLAP)) THEN
                  CALL WRITE_SCALAR_IN_VTP_BIN('Min_overlap', MIN_OVERLAP,PASS)
                  CALL WRITE_SCALAR_IN_VTP_BIN('Max_overlap', MAX_OVERLAP,PASS)
                  CALL WRITE_SCALAR_IN_VTP_BIN('Mean_overlap', MEAN_OVERLAP,PASS)

   ! Non-dimensionalized overlap (by diameter)
                  IF(PASS==WRITE_DATA) THEN
                     WHERE(DES_RADIUS>ZERO) MIN_OVERLAP = MIN_OVERLAP/(2.0D0*DES_RADIUS)
                     WHERE(DES_RADIUS>ZERO) MAX_OVERLAP = MAX_OVERLAP/(2.0D0*DES_RADIUS)
                     WHERE(DES_RADIUS>ZERO) MEAN_OVERLAP = MEAN_OVERLAP/(2.0D0*DES_RADIUS)
                  ENDIF

                  CALL WRITE_SCALAR_IN_VTP_BIN('Min_overlap/Diameter', MIN_OVERLAP, PASS)
                  CALL WRITE_SCALAR_IN_VTP_BIN('Max_overlap/Diameter', MAX_OVERLAP,PASS)
                  CALL WRITE_SCALAR_IN_VTP_BIN('Mean_overlap/Diameter', MEAN_OVERLAP,PASS)
               ENDIF

               ! SuperDEM
               IF(SuperDEM) THEN
   ! Always save quaternions, needed by the GUI to visualize SQP particles
                  IF(ALLOCATED(super_q)) THEN
                     CALL WRITE_VECTOR_IN_VTP_BIN9('Q_ROTATION', super_q(:,1:4),PASS)
                     CALL WRITE_QUATERNION_IN_VTP_BIN('QUATERNION', super_q(:,1:4),PASS)
                  ENDIF
                  ! IF(SQP_POLY.AND.ALLOCATED(DES_RADIUS)) THEN
                  IF(ALLOCATED(DES_RADIUS)) THEN ! Save scale factor. Useful for
                                                ! polydisperse systems and reacting flows with
                                                ! constant solids density (variable size)
                     IF(PASS==WRITE_DATA) THEN
                        ALLOCATE(SQP_POLY_SCALE(MAX_PIP))
                        SQP_POLY_SCALE(:) = 1.0D0
                        DO L = 1, MAX_PIP
                           M = PIJK(L,5)
                           IF(M>0) SQP_POLY_SCALE(L) = 2.0D0*DES_RADIUS(L)/D_P0(M)
                        ENDDO
                     ENDIF

                     CALL WRITE_SCALAR_IN_VTP_BIN('SQP_POLY_SCALE', SQP_POLY_SCALE,PASS)
                     IF(PASS==WRITE_DATA) DEALLOCATE(SQP_POLY_SCALE)
                  ENDIF
               ENDIF

               !GluedSphere
               IF(gsp_explicit) THEN
                  IF(VTK_PART_GLUEID(VTK_REGION) .AND. ALLOCATED(gid_list)) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_ID', DBLE(gid_list) ,PASS)
                  ENDIF
               ENDIF

            ENDIF

         ELSE ! Force chain data or glued particle data

	         IF ( .NOT. VTP_GLUE_DATA) THEN ! Force chain data

               CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_LENGTH',FCHAIN_LENGTH,PASS)
               CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_FN_MAG',FCHAIN_FN_MAG,PASS)
               CALL WRITE_VECTOR_IN_VTP_BIN('FORCE_CHAIN_ORIENTATION',FCHAIN_ORIENT,PASS)
               CALL WRITE_VECTOR_IN_VTP_BIN('FORCE_CHAIN_FN',FCHAIN_FN,PASS)
               CALL WRITE_VECTOR_IN_VTP_BIN('FORCE_CHAIN_FT',FCHAIN_FT,PASS)
               CALL WRITE_VECTOR_IN_VTP_BIN('FORCE_CHAIN_FC',FCHAIN_FC,PASS)
               CALL WRITE_VECTOR_IN_VTP_BIN('FORCE_CHAIN_CONTACT_POINT',FCHAIN_CONTACT_POINT,PASS)
               CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_GLOBAL_ID1',DBLE(FCHAIN_GLOBAL_ID1),PASS)
               CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_GLOBAL_ID2',DBLE(FCHAIN_GLOBAL_ID2),PASS)
               CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_OVERLAP',DBLE(FCHAIN_OVERLAP),PASS)
               IF(VTK_FC_JKR_STATE(VTK_REGION)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_JKR_CONTACT_STATE',FCHAIN_JKR_STATE,PASS)
               IF(VTK_FC_JKR_RADIUS(VTK_REGION)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_JKR_CONTACT_RADIUS',FCHAIN_JKR_RADIUS,PASS)
               IF(VTK_FC_JKR_GAMMA(VTK_REGION)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_JKR_EFFECTIVE_SURFACE_ENERGY',FCHAIN_JKR_GAMMA,PASS)
               IF(VTK_FC_TANG_DISP(VTK_REGION)) &
                  CALL WRITE_VECTOR_IN_VTP_BIN('FORCE_CHAIN_TANGENTIAL_DISPLACEMENT',FCHAIN_TANGENTIAL_DISPLACEMENT,PASS)
               IF(VTK_FC_ROLL_DISP(VTK_REGION)) &
                  CALL WRITE_VECTOR_IN_VTP_BIN('FORCE_CHAIN_ROLLING_DISPLACEMENT',FCHAIN_ROLLING_DISPLACEMENT,PASS)
               IF(VTK_FC_TWIST_ANGLE(VTK_REGION)) &
                  CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_TWIST_ANGLE',FCHAIN_TWIST_ANGLE,PASS)

               ! Normalize a copy: subsequent regions still need dimensional FN.
               ALLOCATE(FCHAIN_FN_MAG_ND(FCHAINC))
               FCHAIN_FN_MAG_ND = FCHAIN_FN_MAG(1:FCHAINC)
               IF(PASS==WRITE_DATA) THEN
                  SUM_FN = ZERO
                  DO L = 1, FCHAINC
                     SUM_FN = SUM_FN + FCHAIN_FN_MAG(L)
                  ENDDO
                  call global_all_sum(SUM_FN, GLOBAL_SUM_FN)
                  call global_all_sum(FCHAINC, GLOBAL_FCHAINC)

                  IF(GLOBAL_SUM_FN>0.AND.GLOBAL_FCHAINC>0) THEN
                     AVG_FN = GLOBAL_SUM_FN/GLOBAL_FCHAINC
                     FCHAIN_FN_MAG_ND = FCHAIN_FN_MAG_ND/AVG_FN
                  ENDIF

               ENDIF
               CALL WRITE_SCALAR_IN_VTP_BIN('FORCE_CHAIN_FN_ND',FCHAIN_FN_MAG_ND,PASS)
               DEALLOCATE(FCHAIN_FN_MAG_ND)

! Add other variables from the vtk keywords. This will write a pair of
! values for each particle in the contact pair.

               IF(ALLOCATED(DES_RADIUS)) &
                  CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Diameter',2.0D0*DES_RADIUS,PASS)

               IF(VTK_PART_VEL(VTK_REGION).AND.ALLOCATED(DES_VEL_NEW)) &
                  CALL ADD_VECTOR_TO_FORCE_CHAIN_VTP_BIN('Velocity',DES_VEL_NEW,PASS)

               IF(CGDEM) THEN
                  IF(VTK_PART_PHYSICAL_DIAMETER(VTK_REGION).AND.ALLOCATED(DES_CGP_RPR)) &
                     CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Physical_diameter',2.0D0*DES_CGP_RPR,PASS)
                  IF(VTK_PART_CGP_STAT_WT(VTK_REGION).AND.ALLOCATED(DES_CGP_STW)) &
                     CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Statistical_weight',DES_CGP_STW,PASS)
               ENDIF

               IF(VTK_PART_ANGULAR_VEL(VTK_REGION).AND.ALLOCATED(OMEGA_NEW)) &
                  CALL ADD_VECTOR_TO_FORCE_CHAIN_VTP_BIN('Angular_velocity', OMEGA_NEW,PASS)

               IF(PARTICLE_ORIENTATION) THEN
                  IF(VTK_PART_ORIENTATION(VTK_REGION).AND.ALLOCATED(ORIENTATION)) &
                     CALL ADD_VECTOR_TO_FORCE_CHAIN_VTP_BIN('Orientation', ORIENTATION,PASS)
               ENDIF

               DO N=1, DES_USR_VAR_SIZE
                  IF(VTK_PART_USR_VAR(VTK_REGION,N).AND.ALLOCATED(DES_USR_VAR)) &
                    CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('User Defined Var '//trim(iVal(N)),DES_USR_VAR(N,:),PASS)
               ENDDO

               IF(ENERGY_EQ.AND.VTK_PART_TEMP(VTK_REGION).AND.ALLOCATED(DES_T_s)) &
                 CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Temperature', DES_T_s,PASS)

               IF(ANY_SOLIDS_SPECIES_EQ) THEN
                  DO N=1, DIMENSION_N_S
                     IF(VTK_PART_X_s(VTK_REGION,N).AND.ALLOCATED(DES_X_s)) THEN
                        M= VTK_PHASE_FOR_DES_X_S(VTK_REGION)
                        IF(M/=0) THEN
                           CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN(trim(SPECIES_ALIAS_s(M,N)), DES_X_s(:,N),PASS)
                        ELSE
                           CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN(trim(iVar('X_s',N)), DES_X_s(:,N),PASS)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF

               IF(SAVE_PART_RRATES) THEN
                  DO N=1, NO_OF_DES_RXNS
                     IF(VTK_PART_RRATE(VTK_REGION,N).AND.ALLOCATED(PART_RRATES_OUT)) THEN
                        CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('RRate_'//trim(DES_RXN_NAME(N)), PART_RRATES_OUT(:,N),PASS)
                     ENDIF
                  ENDDO
               ENDIF

               IF(VTK_PART_DENSITY(VTK_REGION).AND.ALLOCATED(RO_SOL)) &
                  CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Density', RO_SOL,PASS)

               IF(USE_COHESION.AND.VTK_PART_COHESION(VTK_REGION).AND.ALLOCATED(PostCohesive)) &
                  CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Cohesive_force', PostCohesive,PASS)

               IF(VTK_PART_COLL_FORCE(VTK_REGION).AND.ALLOCATED(DES_COLL_FORCE)) &
                  CALL ADD_VECTOR_TO_FORCE_CHAIN_VTP_BIN('Collision_force', DES_COLL_FORCE,PASS)


               IF(VTK_PART_RANK(VTK_REGION)) THEN
                  IF(PASS==WRITE_DATA) THEN
                     ALLOCATE(PARTICLE_RANK(MAX_PIP))
                     DO L = 1, MAX_PIP
                        PARTICLE_RANK(L) = DBLE(MyPE)
                     ENDDO
                  ENDIF

                  CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Particle_rank', PARTICLE_RANK,PASS)
                  IF(PASS==WRITE_DATA) DEALLOCATE(PARTICLE_RANK)
               ENDIF

               IF(VTK_PART_PHASE_ID(VTK_REGION)) THEN
                  IF(PASS==WRITE_DATA) THEN
                     ALLOCATE(PARTICLE_PHASE_ID(MAX_PIP))
                     DO L = 1, MAX_PIP
                        PARTICLE_PHASE_ID(L) = DBLE(PIJK(L,5))
                     ENDDO
                  ENDIF

                  CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Particle_Phase_ID', PARTICLE_PHASE_ID,PASS)
                  IF(PASS==WRITE_DATA) DEALLOCATE(PARTICLE_PHASE_ID)
               ENDIF

               IF(VTK_PART_ID(VTK_REGION).AND.ALLOCATED(iGLOBAL_ID)) &
                  CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Particle_ID', DBLE(iGLOBAL_ID),PASS)

               IF(VTK_PART_RESIDENCE_TIME(VTK_REGION).AND.ALLOCATED(RESIDENCE_TIME)) &
                  CALL ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN('Residence_Time', RESIDENCE_TIME,PASS)

	         ELSE !Glued particle data
               IF(gsp_explicit) THEN
                  IF(VTK_PART_GLUEDIAMETER(VTK_REGION) .AND. ALLOCATED(glueDiameter)) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Diameter', glueDiameter ,PASS)
                  ENDIF
                  IF(VTK_PART_GLUEBOUNDING(VTK_REGION) .AND. ALLOCATED(glueBounding)) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Bounding_Diameter', glueBounding ,PASS)
                  ENDIF
                  IF(VTK_PART_GLUEMASS(VTK_REGION) .AND. ALLOCATED(glueMASS)) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Mass', glueMASS ,PASS)
                  ENDIF
                  IF(VTK_PART_GLUEVOLUME(VTK_REGION) .AND. ALLOCATED(glueVolume)) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Volume', glueVolume ,PASS)
                  ENDIF
                  IF(VTK_PART_GLUEQUAT(VTK_REGION) .AND. ALLOCATED(glueQuat)) THEN
                     CALL WRITE_VECTOR_IN_VTP_BIN9('Glued_Particle_Q_Rotation', glueQuat(:,1:4),PASS)
                     CALL WRITE_QUATERNION_IN_VTP_BIN('Glued_Particle_Quaternion', glueQuat ,PASS)
                  ENDIF
               ENDIF
            ENDIF
	      ENDIF
      ENDDO ! PASS LOOP, EITHER HEADER OR DATA

      IF(GSP_IMPLICIT .and. VTK_GSP_OUTPUT_TYPE(VTK_REGION) == 1) THEN
         TMP_LOCAL_CNT = cur_n_sphere
         TMP_GLOBAL_CNT = total_n_sphere
      ELSE
         TMP_LOCAL_CNT = LOCAL_CNT
         TMP_GLOBAL_CNT = GLOBAL_CNT
      ENDIF

! Write connectivity and offset. This let Paraview define cells (one
! vertex per cell), and apply Clip and Threshold filters on it
      IF(TMP_GLOBAL_CNT>0.AND.MODE==0) THEN
         IF(myPE == PE_IO) THEN
! Write connectivity
            WRITE(UNIT_VTP)  TMP_GLOBAL_CNT*c_sizeof(L)
            DO L=0, TMP_GLOBAL_CNT-1
               WRITE(UNIT_VTP)  L
            ENDDO
! Write offset
            WRITE(UNIT_VTP)  TMP_GLOBAL_CNT*c_sizeof(L)
            DO L=0, TMP_GLOBAL_CNT-1
               WRITE(UNIT_VTP)  L+1
            ENDDO
         ENDIF
      ENDIF

      CALL CLOSE_VTP_FILE_BIN(MODE)

! Only update pvd file when there are particles in vtk region
      IF(TMP_GLOBAL_CNT>0.AND.MODE==0) CALL UPDATE_AND_CLOSE_PVD_FILE(TIMESTAMP)

#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif

! Update Frames
      IF (myPE == PE_IO.AND.TIME_DEPENDENT_FILENAME) THEN
         CALL WRITE_PVD_FRAMES
      ENDIF

      IF (FULL_LOG.AND.myPE == PE_IO) THEN
         IF(TMP_GLOBAL_CNT>0.AND.MODE==0) THEN
           ! WRITE(ERR_MSG,20)' DONE.'
         ELSE
            WRITE(ERR_MSG,20)' VTP file not written (zero particles in vtk region).'
         ENDIF
        CALL LOG_STATUS()
      ENDIF

20    FORMAT(A,1X/)

      CALL END_IO
      RETURN

      END SUBROUTINE WRITE_VTP_FILE

      SUBROUTINE UPDATE_FRAMES
         IMPLICIT NONE
         INTEGER :: L
	 RESET_FRAME_AT_TIME_ZERO = .TRUE.

         IF (TIME_DEPENDENT_FILENAME) THEN
            IF(MYPE==PE_IO) THEN
               CALL READ_PVD_FRAMES
               IF (RESET_FRAME_AT_TIME_ZERO.AND.TIME==ZERO) THEN
                  DO L = 1, DIMENSION_VTK
                     IF(L==VTK_REGION .AND. FRAME(L).ne.ZERO) FRAME(L)=-1
                  ENDDO
                  RESET_FRAME_AT_TIME_ZERO = .FALSE.
               ENDIF
               DO L = 1, DIMENSION_VTK
                  IF(L==VTK_REGION) FRAME(L) = FRAME(L) + 1
               ENDDO
            ENDIF
! Broadcast only required for Distributed IO
            IF(BDIST_IO) CALL BCAST(FRAME(1:DIMENSION_VTK))
         ENDIF
      END SUBROUTINE UPDATE_FRAMES
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_COMP_SPHERE_IMPLICIT                             C
!  Purpose: dump component sphere                                      C
!           when keyword "VTK_GSP_OUTPUT_TYPE == 1" (by default)       C
!                                                                      C
!  Author: Renjie Ke                                  Date: 06-NOV-24  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_COMP_SPHERE_IMPLICIT(PASS, cur_n_sphere, total_n_sphere)
         use constant, only: pi
         use discretelement, only: cspart_data, sphereNumber
         use resize, only: real_grow, real_grow2, real_grow2_reverse
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: PASS
         INTEGER, INTENT(IN) :: cur_n_sphere, total_n_sphere

         INTEGER :: L, JJ, MM, start_index, tmp_index, N,M
         INTEGER :: old_size, new_size
         double precision :: sphere_radius

! @renjieke
! if max_pip increases, NSphereGSP should increase so need to grow arrays
! The problem is we may not tracking the NSphereGSP
! so maybe using size(sphere_radius_out)
         IF(PASS == 2) THEN
         ! @renjieke: INTEGER :: WRITE_DATA = 2
         ! only when it is write_data, we need to do a grow check
         ! and then fill those output arrays
            old_size = size(sphere_radius_out)
            new_size = old_size
! Only when current output array size is not large enough then grow
            IF(old_size < cur_n_sphere) THEN
               DO WHILE (new_size < cur_n_sphere)
                  new_size = 2 * new_size
               ENDDO

! a lot of grow operation
               call real_grow(sphere_radius_out ,new_size)
               call real_grow(sphere_rho_out ,new_size)
               call real_grow(sphere_mass_out ,new_size)
               call real_grow(sphere_vol_out ,new_size)
               call real_grow(sphere_rank_out ,new_size)
               call real_grow(sphere_pid_out ,new_size)
               call real_grow(sphere_id_out ,new_size)
               call real_grow(sphere_res_time_out,new_size)
               call real_grow(sphere_coordination_out,new_size)
               call real_grow(sphere_MIN_OVERLAP_out,new_size)
               call real_grow(sphere_MAX_OVERLAP_out,new_size)
               call real_grow(sphere_MEAN_OVERLAP_out,new_size)

               call real_grow2_reverse(sphere_rot_out,new_size)
               call real_grow2_reverse(sphere_vel_out,new_size)
               call real_grow2_reverse(sphere_orient_out,new_size)

               call real_grow2_reverse(sphere_coll_force_out,new_size)
               call real_grow(sphere_drag_out,new_size)

               call real_grow(sphere_total_mass_out, new_size)
               call real_grow(sphere_bounding_out, new_size)
               call real_grow(sphere_eqv_dia_out, new_size)

               IF(ENERGY_EQ) THEN
                  call real_grow(sphere_temp_out,new_size)
                  call real_grow(sphere_cps_out,new_size)
               ENDIF
               if(DES_USR_VAR_SIZE > 0) &
                  call real_grow2(sphere_usr_var_out,new_size)
               IF(ANY_SOLIDS_SPECIES_EQ) &
                  call real_grow2_reverse(sphere_xs_out,new_size)
               if (SAVE_PART_RRATES) &
                  call real_grow2_reverse(sphere_prates_out,new_size)
            ENDIF

! fill those arrays
            start_index = 0
            DO L = 1, MAX_PIP
               IF(.NOT. BELONGS_TO_VTK_SUBDOMAIN(L)) CYCLE
               MM = PIJK(L,5)

               DO JJ = 1, sphereNumber(MM)
                  sphere_radius = cspart_data(JJ,4,MM) * des_radius(L)
                  sphere_radius_out(start_index + JJ) = sphere_radius
                  sphere_rho_out(start_index + JJ) = RO_SOL(L)
                  sphere_vol_out(start_index + JJ) = 4.0/3.0*PI*sphere_radius**3.0
                  sphere_mass_out(start_index + JJ) = sphere_vol_out(start_index + JJ)*RO_SOL(L)
                  sphere_rank_out(start_index + JJ) = DBLE(mype)
                  sphere_pid_out(start_index + JJ) = DBLE(PIJK(L,5))
                  sphere_id_out(start_index + JJ) = DBLE(iGLOBAL_ID(L))
                  sphere_res_time_out(start_index + JJ) = RESIDENCE_TIME(L)
                  sphere_coordination_out(start_index + JJ) = COORDINATION(L)
                  sphere_MIN_OVERLAP_out(start_index + JJ) = MIN_OVERLAP(L)
                  sphere_MAX_OVERLAP_out(start_index + JJ) = MAX_OVERLAP(L)
                  sphere_MEAN_OVERLAP_out(start_index + JJ) = MEAN_OVERLAP(L)
                  sphere_rot_out(start_index + JJ,:) = omega_new(L,:)
                  sphere_vel_out(start_index + JJ,:) = des_vel_new(L,:)

                  sphere_coll_force_out(start_index + JJ,:) = DES_COLL_FORCE(L,:)
                  sphere_drag_out(start_index + JJ) = F_GP(L)

                  sphere_total_mass_out(start_index + JJ) = PMASS(L)
                  sphere_bounding_out(start_index + JJ) = glueBounding(L)
                  sphere_eqv_dia_out(start_index + JJ) = glueDiameter(L)

                  IF(PARTICLE_ORIENTATION) &
                     sphere_orient_out(start_index + JJ,:) = ORIENTATION(L,:)

                  IF(ENERGY_EQ) THEN
                     sphere_temp_out(start_index + JJ) = des_t_s(L)
                     sphere_cps_out(start_index + JJ) = des_c_Ps(L)
                  ENDIF

                  if(DES_USR_VAR_SIZE > 0) &
                     sphere_usr_var_out(:,start_index + JJ) = DES_USR_VAR(:,L)
                  IF(ANY_SOLIDS_SPECIES_EQ) &
                     sphere_xs_out(start_index + JJ,:) = des_x_s(L,:)
                  if (SAVE_PART_RRATES) &
                     sphere_prates_out(start_index + JJ,:) = PART_RRATES_OUT(L,:)
                  tmp_index = start_index + JJ
               ENDDO
               start_index = tmp_index
            ENDDO

         ENDIF ! end of IF(PASS == 2)
! start writing
! @renjieke WRITE_SCALAR_IN_VTP_BIN: nbytes_scalar = GLOBAL_CNT * c_sizeof(float)
! need to change this line to make sure the size match
         IF(ALLOCATED(sphere_radius_out)) THEN
            CALL WRITE_SCALAR_IN_VTP_BIN('Diameter',2.0D0*sphere_radius_out,PASS,cur_n_sphere,total_n_sphere)
            IF(VTK_PART_GLUEBOUNDING(VTK_REGION)) &
               CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Bounding_Diameter',sphere_bounding_out,PASS,cur_n_sphere,total_n_sphere)
            IF(VTK_PART_GLUEDIAMETER(VTK_REGION)) &
               CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Equivalent_Diameter',sphere_eqv_dia_out,PASS,cur_n_sphere,total_n_sphere)
         ENDIF

         IF(VTK_PART_VEL(VTK_REGION).AND.ALLOCATED(DES_VEL_NEW)) &
            CALL WRITE_VECTOR_IN_VTP_BIN('Velocity',sphere_vel_out,PASS,cur_n_sphere,total_n_sphere)

         IF(VTK_PART_ANGULAR_VEL(VTK_REGION).AND.ALLOCATED(OMEGA_NEW)) &
            CALL WRITE_VECTOR_IN_VTP_BIN('Angular_velocity',sphere_rot_out,PASS,cur_n_sphere,total_n_sphere)

         IF(PARTICLE_ORIENTATION) THEN
            IF(VTK_PART_ORIENTATION(VTK_REGION).AND.ALLOCATED(ORIENTATION)) &
               CALL WRITE_VECTOR_IN_VTP_BIN('Orientation',sphere_orient_out,PASS,cur_n_sphere,total_n_sphere)
         ENDIF

         DO N=1, DES_USR_VAR_SIZE
            IF(VTK_PART_USR_VAR(VTK_REGION,N).AND.ALLOCATED(DES_USR_VAR)) &
               CALL WRITE_SCALAR_IN_VTP_BIN('User Defined Var '//trim(iVal(N)),sphere_usr_var_out(N,:),PASS,cur_n_sphere,total_n_sphere)
         ENDDO

         IF(ENERGY_EQ.AND.VTK_PART_TEMP(VTK_REGION).AND.ALLOCATED(DES_T_s)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Temperature', sphere_temp_out,PASS,cur_n_sphere,total_n_sphere)
         IF(ENERGY_EQ.AND.VTK_PART_C_Ps(VTK_REGION).AND.ALLOCATED(DES_C_Ps)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Specific_heat', sphere_cps_out,PASS,cur_n_sphere,total_n_sphere)

         IF(ANY_SOLIDS_SPECIES_EQ) THEN
            DO N=1, DIMENSION_N_S
               IF(VTK_PART_X_s(VTK_REGION,N).AND.ALLOCATED(DES_X_s)) THEN
                  M = VTK_PHASE_FOR_DES_X_S(VTK_REGION)
                  IF(M/=0) THEN
                     CALL WRITE_SCALAR_IN_VTP_BIN(trim(SPECIES_ALIAS_s(M,N)), sphere_xs_out(:,N),PASS,cur_n_sphere,total_n_sphere)
                  ELSE
                     CALL WRITE_SCALAR_IN_VTP_BIN(trim(iVar('X_s',N)), sphere_xs_out(:,N),PASS,cur_n_sphere,total_n_sphere)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         IF(SAVE_PART_RRATES) THEN
            DO N=1, NO_OF_DES_RXNS
               IF(VTK_PART_RRATE(VTK_REGION,N).AND.ALLOCATED(PART_RRATES_OUT)) THEN
                  CALL WRITE_SCALAR_IN_VTP_BIN('RRate_'//trim(DES_RXN_NAME(N)), sphere_prates_out(:,N),PASS,cur_n_sphere,total_n_sphere)
               ENDIF
            ENDDO
         ENDIF

         IF(VTK_PART_DENSITY(VTK_REGION).AND.ALLOCATED(RO_SOL)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Density', sphere_rho_out,PASS,cur_n_sphere,total_n_sphere)
         IF(VTK_PART_MASS(VTK_REGION).AND.ALLOCATED(PMASS)) THEN
            CALL WRITE_SCALAR_IN_VTP_BIN('Mass', sphere_mass_out,PASS,cur_n_sphere,total_n_sphere)
            IF(VTK_PART_GLUEMASS(VTK_REGION)) &
               CALL WRITE_SCALAR_IN_VTP_BIN('Glued_Particle_Mass', sphere_total_mass_out,PASS,cur_n_sphere,total_n_sphere)
         ENDIF
         IF(VTK_PART_VOLUME(VTK_REGION).AND.ALLOCATED(PVOL)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Volume', sphere_vol_out,PASS,cur_n_sphere,total_n_sphere)

         IF(VTK_PART_COLL_FORCE(VTK_REGION).AND.ALLOCATED(DES_COLL_FORCE)) &
            CALL WRITE_VECTOR_IN_VTP_BIN('Collision_force', sphere_coll_force_out,PASS,cur_n_sphere,total_n_sphere)

         IF(VTK_PART_DRAG_COEFF(VTK_REGION).AND.ALLOCATED(F_GP)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Drag_coefficient', sphere_drag_out,PASS,cur_n_sphere,total_n_sphere)

         ! @renjieke not yet setup arrays
         ! IF(USE_COHESION.AND.VTK_PART_COHESION(VTK_REGION).AND.ALLOCATED(PostCohesive)) &
         !    CALL WRITE_SCALAR_IN_VTP_BIN('CohesiveForce', PostCohesive,PASS)
         ! IF(VTK_PART_COLL_FORCE(VTK_REGION).AND.ALLOCATED(DES_COLL_FORCE)) &
         !    CALL WRITE_VECTOR_IN_VTP_BIN('Collision_force', DES_COLL_FORCE,PASS)
         ! IF(VTK_PART_DRAG_COEFF(VTK_REGION).AND.ALLOCATED(F_GP)) &
         !    CALL WRITE_SCALAR_IN_VTP_BIN('Drag_coefficient', F_GP,PASS)
         ! IF(VTK_PART_RANK(VTK_REGION)) THEN
         !    IF(PASS==WRITE_DATA) THEN
         !       ALLOCATE(PARTICLE_RANK(MAX_PIP))
         !       DO L = 1, MAX_PIP
         !          PARTICLE_RANK(L) = DBLE(MyPE)
         !       ENDDO
         !    ENDIF

         !    CALL WRITE_SCALAR_IN_VTP_BIN('Particle_Rank', PARTICLE_RANK,PASS)
         !    IF(PASS==WRITE_DATA) DEALLOCATE(PARTICLE_RANK)
         ! ENDIF
         IF(VTK_PART_RANK(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Particle_Rank', sphere_rank_out,PASS,cur_n_sphere,total_n_sphere)

         IF(VTK_PART_PHASE_ID(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Particle_Phase_ID', sphere_pid_out,PASS,cur_n_sphere,total_n_sphere)

         IF(VTK_PART_ID(VTK_REGION).AND.ALLOCATED(iGLOBAL_ID)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Particle_ID', sphere_id_out,PASS,cur_n_sphere,total_n_sphere)

         IF(VTK_PART_RESIDENCE_TIME(VTK_REGION).AND.ALLOCATED(RESIDENCE_TIME)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Residence_Time', sphere_res_time_out,PASS,cur_n_sphere,total_n_sphere)

         IF(VTK_PART_COORDINATION(VTK_REGION).AND.ALLOCATED(COORDINATION)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Coordination', sphere_coordination_out,PASS,cur_n_sphere,total_n_sphere)

         IF(VTK_PART_OVERLAP_STATS(VTK_REGION).AND.ALLOCATED(MIN_OVERLAP)) THEN
            CALL WRITE_SCALAR_IN_VTP_BIN('Min_overlap', sphere_MIN_OVERLAP_out,PASS,cur_n_sphere,total_n_sphere)
            CALL WRITE_SCALAR_IN_VTP_BIN('Max_overlap', sphere_MAX_OVERLAP_out,PASS,cur_n_sphere,total_n_sphere)
            CALL WRITE_SCALAR_IN_VTP_BIN('Mean_overlap', sphere_MEAN_OVERLAP_out,PASS,cur_n_sphere,total_n_sphere)

! Non-dimensionalized overlap (by diameter)
               ! @renjieke can't just do this, need new arrays
               ! IF(PASS==WRITE_DATA) THEN
               !    WHERE(DES_RADIUS>ZERO) MIN_OVERLAP = sphere_MIN_OVERLAP_out/(2.0D0*sphere_radius_out)
               !    WHERE(DES_RADIUS>ZERO) MAX_OVERLAP = sphere_MAX_OVERLAP_out/(2.0D0*sphere_radius_out)
               !    WHERE(DES_RADIUS>ZERO) MEAN_OVERLAP = sphere_MEAN_OVERLAP_out/(2.0D0*sphere_radius_out)
               ! ENDIF

               ! CALL WRITE_SCALAR_IN_VTP_BIN('Min_overlap/Diameter', MIN_OVERLAP, PASS)
               ! CALL WRITE_SCALAR_IN_VTP_BIN('Max_overlap/Diameter', MAX_OVERLAP,PASS)
               ! CALL WRITE_SCALAR_IN_VTP_BIN('Mean_overlap/Diameter', MEAN_OVERLAP,PASS)
         ENDIF

         RETURN

      END SUBROUTINE WRITE_COMP_SPHERE_IMPLICIT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_VTP_FILE_BIN                                      C
!  Purpose: Open a vtp file and writes the header                      C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_VTP_FILE_BIN(MODE,TIMESTAMP)

      IMPLICIT NONE
      LOGICAL :: NEED_TO_WRITE_VTP
      INTEGER :: ISTAT
      CHARACTER(256) :: io_message
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)
      DOUBLE PRECISION :: TIMESTAMP
      INTEGER :: M

      IF(BDIST_IO) THEN
         NEED_TO_WRITE_VTP = (LOCAL_CNT>0)
      ELSE
         NEED_TO_WRITE_VTP = (MyPE==0.AND.GLOBAL_CNT>0)
      ENDIF

! Only open the file from head node when not using distributed I/O
      IF (myPE /= PE_IO.AND.(.NOT.BDIST_IO)) THEN
         RETURN
      END IF

      CALL UPDATE_FRAMES

      CALL OPEN_VTP(TIMESTAMP)

   CONTAINS

      SUBROUTINE OPEN_VTP(TIMESTAMP)

      DOUBLE PRECISION :: TIMESTAMP

         IF(.NOT.BDIST_IO) THEN
            VTP_FILENAME = trim(FIND_VTK_BASENAME(MODE)) // ".vtp"
         ELSE
            VTP_FILENAME = trim(FIND_VTK_BASENAME(MODE,MyPE)) // ".vtp"
         ENDIF

         IF(TRIM(VTU_DIR)/='.'.AND.TRIM(VTU_DIR)/='' .AND. myPE == PE_IO) THEN
            CALL CREATE_DIR(trim(VTU_DIR))
         END IF

! Echo
      IF (FULL_LOG) THEN
         IF (.NOT.BDIST_IO) THEN
            WRITE(ERR_MSG,10)' WRITING VTP FILE : ', TRIM(VTP_FILENAME),' .'
            CALL LOG_STATUS()
         ELSE
            WRITE (ERR_MSG, 20) ' WRITING PVTP FILE : ', trim(FIND_VTK_BASENAME(MODE)), '.pvtp (EACH PROCESSOR IS WRITING ITS OWN VTP FILE)'
            CALL LOG_STATUS()
         ENDIF
      ENDIF

! Open File

      IF (NEED_TO_WRITE_VTP) THEN

         OPEN(UNIT     = UNIT_VTP,           &
              FILE     = TRIM(VTP_FILENAME), &
              FORM     = 'UNFORMATTED',      &
              ACCESS   = 'STREAM',           &
              ACTION   = 'WRITE', CONVERT  = 'BIG_ENDIAN', IOSTAT=ISTAT, IOMSG=io_message)

         IF (ISTAT /= 0) THEN
            WRITE(ERR_MSG, "(/1X,/, A, &
               /10X, 'Unable to write to filename:  ',A, &
               /10X, 'Verify that VTU_DIR exists: ', A, /1X)") io_message, VTU_FILENAME, VTU_DIR
            call log_error()
         ENDIF

! Write file Header
         BUFFER='<?xml version="1.0"?>'
         WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

         WRITE(BUFFER,*)'<!-- Time =',TIMESTAMP,' sec. -->'
         WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

         BUFFER='<VTKFile type="PolyData" version="0.1" byte_order="BigEndian" header_type="UInt64">'
         WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

         BUFFER='  <PolyData>'
         WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! Write Superquadric parameters as Field Data
         IF(SuperDEM) THEN
            BUFFER='    <FieldData>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
            DO M=1,DIM_M
               IF(SOLIDS_MODEL(M) == 'SQP') THEN
                  WRITE(BUFFER,100) '      <DataArray type="Float32" Name="SQP_',M,'" NumberOfTuples="5">',&
                                     SQP_a(M),SQP_b(M),SQP_c(M),SQP_m(M),SQP_n(M),&
                                    '</DataArray>'
                  WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
               ENDIF
            ENDDO
            BUFFER='    </FieldData>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
         ENDIF

100 FORMAT(A,I2.2,A,5(E16.8),A)

      ENDIF
! For distributed I/O, open .p))vtp file that combines all *.vtp files for a given FRAME
! this is a simple ASCII file

      IF (myPE == PE_IO.AND.BDIST_IO.AND.GLOBAL_CNT>0) THEN

         PVTU_FILENAME = FIND_VTK_BASENAME(MODE) // ".pvtp"

         OPEN(UNIT = UNIT_PVTP, FILE = TRIM(PVTU_FILENAME))

         WRITE(UNIT_PVTP,105) '<?xml version="1.0"?>'
         WRITE(UNIT_PVTP,110) '<!-- Time =',TIMESTAMP,' sec. -->'
         WRITE(UNIT_PVTP,120) '<VTKFile type="PPolyData"',&
                  ' version="0.1" byte_order="BigEndian">'

         WRITE(UNIT_PVTP,*) '  <PPolyData GhostLevel="0">'
         WRITE(UNIT_PVTP,*) '      <PPoints>'
         WRITE(UNIT_PVTP,*) '        <PDataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
              &format="appended" offset=" 0" />'
         WRITE(UNIT_PVTP,*) '      </PPoints>'
         WRITE(UNIT_PVTP,*) ''
         WRITE(UNIT_PVTP,*) '      <PPointData Scalars="Diameter" Vectors="Velocity">'

      ENDIF

105   FORMAT(A)
110   FORMAT(A,E14.7,A)
120   FORMAT(A,A)
10    FORMAT(/1X,3A)
15    FORMAT(/1X,A)
20    FORMAT(/1X,3A)

      RETURN

   END SUBROUTINE OPEN_VTP

   END SUBROUTINE OPEN_VTP_FILE_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_GEOMETRY_IN_VTP_BIN                              C
!  Purpose: Write Geometry and connectivity in a vtu file              C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number 1                                  Date: 01-Nov-24  C
!  Author: Renjie Ke                                                   C
!  Purpose: Compatible with gsp implicit method                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_GEOMETRY_IN_VTP_BIN(PASS,cur_pe_n_sphere,total_n_sphere)

      use discretelement, only: cspart_data, sphereNumber
      IMPLICIT NONE

      INTEGER, INTENT(IN), OPTIONAL :: cur_pe_n_sphere, total_n_sphere

      REAL(c_float) :: float
      ! INTEGER(c_long) :: long

      INTEGER(KIND=8) ::     nbytes_vector
      INTEGER(KIND=8) ::     offset_xyz

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)  ! local
      DOUBLE PRECISION, ALLOCATABLE :: gtemp_array(:,:)  ! global

      INTEGER :: LB, UB
      INTEGER :: PC, LC1, LC2
      LOGICAL :: gsp_implicit_writing
      INTEGER :: TMP_GLOBAL_CNT, TMP_LOCAL_CNT
      INTEGER :: MM, start_index, lproc
      DOUBLE PRECISION :: xyzloc(3), sphere_pos(3)

      gsp_implicit_writing = present(cur_pe_n_sphere) .and. present(total_n_sphere)

      if(gsp_implicit_writing) then
         TMP_LOCAL_CNT = cur_pe_n_sphere
         TMP_GLOBAL_CNT = total_n_sphere
      else
         TMP_LOCAL_CNT = LOCAL_CNT
         TMP_GLOBAL_CNT = GLOBAL_CNT
      endif

! Loop through all particles and kee a list of particles belonging to a VTK region

! Since the data is appended (i.e., written after all tags), the
! offset, in number of bytes must be specified.  The offset includes
! the size of the data for each field, plus the size of the integer
! that stores the number of bytes.  this is why the offset of a field
! equals the offset of the previous field plus sizeof(int) plus the
! number of bytes of the field.

! Next, the actual data is written for the geometry (PASS=WRITE_DATA)
! The DATA is converted to single precision to save memory.

      IF (.NOT.BDIST_IO) THEN
! The number of points in the pvd file is the global number of particles
! computed from SETUP_VTK_REGION_PARTICLES

         NUMBER_OF_POINTS = TMP_GLOBAL_CNT

! Number of bytes of position field (vector,3 components)
         nbytes_vector = INT(NUMBER_OF_POINTS,8) * 3 * c_sizeof(float)

! Offset of each field
         offset_xyz = 0


         IF(PASS==WRITE_HEADER) THEN
            IF(myPE == PE_IO) THEN

               WRITE(BUFFER,*)'    <Piece NumberOfPoints="',NUMBER_OF_POINTS, &
! JFD: I am turning off the vertex data until a solution is found. The issue is c_sizeof(long) doesn't return the same
! value on Linux and Windows. That makes the vtp file unreadable by Paraview and the GUI.
!                     '"  NumberOfVerts="',NUMBER_OF_POINTS,'" NumberOfLines ="0" NumberOfStrips="0" NumberOfPolys="0" >'
                     '"  NumberOfVerts="', 0,'" NumberOfLines ="0" NumberOfStrips="0" NumberOfPolys="0" >'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Points>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
                                       &format="appended" offset="',offset_xyz,'" />'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      </Points>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'<PointData Scalars="Diameter" Vectors="Velocity"> '!preparing pointData
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! calculate offset for next field
               VTU_offset = offset_xyz + VTK_HEADER_SIZE + nbytes_vector

            ENDIF

         ELSEIF(PASS==WRITE_DATA) THEN

            IF(myPE == PE_IO) THEN

               WRITE(BUFFER,*)'      </PointData>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! Vertex data. This allows Paraview to apply Clip and Threshold filters
               WRITE(BUFFER,*)'      <Verts> </Verts>'
! JFD: I am turning off the vertex data until a solution is found. The issue is c_sizeof(long) doesn't return the same
! value on Linux and Windows. That makes the vtp file unreadable by Paraview and the GUI.
!               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
!               VTU_offset = VTU_offset - c_sizeof(int) + c_sizeof(long) ! Correction because last offset assumed int data
!               WRITE(BUFFER,*)'         <DataArray type="Int64" Name="connectivity" format="appended" RangeMin="" RangeMax="" offset="',VTU_offset,'" />'
!               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
!               VTU_offset = VTU_offset + 1*c_sizeof(long)*NUMBER_OF_POINTS + 1*c_sizeof(long)
!               WRITE(BUFFER,*)'         <DataArray type="Int64" Name="offsets" format="appended" RangeMin="" RangeMax="" offset="', VTU_offset,'" />'
!               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
!               WRITE(BUFFER,*)'      </Verts>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Lines> </Lines>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Strips> </Strips>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Polys> </Polys>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'    </Piece>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'  </PolyData>'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'  <AppendedData encoding="raw">'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! Starting raw binary data with an underscore

               WRITE(BUFFER,*)'_'
               WRITE(UNIT_VTP)TRIM(BUFFER)

! Number of bytes for X,Y,Z coordinates
               WRITE(UNIT_VTP) nbytes_vector
         ENDIF

         IF(VTP_P_DATA) THEN
            LB = LBOUND(DES_POS_NEW,2) ! This should always be 1
            UB = UBOUND(DES_POS_NEW,2) ! This should always be 2
         ELSEIF(VTP_GLUE_DATA) THEN
	         LB = LBOUND(gluePos,2) ! This should always be 1
            UB = UBOUND(gluePos,2) ! This should always be 2
	      ELSE
            LB = LBOUND(FCHAIN_MIDPOINT,2) ! This should always be 1
            UB = UBOUND(FCHAIN_MIDPOINT,2) ! This should always be 2
         ENDIF

         ALLOCATE (dProcBuf(TMP_LOCAL_CNT) )
         ALLOCATE (ltemp_array((UB-LB)+1,TMP_LOCAL_CNT))

         IF(myPE == PE_IO) THEN
            ALLOCATE (dRootBuf(TMP_GLOBAL_CNT))
            ALLOCATE (gtemp_array((UB-LB)+1,TMP_GLOBAL_CNT))
         ELSE
            ALLOCATE (dRootBuf(10))
            ALLOCATE (gtemp_array((UB-LB)+1,10))
         ENDIF

! Pack particle coordinates in a temporary local array
         IF(VTP_P_DATA) THEN ! Regular particle data
            PC = 0
            IF(gsp_implicit_writing) THEN
               DO LC1 = 1, MAX_PIP
                  IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     ! have to use cspart_data to map the gsp to component spheres
                     MM = PIJK(LC1,5)
                     DO LC2 = 1, sphereNumber(MM)
                        xyzloc(1) = cspart_data(LC2,1,MM) * des_radius(LC1) * 2.0
                        xyzloc(2) = cspart_data(LC2,2,MM) * des_radius(LC1) * 2.0
                        xyzloc(3) = cspart_data(LC2,3,MM) * des_radius(LC1) * 2.0

                        sphere_pos(:) = des_pos_new(LC1,:) + xyzloc(1)*glueEX(LC1,:) &
                           + xyzloc(2)*glueEY(LC1,:) + xyzloc(3)*glueEZ(LC1,:)
                        PC = PC + 1
                        ltemp_array(:,PC) = sphere_pos(:)
                     ENDDO
                  ENDIF
                  IF(PC==TMP_LOCAL_CNT) EXIT
               ENDDO
            ELSE
               DO LC1 = 1, MAX_PIP
                  IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     PC =PC + 1
                     DO LC2=LB, UB
                        ltemp_array(LC2,PC) = DES_POS_NEW(LC1,LC2)
                     ENDDO
                  ENDIF
                  IF(PC==TMP_LOCAL_CNT) EXIT
               ENDDO
            ENDIF
	      ELSEIF(VTP_GLUE_DATA) THEN ! glued particle data
	         PC = 0
            DO LC1 = 1, NGluedParticles
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = gluePos(LC1,LC2)
                  ENDDO
               ENDIF
               IF(PC==TMP_LOCAL_CNT) EXIT
            ENDDO
         ELSE ! Force chain data
            PC = 0
            DO LC1 = 1, FCHAINC
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = FCHAIN_MIDPOINT(LC1,LC2)
                  ENDDO
               ENDIF
               IF(PC==TMP_LOCAL_CNT) EXIT
            ENDDO
         ENDIF

! For each coordinate (x,y, and z), gather the local list to global temporary array
         DO LC2 = LB, UB
            dprocbuf(1:TMP_LOCAL_CNT)=ltemp_array(LC2,1:TMP_LOCAL_CNT)
            CALL desmpi_gatherv(ptype=2)
            gtemp_array(LC2,:) = drootbuf(:)
         ENDDO

! Write the list of coordinates
         IF(myPE == PE_IO) THEN
            DO LC1=1, TMP_GLOBAL_CNT
               DO LC2=LB, UB
                  WRITE(UNIT_VTP)  real(gtemp_array(LC2,LC1))
               ENDDO
            ENDDO
         ENDIF

         if(allocated(dProcBuf)) deallocate(dProcBuf)
         if(allocated(dRootBuf)) deallocate(dRootBuf)
         if(allocated(ltemp_array)) deallocate(ltemp_array)
         if(allocated(gtemp_array)) deallocate(gtemp_array)

         ENDIF ! end of pass == write_data

      ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN

         IF(TMP_LOCAL_CNT==0) RETURN
! The number of points in the pvd file is the local number of particles
! computed from SETUP_VTK_REGION_PARTICLES

         NUMBER_OF_POINTS = TMP_LOCAL_CNT

! Number of bytes of position field (vector,3 components)
         nbytes_vector = INT(NUMBER_OF_POINTS,8) * 3 * c_sizeof(float)

! Offset of each field
         offset_xyz = 0

         IF(PASS==WRITE_HEADER) THEN

            WRITE(BUFFER,*)'    <Piece NumberOfPoints="',NUMBER_OF_POINTS, &
                  '"  NumberOfVerts="0" NumberOfLines ="0" NumberOfStrips="0" NumberOfPolys="0" >'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Points>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
                                    &format="appended" offset="',offset_xyz,'" />'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      </Points>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'<PointData Scalars="Diameter" Vectors="Velocity"> '!preparing pointData
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! calculate offset for next field
            VTU_offset = offset_xyz + VTK_HEADER_SIZE + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            WRITE(BUFFER,*)'      </PointData>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Verts> </Verts>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Lines> </Lines>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Strips> </Strips>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Polys> </Polys>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'    </Piece>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'  </PolyData>'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'  <AppendedData encoding="raw">'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! Starting raw binary data with an underscore

            WRITE(BUFFER,*)'_'
            WRITE(UNIT_VTP)TRIM(BUFFER)

! Number of bytes for X,Y,Z coordinates
            WRITE(UNIT_VTP) nbytes_vector

            IF(VTP_P_DATA) THEN
               LB = LBOUND(DES_POS_NEW,2) ! This should always be 1
               UB = UBOUND(DES_POS_NEW,2) ! This should always be 2
            ELSE
               LB = LBOUND(FCHAIN_MIDPOINT,2) ! This should always be 1
               UB = UBOUND(FCHAIN_MIDPOINT,2) ! This should always be 2
            ENDIF

            ALLOCATE (ltemp_array((UB-LB)+1,TMP_LOCAL_CNT))

! Pack particle coordinates in a temporary local array
            IF(VTP_P_DATA) THEN ! Regular particle data
               PC = 0

               IF(gsp_implicit_writing) THEN
                  DO LC1 = 1, MAX_PIP
                     IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                        ! have to use cspart_data to map the gsp to component spheres
                        MM = PIJK(LC1,5)
                        DO LC2 = 1, sphereNumber(MM)
                           xyzloc(1) = cspart_data(LC2,1,MM) * des_radius(LC1) * 2.0
                           xyzloc(2) = cspart_data(LC2,2,MM) * des_radius(LC1) * 2.0
                           xyzloc(3) = cspart_data(LC2,3,MM) * des_radius(LC1) * 2.0

                           sphere_pos(:) = des_pos_new(LC1,:) + xyzloc(1)*glueEX(LC1,:) &
                              + xyzloc(2)*glueEY(LC1,:) + xyzloc(3)*glueEZ(LC1,:)

                           PC = PC + 1
                           ltemp_array(:,PC) = sphere_pos(:)
                        ENDDO
                     ENDIF
                     IF(PC==TMP_LOCAL_CNT) EXIT
                  ENDDO
               ELSE
                  DO LC1 = 1, MAX_PIP
                     IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                        PC =PC + 1
                        DO LC2=LB, UB
                           ltemp_array(LC2,PC) = DES_POS_NEW(LC1,LC2)
                        ENDDO
                     ENDIF
                     IF(PC==TMP_LOCAL_CNT) EXIT
                  ENDDO
               ENDIF
	         ELSEIF(VTP_GLUE_DATA) THEN ! glued particle data
               PC = 0
               DO LC1 = 1, NGluedParticles
                  IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     PC =PC + 1
                     DO LC2=LB, UB
                        ltemp_array(LC2,PC) = gluePos(LC1,LC2)
                     ENDDO
                  ENDIF
                  IF(PC==TMP_LOCAL_CNT) EXIT
               ENDDO
            ELSE ! Force chain data
               PC = 0
               DO LC1 = 1, FCHAINC
                  IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     PC =PC + 1
                     DO LC2=LB, UB
                        ltemp_array(LC2,PC) = FCHAIN_MIDPOINT(LC1,LC2)
                     ENDDO
                  ENDIF
                  IF(PC==TMP_LOCAL_CNT) EXIT
               ENDDO
            ENDIF

! Write the list of coordinates
            DO LC1=1, TMP_LOCAL_CNT
               DO LC2=LB, UB
                  WRITE(UNIT_VTP)  real(ltemp_array(LC2,LC1))
               ENDDO
            ENDDO

            deallocate (ltemp_array)

         ENDIF

      ENDIF

      RETURN

      END SUBROUTINE WRITE_GEOMETRY_IN_VTP_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SCALAR_IN_VTP_BIN                                C
!  Purpose: Write Scalar variable in a vtp file                        C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_SCALAR_IN_VTP_BIN(VAR_NAME,VAR,PASS,cur_pe_n_sphere,total_n_sphere)

      IMPLICIT NONE
      INTEGER :: I,LC1,PC

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, INTENT(in) :: VAR(:)
      INTEGER, INTENT(IN), OPTIONAL :: cur_pe_n_sphere
      INTEGER, INTENT(IN), OPTIONAL :: total_n_sphere

      REAL(c_float) :: float

      INTEGER(KIND=8) :: nbytes_scalar

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2
      LOGICAL :: gsp_implicit_writing
      gsp_implicit_writing = present(cur_pe_n_sphere) .and. present(total_n_sphere)
      IF (.NOT.BDIST_IO) THEN
! Number of bytes for each scalar field
! @renjieke, if global_cnt can be translated to total_n_sphere / total_n_sphere should be == global_cnt * sphereNumber(MM)/
         nbytes_scalar = INT(GLOBAL_CNT,8) * c_sizeof(float)
         if(gsp_implicit_writing) nbytes_scalar = INT(total_n_sphere,8) * c_sizeof(float)

         IF(PASS==WRITE_HEADER) THEN

! Remove possible white space with underscore
            DO I = 1,LEN_TRIM(VAR_NAME)
               IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
            ENDDO

            IF (myPE == PE_IO) THEN
! For each scalar, write a tag, with corresponding offset
               WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                    TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
            ENDIF

! Prepare the offset for the next field
            VTU_offset = VTU_offset + VTK_HEADER_SIZE + nbytes_scalar


         ELSEIF(PASS==WRITE_DATA) THEN
            IF(gsp_implicit_writing) THEN
               allocate (dProcBuf(cur_pe_n_sphere) )
               allocate (dRootBuf(total_n_sphere))
            ELSE
               allocate (dProcBuf(LOCAL_CNT) )
               allocate (dRootBuf(GLOBAL_CNT))
            ENDIF

! Pack scalar list in a local buffer before gathering to root
            PC = 0

            IF(gsp_implicit_writing) THEN
               ! @renjieke BELONGS_TO_VTK_SUBDOMAIN(LC1)) already filter the arrays, so no need to put any if condition here
               DO LC1 = 1, size(VAR(:),1)
                  ! if no particles in current rank, no need to append data
                  IF(size(dProcBuf) == 0) CYCLE
                  PC = PC + 1
                  dProcBuf(PC) = VAR(LC1)
                  IF(PC == cur_pe_n_sphere) EXIT
               ENDDO
            ELSE
               DO LC1 = 1, VTP_LUB !MAX_PIP for particles, FCHAINC for Force chain, NGluedParticles for Glued particles
                  IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     PC =PC + 1
                     dProcBuf(PC) = VAR(LC1)
                  ENDIF
                  IF(PC==LOCAL_CNT) EXIT
               ENDDO
            ENDIF

! Gather local buffer to root
            CALL desmpi_gatherv(ptype=2)

! Write the data, always preceded by its size in number of bytes
! Write root buffer to file

            IF(myPE == PE_IO) THEN
               WRITE(UNIT_VTP) nbytes_scalar

               IF(gsp_implicit_writing) THEN
                  DO LC1=1, total_n_sphere
                     WRITE(UNIT_VTP) real(drootBuf(LC1))
                  ENDDO
               ELSE
                  DO LC1=1, GLOBAL_CNT
                     WRITE(UNIT_VTP) real(drootBuf(LC1))
                  ENDDO
               ENDIF
            ENDIF

            deallocate (dProcBuf, dRootBuf)

         ENDIF ! end of PASS==WRITE_HEADER

      ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN

! Number of bytes for each scalar field
         nbytes_scalar = INT(LOCAL_CNT,8) * c_sizeof(float)
         IF(gsp_implicit_writing) nbytes_scalar = INT(cur_pe_n_sphere,8) * c_sizeof(float)
! Remove possible white space with underscore
         DO I = 1,LEN_TRIM(VAR_NAME)
            IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
         ENDDO

         IF(PASS==WRITE_HEADER) THEN

! For each scalar, write a tag, with corresponding offset
            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! Prepare the offset for the next field
            VTU_offset = VTU_offset + VTK_HEADER_SIZE + nbytes_scalar


         ELSEIF(PASS==WRITE_DATA) THEN

            IF(gsp_implicit_writing) THEN
               allocate (dProcBuf(cur_pe_n_sphere))
            ELSE
               allocate (dProcBuf(LOCAL_CNT))
            ENDIF

! Pack scalar list in a local buffer before writing in file
            PC = 0
            IF(gsp_implicit_writing) THEN
               ! @renjieke BELONGS_TO_VTK_SUBDOMAIN(LC1)) already filter the arrays, so no need to put any if condition here
               DO LC1 = 1, size(VAR(:),1)
                  IF(size(dProcBuf) == 0) CYCLE
                  dProcBuf(LC1) = VAR(LC1)
                  IF(LC1 == cur_pe_n_sphere) EXIT
               ENDDO
            ELSE
               DO LC1 = 1, VTP_LUB !MAX_PIP for particles, FCHAINC for Force chain, NGluedParticles for Glued particles
                  IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     PC =PC + 1
                     dProcBuf(PC) = VAR(LC1)
                  ENDIF
                  IF(PC==LOCAL_CNT) EXIT
               ENDDO
            ENDIF

! Write the data, always preceded by its size in number of bytes
! Write root buffer to file
            WRITE(UNIT_VTP) nbytes_scalar
            ! no gather operation so write dProcBuf only
            IF(gsp_implicit_writing) THEN
               DO LC1=1, cur_pe_n_sphere
                  IF(size(dProcBuf) == 0) CYCLE
                  WRITE(UNIT_VTP)  real(dProcBuf(LC1))
               ENDDO
            ELSE
               DO LC1=1, LOCAL_CNT
                  WRITE(UNIT_VTP)  real(dProcBuf(LC1))
               ENDDO
            ENDIF

            deallocate (dProcBuf)

         ENDIF

      ENDIF

! Update pvtu file with variable name
      IF(BDIST_IO.AND.MyPE==PE_IO.AND.PASS==WRITE_DATA) THEN
         WRITE(UNIT_PVTP,100) '        <PointArray type="Float32" Name="', &
                              TRIM(VAR_NAME),'" format="appended"  />'
      ENDIF

      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10,ADVANCE='NO')'.'

10    FORMAT(A)
90    FORMAT(A,A,A,I20,A)
100   FORMAT(A,A,A)

      RETURN

      END SUBROUTINE WRITE_SCALAR_IN_VTP_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VECTOR_IN_VTP                                    C
!  Purpose: Write Vector variable in a vtp file                        C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VECTOR_IN_VTP_BIN(VAR_NAME,VAR,PASS,cur_pe_n_sphere,total_n_sphere)
   ! @renjieke, local_cnt and global_cnt --> they are cur_n_sphere and total_n_sphere?

      IMPLICIT NONE

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, INTENT(in) :: VAR(:,:)

      REAL(c_float) :: float

      INTEGER(KIND=8) :: nbytes_vector

      INTEGER :: PASS
      INTEGER, INTENT(IN),OPTIONAL :: total_n_sphere,cur_pe_n_sphere
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)  ! local
      DOUBLE PRECISION, ALLOCATABLE :: gtemp_array(:,:)  ! global

      INTEGER :: LB, UB
      INTEGER :: PC, LC1, LC2
      LOGICAL :: gsp_implicit_writing
      INTEGER :: TMP_LOCAL_CNT, TMP_GLOBAL_CNT

      gsp_implicit_writing = present(cur_pe_n_sphere) .and. present(total_n_sphere)

      IF(gsp_implicit_writing) THEN
         TMP_LOCAL_CNT = cur_pe_n_sphere
         TMP_GLOBAL_CNT = total_n_sphere
      ELSE
         TMP_LOCAL_CNT = LOCAL_CNT
         TMP_GLOBAL_CNT = GLOBAL_CNT
      ENDIF

      IF (.NOT.BDIST_IO) THEN
! Number of bytes for each vector field
         nbytes_vector = INT(GLOBAL_CNT,8) * 3 * c_sizeof(float)
         IF(present(total_n_sphere)) nbytes_vector = INT(total_n_sphere,8) * 3 * c_sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
! For each vector, write a tag, with corresponding offset

            IF (myPE == PE_IO) THEN
               WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                    TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
            ENDIF

! Prepare the offset for the next field
            VTU_offset = VTU_offset + VTK_HEADER_SIZE + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            LB = LBOUND(VAR,2) ! This should always be 1
            UB = UBOUND(VAR,2) ! This should always be 2

            ALLOCATE (dProcBuf(TMP_LOCAL_CNT) )
            ALLOCATE (ltemp_array((UB-LB)+1,TMP_LOCAL_CNT))

            IF(myPE == PE_IO) THEN
               ALLOCATE (dRootBuf(TMP_GLOBAL_CNT))
               ALLOCATE (gtemp_array((UB-LB)+1,TMP_GLOBAL_CNT))
            ELSE
               ALLOCATE (dRootBuf(10))
               ALLOCATE (gtemp_array((UB-LB)+1,10))
            ENDIF

! For each vector component, pack component list in a local array
            PC = 0
            IF(gsp_implicit_writing) THEN
               DO LC1 = 1, size(VAR,1)
                  IF(size(ltemp_array,2) == 0) CYCLE
                  DO LC2 = LB, UB
                     ltemp_array(LC2,LC1) = VAR(LC1,LC2)
                  ENDDO
                  IF(LC1 == TMP_LOCAL_CNT) EXIT
               ENDDO
            ELSE
               DO LC1 = 1, VTP_LUB !MAX_PIP for particles, FCHAINC for Force chain, NGluedParticles for Glued particles
                  IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     PC =PC + 1
                     DO LC2=LB, UB
                        ltemp_array(LC2,PC) = VAR(LC1,LC2)
                     ENDDO
                  ENDIF
                  IF(PC==TMP_LOCAL_CNT) EXIT
               ENDDO
            ENDIF

! For each component, gather the local list to global temporary array
         DO LC1 = LB, UB
            dprocbuf(1:TMP_LOCAL_CNT)=ltemp_array(LC1,1:TMP_LOCAL_CNT)
            CALL desmpi_gatherv(ptype=2)
            gtemp_array(LC1,:) = drootbuf(:)
         ENDDO

! Write the data, always preceded by its size in number of bytes
         IF(myPE == PE_IO) THEN
            WRITE(UNIT_VTP) nbytes_vector
            DO LC1=1, TMP_GLOBAL_CNT
               DO LC2=LB, UB
                  WRITE(UNIT_VTP) real(gtemp_array(LC2,LC1))
               ENDDO
            ENDDO
         ENDIF

         deallocate (dProcBuf, dRootBuf, ltemp_array,gtemp_array)

         ENDIF

      ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN

! Number of bytes for each vector field
         nbytes_vector = INT(TMP_LOCAL_CNT,8) * 3 * c_sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
! For each vector, write a tag, with corresponding offset

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! Prepare the offset for the next field
            VTU_offset = VTU_offset + VTK_HEADER_SIZE + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            LB = LBOUND(VAR,2) ! This should always be 1
            UB = UBOUND(VAR,2) ! This should always be 2

            ALLOCATE (ltemp_array((UB-LB)+1,TMP_LOCAL_CNT))

! For each vector component, pack component list in a local array
            PC = 0
            IF(gsp_implicit_writing) THEN
               DO LC1 = 1, size(VAR,1)
                  IF(size(ltemp_array,2) == 0) CYCLE
                  DO LC2 = LB, UB
                     ltemp_array(LC2,LC1) = VAR(LC1,LC2)
                  ENDDO
                  IF(LC1 == TMP_LOCAL_CNT) EXIT
               ENDDO
            ELSE
               DO LC1 = 1, VTP_LUB !MAX_PIP for particles, FCHAINC for Force chain, NGluedParticles for Glued particles
                  IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     PC =PC + 1
                     DO LC2=LB, UB
                        ltemp_array(LC2,PC) = VAR(LC1,LC2)
                     ENDDO
                  ENDIF
                  IF(PC==TMP_LOCAL_CNT) EXIT
               ENDDO
            ENDIF


! Write the data, always preceded by its size in number of bytes
            WRITE(UNIT_VTP) nbytes_vector
            DO LC1=1, TMP_LOCAL_CNT
               DO LC2=LB, UB
                  WRITE(UNIT_VTP)  real(ltemp_array(LC2,LC1))
               ENDDO
            ENDDO

            deallocate (ltemp_array)

         ENDIF

      ENDIF

! Update pvtu file with variable name
      IF(BDIST_IO.AND.MyPE==PE_IO.AND.PASS==WRITE_DATA) THEN
         WRITE(UNIT_PVTP,100) '        <PointArray type="Float32" Name="', &
                              TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended"  />'
      ENDIF

      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10,ADVANCE='NO')'.'

10    FORMAT(A)
90    FORMAT(A,A,A,I20,A)
100   FORMAT(A,A,A)

      RETURN

      END SUBROUTINE WRITE_VECTOR_IN_VTP_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_QUATERNION_IN_VTP                                C
!  Purpose: Write quaternion in a vtp file                             C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 13-Dec-22  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_QUATERNION_IN_VTP_BIN(VAR_NAME,VAR,PASS)

      IMPLICIT NONE

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, INTENT(in) :: VAR(:,:)

      REAL(c_float) :: float

      INTEGER(KIND=8) :: nbytes_vector

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)  ! local
      DOUBLE PRECISION, ALLOCATABLE :: gtemp_array(:,:)  ! global

      INTEGER :: LB, UB
      INTEGER :: PC, LC1, LC2

      IF (.NOT.BDIST_IO) THEN

! Number of bytes for each vector field
         nbytes_vector = INT(GLOBAL_CNT,8) * 4 * c_sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
! For each vector, write a tag, with corresponding offset

            IF (myPE == PE_IO) THEN
               WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                    TRIM(VAR_NAME),'"  NumberOfComponents="4" format="appended" offset="',VTU_offset,'" />'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
            ENDIF

! Prepare the offset for the next field
            VTU_offset = VTU_offset + VTK_HEADER_SIZE + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            LB = LBOUND(VAR,2) ! This should always be 1
            UB = UBOUND(VAR,2) ! This should always be 4

            ALLOCATE (dProcBuf(LOCAL_CNT) )
            ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))

            IF(myPE == PE_IO) THEN
               ALLOCATE (dRootBuf(GLOBAL_CNT))
               ALLOCATE (gtemp_array((UB-LB)+1,GLOBAL_CNT))
            ELSE
               ALLOCATE (dRootBuf(10))
               ALLOCATE (gtemp_array((UB-LB)+1,10))
            ENDIF

! For each vector component, pack component list in a local array
            PC = 0
            DO LC1 = 1, VTP_LUB !MAX_PIP for particles, FCHAINC for Force chain, NGluedParticles for Glued particles
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = VAR(LC1,LC2)
                  ENDDO
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO


! For each component, gather the local list to global temporary array
         DO LC1 = LB, UB
            dprocbuf(1:LOCAL_CNT)=ltemp_array(LC1,1:LOCAL_CNT)
            CALL desmpi_gatherv(ptype=2)
            gtemp_array(LC1,:) = drootbuf(:)
         ENDDO

! Write the data, always preceded by its size in number of bytes
         IF(myPE == PE_IO) THEN
            WRITE(UNIT_VTP) nbytes_vector
            DO LC1=1, GLOBAL_CNT
               DO LC2=LB, UB
                  WRITE(UNIT_VTP)  real(gtemp_array(LC2,LC1))
               ENDDO
            ENDDO
         ENDIF

         deallocate (dProcBuf, dRootBuf, ltemp_array,gtemp_array)

         ENDIF

      ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN

! Number of bytes for each vector field
         nbytes_vector = INT(LOCAL_CNT,8) * 4 * c_sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
! For each vector, write a tag, with corresponding offset

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="4" format="appended" offset="',VTU_offset,'" />'
            WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

! Prepare the offset for the next field
            VTU_offset = VTU_offset + VTK_HEADER_SIZE + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            LB = LBOUND(VAR,2) ! This should always be 1
            UB = UBOUND(VAR,2) ! This should always be 2

            ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))

! For each vector component, pack component list in a local array
            PC = 0
            DO LC1 = 1, VTP_LUB !MAX_PIP for particles, FCHAINC for Force chain, NGluedParticles for Glued particles
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = VAR(LC1,LC2)
                  ENDDO
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO


! Write the data, always preceded by its size in number of bytes
            WRITE(UNIT_VTP) nbytes_vector
            DO LC1=1, LOCAL_CNT
               DO LC2=LB, UB
                  WRITE(UNIT_VTP)  real(ltemp_array(LC2,LC1))
               ENDDO
            ENDDO

            deallocate (ltemp_array)


         ENDIF

      ENDIF

! Update pvtu file with variable name
      IF(BDIST_IO.AND.MyPE==PE_IO.AND.PASS==WRITE_DATA) THEN
         WRITE(UNIT_PVTP,100) '        <PointArray type="Float32" Name="', &
                              TRIM(VAR_NAME),'"  NumberOfComponents="4" format="appended"  />'
      ENDIF

      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10,ADVANCE='NO')'.'

10    FORMAT(A)
90    FORMAT(A,A,A,I20,A)
100   FORMAT(A,A,A)

      RETURN

      END SUBROUTINE WRITE_QUATERNION_IN_VTP_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VECTOR_IN_VTP_BIN9                               C
!  Purpose: Write SuperDEM and GluedSphere transformation matrix in a  C
!           vtp file                                                   C
!           Binary format                                              C
!                                                                      C
!  Author: Xi Gao                                     Date: 25-Feb-19  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

  SUBROUTINE WRITE_VECTOR_IN_VTP_BIN9(VAR_NAME,VAR,PASS)

      IMPLICIT NONE

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, INTENT(in) :: VAR(:,:)

      REAL(c_float) :: float

      INTEGER(KIND=8) :: nbytes_vector

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)  ! local
      DOUBLE PRECISION, ALLOCATABLE :: gtemp_array(:,:)  ! global
      DOUBLE PRECISION :: qc(4), qrot(3,3),qrott(9)
      INTEGER :: LB, UB,UBB
      INTEGER :: PC, LC1, LC2

      UBB=9

      IF (.NOT.BDIST_IO) THEN

! Number of bytes for each vector field
         nbytes_vector = INT(GLOBAL_CNT,8) * 9 * c_sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
! For each vector, write a tag, with corresponding offset

            IF (myPE == PE_IO) THEN
               WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                    TRIM(VAR_NAME),'"  NumberOfComponents="9" format="appended" offset="',VTU_offset,'" />'
               WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC
            ENDIF

! Prepare the offset for the next field
            VTU_offset = VTU_offset + VTK_HEADER_SIZE + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            LB = LBOUND(VAR,2) ! This should always be 1
            UB = UBOUND(VAR,2) ! This should always be 2

            ALLOCATE (dProcBuf(LOCAL_CNT) )
            ALLOCATE (dRootBuf(GLOBAL_CNT))
            ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))
            ALLOCATE (gtemp_array((UB-LB)+1,GLOBAL_CNT))

! For each vector component, pack component list in a local array
            PC = 0
            DO LC1 = 1, MAX_PIP
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = VAR(LC1,LC2)
                  ENDDO
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO


! For each component, gather the local list to global temporary array
         DO LC1 = LB, UB
            dprocbuf(1:LOCAL_CNT)=ltemp_array(LC1,1:LOCAL_CNT)
            CALL desmpi_gatherv(ptype=2)
            gtemp_array(LC1,:) = drootbuf(:)
         ENDDO

! Write the data, always preceded by its size in number of bytes
         IF(myPE == PE_IO) THEN
            WRITE(UNIT_VTP) nbytes_vector
            DO LC1=1, GLOBAL_CNT
!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 Qc(1)=gtemp_array(1,LC1)
                 Qc(2)=gtemp_array(2,LC1)
                 Qc(3)=gtemp_array(3,LC1)
                 Qc(4)=gtemp_array(4,LC1)
                 call  QROTATIONMATRIX(Qc,Qrot)

                 Qrott(1)=2.0D0*Qrot(1,1)
                 Qrott(2)=2.0D0*Qrot(1,2)
                 Qrott(3)=2.0D0*Qrot(1,3)
                 Qrott(4)=2.0D0*Qrot(2,1)
                 Qrott(5)=2.0D0*Qrot(2,2)
                 Qrott(6)=2.0D0*Qrot(2,3)
                 Qrott(7)=2.0D0*Qrot(3,1)
                 Qrott(8)=2.0D0*Qrot(3,2)
                 Qrott(9)=2.0D0*Qrot(3,3)

                 DO LC2=LB, UBB
                  WRITE(UNIT_VTP)  real(Qrott(LC2))
               ENDDO
            ENDDO
         ENDIF

         deallocate (dProcBuf, dRootBuf, ltemp_array,gtemp_array)

         ENDIF

         ! not implemented yet!
     ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN
         ! Number of bytes for each vector field
         nbytes_vector = INT(LOCAL_CNT,8) * 3 * c_sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
             ! For each vector, write a tag, with corresponding offset

             WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
             WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

             ! Prepare the offset for the next field
             VTU_offset = VTU_offset + VTK_HEADER_SIZE + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

             LB = LBOUND(VAR,1) ! This should always be 1
             UB = UBOUND(VAR,1) ! This should always be 2

             ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))

             ! For each vector component, pack component list in a local array
             PC = 0
             DO LC1 = 1, MAX_PIP
                 IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                     PC =PC + 1
                     DO LC2=LB, UB
                         ltemp_array(LC2,PC) = VAR(LC2,LC1)
                     ENDDO
                 ENDIF
                 IF(PC==LOCAL_CNT) EXIT
             ENDDO


             ! Write the data, always preceded by its size in number of bytes
             WRITE(UNIT_VTP) nbytes_vector
            DO LC1=1, LOCAL_CNT
                DO LC2=LB, UB
                    WRITE(UNIT_VTP)  real(ltemp_array(LC2,LC1))
                ENDDO
            ENDDO

            deallocate (ltemp_array)

            IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
                WRITE(UNIT_PVTP,90)'        <PointArray type="Float32" Name="', &
                    TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            ENDIF

        ENDIF

    ENDIF

    IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10,ADVANCE='NO')'.'

10  FORMAT(A)
90  FORMAT(A,A,A,I20,A)

    RETURN

END SUBROUTINE WRITE_VECTOR_IN_VTP_BIN9


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLOSE_VTP_FILE_BIN                                     C
!  Purpose: Close a vtp file                                           C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLOSE_VTP_FILE_BIN(MODE)

      IMPLICIT NONE

      INTEGER:: N
      CHARACTER (LEN=256) :: VTU_NAME
      INTEGER, DIMENSION(0:numPEs-1) :: ALL_PART_CNT
      INTEGER :: IERR
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)

      IF((myPE == PE_IO.AND.(.NOT.BDIST_IO)).OR.(BDIST_IO.AND.LOCAL_CNT>0)) THEN

! Write last tags and close the vtp file
      WRITE(BUFFER, *)'  </AppendedData>'
      WRITE(UNIT_VTP)END_REC//TRIM(BUFFER)//END_REC

      WRITE(BUFFER, *)'</VTKFile>'
      WRITE(UNIT_VTP)TRIM(BUFFER)//END_REC

      CLOSE(UNIT_VTP)
      IF (FULL_LOG) WRITE(*,10)''

10    FORMAT(A)

      ENDIF

! Update pvtu file and close

      IF(BDIST_IO)  THEN
         CALL allgather_1i (LOCAL_CNT,ALL_PART_CNT,IERR)

         IF (myPE == PE_IO.AND.GLOBAL_CNT>0) THEN
            WRITE(UNIT_PVTP, *) '      </PPointData>'

            DO N = 0,NumPEs-1
               IF(ALL_PART_CNT(N)>0) THEN
! The pvtp and vtp files are in the same directory (either project dir or VTU_DIR)
! The VTU_DIR should not be prepended to the basename (see .FALSE. as last argument)
                  VTU_NAME = FIND_VTK_BASENAME(MODE, N, .FALSE.) // ".vtp"
                  WRITE(UNIT_PVTP,110) '      <Piece Source="',TRIM(VTU_NAME),'"/>'
               ENDIF
            ENDDO

            WRITE(UNIT_PVTP, *) '  </PPolyData>'
            WRITE(UNIT_PVTP, *) '</VTKFile>'
            CLOSE(UNIT_PVTP)
         ENDIF
      ENDIF

110   FORMAT(A,A,A)

      RETURN

      END SUBROUTINE CLOSE_VTP_FILE_BIN


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SETUP_VTK_REGION_PARTICLES                             C
!                                                                      C
!  Purpose: Filter the particles  based on the VTK region bounds and   C
!           set the flag BELONGS_TO_VTK_SUBDOMAIN to .TRUE.            C
!           to keep the particle.                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SETUP_VTK_REGION_PARTICLES

      IMPLICIT NONE

      INTEGER :: PC,LC1,M
      INTEGER :: NXS,NYS,NZS,NS
      DOUBLE PRECISION :: X_SLICE(DIM_I),Y_SLICE(DIM_J),Z_SLICE(DIM_K)
      DOUBLE PRECISION :: XE,XW,YS,YN,ZB,ZT
      DOUBLE PRECISION :: XP,YP,ZP,XP1,YP1,ZP1,XP2,YP2,ZP2,R

      DOUBLE PRECISION :: SLICE_TOL
      LOGICAL :: KEEP_XDIR,KEEP_YDIR,KEEP_ZDIR

! Variables related to gather
      integer lgathercnts(0:numpes-1), lproc

      CHARACTER(LEN=1) :: SELECT_PARTICLE_BY

! Get VTK region bounds
      XE = VTK_X_E(VTK_REGION)
      XW = VTK_X_W(VTK_REGION)
      YS = VTK_Y_S(VTK_REGION)
      YN = VTK_Y_N(VTK_REGION)
      ZB = VTK_Z_B(VTK_REGION)
      ZT = VTK_Z_T(VTK_REGION)

      NXS = VTK_NXS(VTK_REGION)
      NYS = VTK_NYS(VTK_REGION)
      NZS = VTK_NZS(VTK_REGION)

      SLICE_TOL = VTK_SLICE_TOL(VTK_REGION)

      SELECT_PARTICLE_BY = VTK_SELECT_MODE(VTK_REGION)

! get slice(s) location
      DO NS = 1,NXS
         X_SLICE(NS) = XW + (XE-XW)/FLOAT(NXS-1)*FLOAT(NS-1)
      ENDDO

      DO NS = 1,NYS
         Y_SLICE(NS) = YS + (YN-YS)/FLOAT(NYS-1)*FLOAT(NS-1)
      ENDDO

      DO NS = 1,NZS
         Z_SLICE(NS) = ZB + (ZT-ZB)/FLOAT(NZS-1)*FLOAT(NS-1)
      ENDDO


! Loop through all particles on local rank and keep a list of particles
! belonging to VTK region

      IF(ALLOCATED(BELONGS_TO_VTK_SUBDOMAIN)) DEALLOCATE(BELONGS_TO_VTK_SUBDOMAIN)
      ALLOCATE(BELONGS_TO_VTK_SUBDOMAIN(VTP_LUB))

      BELONGS_TO_VTK_SUBDOMAIN = .FALSE.

      LOCAL_CNT = 0
      PC = 1
      DO LC1 = 1, VTP_LUB !MAX_PIP for particles, FCHAINC for Force chain, NGluedParticles for Glued particles

         IF(VTP_P_DATA) THEN ! Particle data
            ! IF(.NOT.IS_NORMAL(LC1)) CYCLE
            IF(.NOT.(IS_NORMAL(LC1).OR.IS_RIGID_MOTION(LC1))) CYCLE
            IF(PC > PIP) EXIT
            ! IF(IS_NONEXISTENT(LC1)) CYCLE
            PC = PC+1
            ! IF(IS_GHOST(LC1) .OR. IS_ENTERING_GHOST(LC1) .OR. IS_EXITING_GHOST(LC1)) CYCLE

            M = PIJK(LC1,5)

            IF(.NOT.VTK_PART_PHASE(VTK_REGION,M)) CYCLE

               XP = DES_POS_NEW(LC1,1)
               YP = DES_POS_NEW(LC1,2)
               ZP = DES_POS_NEW(LC1,3)

               R = DES_RADIUS(LC1)

               XP1 = DES_POS_NEW(LC1,1) - R
               YP1 = DES_POS_NEW(LC1,2) - R
               ZP1 = DES_POS_NEW(LC1,3) - R

               XP2 = DES_POS_NEW(LC1,1) + R
               YP2 = DES_POS_NEW(LC1,2) + R
               ZP2 = DES_POS_NEW(LC1,3) + R

         ELSEIF (VTP_GLUE_DATA) THEN ! Glued Particle
	            XP = gluePos(LC1,1)
               YP = gluePos(LC1,2)
               ZP = gluePos(LC1,3)

               R = glueBounding(LC1)/2.0

               XP1 = gluePos(LC1,1) - R
               YP1 = gluePos(LC1,2) - R
               ZP1 = gluePos(LC1,3) - R

               XP2 = gluePos(LC1,1) + R
               YP2 = gluePos(LC1,2) + R
               ZP2 = gluePos(LC1,3) + R
	 ELSE ! Force chain
              ! Only selection mode available is 'C' (Force chain midpoint must be inside vtk region)
               XP = FCHAIN_MIDPOINT(LC1,1)
               YP = FCHAIN_MIDPOINT(LC1,2)
               ZP = FCHAIN_MIDPOINT(LC1,3)

               IF(VTK_FCMAX_ONLY(VTK_REGION)) THEN ! Only keep the max FCHAIN
                                                   ! setting the x coordinate to
                                                   ! -UNDEFINED will exclude it
                                                   ! in the next step
                  IF(.NOT.FCHAIN_FCMAX(LC1)) XP = -UNDEFINED
               ENDIF

         ENDIF

         SELECT CASE(SELECT_PARTICLE_BY)
            CASE('C')  ! Particle center must be inside vtk region

! X-direction
               KEEP_XDIR=.FALSE.
               IF(NXS==0) THEN
                  IF(XW<=XP.AND.XP<=XE) KEEP_XDIR=.TRUE.
               ELSE
                  DO NS = 1,NXS
                     IF((X_SLICE(NS)-SLICE_TOL)<=XP.AND.XP<=(X_SLICE(NS)+SLICE_TOL)) KEEP_XDIR=.TRUE.
                  ENDDO
               ENDIF

! Y-direction
               KEEP_YDIR=.FALSE.
               IF(NYS==0) THEN
                  IF(YS<=YP.AND.YP<=YN) KEEP_YDIR=.TRUE.
               ELSE
                  DO NS = 1,NYS
                     IF((Y_SLICE(NS)-SLICE_TOL)<=YP.AND.YP<=(Y_SLICE(NS)+SLICE_TOL)) KEEP_YDIR=.TRUE.
                  ENDDO
               ENDIF

! Z-direction
               KEEP_ZDIR=.FALSE.
               IF(NZS==0) THEN
                  IF(ZB<=ZP.AND.ZP<=ZT) KEEP_ZDIR=.TRUE.
               ELSE
                  DO NS = 1,NZS
                     IF((Z_SLICE(NS)-SLICE_TOL)<=ZP.AND.ZP<=(Z_SLICE(NS)+SLICE_TOL)) KEEP_ZDIR=.TRUE.
                  ENDDO
               ENDIF


            CASE('P')  ! Entire particle must be inside vtk region

! X-direction
               KEEP_XDIR=.FALSE.
               IF(NXS==0) THEN
                  IF(XW<=XP1.AND.XP2<=XE) KEEP_XDIR=.TRUE.
               ELSE
                  DO NS = 1,NXS
                     IF((X_SLICE(NS)-SLICE_TOL)<=XP1.AND.XP2<=(X_SLICE(NS)+SLICE_TOL)) KEEP_XDIR=.TRUE.
                  ENDDO
               ENDIF

! Y-direction
               KEEP_YDIR=.FALSE.
               IF(NYS==0) THEN
                  IF(YS<=YP1.AND.YP2<=YN) KEEP_YDIR=.TRUE.
               ELSE
                  DO NS = 1,NYS
                     IF((Y_SLICE(NS)-SLICE_TOL)<=YP1.AND.YP2<=(Y_SLICE(NS)+SLICE_TOL)) KEEP_YDIR=.TRUE.
                  ENDDO
               ENDIF

! Z-direction
               KEEP_ZDIR=.FALSE.
               IF(NZS==0) THEN
                  IF(ZB<=ZP1.AND.ZP2<=ZT) KEEP_ZDIR=.TRUE.
               ELSE
                  DO NS = 1,NZS
                     IF((Z_SLICE(NS)-SLICE_TOL)<=ZP1.AND.ZP2<=(Z_SLICE(NS)+SLICE_TOL)) KEEP_ZDIR=.TRUE.
                  ENDDO
               ENDIF


            CASE('I')  ! Particle must be inside or intersect the edge of the vtk region

! X-direction
               KEEP_XDIR=.FALSE.
               IF(NXS==0) THEN
                  IF(.NOT.(XE<=XP1.OR.XP2<=XW)) KEEP_XDIR=.TRUE.
               ELSE
                  DO NS = 1,NXS
                     IF(.NOT.((X_SLICE(NS)+SLICE_TOL)<=XP1.OR.XP2<=(X_SLICE(NS)-SLICE_TOL))) KEEP_XDIR=.TRUE.
                  ENDDO
               ENDIF

! Y-direction
               KEEP_YDIR=.FALSE.
               IF(NYS==0) THEN
                  IF(.NOT.(YN<=YP1.OR.YP2<=YS)) KEEP_YDIR=.TRUE.
               ELSE
                  DO NS = 1,NYS
                     IF(.NOT.((Y_SLICE(NS)+SLICE_TOL)<=YP1.OR.YP2<=(Y_SLICE(NS)-SLICE_TOL))) KEEP_YDIR=.TRUE.
                  ENDDO
               ENDIF

! Z-direction
               KEEP_ZDIR=.FALSE.
               IF(NZS==0) THEN
                  IF(.NOT.(ZT<=ZP1.OR.ZP2<=ZB)) KEEP_ZDIR=.TRUE.
               ELSE
                  DO NS = 1,NZS
                     IF(.NOT.((Z_SLICE(NS)+SLICE_TOL)<=ZP1.OR.ZP2<=(Z_SLICE(NS)-SLICE_TOL))) KEEP_ZDIR=.TRUE.
                  ENDDO
               ENDIF

            CASE DEFAULT
               print*,'should not be here (select particle by)'
         END SELECT

! Now combine
         IF(KEEP_XDIR.AND.KEEP_YDIR.AND.KEEP_ZDIR) THEN
            BELONGS_TO_VTK_SUBDOMAIN(LC1) = .TRUE.
            LOCAL_CNT = LOCAL_CNT + 1
         ENDIF
      ENDDO ! particle loop

! Calculate the total number of particles system-wide.
      GLOBAL_CNT=10
      call global_sum(LOCAL_CNT, GLOBAL_CNT)

! No need to set the send/reccv when using distributed IO
      IF (BDIST_IO) RETURN
! Set the send count from the local process.
      igath_sendcnt = LOCAL_CNT

! Collect the number of particles on each rank.all ranks.
      lgathercnts = 0
      lgathercnts(myPE) = LOCAL_CNT
      call global_sum(lgathercnts,igathercnts)

! Calculate the rank displacements.
      idispls(0) = 0
      DO lPROC = 1,NUMPEs-1
         idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
      ENDDO

      RETURN

   END SUBROUTINE SETUP_VTK_REGION_PARTICLES

   SUBROUTINE  ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN(VAR_NAME, VAR, PASS)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: VAR_NAME
      DOUBLE PRECISION, INTENT(in) :: VAR(:)
      INTEGER, INTENT(in) :: PASS

      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2
      INTEGER :: L, P1, P2
      CHARACTER(LEN=128) :: VAR_NAME1, VAR_NAME2
      DOUBLE PRECISION, ALLOCATABLE :: VAR1(:), VAR2(:)

! Below 1 and 2 refers to first and second particles in a contact pair.
! 1 represents particle LL , and 2 represents particle I in calc_force_dem

! Append index to variable name
      VAR_NAME1 = TRIM(VAR_NAME)//'_1'
      VAR_NAME2 = TRIM(VAR_NAME)//'_2'

! Get data for pair of particles (P1, P2)
      IF(PASS==WRITE_DATA) THEN
         ALLOCATE(VAR1(FCHAINC))
         ALLOCATE(VAR2(FCHAINC))
         DO L = 1, FCHAINC
            P1 = FCHAIN_LOCAL_ID1(L)
            P2 = FCHAIN_LOCAL_ID2(L)
            VAR1(L) = VAR(P1)
            VAR2(L) = VAR(P2)
         ENDDO
      ENDIF

! Write data for pair of particles (VAR1, VAR2)
      CALL WRITE_SCALAR_IN_VTP_BIN(VAR_NAME1,VAR1,PASS)
      CALL WRITE_SCALAR_IN_VTP_BIN(VAR_NAME2,VAR2,PASS)

      IF(PASS==WRITE_DATA) THEN
         DEALLOCATE(VAR1)
         DEALLOCATE(VAR2)
      ENDIF

      RETURN
   END SUBROUTINE  ADD_SCALAR_TO_FORCE_CHAIN_VTP_BIN

   SUBROUTINE  ADD_VECTOR_TO_FORCE_CHAIN_VTP_BIN(VAR_NAME, VAR, PASS)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in) :: VAR_NAME
      DOUBLE PRECISION, INTENT(in) :: VAR(:,:)
      INTEGER, INTENT(in) :: PASS

      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2
      INTEGER :: L, P1, P2
      CHARACTER(LEN=128) :: VAR_NAME1, VAR_NAME2
      DOUBLE PRECISION, ALLOCATABLE :: VAR1(:,:), VAR2(:,:)

! Below 1 and 2 refers to first and second particles in a contact pair.
! 1 represents particle LL , and 2 represents particle I in calc_force_dem

! Append index to variable name
      VAR_NAME1 = TRIM(VAR_NAME)//'_1'
      VAR_NAME2 = TRIM(VAR_NAME)//'_2'

! Get data for pair of particles (P1, P2)
      IF(PASS==WRITE_DATA) THEN
         ALLOCATE(VAR1(FCHAINC,3))
         ALLOCATE(VAR2(FCHAINC,3))
         DO L = 1, FCHAINC
            P1 = FCHAIN_LOCAL_ID1(L)
            P2 = FCHAIN_LOCAL_ID2(L)
            VAR1(L,:) = VAR(P1,:)
            VAR2(L,:) = VAR(P2,:)
         ENDDO
      ENDIF

! Write data for pair of particles (VAR1, VAR2)
      CALL WRITE_VECTOR_IN_VTP_BIN(VAR_NAME1,VAR1,PASS)
      CALL WRITE_VECTOR_IN_VTP_BIN(VAR_NAME2,VAR2,PASS)

      IF(PASS==WRITE_DATA) THEN
         DEALLOCATE(VAR1)
         DEALLOCATE(VAR2)
      ENDIF

      RETURN
   END SUBROUTINE  ADD_VECTOR_TO_FORCE_CHAIN_VTP_BIN

END MODULE VTP
