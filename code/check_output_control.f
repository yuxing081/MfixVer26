#include "error.inc"

MODULE CHECK_OUTPUT_CONTROL_MOD

   use error_manager

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_OUTPUT_CONTROL                                    !
!  Purpose: Check the output control namelist section                  !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_OUTPUT_CONTROL

! Global Variables:
!---------------------------------------------------------------------//
! Time intervalue between updating the RES and SPx files.
      use output, only: RES_DT
! Time-step intervalue between updating the .LOG file.
      use output, only: NLOG
! Screen verbosity
      use output, only: FULL_LOG
! VTK
      use vtk

! Global Parameters:
!---------------------------------------------------------------------//
! Number aliases
      use param1, only: UNDEFINED, ZERO

! Skip data check when doing preprocessing only
      USE run, only:ppo

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//

!......................................................................!

      IF(PPO) RETURN

! Check the values specified for the RES file.
      IF (RES_DT==UNDEFINED)THEN
         WRITE(ERR_MSG,1000) 'RES_DT'
         CALL LOG_ERROR()
      ELSEIF(RES_DT <= ZERO) THEN
         WRITE(ERR_MSG,1002) 'RES_DT', RES_DT
         CALL LOG_ERROR()
      ENDIF

!     Verify that the LOG frequency is valid.
      IF(NLOG <= 0) THEN
         WRITE(ERR_MSG,1003) 'NLOG', NLOG
         CALL LOG_ERROR()
      ENDIF


! Screen verbosity

      CALL CHECK_OUTPUT_SPX()
      CALL CHECK_OUTPUT_VTK()

! Monitor data check is now called in INITIALIZE subroutine (main.f)
! Monitor data check needs to be called after variable allocation so
! invalid monitors can be turned off when variables are not used in a simulation
!      CALL CHECK_OUTPUT_MONITOR()

! New data check for reaction rates only
      CALL CHECK_MONITOR_RRATES

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A)

 1002 FORMAT('Error 1002: Invalid or unknown input: ',A,' = ',E14.6)

 1003 FORMAT('Error 1003: Invalid or unknown input: ',A,' = ',I4)

 2010 FORMAT('Error 2010: Invalid or unknown input: ',A,I4)
      END SUBROUTINE CHECK_OUTPUT_CONTROL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_OUTPUT_CONTROL                                    !
!  Purpose: Check the output control namelist section                  !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_OUTPUT_SPX

! Global Variables:
!---------------------------------------------------------------------//
! Time intervalue between updating the RES and SPx files.
! Flag: Use the K-Epsilon model
      use turb, only: K_EPSILON
      use output, only: SPX_DT
! Number of arrays to store in SPA
      use rxns, only: nRR

! Global Parameters:
!---------------------------------------------------------------------//
! Number aliases
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, LARGE_NUMBER
! Number of SPx files.
      USE param1, only: N_SPX

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: LC

! Check the SPx Files
      SPx_LP: DO LC = 1, N_SPX

! Disable writing the .SPA file if nRR is unspecified.
         IF(LC == 10) THEN
            IF(nRR == 0) THEN
               IF (SPX_DT(LC) == UNDEFINED) SPX_DT(LC) = LARGE_NUMBER
               CYCLE SPx_LP
            ENDIF

! Disable writing the .SPB file if K-Epsilon is unspecified.
         ELSEIF(LC == 11) THEN
            IF(.NOT.K_Epsilon) THEN
               IF (SPX_DT(LC)==UNDEFINED)SPX_DT(LC) = LARGE_NUMBER
               CYCLE SPx_LP
            ENDIF

! Verify the remaining SPx files.
         ELSE
            IF(SPX_DT(LC) <= ZERO) THEN
               WRITE(ERR_MSG,1001) iVar('SPX_DT',LC), SPX_DT(LC)
               CALL LOG_ERROR()
            ENDIF
         ENDIF
      ENDDO SPx_LP

      RETURN

 1001 FORMAT('Error 1001: Invalid or unknown input: ',A,' = ',A)

      END SUBROUTINE CHECK_OUTPUT_SPX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_OUTPUT_VTK                                        !
!  Purpose: Check the output control namelist section                  !
!                                                                      !
!  Author: P. Nicoletti                               Date: 27-NOV-91  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_OUTPUT_VTK

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Use the K-Epsilon model
      use turb, only: K_EPSILON
! Number of arrays to store in SPA
      use rxns, only: nRR
! VTK
      use vtk
      USE run, only: RUN_NAME, REINITIALIZING, DT, TIME
      USE physprop, only: MMAX
      USE scalars, only :NSCALAR
      USE DISCRETELEMENT, only:DISCRETE_ELEMENT
      USE DISCRETELEMENT, only: PARTICLE_ORIENTATION
      USE DISCRETELEMENT, only: DES_USR_VAR_SIZE
      USE DISCRETELEMENT, ONLY: DES_MMAX
      USE DISCRETELEMENT, ONLY: WRITE_FORCE_CHAIN
      USE cutcell, only: USE_STL
      use des_rxns, only: NO_OF_DES_RXNS
      use rxns, only: NO_OF_RXNS
      use param, only: dimension_3
      use parallel, only: is_smp

! Global Parameters:
!---------------------------------------------------------------------//
! Number aliases
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, LARGE_NUMBER
! Number of SPx files.
      USE param1, only: N_SPX
! Geometry bounds
      use geometry, only: X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//

      INTEGER :: L,M,N,LV,N_VTK_REGIONS,R, M_IN_VTK

      CHARACTER(LEN=32) :: SUBM

      LOGICAL :: SMP_WARNING_PRINTED = .FALSE.
!......................................................................!

! Check VTK regions

      IF(FRAME(1)<-1) THEN
         WRITE(ERR_MSG, 2000) trim(iVAL(FRAME(1)))
         CALL LOG_ERROR()
      ENDIF

 2000 FORMAT('Error 2000: Invalid value for FRAME = ',A, '. Valid ',&
         'values',/'are integers >= -1.')

      IF(VTK_DT(1)<ZERO) THEN
          WRITE(ERR_MSG,2001) trim(iVal(VTK_DT(1)))
         CALL LOG_ERROR()
      ENDIF

 2001 FORMAT('Error 2001: Invalid value for VTK_DT = ',A, '. Valid',&
         ' values',/'are positive numbers (e.g., 0.1).')

      N_VTK_REGIONS = 0

      SAVE_FLUID_RRATES = .FALSE.
      SAVE_DES_RRATES = .FALSE.
      SAVE_PART_RRATES = .FALSE.

      DO L = 1, DIMENSION_VTK
         VTK_DEFINED(L) = .FALSE.
         IF (VTK_X_W(L) /= -UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_X_E(L) /=  UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_Y_S(L) /= -UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_Y_N(L) /=  UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_Z_B(L) /= -UNDEFINED)   VTK_DEFINED(L) = .TRUE.
         IF (VTK_Z_T(L) /=  UNDEFINED)   VTK_DEFINED(L) = .TRUE.

         IF (VTK_DATA(L) == 'G')   VTK_DEFINED(L) = .TRUE.

         IF (VTK_DOMAIN_DECOMPOSITION(L)) VTK_DEFINED(L) = .TRUE.

         IF(VTK_DEFINED(L)) THEN
            N_VTK_REGIONS =  N_VTK_REGIONS + 1
         ELSE IF(REINITIALIZING) THEN
            VTK_TIME(L) = UNDEFINED
         ENDIF
         ! If a new vtk region was created while the solver is paused, need
         ! to set the vtk_time.
         IF(VTK_DEFINED(L).AND.VTK_TIME(L) == UNDEFINED) THEN
           VTK_TIME(L) = (INT((TIME + 0.1d0*DT)/VTK_DT(L))+1)*VTK_DT(L)
         ENDIF
      ENDDO   ! end loop over (l = 1,dimension_vtk)

! There must be at least one VTK region defined
! If this is not the case, define the entire domain as default region
      IF(WRITE_VTK_FILES.AND.N_VTK_REGIONS==0) THEN
         VTK_DEFINED(1) = .TRUE.
         VTK_X_W(1) = X_MIN
         VTK_X_E(1) = X_MAX
         VTK_Y_S(1) = Y_MIN
         VTK_Y_N(1) = Y_MAX
         VTK_Z_B(1) = Z_MIN
         VTK_Z_T(1) = Z_MAX
         VTK_FILEBASE(1) = RUN_NAME
      ENDIF

! If VTK_VAR is defined, fill-up the variable list
! for the vtk subdomains
      DO L = 1, DIM_VTK_VAR
         IF(VTK_VAR(L)/=UNDEFINED_I) VTK_VARLIST(:,L) = VTK_VAR(L)
      ENDDO

      WRITE_FORCE_CHAIN = .FALSE.

      DO L = 1, DIMENSION_VTK

         IF(.NOT.VTK_DEFINED(L)) CYCLE

         DO LV = 1, DIM_VTK_VAR

            SELECT CASE (VTK_VARLIST(L,LV))

               CASE (1)
                  VTK_EP_g(L) = .TRUE.

               CASE (2)
                  VTK_P_g(L)        = .TRUE.
                  VTK_P_s(L,1:MMAX) = .TRUE.

               CASE (3)
                  VTK_VEL_G(L) = .TRUE.

               CASE (4)
                  DO M = 1,MMAX
                     VTK_VEL_S(L,M) = .TRUE.
                  END DO

               CASE (5)
                  DO M = 1,MMAX
                     VTK_ROP_s(L,M) = .TRUE.
                  END DO

               CASE (6)
                  VTK_T_g(L) = .TRUE.
                  DO M = 1,MMAX
                     VTK_T_s(L,M) = .TRUE.
                  END DO

               CASE (7)
                 !DO N = 1,NMAX(0)
                    VTK_X_g(L,:) = .TRUE.
                 !END DO

                  DO M = 1, MMAX
                    !DO N = 1,NMAX(M)
                        VTK_X_s(L,M,:) = .TRUE.
                    !END DO
                  END DO

               CASE (8)
                  DO M = 1,MMAX
                     VTK_Theta_m(L,M) = .TRUE.
                  END DO

               CASE (9)
                  DO N = 1,NSCALAR
                     VTK_Scalar(L,N) =.TRUE.
                  END DO

               CASE (10)
                  DO R = 1,nRR
                     VTK_RRate(L,R) = .TRUE.
                  END DO

               CASE (11)
                  IF(K_EPSILON) THEN
                     VTK_K_Turb_G(L) = .TRUE.
                     VTK_E_Turb_G(L) = .TRUE.
                  ENDIF

               CASE (12)
                  VTK_VORTICITY(L) = .TRUE.
                  VTK_LAMBDA_2(L)  = .TRUE.

               CASE (100)
                  VTK_PARTITION(L) = .TRUE.

               CASE (101)
                  VTK_BC_ID(L) = .TRUE.

               CASE (102)
                  VTK_DWALL(L) = .TRUE.

               CASE (103)
                  IF(DISCRETE_ELEMENT.AND.USE_STL) THEN
                     VTK_FACET_COUNT_DES(L) = .TRUE.
                  ENDIF

               CASE (104)
                  IF(DISCRETE_ELEMENT.AND.USE_STL) THEN
                     VTK_NB_FACET_DES(L) = .TRUE.
                  ENDIF

               CASE(999)
                  VTK_IJK(L) = .TRUE.

               CASE(1000)
                  VTK_NORMAL(L) = .TRUE.

               CASE (1001)
                  VTK_DEBUG(L,1) = .TRUE.

               CASE (1002)
                  VTK_DEBUG(L,2) = .TRUE.

               CASE (1003)
                  VTK_DEBUG(L,3) = .TRUE.

               CASE (1004)
                  VTK_DEBUG(L,4) = .TRUE.

               CASE (1005)
                  VTK_DEBUG(L,5) = .TRUE.

               CASE (1006)
                  VTK_DEBUG(L,6) = .TRUE.

               CASE (1007)
                  VTK_DEBUG(L,7) = .TRUE.

               CASE (1008)
                  VTK_DEBUG(L,8) = .TRUE.

               CASE (1009)
                  VTK_DEBUG(L,9) = .TRUE.

               CASE (1010)
                  VTK_DEBUG(L,10) = .TRUE.

               CASE (1011)
                  VTK_DEBUG(L,11) = .TRUE.

               CASE (1012)
                  VTK_DEBUG(L,12) = .TRUE.

               CASE (1013)
                  VTK_DEBUG(L,13) = .TRUE.

               CASE (1014)
                  VTK_DEBUG(L,14) = .TRUE.

               CASE (1015)
                  VTK_DEBUG(L,15) = .TRUE.


               CASE (0) ! do nothing

               CASE (UNDEFINED_I) ! do nothing

               CASE DEFAULT
                  WRITE(ERR_MSG,2100) trim(iVal(L)),                   &
                     trim(iVal(VTK_VAR(L)))
                  CALL LOG_ERROR()
            END SELECT

 2100 FORMAT(' Error 2100: Unknown VTK variable flag ',A,':',A,       /&
         'Available flags are:',                                      /&
         '  1 : Void fraction (EP_g)',                                /&
         '  2 : Gas pressure, solids pressure (P_g, P_star)',         /&
         '  3 : Gas velocity (U_g, V_g, W_g)',                        /&
         '  4 : Solids velocity (U_s, V_s, W_s)',                     /&
         '  5 : Solids density (ROP_s)',                              /&
         '  6 : Gas and solids temperature (T_g, T_s1, T_s2)',        /&
         '  7 : Gas and solids mass fractions (X_g, X-s)',            /&
         '  8 : Granular temperature (G)',                            /&
         '  9 : User defined scalars',                                /&
         ' 10 : Reaction Rates',                                      /&
         ' 11 : Turbulence quantities (k and ε)',                     /&
         ' 12 : Gas Vorticity magn and Lambda_2(VORTICITY,LAMBDA_2)', /&
         '100: Processor assigned to scalar cell (Partition)',        /&
         '101: Boundary condition flag for scalar cell (BC_ID)')

         ENDDO

! Activate particle orientation calculation if one vtk region needs it.
         IF(VTK_PART_ORIENTATION(L)) PARTICLE_ORIENTATION = .TRUE.


! If there is only one DES solids phase, use species aliases when
! saving DES_X_S in the vtp file
! Set VTK_PHASE_FOR_DES_X_S(VTK_REGION) equal to the DES phase ID
! if there is only one des phase in the vtk region, set it to zero
! otherwise
! VTK_PART_GROUP_STRING keeps a list of all phases present in the vtp file. Used
! as meta data in the pvd file.

         VTK_PART_GROUP_STRING(L) = ''
         M_IN_VTK=0
         VTK_PHASE_FOR_DES_X_S(L) = 0
         DO M = MMAX+1,MMAX+DES_MMAX
            IF(VTK_PART_PHASE(L,M)) THEN
               M_IN_VTK = M_IN_VTK + 1
               VTK_PHASE_FOR_DES_X_S(L) = M
               WRITE(SUBM,*)M
               VTK_PART_GROUP_STRING(L) = TRIM(VTK_PART_GROUP_STRING(L)) // ', ' //ADJUSTL(SUBM)
            ENDIF
         ENDDO

         IF(VTK_PART_GROUP_STRING(L)(1:1)==',') VTK_PART_GROUP_STRING(L)(1:1)=''

         IF(M_IN_VTK/=1) VTK_PHASE_FOR_DES_X_S(L) = 0

! Force chain

         IF(VTK_DATA(L)=='F') WRITE_FORCE_CHAIN = .TRUE.
         IF(VTK_PART_COORDINATION(L)) WRITE_FORCE_CHAIN = .TRUE.
         IF(VTK_FC_JKR_STATE(L)) WRITE_FORCE_CHAIN = .TRUE.
         IF(VTK_FC_JKR_RADIUS(L)) WRITE_FORCE_CHAIN = .TRUE.
         IF(VTK_FC_JKR_GAMMA(L)) WRITE_FORCE_CHAIN = .TRUE.
         IF(VTK_FC_TANG_DISP(L)) WRITE_FORCE_CHAIN = .TRUE.
         IF(VTK_FC_ROLL_DISP(L)) WRITE_FORCE_CHAIN = .TRUE.
         IF(VTK_FC_TWIST_ANGLE(L)) WRITE_FORCE_CHAIN = .TRUE.
         IF(VTK_PART_OVERLAP_STATS(L)) WRITE_FORCE_CHAIN = .TRUE.

         IF(IS_SMP .and. WRITE_FORCE_CHAIN) THEN
            IF (.NOT. SMP_WARNING_PRINTED) THEN
               WRITE(ERR_MSG, *) 'WARNING: Force chain data is not supported in SMP runs, and will be turned off.'
               CALL LOG_WARNING()
               SMP_WARNING_PRINTED = .TRUE.
            ENDIF
            WRITE_FORCE_CHAIN = .FALSE.
         ENDIF

! Set up flags to save reaction rates in vtu/vtp files
! Cell data (vtu files): Fluid (gas phase or gas and solids phase for TFM) rates
!                        and particle rates, accumulated at the cell level.
! Particle data (DEM, CGDEM, PIC): each particle store its own rates.
         IF(VTK_DATA(L)=='C') THEN
            DO R = 1,NO_OF_RXNS
               IF(VTK_FLUID_RRate(L,R)) SAVE_FLUID_RRATES = .TRUE.
            END DO
            DO R = 1,NO_OF_DES_RXNS
               IF(VTK_DES_RRate(L,R)) SAVE_DES_RRATES = .TRUE.
            END DO
         ELSEIF(VTK_DATA(L)=='P') THEN
            DO R = 1,NO_OF_DES_RXNS
               IF(VTK_PART_RRate(L,R)) SAVE_PART_RRATES = .TRUE.
            END DO
         ENDIF


      ENDDO   ! end loop over (l = 1,dimension_vtk)


      IF(DES_USR_VAR_SIZE>VTK_PART_USRmax) THEN
         WRITE(ERR_MSG,1005) DES_USR_VAR_SIZE,VTK_PART_USRmax
         CALL LOG_ERROR()
      ENDIF

 1005 FORMAT('Error 1005: Invalid input: DES_USR_VAR_SIZE>&
      &VTK_PART_USRmax:',/ 'DES_USR_VAR_SIZE = ',I6,/ &
      'VTK_PART_USRmax  = ',I6,/ &
      'Correct DES_USR_VAR_SIZE or increase ', / &
      'VTK_PART_USRmax in model/cartesian_grid/vtk_mod.f')

      RETURN

      END SUBROUTINE CHECK_OUTPUT_VTK


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_OUTPUT_MONITOR                                    !
!  Purpose: Check the output control namelist section                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_OUTPUT_MONITOR

! Global Variables:
!---------------------------------------------------------------------//
      use calc_cell_mod, only: calc_cell, calc_loc

      use monitor

      use fldvar

! Global Parameters:
!---------------------------------------------------------------------//
! Basic grid information
      use geometry, only: NO_I, XLENGTH, DX, IMAX, IMAX2
      use geometry, only: NO_J, YLENGTH, DY, JMAX, JMAX2
      use geometry, only: NO_K, ZLENGTH, DZ, KMAX, KMAX2
      use geometry, only: X_MIN, Y_MIN, Z_MIN
! Number aliases
      use param1, only: UNDEFINED, ZERO
      USE physprop, only: mmax, nmax, mw_mix_g
      USE physprop, only: mmax, nmax
      use rxns,     only: nrr, reactionrates, no_of_rxns, nrrmax
      use scalars, only : nscalar
      use discretelement, only: des_radius, pmass, pvol, ro_sol
      use discretelement, only: des_vel_new, des_pos_new, omega_new
      use discretelement, only: des_usr_var_size, des_usr_var
      use discretelement, only: residence_time
      use discretelement, only: max_pip
      use des_thermo,     only: des_t_s
      use des_rxns,       only: des_x_s, no_of_des_rxns, part_rrates_out


      use run, only: reinitializing, time
      use compar, only: mype,pe_io

      use monitor_functions, only: monitor_open_file
      use monitor_functions, only: monitor_write_header
      use vtk
      use param, only: dimension_3

! Function to compare two values
      use toleranc, only: COMPARE
      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      integer :: lc1, lc2, lc3, EEcnt
      integer :: wcnt, funit, funit_tavg

! Calculated indices of the wall boundary
      INTEGER :: I_w , I_e , J_s , J_n , K_b , K_t
! Integer error flag
      INTEGER :: IER
!......................................................................!

      do lc1 = 1, dimension_monitor

         monitor_defined(lc1) = .false.
         if (monitor_dt(lc1)  /=  undefined) monitor_defined(lc1) = .true.
         if (monitor_x_w(lc1) /= -undefined) monitor_defined(lc1) = .true.
         if (monitor_x_e(lc1) /=  undefined) monitor_defined(lc1) = .true.
         if (monitor_y_s(lc1) /= -undefined) monitor_defined(lc1) = .true.
         if (monitor_y_n(lc1) /=  undefined) monitor_defined(lc1) = .true.
         if (monitor_z_b(lc1) /= -undefined) monitor_defined(lc1) = .true.
         if (monitor_z_t(lc1) /=  undefined) monitor_defined(lc1) = .true.

         if(monitor_defined(lc1)) then
            if(monitor_dt(lc1)<zero) then
               write(err_msg,2001) trim(ival(monitor_dt(lc1)))
               call log_error()
            endif
         else if(reinitializing) then
            monitor_next_time(lc1) = undefined
         endif

 2001 FORMAT('Error 2001: Invalid value for MONITOR_DT = ',A, '. &
         &Acceptable values',/'are positive numbers (e.g., 0.1).')

         if(monitor_defined(lc1)) then

! Convert monitor coordinates to indices


            IF(MONITOR_X_W(lc1)/=-UNDEFINED .AND. MONITOR_X_E(lc1)/=UNDEFINED) THEN

! setting indices to 1 if there is no variation in the i (x) direction
               IF (NO_I) THEN
                  I_W = 1
                  I_E = 1
               ELSE
                  CALL CALC_CELL (X_MIN, MONITOR_X_W(lc1), DX, IMAX, I_W)
                  I_W = I_W + 1
                  CALL CALC_CELL (X_MIN, MONITOR_X_E(lc1), DX, IMAX, I_E)
! Monitor along zy plane, checking if far west or far east of domain
                  IF(MONITOR_X_W(lc1) == MONITOR_X_E(lc1)) THEN
                     IF(COMPARE(MONITOR_X_W(lc1),X_MIN)) THEN
                        I_W = 1
                        I_E = 1
                     ELSEIF(COMPARE(MONITOR_X_W(lc1),X_MIN+XLENGTH)) THEN
                        I_W = IMAX2
                        I_E = IMAX2
                     ELSE
                        I_W = I_E
                     ENDIF
                  ENDIF
               ENDIF
               MONITOR_I_W(lc1) = I_W
               MONITOR_I_E(lc1) = I_E
            ENDIF

            IF(MONITOR_Y_S(lc1)/=-UNDEFINED .AND. MONITOR_Y_N(lc1)/=UNDEFINED) THEN
! setting indices to 1 if there is no variation in the j (y) direction
               IF(NO_J) THEN
                  J_S = 1
                  J_N = 1
               ELSE
                  CALL CALC_CELL (Y_MIN, MONITOR_Y_S(lc1), DY, JMAX, J_S)
                  J_S = J_S + 1
                  CALL CALC_CELL (Y_MIN, MONITOR_Y_N(lc1), DY, JMAX, J_N)
! Monitor along xz plane, checking if far south or far north of domain
                  IF(MONITOR_Y_S(lc1) == MONITOR_Y_N(lc1)) THEN
                     IF(COMPARE(MONITOR_Y_S(lc1),Y_MIN)) THEN
                        J_S = 1
                        J_N = 1
                     ELSE IF (COMPARE(MONITOR_Y_S(lc1),Y_MIN+YLENGTH)) THEN
                        J_S = JMAX2
                        J_N = JMAX2
                     ELSE
                        J_S = J_N
                     ENDIF
                  ENDIF
               ENDIF
               MONITOR_J_S(lc1) = J_S
               MONITOR_J_N(lc1) = J_N
            ENDIF

            IF(MONITOR_Z_B(lc1)/=-UNDEFINED .AND. MONITOR_Z_T(lc1)/=UNDEFINED.OR.NO_K) THEN
! setting indices to 1 if there is no variation in the k (z) direction
               IF(NO_K)THEN
                  K_B = 1
                  K_T = 1
               ELSE
                  CALL CALC_CELL (Z_MIN, MONITOR_Z_B(lc1), DZ, KMAX, K_B)
                  K_B = K_B + 1
                  CALL CALC_CELL (Z_MIN, MONITOR_Z_T(lc1), DZ, KMAX, K_T)
! Monitor along xy plane, checking if far bottom or far top of domain
                  IF(MONITOR_Z_B(lc1) == MONITOR_Z_T(lc1)) THEN
                     IF(COMPARE(MONITOR_Z_B(lc1),Z_MIN)) THEN
                        K_B = 1
                        K_T = 1
                     ELSEIF(COMPARE(MONITOR_Z_B(lc1),Z_MIN+ZLENGTH)) THEN
                        K_B = KMAX2
                        K_T = KMAX2
                     ELSE
                        K_B = K_T
                     ENDIF
                  ENDIF
               ENDIF
               MONITOR_K_B(lc1) = K_B
               MONITOR_K_T(lc1) = K_T
            ENDIF


            IF(MONITOR_TYPE(LC1)<200) THEN
! CHECK FOR VALID VALUES
            IER = 0
            IF (MONITOR_K_B(lc1)<1 .OR. MONITOR_K_B(lc1)>KMAX2) IER = 1
            IF (MONITOR_J_S(lc1)<1 .OR. MONITOR_J_S(lc1)>JMAX2) IER = 1
            IF (MONITOR_I_W(lc1)<1 .OR. MONITOR_I_W(lc1)>IMAX2) IER = 1
            IF (MONITOR_K_T(lc1)<1 .OR. MONITOR_K_T(lc1)>KMAX2) IER = 1
            IF (MONITOR_J_N(lc1)<1 .OR. MONITOR_J_N(lc1)>JMAX2) IER = 1
            IF (MONITOR_I_E(lc1)<1 .OR. MONITOR_I_E(lc1)>IMAX2) IER = 1
            IF (MONITOR_K_B(lc1) > MONITOR_K_T(lc1)) IER = 1
            IF (MONITOR_J_S(lc1) > MONITOR_J_N(lc1)) IER = 1
            IF (MONITOR_I_W(lc1) > MONITOR_I_E(lc1)) IER = 1

            IF(IER /= 0)THEN
               WRITE(ERR_MSG,1100) lc1,                                      &
                  'X', MONITOR_X_W(lc1), MONITOR_X_E(lc1),'I',MONITOR_I_W(lc1),MONITOR_I_E(lc1), &
                  'Y', MONITOR_Y_S(lc1), MONITOR_Y_N(lc1),'J',MONITOR_J_S(lc1),MONITOR_J_N(lc1), &
                  'Z', MONITOR_Z_B(lc1), MONITOR_Z_T(lc1),'K',MONITOR_K_B(lc1),MONITOR_K_T(lc1)
               CALL LOG_ERROR()
            ENDIF

            ENDIF

 1100 FORMAT('Error 1100: Invalid location specified for Monitor ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),&
         /3x,'This usually occurs when the Monitor region is outside the fluid region,', &
         /3x,'or when the Monitor region is smaller than the grid spacing.')

   ! Get the monitor geometry type (point, plane, volume)
      call get_monitor_geo_type(lc1)

! Get the type
            select case (monitor_type(lc1))
               case(0)
                  monitor_typename(lc1) = 'Point value'
               case(1)
                  monitor_typename(lc1) = 'Sum over region'
               case(2)
                  monitor_typename(lc1) = 'Minimum over region'
               case(3)
                  monitor_typename(lc1) = 'Maximum over region'
               case(4)
                  monitor_typename(lc1) = 'Arithmetic average over region'
               case(5)
                  monitor_typename(lc1) = 'Standard deviation over region'
               case(6)
                  monitor_typename(lc1) = 'Area-weighted average over surface'
               case(7)
                  monitor_typename(lc1) = 'Flow rate across surface'
               case(8)
                  monitor_typename(lc1) = 'Mass flow rate across surface'
               case(9)
                  monitor_typename(lc1) = 'Mass-weighted average over surface'
               case(10)
                  monitor_typename(lc1) = 'Volumetric flow rate over surface'
               case(11)
                  monitor_typename(lc1) = 'Volume integral'
               case(12)
                  monitor_typename(lc1) = 'Volume-weighted average'
               case(13)
                  monitor_typename(lc1) = 'Mass-weighted volume integral'
               case(14)
                  monitor_typename(lc1) = 'Mass-weighted volume average'
               case(101)
                  monitor_typename(lc1) = 'Sum over region'
               case(102)
                  monitor_typename(lc1) = 'Minimum over region'
               case(103)
                  monitor_typename(lc1) = 'Maximum over region'
               case(104)
                  monitor_typename(lc1) = 'Arithmetic average over region'
               case(105)
                  monitor_typename(lc1) = 'Standard deviation over region'
               case(106)
                  monitor_typename(lc1) = 'Mass-weighted average over region'
               case(107)
                  monitor_typename(lc1) = 'Volume-weighted average over region'
               case(108)
                  monitor_typename(lc1) = 'Flow rate across surface'
               case(109)
                  monitor_typename(lc1) = 'Mass-weighted flow rate across surface'
               case(110)
                  monitor_typename(lc1) = 'Volume-weighted flow rate across surface'
               case(201)
                  monitor_typename(lc1) = 'Sum of forces acting on STL geometry'
            end select


! Get the list of variables in the monitor
! Turn off monitor variable if variable is not allocated
            wcnt = 0
            if(.not.allocated(ep_g)) monitor_ep_g(lc1) = .False.
            if(monitor_ep_g(lc1)) call increment_monitor_var(lc1,'Void fraction (monitor_ep_g)',wcnt)

            if(.not.allocated(ro_g)) monitor_ro_g(lc1) = .False.
            if(monitor_ro_g(lc1)) call increment_monitor_var(lc1,'Gas density (monitor_ro_g)',wcnt)

            if(.not.allocated(p_g)) monitor_p_g(lc1) = .False.
            if(monitor_p_g (lc1)) call increment_monitor_var(lc1,'Gas pressure (monitor_p_g)',wcnt)

            if(.not.allocated(u_g)) monitor_u_g(lc1) = .False.
            if(monitor_u_g (lc1)) call increment_monitor_var(lc1,'Gas u-velocity (monitor_u_g)',wcnt)

            if(.not.allocated(v_g)) monitor_v_g(lc1) = .False.
            if(monitor_v_g (lc1)) call increment_monitor_var(lc1,'Gas v-velocity (monitor_v_g)',wcnt)

            if(.not.allocated(w_g)) monitor_w_g(lc1) = .False.
            if(monitor_w_g (lc1)) call increment_monitor_var(lc1,'Gas w-velocity (monitor_w_g)',wcnt)

            if(.not.allocated(t_g)) monitor_t_g(lc1) = .False.
            if(monitor_t_g (lc1)) call increment_monitor_var(lc1,'Gas temperature (monitor_t_g)',wcnt)

            if(.not.allocated(mw_mix_g)) monitor_mw_mix_g(lc1) = .False.
            if(monitor_mw_mix_g (lc1)) call increment_monitor_var(lc1,'Gas molecular weight (monitor_mw_mix)',wcnt)

            if(.not.allocated(x_g)) monitor_x_g(lc1,:) = .False.
            if(.not.allocated(x_g)) monitor_y_g(lc1,:) = .False.  ! molar fraction is computed from mass fraction, y_g is not a field variable
            do lc2 = 1, nmax(0)
               if(monitor_x_g(lc1,lc2)) call increment_monitor_var(lc1,'Gas species mass fraction (monitor_x_g)',wcnt,', species:', lc2)
               if(monitor_y_g(lc1,lc2)) call increment_monitor_var(lc1,'Gas species molar fraction (monitor_y_g)',wcnt,', species:', lc2)
            enddo

            if(.not.allocated(k_turb_g)) monitor_k_turb_g(lc1) = .False.
            if(monitor_k_turb_g(lc1)) call increment_monitor_var(lc1,'Turbulent kinetic energy (monitor_k_turb_g)',wcnt)

            if(.not.allocated(e_turb_g)) monitor_e_turb_g(lc1) = .False.
            if(monitor_e_turb_g(lc1)) call increment_monitor_var(lc1,'Turbulent dissipation rate (monitor_e_turb_g)',wcnt)

            if(.not.allocated(p_star)) monitor_p_star(lc1) = .False.
            if(monitor_p_star (lc1)) call increment_monitor_var(lc1,'Solids pressure preventing overpacking (monitor_p_star)',wcnt)

            do lc2 = 1, mmax
               if(.not.allocated(ro_s)) monitor_ep_s(lc1,lc2) = .False. ! ep_s is not a field variable, it is computed from solids density
               if(monitor_ep_s   (lc1,lc2)) call increment_monitor_var(lc1,'Solids volume fraction (monitor_ep_s)',wcnt,', phase:', lc2)

               if(.not.allocated(rop_s)) monitor_rop_s(lc1,lc2) = .False.
               if(monitor_rop_s  (lc1,lc2)) call increment_monitor_var(lc1,'Solids bulk density (monitor_rop_s)',wcnt,', phase:', lc2)

               if(.not.allocated(ro_s)) monitor_ro_s(lc1,lc2) = .False.
               if(monitor_ro_s   (lc1,lc2)) call increment_monitor_var(lc1,'Solids density (monitor_rop_s)',wcnt,', phase:', lc2)

               if(.not.allocated(u_s)) monitor_u_s(lc1,lc2) = .False.
               if(monitor_u_s    (lc1,lc2)) call increment_monitor_var(lc1,'Solids u-velocity (monitor_u_s)',wcnt,', phase:', lc2)

               if(.not.allocated(v_s)) monitor_v_s(lc1,lc2) = .False.
               if(monitor_v_s    (lc1,lc2)) call increment_monitor_var(lc1,'Solids v-velocity (monitor_v_s)',wcnt,', phase:', lc2)

               if(.not.allocated(w_s)) monitor_w_s(lc1,lc2) = .False.
               if(monitor_w_s    (lc1,lc2)) call increment_monitor_var(lc1,'Solids w-velocity (monitor_w_s)',wcnt,', phase:', lc2)

               if(.not.allocated(p_s)) monitor_p_s(lc1,lc2) = .False.
               if(monitor_p_s    (lc1,lc2)) call increment_monitor_var(lc1,'Solids pressure (monitor_p_s)',wcnt,', phase:', lc2)

               if(.not.allocated(t_s)) monitor_t_s(lc1,lc2) = .False.
               if(monitor_t_s    (lc1,lc2)) call increment_monitor_var(lc1,'Solids temperature (monitor_t_s)',wcnt,', phase:', lc2)

               if(.not.allocated(theta_m)) monitor_theta_m(lc1,lc2) = .False.
               if(monitor_theta_m(lc1,lc2)) call increment_monitor_var(lc1,'Solids granular temperature (monitor_theta_m)',wcnt,', phase:', lc2)


               if(.not.allocated(x_s)) monitor_x_s(lc1,lc2,:) = .False.
               do lc3 = 1, nmax(lc2)
                  if(monitor_x_s(lc1,lc2,lc3)) call increment_monitor_var(lc1,'Solids species mass fraction (monitor_x_s)',wcnt,', phase:', lc2, ', species:', lc3)
               enddo
            enddo

            if(.not.allocated(scalar)) monitor_scalar(lc1,:) = .False.
            do lc2 = 1, nscalar
               if(monitor_scalar(lc1,lc2)) call increment_monitor_var(lc1,'Scalar (monitor_scalar)',wcnt,', scalar:', lc2)
            enddo

            if(.not.allocated(reactionrates)) monitor_rrate(lc1,:) = .False.
            do lc2 = 1, min(nrrmax,nrr)
               if(monitor_rrate(lc1,lc2)) call increment_monitor_var(lc1,'Rrate array, cell data (monitor_rrate)',wcnt,', index:', lc2)
            enddo

            if(.not.allocated(fluid_rrates_out)) monitor_fluid_rrate(lc1,:) = .False.
            do lc2 = 1, NO_OF_RXNS
               if(monitor_fluid_rrate(lc1,lc2)) then
                  call increment_monitor_var(lc1,'Fluid reaction rate (monitor_fluid_rrate)',wcnt,', reaction:', lc2)
               endif
            enddo

            if(.not.allocated(Des_rrates_out)) monitor_des_rrate(lc1,:) = .False.
            do lc2 = 1, NO_OF_DES_RXNS
               if(monitor_des_rrate(lc1,lc2)) then
                  call increment_monitor_var(lc1,'DES reaction rate (monitor_des_rrate)',wcnt,', reaction:', lc2)
               endif
            enddo

            ! Save the number of Eulerian monitors
            EEcnt = wcnt

            if(.not.allocated(des_radius)) monitor_radius(lc1) = .False.
            if(monitor_radius(lc1)) call increment_monitor_var(lc1,'Particle radius (monitor_radius)',wcnt)

            if(.not.allocated(pmass)) monitor_pmass(lc1) = .False.
            if(monitor_pmass(lc1)) call increment_monitor_var(lc1,'Particle mass (monitor_pmass)',wcnt)

            if(.not.allocated(pmass)) monitor_pcount(lc1) = .False. ! pcount is computed directly from pmass
            if(monitor_pcount(lc1)) call increment_monitor_var(lc1,'Particle count (monitor_pcount)',wcnt)

            if(.not.allocated(pvol)) monitor_pvol(lc1) = .False.
            if(monitor_pvol(lc1)) call increment_monitor_var(lc1,'Particle volume (monitor_pvol)',wcnt)

            if(.not.allocated(ro_sol)) monitor_ro_p(lc1) = .False.
            if(monitor_ro_p(lc1)) call increment_monitor_var(lc1,'Particle density (monitor_ro_p)',wcnt)

            if(.not.allocated(des_vel_new)) monitor_vel_x(lc1) = .False.
            if(monitor_vel_x(lc1)) call increment_monitor_var(lc1,'Particle x-velocity (monitor_vel_x)',wcnt)

            if(.not.allocated(des_vel_new)) monitor_vel_y(lc1) = .False.
            if(monitor_vel_y(lc1)) call increment_monitor_var(lc1,'Particle y-velocity (monitor_vel_y)',wcnt)

            if(.not.allocated(des_vel_new)) monitor_vel_z(lc1) = .False.
            if(monitor_vel_z(lc1)) call increment_monitor_var(lc1,'Particle z-velocity (monitor_vel_z)',wcnt)

            if(.not.allocated(des_pos_new)) monitor_pos_x(lc1) = .False.
            if(monitor_pos_x(lc1)) call increment_monitor_var(lc1,'Particle x-position (monitor_pos_x)',wcnt)

            if(.not.allocated(des_pos_new)) monitor_pos_y(lc1) = .False.
            if(monitor_pos_y(lc1)) call increment_monitor_var(lc1,'Particle y-position (monitor_pos_y)',wcnt)

            if(.not.allocated(des_pos_new)) monitor_pos_z(lc1) = .False.
            if(monitor_pos_z(lc1)) call increment_monitor_var(lc1,'Particle z-position (monitor_pos_z)',wcnt)

            if(.not.allocated(omega_new)) monitor_rot_x(lc1) = .False.
            if(monitor_rot_x(lc1)) call increment_monitor_var(lc1,'Particle angular x-velocity (monitor_rot_x)',wcnt)

            if(.not.allocated(omega_new)) monitor_rot_y(lc1) = .False.
            if(monitor_rot_y(lc1)) call increment_monitor_var(lc1,'Particle angular y-velocity (monitor_rot_y)',wcnt)

            if(.not.allocated(omega_new)) monitor_rot_z(lc1) = .False.
            if(monitor_rot_z(lc1)) call increment_monitor_var(lc1,'Particle angular z-velocity (monitor_rot_z)',wcnt)

            if(.not.allocated(des_t_s)) monitor_t_p(lc1) = .False.
            if(monitor_t_p(lc1)) call increment_monitor_var(lc1,'Particle temperature (monitor_t_p)',wcnt)

            if(.not.allocated(des_x_s)) monitor_x_p(lc1,:) = .False.
            do lc2 = 1, dimension_n_s
               if(monitor_x_p(lc1, lc2)) call increment_monitor_var(lc1,'Particle species mass fraction (monitor_x_p)',wcnt,', species:', lc2)
            end do

            if(.not.allocated(des_usr_var)) monitor_des_usr_var(lc1,:) = .False.
            do lc2 = 1, des_usr_var_size
               if(monitor_des_usr_var(lc1, lc2)) call increment_monitor_var(lc1,'Particle user-defined variable (monitor_des_usr_var)',wcnt,', index:', lc2)
            end do

            do lc2 = 1, NO_OF_DES_RXNS
               if(monitor_part_rrate(lc1,lc2)) then
                  call increment_monitor_var(lc1,'Particle reaction rate (monitor_part_rrate)',wcnt,', reaction:', lc2)
               endif
            enddo

            if(.not.allocated(residence_time)) monitor_part_residence_time(lc1) = .False.
            if(monitor_part_residence_time(lc1)) call increment_monitor_var(lc1,'Particle residence time (monitor_residence_time)',wcnt)

            if(monitor_coll_count (lc1)) call increment_monitor_var(lc1,'STL collision count (monitor_coll_count)',wcnt)

            if(monitor_wear (lc1)) call increment_monitor_var(lc1,'STL wear (monitor_collision_wear)',wcnt)

            if(monitor_normal_force (lc1)) then
               call increment_monitor_var(lc1,'STL normal force x (monitor_normal_force)',wcnt)
               call increment_monitor_var(lc1,'STL normal force y (monitor_normal_force)',wcnt)
               call increment_monitor_var(lc1,'STL normal force z (monitor_normal_force)',wcnt)
            endif

            if(monitor_tang_force (lc1)) then
               call increment_monitor_var(lc1,'STL tang force x (monitor_tang_force)',wcnt)
               call increment_monitor_var(lc1,'STL tang force y (monitor_tang_force)',wcnt)
               call increment_monitor_var(lc1,'STL tang force z (monitor_tang_force)',wcnt)
            endif

            if(monitor_force (lc1)) then
               call increment_monitor_var(lc1,'STL force x (monitor_force)',wcnt)
               call increment_monitor_var(lc1,'STL force y (monitor_force)',wcnt)
               call increment_monitor_var(lc1,'STL force z (monitor_force)',wcnt)
            endif

            if(monitor_torque (lc1)) then
               call increment_monitor_var(lc1,'STL torque x (monitor_torque)',wcnt)
               call increment_monitor_var(lc1,'STL torque y (monitor_torque)',wcnt)
               call increment_monitor_var(lc1,'STL torque z (monitor_torque)',wcnt)
            endif

            if(monitor_q(lc1)) call increment_monitor_var(lc1,'STL heat transfer (monitor_q)',wcnt)


            if(wcnt == 0) then
               monitor_defined(lc1) = .false.
               write(err_msg, 4455) lc1
               call log_warning()

            else if( EEcnt /= 0 .and. EEcnt /= wcnt) then
               write(err_msg, 4460) lc1
               call log_error()

            else
               monitor_var_count(lc1) = wcnt
            endif

            if(reinitializing .and. monitor_next_time(lc1) == undefined) &
               monitor_next_time(lc1) = time

! Time-averaged data
            if(monitor_tavg_dt(lc1) <= 0.0D0) then
               write(err_msg, 4470) lc1
               call log_error()
            endif

            if(monitor_tavg_reset_dt(lc1) <= 0.0D0) then
               write(err_msg, 4480) lc1
               call log_error()
            endif


            monitor_tavg(lc1) = (monitor_tavg_dt(lc1) /= undefined)

            if(monitor_tavg(lc1)) then
               monitor_tavg_vars(lc1,:) = zero
               monitor_tavg_count(lc1) = zero
            endif

! Header is now written in write_monitor. Change was needed
! because we now write the monitors's area or volume in the header,
! and this can only be computed after mesh generation.
            ! if(mype == pe_io) then
            !    call monitor_open_file(lc1, funit,.false.)
            !    call monitor_write_header(lc1, funit)
            !    close(funit)
            !    if(monitor_tavg(lc1)) then
            !       call monitor_open_file(lc1, funit_tavg,.true.)
            !       call monitor_write_header(lc1, funit_tavg)
            !       close(funit_tavg)
            !    endif
            ! endif

         endif ! monitor is defined
      enddo ! loop over monitors (lc1)


4455  format('Warning 4455: Monitor region ',i3,&
           ' has no variables selected.')

4460  format('Monitors cannot collect both Eulerian and Lagrangian data.',&
           &/'Correct the variable list for monitor:',I3)

4470  format('Sampling interval (MONITOR_TAVG_DT) is negative or zero, must be strictly positive.',&
           &/'Correct the sampling interval for monitor:',I3)

4480  format('Reset interval (MONITOR_TAVG_RESET_DT) is negative or zero, must be strictly positive.',&
           &/'Correct the reset interval for monitor:',I3)

      RETURN

      END SUBROUTINE CHECK_OUTPUT_MONITOR

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CHECK_MONITOR_RRATES                                    !
!  Purpose: Check if reaction rates need to be written in monitors     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_MONITOR_RRATES

! Global Variables:
!---------------------------------------------------------------------//
      use monitor
      use rxns,     only: no_of_rxns
      use rxns,     only: SAVE_FLUID_RRATES, SAVE_DES_RRATES
      use rxns,     only: FLUID_RRATES_OUT, DES_RRATES_OUT
      use des_rxns, only: no_of_des_rxns
      use des_rxns, only: SAVE_PART_RRATES
      use des_rxns, only: part_rrates_out

      use vtk
      use param, only: dimension_3

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      integer :: lc1, lc2

!......................................................................!

      do lc1 = 1, dimension_monitor

! Set flags to save fluid reaction rates (cell data)
         do lc2 = 1, NO_OF_RXNS
            if(monitor_fluid_rrate(lc1,lc2)) SAVE_FLUID_RRATES = .TRUE.
         enddo

! Set flags to save des reaction rates (will be converted to cell data)
         do lc2 = 1, NO_OF_DES_RXNS
            if(monitor_des_rrate(lc1,lc2)) SAVE_DES_RRATES = .TRUE.
         enddo

! Set flags to save des reaction rates (particle data)
         do lc2 = 1, NO_OF_DES_RXNS
            if(monitor_part_rrate(lc1,lc2)) SAVE_PART_RRATES = .TRUE.
         enddo

! Allocate reaction rates arrays used to write reaction rates in vtk files or monitors
      if (SAVE_FLUID_RRATES) then
         if(allocated(Fluid_RRates_out)) deallocate(Fluid_RRates_out)
         allocate(Fluid_RRates_out(DIMENSION_3,NO_OF_RXNS))
         Fluid_RRates_out(:,:) = ZERO
      endif

      if (SAVE_DES_RRATES) then
         if(allocated(Des_RRates_out)) deallocate(Des_RRates_out)
         allocate(Des_RRates_out(DIMENSION_3,NO_OF_DES_RXNS))
         Des_RRates_out(:,:) = ZERO
      endif

      enddo ! loop over monitors (lc1)


      RETURN

      END SUBROUTINE CHECK_MONITOR_RRATES

      subroutine increment_monitor_var(lid, label, wcnt, i2_label, i2, i3_label, i3)

         use monitor

         implicit none
         integer, intent(in) :: lid
         character(len=*), intent(in) :: label
         character(len=*), intent(in), optional :: i2_label
         integer, intent(in), optional :: i2
         character(len=*), intent(in), optional :: i3_label
         integer, intent(in), optional :: i3
         integer, intent(inout) :: wcnt


         character (len=6) :: sub2,sub3
         character (len=256) :: label2, label3
         character (len=284) :: var_name

         wcnt = wcnt + 1

         monitor_varname_list(lid,wcnt) = '- '//trim(label)

         if(present(i2_label).and.present(i2)) then
            write(label2,'(A,I3)') trim(i2_label),i2
            monitor_varname_list(lid,wcnt) = trim(monitor_varname_list(lid,wcnt))//adjustl(label2)
         endif

         if(present(i3_label).and.present(i3)) then
            write(label3,'(A,I3)') trim(i3_label),i3
            monitor_varname_list(lid,wcnt) = trim(monitor_varname_list(lid,wcnt))//adjustl(label3)
         endif

         return

      end subroutine increment_monitor_var

      subroutine get_monitor_geo_type(lid)

         use monitor
         use geometry, only: no_k

         implicit none
         integer, intent(in) :: lid
         double precision :: xw,xe,ys,yn,zb,zt
         logical :: same_x, same_y, same_z

         xw = monitor_x_w(lid)
         xe = monitor_x_e(lid)
         ys = monitor_y_s(lid)
         yn = monitor_y_n(lid)
         zb = monitor_z_b(lid)
         zt = monitor_z_t(lid)

         same_x = (xw==xe)
         same_y = (ys==yn)
         same_z = (zb==zt)

         if(no_k) then

            if(same_x.and.same_y) then
               monitor_geotype(lid) = 'PT'
            elseif(same_x.and.(.not.same_y)) then
               monitor_geotype(lid) = 'YZ'
            elseif(same_y.and.(.not.same_x)) then
               monitor_geotype(lid) = 'XZ'
            elseif((.not.same_x).and.(.not.same_y)) then
               monitor_geotype(lid) = 'XY'
            else
               WRITE(ERR_MSG,1000) lid,xw,xe,ys,yn,zb,zt
               CALL LOG_ERROR()
            endif

         else
            if(same_x.and.same_y.and.same_z) then
               monitor_geotype(lid) = 'PT'
            elseif(same_x.and.(.not.same_y).and.(.not.same_z)) then
               monitor_geotype(lid) = 'YZ'
            elseif(same_y.and.(.not.same_x).and.(.not.same_z)) then
               monitor_geotype(lid) = 'XZ'
            elseif(same_z.and.(.not.same_x).and.(.not.same_y)) then
               monitor_geotype(lid) = 'XY'
            elseif((.not.same_x).and.(.not.same_y).and.(.not.same_z)) then
               monitor_geotype(lid) = '3D'
            else
               WRITE(ERR_MSG,1000) lid,xw,xe,ys,yn,zb,zt
               CALL LOG_ERROR()
            endif
         endif


1000 FORMAT('Error 1000: Unable to determine monitor geometry type (point, plane, volume) for monitor: ',I4, &
           /'West coordinate   (monitor_x_w)= ', G12.5,&
           /'East coordinate   (monitor_x_e)= ', G12.5,&
           /'South coordinate  (monitor_y_s)= ', G12.5,&
           /'North coordinate  (monitor_y_n)= ', G12.5,&
           /'Bottom coordinate (monitor_z_b)= ', G12.5,&
           /'Top coordinate    (monitor_z_t)= ', G12.5 &
           )

         return

      end subroutine get_monitor_geo_type

END MODULE CHECK_OUTPUT_CONTROL_MOD
