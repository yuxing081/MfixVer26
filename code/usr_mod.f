MODULE usr
! MFIX user module: user-defined namelist variables and other user data
! may be declared here and accessed with USE usr. Allocate user arrays
! in usr0 when using the CALL_USR hooks.
!
! These stateless JKR kernels are called directly by the contact routines,
! independently of CALL_USR. Contact histories remain in discretelement.
   IMPLICIT NONE
   PUBLIC :: JKR_TWIST_UPDATE, JKR_ROLL_UPDATE, JKR_TRANSPORT_HISTORY
CONTAINS
! Transport both tangent histories with one contact-frame rotation. FRAME
! stores the owner-directed normal and the time of the last force evaluation.
! CFNEWVALUES rotates spheres with its updated OMEGA_NEW; at the next force
! evaluation that angular velocity spans TIME_NOW-FRAME(4), not DTSOLID.
   SUBROUTINE JKR_TRANSPORT_HISTORY(NORMAL, MEAN_OMEGA, TIME_NOW, &
      FRAME, PFT, PFR, VALID)
      USE GSP_MATH_MOD, ONLY: ROTMAT_FROM_NORMALS
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      DOUBLE PRECISION, INTENT(IN) :: NORMAL(3), MEAN_OMEGA(3), TIME_NOW
      DOUBLE PRECISION, INTENT(INOUT) :: FRAME(4), PFT(3), PFR(3)
      LOGICAL, INTENT(OUT) :: VALID
      DOUBLE PRECISION :: N(3), OLD_N(3), R(3,3), SPIN(3,3), Q(3,3)
      DOUBLE PRECISION :: NORM, OLD_NORM, ANGLE, C, S, DT, H(3,2)
      DOUBLE PRECISION :: RESIDUAL
      INTEGER :: I, J

      VALID = .FALSE.
      IF(.NOT.ALL(IEEE_IS_FINITE(NORMAL))) RETURN
      IF(.NOT.ALL(IEEE_IS_FINITE(MEAN_OMEGA))) RETURN
      IF(.NOT.ALL(IEEE_IS_FINITE(FRAME))) RETURN
      IF(.NOT.ALL(IEEE_IS_FINITE(PFT))) RETURN
      IF(.NOT.ALL(IEEE_IS_FINITE(PFR))) RETURN
      IF(.NOT.IEEE_IS_FINITE(TIME_NOW)) RETURN
      NORM = SQRT(DOT_PRODUCT(NORMAL,NORMAL))
      IF(NORM <= TINY(1.D0)) RETURN
      N = NORMAL/NORM
      OLD_NORM = SQRT(DOT_PRODUCT(FRAME(1:3),FRAME(1:3)))
      H(:,1) = PFT
      H(:,2) = PFR
      IF(OLD_NORM > TINY(1.D0)) THEN
         DT = TIME_NOW-FRAME(4)
         IF(DT < 0.D0) RETURN
         OLD_N = FRAME(1:3)/OLD_NORM
! The native antiparallel branch writes the diagonal before adding uu^T;
! initialize off-diagonal entries as well. Never pass a zero normal.
         R = 0.D0
         CALL ROTMAT_FROM_NORMALS(OLD_N,N,R,EPS=1.D-14)
         ANGLE = DT*DOT_PRODUCT(MEAN_OMEGA,N)
         C = COS(ANGLE)
         S = SIN(ANGLE)
         DO I=1,3
            DO J=1,3
               SPIN(I,J) = (1.D0-C)*N(I)*N(J)
            ENDDO
            SPIN(I,I) = SPIN(I,I)+C
         ENDDO
         SPIN(1,2) = SPIN(1,2)-S*N(3)
         SPIN(1,3) = SPIN(1,3)+S*N(2)
         SPIN(2,1) = SPIN(2,1)+S*N(3)
         SPIN(2,3) = SPIN(2,3)-S*N(1)
         SPIN(3,1) = SPIN(3,1)-S*N(2)
         SPIN(3,2) = SPIN(3,2)+S*N(1)
         Q = MATMUL(SPIN,R)
         H = MATMUL(Q,H)
      ENDIF
! Remove roundoff normal components after transport. A relative residual
! divided by a nearly cancelled history is not a valid stopping criterion.
      DO J=1,2
         RESIDUAL = DOT_PRODUCT(N,H(:,J))
         H(:,J) = H(:,J)-RESIDUAL*N
      ENDDO
      PFT = H(:,1)
      PFR = H(:,2)
      FRAME(1:3) = N
      FRAME(4) = TIME_NOW
      VALID = .TRUE.
   END SUBROUTINE JKR_TRANSPORT_HISTORY

   PURE SUBROUTINE JKR_TWIST_UPDATE(DT, RATE, A, KT, ETAT, FSLIDE, &
      ANGLE, MOMENT)
      DOUBLE PRECISION, INTENT(IN) :: DT, RATE, A, KT, ETAT, FSLIDE
      DOUBLE PRECISION, INTENT(INOUT) :: ANGLE
      DOUBLE PRECISION, INTENT(OUT) :: MOMENT
      DOUBLE PRECISION :: KQ, CQ, LIMIT, TRIAL
      KQ = KT*A*A/2.D0
      CQ = ETAT*A*A/2.D0
      LIMIT = 3.D0*ACOS(-1.D0)*A*FSLIDE/16.D0
      IF(KQ <= 0.D0 .OR. LIMIT <= 0.D0) THEN
         ANGLE = 0.D0
         MOMENT = 0.D0
         RETURN
      ENDIF
      ANGLE = ANGLE + DT*RATE
      TRIAL = -KQ*ANGLE - CQ*RATE
      MOMENT = MAX(-LIMIT,MIN(LIMIT,TRIAL))
! Return mapping keeps the *total* spring/dashpot torque on the yield surface.
      IF(ABS(TRIAL) > LIMIT) ANGLE = -(MOMENT+CQ*RATE)/KQ
   END SUBROUTINE JKR_TWIST_UPDATE

   PURE SUBROUTINE JKR_ROLL_UPDATE(DT, RATE, STIFFNESS, XI_CRIT, &
      DAMPING, XI, RESISTANCE)
      DOUBLE PRECISION, INTENT(IN) :: DT, RATE(3), STIFFNESS, XI_CRIT
      DOUBLE PRECISION, INTENT(IN) :: DAMPING
      DOUBLE PRECISION, INTENT(INOUT) :: XI(3)
      DOUBLE PRECISION, INTENT(OUT) :: RESISTANCE(3)
      DOUBLE PRECISION :: MAGNITUDE, LIMIT
      XI = XI + DT*RATE
      IF(DAMPING == 0.D0) THEN
! Preserve the legacy elastic spring-slider update, including operation order.
         MAGNITUDE = SQRT(DOT_PRODUCT(XI,XI))
         IF(MAGNITUDE > XI_CRIT .AND. MAGNITUDE > 0.D0) &
            XI = XI_CRIT*XI/MAGNITUDE
         RESISTANCE = STIFFNESS*XI
         RETURN
      ENDIF
      IF(STIFFNESS <= 0.D0) THEN
         XI = 0.D0
         RESISTANCE = 0.D0
         RETURN
      ENDIF
      RESISTANCE = STIFFNESS*XI + DAMPING*RATE
      LIMIT = STIFFNESS*XI_CRIT
      MAGNITUDE = SQRT(DOT_PRODUCT(RESISTANCE,RESISTANCE))
      IF(MAGNITUDE > LIMIT .AND. MAGNITUDE > 0.D0) THEN
         RESISTANCE = (LIMIT/MAGNITUDE)*RESISTANCE
         XI = (RESISTANCE-DAMPING*RATE)/STIFFNESS
      ENDIF
   END SUBROUTINE JKR_ROLL_UPDATE
END MODULE usr
