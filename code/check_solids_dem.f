#include "error.inc"

MODULE CHECK_SOLIDS_DEM_MOD

   USE constant, only: pi

   use des_thermo, only: SB_CONST
   use des_thermo, only: DES_Em
   USE des_thermo, only: CALC_COND_DES
   USE des_thermo, only: CALC_RADT_DES
   USE des_thermo, only: DES_MIN_COND_DIST
   USE des_thermo_cond, only: DO_AREA_CORRECTION
   USE discretelement
   USE param1, only: zero, one, half, undefined

   USE physprop, only: mmax, d_p0, ro_s0
   use physprop, only: K_S0
   use run, only: ks_model, ks_model_enum, ks_none, ks_musser
   USE run, only: time, units, solids_model, energy_eq
   USE SQ_PROPERTIES_MOD

   USE error_manager


CONTAINS

! Marshall (2009), Eq. (53): restitution-to-damping fit used by Chen JKR.
   DOUBLE PRECISION FUNCTION JKR_DAMPING_ALPHA(RESTITUTION)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: RESTITUTION

      IF(RESTITUTION >= ONE) THEN
         JKR_DAMPING_ALPHA = ZERO
      ELSE
         JKR_DAMPING_ALPHA = 1.2728D0 - 4.2783D0*RESTITUTION + &
            11.087D0*RESTITUTION**2 - 22.348D0*RESTITUTION**3 + &
            27.467D0*RESTITUTION**4 - 18.022D0*RESTITUTION**5 + &
            4.8218D0*RESTITUTION**6
      ENDIF
   END FUNCTION JKR_DAMPING_ALPHA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_DEM                                        !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CHECK_SOLIDS_DEM

      implicit none

! Skip the check when using DLB and the partition is adjusted
      IF(ADJUST_PARTITION) return

! Get reference radius and mass for each phase
      CALL GET_REF_PARTICLE_SIZE_AND_MASS
! Particle-particle collision parameters.
      CALL CHECK_SOLIDS_DEM_COLLISION
! DES cohesion model parameters.
      CALL CHECK_SOLIDS_DEM_COHESION
! Particle-particle conduction model parameters.
      IF (ENERGY_EQ) CALL CHECK_SOLIDS_DEM_ENERGY

      RETURN

      END SUBROUTINE CHECK_SOLIDS_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE CHECK_SOLIDS_DEM_ENERGY                                 !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_ENERGY


      IMPLICIT NONE
!......................................................................!
      INTEGER :: M, L  ! Loop indices for DEM solids
! Calculate phase index offset for certain inputs until it can be
! addressed in other ways. should not matter unlesss hybrid
      CHARACTER(len=64) :: MSG
      INTEGER :: MMAX_TOT
      LOGICAL :: ANY_CONDUCTION = .FALSE.

      DOUBLE PRECISION :: MASS_L, MASS_M, MASS_EFF
      DOUBLE PRECISION :: R_EFF, E_EFF
      INTEGER :: zero_ks0, undef_ks0, no_dem


! Radiation Equation:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Set flags for radiation in dem specific routine
      DO M = MMAX+1, MMAX+DES_MMAX
! GSP stands for glued sphere model
         IF(.NOT.(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP' .OR.SOLIDS_MODEL(M) == 'SQP' .OR.SOLIDS_MODEL(M) == 'GSP')) CYCLE
         IF(DES_Em(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_Em',M))
            CALL LOG_ERROR()
         ELSEIF(DES_Em(M) > ONE .OR. DES_Em(M) < ZERO) THEN
            WRITE(ERR_MSG, 1001) trim(iVar('DES_Em',M)), &
               ival(DES_EM(M))
            CALL LOG_ERROR()
         ELSEIF(DES_Em(M) > ZERO) THEN
! Flag to calculate radiation.
            CALC_RADT_DES(M) = .TRUE.
         ENDIF
      ENDDO

! Set the value of the Stefan-Boltzman Constant based on the units
      IF(UNITS == 'SI')THEN
         SB_CONST = 5.6704d0*(10.0d0**(-8)) ! W/((m^2).K^4)
      ELSE
         SB_CONST = 1.355282d0*(10.0d0**(-12)) ! cal/((cm^2).sec.K^4)
      ENDIF


! Conduction Equations:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

! At this point, ks_model must be defined for all dem solids
      zero_ks0 = 0
      undef_ks0 = 0
      no_dem = 0   ! distinguish different types of discrete phases
      ANY_CONDUCTION = .TRUE.

      DO M = MMAX+1, MMAX+DES_MMAX
         IF(.NOT.(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP' .OR.SOLIDS_MODEL(M) == 'SQP'.OR.SOLIDS_MODEL(M) == 'GSP')) CYCLE

         no_dem = no_dem + 1
         SELECT CASE (trim(adjustl(KS_MODEL(M))))
         CASE ('MUSSER')
            KS_MODEL_ENUM(M) = ks_musser
            IF (K_S0(M) == UNDEFINED .OR. K_S0(M) < ZERO) THEN
! while k_s0=0 is permissible, it would involve an unnecessary computationally
! intensive for a result that would be zero for any/all conductive heat
! transfer involving a particle of that phase. however, one could permit this
! to allow a case where some particles have ks=0 and others ks>0
               WRITE(ERR_MSG, 1001) trim(iVar('K_s0',M)), iVal(K_s0(M))
               CALL LOG_ERROR()
            ELSEIF (K_S0(M) == ZERO) THEN
               zero_ks0 = zero_ks0+1
            ENDIF
! Flag to calculate conduction.
            CALC_COND_DES(M) = .TRUE.

         CASE ('NONE')
! in this case all particle conductive heat transfer must be neglected. given
! the way particle-particle conduction is resolved (through the neighbor list)
! if any particle conduction is skipped then for consistency all particle
! conduction should be skipped
            KS_MODEL_ENUM(M) = ks_none
            IF (K_S0(M) /= UNDEFINED .AND. K_S0(M) /= ZERO) THEN
! flag warning k_s0 will be overwritten to 0
               WRITE(ERR_MSG, 1201) trim(iVar('K_S0',M)), K_s0(M), &
                  trim(iVar('KS_MODEL',M)), trim(adjustl(KS_MODEL(M)))
               CALL LOG_WARNING()
               K_S0(M) = zero
            ENDIF
            undef_ks0 = undef_ks0 + 1
            CALC_COND_DES(M) = .FALSE.
            ANY_CONDUCTION = .FALSE.
         CASE DEFAULT
            WRITE(ERR_MSG, 1001) trim(adjustl(KS_MODEL(M))), &
               trim(iVar('KS_MODEL',M))
            CALL LOG_ERROR()
         END SELECT
      ENDDO

! check only against dem discrete phases
      IF (zero_ks0>0 .AND. zero_ks0 == no_dem) THEN
          WRITE(ERR_MSG, 1202)
          CALL LOG_WARNING()
      ENDIF
      IF (undef_ks0>0 .AND. undef_ks0 /= no_dem) THEN
         WRITE(ERR_MSG, 1203)
         CALL LOG_ERROR()
      ENDIF


      IF (.NOT.ANY_CONDUCTION) RETURN

! Set the default value for the minimum distance separating particles'
! surfaces.
      IF(DES_MIN_COND_DIST == UNDEFINED)THEN
         DES_MIN_COND_DIST = 1.0D-04 ! cm
         IF (UNITS == 'SI') DES_MIN_COND_DIST = &
            DES_MIN_COND_DIST/100.0  ! m
      ENDIF

! Setup code for conduction correction terms for artificial softening
      ! Calculate masses used for collision calculations.

! Shift the phase index for certain inputs to match the global phase
! index until this matter can be addressed otherwise (i.e., require
! the user specify correct indexing in project settings). This should have no
! impact if not running a hybrid case
      MMAX_TOT = DES_MMAX+MMAX
      e_young_actual((MMAX+1):MMAX_TOT) = e_young_actual(1:DES_MMAX)
      v_poisson_actual((MMAX+1):MMAX_TOT) = v_poisson_actual(1:DES_MMAX)
      DO_AREA_CORRECTION = .TRUE.
      DO M=MMAX+1,MMAX_TOT
         IF(.NOT.(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP'.OR.SOLIDS_MODEL(M) == 'SQP'.OR.SOLIDS_MODEL(M) == 'GSP')) CYCLE
         IF(E_YOUNG_ACTUAL(M) == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Phase ',I2,' Actual EYoungs')") M
            WRITE(ERR_MSG,2002) 'E_YOUNG_ACTUAL', MSG
            DO_AREA_CORRECTION = .FALSE.
            CALL LOG_WARNING()
         ENDIF
         IF(V_POISSON_ACTUAL(M) == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Phase ',I2,' Actual Poisson''s ratio')") M
            WRITE(ERR_MSG,2002) 'V_POISSON_ACTUAL', MSG
            DO_AREA_CORRECTION = .FALSE.
            CALL LOG_WARNING()
         ELSEIF(V_POISSON_ACTUAL(M) > 0.5d0 .OR. &
                V_POISSON_ACTUAL(M) <= -ONE) THEN
            WRITE(ERR_MSG,1001) trim(iVar('V_POISSON_ACTUAL',M)),  &
               iVal(V_POISSON_ACTUAL(M))
            CALL LOG_ERROR()
         ENDIF
      ENDDO

      IF(EW_YOUNG_ACTUAL == UNDEFINED) THEN
         MSG=''; WRITE(MSG,"('Actual EYoungs (wall)')")
         WRITE(ERR_MSG,2002) 'EW_YOUNG_ACTUAL', MSG
         DO_AREA_CORRECTION = .FALSE.
         CALL LOG_WARNING()
      ENDIF
      IF(VW_POISSON_ACTUAL == UNDEFINED) THEN
         MSG=''; WRITE(MSG,"('Actual Poisson''s ratio (wall)')")
         WRITE(ERR_MSG,2002) 'VW_POISSON_ACTUAL', MSG
         DO_AREA_CORRECTION = .FALSE.
         CALL LOG_WARNING()
      ENDIF

      DO M=MMAX+1,MMAX_TOT
         IF(.NOT.(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP'.OR.SOLIDS_MODEL(M) == 'SQP'.OR.SOLIDS_MODEL(M) == 'GSP')) CYCLE
         IF(.NOT.DO_AREA_CORRECTION)CYCLE
! Calculate the mass of a phase M particle.
         ! MASS_M = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)
         MASS_M = REF_MASS(M)

! Particle-Particle Collision Parameters ------------------------------>
         DO L=M,MMAX_TOT
            ! MASS_L = (PI/6.d0)*(D_P0(L)**3)*RO_S0(L)
            MASS_L = REF_MASS(L)
            MASS_EFF = (MASS_M*MASS_L)/(MASS_M+MASS_L)

! Calculate the effective radius, Young's modulus, and shear modulus.
            ! R_EFF = 0.5d0*(D_P0(M)*D_P0(L)/(D_P0(M) + D_P0(L)))
            R_EFF = REF_RADIUS(M)*REF_RADIUS(L)/(REF_RADIUS(M)+REF_RADIUS(L))
            E_EFF = E_YOUNG_ACTUAL(M)*E_YOUNG_ACTUAL(L) /  &
            &  (E_YOUNG_ACTUAL(M)*(1.d0 - V_POISSON_ACTUAL(L)**2) + &
            &   E_YOUNG_ACTUAL(L)*(1.d0 - V_POISSON_ACTUAL(M)**2))

! Calculate the spring properties and store in symmetric matrix format.
            HERT_KN_ACTUAL(M,L)=(4.d0/3.d0)*SQRT(R_EFF)*E_EFF

! Compute baseline for Hertzian collision time
            TAU_C_BASE_ACTUAL(M,L)=3.21D0*(MASS_Eff/HERT_KN_ACTUAL(M,L))**0.4
            ! Can compute actual collision time via:
            !    TAU_C_ACTUAL = TAU_C_BASE_ACTUAL * (1/ImpactVel)^0.2

! Compute base for simulated collision time.  If Hertzian, only include
! base (without impact vel component).
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN)THEN
               TAU_C_BASE_SIM(M,L)=3.21D0*(MASS_Eff/(HERT_KN(M,L)*SQRT(R_EFF)))**0.4
            ELSE
               TAU_C_BASE_SIM(M,L)=PI/SQRT(KN/MASS_EFF - &
               ((DES_ETAN(M,L)/SQRT(MASS_EFF))**2)/4.d0)
            ENDIF
         ENDDO

! Do particle-wall calculations
         MASS_EFF = MASS_M
         ! R_EFF = 0.5d0*D_P0(M)
         R_EFF = REF_RADIUS(M)
         E_EFF = E_YOUNG_ACTUAL(M)*EW_YOUNG_ACTUAL /  &
         &  (E_YOUNG_ACTUAL(M)*(1.d0 - VW_POISSON_ACTUAL**2) + &
         &   EW_YOUNG_ACTUAL*(1.d0 - V_POISSON_ACTUAL(M)**2))

         HERT_KWN_ACTUAL(M) = (4.d0/3.d0)*SQRT(R_EFF)*E_EFF
         TAUW_C_BASE_ACTUAL(M) = 3.21D0 * (MASS_Eff/HERT_KWN_ACTUAL(M))**0.4

         IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN)THEN
            TAUW_C_BASE_SIM(M)=3.21D0*(MASS_Eff/(HERT_KWN(M)*SQRT(R_EFF)))**0.4
         ELSE
            TAUW_C_BASE_SIM(M)=PI/SQRT(KN_w/MASS_EFF - &
            ((DES_ETAN_WALL(M)/SQRT(MASS_EFF))**2)/4.d0)
         ENDIF

      ENDDO

      RETURN

 1000 FORMAT('Required input not specified: ',A)

 1001 FORMAT('Illegal or unknown input: ',A,' = ',A)

 2002 FORMAT('Recommended input not specified: ',A,/       &
         'Description: ',A,/'Not correcting contact area.')

 1201 FORMAT('WARNING 1201: Inconsistent settings: ',A,' = ',g14.4,/,&
        'with ',A,' = ',A, /,'Overwriting this value with zero.')

 1202 FORMAT('WARNING 1202: All K_s0(M)=0 but have requested a KS_MODEL(M) ',&
        'other than',/,'''NONE''. Suggest revising project settings so that ',&
        'KS_MODEL(M)=''NONE'' for ',/,'all M.')

 1203 FORMAT('ERROR 1203: If any KS_MODEL(M)=''NONE'', all must be NONE.')



      END SUBROUTINE CHECK_SOLIDS_DEM_ENERGY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_SOLIDS_DEM_COHESION                              !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check/set parameters for DES cohesion models.              !
!                                                                      !
!  Comments: Original code moved from CHECK_DES_DATA                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_COHESION

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Neighborhood size for Van der Waals force.
      DOUBLE PRECISION :: VDW_NEIGHBORHOOD
!......................................................................!


! Override the following settings if cohesion not used.
      IF(.NOT.USE_COHESION) THEN
!No more square well in the code
         SQUARE_WELL = .FALSE.
         VAN_DER_WAALS = .FALSE.
         WALL_VDW_OUTER_CUTOFF = ZERO
         RETURN
      ENDIF

! Verify that only one cohesion model is specified.
      IF (SQUARE_WELL .AND. VAN_DER_WAALS) THEN
         WRITE(ERR_MSG,1002)
         CALL LOG_ERROR()

! Verify that at a cohesion model is specified.
      ELSEIF(.NOT.SQUARE_WELL .AND. .NOT.VAN_DER_WAALS) THEN
         WRITE(ERR_MSG,1003)
         CALL LOG_ERROR()
      ENDIF


! Van der Waals model checks.
      IF (VAN_DER_WAALS) THEN

         IF (VDW_INNER_CUTOFF .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'VDW_INNER_CUTOFF'
            CALL LOG_ERROR()
         ENDIF

         IF(VDW_OUTER_CUTOFF .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'VDW_OUTER_CUTOFF'
            CALL LOG_ERROR()
         ENDIF

         IF(HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'HAMAKER_CONSTANT'
            CALL LOG_ERROR()
         ENDIF

         IF (WALL_VDW_INNER_CUTOFF .EQ. UNDEFINED)THEN
            WRITE(ERR_MSG,1201) 'WALL_VDW_INNER_CUTOFF'
            CALL LOG_ERROR()
         ENDIF

         IF (WALL_VDW_OUTER_CUTOFF .EQ. UNDEFINED)THEN
            WRITE(ERR_MSG,1201) 'WALL_VDW_OUTER_CUTOFF'
            CALL LOG_ERROR()
         ENDIF

         IF(WALL_HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'WALL_HAMAKER_CONSTANT'
            CALL LOG_ERROR()
         ENDIF

         VDW_NEIGHBORHOOD = 1.0d0 + (VDW_OUTER_CUTOFF/(2.d0*MAX_RADIUS))
         IF (FACTOR_RLM < VDW_NEIGHBORHOOD) THEN
            WRITE(ERR_MSG,1202)
            CALL LOG_ERROR()
         ENDIF

         IF (ASPERITIES < ZERO) THEN
            WRITE(ERR_MSG,1001) 'ASPERITIES', trim(iVal(ASPERITIES))
            CALL LOG_ERROR()
         ENDIF

         SURFACE_ENERGY=HAMAKER_CONSTANT/&
            (24.d0*Pi*VDW_INNER_CUTOFF**2)

         WALL_SURFACE_ENERGY=WALL_HAMAKER_CONSTANT/&
            (24.d0*Pi*WALL_VDW_INNER_CUTOFF**2)

      ENDIF

      RETURN

1001  FORMAT('Illegal or unknown input: ',A, ' = ',A)

1002  FORMAT('Cannot use SQUARE_WELL and VAN_DER_WAALS ',  &
         'cohesion',/'models simultaneously.')

1003  FORMAT('A cohesion model was not selected. Specify ',&
         'one of the available models in the project settings.')


!<------------------- Van der Waals model messages. ----------------->!

1201  FORMAT('Missing input data for Van der Waals ',      &
         'cohesion model.',/'Input parameter ',A,' is UNDEFINED.')

1202  FORMAT('VDW_OUTER_CUTOFF outside of the neighbor ',  &
         'search distance.',/'Increase FACTOR_RLM to increase the ',   &
         'search distance.')

      END SUBROUTINE CHECK_SOLIDS_DEM_COHESION


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_DEM_COLLISION                              !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for DES collision calculations.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_COLLISION

      IMPLICIT NONE

! Check coefficient friction
      IF(MEW == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'MEW'
         CALL LOG_ERROR()
      ELSEIF (MEW < ZERO .OR. MEW_W > ONE) THEN
         WRITE(ERR_MSG,1001) 'MEW', trim(iVal(MEW))
         CALL LOG_ERROR()
      ENDIF

      IF(MEW_W == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'MEW_W'
         CALL LOG_ERROR()
      ELSEIF(MEW_w < ZERO .OR. MEW_W > ONE) THEN
         WRITE(ERR_MSG,1001) 'MEW_W', trim(iVal(MEW_W))
         CALL LOG_ERROR()
      ENDIF

! Check rolling coefficient friction
      IF(MEW_R == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'MEW_R'
         CALL LOG_ERROR()
      ELSEIF (MEW_R < ZERO ) THEN
         WRITE(ERR_MSG,1001) 'MEW_R', trim(iVal(MEW_R))
         CALL LOG_ERROR()
      ENDIF

      IF(MEW_RW == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'MEW_RW'
         CALL LOG_ERROR()
      ELSEIF(MEW_RW < ZERO ) THEN
         WRITE(ERR_MSG,1001) 'MEW_RW', trim(iVal(MEW_RW))
         CALL LOG_ERROR()
      ENDIF

      IF(DTSOLID_INIT == UNDEFINED) THEN
         IF(DTSOLID_FACTOR == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'DTSOLID_FACTOR'
            CALL LOG_ERROR()
         ELSEIF(DTSOLID_FACTOR < 10.0D0) THEN
            WRITE(ERR_MSG,1001) 'DT_SOLID_FACTOR', trim(iVal(DTSOLID_FACTOR))
            CALL LOG_ERROR()
         ENDIF
      ELSEIF(DTSOLID_INIT<=ZERO) THEN
         WRITE(ERR_MSG,1002) 'DTSOLID', DTSOLID_INIT
         CALL LOG_ERROR()
      ENDIF

! Check collision model specific parameters.

! JKR is implemented only for ordinary spherical Hertzian DEM particles.
      IF(USE_JKR .AND. USE_COHESION) THEN
         WRITE(ERR_MSG,'(A)') &
            'USE_JKR and USE_COHESION cannot both be enabled.'
         CALL LOG_ERROR()
      ENDIF

      IF(USE_JKR .AND. TRIM(DES_COLL_MODEL) /= 'HERTZIAN') THEN
         WRITE(ERR_MSG,'(A)') &
            'USE_JKR requires DES_COLL_MODEL = "HERTZIAN".'
         CALL LOG_ERROR()
      ENDIF

      IF(USE_JKR .AND. (SuperDEM .OR. GSP_EXPLICIT .OR. GSP_IMPLICIT)) THEN
         WRITE(ERR_MSG,'(A)') 'USE_JKR supports spherical DEM only.'
         CALL LOG_ERROR()
      ENDIF

      SELECT CASE (trim(DES_COLL_MODEL))
! Linear spring-dashpot model.
      CASE('LSD')
        DES_COLL_MODEL_ENUM = LSD
        CALL CHECK_SOLIDS_DEM_COLL_LSD
! Hertzian collision model.
      CASE('HERTZIAN')
         DES_COLL_MODEL_ENUM = HERTZIAN
         CALL CHECK_SOLIDS_DEM_COLL_HERTZ
! Unknown collision model.
      CASE DEFAULT
         WRITE(ERR_MSG,2000) TRIM(DES_COLL_MODEL)
         CALL LOG_ERROR()
      END SELECT

2000  FORMAT('Invalid particle-particle collision model: ',A)

      RETURN

1000  FORMAT('Required input not specified: ',A)

1001  FORMAT('Illegal or unknown input: ',A,' = ',A)

1002  FORMAT('Invalid (negative or zero) solids time step: ',A,' = ',E14.8)

      END SUBROUTINE CHECK_SOLIDS_DEM_COLLISION


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_DEM_COLL_LSD                               !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for DES collision calculations.      !
!                                                                      !
!  References:                                                         !
!   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
!   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
!      65-149 (pages 94-95)                                            !
!   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_COLL_LSD

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M, L, LC, MMAX_TOT
! Calculate phase index offset for certain inputs until it can be
! addressed in other ways. should not matter unlesss hybrid
      INTEGER :: lent, lend, lenc
! Flag to warn user.
      LOGICAL :: FLAG_WARN1, FLAG_WARN2
! Collision length scale.
      DOUBLE PRECISION :: TCOLL, TCOLL_TMP
! Collision length scale.
      DOUBLE PRECISION :: MASS_M, MASS_L, MASS_EFF
! Alias for coefficient restitution
      DOUBLE PRECISION :: EN
! SQP volume
      DOUBLE PRECISION :: aaxi(3),exp_m,exp_n,svolume
! critical dtsolid from virtual bond
      DOUBLE PRECISION :: TBOND, TMP_DTSOLID
!......................................................................!

! Initialize.
      TCOLL = UNDEFINED
      TBOND = UNDEFINED
      TMP_DTSOLID = 1000000.0

! Check for particle-particle normal spring constants.
      IF(KN == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'KN'
         CALL LOG_ERROR()
      ENDIF

! Check for particle-particle tangential spring constant factors.
      IF(KT_FAC == UNDEFINED) THEN
         WRITE (ERR_MSG, 2100) 'KT_FAC'
         CALL LOG_WARNING()
         KT_FAC = 2.0d0/7.0d0
      ELSEIF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
         WRITE(ERR_MSG,1001) 'KT_FAC', trim(iVal(KT_FAC))
         CALL LOG_ERROR()
      ENDIF
! Calculate the particle-particle tangential spring factor.
      KT = KT_FAC*KN

! Check for particle-wall normal spring constants.
      IF(KN_W == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'KN_W'
         CALL LOG_ERROR()
      ENDIF

! Check for particle-wall tangential spring constant factors.
      IF(KT_W_FAC == UNDEFINED) THEN
         WRITE (ERR_MSG, 2100) 'KT_W_FAC'
         CALL LOG_WARNING()
         KT_W_FAC = 2.0d0/7.0d0
      ELSEIF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
         WRITE(ERR_MSG,1001) 'KT_W_FAC', trim(iVal(KT_W_FAC))
         CALL LOG_ERROR()
      ENDIF
! Calculate the particle-wall tangential spring factor.
      KT_W = KT_W_FAC*KN_W

2100  FORMAT('Tangential spring factor ',A,' not ',      &
         'specified in project settings.',/'Setting to default: (2/7).')

! Check for particle-particle tangential damping coefficients
      IF(DES_ETAT_FAC == UNDEFINED) THEN
         WRITE (ERR_MSG, 2101) 'DES_ETAT_FAC'
         CALL LOG_WARNING()
         DES_ETAT_FAC = HALF
      ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
         WRITE(ERR_MSG,1001) 'DES_ETAT_FAC', iVal(DES_ETAT_FAC)
         CALL LOG_ERROR()
      ENDIF

! Check for particle-wall tangential damping coefficients
      IF(DES_ETAT_W_FAC == UNDEFINED) THEN
         WRITE (ERR_MSG, 2101) 'DES_ETAT_W_FAC'
         CALL LOG_WARNING()
         DES_ETAT_W_FAC = HALF
      ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
         WRITE(ERR_MSG,1001) 'DES_ETAT_W_FAC', iVal(DES_ETAT_W_FAC)
         CALL LOG_ERROR()
      ENDIF

2101  FORMAT('Tangential damping factor ',A,' not ', &
         'specified, setting to default: (1/2).')

! Shift the phase index for certain inputs to match the global phase
! index until this matter can be addressed otherwise (i.e., require
! the user specify correct indexing in project settings). This should have no
! impact if not running a hybrid case
      MMAX_TOT = DES_MMAX+MMAX
      des_en_wall_input((MMAX+1):MMAX_TOT) = des_en_wall_input(1:DES_MMAX)
      des_et_wall_input((MMAX+1):MMAX_TOT) = des_et_wall_input(1:DES_MMAX)
      lent = MMAX_TOT+MMAX_TOT*(MMAX_TOT-1)/2
      lend = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
      lenc = lent-lend
      des_en_input((lenc+1):lent) = des_en_input(1:lend)
      des_et_input((lenc+1):lent) = des_et_input(1:lend)
      LC = lenc

! if below indicated quantities are assigned warn user they are discarded
      FLAG_WARN1 = .FALSE.
      FLAG_WARN2 = .FALSE.

      DO M = MMAX+1, MMAX_TOT
         IF(.NOT.(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP'.OR.SOLIDS_MODEL(M) == 'SQP'.OR.SOLIDS_MODEL(M) == 'GSP')) CYCLE

         if(flexible_gsp .and. allocated(glueDiameter)) TBOND = MIN(TBOND, MEFF_OVER_K(M))

! Calculate the mass of a phase M particle, need future work here
         IF(SuperDEM) THEN
             aaxi(1)= SQP_a(M)
             aaxi(2)= SQP_b(M)
             aaxi(3)= SQP_c(M)
             exp_m=SQP_m(M)
             exp_n=SQP_n(M)
             CALL SQ_VOLUME(aaxi,exp_m,exp_n,svolume)
             MASS_M = (svolume*RO_S0(M))
             sqp_vol(M) = svolume
         else
! Calculate the mass of a phase M particle.
             ! MASS_M = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)
             MASS_M = REF_MASS(M)
         endif

! Particle-Particle Collision Parameters ------------------------------>
         DO L = M, MMAX_TOT

! not necessary to cycle here given requirements on ordering of solids input
! (all dem phase indices are contiguous)
            LC = LC+1

! Check particle-particle normal restitution coefficient
            IF(DES_EN_INPUT(LC) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('DES_EN_INPUT',LC))
               CALL LOG_ERROR()
            ELSEIF(DES_EN_INPUT(LC) > ONE .OR. &
                   DES_EN_INPUT(LC) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('DES_EN_INPUT',LC)), &
                  trim(iVal(DES_EN_INPUT(LC)))
               CALL LOG_ERROR()
            ENDIF
            EN = DES_EN_INPUT(LC)

! warn user that particle-particle tangential restitution coefficient
! is discarded for this model
            IF(DES_ET_INPUT(LC) .NE. UNDEFINED) FLAG_WARN1 = .TRUE.

! Calculate masses used for collision calculations.
            ! MASS_L = (PI/6.d0)*(D_P0(L)**3)*RO_S0(L)
            MASS_L = REF_MASS(L)
            MASS_EFF = MASS_M*MASS_L/(MASS_M+MASS_L)

! Calculate the M-L normal and tangential damping coefficients.
            IF(EN .NE. ZERO) THEN
               DES_ETAN(M,L) = 2.0D0*SQRT(KN*MASS_EFF) * ABS(LOG(EN))
               DES_ETAN(M,L) = DES_ETAN(M,L)/SQRT(PI*PI + (LOG(EN)**2))
            ELSE
               DES_ETAN(M,L) = 2.0D0*SQRT(KN*MASS_EFF)
            ENDIF
!            DES_ETAT(M,L) = DES_ETAT_FAC*DES_ETAN(M,L)

!{SBAN
! Calculate the collision time scale (added conditional).
            IF(EN .NE. ZERO) THEN
               TCOLL_TMP = PI/SQRT(KN/MASS_EFF - &
                  ((DES_ETAN(M,L)/MASS_EFF)**2)/4.d0)
            ELSE
               TCOLL_TMP = PI/SQRT(KN/MASS_EFF)
            ENDIF
            TCOLL = MIN(TCOLL_TMP, TCOLL)

! Convert the entries to particle-agnostic form.
            DES_ETAN(M,L) = DES_ETAN(M,L) / SQRT(MASS_EFF)
            DES_ETAT(M,L) = DES_ETAT_FAC * DES_ETAN(M,L)

! Store the entries in the symmetric matrix in particle-agnostic form.
            DES_ETAN(L,M) = DES_ETAN(M,L)
            DES_ETAT(L,M) = DES_ETAT(M,L)
!}
         ENDDO


! Particle-Wall Collision Parameters ---------------------------------->
! Check particle-wall normal restitution coefficient.
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_EN_WALL_INPUT',M))
            CALL LOG_ERROR()
         ELSEIF(DES_EN_WALL_INPUT(M) > ONE .OR. &
            DES_EN_WALL_INPUT(M) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('DES_EN_WALL_INPUT',M)), &
               trim(iVal(DES_EN_WALL_INPUT(M)))
            CALL LOG_ERROR()
         ENDIF
         EN = DES_EN_WALL_INPUT(M)

         IF(DES_ET_WALL_INPUT(M) .NE. UNDEFINED) FLAG_WARN2 = .TRUE.

! Calculate masses used for collision calculations.
         MASS_EFF = MASS_M

! Calculate the M-Wall normal and tangential damping coefficients.
         IF(EN .NE. ZERO) THEN
            DES_ETAN_WALL(M) = 2.d0*SQRT(KN_W*MASS_EFF)*ABS(LOG(EN))
            DES_ETAN_WALL(M) = DES_ETAN_WALL(M)/SQRT(PI*PI+(LOG(EN))**2)
         ELSE
            DES_ETAN_WALL(M) = 2.D0*SQRT(KN_W*MASS_EFF)
         ENDIF
!         DES_ETAT_WALL(M) = DES_ETAT_W_FAC*DES_ETAN_WALL(M)

!{SBAN
! Calculate the collision time scale (added conditional).
            IF(EN .NE. ZERO) THEN
               TCOLL_TMP = PI/SQRT(KN_W/MASS_EFF - &
                  ((DES_ETAN_WALL(M)/MASS_EFF)**2.d0)/4.d0)
            ELSE
               TCOLL_TMP = PI/SQRT(KN_W/MASS_EFF)
            ENDIF
!            TCOLL = MIN(TCOLL_TMP, TCOLL)

! Convert the entries to particle-agnostic form.
            DES_ETAN_WALL(M) = DES_ETAN_WALL(M) / SQRT(MASS_EFF)
            DES_ETAT_WALL(M) = DES_ETAT_W_FAC * DES_ETAN_WALL(M)
!}
      ENDDO

      IF (FLAG_WARN1) THEN
         WRITE(ERR_MSG,2102) 'DES_ET_INPUT'
         CALL LOG_WARNING()
      ENDIF
      IF (FLAG_WARN2)THEN
         WRITE(ERR_MSG,2102) 'DES_ET_WALL_INPUT'
         CALL LOG_WARNING()
      ENDIF

2102  FORMAT(A,' values are not used with the',  &
         ' linear spring-dashpot collision model.')

! DEM solids time step is computed based on collision time scale.

! When DTSOLID is defined:
!  - Use DTSOLID for the entire simulation if DTSOLID_UPDTATE_DT is undefined.
!  - Use DTSOLID initially, then update it every DTSOLID_UPDATE_DT if
!    DTSOLID_UPDATE_DT is defined.
!
! When DTSOLID is undefined:
!  - Compute DTSOLID at the beginning of the simulation
!  - Recompute DTSOLID every DTSOLID_UPDATE_DT if
!    DTSOLID_UPDATE_DT is defined.

      IF(DTSOLID_INIT==UNDEFINED) THEN
         DTSOLID = TCOLL/DTSOLID_FACTOR

         IF(gsp_explicit .and. allocated(glueDiameter)) THEN
            IF(flexible_gsp) THEN
               ! critical tbond = 2*sqrt(m/k)
               TBOND = 2.0D0 * SQRT(TBOND)
               TBOND = TBOND/DTSOLID_FACTOR ! safety factor
               TMP_DTSOLID = MIN(TMP_DTSOLID, TBOND)
            ENDIF

            DTSOLID = MIN(TMP_DTSOLID, DTSOLID)
         ENDIF

         WRITE(ERR_MSG,1100) TIME, DTSOLID
         CALL LOG_INFO()
      ELSE
         IF(TIME<DTSOLID_UPDATE_DT) THEN

            DTSOLID = DTSOLID_INIT
            WRITE(ERR_MSG,1200) TIME, DTSOLID
            CALL LOG_INFO()

         ELSE
            DTSOLID = TCOLL/DTSOLID_FACTOR

            IF(gsp_explicit .and. allocated(glueDiameter)) THEN
               IF(flexible_gsp) THEN
                  ! critical tbond = 2*sqrt(m/k)
                  TBOND = 2.0D0 * SQRT(TBOND)
                  TBOND = TBOND/DTSOLID_FACTOR ! safety factor
                  TMP_DTSOLID = MIN(TMP_DTSOLID, TBOND)
               ENDIF

               DTSOLID = MIN(TMP_DTSOLID, DTSOLID)
            ENDIF

            WRITE(ERR_MSG,1300) TIME, DTSOLID
            CALL LOG_INFO()
         ENDIF
      ENDIF

      RETURN

1000  FORMAT('Required input not specified: ',A)

1001  FORMAT('Illegal or unknown input: ',A,' = ',A)

1100  FORMAT('Info: At time = ',E14.6, ' sec., DEM solids time step is computed based on particle properties,',/,'DTSOLID (sec) = ', E14.6)

1200  FORMAT('Info: At time = ',E14.6, ' sec., User-defined DEM solids time step is ',/,'DTSOLID (sec) = ', E14.6)

1300  FORMAT('Info: At time = ',E14.6, ' sec., DEM solids time step is updated based on particle properties,',/,'DTSOLID (sec) = ', E14.6)

      END SUBROUTINE CHECK_SOLIDS_DEM_COLL_LSD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_DEM_COLL_HERTZ                             !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for Hertzian collisions.             !
!                                                                      !
!  References:                                                         !
!   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
!   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
!      65-149 (pages 94-95)                                            !
!   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_COLL_HERTZ

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M, L, LC, MMAX_TOT
! Calculate phase index offset for certain inputs until it can be
! addressed in other ways. should not matter unlesss hybrid
      INTEGER :: lent, lend, lenc
! Message for formatted output.
      CHARACTER(len=64) :: MSG
! Collision length scale.
      DOUBLE PRECISION :: TCOLL, TCOLL_TMP
! Particle and effective mass.
      DOUBLE PRECISION :: MASS_M, MASS_L, MASS_EFF
! Effective physical quantities. Radius, Young's, Shear
      DOUBLE PRECISION :: R_EFF, E_EFF, E_EFF_ACTUAL
      DOUBLE PRECISION :: G_MOD_EFF, RED_MASS_EFF
! Alias for coefficient restitution
      DOUBLE PRECISION :: EN, ET
! shear modulus
      DOUBLE PRECISION :: G_MOD(DIM_M)
! Shear modulus for wall
      DOUBLE PRECISION :: G_MOD_WALL
! SQP volume
      DOUBLE PRECISION :: aaxi(3),exp_m,exp_n,svolume
! critical dtsolid from virtual bond
      DOUBLE PRECISION :: TBOND, TMP_DTSOLID
!......................................................................!

! Initialize.
      TCOLL = UNDEFINED
      TBOND = UNDEFINED
      TMP_DTSOLID = 1000000.0
      MMAX_TOT = DES_MMAX+MMAX

      IF(USE_JKR) THEN
         IF(EW_YOUNG_ACTUAL == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'EW_YOUNG_ACTUAL'
            CALL LOG_ERROR()
         ELSEIF(EW_YOUNG_ACTUAL <= ZERO .OR. &
            EW_YOUNG_ACTUAL < EW_YOUNG) THEN
            WRITE(ERR_MSG,1001) 'EW_YOUNG_ACTUAL', &
               iVal(EW_YOUNG_ACTUAL)
            CALL LOG_ERROR()
         ENDIF

         IF(DES_SURFACE_ENERGY_WALL == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'DES_SURFACE_ENERGY_WALL'
            CALL LOG_ERROR()
         ELSEIF(.NOT.(DES_SURFACE_ENERGY_WALL > ZERO)) THEN
            WRITE(ERR_MSG,1001) 'DES_SURFACE_ENERGY_WALL', &
               iVal(DES_SURFACE_ENERGY_WALL)
            CALL LOG_ERROR()
         ENDIF

         DO M=1,DES_MMAX
            IF(E_YOUNG_ACTUAL(M) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('E_YOUNG_ACTUAL',M))
               CALL LOG_ERROR()
            ELSEIF(E_YOUNG_ACTUAL(M) <= ZERO .OR. &
               E_YOUNG_ACTUAL(M) < E_YOUNG(M)) THEN
               WRITE(ERR_MSG,1001) trim(iVar('E_YOUNG_ACTUAL',M)), &
                  iVal(E_YOUNG_ACTUAL(M))
               CALL LOG_ERROR()
            ENDIF

            IF(DES_SURFACE_ENERGY(M) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) 'DES_SURFACE_ENERGY'
               CALL LOG_ERROR()
            ELSEIF(.NOT.(DES_SURFACE_ENERGY(M) > ZERO)) THEN
               WRITE(ERR_MSG,1001) 'DES_SURFACE_ENERGY', &
                  iVal(DES_SURFACE_ENERGY(M))
               CALL LOG_ERROR()
            ENDIF
         ENDDO

         IF(.NOT.(JKR_ROLL_DAMP_COEFF >= ZERO) .OR. &
            JKR_ROLL_DAMP_COEFF > HUGE(ONE)) THEN
            WRITE(ERR_MSG,1001) 'JKR_ROLL_DAMP_COEFF', &
               iVal(JKR_ROLL_DAMP_COEFF)
            CALL LOG_ERROR()
         ENDIF
         IF(JKR_THETAR_CRIT < 1.0D-4 .OR. &
            JKR_THETAR_CRIT > 1.0D-1) THEN
            WRITE(ERR_MSG,1001) 'JKR_THETAR_CRIT', iVal(JKR_THETAR_CRIT)
            CALL LOG_ERROR()
         ENDIF
         IF(WALL_JKR_THETAR_CRIT < 1.0D-4 .OR. &
            WALL_JKR_THETAR_CRIT > 1.0D-1) THEN
            WRITE(ERR_MSG,1001) 'WALL_JKR_THETAR_CRIT', &
               iVal(WALL_JKR_THETAR_CRIT)
            CALL LOG_ERROR()
         ENDIF
         IF(JKR_ALPHA_SAFETY < 0.7D0 .OR. JKR_ALPHA_SAFETY > ONE) THEN
            WRITE(ERR_MSG,1001) 'JKR_ALPHA_SAFETY', iVal(JKR_ALPHA_SAFETY)
            CALL LOG_ERROR()
         ENDIF

         E_YOUNG_ACTUAL((MMAX+1):MMAX_TOT) = &
            E_YOUNG_ACTUAL(1:DES_MMAX)
      ENDIF


! check young's modulus and poisson ratio
      IF(Ew_YOUNG == UNDEFINED ) THEN
         MSG="Wall value for Young's modulus"
         WRITE(ERR_MSG,1002) 'Ew_YOUNG', MSG
         CALL LOG_ERROR()
      ENDIF

      IF(Vw_POISSON == UNDEFINED) THEN
         MSG="Wall value for Poisson's ratio"
         WRITE(ERR_MSG,1002) 'Vw_POISSON', MSG
         CALL LOG_ERROR()
      ELSEIF (Vw_POISSON > 0.5d0 .OR. Vw_POISSON <= -ONE) THEN
         WRITE(ERR_MSG,1001) 'Vw_POISSON',iVal(Vw_POISSON)
         CALL LOG_ERROR()
      ENDIF

      G_MOD_WALL = 0.5d0*Ew_YOUNG/(1.d0+Vw_POISSON)

! Shift the phase index for certain inputs to match the global phase
! index until this matter can be addressed otherwise (i.e., require
! the user specify correct indexing in project settings). This should have no
! impact if not running a hybrid case
      e_young((MMAX+1):MMAX_TOT) = e_young(1:DES_MMAX)
      v_poisson((MMAX+1):MMAX_TOT) = v_poisson(1:DES_MMAX)
      des_en_wall_input((MMAX+1):MMAX_TOT) = des_en_wall_input(1:DES_MMAX)
      des_et_wall_input((MMAX+1):MMAX_TOT) = des_et_wall_input(1:DES_MMAX)
      IF(USE_JKR) DES_SURFACE_ENERGY((MMAX+1):MMAX_TOT) = &
         DES_SURFACE_ENERGY(1:DES_MMAX)
      lent = MMAX_TOT+MMAX_TOT*(MMAX_TOT-1)/2
      lend = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
      lenc = lent-lend
      des_en_input((lenc+1):lent) = des_en_input(1:lend)
      des_et_input((lenc+1):lent) = des_et_input(1:lend)

      DO M=MMAX+1,MMAX_TOT
         IF(.NOT.(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP'.OR.SOLIDS_MODEL(M) == 'SQP'.OR.SOLIDS_MODEL(M) == 'GSP')) CYCLE

         IF(E_YOUNG(M) == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Phase ',I2,' Young''s modulus')") M
            WRITE(ERR_MSG,1002) 'E_YOUNG', MSG
            CALL LOG_ERROR()
         ENDIF
         IF(V_POISSON(M) == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Phase ',I2,' Poisson''s ratio')") M
            WRITE(ERR_MSG,1002) 'V_POISSON', MSG
            CALL LOG_ERROR()
         ELSEIF(V_POISSON(M) > 0.5d0 .OR. &
                V_POISSON(M) <= -ONE) THEN
            WRITE(ERR_MSG,1001) trim(iVar('V_POISSON',M)),  &
               iVal(V_POISSON(M))
            CALL LOG_ERROR()
         ENDIF
! Calculate the shear modulus for phase M.
         G_MOD(M) = 0.5d0*E_YOUNG(M)/(1.d0+V_POISSON(M))
      ENDDO

! see above index shift
      LC = LENC
      DO M=MMAX+1,MMAX_TOT
         IF(.NOT.(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP'.OR.SOLIDS_MODEL(M) == 'SQP'.OR.SOLIDS_MODEL(M) == 'GSP')) CYCLE

! Calculate the mass of a phase M particle, need future work here
         IF(SuperDEM) THEN
             aaxi(1)= SQP_a(M)
             aaxi(2)= SQP_b(M)
             aaxi(3)= SQP_c(M)
             exp_m=SQP_m(M)
             exp_n=SQP_n(M)
             CALL SQ_VOLUME(aaxi,exp_m,exp_n,svolume)
             MASS_M = (svolume*RO_S0(M))
             sqp_vol(M) = svolume
         else
! Calculate the mass of a phase M particle.
             ! MASS_M = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)
             MASS_M = REF_MASS(M)
         endif
! Particle-Particle Collision Parameters ------------------------------>
         DO L=M,MMAX_TOT
! not necessary to cycle here given requirements on ordering of solids input
            LC = LC+1

! Check particle-particle normal restitution coefficient
            IF(DES_EN_INPUT(LC) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('DES_EN_INPUT',LC))
               CALL LOG_ERROR()
            ELSEIF(DES_EN_INPUT(LC) > ONE .OR. &
               DES_EN_INPUT(LC) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('DES_EN_INPUT',LC)), &
                  trim(iVal(DES_EN_INPUT(LC)))
               CALL LOG_ERROR()
            ENDIF
            EN = DES_EN_INPUT(LC)

! Check particle-particle tangential restitution coefficient
            IF(DES_ET_INPUT(LC) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('DES_ET_INPUT',LC))
               CALL LOG_ERROR()
            ELSEIF(DES_ET_INPUT(LC) > ONE .OR. &
                   DES_ET_INPUT(LC) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('DES_ET_INPUT',LC)), &
                  iVal(DES_ET_INPUT(LC))
               CALL LOG_ERROR()
            ENDIF
            ET = DES_ET_INPUT(LC)

! Calculate masses used for collision calculations.
         IF(SuperDEM) THEN
              aaxi(1)= SQP_a(L)
              aaxi(2)= SQP_b(L)
              aaxi(3)= SQP_c(L)
              exp_m=SQP_m(L)
              exp_n=SQP_n(L)
             CALL SQ_VOLUME(aaxi,exp_m,exp_n,svolume)
             MASS_L = (svolume*RO_S0(L))
             MASS_EFF = MASS_M*MASS_L/(MASS_M+MASS_L)
          ELSE
! Calculate masses used for collision calculations.
            ! MASS_L = (PI/6.d0)*(D_P0(L)**3)*RO_S0(L)
            MASS_L = REF_MASS(L)
            MASS_EFF = (MASS_M*MASS_L)/(MASS_M+MASS_L)
          ENDIF
            RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF
! Calculate the effective radius, Young's modulus, and shear modulus.


! note that R_EFF is corrected for superquadric particle in calc_for_dem.f
!            R_EFF = 0.5d0*(D_P0(M)*D_P0(L)/(D_P0(M) + D_P0(L)))
            R_EFF = REF_RADIUS(M)*REF_RADIUS(L)/(REF_RADIUS(M)+REF_RADIUS(L))
            E_EFF = E_YOUNG(M)*E_YOUNG(L) /  &
               (E_YOUNG(M)*(1.d0 - V_POISSON(L)**2) + &
                E_YOUNG(L)*(1.d0 - V_POISSON(M)**2))
            G_MOD_EFF = G_MOD(M)*G_MOD(L)/ &
               (G_MOD(M)*(2.d0 - V_POISSON(L)) + &
                G_MOD(L)*(2.d0 - V_POISSON(M)))

! Calculate the spring properties and store in symmetric matrix format.
            HERT_KN(M,L)=(4.d0/3.d0)*SQRT(R_EFF)*E_EFF
            HERT_KT(M,L)= 8.d0*SQRT(R_EFF)*G_MOD_EFF

            HERT_KN(L,M) = HERT_KN(M,L)
            HERT_KT(L,M) = HERT_KT(M,L)

! Calculate the normal damping coefficient.
            IF(EN .NE. ZERO) THEN
               DES_ETAN(M,L) = SQRT(5.d0)*SQRT(HERT_KN(M,L)*MASS_EFF)* &
                  ABS(LOG(EN))
               DES_ETAN(M,L) = DES_ETAN(M,L)/ &
                  SQRT(PI*PI + (LOG(EN))**2)
            ELSE
               DES_ETAN(M,L) = SQRT(5.d0)*SQRT(HERT_KN(M,L)*MASS_EFF)
            ENDIF
            DES_ETAN(L,M) = DES_ETAN(M,L)

! Calculate the tangential coefficients.
            IF(ET .NE. ZERO) THEN
               DES_ETAT(M,L) = SQRT(5.d0)*SQRT(HERT_KT(M,L)*RED_MASS_EFF)* &
                  ABS(LOG(ET))
               DES_ETAT(M,L) = DES_ETAT(M,L)/ SQRT(PI*PI+(LOG(ET))**2)
            ELSE
               DES_ETAT(M,L) = SQRT(5.d0)*SQRT(HERT_KT(M,L)*RED_MASS_EFF)
            ENDIF
            DES_ETAT(L,M) = DES_ETAT(M,L)

! JKR coefficients retain their own contact-radius runtime contract:
! KN/KT are multiplied by a, while damping is multiplied by sqrt(a).
            IF(USE_JKR) THEN
               JKR_GAMAN(M,L) = SQRT(DES_SURFACE_ENERGY(M)* &
                  DES_SURFACE_ENERGY(L))
               JKR_AoN(M,L) = (9.D0*PI*JKR_GAMAN(M,L)*R_EFF**2/ &
                  E_EFF)**(ONE/3.D0)
               JKR_DeltacN(M,L) = JKR_AoN(M,L)**2 / &
                  (2.D0*6.D0**(ONE/3.D0)*R_EFF)
               JKR_FcN(M,L) = 3.D0*PI*R_EFF*JKR_GAMAN(M,L)

               E_EFF_ACTUAL = E_YOUNG_ACTUAL(M)*E_YOUNG_ACTUAL(L) / &
                  (E_YOUNG_ACTUAL(M)*(ONE - V_POISSON(L)**2) + &
                   E_YOUNG_ACTUAL(L)*(ONE - V_POISSON(M)**2))
               JKR_GAMAN_ACTUAL(M,L) = JKR_GAMAN(M,L)* &
                  (E_EFF_ACTUAL/E_EFF)**(2.D0/5.D0)
               JKR_FcN_ACTUAL(M,L) = 3.D0*PI*R_EFF* &
                  JKR_GAMAN_ACTUAL(M,L)

               HERT_KN_JKR(M,L) = (4.D0/3.D0)*E_EFF
               HERT_KT_JKR(M,L) = 8.D0*G_MOD_EFF
               DES_ETAN_JKR(M,L) = JKR_DAMPING_ALPHA(EN)* &
                  SQRT(HERT_KN_JKR(M,L)*MASS_EFF)
               DES_ETAT_JKR(M,L) = DES_ETAN_JKR(M,L)

               JKR_GAMAN(L,M) = JKR_GAMAN(M,L)
               JKR_AoN(L,M) = JKR_AoN(M,L)
               JKR_DeltacN(L,M) = JKR_DeltacN(M,L)
               JKR_FcN(L,M) = JKR_FcN(M,L)
               JKR_GAMAN_ACTUAL(L,M) = JKR_GAMAN_ACTUAL(M,L)
               JKR_FcN_ACTUAL(L,M) = JKR_FcN_ACTUAL(M,L)
               HERT_KN_JKR(L,M) = HERT_KN_JKR(M,L)
               HERT_KT_JKR(L,M) = HERT_KT_JKR(M,L)
               DES_ETAN_JKR(L,M) = DES_ETAN_JKR(M,L)
               DES_ETAT_JKR(L,M) = DES_ETAT_JKR(M,L)
            ENDIF

!{SBAN
! Calculate the collision time scale (added conditional).
            IF(EN .NE. ZERO) THEN
               TCOLL_TMP = PI/SQRT(HERT_KN(M,L)/MASS_EFF - &
               ((DES_ETAN(M,L)/MASS_EFF)**2)/5.d0)
            ELSE
               TCOLL_TMP = PI/SQRT(HERT_KN(M,L)/MASS_EFF)
            ENDIF
            TCOLL = MIN(TCOLL_TMP, TCOLL)

! Convert the entries to particle-agnostic form.
            DES_ETAN(M,L) = DES_ETAN(M,L) / SQRT(HERT_KN(M,L)*MASS_EFF)
            DES_ETAT(M,L) = DES_ETAT(M,L) / SQRT(HERT_KT(M,L)*MASS_EFF)
            HERT_KN(M,L) = HERT_KN(M,L) / SQRT(R_EFF)
            HERT_KT(M,L) = HERT_KT(M,L) / SQRT(R_EFF)

! Store the entries in the symmetric matrix in particle-agnostic form.
            HERT_KN(L,M) = HERT_KN(M,L)
            HERT_KT(L,M) = HERT_KT(M,L)
            DES_ETAN(L,M) = DES_ETAN(M,L)
            DES_ETAT(L,M) = DES_ETAT(M,L)
!}
         ENDDO

! Particle-Wall Collision Parameters ---------------------------------->
! Check particle-wall normal restitution coefficient.
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_EN_WALL_INPUT',M))
            CALL LOG_ERROR()
         ELSEIF(DES_EN_WALL_INPUT(M) > ONE .OR. &
                DES_EN_WALL_INPUT(M) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('DES_EN_WALL_INPUT',M)), &
               trim(iVal(DES_EN_WALL_INPUT(M)))
            CALL LOG_ERROR()
         ENDIF
         EN = DES_EN_WALL_INPUT(M)

! Check particle-wall tangential restitution coefficient
         IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_ET_WALL_INPUT',M))
            CALL LOG_ERROR()
         ELSEIF(DES_ET_WALL_INPUT(M) > ONE .OR. &
                DES_ET_WALL_INPUT(M) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('DES_ET_WALL_INPUT',M)), &
               trim(iVal(DES_ET_WALL_INPUT(M)))
            CALL LOG_ERROR()
         ENDIF
         ET = DES_ET_WALL_INPUT(M)

! Calculate masses used for collision calculations.
         MASS_EFF = MASS_M
         RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF
! Calculate the effective radius, Young's modulus, and shear modulus.
         ! R_EFF = 0.5d0*D_P0(M)
         R_EFF = REF_RADIUS(M)
         E_EFF = E_YOUNG(M)*Ew_YOUNG /                                 &
            (E_YOUNG(M)*(1.d0-Vw_POISSON**2) +                         &
             Ew_YOUNG  *(1.d0-V_POISSON(M)**2))
         G_MOD_EFF = G_MOD(M)*G_MOD_WALL /                             &
            (G_MOD(M)*(2.d0 - Vw_POISSON) +                            &
             G_MOD_WALL*(2.d0 - V_POISSON(M)))

! Calculate the spring properties.
         HERT_Kwn(M) = (4.d0/3.d0)*SQRT(R_EFF)*E_EFF
         HERT_Kwt(M) = 8.0*SQRT(R_EFF)*G_MOD_EFF

! Calculate the tangential coefficients.
         IF(EN /= ZERO) THEN
            DES_ETAN_WALL(M) = SQRT(5.d0)*SQRT(HERT_Kwn(M)*MASS_EFF)*&
               ABS(LOG(EN))
            DES_ETAN_WALL(M) = DES_ETAN_WALL(M)/&
               SQRT(PI*PI + (LOG(EN))**2)
         ELSE
            DES_ETAN_WALL(M) = SQRT(5.d0)*SQRT(HERT_Kwn(M)*MASS_EFF)
         ENDIF

         IF(ET /= ZERO) THEN
            DES_ETAT_WALL(M) = SQRT(5.d0)*SQRT(HERT_Kwt(M)*RED_MASS_EFF)*    &
                ABS(LOG(ET))
            DES_ETAT_WALL(M) = DES_ETAT_WALL(M)/SQRT(PI*PI+(LOG(ET))**2)
         ELSE
            DES_ETAT_WALL(M) = SQRT(5.d0)*SQRT(HERT_Kwt(M)*RED_MASS_EFF)
         ENDIF

! Particle-wall JKR coefficients use the same contact-radius contract.
         IF(USE_JKR) THEN
            JKR_GAMAWN(M) = SQRT(DES_SURFACE_ENERGY(M)* &
               DES_SURFACE_ENERGY_WALL)
            JKR_AoWN(M) = (9.D0*PI*JKR_GAMAWN(M)*R_EFF**2/ &
               E_EFF)**(ONE/3.D0)
            JKR_DeltacWN(M) = JKR_AoWN(M)**2 / &
               (2.D0*6.D0**(ONE/3.D0)*R_EFF)
            JKR_FcWN(M) = 3.D0*PI*R_EFF*JKR_GAMAWN(M)
            E_EFF_ACTUAL = E_YOUNG_ACTUAL(M)*EW_YOUNG_ACTUAL / &
               (E_YOUNG_ACTUAL(M)*(ONE-VW_POISSON**2) + &
                EW_YOUNG_ACTUAL*(ONE-V_POISSON(M)**2))
            JKR_GAMAWN_ACTUAL(M) = JKR_GAMAWN(M)* &
               (E_EFF_ACTUAL/E_EFF)**(2.D0/5.D0)
            JKR_FcWN_ACTUAL(M) = 3.D0*PI*R_EFF* &
               JKR_GAMAWN_ACTUAL(M)
            HERT_KWN_JKR(M) = (4.D0/3.D0)*E_EFF
            HERT_KWT_JKR(M) = 8.D0*G_MOD_EFF
            DES_ETAN_WALL_JKR(M) = JKR_DAMPING_ALPHA(EN)* &
               SQRT(HERT_KWN_JKR(M)*MASS_EFF)
            DES_ETAT_WALL_JKR(M) = DES_ETAN_WALL_JKR(M)
         ENDIF

!{SBAN
! Calculate the collision time scale (added conditional).
            IF(EN .NE. ZERO) THEN
               TCOLL_TMP = PI/SQRT(HERT_Kwn(M)/MASS_EFF - &
               ((DES_ETAN_WALL(M)/MASS_EFF)**2)/5.d0)
            ELSE
               TCOLL_TMP = PI/SQRT(HERT_KN(M,L)/MASS_EFF)
            ENDIF

! Convert the entries to particle-agnostic form.
            DES_ETAN_WALL(M) = DES_ETAN_WALL(M) / SQRT(HERT_Kwn(M)*MASS_EFF)
            DES_ETAT_WALL(M) = DES_ETAT_WALL(M) / SQRT(HERT_Kwt(M)*MASS_EFF)
            HERT_Kwn(M) = HERT_Kwn(M) / SQRT(R_EFF)
            HERT_Kwt(M) = HERT_Kwt(M) / SQRT(R_EFF)
!}
      ENDDO


! Store the smalled calculated collision time scale. This value is used
! in time-marching the DEM solids.

! When DTSOLID is defined:
!  - Use DTSOLID for the entire simulation if DTSOLID_UPDTATE_DT is undefined.
!  - Use DTSOLID initially, then update it every DTSOLID_UPDATE_DT if
!    DTSOLID_UPDATE_DT is defined.
!
! When DTSOLID is undefined:
!  - Compute DTSOLID at the beginning of the simulation
!  - Recompute DTSOLID every DTSOLID_UPDATE_DT if
!    DTSOLID_UPDATE_DT is defined.

      IF(DTSOLID_INIT==UNDEFINED) THEN
         DTSOLID = TCOLL/DTSOLID_FACTOR

         IF(gsp_explicit .and. allocated(glueDiameter)) THEN
            IF(flexible_gsp) THEN
               ! critical tbond = 2*sqrt(m/k)
               TBOND = 2.0D0 * SQRT(TBOND)
               TBOND = TBOND/DTSOLID_FACTOR ! safety factor
               TMP_DTSOLID = MIN(TMP_DTSOLID, TBOND)
            ENDIF

            DTSOLID = MIN(TMP_DTSOLID, DTSOLID)
         ENDIF

         WRITE(ERR_MSG,1100) TIME, DTSOLID
         CALL LOG_INFO()
      ELSE
         IF(TIME<DTSOLID_UPDATE_DT) THEN
            DTSOLID = DTSOLID_INIT
            WRITE(ERR_MSG,1200) TIME, DTSOLID
            CALL LOG_INFO()
         ELSE
            DTSOLID = TCOLL/DTSOLID_FACTOR

            IF(gsp_explicit .and. allocated(glueDiameter)) THEN
               IF(flexible_gsp) THEN
                  ! critical tbond = 2*sqrt(m/k)
                  TBOND = 2.0D0 * SQRT(TBOND)
                  TBOND = TBOND/DTSOLID_FACTOR ! safety factor
                  TMP_DTSOLID = MIN(TMP_DTSOLID, TBOND)
               ENDIF

               DTSOLID = MIN(TMP_DTSOLID, DTSOLID)
            ENDIF

            WRITE(ERR_MSG,1300) TIME, DTSOLID
            CALL LOG_INFO()
         ENDIF
      ENDIF

      RETURN

1000  FORMAT('Required input not specified: ',A)

1001  FORMAT('Illegal or unknown input: ',A,' = ',A)

1002  FORMAT('Required input not specified: ',A,/          &
         'Description: ',A)

1100  FORMAT('Info: At time = ',E14.6, ' sec., DEM solids time step is computed based on particle properties,',/,'DTSOLID (sec) = ', E14.6)

1200  FORMAT('Info: At time = ',E14.6, ' sec., DEM solids time step is set to',/,'DTSOLID (sec) = ', E14.6)

1300  FORMAT('Info: At time = ',E14.6, ' sec., DEM solids time step is updated based on particle properties,',/,'DTSOLID (sec) = ', E14.6)

   END SUBROUTINE CHECK_SOLIDS_DEM_COLL_HERTZ

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GET_REF_PARTICLE_SIZE_AND_MASS                          !
!  Author: J.-F. Dietiker                             Date: 11-Sep-20  !
!                                                                      !
!  Purpose: Get reference particle size and mass for each DEM phase    !
!           This will be used to compute the DTSOLID at the beginning  !
!           of the simulation and at regular intervals (when particle  !
!           size/mass changes in time.                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GET_REF_PARTICLE_SIZE_AND_MASS

      use functions, only: IS_NORMAL
      USE mpi_utility, only: global_all_min
      use run, only: RUN_TYPE, TIME

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M, L, LC, MMAX_TOT, NP, gid
! Reference radius and mass used to compute solids time step
      DOUBLE PRECISION :: lradius_min(DIM_M), lmass_min(DIM_M), cur_gsp_radius
! Message for formatted output.
      CHARACTER(len=64) :: MSG
! Math
      DOUBLE PRECISION :: THIRD = 1.0/3.0
!......................................................................!

      INTEGER :: CC_START, CC_END, CC, pII, nbid
      DOUBLE PRECISION :: meff_k
      DOUBLE PRECISION :: lmok_min(DIM_M)
! Initialize with baseline values.

      MMAX_TOT = DES_MMAX+MMAX
      MEFF_OVER_K(:) = undefined

      DO M=MMAX+1,MMAX_TOT
         IF(.NOT.(SOLIDS_MODEL(M) == 'DEM'.OR.SOLIDS_MODEL(M) == 'CGP'.OR.SOLIDS_MODEL(M) == 'SQP'.OR.SOLIDS_MODEL(M) == 'GSP')) CYCLE

         REF_RADIUS(M) = HALF * D_P0(M)
         REF_MASS(M)   = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)

         IF(gsp_explicit .AND. ALLOCATED(glueDiameter)) THEN
            ! loop on components and find the attached gid
            DO NP = 1, MAX_PIP
               IF ((PIJK(NP,5) .NE. M) .OR. .NOT.IS_NORMAL(NP)) CYCLE

               ! INTEGER :: CC_START, CC_END, CC, pII, nbid

               IF(flexible_gsp) THEN

                  CC_START = 1
                  IF (NP.gt.1) CC_START = NEIGHBOR_INDEX(NP-1)
                  CC_END   = NEIGHBOR_INDEX(NP)
                  cc_loop: DO CC = CC_START, CC_END-1
                     pII  = NEIGHBORS(CC)
                     IF(gid_list(NP) .ne. gid_list(pII)) CYCLE cc_loop
                     DO nbid = 1, 6
                        IF(iglobal_id(pII) == gp_neighsid(NP,nbid)) THEN
                           meff_k = 0.0
                           CALL DETERMINE_BOND_DT(NP, pII, M, meff_k)
                           ! find the minimum m_eff/bond_stiffness for each phase
                           MEFF_OVER_K(M)= min(MEFF_OVER_K(M), meff_k)
                        ENDIF
                     ENDDO
                  ENDDO cc_loop
               ENDIF ! end of flexible GSP

               gid = gid_list(NP)
               cur_gsp_radius = HALF * glueDiameter(gid)
               ! reference radius is determined by gsp equivalent diameter instead of bounding diameter per phase.
               REF_RADIUS(M) = min(REF_RADIUS(M), cur_gsp_radius)
            ENDDO
            REF_MASS(M) = (4.d0 * PI/3.d0)*(REF_RADIUS(M)**3)*RO_S0(M)
         ENDIF

         IF(gsp_implicit .AND. ALLOCATED(cspart_data)) THEN
            ! loop on components and find the attached gid
            DO NP = 1, MAX_PIP
               IF ((PIJK(NP,5) .NE. M) .OR. .NOT.IS_NORMAL(NP)) CYCLE
               cur_gsp_radius = HALF * (6.0D0*PVOL(NP)/PI)**THIRD
               ! reference radius is determined by gsp equivalent diameter instead of bounding diameter per phase.
               REF_RADIUS(M) = min(REF_RADIUS(M), cur_gsp_radius)
            ENDDO
            REF_MASS(M) = (4.d0 * PI/3.d0)*(REF_RADIUS(M)**3)*RO_S0(M)
         ENDIF

      ENDDO

! Get the min particle size for each solids phase
      lradius_min(MMAX+1: MMAX_TOT) = REF_RADIUS(MMAX+1: MMAX_TOT) !UNDEFINED
      lmass_min(MMAX+1: MMAX_TOT) = REF_MASS(MMAX+1: MMAX_TOT) !UNDEFINED
      if(flexible_gsp) lmok_min(MMAX+1: MMAX_TOT) = MEFF_OVER_K(MMAX+1: MMAX_TOT)

      IF(GluedSphereDEM) THEN
         ! call MPI to share reference values
         call global_all_min(lradius_min(MMAX+1: MMAX_TOT), REF_RADIUS(MMAX+1: MMAX_TOT))
         call global_all_min(lmass_min(MMAX+1: MMAX_TOT), REF_MASS(MMAX+1: MMAX_TOT))
         if(flexible_gsp) call global_all_min(lmok_min(MMAX+1: MMAX_TOT), MEFF_OVER_K(MMAX+1: MMAX_TOT))
         RETURN
         ! force to return, do not access the component loop
         ! so for gsp model the solid time step is always depending on glueDiameter
         ! call check_solid_dem in make_arrays_dem.f to make sure update the DTSOLID based on GSP equivalent diameter
      ENDIF

      IF(PIP>0) THEN

          DO NP=1,MAX_PIP
             IF(.NOT.IS_NORMAL(NP)) CYCLE
! Solids phase
                M              = PIJK(NP,5)
                lradius_min(M) = MIN(lradius_min(M),DES_RADIUS(NP))
                lmass_min(M)   = MIN(lmass_min(M),PMASS(NP))
          ENDDO

      ENDIF

      call global_all_min(lradius_min(MMAX+1: MMAX_TOT), REF_RADIUS(MMAX+1: MMAX_TOT))
      call global_all_min(lmass_min(MMAX+1: MMAX_TOT), REF_MASS(MMAX+1: MMAX_TOT))

      RETURN

1000  FORMAT('Required input not specified: ',A)

   END SUBROUTINE GET_REF_PARTICLE_SIZE_AND_MASS


   SUBROUTINE UPDATE_DTSOLID

   use run, only: time
   use des_bc

   IMPLICIT NONE

   DOUBLE PRECISION :: dtsolid_copy
   INTEGER :: BCV_I

   dtsolid_copy = dtsolid
   IF(DTSOLID_UPDATE_DT<=ZERO.OR.DTSOLID_UPDATE_DT==UNDEFINED) RETURN
   IF(DTSOLID_UPDATE_TIME<=ZERO) RETURN

! Update DEM time step
   IF(TIME>=DTSOLID_UPDATE_TIME) THEN
      CALL CHECK_SOLIDS_DEM

! Set time for the next update
      DTSOLID_UPDATE_TIME = DTSOLID_UPDATE_TIME + DTSOLID_UPDATE_DT

! The number of particles injected into the system from BCs per DTSOLID need to be modified
! Loop over BCs that flagged for DEM mass inflow.
      DO BCV_I = 1, DEM_BCMI
         IF(PI_COUNT(BCV_I) .EQ. 1) THEN
            PI_FACTOR(BCV_I) = real(PI_FACTOR(BCV_I)*DTSOLID_COPY/DTSOLID)
         ELSEIF(PI_COUNT(BCV_I) .GT. 1) THEN
            PI_COUNT(BCV_I) = CEILING(real(PI_COUNT(BCV_I)/DTSOLID_COPY*DTSOLID))
         ENDIF
      ENDDO

   ENDIF


   RETURN
   END SUBROUTINE UPDATE_DTSOLID

   SUBROUTINE DETERMINE_BOND_DT(NP, pII, PID, Meff_K)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NP, pII
   INTEGER, INTENT(IN) :: PID
   DOUBLE PRECISION, INTENT(INOUT) :: Meff_K

   DOUBLE PRECISION :: length_vector(3)
   DOUBLE PRECISION :: initial_bond_length, rad_eff, area_bond, mass_eff
   DOUBLE PRECISION :: G_bond
   DOUBLE PRECISION :: current_min
   DOUBLE PRECISION :: moi_bond, pmoi_bond, stiff_bond_bending, stiff_bond_norm, stiff_bond_tan, stiff_bond_twisting

      length_vector(:) = sc2gpc_vec(pII,:) - sc2gpc_vec(NP,:)
      initial_bond_length = dot_product(length_vector,length_vector)
      initial_bond_length = sqrt(initial_bond_length)

      IF(initial_bond_length == 0.0) THEN
         WRITE(*,*) "Error: Initially overlapping components"
         ERROR STOP
      ENDIF

      current_min = undefined

      ! calculate the effective radius
      rad_eff = 2.0D0 * des_radius(pII)*des_radius(NP)/(des_radius(pII) + des_radius(NP))
      ! cross section of the bond, the bond is assumed to be a cylinder
      area_bond = pi*rad_eff**(2.0)
      ! calculate the effective mass
      mass_eff = 2.0D0 * pmass(pII)*pmass(NP)/(pmass(pII) + pmass(NP))

      ! calculate the bond normal stiffness, this is based on provided young's modulus
      stiff_bond_norm = E_YOUNG_BOND(PID)*area_bond/initial_bond_length
      current_min = min(current_min, mass_eff/stiff_bond_norm)

      ! shear modulus is derived from provided young's modulus
      G_bond = E_YOUNG_BOND(PID)/(2.0d0 * (1.0d0 + V_POISSON_BOND(PID)))
      ! calculate the bond shear stiffness based on the shear modulus
      stiff_bond_tan = G_bond*area_bond/initial_bond_length
      current_min = min(current_min, mass_eff/stiff_bond_tan)

      ! area moment of inertia, bending module
      moi_bond = pi*rad_eff**(4.0d0)/4.0d0
      stiff_bond_bending = E_YOUNG_BOND(PID)*moi_bond/initial_bond_length
      current_min = min(current_min, mass_eff/stiff_bond_bending)

      ! polar moment of inertia, twisting module
      pmoi_bond = pi*rad_eff**(4.0d0)/2.0d0
      stiff_bond_twisting = G_bond*pmoi_bond/initial_bond_length
      current_min = min(current_min, mass_eff/stiff_bond_twisting)

      Meff_K = current_min


   END SUBROUTINE DETERMINE_BOND_DT

END MODULE CHECK_SOLIDS_DEM_MOD
