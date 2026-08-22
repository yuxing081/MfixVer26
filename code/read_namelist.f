#include "error.inc"

MODULE READ_NAMELIST_MOD

      USE bc
      USE cdist
      USE compar
      USE constant
      USE cutcell
      USE dashboard
      USE deprecated_mod, only: deprecated_or_unknown
      USE des_bc
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE error_manager
      USE fldvar
      USE fpe
      USE funits
      USE geometry
      USE ic
      USE indices
      USE is
      USE iterate, only: max_nit
      USE keyframe_mod
      USE leqsol
      USE make_upper_case_mod, only: make_upper_case, replace_tab
      USE mfix_pic
      USE monitor
      USE output
      USE parallel
      USE param1, only: undefined
      USE parse, only: arrhenius_rrates_fluid, arrhenius_rrates_des
      USE parse_line_mod, only: parse_line
      USE particle_filter
      USE physprop
      USE polygon
      USE ps
      USE qmom_kinetic_equation
      USE quadric
      USE randomno, only : repeatable
      USE remove_comment_mod, only: remove_comment
      USE remove_comment_mod, only: remove_par_blanks
      USE residual
      USE residual_pub
      USE run
      USE rxns
      USE scalars
      USE scales
      USE stiff_chem
      USE toleranc
      USE ur_facs
      USE usr
      USE utilities
      USE vtk
      Use stl
      use turb, only: turbulence_model
      use usr_prop, only: usr_fgs, usr_fss, usr_gama
      use usr_prop, only: usr_rog, usr_cpg, usr_mug, usr_kg, usr_difg
      use usr_prop, only: usr_ros, usr_cps, usr_mus, usr_ks, usr_difs
      use usr_src, only: call_usr_source
      use visc_g, only: mu_gmax

! Coefficient of restitution (old symbol)
   DOUBLE PRECISION e

! Indicates whether currently reading rxns or rate
   LOGICAL :: RXN_FLAG

! Flags restricting what keywords to process
   LOGICAL :: READ_LOCKED, READ_FULL

! Holds one line in the input
   CHARACTER(LEN=1024) :: LINE_STRING
! Length of noncomment string
   INTEGER :: LINE_LEN
! Line number
   INTEGER :: LINE_NO

   INTEGER :: IOS, II

! enumerator for READ_KEYWORD_NAMELIST to read certain keywords ahead of the others (to initialize error logs)
   enum, bind(c)
      enumerator ENUM_ENABLE_DMP_LOG
      enumerator ENUM_FULL_LOG
      enumerator ENUM_RUN_NAME
   end enum

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: READ_NAMELIST(POST)                                 !
!     Author: P. Nicoletti                            Date: 25-NOV-91  !
!                                                                      !
!     Purpose: Read in the NAMELIST variables                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE READ_NAMELIST(READ_ACTION, MFIX_DAT, CMD_LINE_ARGS, CMD_LINE_ARGS_COUNT)

      IMPLICIT NONE

! Dummy Arguments:
!------------------------------------------------------------------------//
! Specify how much of the input to process.
      INTEGER, INTENT(IN) :: READ_ACTION

! Filename of the input file
      CHARACTER(LEN=*), INTENT(IN) :: MFIX_DAT

      CHARACTER(LEN=80), DIMENSION(100), INTENT(IN), OPTIONAL :: CMD_LINE_ARGS
      INTEGER, INTENT(IN), OPTIONAL :: CMD_LINE_ARGS_COUNT

! Local Parameters:
!---------------------------------------------------------------------//
      INTEGER, PARAMETER :: READ_MFIX = 0
      INTEGER, PARAMETER :: READ_POST = 1
      INTEGER, PARAMETER :: READ_REINIT = 2

      E = UNDEFINED
      RXN_FLAG = .FALSE.
      NO_OF_RXNS = 0

      SELECT CASE(READ_ACTION)
      CASE(READ_MFIX)
         READ_LOCKED = .TRUE.
         READ_FULL = .TRUE.
      CASE(READ_POST)
         READ_LOCKED = .TRUE.
         READ_FULL = .FALSE.
      CASE(READ_REINIT)
         READ_LOCKED = .FALSE.
         READ_FULL = .TRUE.
      END SELECT

      CALL READ_KEYWORD_NAMELIST(MFIX_DAT)

      IF (PRESENT(CMD_LINE_ARGS_COUNT) .AND. PRESENT(CMD_LINE_ARGS)) THEN
         CALL READ_CMDLINE_KEYWORDS(CMD_LINE_ARGS, CMD_LINE_ARGS_COUNT)
      ENDIF

      IF (E /= UNDEFINED) C_E = E

      RETURN

   END SUBROUTINE READ_NAMELIST


   !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
   !                                                                      !
   ! Subroutine: READ_KEYWORD_NAMELIST                                    !
   ! Author: Mark Meredith                               Date: 31-Jan-18  !
   !                                                                      !
   ! Purpose: Process MFiX keywords from an input file.                   !
   !                                                                      !
   !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE READ_KEYWORD_NAMELIST(MFIX_DAT, enum_kw)
      CHARACTER(LEN=*), INTENT(IN) :: MFIX_DAT
      INTEGER, INTENT(IN), OPTIONAL :: ENUM_KW
      INTEGER :: START_INDEX, END_INDEX, PRE_KEYWORD

! LINE_STRING(1:MAXCOL) has valid input data
      INTEGER, PARAMETER :: MAXCOL = 300

! Error flag
      LOGICAL :: ERROR

! temp string
      CHARACTER(len=256) :: STRING

      NAMELIST / RUN_NAME_NAMELIST / RUN_NAME
      NAMELIST / OUTPUT_CONTROL_UNLOCKED / FULL_LOG
      NAMELIST / OUTPUT_CONTROL_LOCKED / ENABLE_DMP_LOG

      IF(PRESENT(ENUM_KW)) THEN
         pre_keyword = enum_kw
      ELSE
         pre_keyword = -1
      ENDIF

      LINE_NO = 0

! Loop through the lines of MFIX_DAT and process the input data.
      START_INDEX = 1
      READ_LP: DO
         END_INDEX = INDEX(MFIX_DAT(START_INDEX:), NEW_LINE(''))
         IF (END_INDEX == 0) THEN
            LINE_STRING = MFIX_DAT(START_INDEX:)
         ELSE
            LINE_STRING = MFIX_DAT(START_INDEX:START_INDEX+END_INDEX-2)
         END IF
         START_INDEX = START_INDEX + END_INDEX
         IF (END_INDEX == 0) EXIT READ_LP ! end of file

         LINE_NO = LINE_NO + 1

         LINE_LEN = SEEK_COMMENT(LINE_STRING,LEN(trim(LINE_STRING))) - 1
         CALL REMOVE_COMMENT(LINE_STRING, LINE_LEN+1, LEN(trim(LINE_STRING)))

         IF(LINE_LEN <= 0) CYCLE READ_LP           ! comment line
         IF(BLANK_LINE(LINE_STRING)) CYCLE READ_LP ! blank line

         IF(LINE_TOO_BIG(LINE_STRING,LINE_LEN,MAXCOL) > 0) THEN
            WRITE (ERR_MSG, 1100) trim(iVAL(LINE_NO)), trim(ival(MAXCOL)), &
               LINE_STRING(1:MAXCOL)
            CALL LOG_ERROR()
         ENDIF

1100     FORMAT(//1X,/1x,'From: READ_NAMELIST',/1x,'Error 1100: ', &
            'Line ',A,' in input file is too long. Input lines should', &
            /1x,'not pass column ',A,'.',2/3x,A,2/1x)

         SELECT CASE (pre_keyword)

         CASE (ENUM_RUN_NAME)
            CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN)
            IF(INDEX(LINE_STRING, 'RUN_NAME') > 0) THEN
               STRING=''; STRING = '&RUN_NAME_NAMELIST '//&
                  trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
               READ(STRING, NML=RUN_NAME_NAMELIST, IOSTAT=IOS)
               IF(IOS == 0)  RETURN
            END IF

         CASE (ENUM_FULL_LOG)
            CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN)
            IF(INDEX(LINE_STRING, 'FULL_LOG') > 0) THEN
               STRING=''; STRING = '&OUTPUT_CONTROL_UNLOCKED '//&
                  trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
               READ(STRING, NML=OUTPUT_CONTROL_UNLOCKED, IOSTAT=IOS)
               IF(IOS == 0)  RETURN
            END IF

         CASE (ENUM_ENABLE_DMP_LOG)
            CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN)
            IF(INDEX(LINE_STRING, 'ENABLE_DMP') > 0) THEN
               STRING=''; STRING = '&OUTPUT_CONTROL_LOCKED '//&
                  trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
               READ(STRING, NML=OUTPUT_CONTROL_LOCKED, IOSTAT=IOS)
               IF(IOS == 0)  RETURN
            END IF

         CASE DEFAULT

! All subsequent lines are thermochemical data
            IF(LINE_STRING(1:11) == 'THERMO DATA') EXIT READ_LP

            CALL SET_KEYWORD(ERROR)
            IF (ERROR) THEN
! At this point, the keyword was not identified therefore it is
! either deprecated or unknown.
               CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
            ENDIF

         END SELECT

      ENDDO READ_LP

   END SUBROUTINE READ_KEYWORD_NAMELIST


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: READ_CMDLINE_KEYWORDS                                    !
! Author: Mark Meredith                               Date: 20-NOV-17  !
!                                                                      !
! Purpose: Process MFiX keywords from the command line.                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
! is this used?  can we remove?  cgw 2020-05-18
   SUBROUTINE READ_CMDLINE_KEYWORDS(CMD_LINE_ARGS, CMD_LINE_ARGS_COUNT)
      CHARACTER(LEN=80), DIMENSION(100), INTENT(IN) :: CMD_LINE_ARGS
      INTEGER, INTENT(IN) :: CMD_LINE_ARGS_COUNT
! Error flag
      LOGICAL :: ERROR

      DO II=1, CMD_LINE_ARGS_COUNT
         LINE_STRING = CMD_LINE_ARGS(ii)
         LINE_LEN = len(line_string)
         CALL SET_KEYWORD(ERROR)
         IF (ERROR) THEN
            CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
         ENDIF
      ENDDO

   END SUBROUTINE READ_CMDLINE_KEYWORDS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: SET_KEYWORD(ERROR)                                       !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Process LINE_STRING for MFIX keyword data.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE SET_KEYWORD(ERROR)

      IMPLICIT NONE

      LOGICAL, INTENT(OUT) ::ERROR

!Temporary variable for checking keywords
      CHARACTER(LEN=512) :: LINE_STRING_TMP

! Indicate whether to do a namelist read on the line
      LOGICAL :: READ_FLAG

      CHARACTER(len=256) :: STRING

! External namelist files:
!---------------------------------------------------------------------//
      INCLUDE 'run_control.inc'
      INCLUDE 'physical_params.inc'
      INCLUDE 'numerical_params.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'gas_phase.inc'
      INCLUDE 'solids_phase.inc'
      INCLUDE 'tfm_solids.inc'
      INCLUDE 'initial_conditions.inc'
      INCLUDE 'boundary_conditions.inc'
      INCLUDE 'internal_surfaces.inc'
      INCLUDE 'point_sources.inc'
      INCLUDE 'output_control.inc'
      INCLUDE 'usr_hooks.inc'
      INCLUDE 'chem_equations.inc'
      INCLUDE 'dmp_batch_control.inc'
      INCLUDE 'desnamelist.inc'
      INCLUDE 'cartesian_grid_namelist.inc'
      INCLUDE 'qmomknamelist.inc'
      INCLUDE 'legacy.inc'
      INCLUDE 'monitors.inc'

      ERROR = .FALSE.

! Make upper case every line except species names
      LINE_STRING_TMP = LINE_STRING
      CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN)

      IF(INDEX(LINE_STRING,'SPECIES_NAME') /= 0 .OR. &
         INDEX(LINE_STRING,'SPECIES_G') /= 0 .OR.    &
         INDEX(LINE_STRING,'SPECIES_S') /= 0) THEN
         LINE_STRING = LINE_STRING_TMP
      ENDIF

      CALL REPLACE_TAB (LINE_STRING, LINE_LEN)
      CALL REMOVE_PAR_BLANKS(LINE_STRING)

! Complete arithmetic operations and expand line
      CALL PARSE_LINE (LINE_STRING, LINE_LEN, RXN_FLAG, READ_FLAG)

! Write the current line to a scratch file
! and read the scratch file in NAMELIST format
      IF(.NOT.READ_FLAG) RETURN


! Run control keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&RUN_CONTROL_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=RUN_CONTROL_LOCKED,  IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&RUN_CONTROL_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=RUN_CONTROL_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Physical parameter keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&PHYSICAL_PARAM_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=PHYSICAL_PARAM_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&PHYSICAL_PARAM_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=PHYSICAL_PARAM_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Numerical parameter keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&NUMERICAL_PARAM_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=NUMERICAL_PARAM_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&NUMERICAL_PARAM_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=NUMERICAL_PARAM_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Geometry and discretization keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&GEOMETRY_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=GEOMETRY_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&GEOMETRY_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GEOMETRY_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Gas phase keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&GAS_PHASE_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=GAS_PHASE_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&GAS_PHASE_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GAS_PHASE_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Solids phase keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&SOLIDS_PHASE_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=SOLIDS_PHASE_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&SOLIDS_PHASE_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=SOLIDS_PHASE_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Two-fluid solids keywords
      STRING=''; STRING = '&TFM_SOLIDS_LOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=TFM_SOLIDS_LOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN

      STRING=''; STRING = '&TFM_SOLIDS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=TFM_SOLIDS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Initial condition keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&INITIAL_CONDITIONS_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=INITIAL_CONDITIONS_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

!      STRING=''; STRING = '&INITIAL_CONDITIONS_UNLOCKED '//&
!         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
!      READ(STRING, NML=INITIAL_CONDITIONS_UNLOCKED, IOSTAT=IOS)
!      IF(IOS == 0)  RETURN


! Boundary condition keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&BOUNDARY_CONDITIONS_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=BOUNDARY_CONDITIONS_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&BOUNDARY_CONDITIONS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=BOUNDARY_CONDITIONS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Internal surface keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&INTERNAL_SURFACES_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=INTERNAL_SURFACES_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&INTERNAL_SURFACES_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=INTERNAL_SURFACES_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Point source keywords
      STRING=''; STRING = '&POINT_SOURCES_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=POINT_SOURCES_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Output control keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&OUTPUT_CONTROL_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=OUTPUT_CONTROL_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&OUTPUT_CONTROL_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=OUTPUT_CONTROL_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! User hook keywords
      STRING=''; STRING = '&USER_HOOKS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=USER_HOOKS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Chemical equation keywords
!      IF(READ_LOCKED) THEN
!         STRING=''; STRING = '&CHEM_EQUATIONS_LOCKED '//&
!            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
!        READ(STRING, NML=CHEM_EQUATIONS_LOCKED, IOSTAT=IOS)
!         IF(IOS == 0)  RETURN
!     ENDIF

      STRING=''; STRING = '&CHEM_EQUATIONS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=CHEM_EQUATIONS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! DMP and Batch Queue control keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&DMP_BATCH_CONTROL_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=DMP_BATCH_CONTROL_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&DMP_BATCH_CONTROL_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=DMP_BATCH_CONTROL_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Legacy keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&LEGACY_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=LEGACY_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN

      ENDIF

!     Stop processing keyword inputs if running POSTMFIX
      IF(.NOT.READ_FULL) RETURN

! Monitor keywords - all unlocked

      STRING=''; STRING = '&MONITORS_UNLOCKED '//&
      trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=MONITORS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0) RETURN


      IF(READ_LOCKED) THEN

! Discrete Element model input parameters - all locked
         STRING=''; STRING = '&DES_INPUT_DATA_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=DES_INPUT_DATA_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN

         CALL USR_READ_NAMELIST(STRING, LINE_STRING, LINE_LEN, IOS)
         IF(IOS == 0)  RETURN

! Cartesian grid cut-cell input parameters - all locked
         STRING=''; STRING = '&CARTESIAN_GRID_INPUT_DATA_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=CARTESIAN_GRID_INPUT_DATA_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN


! QMOMK input parameters - all locked
         STRING=''; STRING = '&QMOMK_INPUT_DATA_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=QMOMK_INPUT_DATA_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      IF(READ_LOCKED) ERROR = .TRUE.

      RETURN
   END SUBROUTINE SET_KEYWORD

END MODULE READ_NAMELIST_MOD
