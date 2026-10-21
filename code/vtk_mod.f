MODULE vtk

   use cdist, only: bdist_io
   USE compar, only: myPE, pe_io
   Use param1
! Maximum number of solids phases.
   use param, only: DIM_M
! Maximum number of gas phase species
   use param, only: DIM_N_g
! Maximum number of solids phase species
   use param, only: DIM_N_s
! Maximum number of scalar equations
   use param, only: DIM_Scalar
! Maximum number of DEM solids phase species
   use param, only: DIMENSION_N_S
! Number of reactions rates
   use rxns, only: NO_OF_RXNS, SAVE_FLUID_RRATES, SAVE_DES_RRATES, nrrmax
   use rxns, only: FLUID_RRATES_OUT, DES_RRATES_OUT
   use des_rxns, only: NO_OF_DES_RXNS, SAVE_PART_RRATES, PART_RRATES_OUT

   INTEGER NUMBER_OF_CELLS
   INTEGER NUMBER_OF_CUT_CELLS
   INTEGER NUMBER_OF_BLOCKED_CELLS
   INTEGER NUMBER_OF_STANDARD_CELLS
   INTEGER NUMBER_OF_VTK_CELLS
   INTEGER NUMBER_OF_VTK_PARTS

   LOGICAL             :: WRITE_VTK_FILES
   LOGICAL             :: TIME_DEPENDENT_FILENAME
   LOGICAL             :: RESET_FRAME_AT_TIME_ZERO = .TRUE.
   CHARACTER(LEN=255)  :: VTU_DIR
   CHARACTER(LEN=255)  :: VTK_FILENAME, FRAME_CHAR, VTU_FILENAME, PVD_FILENAME, PVTU_FILENAME
   CHARACTER(LEN=255)  :: VTP_FILENAME, PVTP_FILENAME
   CHARACTER(LEN=255)  :: VTU_FRAME_FILENAME = 'VTU_FRAME_INDEX.TXT'
   CHARACTER(LEN=512) :: BUFFER

   CHARACTER(LEN=1), PARAMETER  :: END_REC = CHAR(10)

   INTEGER, PARAMETER :: DIM_VTK_VAR = 20
   INTEGER, DIMENSION(DIM_VTK_VAR) :: VTK_VAR

   INTEGER :: POLY_COUNTER, NUMBER_OF_POINTS

   integer, allocatable :: CLEANED_CONNECTIVITY(:, :)
   REAL, allocatable :: COORDS_OF_POINTS(:, :)

   integer, allocatable :: GLOBAL_I_OF(:)
   integer, allocatable :: GLOBAL_J_OF(:)
   integer, allocatable :: GLOBAL_K_OF(:)
   integer, allocatable :: GLOBAL_CONNECTIVITY(:, :)
   integer, allocatable :: GLOBAL_CLEANED_CONNECTIVITY(:, :)
   integer, allocatable :: GLOBAL_NUMBER_OF_NODES(:)

   REAL, allocatable :: GLOBAL_COORDS_OF_POINTS(:, :)

   LOGICAL, allocatable :: GLOBAL_INTERIOR_CELL_AT(:)
   LOGICAL, allocatable :: GLOBAL_BLOCKED_CELL_AT(:)
   LOGICAL, allocatable :: GLOBAL_STANDARD_CELL_AT(:)
   LOGICAL, allocatable :: GLOBAL_CUT_CELL_AT(:)
   LOGICAL, allocatable :: GLOBAL_SMALL_CELL_AT(:)
   LOGICAL, allocatable :: GLOBAL_SNAP(:)
   DOUBLE PRECISION, allocatable :: GLOBAL_F_AT(:)
   INTEGER, allocatable :: GLOBAL_BC_ID(:)
   INTEGER, allocatable :: GLOBAL_FLAG(:)


   double precision, allocatable :: GLOBAL_X_NEW_POINT(:)
   double precision, allocatable :: GLOBAL_Y_NEW_POINT(:)
   double precision, allocatable :: GLOBAL_Z_NEW_POINT(:)

   INTEGER :: GLOBAL_NUMBER_OF_NEW_POINTS

   LOGICAL :: GLOBAL_VAR_ALLOCATED

   LOGICAL :: GRID_INFO_PRINTED_ON_SCREEN

   LOGICAL :: WRITE_ANI_CUTCELL

   INTEGER(KIND=8) :: VTU_offset
   INTEGER(KIND=8), PARAMETER :: VTK_HEADER_SIZE = 8  ! UInt64 appended-block size-prefix width (bytes)

   LOGICAL, allocatable :: BELONGS_TO_VTK_SUBDOMAIN(:)
! LOGICAL, allocatable :: PART_BELONGS_TO_VTK_SUBDOMAIN(:)

   INTEGER, PARAMETER :: DIMENSION_VTK = 100
! Max size for VTK_PART_USR_VAR
   INTEGER, PARAMETER :: VTK_PART_USRmax = 100

! Current VTK region
   INTEGER :: VTK_REGION

! Time interval at which vtk files are saved
   DOUBLE PRECISION :: VTK_DT(0:DIMENSION_VTK)

! Type of GSP output shape in VTK files.
   INTEGER :: VTK_GSP_OUTPUT_TYPE(0:DIMENSION_VTK)

! Starting time at which vtk files are saved
   DOUBLE PRECISION :: VTK_T_START(0:DIMENSION_VTK)

! Ending time at which vtk files are saved
   DOUBLE PRECISION :: VTK_T_END(0:DIMENSION_VTK)

! Flag to distinguish between regular file (VTK_DBG_FILE=.FALSE., default value)
!                        and debug file (VTK_DBG_FILE=.TRUE.)
   LOGICAL :: VTK_DBG_FILE(0:DIMENSION_VTK)

! Current vtk time
   DOUBLE PRECISION :: VTK_TIME(0:DIMENSION_VTK)

! Type of data in vtk region: 'C':cell data, 'P': particle data, 'G': Geometry data,
! 'F': Force chain, 'N': Non-spherical (Glued) Particle
   CHARACTER(LEN=1) :: VTK_DATA(0:DIMENSION_VTK)

! Max size  of geometry files in a single VTK region
   INTEGER, PARAMETER :: VTK_GEO_max = 100
! Geometry fine name when using VTK_DATA='G'
   CHARACTER(LEN=20) :: VTK_GEO(0:DIMENSION_VTK, VTK_GEO_max)

! STL GROUP ID in the current VTK REGION
   INTEGER :: VTK_GEO_GROUP

! FRAME index of vtk file
   INTEGER :: FRAME(0:DIMENSION_VTK)

! PVD file initialization flag
   LOGICAL :: PVD_FILE_INITIALIZED(0:DIMENSION_VTK) = .FALSE.

! Logical variable to determine whether an vtk region is defined
   LOGICAL :: VTK_DEFINED(0:DIMENSION_VTK)

! VTK region West face, X-coordinate
   DOUBLE PRECISION :: VTK_X_w(0:DIMENSION_VTK)

! VTK region East face, X-coordinate
   DOUBLE PRECISION :: VTK_X_e(0:DIMENSION_VTK)

! VTK region South face, Y-coordinate
   DOUBLE PRECISION :: VTK_Y_s(0:DIMENSION_VTK)

! VTK region North face, Y-coordinate
   DOUBLE PRECISION :: VTK_Y_n(0:DIMENSION_VTK)

! VTK region Bottom face, Z-coordinate
   DOUBLE PRECISION :: VTK_Z_b(0:DIMENSION_VTK)

! VTK region Top face, Z-coordinate
   DOUBLE PRECISION :: VTK_Z_t(0:DIMENSION_VTK)

! VTK number of slices in x-direction
   INTEGER :: VTK_NXS(0:DIMENSION_VTK)

! VTK number of slices in y-direction
   INTEGER :: VTK_NYS(0:DIMENSION_VTK)

! VTK number of slices in z-direction
   INTEGER :: VTK_NZS(0:DIMENSION_VTK)

! VTK slice tolerance
   DOUBLE PRECISION :: VTK_SLICE_TOL(0:DIMENSION_VTK)

! Flag to write only cut cell data in VTK file
   LOGICAL :: VTK_CUTCELL_ONLY(0:DIMENSION_VTK)

! VTK filename base
   CHARACTER(LEN=255) :: VTK_FILEBASE(0:DIMENSION_VTK)

! Gas phase volume fraction
   LOGICAL :: VTK_EP_g(0:DIMENSION_VTK)

! Gas density
   LOGICAL :: VTK_RO_G(0:DIMENSION_VTK)

! Gas viscosity
   LOGICAL :: VTK_MU_G(0:DIMENSION_VTK)

! Gas thermal conductivity
   LOGICAL :: VTK_K_G(0:DIMENSION_VTK)

! Gas pressure
   LOGICAL :: VTK_P_g(0:DIMENSION_VTK)

! X-component of gas velocity
   LOGICAL :: VTK_U_g(0:DIMENSION_VTK)

! X-component of solids phase velocity
   LOGICAL :: VTK_U_s(0:DIMENSION_VTK, DIM_M)

! Y-component of gas velocity
   LOGICAL :: VTK_V_g(0:DIMENSION_VTK)

! Y-component of solids phase velocity
   LOGICAL :: VTK_V_s(0:DIMENSION_VTK, DIM_M)

! Z-component of gas velocity
   LOGICAL :: VTK_W_g(0:DIMENSION_VTK)

! Z-component of solids phase velocity
   LOGICAL :: VTK_W_s(0:DIMENSION_VTK, DIM_M)

! Gas velocity vector
   LOGICAL :: VTK_VEL_g(0:DIMENSION_VTK)

! Solids pressure preventing overpacking
   LOGICAL :: VTK_P_star(0:DIMENSION_VTK)

! Solids pressure as a result of granular motion
   LOGICAL :: VTK_P_s(0:DIMENSION_VTK, DIM_M)

! Solids velocity vector
   LOGICAL :: VTK_VEL_s(0:DIMENSION_VTK, DIM_M)

! Macroscopic density of solids phases
   LOGICAL :: VTK_ROP_s(0:DIMENSION_VTK, DIM_M)

! Density of solids phases
   LOGICAL :: VTK_RO_s(0:DIMENSION_VTK, DIM_M)

! Solids phase volume fraction
   LOGICAL :: VTK_EP_s(0:DIMENSION_VTK, DIM_M)

! Gas temperature
   LOGICAL :: VTK_T_g(0:DIMENSION_VTK)

! Solids temperature
   LOGICAL :: VTK_T_s(0:DIMENSION_VTK, DIM_M)

! Gas specific heat
   LOGICAL :: VTK_C_Pg(0:DIMENSION_VTK)

! Solids specific heat
   LOGICAL :: VTK_C_Ps(0:DIMENSION_VTK, DIM_M)

! Gas mixture molecular weight
   LOGICAL :: VTK_MW_MIX_g(0:DIMENSION_VTK)

! Gas species mass fractions
   LOGICAL :: VTK_X_g(0:DIMENSION_VTK, DIM_N_g)

! Gas species molar fractions
   LOGICAL :: VTK_Y_g(0:DIMENSION_VTK, DIM_N_g)

! Solids species mass fractions
   LOGICAL :: VTK_X_s(0:DIMENSION_VTK, DIM_M, DIM_N_s)

! Granular temperature
   LOGICAL :: VTK_Theta_m(0:DIMENSION_VTK, DIM_M)

! Scalar value
   LOGICAL :: VTK_Scalar(0:DIMENSION_VTK, DIM_scalar)

! Reaction rates
   LOGICAL :: VTK_RRate(0:DIMENSION_VTK, nRRmax)
   LOGICAL :: VTK_FLUID_RRate(0:DIMENSION_VTK, nRRmax)
   LOGICAL :: VTK_DES_RRate(0:DIMENSION_VTK, nRRmax)

! Reaction rates labels
   CHARACTER(LEN=255) :: VTK_RRate_label(0:DIMENSION_VTK, nRRmax)

! K & Epsilon values
   LOGICAL :: VTK_K_Turb_G(0:DIMENSION_VTK)
   LOGICAL :: VTK_E_Turb_G(0:DIMENSION_VTK)

! Vorticity magnitude
   LOGICAL :: VTK_VORTICITY(0:DIMENSION_VTK)

! Lambda_2 (vortex core indicator)
   LOGICAL :: VTK_LAMBDA_2(0:DIMENSION_VTK)

! Grid partition
   LOGICAL :: VTK_PARTITION(0:DIMENSION_VTK)

! Boundary condition ID
   LOGICAL :: VTK_BC_ID(0:DIMENSION_VTK)

! Wall distance
   LOGICAL :: VTK_DWALL(0:DIMENSION_VTK)

! Facet count (DES)
   LOGICAL :: VTK_FACET_COUNT_DES(0:DIMENSION_VTK)

! Neighboring facets (DES)
   LOGICAL :: VTK_NB_FACET_DES(0:DIMENSION_VTK)

! Cell IJK index
   LOGICAL :: VTK_IJK(0:DIMENSION_VTK)

! Scalar cell center coordinates
   LOGICAL :: VTK_SCALAR_CC(0:DIMENSION_VTK)

! Cut face normal vector
   LOGICAL :: VTK_NORMAL(0:DIMENSION_VTK)

! Debug variable
   LOGICAL :: VTK_DEBUG(0:DIMENSION_VTK, 15)

   INTEGER, DIMENSION(0:DIMENSION_VTK, DIM_VTK_VAR) :: VTK_VARLIST

! VTK Particle selection: 'C': center, 'P': entire particle, 'I': particle intersect
   CHARACTER(LEN=1) :: VTK_SELECT_MODE(0:DIMENSION_VTK)

! Phase ID to use as labels in vtp files when there is only
! one phase in the vtk region
   INTEGER :: VTK_PHASE_FOR_DES_X_S(0:DIMENSION_VTK)

! Particle phase

   LOGICAL :: VTK_PART_PHASE(0:DIMENSION_VTK, DIM_M)

! Particle diameter (Parcel diameter for PIC, Coarse grain diameter for CGDEM)
   LOGICAL :: VTK_PART_DIAMETER(0:DIMENSION_VTK)

! Particle physical diameter for CGDEM
   LOGICAL :: VTK_PART_PHYSICAL_DIAMETER(0:DIMENSION_VTK)

! Particle statistical weight for CGDEM
   LOGICAL :: VTK_PART_CGP_STAT_WT(0:DIMENSION_VTK)

! Particle velocity
   LOGICAL :: VTK_PART_VEL(0:DIMENSION_VTK)

! Particle angular velocity
   LOGICAL :: VTK_PART_ANGULAR_VEL(0:DIMENSION_VTK)

! Particle orientation
   LOGICAL :: VTK_PART_ORIENTATION(0:DIMENSION_VTK)

! Particle user-defined variable
   LOGICAL :: VTK_PART_USR_VAR(0:DIMENSION_VTK, VTK_PART_USRmax)

! Particle temperature
   LOGICAL :: VTK_PART_TEMP(0:DIMENSION_VTK)

! Particle specific heat
   LOGICAL :: VTK_PART_C_Ps(0:DIMENSION_VTK)

! Particle species mass fraction
   LOGICAL :: VTK_PART_X_s(0:DIMENSION_VTK, 100)

! Particle density
   LOGICAL :: VTK_PART_DENSITY(0:DIMENSION_VTK)

! Particle mass
   LOGICAL :: VTK_PART_MASS(0:DIMENSION_VTK)

! Particle volume
   LOGICAL :: VTK_PART_VOLUME(0:DIMENSION_VTK)

! Particle drag coefficient
   LOGICAL :: VTK_PART_DRAG_COEFF(0:DIMENSION_VTK)

! Particle cohesion
   LOGICAL :: VTK_PART_COHESION(0:DIMENSION_VTK)

! Particle rank
   LOGICAL :: VTK_PART_RANK(0:DIMENSION_VTK)

! Particle ID
   LOGICAL :: VTK_PART_ID(0:DIMENSION_VTK)

! Residence time
   LOGICAL :: VTK_PART_RESIDENCE_TIME(0:DIMENSION_VTK)

! Particle Phase ID
   LOGICAL :: VTK_PART_PHASE_ID(0:DIMENSION_VTK)

! Domain decomposition
   LOGICAL :: VTK_DOMAIN_DECOMPOSITION(0:DIMENSION_VTK)

! Reaction rate
   LOGICAL :: VTK_PART_RRate(0:DIMENSION_VTK, nRRmax)

! Reaction rate
   CHARACTER(LEN=64) :: VTK_PART_GROUP_STRING(0:DIMENSION_VTK)

! Maximum force chain
   LOGICAL :: VTK_FCMAX_ONLY(0:DIMENSION_VTK)
! JKR force-chain fields
   LOGICAL :: VTK_FC_JKR_STATE(0:DIMENSION_VTK)
   LOGICAL :: VTK_FC_JKR_RADIUS(0:DIMENSION_VTK)
   LOGICAL :: VTK_FC_JKR_GAMMA(0:DIMENSION_VTK)
   LOGICAL :: VTK_FC_TANG_DISP(0:DIMENSION_VTK)
   LOGICAL :: VTK_FC_ROLL_DISP(0:DIMENSION_VTK)
   LOGICAL :: VTK_FC_TWIST_ANGLE(0:DIMENSION_VTK)

! Collision force
   LOGICAL :: VTK_PART_COLL_FORCE(0:DIMENSION_VTK)

! Coordination number
   LOGICAL :: VTK_PART_COORDINATION(0:DIMENSION_VTK)

! Overlap statistics (min, max, mean)
   LOGICAL :: VTK_PART_OVERLAP_STATS(0:DIMENSION_VTK)

! Index of glued partices
   LOGICAL :: VTK_PART_GLUEID(0:DIMENSION_VTK)

! Diameter of glued partices
   LOGICAL :: VTK_PART_GLUEDIAMETER(0:DIMENSION_VTK)

! Bounding Diameter of glued partices
   LOGICAL :: VTK_PART_GLUEBOUNDING(0:DIMENSION_VTK)

! Mass of glued partices
   LOGICAL :: VTK_PART_GLUEMASS(0:DIMENSION_VTK)

! Actual volume of glued partices
   LOGICAL :: VTK_PART_GLUEVOLUME(0:DIMENSION_VTK)

! Quaternion and rotation of glued partices
   LOGICAL :: VTK_PART_GLUEQUAT(0:DIMENSION_VTK)

! STL Group
   LOGICAL :: VTK_STL_GROUP_ID(0:DIMENSION_VTK)

! STL Collision count
   LOGICAL :: VTK_STL_COLL_COUNT(0:DIMENSION_VTK)

! STL heat transfer
   LOGICAL :: VTK_STL_Q(0:DIMENSION_VTK)

! STL surface area
   LOGICAL :: VTK_STL_AREA(0:DIMENSION_VTK)

! STL wear
   LOGICAL :: VTK_STL_WEAR(0:DIMENSION_VTK)

! STL nornal force
   LOGICAL :: VTK_STL_NORMAL_FORCE(0:DIMENSION_VTK)

! STL tangential force
   LOGICAL :: VTK_STL_TANG_FORCE(0:DIMENSION_VTK)

! STL force
   LOGICAL :: VTK_STL_FORCE(0:DIMENSION_VTK)

! STL torque
   LOGICAL :: VTK_STL_TORQUE(0:DIMENSION_VTK)

! STL normal velocity
   LOGICAL :: VTK_STL_NORMAL_VEL(0:DIMENSION_VTK)

! STL tangential velocity
   LOGICAL :: VTK_STL_TANG_VEL(0:DIMENSION_VTK)

! STL velocity
   LOGICAL :: VTK_STL_VEL(0:DIMENSION_VTK)


CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: STRIP                                              C
!  Purpose: Remove trailing *and* leading whitespace from a string     C
!                                                                      C
!  Author: Mark Meredith                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   FUNCTION STRIP(str) result(stripped)
      IMPLICIT NONE

      CHARACTER(:), allocatable :: stripped
      CHARACTER(LEN=*) :: str

      stripped = trim(adjustl(str))
   END FUNCTION STRIP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: ADD_VTU_DIRECTORY                                  C
!  Purpose: Prepend VTU_DIR to a filename (if VTU_DIR exists)          C
!                                                                      C
!  Author: Mark Meredith                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   FUNCTION ADD_VTU_DIRECTORY(str) result(path)
      IMPLICIT NONE
      CHARACTER(LEN=*), intent(in) :: str
      CHARACTER(:), allocatable  :: path
      CHARACTER(256) :: tmp

      path = str
      IF (TRIM(VTU_DIR) /= '.' .AND. TRIM(VTU_DIR) /= '') THEN
         tmp = './'//strip(VTU_DIR)//'/'//str
         path = strip(tmp)
      ENDIF
   END FUNCTION ADD_VTU_DIRECTORY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: FIND_VTK_BASENAME                                  C
!  Purpose: Return the filename without extension for .vtu .vtp .pvtp  C
!           .pvtu files. If second argument is present, append for     C
!           distributed I/O.                                           C
!  Author: Mark Meredith                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   FUNCTION FIND_VTK_BASENAME(MODE, pe, add_vtu_dir) RESULT(FILENAME)
      IMPLICIT NONE

      CHARACTER(:), allocatable :: FILENAME
      INTEGER, intent(in) :: MODE
      INTEGER, intent(in), optional :: pe
      LOGICAL, intent(in), optional :: add_vtu_dir

      INTEGER :: suffix
      CHARACTER(256) :: tmp, frame_c, suffix_c

      tmp = strip(VTK_FILEBASE(VTK_REGION))

      IF (TIME_DEPENDENT_FILENAME .AND. MODE == 0) THEN
         IF (FRAME(VTK_REGION) < 10000) THEN
            WRITE (frame_c, "(I4.4)") FRAME(VTK_REGION)
         ELSE
            WRITE (frame_c, *) FRAME(VTK_REGION)
         ENDIF
         tmp = trim(tmp) // "_" // strip(frame_c)
      END IF


! For distributed I/O, define the file name for each processor that owns particles
      IF (BDIST_IO) THEN
         IF (PRESENT(PE)) THEN
            suffix = pe
            write (suffix_c, "(i4.4)") suffix
            tmp = trim(tmp) // "_" // strip(suffix_c)
         ENDIF

      END IF

! Add the VTU directory path if necessary
! When using distributed I/O along with a VTU_DIR, the pvtu and pvtp files
! are written inside the VTU_DIR, and the VTU_DIR should not be added
! to the vtu and vtp files.
      IF (PRESENT(add_vtu_dir)) THEN
         IF(add_vtu_dir) THEN
            FILENAME = ADD_VTU_DIRECTORY(trim(tmp))
         ELSE
            FILENAME = trim(tmp)
         ENDIF
      ELSE
         FILENAME = ADD_VTU_DIRECTORY(trim(tmp))
      ENDIF

   END FUNCTION FIND_VTK_BASENAME


! Save fluid reaction rates at each cell
   SUBROUTINE SAVE_FLUID_RRATES_CELL_DATA(IJK,RATES)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IJK
   DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(IN) :: RATES
   INTEGER :: L

   IF(SAVE_FLUID_RRATES) THEN
      DO L=1, NO_OF_RXNS
         Fluid_RRates_out(IJK,L) = RATES(L)
      ENDDO
   ENDIF

   RETURN
   END SUBROUTINE SAVE_FLUID_RRATES_CELL_DATA


! Save particle/parcel reaction rates at each cell and each particle/parcel
   SUBROUTINE SAVE_DES_RRATES_DATA(NP,IJK,DES_RATES)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NP  ! Global index of particle
   INTEGER, INTENT(IN) :: IJK ! Fluid cell index containing NP
   DOUBLE PRECISION, INTENT(IN) :: DES_RATES(NO_OF_DES_RXNS)
   INTEGER :: L

      IF(SAVE_DES_RRATES) THEN
         DO L=1, NO_OF_DES_RXNS
            Des_RRates_out(IJK,L) = Des_RRates_out(IJK,L) + DES_RATES(L)
         ENDDO
      ENDIF

      IF(SAVE_PART_RRATES) THEN
         DO L=1, NO_OF_DES_RXNS
            Part_RRates_out(NP,L) = DES_RATES(L)
         ENDDO
      ENDIF
      RETURN

      END SUBROUTINE SAVE_DES_RRATES_DATA

END MODULE vtk
