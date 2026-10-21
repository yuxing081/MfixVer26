!----------------------------------------------------------------------!
! Project-local JKR keyword metadata for MFIX-GUI.                     !
!                                                                      !
! These local defaults intentionally match DES_INIT_NAMELIST. The      !
! solver declares, initializes, and reads the real variables through   !
! the DES namelist; this metadata routine only lets the GUI        !
! discover them in a custom project.                                   !
!----------------------------------------------------------------------!
      SUBROUTINE USR_INIT_NAMELIST

      USE vtk, ONLY: DIMENSION_VTK
      IMPLICIT NONE

      INTEGER, PARAMETER :: GUI_MAX_PHASES = 10
      DOUBLE PRECISION, PARAMETER :: UNDEFINED = -9.87654321D30
      LOGICAL :: USE_JKR
      LOGICAL :: VTK_FC_JKR_STATE(0:DIMENSION_VTK)
      LOGICAL :: VTK_FC_JKR_RADIUS(0:DIMENSION_VTK)
      LOGICAL :: VTK_FC_JKR_GAMMA(0:DIMENSION_VTK)
      LOGICAL :: VTK_FC_TANG_DISP(0:DIMENSION_VTK)
      LOGICAL :: VTK_FC_ROLL_DISP(0:DIMENSION_VTK)
      LOGICAL :: VTK_FC_TWIST_ANGLE(0:DIMENSION_VTK)
      DOUBLE PRECISION :: DES_SURFACE_ENERGY(GUI_MAX_PHASES)
      DOUBLE PRECISION :: DES_SURFACE_ENERGY_WALL
      DOUBLE PRECISION :: JKR_THETAR_CRIT
      DOUBLE PRECISION :: WALL_JKR_THETAR_CRIT
      DOUBLE PRECISION :: JKR_ALPHA_SAFETY
      DOUBLE PRECISION :: JKR_ROLL_DAMP_COEFF

!<keyword dtype="LOGICAL" category="Discrete Element Model" required="false" locked="true">
!  <description>Enable spherical-particle JKR adhesion and Chen rolling.</description>
!  <conflict keyword="USE_COHESION" value=".TRUE."/>
      USE_JKR = .FALSE.
!</keyword>

!<keyword dtype="REAL" category="Discrete Element Model" required="false" locked="true">
!  <description>Positive per-phase particle surface energy used by JKR, in J/m^2.</description>
!  <dependent keyword="USE_JKR" value=".TRUE."/>
!  <arg index="1" id="Phase" min="1" max="MMAX"/>
!  <range min="+0"/>
      DES_SURFACE_ENERGY(:) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Discrete Element Model" required="false" locked="true">
!  <description>Positive wall surface energy used by JKR, in J/m^2.</description>
!  <dependent keyword="USE_JKR" value=".TRUE."/>
!  <range min="+0"/>
      DES_SURFACE_ENERGY_WALL = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Discrete Element Model" required="false" locked="true">
!  <description>Critical Chen rolling angle for particle-particle contacts.</description>
!  <dependent keyword="USE_JKR" value=".TRUE."/>
!  <range min="1.0e-4" max="1.0e-1"/>
      JKR_THETAR_CRIT = 1.0D-2
!</keyword>

!<keyword dtype="REAL" category="Discrete Element Model" required="false" locked="true">
!  <description>Critical Chen rolling angle for particle-wall contacts.</description>
!  <dependent keyword="USE_JKR" value=".TRUE."/>
!  <range min="1.0e-4" max="1.0e-1"/>
      WALL_JKR_THETAR_CRIT = 1.0D-2
!</keyword>

!<keyword dtype="REAL" category="Discrete Element Model" required="false" locked="true">
!  <description>Safety factor applied to the tensile JKR separation cutoff.</description>
!  <dependent keyword="USE_JKR" value=".TRUE."/>
!  <range min="0.7" max="1.0"/>
      JKR_ALPHA_SAFETY = 1.0D0
!</keyword>

!<keyword dtype="REAL" category="Discrete Element Model" required="false" locked="true">
!  <description>Rolling damping time mu_R in seconds; eta_R=mu_R*abs(F_ne). Zero preserves elastic rolling.</description>
!  <dependent keyword="USE_JKR" value=".TRUE."/>
!  <range min="0"/>
      JKR_ROLL_DAMP_COEFF = 0.0D0
!</keyword>


!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write JKR contact state in particle-particle force-chain data: 1 for an active JKR contact, 0 otherwise.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <dependent keyword="VTK_DATA" value="F"/>
      VTK_FC_JKR_STATE(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write simulated JKR contact radius in meters in particle-particle force-chain data.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <dependent keyword="VTK_DATA" value="F"/>
      VTK_FC_JKR_RADIUS(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write effective JKR surface energy in J/m^2 in particle-particle force-chain data.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <dependent keyword="VTK_DATA" value="F"/>
      VTK_FC_JKR_GAMMA(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write updated tangential displacement history in meters in particle-particle force-chain data.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <dependent keyword="VTK_DATA" value="F"/>
      VTK_FC_TANG_DISP(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write updated JKR rolling displacement history in meters in particle-particle force-chain data, after the rolling limit.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <dependent keyword="VTK_DATA" value="F"/>
      VTK_FC_ROLL_DISP(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write signed JKR elastic twist angle in radians in particle-particle force-chain data; after torque-limit return mapping, not cumulative rotation.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <dependent keyword="VTK_DATA" value="F"/>
      VTK_FC_TWIST_ANGLE(0:DIMENSION_VTK) = .FALSE.
!</keyword>

      RETURN
      END SUBROUTINE USR_INIT_NAMELIST
