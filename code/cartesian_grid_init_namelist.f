MODULE CG_INIT_NAMELIST
   CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C
!     Module name: CARTESIAN_GRID_INIT_NAMELIST                           C
!     Purpose: initialize the cartesian_grid-namelist                     C
!                                                                         C
!                                                                         C
!     Author: Jeff Dietiker                              Date: 26-Aug-08  C
!     Reviewer:                                          Date:            C
!     Comments:                                                           C
!                                                                         C
!                                                                         C
!  Keyword Documentation Format:                                          C
!<keyword category="category name" required="true/false"                  C
!                                    legacy="true/false">                 C
!  <description></description>                                            C
!  <arg index="" id="" max="" min=""/>                                    C
!  <dependent keyword="" value="DEFINED"/>                                C
!  <conflict keyword="" value="DEFINED"/>                                 C
!  <valid value="" note="" alias=""/>                                     C
!  <range min="" max="" />                                                C
!  MFIX_KEYWORD=INIT_VALUE                                                C
!</keyword>                                                               C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CARTESIAN_GRID_INIT_NAMELIST

      USE param1
      USE quadric
      USE cutcell
      USE polygon
      USE vtk
      USE progress_bar
      USE dashboard
      Use stl

      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!-----------------------------------------------
!
!
      INCLUDE 'cartesian_grid_namelist.inc'

!<keyword dtype="LOGICAL" category="cartesian grid" required="false" locked="true">
!  <description>Activate Cartesian grid cut cell technique.</description>
!  <dependent keyword="COORDINATES" value="CARTESIAN"/>
!  <conflict keyword="COORDINATES" value="CYLINDRICAL"/>
!  <valid value=".FALSE." note="Do not use Cartesian grid cut cell technique."/>
!  <valid value=".TRUE." note="Use Cartesian grid cut cell technique."/>
      CARTESIAN_GRID = .FALSE.
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Number of quadric surfaces defining the boundaries (&lt;=100).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <range min="0" max="100" />
      N_QUADRIC = 0
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>Use STL file to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value=".FALSE." note="Do not use STL file."/>
!  <valid value=".TRUE." note="Read triangulated geometry (for 3d geometry only) from geometry_####.stl."/>
      USE_STL = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>Use .msh file to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value=".FALSE." note="Do not use .msh file."/>
!  <valid value=".TRUE." note="Read geometry (for 3d geometry only) from geometry.msh."/>
      USE_MSH = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>Use polygons to describe geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value=".FALSE." note="Do not use polygons."/>
!  <valid value=".TRUE." note="Read polygon data (for 2d geometry only) from poly.dat."/>
      USE_POLYGON = .FALSE.
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Number of user-defined functions (currently limited to
!  0 or 1). If set to 1, the geometry is defined in the user
!  subroutine eval_usr_fct.f.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value="0" note="Do not use user-defined function" alias=""/>
!  <valid value="1" note="Use one user-defined function" alias=""/>
!  <range min="0" max="1" />
      N_USR_DEF = 0
!</keyword>

!<keyword dtype="CHARACTER" category="Cartesian grid" required="false" locked="true">
!  <description>Form of the quadric surface equation.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value="normal" note="Use normal form, as defined in equation
!    (1). The LAMDBAs and D must be defined"/>
!  <valid value="plane" note="Plane. Needs to define N_X,N_Y,N_Z
!    (unit normal vector pointing away from fluid cells)."/>
!  <valid value="x_cyl_int" note="Cylinder aligned with x-axis,
!    internal flow. Needs to define RADIUS(QID)."/>
!  <valid value="x_cyl_ext" note="Cylinder aligned with x-axis,
!    external flow. Needs to define RADIUS(QID)."/>
!  <valid value="y_cyl_int" note="Cylinder aligned with y-axis,
!    internal flow. Needs to define RADIUS(QID)."/>
!  <valid value="y_cyl_ext" note="Cylinder aligned with y-axis,
!    external flow. Needs to define RADIUS(QID)."/>
!  <valid value="z_cyl_int" note="Cylinder aligned with z-axis,
!    internal flow. Needs to define RADIUS(QID)."/>
!  <valid value="z_cyl_ext" note="Cylinder aligned with z-axis,
!    external flow. Needs to define RADIUS(QID)."/>
!  <valid value="x_cone" note="Cone aligned with x-axis, internal flow.
!    Needs to define HALF_ANGLE(QID)."/>
!  <valid value="y_cone" note="Cone aligned with y-axis, internal flow.
!    Needs to define HALF_ANGLE(QID)."/>
!  <valid value="z_cone" note="Cone aligned with z-axis, internal flow.
!    Needs to define HALF_ANGLE(QID)."/>
!  <valid value="sphere_int" note="Sphere, internal flow.
!    Needs to define RADIUS(QID)."/>
!  <valid value="sphere_ext" note="Sphere, external flow.
!    Needs to define RADIUS(QID)."/>
!  <valid value="C2C" note="Cylinder-to-cylinder conical junction,
!    internal flow. Needs to be defined between two cylinders."/>
!  <valid value="Torus_int" note="Torus, internal flow. Needs to
!    define TORUS_R1(QID) and TORUS_R2(QID).A torus is not a quadric
!    surface but is defined as a basic shape."/>
!  <valid value="Torus_ext" note="Torus, external flow. Needs to define
!    TORUS_R1(QID) and TORUS_R2(QID)."/>
!  <valid value="Y_UCOIL_EXT" note="Pair of parallel cylinders (y-direction),
!     capped at both ends by half a torus, to create a U-shaped coil. Needs
!     UCOIL_R1, UCOIL_R2, UCOIL_Y1, UCOIL_Y2."/>
!  <valid value="Y_UCOIL_EXT" note="Pair of parallel cylinders (y-direction),
!     capped at both ends by a cylinder at 90 degree angle
!     to create a U-shaped coil. Needs
!     UCOIL_R1, UCOIL_R2, UCOIL_Y1, UCOIL_Y2."/>
!  <valid value="XY_BEND_INT" note="Bend between
!     two cylinders in the XY plane, Needs
!     BEND_R1,BEND_R2,BEND_THETA1,BEND_THETA2."/>
!  <valid value="Y_C2C_INT" note="connects two vertical cylinders by a
!     conical section. Needs C2C_R1,C2C_R2,C2C_Y1,C2C_Y2."/>
!  <valid value="REACTOR1" note="Reactor, made of two vertical cylinders,
!     connected by a conical section. Each cylinder is rounded and closed
!     by a conical cap. Needs REACTOR1_R1,REACTOR1_R2,REACTOR1_Y1,REACTOR1_Y2,
!     REACTOR1_YR1,REACTOR1_YR2,REACTOR1_RR1,REACTOR1_RR2,
!     REACTOR1_THETA1,REACTOR1_THETA2."/>
      quadric_form(1:DIM_QUADRIC) = 'NORMAL'
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Scaling factor, applied to all quadric geometry parameters. Must
!    be a positive number.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <range min="0.0" max="" />
      quadric_scale = ONE
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Coefficient LAMBDA_X in equation (1) ('NORMAL' form) or
!    x-component of normal vector defining plane in equation (5)
!    ('DEGENERATE' form).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      lambda_x(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Coefficient LAMBDA_Y in equation (1) ('NORMAL' form) or
!    y-component of normal vector defining plane in equation (5)
!    ('DEGENERATE' form).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      lambda_y(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Coefficient LAMBDA_Z in equation (1) ('NORMAL' form) or
!    z-component of normal vector defining plane in equation (5)
!    ('DEGENERATE' form).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      lambda_z(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Coefficient D in equation (1).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      dquadric(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Rotation angle with respect to x-axis (degrees).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      theta_x(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Rotation angle with respect to y-axis (degrees).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      theta_y(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Rotation angle with respect to z-axis (degrees).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      theta_z(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Cylinder radius (used when QUADRIC_FORM = *_CYL_***)
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      Radius(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Cone half angle, expressed in degrees (used when
!    QUADRIC_FORM = *_CONE)
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      Half_angle(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Torus Radius 1 (used when QUADRIC_FORM = TORUS_*),
!    R1>R2 for a ring.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      Torus_R1(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Torus Radius 2 (used when QUADRIC_FORM = TORUS_*),
!    R1>R2 for a ring.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      Torus_R2(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    U-shaped coil Radius 1 (used when QUADRIC_FORM = UCOIL*),
!    UCOIL_R1>UCOIL_R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      UCOIL_R1(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    U-shaped coil Radius 2 (used when QUADRIC_FORM = UCOIL*),
!    UCOIL_R1>UCOIL_R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      UCOIL_R2(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    U-shaped coil ymax (used when QUADRIC_FORM = UCOIL*),
!    UCOIL_Y2>UCOIL_Y1.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      UCOIL_Y1(1:DIM_QUADRIC) = -UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    U-shaped coil ymin (used when QUADRIC_FORM = UCOIL*),
!    UCOIL_Y2>UCOIL_Y1.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      UCOIL_Y2(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Bend Radius 1 (used when QUADRIC_FORM = BEND*),
!    BEND_R1>BEND_R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BEND_R1(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Bend Radius 2 (used when QUADRIC_FORM = BEND*),
!    BEND_R1>BEND_R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BEND_R2(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Bend start angle, in degrees (used when QUADRIC_FORM = BEND*).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BEND_THETA1(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Bend end angle, in degrees (used when QUADRIC_FORM = BEND*).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BEND_THETA2(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Cylinder-cone_cylinder Radius 1 (used when QUADRIC_FORM = C2C*).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      C2C_R1(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Cylinder-cone_cylinder Radius 2 (used when QUADRIC_FORM = C2C*).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      C2C_R2(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Cylinder-cone_cylinder Y1 (used when QUADRIC_FORM = C2C*).
!    If Y1=Y2, then R1=R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      C2C_Y1(1:DIM_QUADRIC) = -UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Cylinder-cone_cylinder Y2 (used when QUADRIC_FORM = C2C*).
!    If Y1=Y2, then R1=R2.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      C2C_Y2(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Reactor 1, lower cylinder radius.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_R1(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Reactor 1, upper cylinder radius.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_R2(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Reactor 1, lower conical transition between cylinders.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_Y1(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Reactor 1, upper conical transition between cylinders.
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_Y2(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Reactor 1, lower rounding below cylinder.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_YR1(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Reactor 1, upper rounding above cylinder.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_YR2(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Reactor 1, lower rounding radius.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_RR1(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Reactor 1, upper rounding radius.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_RR2(1:DIM_QUADRIC) = UNDEFINED
!</keyword>


!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Reactor 1, lower rounding angle (degrees).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_THETA1(1:DIM_QUADRIC) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Reactor 1, upper rounding angle (degrees).</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REACTOR1_THETA2(1:DIM_QUADRIC) = UNDEFINED
!</keyword>


!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    X-component of normal vector defining the plane (used when
!    QUADRIC_FORM = PLANE).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      n_x(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Y-component of normal vector defining the plane (used when
!    QUADRIC_FORM = PLANE).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      n_y(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Z-component of normal vector defining the plane (used when
!    QUADRIC_FORM = PLANE).
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      n_z(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Translation in x-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      t_x(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Translation in y-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      t_y(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Translation in z-direction.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      t_z(1:DIM_QUADRIC) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Lower x-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_xmin(1:DIM_QUADRIC) = - LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Upper x-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_xmax(1:DIM_QUADRIC) =   LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Lower y-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_ymin(1:DIM_QUADRIC) = - LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Upper y-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_ymax(1:DIM_QUADRIC) =   LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Lower z-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_zmin(1:DIM_QUADRIC) = - LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Upper z-limit where the quadric is defined.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      clip_zmax(1:DIM_QUADRIC) =   LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Lower x-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_xmin(1:DIM_QUADRIC) = - LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Upper x-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_xmax(1:DIM_QUADRIC) =   LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Lower y-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_ymin(1:DIM_QUADRIC) = - LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Upper y-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_ymax(1:DIM_QUADRIC) =   LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Lower z-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_zmin(1:DIM_QUADRIC) = - LARGE_NUMBER
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Upper z-limit where the quadric is defined in a piecewise group.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      piece_zmax(1:DIM_QUADRIC) =   LARGE_NUMBER
!</keyword>


!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Flag defining the type of cells that are outside of the zone
!    defined by [CLIP_XMIN; CLIP_XMAX], [CLIP_YMIN; CLIP_YMAX],
!    [CLIP_ZMIN; CLIP_ZMAX].
!  </description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <valid value=".FALSE." note="Remove cells from computational domain."/>
!  <valid value=".TRUE." note="Treat cells as fluid cells."/>
      FLUID_IN_CLIPPED_REGION(1:DIM_QUADRIC) = .TRUE.
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Boundary condition flag.</description>
!  <arg index="1" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BC_ID_Q(1:DIM_QUADRIC) = UNDEFINED_I
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Number of group(s) of quadrics (&lt;=50).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      N_GROUP = 1
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Number of quadrics in the group.</description>
!  <arg index="1" id="Group ID" min="1" max="DIM_GROUP"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      GROUP_SIZE(1:DIM_GROUP) = 0
!</keyword>
      GROUP_SIZE(1) = 1 ! Default is one quadric in group 1

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Quadric ID assigned to a group.</description>
!  <arg index="1" id="Group ID" min="1" max="DIM_GROUP"/>
!  <arg index="2" id="Quadric ID" min="1" max="DIM_QUADRIC"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      GROUP_Q(DIM_GROUP,DIM_QUADRIC) = 0
!</keyword>

      GROUP_Q(1,1) = 1 ! Default is 1st quadric of group 1 is quadric 1

!<keyword dtype="CHARACTER" category="Cartesian grid" required="false" locked="true">
!  <description>Relation among quadrics of a same group.</description>
!  <valid value="OR" note="A point belongs to the computational domain
!    if at least one of f(x,y,z) among all quadrics is negative."/>
!  <valid value="AND" note="A point belongs to the computational domain
!   if all of f(x,y,z) among all quadrics are negative."/>
!  <valid value="PIECEWISE" note="When quadrics intersect along planes
!   that are perpendicular to either the x, y, or z-axis, quadrics can
!   be smoothly combined in a piecewise manner."/>
!  <arg index="1" id="Group ID" min="1" max="DIM_GROUP"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      GROUP_RELATION(1:DIM_GROUP) = 'OR'
!</keyword>

!<keyword dtype="CHARACTER" category="Cartesian grid" required="false" locked="true">

!  <description>
!    Relation between current group and combination of all previous groups.
!  </description>
!  <valid value="OR" note="A point belongs to the computational domain
!    if f-value for the current group or f-value for the combination of
!    previous groups is negative."/>
!  <valid value="AND" note="A point belongs to the computational
!    domain if f-value for the current group and f-value for the
!    combination of previous groups is negative."/>
!  <arg index="1" id="Group ID" min="1" max="DIM_GROUP"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      RELATION_WITH_PREVIOUS(1:DIM_GROUP) = 'OR'
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Tolerance used to snap an intersection point onto an
!  existing cell corner (expressed as a fraction of edge length,
!  between 0.0 and 0.5). For stretched grids, three values can be
!  entered in the x, y and z directions.</description>
!  <range min="0.0" max="0.5" />
!  <arg index="1" id="Direction" min="1" max="3"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_SNAP(1)    = 0.00D0  ! 0% of original edge length
!</keyword>
      TOL_SNAP(2)    = UNDEFINED
      TOL_SNAP(3)    = UNDEFINED

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Tolerance used to limit acceptable values of normal
!  distance to the wall (expressed as a fraction of cell diagonal,
!  between 0.0 and 1.0).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_DELH       = 0.00D0  ! 0% of original Diagonal
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Tolerance used to detect small cells (expressed as a
!  fraction of cell volume, between 0.0 and 1.0).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_SMALL_CELL = 0.01D0  ! 1% of original cell volume
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>Option to distribute small cut cells among their neighbors.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      DISTRIBUTE_SMALL_CUT_CELLS = .FALSE.
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Tolerance used to remove duplicate nodes (expressed as
!  a fraction of cell diagonal, between 0.0 and 1.0).</description>
      TOL_MERGE      = 1.0D-12 ! fraction of original cell diagonal
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Tolerance used to detect small faces (expressed as a
!  fraction of original face area, between 0.0 and 1.0).</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_SMALL_AREA = 0.01D0  ! 1% of original face area
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Maximum acceptable value of interpolation correction factor.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      ALPHA_MAX      = ONE
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Tolerance used to find intersection of quadric surfaces or user-defined function with background grid.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_F     = 1.0D-9
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Tolerance used to find intersection of polygon with background grid.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_POLY  = 1.0D-9
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Maximum number of iterations used to find intersection points.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      ITERMAX_INT = 10000
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Tolerance used to find intersection of STL triangles with background grid.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_STL = 1.0D-6        ! Settings for STL file
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Smallest angle accepted for valid STL triangles (in
!  degrees). Triangles having an angle smaller that this value will be
!  ignored.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      STL_SMALL_ANGLE = 0.1  ! Degrees
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Angle (degrees) to detect if two STL facets are nearly aligned.
!  This is used to ignore edge collision in DEM simulations.
!  Particle-edge collision will be ignored if the angle between two facets is smaller than STL_NB_ANGLE.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      STL_NB_ANGLE = 5.0  ! Degrees
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Dot product tolerance when determining if a point lies in a facet.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_STL_DP = 1.0D-3        ! Settings for STL file
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Maximum number of STL facets per cell.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      DIM_FACETS_PER_CELL  = 10        !
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Defines value of F_STL outside of the STL geometry. A
!  value of 1.0 means the domain outside of the STL geometry is
!  excluded from computation, i.e., an internal flow is
!  computed. Note: This depends on the direction of the facet normal vectors.
!  </description>
!  <valid value="-1.0" note="model an external flow"/>
!  <valid value="1.0" note="model an internal flow"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      OUT_STL_VALUE = 1.0
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>Option to flip STL facet normals.
!               The index corresponds to the BC ID the STL
!               file is applied to.
!  </description>
!  <valid value=".True." note="Flip normals"/>
!  <valid value=".False." note="Do not flip normals"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="BC" min="1" max="DIMENSION_BC"/>
      FLIP_STL_NORMALS(1:DIMENSION_BC) = .False.
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Boundary condition flag for the STL geometry</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      STL_BC_ID = UNDEFINED_I
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Translation in x-direction, applied to the STL geometry.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TX_STL = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Translation in y-direction, applied to the STL geometry.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TY_STL = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Translation in z-direction, applied to the STL geometry.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TZ_STL = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Scaling factor, applied to the STL geometry. Note that translation
!    occurs after scaling.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      SCALE_STL = ONE
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Tolerance used to find intersection of .msh file with background grid.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TOL_MSH = 1.0D-6        ! Settings for MSH file
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Defines value of f outside of the .msh geometry. a
!  value of 1.0 means the domain outside of the .msh geometry is
!  excluded from computation, i.e., an internal flow is
!  computed.</description>
!  <valid value="-1.0" note="model an external flow"/>
!  <valid value="1.0" note="model an internal flow"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      OUT_MSH_VALUE = 1.0
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Translation in x-direction, applied to the .msh geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TX_MSH = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Translation in y-direction, applied to the .msh geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TY_MSH = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Translation in z-direction, applied to the .msh geometry.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TZ_MSH = ZERO
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>Scaling factor, applied to the .msh geometry. Note that translation occurs after scaling.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      SCALE_MSH = ONE
!</keyword>

!<keyword dtype="CHARACTER" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Ray propagation order used to determine whether any point is
!    located inside or outside of the STL surface.
!  </description>
!  <valid value="   " note="Propagation occurs in unique cell index order,
!    from IJK=IJKSTART3 to IJKEND3, one neighbor at a time (West, East,
!    South, North, Bottom, Top)"/>
!  <valid value="ijk" note="Propagation occurs in the
!    I, followed by J, and K directions"/>
!  <valid value="jki" note="Propagation occurs in the
!    J, followed by K, and I directions"/>
!  <valid value="kij" note="Propagation occurs in the
!    K, followed by I, and J directions"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      CAD_PROPAGATE_ORDER = '   '
!</keyword>

!<keyword dtype="CHARACTER" category="Cartesian grid" required="false" locked="true">
!  <description>Ray direction when propagating CAD value</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      RAY_DIR = 'X-'
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Flag to detect and treat corner cells the same way as
!    in the original MFiX version (i.e. without cut cells).
!  </description>
!  <valid value=".TRUE." note="Some cut cells may be treated as corner cells."/>
!  <valid value=".FALSE." note="Do not treat cut cells as corner cells."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      SET_CORNER_CELLS = .FALSE.
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Factor used to allocate cut cell arrays (expressed as a fraction
!    of DIMENSION_3G).
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      FAC_DIM_MAX_CUT_CELL = 1.0
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="true">
!  <description>Write VTK files at regular intervals.</description>
!  <valid value=".FALSE." note="Do not write VTK files. If there are
!  cut cells, they will not be displayed from the usual .res file"/>
!  <valid value=".TRUE." note="Valid only if Cartesian_grid = .TRUE."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      WRITE_VTK_FILES = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="true">
!  <description>Use time-dependent VTK file names</description>
!  <valid value=".FALSE." note="The VTK file overwrites the previous
!    file (recommended for steady-state computation)."/>
!  <valid value=".TRUE." note="A sequential integer is appended to the
!    VTK filenames as they are written to create a series of files
!    (recommended for transient computation)."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      TIME_DEPENDENT_FILENAME = .TRUE.
!</keyword>

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>
!    Time interval (expressed in seconds of simulation time) at which VTK
!    files are written.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_DT(0:DIMENSION_VTK)     = UNDEFINED
!</keyword>

!<keyword dtype="INTEGER" category="Output Control" required="false" locked="false">
!  <description>
!    Type of GSP output shape in VTK files.
!  </description>
!  <valid value="1" note="Use actual glued particle shape (show all component spheres)."/>
!  <valid value="2" note="Use bounding sphere to represent entire GSP."/>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_GSP_OUTPUT_TYPE(0:DIMENSION_VTK) = 1
!</keyword>

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>
!    Starting time (expressed in seconds of simulation time) at which VTK
!    files are written.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_T_START(0:DIMENSION_VTK)     = -UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>
!    Ending time (expressed in seconds of simulation time) at which VTK
!    files are written.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_T_END(0:DIMENSION_VTK)     = UNDEFINED
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>When VTK_DBG_FILE is .TRUE., the VTK region file is only written
!   when the subroutine WRITE_DBG_VTU_AND_VTP_FILES is called, typically in a UDF.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTK_DBG_FILE(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="INTEGER" category="Output Control" required="false" locked="false">
!  <description>List of variables written in the VTK files.</description>
!  <valid value="1" note="Void fraction (EP_g)"/>
!  <valid value="2" note="Gas pressure, solids pressure (P_g, P_star)"/>
!  <valid value="3" note="Gas velocity (U_g, V_g, W_g)"/>
!  <valid value="4" note="Solids velocity (U_s, V_s, W_s)"/>
!  <valid value="5" note="Solids density (ROP_s)"/>
!  <valid value="6" note="Gas and solids temperature (T_g, T_s)"/>
!  <valid value="7" note="Gas and solids mass fractions (X_g, X_s)"/>
!  <valid value="8" note="Granular temperature (Theta_m)"/>
!  <valid value="9" note="Scalar"/>
!  <valid value="10" note="Reaction rates"/>
!  <valid value="11" note="K and Epsilon"/>
!  <valid value="12" note="Vorticity magnitude and lambda_2"/>
!  <valid value="100" note="Grid Partition"/>
!  <valid value="101" note="Boundary Condition ID"/>
!  <valid value="102" note="Distance to wall"/>
!  <valid value="103" note="DEM facet count"/>
!  <valid value="104" note="DEM Neighboring facets"/>
!  <valid value="999" note="Cell IJK index"/>
!  <valid value="1000" note="Cut face normal vector"/>
!  <arg index="1" id="IDX" min="1" max="20"/>
      VTK_VAR(1) = UNDEFINED_I
!</keyword>
      VTK_VAR(2:20) = UNDEFINED_I

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>X coordinate of the VTK region west face.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_X_W(0:DIMENSION_VTK) = -UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>X coordinate of the VTK region east face.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_X_E(0:DIMENSION_VTK) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>Y coordinate of the VTK region south face.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_Y_S(0:DIMENSION_VTK) = -UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>Y coordinate of the VTK region north face.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_Y_N(0:DIMENSION_VTK) = UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>Z coordinate of the VTK region bottom face.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_Z_B(0:DIMENSION_VTK) = -UNDEFINED
!</keyword>

!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>Z coordinate of the VTK region top face.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_Z_T(0:DIMENSION_VTK) = UNDEFINED
!</keyword>

!<keyword dtype="CHARACTER" category="Output Control" required="false" locked="false">
!  <description>VTK region output file name base.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_FILEBASE(0:DIMENSION_VTK) = UNDEFINED_C
!</keyword>

!<keyword dtype="CHARACTER" category="Output Control" required="false" locked="false">
!  <description>Type of data to write in the VTK file.</description>
!  <valid value="C" note="Cell data (VTU file)."/>
!  <valid value="P" note="Particle data (VTP file)."/>
!  <valid value="G" note="Geometry data (STL file)."/>
!  <valid value="F" note="Force chain data (VTP file)."/>
!  <valid value="N" note="Glued (non-spherical) particle data (VTP file)."/>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_DATA(0:DIMENSION_VTK) = 'C'
!</keyword>

!<keyword dtype="CHARACTER" category="Output Control" required="false" locked="false">
!  <description>Geometry source for VTK geometry output.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="GEO" min="1" max="VTK_GEO_MAX"/>
      VTK_GEO(0:DIMENSION_VTK,1:VTK_GEO_MAX) = ''
!</keyword>

!<keyword dtype="CHARACTER" category="Output Control" required="false" locked="false">
!  <description>Particle selection mode in a VTK region.</description>
!  <valid value="C" note="Select particles with centers inside
!    the VTK region."/>
!  <valid value="P" note="Select particles that are entirely inside
!    the VTK region."/>
!  <valid value="I" note="Select particles that are inside or
!    intersect the edges of the VTK region."/>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_SELECT_MODE(0:DIMENSION_VTK) = 'C'
!</keyword>

!<keyword dtype="INTEGER" category="Output Control" required="false" locked="false">
!  <description>
!    Specifies the number of subdivisions in the x-axial direction
!    to decompose a VTK region. Leave undefined to export the full
!    region. (Slice a volume into planes; cut a plane into lines;
!    or break a line into points.)
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_NXS(0:DIMENSION_VTK) = 0
!</keyword>

!<keyword dtype="INTEGER" category="Output Control" required="false" locked="false">
!  <description>
!    Specifies the number of subdivisions in the y-axial direction
!    to decompose a VTK region. Leave undefined to export the full
!    region. (Slice a volume into planes; cut a plane into lines;
!    or break a line into points.)
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_NYS(0:DIMENSION_VTK) = 0
!</keyword>

!<keyword dtype="INTEGER" category="Output Control" required="false" locked="false">
!  <description>
!    Specifies the number of subdivisions in the z-axial direction
!    to decompose a VTK region. Leave undefined to export the full
!    region. (Slice a volume into planes; cut a plane into lines;
!    or break a line into points.)
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_NZS(0:DIMENSION_VTK) = 0
!</keyword>


!<keyword dtype="REAL" category="Output Control" required="false" locked="false">
!  <description>Tolerance to detect particles in a VTK region's slice.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_SLICE_TOL(0:DIMENSION_VTK) = ZERO
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write the void fraction in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_EP_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write the gas density in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_RO_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write the gas pressure in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_P_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write the gas velocity vector in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_VEL_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>
!    Write x-component of gas velocity vector in region.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_U_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>
!    Write y-component of gas velocity vector in region.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_V_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>
!    Write z-component of gas velocity vector in region.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_W_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write the solids pressure that prevents overpacking in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_P_star(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write the solids pressure as a result of granular motion in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_P_s(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write solids velocity vector in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_VEL_S(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>
!     Write x-component of solids velocity vector in region.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_U_s(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>
!    Write y-component of solids velocity vector in region.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_V_s(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>
!    Write z-component of solids velocity vector in region.
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_W_s(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write solids bulk density in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_ROP_s(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write solids density in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_RO_s(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write solids volume fraction in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_EP_s(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write gas temperature in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_T_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write solids temperature in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_T_s(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write gas specific heat in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_C_Pg(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write solids specific heat in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_C_Ps(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write gas phase viscosity in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_Mu_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write gas phase thermal conductivity in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_K_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>


!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write gas phase mixture molecular weight in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_MW_MIX_g(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write gas phase species mass fraction in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Gas Species" min="1" max="DIM_N_g"/>
      VTK_X_g(0:DIMENSION_VTK,1:DIM_N_g) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write gas phase species molar fraction in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Gas Species" min="1" max="DIM_N_g"/>
      VTK_Y_g(0:DIMENSION_VTK,1:DIM_N_g) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write solids phase species mass fraction in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
!  <arg index="3" id="Solids Species" min="1" max="DIM_N_s"/>
      VTK_X_s(0:DIMENSION_VTK,1:DIM_M,1:DIM_N_s) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write granular temperature in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_Theta_m(0:DIMENSION_VTK,1:DIM_M) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write scalar in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Scalar" min="1" max="DIM_scalar"/>
      VTK_Scalar(0:DIMENSION_VTK,1:DIM_scalar) =.FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write reaction rate in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="rate" min="1" max="nRRmax"/>
      VTK_RRate(0:DIMENSION_VTK,1:nRRmax) = .FALSE.
!</keyword>

!<keyword dtype="CHARACTER" category="Output Control" required="false" locked="false">
!  <description>Reaction rate custom label in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="rate" min="1" max="nRRmax"/>
      VTK_RRate_label(0:DIMENSION_VTK,1:nRRmax) = ''
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write fluid reaction rates in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="rate" min="1" max="nRRmax"/>
      VTK_FLUID_RRate(0:DIMENSION_VTK,1:nRRmax) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write Lagrangian (DEM, CGDEM, PIC) reaction rates in region. Particle/parcel data is
!  accumulated in fluid cells.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="rate" min="1" max="nRRmax"/>
      VTK_DES_RRate(0:DIMENSION_VTK,1:nRRmax) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write turbulent kinetic energy in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_K_Turb_G(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write turbulent dissipation rate in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_E_Turb_G(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write vorticity magnitude in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_VORTICITY(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write lambda_2 in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_LAMBDA_2(0:DIMENSION_VTK)  = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write void grid partition in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PARTITION(0:DIMENSION_VTK) = .FALSE.
!</keyword>

      VTK_DOMAIN_DECOMPOSITION(0:DIMENSION_VTK) = .FALSE.

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write boundary condition ID in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_BC_ID(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write wall distance in region.</description>
      VTK_DWALL(0:DIMENSION_VTK) = .FALSE.
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write STL facet count for DES in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_FACET_COUNT_DES(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write neighboring facets in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_NB_FACET_DES(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write cell IJK index in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_IJK(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write scalar cell center coordinates in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_SCALAR_CC(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write cut face normal vector in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_NORMAL(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write debug variable in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="IDX" min="1" max="15"/>
      VTK_DEBUG(0:DIMENSION_VTK,15) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle diameter in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_DIAMETER(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write Coarse Grain physical particle diameter in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_PHYSICAL_DIAMETER(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write Coarse Grain particle statistical weight in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_CGP_STAT_WT(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle velocity in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_VEL(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle angular velocity in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_ANGULAR_VEL(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle orientation in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_ORIENTATION(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle user-defined variable in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="USR_VAR" min="1" max="VTK_PART_USRmax"/>
      VTK_PART_USR_VAR(0:DIMENSION_VTK,1:VTK_PART_USRmax)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle temperature in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_TEMP(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle specific heat in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_C_Ps(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle species mass fraction in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Des_species" min="1" max="100"/>
      VTK_PART_X_s(0:DIMENSION_VTK,100)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle/parcel (DEM, CGDEM, PIC) reaction rates in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="rate" min="1" max="nRRmax"/>
      VTK_PART_RRate(0:DIMENSION_VTK,1:nRRmax) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle density in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_DENSITY(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle mass in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_MASS(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle volume in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_VOLUME(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle drag coefficient (F_GP) in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_DRAG_COEFF(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle cohesion in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_COHESION(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle rank in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_RANK(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle ID in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_ID(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle residence time in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_RESIDENCE_TIME(0:DIMENSION_VTK) = .FALSE.
!</keyword>
!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle phase ID in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_PHASE_ID(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Option to save or not save particles belonging to solids phases in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <arg index="2" id="Phase" min="1" max="DIM_M"/>
      VTK_PART_PHASE(0:DIMENSION_VTK,1:DIM_M) = .TRUE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Option to only save the maximum force chain data per particle.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_FCMAX_ONLY(0:DIMENSION_VTK) = .FALSE.

!</keyword>
!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write JKR contact state in force-chain data.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_FC_JKR_STATE(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write JKR contact radius in force-chain data.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_FC_JKR_RADIUS(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write JKR effective surface energy in force-chain data.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_FC_JKR_GAMMA(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write tangential displacement in force-chain data.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_FC_TANG_DISP(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write rolling displacement in force-chain data.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_FC_ROLL_DISP(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write signed JKR elastic twist angle in radians in particle-particle force-chain data; after torque-limit return mapping, not cumulative rotation.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
!  <dependent keyword="VTK_DATA" value="F"/>
      VTK_FC_TWIST_ANGLE(0:DIMENSION_VTK) = .FALSE.
!</keyword>


!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle collision force in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_COLL_FORCE(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle coordination number in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_COORDINATION(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write particle overlap statistics (min, max, mean) in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_OVERLAP_STATS(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write cut-cell data only in region.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_CUTCELL_ONLY(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write Index of glued particles of each sphere in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_GLUEID(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write diameter of glued particles in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_GLUEDIAMETER(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write bounding diameter of glued particles in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_GLUEBOUNDING(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write mass of glued particles in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_GLUEMASS(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write actual volume of glued particles in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_GLUEVOLUME(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write rotation and quaternion of glued particles in region.</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_PART_GLUEQUAT(0:DIMENSION_VTK)= .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write STL facet group ID in region (geometry data).</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_GROUP_ID(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write collision count between particles and STL facet in region (geometry data).</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_COLL_COUNT(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write heat transfer between particles and STL facet in region (geometry data).</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_Q(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write STL facet surface area in region (geometry data).</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_AREA(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write a measure of the wear on the STL facet in region (geometry data).</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_WEAR(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write normal force exerted by particles onto the STL facet in region (geometry data).</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_NORMAL_FORCE(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write tangential force exerted by particles onto the STL facet in region (geometry data).</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_TANG_FORCE(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write total force exerted by particles onto the STL facet in region (geometry data).</description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_FORCE(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write torque exerted by particles onto the STL facet in region (geometry data).
!     The torque pivot point must be defined for each geometry (STL BC, regular BC or default BC).
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_TORQUE(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write normal velocity of a moving STL facet in region (geometry data).
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_NORMAL_VEL(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write tangential velocity of a moving STL facet in region (geometry data).
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_TANG_VEL(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Output Control" required="false" locked="false">
!  <description>Write velocity of a moving STL facet in region (geometry data).
!  </description>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      VTK_STL_VEL(0:DIMENSION_VTK) = .FALSE.
!</keyword>

!<keyword dtype="INTEGER" category="Output Control" required="false" locked="false">
!  <description>Starting Index appended to VTU regions.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
!  <arg index="1" id="VTK" min="1" max="DIMENSION_VTK"/>
      FRAME(0:DIMENSION_VTK) = -1
!</keyword>

!<keyword dtype="CHARACTER" category="Output Control" required="false" locked="false">
!  <description>Directory where VTK files are stored. By default,
!   VTK files are written in the run directory.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      VTU_DIR = ''
!</keyword>


!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Option for pressure gradient computation in cut cells.
!  </description>
!  <valid value="1" note="Use maximum of (east/west), (north/south),
!    and (top/bottom) pairs of velocity cells."/>
!  <valid value="2" note="Use both (east/west), (north/south), and
!    (top/bottom) areas of velocity cells."/>
!  <valid value="0" note="Use east, north and top areas of pressure cell (same as standard cells)."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      PG_OPTION = 0
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>Run code in safe mode.</description>
!  <valid value="1" note="Performs initial preprocessing but use all
!    original MFiX subroutines during flow solution (using only cell
!    volumes and areas of cut cells)."/>
!  <valid value="0" note="Runs the code with modified subroutines for
!    cut cell treatment."/>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      CG_SAFE_MODE = 0
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Prints any warning message encountered during preprocessing
!    on the screen.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      PRINT_WARNINGS = .FALSE.
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Underrelaxation factor used in cut cells (only CG_UR_FAC(2) is used).
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      CG_UR_FAC = 1.0
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Print a progress bar during each major step of preprocessing stage.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      PRINT_PROGRESS_BAR = .FALSE.
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Width of the progress bar (complete status), expressed in number
!    of characters (between 10 and 80).
!  </description>
!  <range min="10" max="80" />
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BAR_WIDTH = 50
!</keyword>

!<keyword dtype="CHARACTER" category="Cartesian grid" required="false" locked="true">
!  <description>Character used to create the progress bar.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BAR_CHAR  = '='
!</keyword>

!<keyword dtype="REAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Update frequency of progress bar, expressed in percent of total
!    length (between 1.0 and 100.0).
!  </description>
!  <range min="1.0" max="100.0" />
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      BAR_RESOLUTION = 5.0
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Writes the file dashboard.txt at regular intervals. The file shows
!    a summary of the simulation progress.
!  </description>
      WRITE_DASHBOARD = .FALSE.
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Frequency, expressed in terms of iterations, at which the dashboard
!    is updated.
!  </description>
      F_DASHBOARD = 1
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>Location of control points in x-direction.</description>
!  <arg index="1" id="CTRL" min="0" max="MAX_CP" />
      CPX(0:MAX_CP) = -UNDEFINED
!</keyword>

!<keyword dtype="INTEGER" category="Geometry and Discretization" required="false" locked="true">
!  <description>Number of cells within a segment (x-direction).</description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      NCX(1:MAX_CP)      = 0
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Expansion ratio (last DX/first DX) in a segment (x-direction).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      ERX(1:MAX_CP)      = ONE
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>Value of first DX in a segment (x-direction). A negative
!      value will copy DX from previous segment (if available).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      FIRST_DX(1:MAX_CP) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>Value of last DX in a segment (x-direction). A
!  negative value will copy DX from next segment (if
!  available).</description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      LAST_DX(1:MAX_CP)  = ZERO
!</keyword>


!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>Location of control points in y-direction.</description>
!  <arg index="1" id="CTRL" min="0" max="MAX_CP" />
      CPY(0:MAX_CP) = -UNDEFINED
!</keyword>

!<keyword dtype="INTEGER" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Number of cells within a segment (y-direction).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      NCY(1:MAX_CP) = 0
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Expansion ratio (last DY/first DY) in a segment (y-direction).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      ERY(1:MAX_CP) = ONE
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Value of first DY in a segment (y-direction). A negative value will
!    copy DY from previous segment (if available).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      FIRST_DY(1:MAX_CP) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Value of last DY in a segment (y-direction). A negative value will
!    copy DY from next segment (if available).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      LAST_DY(1:MAX_CP)  = ZERO
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>Location of control points in z-direction.</description>
!  <arg index="1" id="CTRL" min="0" max="MAX_CP" />
      CPZ(0:MAX_CP) = -UNDEFINED
!</keyword>

!<keyword dtype="INTEGER" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Number of cells within a segment (z-direction).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      NCZ(1:MAX_CP) = 0
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Expansion ratio (last DZ/first DZ) in a segment (z-direction).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      ERZ(1:MAX_CP)      = ONE
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Value of first DZ in a segment (z-direction). A negative value will
!    copy DZ from previous segment (if available).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      FIRST_DZ(1:MAX_CP) = ZERO
!</keyword>

!<keyword dtype="REAL" category="Geometry and Discretization" required="false" locked="true">
!  <description>
!    Value of last DZ in a segment (z-direction). A negative value will
!    copy DZ from next segment (if available).
!  </description>
!  <arg index="1" id="CTRL" min="1" max="MAX_CP" />
      LAST_DZ(1:MAX_CP)  = ZERO
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Turns on the re-indexing of cells. When true, inactive (dead)
!    cells are removed from computational domain.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      RE_INDEXING = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Attempts to adjust grid partition. Each processor will be assigned
!    its own size to minimize load imbalance.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      ADJUST_PROC_DOMAIN_SIZE = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Attempts to adjust grid partition. Each processor will
!    be assigned its own size to minimize load imbalance.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      REPORT_BEST_DOMAIN_SIZE = .FALSE.
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Temporary setting used in serial run to report best domain size
!    for parallel run.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      NODESI_REPORT = 1
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Temporary setting used in serial run to report best domain size
!    for parallel run.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      NODESJ_REPORT = 1
!</keyword>

!<keyword dtype="INTEGER" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Temporary setting used in serial run to report best domain size
!    for parallel run.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      NODESK_REPORT = 1
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>
!    Attempts to minimize the size of the send/receive layers.
!  </description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      MINIMIZE_SEND_RECV = .TRUE.
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>Brute force calculation of wall distance.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      DWALL_BRUTE_FORCE = .FALSE.
!</keyword>

!<keyword dtype="LOGICAL" category="Cartesian grid" required="false" locked="true">
!  <description>Flag to turn off non-orthogonality terms in phi diffusion derivatives.</description>
!  <dependent keyword="CARTESIAN_GRID" value=".TRUE."/>
      NOC_PHI_DIF = .TRUE.
!</keyword>
      RETURN
      END SUBROUTINE CARTESIAN_GRID_INIT_NAMELIST
END MODULE CG_INIT_NAMELIST
