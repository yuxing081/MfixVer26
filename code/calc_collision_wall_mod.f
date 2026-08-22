#include "error.inc"

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COLLISION_WALL                                    C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!                                                                      C
!  Purpose: subroutines for particle-wall collisions when cutcell is   C
!           used. Also contains rehack of routines for cfslide and     C
!           cffctow which might be different from the stand alone      C
!           routines. Eventually all the DEM routines will be          C
!           consolidated.                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
MODULE CALC_COLLISION_WALL

   USE bc
   USE compar, only: istart3, iend3, jstart3, jend3, kstart3, kend3
   USE compar, only: pe_io, mype, numPEs
   USE constant, only: pi
   USE cutcell, only: area_cut, blocked_cell_at, cartesian_grid
   USE des_thermo, only: flpc, calc_cond_des, des_t_s, q_source
   USE des_thermo_cond, only: des_conduction_wall, des_qw_cond
   USE discretelement
   USE error_manager
   USE exit, only: mfix_exit
   USE functions, only: funijk, fluid_at, is_normal
   USE geometry, only: do_k, imax1, imax2, imin1, imin2
   USE geometry, only: dx, dy, dz
   USE geometry, only: kmax1, kmax2, kmin1, kmin2, zlength, xlength, ylength
   USE geometry, only: no_k, jmax1, jmax2, jmin1, jmin2
   USE geometry, only: cyclic_x, cyclic_y, cyclic_z
   USE param, only: dimension_i, dimension_j, dimension_k
   USE param1, only: half, one, small_number, undefined
   USE particles_in_cell_mod, only: particles_in_cell
   USE physprop, only: k_g0, k_s0
   USE run, only: units, time
   USE stl
   USE stl_dbg_des, only: write_stls_this_dg
   USE stl_dbg_des, only: write_this_stl
   USE stl_functions_des, only: closestptpointtriangle

   use sq_contact_wall
   use SQ_EQUIVALENT_RADIUS_MOD
   use SQ_PROPERTIES_MOD

      use is

   PRIVATE
   PUBLIC :: CALC_DEM_FORCE_WITH_WALL_STL,&
      & CALC_DEM_THERMO_WITH_WALL_STL

CONTAINS

   SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

      Implicit none

      INTEGER :: LL
      INTEGER :: NF, NG
      DOUBLE PRECISION ::OVERLAP_N, SQRT_OVERLAP

      DOUBLE PRECISION :: V_REL_TRANS_NORM, DISTSQ, RADSQ, CLOSEST_PT(DIMN)
! local normal and tangential forces
      DOUBLE PRECISION :: NORMAL(DIMN), VREL_T(DIMN), DIST(DIMN), DISTMOD
      DOUBLE PRECISION, DIMENSION(DIMN) :: FTAN, FNORM, OVERLAP_T
      DOUBLE PRECISION :: ARM(3)

      LOGICAL :: DES_LOC_DEBUG
      INTEGER :: CELL_ID, cell_count
      INTEGER :: PHASELL

      DOUBLE PRECISION :: TANGENT(DIMN)
      DOUBLE PRECISION :: FTMD, FNMD
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W

      double precision :: MAG_OVERLAP_T
      double precision :: vel_mag

      double precision :: line_t
! flag to tell if the orthogonal projection of sphere center to
! extended plane detects an overlap

      DOUBLE PRECISION :: MAX_DISTSQ, DISTAPART, FORCE_COH, R_LM
      INTEGER :: MAX_NF, axis
      DOUBLE PRECISION, DIMENSION(3) :: PARTICLE_MIN, PARTICLE_MAX, POS_TMP
      DOUBLE PRECISION, DIMENSION(3) :: NORM_PREVIOUS_FACET_COH
      DOUBLE PRECISION, DIMENSION(3) :: COH_FORCE_PREVIOUS_FACET
      DOUBLE PRECISION :: COHMAG_PREVIOUS_FACET
      DOUBLE PRECISION :: NORMTEST
! Additional relative translational motion due to rotation
      DOUBLE PRECISION, DIMENSION(DIMN) :: V_ROT
! Wall tangential velocity at contact point
      double precision, dimension(3) :: wall_tang_vel

! Flag to keep only cohesion force with one STL facet
      LOGICAL :: COHESION_WITH_STL_FACET
! Flag to distinguish point or edge intersection
      logical :: point_or_edge_int, moving_wall
      integer :: edge

! GluedSphereDEM model
      INTEGER :: gid, num_contact
      DOUBLE PRECISION, DIMENSION(3) :: DIST2CENTERVEC, ftmp, ttmp, VRELTRANS, tmp_TNORM, tmp_TTAN
      DOUBLE PRECISION, DIMENSION(3) :: CONTACT_GLOBAL

! Superquadric DEM model
! Shape exponents of superquadric surface
      DOUBLE PRECISION :: m,n
! Euler parameters of superquadric surface, also called quaternions
      DOUBLE PRECISION :: euler0, euler1, euler2, euler3
! P is the location of superquadric center, axi is the semi-axes
      DOUBLE PRECISION :: p(3), axi(3),X(3)
      DOUBLE PRECISION :: contact_point_global(3),contact_point_wall_global(3)
      DOUBLE PRECISION :: norm_wall(3),vertex_wall(3),norm_wallsq,DIST2(3)
      DOUBLE PRECISION :: DIST3(3),dist_t
      integer :: contact_test_wall
! Spherical Chen-JKR wall contact and rolling-history state
      DOUBLE PRECISION :: A_OVER_A0, A_JKR
      INTEGER :: CC_JKR
      DOUBLE PRECISION :: JKR_NORM(DIMN)
      DOUBLE PRECISION :: JKR_A0_W_LOCAL, JKR_DELTAC_W_LOCAL
      DOUBLE PRECISION :: JKR_FC_W_LOCAL, JKR_GAMMA_W_LOCAL
      DOUBLE PRECISION :: JKR_E_EFF_W_LOCAL
      DOUBLE PRECISION :: PFT_OLD(DIMN)
      DOUBLE PRECISION :: V_ROLL_W(3), XI_MAG_W, XI_CRIT_W, KR_JKR_W
      DOUBLE PRECISION :: MR_TORQUE_W(3)
! equivalent rad of local contact curvature
      DOUBLE PRECISION :: eq_r,r_eff_bs,r_eff_sq
      DOUBLE PRECISION :: HERT_KwN_SUPER(DIM_M),HERT_KwT_SUPER(DIM_M)
      DOUBLE PRECISION :: DES_ETAN_WALL_SUPER(DIM_M),DES_ETAT_WALL_SUPER(DIM_M)
      INTEGER :: shape_ll
      logical, parameter :: ldebug = .false.
! Rolling friction
      DOUBLE PRECISION :: RFTOW(3), omega_mag
! Temporary variable
      DOUBLE PRECISION, DIMENSION(3) :: TMP
      ! Boundary Condition ID associated with STL triangle
      INTEGER :: bc_id
! IOSTAT
      INTEGER :: io_err

      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false.
      FOCUS_PARTICLE = -1

!--------------------------------------------------------------------------
! Glued sphere model: reset force, torque, glued
      IF (GSP_EXPLICIT) THEN
         glueForce(:,:) = ZERO
         glueTorque(:,:) =ZERO
      ENDIF

! Reset Particle-STL collision counter
      STL_COLL_COUNT(:) = ZERO
      STL_FNORM(:,:) = ZERO
      STL_FTAN(:,:) = ZERO
      STL_FORCE(:,:) = ZERO
      STL_TORQUE(:,:) = ZERO
      STL_CURRENT_WEAR(:) = ZERO

!$omp parallel default(none) private(gid,LL,MAG_OVERLAP_T,                &
!$omp    cell_id,radsq,particle_max,particle_min,tangent,                 &
!$omp    axis,nf,ng,closest_pt,dist,r_lm,distapart,force_coh,distsq,      &
!$omp    DIST2CENTERVEC, tmp_TNORM, tmp_TTAN,                             &
!$omp    line_t,max_distsq,max_nf,normal,distmod,overlap_n,VREL_T,        &
!$omp    v_rel_trans_norm,phaseLL,sqrt_overlap,kn_des_w,kt_des_w,         &
!$omp    cohesion_with_stl_facet, point_or_edge_int, edge, des_usr_var,   &
!$omp    norm_previous_facet_coh, cohmag_previous_facet,                  &
!$omp    coh_force_previous_facet, normtest, vel_mag,                     &
!$omp    etan_des_w,etat_des_w,fnorm,overlap_t,ftan,ftmd,fnmd,pos_tmp,    &
!$omp    rftow, omega_mag,euler0,euler1,euler2,euler3,cell_count,         &
!$omp    shape_ll,p,axi,m,n,norm_wall,norm_wallsq,vertex_wall,VRELTRANS,  &
!$omp    contact_test_wall,contact_point_wall_global,contact_point_global,&
!$omp    dist_t,eq_r,r_eff_sq,r_eff_bs,dist3, tmp,                        &
!$omp    bc_id, v_rot, num_contact, contact_global, arm,                 &
!$omp    a_over_a0,a_jkr,cc_jkr,jkr_norm,jkr_a0_w_local,                  &
!$omp    jkr_deltac_w_local,jkr_fc_w_local,jkr_gamma_w_local,             &
!$omp    jkr_e_eff_w_local,pft_old,v_roll_w,xi_mag_w,xi_crit_w,           &
!$omp    kr_jkr_w,mr_torque_w)                                            &
!$omp shared(max_pip,focus_particle,debug_des,                            &
!$omp    pijk,dg_pijk,des_pos_new,omega_new,pmass,                        &
!$omp    des_radius,facets_at_dg,vertex,mew_rw,                           &
!$omp    ignore_stl_edge_intersection, cos_stl_nb_angle,                  &
!$omp    wall_surface_energy,                                             &
!$omp    hert_kwn,hert_kwt,kn_w,kt_w,des_coll_model_enum,mew_w,tow,       &
!$omp    gluePos,des_vel_new,glueomg,                                     &
!$omp    glueforce,gluetorque,gid_list,iglobal_id,particle_state,         &
!$omp    des_etan_wall,des_etat_wall,dtsolid,fc,norm_face,                &
!$omp    wall_collision_facet_id,wall_collision_PFT,wall_collision_PFR,   &
!$omp    jkr_wall_contact_state,use_jkr,jkr_thetar_crit,jkr_alpha_safety,  &
!$omp    wall_jkr_thetar_crit,des_surface_energy_wall,                    &
!$omp    hert_kwn_jkr,hert_kwt_jkr,des_etan_wall_jkr,des_etat_wall_jkr,   &
!$omp    use_cohesion,                                                     &
!$omp    van_der_waals,wall_hamaker_constant,wall_vdw_outer_cutoff,       &
!$omp    wall_vdw_inner_cutoff,asperities,surface_energy,                 &
!$omp    super_q, super_mn, super_r, SuperDEM,GluedSphereDEM,             &
!$omp    gsp_implicit, gsp_explicit, rigid_gsp, flexible_gsp,             &
!$omp    des_etan_wall_super,des_etat_wall_super,                         &
!$omp    hert_kwn_super,hert_kwt_super,bc_id_stl_face,                    &
!$omp    stl_coll_count, stl_current_wear, stl_fnorm, stl_ftan,           &
!$omp    stl_force, stl_q, stl_torque, stl_group, stl_group_id)

!$omp do
      DO LL = 1, MAX_PIP
         IF(LL.EQ.FOCUS_PARTICLE) DEBUG_DES = .TRUE.
! skipping non-existent particles or ghost particles
! make sure the particle is not classified as a new 'entering' particle
! or is already marked as a potential exiting particle
         IF(.NOT.IS_NORMAL(LL)) CYCLE

         IF (GSP_EXPLICIT) THEN
             gid = gid_list(LL)
         ENDIF

         CELL_ID = DG_PIJK(LL)

! If no neighboring facet in the surrounding 26 cells, then exit
         IF(facets_at_dg(CELL_ID)%COUNT < 1) THEN
            WALL_COLLISION_FACET_ID(:,LL) = -1
            WALL_COLLISION_PFT(:,:,LL) = 0.0d0
            WALL_COLLISION_PFR(:,:,LL) = ZERO
            JKR_WALL_CONTACT_STATE(:,LL) = .FALSE.
            CYCLE
         ENDIF

         if (SuperDEM) then
! Quaternions for LL particle
             euler0= super_q(ll,1)
             euler1= super_q(ll,2)
             euler2= super_q(ll,3)
             euler3= super_q(ll,4)
! position of LL particle
             p(1)=DES_POS_NEW(LL,1)
             p(2)=DES_POS_NEW(LL,2)
             p(3)=DES_POS_NEW(LL,3)
! Semi-axis and roundness of LL particle
             axi(1)=super_r(ll,1)
             axi(2)=super_r(ll,2)
             axi(3)=super_r(ll,3)
             m=super_mn(ll,1)
             n=super_mn(ll,2)

! Check the shape of ll particle
             call Sq_shapes(axi,m,n,shape_ll)
         endif

! Check particle LL for wall contacts
         phaseLL = PIJK(LL,5)
         RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)
         JKR_FC_W_LOCAL = ZERO
         JKR_A0_W_LOCAL = ZERO
         JKR_DELTAC_W_LOCAL = ZERO
         IF(USE_JKR) THEN
            JKR_GAMMA_W_LOCAL = SQRT(GET_DES_SURFACE_ENERGY(LL)* &
               DES_SURFACE_ENERGY_WALL)
            JKR_E_EFF_W_LOCAL = 0.75D0*HERT_KWN_JKR(phaseLL)
            JKR_A0_W_LOCAL = (9.D0*PI*JKR_GAMMA_W_LOCAL* &
               DES_RADIUS(LL)**2/JKR_E_EFF_W_LOCAL)**(ONE/3.D0)
            JKR_DELTAC_W_LOCAL = JKR_A0_W_LOCAL**2 / &
               (2.D0*6.D0**(ONE/3.D0)*DES_RADIUS(LL))
            JKR_FC_W_LOCAL = 3.D0*PI*DES_RADIUS(LL)*JKR_GAMMA_W_LOCAL
         ENDIF

         particle_max(:) = des_pos_new( LL,:) + des_radius(LL)
         particle_min(:) = des_pos_new( LL,:) - des_radius(LL)
         IF(USE_JKR .AND. JKR_FC_W_LOCAL > ZERO) THEN
            particle_max(:) = particle_max(:) + JKR_DELTAC_W_LOCAL
            particle_min(:) = particle_min(:) - JKR_DELTAC_W_LOCAL
         ENDIF

         COHESION_WITH_STL_FACET  = .FALSE.

         DO CELL_COUNT = 1, facets_at_dg(cell_id)%count

            axis = facets_at_dg(cell_id)%dir(cell_count)

            NF = facets_at_dg(cell_id)%id(cell_count)

            BC_ID = BC_ID_STL_FACE(NF)

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN

               CALL ClosestPtPointTriangle(NF, DES_POS_NEW(LL,:),          &
                  VERTEX(:,:,NF), CLOSEST_PT(:),point_or_edge_int, edge)

! To avoid counting cohesion forces multiple times when several
! triangles define a plane (or near planar) surface, edge or point
! interaction is ignored.
              IF(.NOT.(point_or_edge_int.AND.IGNORE_STL_EDGE_INTERSECTION(NF, EDGE))) THEN

                  DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(LL,:)
                  DISTSQ = DOT_PRODUCT(DIST, DIST)
                  R_LM = 1*DES_RADIUS(LL)

                  IF(DISTSQ < (R_LM+WALL_VDW_OUTER_CUTOFF)**2) THEN
                     IF(DISTSQ > (WALL_VDW_INNER_CUTOFF+R_LM)**2) THEN
                        DistApart = (SQRT(DISTSQ)-R_LM)
                        FORCE_COH = WALL_HAMAKER_CONSTANT*DES_RADIUS(LL) /&
                           (6.0d0*DistApart**2)*(Asperities/(Asperities +  &
                           DES_RADIUS(LL)) + ONE/(ONE+Asperities/         &
                           DistApart)**2)
                     ELSE
                        FORCE_COH = 4.0d0*PI*WALL_SURFACE_ENERGY*DES_RADIUS(LL)* &
                           (Asperities/(Asperities+DES_RADIUS(LL)) + ONE/ &
                           (ONE+Asperities/WALL_VDW_INNER_CUTOFF)**2 )
                     ENDIF


! The logic below is also aimed at avoiding counting the contribution of
! a wall multiple times. Once the contribution of a facet is accounted
! for, the contribution of the next facet is either:
! - ignored if the two normal vectors are close to each other (angle below STL_NB_ANGLE defined in
!   stl_mod.f), and its cohesion force magnitude is smaller (because only the closest
!   point should be kept)
! - kept, if the normals are close to each other and the cohesion force
!   magnitude is larger. In that case the force from the previous facet is
!   removed.
! If the facet is not part tof the same plane (say near a corner), the
! its contribution is added
                     IF(.NOT.COHESION_WITH_STL_FACET) THEN ! Record contribution of the first active facet
                       COHESION_WITH_STL_FACET = .TRUE.
                       NORM_PREVIOUS_FACET_COH = NORM_FACE(:,NF)
                       COHMAG_PREVIOUS_FACET = FORCE_COH/SQRT(DISTSQ)
                       COH_FORCE_PREVIOUS_FACET = DIST(:)*FORCE_COH/SQRT(DISTSQ)
                     ELSE
                        NORMTEST = DOT_PRODUCT(NORM_PREVIOUS_FACET_COH,NORM_FACE(:,NF))
                        IF(NORMTEST>COS_STL_NB_ANGLE) THEN
                           IF((FORCE_COH/SQRT(DISTSQ))>COHMAG_PREVIOUS_FACET) THEN
                              FC(LL,:) = FC(LL,:) - COH_FORCE_PREVIOUS_FACET ! Remove contribution from previous facet
                              COHESION_WITH_STL_FACET = .TRUE. ! and record contribution of current facet
                              NORM_PREVIOUS_FACET_COH = NORM_FACE(:,NF)
                              COHMAG_PREVIOUS_FACET = FORCE_COH/SQRT(DISTSQ)
                              COH_FORCE_PREVIOUS_FACET = DIST(:)*FORCE_COH/SQRT(DISTSQ)
                           ELSE
                              COHESION_WITH_STL_FACET = .FALSE.
                           ENDIF
                        ELSE  ! Keep adding contribution of current facet because it is not is he same plane (e.g, near a corner)
                           COHESION_WITH_STL_FACET = .TRUE.
                           NORM_PREVIOUS_FACET_COH = NORM_FACE(:,NF)
                           COHMAG_PREVIOUS_FACET = FORCE_COH/SQRT(DISTSQ)
                           COH_FORCE_PREVIOUS_FACET = DIST(:)*FORCE_COH/SQRT(DISTSQ)
                        ENDIF
                     ENDIF

                     IF(COHESION_WITH_STL_FACET) FC(LL,:) = FC(LL,:) + DIST(:)*FORCE_COH/SQRT(DISTSQ)


                  ENDIF
               ENDIF
            ENDIF

            if (facets_at_dg(cell_id)%min(cell_count) >    &
               particle_max(axis)) then
               call remove_collision(LL, nf, wall_collision_facet_id)
               cycle
            endif

            if (facets_at_dg(cell_id)%max(cell_count) <    &
               particle_min(axis)) then
               call remove_collision(LL, nf, wall_collision_facet_id)
               cycle
            endif

! Checking all the facets is time consuming due to the expensive
! separating axis test. Remove this facet from contention based on
! a simple orthogonal projection test.

! Parametrize a line as p = p_0 + t normal and intersect with the
! triangular plane. If t>0, then point is on the non-fluid side of
! the plane. If the plane normal is assumed to point toward the fluid.

! -undefined, because non zero values will imply the sphere center
! is on the non-fluid side of the plane. Since the testing
! is with extended plane, this could very well happen even
! when the particle is well inside the domain (assuming the plane
! normal points toward the fluid). See the pic below. So check
! only when line_t is negative

!                 \   Solid  /
!                  \  Side  /
!                   \      /
!                    \    /
! Wall 1, fluid side  \  /  Wall 2, fluid side
!                      \/
!                        o particle
!
! line_t will be positive for wall 1 (incorrectly indicating center
! is outside the domain) and line_t will be negative for wall 2.
!
! Therefore, only stick with this test when line_t is negative and let
! the separating axis test take care of the other cases.

! Since this is for checking static config, line's direction is the
! same as plane's normal. For moving particles, the line's normal will
! be along the point joining new and old positions.

            line_t = DOT_PRODUCT(VERTEX(1,:,NF) - des_pos_new(LL,:),&
               NORM_FACE(:,NF))

! k - rad >= tol_orth, where k = -line_t, then orthogonal
! projection is false. Substituting for k
! => line_t + rad <= -tol_orth
! choosing tol_orth = 0.01% of des_radius = 0.0001*des_radius

! Orthogonal projection will detect false positives even
! when the particle does not overlap the triangle.
! However, if the orthogonal projection shows no overlap, then
! that is a big fat negative and overlaps are not possible.

            IF(USE_JKR .AND. JKR_FC_W_LOCAL > ZERO) THEN
               IF(line_t <= -(DES_RADIUS(LL)+JKR_DELTAC_W_LOCAL)) THEN
                  CALL REMOVE_COLLISION(LL,NF,WALL_COLLISION_FACET_ID)
                  CYCLE
               ENDIF
            ELSEIF(line_t <= -1.0001D0*DES_RADIUS(LL)) THEN
               CALL REMOVE_COLLISION(LL,NF,WALL_COLLISION_FACET_ID)
               CYCLE
            ENDIF

! SuperDEM
            IF(SuperDEM .and. shape_ll .gt. 0) THEN
! Facet wall normal point to the fluid domain side
                norm_wall(:)=NORM_FACE(:,NF)
                norm_wallsq=(norm_wall(1))**2+(norm_wall(2))**2+(norm_wall(3))**2
                norm_wall(:)=norm_wall(:)/dsqrt(norm_wallsq)
                vertex_wall(:)=VERTEX(1,:,NF)

                call SQ_CONTACT_WITH_STL(p,axi,m,n,&
                    euler0,euler1, euler2, euler3,NORM_wall,&
                    vertex_wall,contact_point_global,&
                    contact_point_wall_global,contact_test_wall)

                pos_tmp(:) = contact_point_global(:)
                CALL ClosestPtPointTriangle(NF, pos_tmp,             &
                    VERTEX(:,:,NF), CLOSEST_PT(:),point_or_edge_int, edge)

                DIST(:) = contact_point_global(:)-contact_point_wall_global(:)
               dist_t=dot_product(dist(:), norm_wall(:))

               DISTSQ = dist_t**2.0


               if (contact_test_wall > 0) then  !not contact

                   call remove_collision(LL,nf,wall_collision_facet_id)
                   CYCLE
               endif

            ELSE ! Sphere (including special case when SQP is a sphere)

               pos_tmp = DES_POS_NEW(LL,:)
               CALL ClosestPtPointTriangle(NF, pos_tmp,             &
                   VERTEX(:,:,NF), CLOSEST_PT(:),point_or_edge_int, edge)

               DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(LL,:)
               DISTSQ = DOT_PRODUCT(DIST, DIST)

               IF(USE_JKR .AND. JKR_FC_W_LOCAL > ZERO) THEN
                  IF(DISTSQ > (DES_RADIUS(LL)+JKR_ALPHA_SAFETY* &
                     JKR_DELTAC_W_LOCAL)**2) THEN
                     CALL REMOVE_COLLISION(LL,NF,WALL_COLLISION_FACET_ID)
                     CYCLE
                  ENDIF
               ELSEIF(DISTSQ .GE. RADSQ - SMALL_NUMBER) THEN !No overlap exists
                  CALL REMOVE_COLLISION(LL,NF,WALL_COLLISION_FACET_ID)
                  CYCLE
               ENDIF

            ENDIF !SuperDEM
! If angle between two facets is small (below STL_NB_ANGLE defined in
! stl_mod.f), the collision with edge is ignored. This is done to avoid
! unwanted collision with all triangle edges making up a large plane or
! surface with low curvature
            IF(point_or_edge_int.AND.IGNORE_STL_EDGE_INTERSECTION(NF, EDGE)) CYCLE

            MAX_DISTSQ = DISTSQ
            MAX_NF = NF

! Assign the collision normal based on the facet with the
! largest overlap.

! Facet's normal is correct normal only when the intersection is with
! the face. When the intersection is with edge or vertex, then the
! normal is based on closest pt and sphere center. The definition above
! of the normal is generic enough to account for differences between
! vertex, edge, and facet.

! Calculate the particle/wall overlap and V_ROT.
! V_ROT is the additional relative translational motion due to rotation.
! It is computed here because it will be passed to CFRELVEL_WALL below.
! It is computed differently for Spheres, SQP and GSP.
            DISTMOD = SQRT(MAX_DISTSQ)
            IF(SUPERDEM .AND. SHAPE_LL .GT. 0) THEN
               OVERLAP_N = DISTMOD
               NORMAL(:)= - NORM_WALL(:)
               TMP = OMEGA_NEW(LL,:)
               V_ROT(:) = CROSS(TMP, POS_TMP(:) - P(:))
            ELSEIF(rigid_gsp) THEN
               NORMAL(:) = DIST(:)/SQRT(DISTSQ)
               OVERLAP_N = DES_RADIUS(LL) - DISTMOD

               DIST2CENTERVEC(:) = 0.0
               CALL GSP_MIN_IMAGE_VECTOR(CLOSEST_PT, GLUEPOS(GID,:), DIST2CENTERVEC)
               IF (DOT_PRODUCT(DIST, DIST2CENTERVEC) < 0) THEN
                  NORMAL(:) = -NORMAL(:)
                  OVERLAP_N = DES_RADIUS(LL) + DISTMOD
               ENDIF
               TMP = GLUEOMG(GID,:)
               V_ROT(:) =  CROSS(TMP,DIST2CENTERVEC)
            ELSEIF(GSP_IMPLICIT) THEN
               CALL GSP_WALL_CONTACT_IMPLICIT(LL, NF, CLOSEST_PT, OVERLAP_N, NUM_CONTACT, NORMAL, CONTACT_GLOBAL)
               IF(NUM_CONTACT > 0) THEN
                  ! OVERLAP_N = OVERLAP_N/NUM_CONTACT
                  NORMAL = NORMAL/NUM_CONTACT
                  NORMAL = -NORMAL/DSQRT(NORMAL(1)**2+NORMAL(2)**2+NORMAL(3)**2)
                  CONTACT_GLOBAL = CONTACT_GLOBAL/NUM_CONTACT
                  DIST2CENTERVEC = CONTACT_GLOBAL - DES_POS_NEW(LL,:)
               ELSE
                  CALL REMOVE_COLLISION(LL,NF,WALL_COLLISION_FACET_ID)
                  CYCLE
               ENDIF
               TMP = OMEGA_NEW(LL,:)
               V_ROT(:) =  CROSS(TMP,DIST2CENTERVEC)
            ELSE
               NORMAL(:) = DIST(:)
               if(DISTSQ > 0) then
                  NORMAL(:) = DIST(:)/SQRT(DISTSQ)
               endif
               OVERLAP_N = DES_RADIUS(LL) - DISTMOD
               TMP = OMEGA_NEW(LL,:)
               V_ROT(:)  = DISTMOD*CROSS(TMP, NORMAL)
            ENDIF

! Calculate the translational relative velocity
! The same routine works with DEM, SQP, CGP and GSP due to the computation of
! V_ROT above
            CALL CFRELVEL_WALL(BC_ID, NF, LL, CLOSEST_PT, NORMAL, V_ROT, V_REL_TRANS_NORM, VREL_T)

! Obtain the facet slot before evaluating tensile JKR contact history.
            PFT_OLD(:) = GET_COLLISION(LL,NF,WALL_COLLISION_FACET_ID, &
               WALL_COLLISION_PFT,CC_JKR)

! Calculate the spring model parameters.
            phaseLL = PIJK(LL,5)

! Hertz vs linear spring-dashpot contact model
            IF(USE_JKR .AND. JKR_FC_W_LOCAL > ZERO) THEN
               IF(OVERLAP_N >= ZERO) THEN
                  JKR_WALL_CONTACT_STATE(CC_JKR,LL) = .TRUE.
               ELSEIF(.NOT.JKR_WALL_CONTACT_STATE(CC_JKR,LL)) THEN
                  CALL REMOVE_COLLISION(LL,NF,WALL_COLLISION_FACET_ID)
                  CYCLE
               ENDIF
               A_OVER_A0 = NEWTON_RAPHSON_METHOD(OVERLAP_N / &
                  (6.D0**(ONE/3.D0)*JKR_DELTAC_W_LOCAL))
               A_JKR = A_OVER_A0*JKR_A0_W_LOCAL
               KN_DES_W = HERT_KWN_JKR(phaseLL)*A_JKR
               KT_DES_W = HERT_KWT_JKR(phaseLL)*A_JKR
               ETAN_DES_W = DES_ETAN_WALL_JKR(phaseLL)*SQRT(A_JKR)
               ETAT_DES_W = DES_ETAT_WALL_JKR(phaseLL)*SQRT(A_JKR)
               JKR_NORM(:) = -4.D0*JKR_FC_W_LOCAL* &
                  (A_OVER_A0**3-A_OVER_A0**1.5D0)*NORMAL(:)
            ELSE
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = 0.0D0
               if (overlap_n > 0) sqrt_overlap = SQRT(OVERLAP_N)

               IF(SuperDEM .and. shape_ll .gt. 0) THEN
! Equivalent radius of particle A at contact point
                  call EQUIVALENT_RADIUS(p,contact_point_global,axi,m,n,&
                       euler0,euler1,euler2,euler3, &
                       eq_r)

! Update the value of R_EFF for superquadric particle
                  R_EFF_sq = eq_r
                  R_EFF_bs=DES_RADIUS(LL)

! Update some parameters for superquadric particle
!{SBAN
                  hert_kWn_SUPER(phaseLL)=hert_kWn(phaseLL)*SQRT(R_eff_sq)
                  hert_kWt_SUPER(phaseLL)=hert_kWt(phaseLL)*SQRT(R_eff_sq)

                  DES_ETAN_wall_SUPER(phaseLL)=DES_ETAN_WALL(phaseLL)*&
                     SQRT(HERT_KwN_SUPER(phaseLL))*sqrt(PMASS(LL))
                  DES_ETAT_wall_SUPER(phaseLL)=DES_ETAT_WALL(phaseLL)*&
                     SQRT(HERT_KwT_SUPER(phaseLL))*sqrt(PMASS(LL))
!}

                  KN_DES_W = hert_kWn_super(phaseLL)*sqrt_overlap
                  KT_DES_W = hert_kWt_super(phaseLL)*sqrt_overlap
                  ETAN_DES_W = DES_ETAN_WALL_super(phaseLL)*sqrt(sqrt_overlap)
                  ETAT_DES_W = DES_ETAT_WALL_super(phaseLL)*sqrt(sqrt_overlap)

               ELSE  ! for sphere, hertz
!{SBAN
                  KN_DES_W = hert_kwn(phaseLL)*sqrt_overlap * sqrt(DES_RADIUS(LL))
                  KT_DES_W = hert_kwt(phaseLL)*sqrt_overlap * sqrt(DES_RADIUS(LL))
                  ETAN_DES_W = DES_ETAN_WALL(phaseLL)*sqrt(sqrt_overlap) * sqrt(hert_kwn(phaseLL)*sqrt(DES_RADIUS(LL))*PMASS(LL))
                  ETAT_DES_W = DES_ETAT_WALL(phaseLL)*sqrt(sqrt_overlap) * sqrt(hert_kwt(phaseLL)*sqrt(DES_RADIUS(LL))*PMASS(LL))
               ENDIF ! end of SuperDEM .and. shape_ll .gt. 0

            ELSE ! for sphere, lsd
               KN_DES_W = KN_W
               KT_DES_W = KT_W
               ETAN_DES_W = DES_ETAN_WALL(phaseLL) * sqrt(PMASS(LL))
               ETAT_DES_W = DES_ETAT_WALL(phaseLL) * sqrt(PMASS(LL))
!}
            ENDIF ! end of DES_COLL_MODEL_ENUM .EQ. HERTZIAN
            ENDIF

! Calculate the normal contact force
            IF(USE_JKR .AND. JKR_FC_W_LOCAL > ZERO) THEN
               FNORM(:) = JKR_NORM(:) - &
                  ETAN_DES_W*V_REL_TRANS_NORM*NORMAL(:)
            ELSE
               JKR_WALL_CONTACT_STATE(CC_JKR,LL) = .FALSE.
               WALL_COLLISION_PFR(:,CC_JKR,LL) = ZERO
               FNORM(:) = -(KN_DES_W * OVERLAP_N * NORMAL(:) + &
                  ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:))
            ENDIF
            if (ldebug) then
               WRITE(*,*) 'norm force and v_rel_trans_norm'
               WRITE(*,'(3(f18.5))')   FNORM(:)
               WRITE(*,'(3(f18.5))') V_REL_TRANS_NORM,overlap_n
            endif
! Calculate the tangential displacement.
            OVERLAP_T(:) = DTSOLID*VREL_T(:) + PFT_OLD(:)
            MAG_OVERLAP_T = sqrt(DOT_PRODUCT(OVERLAP_T, OVERLAP_T))
            if (ldebug) then
               WRITE(*,*) 'overlaps_t and mag_overlap_t'
               WRITE(*,'(3(f18.5))')   overlap_t
               WRITE(*,'(3(f18.5))') mag_overlap_t
            endif

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Calculate the tangential contact force.
               FTAN = -KT_DES_W*OVERLAP_T - ETAT_DES_W*VREL_T
               FTMD = sqrt(dot_product(FTAN,FTAN))
! Max force before the on set of frictional slip.
               IF(USE_JKR .AND. JKR_FC_W_LOCAL > ZERO) THEN
                  FNMD = MEW_W*ABS(-DOT_PRODUCT(JKR_NORM,NORMAL) + &
                     2.D0*JKR_FC_W_LOCAL)
               ELSE
                  FNMD = MEW_W*sqrt(dot_product(FNORM,FNORM))
               ENDIF
! Frictional slip
               IF(FTMD > FNMD) THEN
! Direction of tangential force.
                  IF(USE_JKR .AND. JKR_FC_W_LOCAL > ZERO) THEN
                     FTAN = (FNMD/FTMD)*FTAN
                     OVERLAP_T = -(FTAN+ETAT_DES_W*VREL_T)/KT_DES_W
                  ELSE
                     TANGENT = OVERLAP_T/MAG_OVERLAP_T
                     FTAN = -FNMD * TANGENT
                     OVERLAP_T = (FNMD/KT_DES_W) * TANGENT
                  ENDIF
               ENDIF
            ELSE
               FTAN = 0.0
            ENDIF

! Store STL forces and torque
            STL_COLL_COUNT(NF) = STL_COLL_COUNT(NF) + 1.0
            VEL_MAG = DOT_PRODUCT(DES_VEL_NEW(LL,:), DES_VEL_NEW(LL,:))
! As a simple measure of wear on an STL surface, record the cumulative
! overlap*dtsolid, as long as particles are moving above some arbitrary
! threshold (so that the wear doesn't keep on increasing when a bed of
! particles is static).
            IF(VEL_MAG>1.0D-08) STL_CURRENT_WEAR(NF) = STL_CURRENT_WEAR(NF) + OVERLAP_N * DTSOLID
            STL_FNORM(NF,:) = STL_FNORM(NF,:) - FNORM(:)
            STL_FTAN(NF,:) = STL_FTAN(NF,:) - FTAN(:)
            STL_FORCE(NF,:) = STL_FORCE(NF,:) - (FNORM(:) + FTAN(:))
            NG = STL_GROUP_ID(NF)
            ARM = CLOSEST_PT(:) - STL_GROUP(NG)%TORQUE_PIVOT(:)
            STL_TORQUE(NF,:) = STL_TORQUE(NF,:) - CROSS(ARM,(FTAN+FNORM))


! Save the tangential displacement.
            CALL UPDATE_COLLISION(OVERLAP_T, LL, NF,                   &
               WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)

            IF (GSP_EXPLICIT .OR. GSP_IMPLICIT) THEN
!$omp atomic
               FC(LL,1) = FC(LL,1) + FNORM(1) + FTAN(1)
!$omp atomic
               FC(LL,2) = FC(LL,2) + FNORM(2) + FTAN(2)
!$omp atomic
               FC(LL,3) = FC(LL,3) + FNORM(3) + FTAN(3)
            ELSE
               FC(LL,:) = FC(LL,:) + FNORM(:) + FTAN(:)
            ENDIF

! Add the torque force to the total torque acting on the particle.
            IF (SuperDEM .and. shape_ll .gt. 0) THEN
               DIST3(:) = contact_point_wall_global(:) - DES_POS_NEW(LL,:)
               TOW(LL,:) = TOW(LL,:) + CROSS(DIST3,(FTAN+FNORM))
            ! @renjieke, for flexible GSP, should the torque in the component be calculated
            ! based on the center vector?
            ELSEIF (rigid_gsp .OR. GSP_IMPLICIT) THEN
               tmp_TNORM(:) = CROSS(DIST2CENTERVEC,FNORM)
               tmp_TTAN(:) = CROSS(DIST2CENTERVEC,FTAN)
!$omp atomic
               TOW(LL,1) = TOW(LL,1) + tmp_TNORM(1) + tmp_TTAN(1)
!$omp atomic
               TOW(LL,2) = TOW(LL,2) + tmp_TNORM(2) + tmp_TTAN(2)
!$omp atomic
               TOW(LL,3) = TOW(LL,3) + tmp_TNORM(3) + tmp_TTAN(3)
            ELSE
               TOW(LL,:) = TOW(LL,:) + DISTMOD*CROSS(NORMAL,FTAN)
            ENDIF

!Rolling friction
! JFD: Here the rolling friction coefficient is non-dimensional. This is the equivalent to
! the definition from Zhou (Physica A 269 (1999) 536-553) divided by the
! particle diameter. Therefore here mew_rw is multiplied by the diameter.

! ---------------------------------------
! OMEGA_NEW(LL,:) is calculated in cfnewvalues.f subroutine CFNEWVALUES
! so Lu did not call CFNEWVALUES, so he did not have OMEGA_NEW(LL,:) calculated.
! have to find another way to adding this rolling friction
! may use glueOmg instead
! ---------------------------------------
            IF(USE_JKR .AND. JKR_FC_W_LOCAL > ZERO) THEN
               V_ROLL_W(:) = -DES_RADIUS(LL)* &
                  CROSS(OMEGA_NEW(LL,:),NORMAL(:)) - 0.5D0*VREL_T(:)
               WALL_COLLISION_PFR(:,CC_JKR,LL) = &
                  WALL_COLLISION_PFR(:,CC_JKR,LL) + DTSOLID*V_ROLL_W(:)
               XI_MAG_W = SQRT(DOT_PRODUCT(WALL_COLLISION_PFR(:,CC_JKR,LL), &
                  WALL_COLLISION_PFR(:,CC_JKR,LL)))
               XI_CRIT_W = WALL_JKR_THETAR_CRIT*DES_RADIUS(LL)
               IF(XI_MAG_W > XI_CRIT_W .AND. XI_MAG_W > ZERO) &
                  WALL_COLLISION_PFR(:,CC_JKR,LL) = XI_CRIT_W* &
                  WALL_COLLISION_PFR(:,CC_JKR,LL)/XI_MAG_W
               KR_JKR_W = 4.D0*JKR_FC_W_LOCAL*A_OVER_A0**1.5D0
               MR_TORQUE_W(:) = -KR_JKR_W* &
                  CROSS(NORMAL(:),WALL_COLLISION_PFR(:,CC_JKR,LL))
               TOW(LL,:) = TOW(LL,:) + MR_TORQUE_W(:)
            ELSE
               IF(MEW_RW>ZERO .AND. .NOT. GluedSphereDEM) THEN
                  OMEGA_MAG = sqrt(dot_product(OMEGA_NEW(LL,:),OMEGA_NEW(LL,:)))
                  if(OMEGA_MAG > ZERO) then
                    RFTOW = -MEW_RW*(2.0*DES_RADIUS(LL))*KN_DES_W * OVERLAP_N*OMEGA_NEW(LL,:)/OMEGA_MAG
                    TOW(LL,:) = TOW(LL,:) + RFTOW
                  endif
               ENDIF
            ENDIF

       ENDDO ! end of wall search iteration
   ENDDO ! end of sphere iteration
!$omp end do
!$omp end parallel

   IF(rigid_gsp) THEN
      DO LL = 1, MAX_PIP
         if(particle_state(LL) .NE. 1) CYCLE
         gid = gid_list(LL)
         glueForce(gid,:) = glueForce(gid,:) + FC(LL,:)
         glueTorque(gid,:) = glueTorque(gid,:) + TOW(LL,:)
      ENDDO
   ENDIF

! Sum across ranks
   call global_all_sum(stl_coll_count)
   call global_all_sum(stl_fnorm)
   call global_all_sum(stl_ftan)
   call global_all_sum(stl_force)
   call global_all_sum(stl_torque)
   call global_all_sum(stl_current_wear)

   stl_wear(:) = stl_wear(:) + stl_current_wear(:)

   RETURN

contains

      SUBROUTINE GSP_MIN_IMAGE_VECTOR(ref_pos, center_pos, out_vec)

      DOUBLE PRECISION, INTENT(IN) :: ref_pos(3), center_pos(3)
      DOUBLE PRECISION, INTENT(OUT) :: out_vec(3)

      out_vec(:) = ref_pos(:) - center_pos(:)

      if (cyclic_x .and. xlength > zero) then
         do while (out_vec(1) >  half*xlength)
            out_vec(1) = out_vec(1) - xlength
         end do
         do while (out_vec(1) < -half*xlength)
            out_vec(1) = out_vec(1) + xlength
         end do
      endif

      if (cyclic_y .and. ylength > zero) then
         do while (out_vec(2) >  half*ylength)
            out_vec(2) = out_vec(2) - ylength
         end do
         do while (out_vec(2) < -half*ylength)
            out_vec(2) = out_vec(2) + ylength
         end do
      endif

      if (cyclic_z .and. zlength > zero) then
         do while (out_vec(3) >  half*zlength)
            out_vec(3) = out_vec(3) - zlength
         end do
         do while (out_vec(3) < -half*zlength)
            out_vec(3) = out_vec(3) + zlength
         end do
      endif

      END SUBROUTINE GSP_MIN_IMAGE_VECTOR

! Return the stable JKR contact-radius branch for a normalized overlap.
      DOUBLE PRECISION FUNCTION NEWTON_RAPHSON_METHOD(Y)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: Y
      DOUBLE PRECISION :: X_OLD, X_NEW, FX, DFX
      INTEGER :: ITER

      X_OLD = ONE
      X_NEW = X_OLD
      DO ITER=1,1000
         FX = 2.D0*X_NEW**2 - (4.D0/3.D0)*SQRT(X_NEW) - Y
         DFX = 4.D0*X_NEW - (2.D0/3.D0)/SQRT(X_NEW)
         IF(ABS(DFX) < 1.D-15) THEN
!$omp critical (jkr_nr_warn_wall)
            WRITE(*,'(A,ES12.5,A,ES12.5,A,I0)') &
               'WARNING: JKR wall Newton derivative small; y=',Y, &
               ' x=',X_NEW,' iter=',ITER
!$omp end critical (jkr_nr_warn_wall)
            NEWTON_RAPHSON_METHOD = X_NEW
            RETURN
         ENDIF
         X_NEW = MAX(X_OLD-FX/DFX,1.D-12)
         IF(ABS((X_NEW-X_OLD)/X_NEW) < 1.D-8) THEN
            NEWTON_RAPHSON_METHOD = X_NEW
            RETURN
         ENDIF
         X_OLD = X_NEW
      ENDDO
!$omp critical (jkr_nr_warn_wall)
      WRITE(*,'(A,ES12.5,A,ES12.5)') &
         'WARNING: JKR wall Newton did not converge; y=',Y,' x=',X_NEW
!$omp end critical (jkr_nr_warn_wall)
      NEWTON_RAPHSON_METHOD = X_NEW
      END FUNCTION NEWTON_RAPHSON_METHOD

!......................................................................!
!  Function: GET_COLLISION                                             !
!                                                                      !
!  Purpose: Return the integrated (t0->t) tangential displacement.     !
!......................................................................!
     FUNCTION GET_COLLISION(LLL,FACET_ID,WALL_COLLISION_FACET_ID,     &
         WALL_COLLISION_PFT,CC_OUT)

          IMPLICIT NONE

      DOUBLE PRECISION :: GET_COLLISION(DIMN)
      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(OUT) :: CC_OUT
      INTEGER, allocatable, INTENT(INOUT) :: WALL_COLLISION_FACET_ID(:,:)
      DOUBLE PRECISION, allocatable, INTENT(INOUT) :: WALL_COLLISION_PFT(:,:,:)
      INTEGER :: CC, FREE_INDEX, LC, dgIJK


      free_index = -1

      do cc = 1, COLLISION_ARRAY_MAX
         if (facet_id == wall_collision_facet_id(cc,LLL)) then
            CC_OUT = CC
            get_collision(:) = wall_collision_PFT(:,cc,LLL)
            return
         else if (-1 == wall_collision_facet_id(cc,LLL)) then
            free_index = cc
         endif
      enddo

! Overwrite old data. This is needed because a particle moving from
! one dg cell to another may no longer 'see' an STL before it moved
! out of contact range. Therefore, the 'remove_collision' function
! does not get called to cleanup the stale data.
      if(-1 == free_index) then
         dgIJK=DG_PIJK(LLL)
         cc_lp: do cc=1, COLLISION_ARRAY_MAX
            do lc=1, facets_at_dg(dgIJK)%count
               if(wall_collision_facet_id(cc,LLL) == &
                  facets_at_dg(dgIJK)%id(LC))  cycle cc_lp
            enddo
            free_index = cc
            exit cc_lp
         enddo cc_lp
      endif

! Last resort... grow the collision array
      if(-1 == free_index) then
! The DMP particle packet is sized from COLLISION_ARRAY_MAX during
! initialization. A rank-local resize would make packet layouts differ
! across ranks, so fail explicitly instead of corrupting wall history.
         IF(USE_JKR .AND. NUMPES > 1) THEN
            WRITE(ERR_MSG, 1090) COLLISION_ARRAY_MAX
            CALL LOG_ERROR()
         ENDIF
         free_index=COLLISION_ARRAY_MAX+1
         COLLISION_ARRAY_MAX = 2*COLLISION_ARRAY_MAX
         CALL GROW_WALL_COLLISION(COLLISION_ARRAY_MAX)
      endif

      CC_OUT = FREE_INDEX
      wall_collision_facet_id(free_index,LLL) = facet_id
      wall_collision_PFT(:,free_index,LLL) = ZERO
      wall_collision_PFR(:,free_index,LLL) = ZERO
      JKR_WALL_CONTACT_STATE(free_index,LLL) = .FALSE.
      get_collision(:) = wall_collision_PFT(:,free_index,LLL)
      return

 1090 FORMAT('JKR DMP wall-history capacity exhausted at ',I0, &
         ' facet slots per particle. Increase the initial capacity ', &
         'before MPI initialization.')

      END FUNCTION GET_COLLISION


!......................................................................!
!  Subroutine: GROW_WALL_COLLISION                                     !
!                                                                      !
!  Purpose: Return the integrated (t0->t) tangential displacement.     !
!......................................................................!
      SUBROUTINE GROW_WALL_COLLISION(NEW_SIZE)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NEW_SIZE
      INTEGER :: lSIZE1, lSIZE2, lSIZE3
      INTEGER, ALLOCATABLE :: tmpI2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: tmpR3(:,:,:)
      LOGICAL, ALLOCATABLE :: tmpL2(:,:)

      lSIZE1 = size(wall_collision_facet_id,1)
      lSIZE2 = size(wall_collision_facet_id,2)

      allocate(tmpI2(NEW_SIZE, lSIZE2))
      tmpI2(1:lSIZE1,:) = WALL_COLLISION_FACET_ID(1:lSIZE1,:)
      call move_alloc(tmpI2, WALL_COLLISION_FACET_ID)
      WALL_COLLISION_FACET_ID(lSIZE1+1:NEW_SIZE,:) = -1

      lSIZE1 = size(wall_collision_pft,1)
      lSIZE2 = size(wall_collision_pft,2)
      lSIZE3 = size(wall_collision_pft,3)

      allocate(tmpR3(lSIZE1, NEW_SIZE, lSIZE3))
      tmpR3(:,1:lSIZE2,:) = WALL_COLLISION_PFT(:,1:lSIZE2,:)
      tmpR3(:,lSIZE2+1:NEW_SIZE,:) = ZERO
      call move_alloc(tmpR3, WALL_COLLISION_PFT)

      allocate(tmpR3(lSIZE1, NEW_SIZE, lSIZE3))
      tmpR3(:,1:lSIZE2,:) = WALL_COLLISION_PFR(:,1:lSIZE2,:)
      tmpR3(:,lSIZE2+1:NEW_SIZE,:) = ZERO
      call move_alloc(tmpR3, WALL_COLLISION_PFR)

      allocate(tmpL2(NEW_SIZE,lSIZE3))
      tmpL2(1:lSIZE2,:) = JKR_WALL_CONTACT_STATE(1:lSIZE2,:)
      tmpL2(lSIZE2+1:NEW_SIZE,:) = .FALSE.
      call move_alloc(tmpL2,JKR_WALL_CONTACT_STATE)

      RETURN
      END SUBROUTINE GROW_WALL_COLLISION




!......................................................................!
!  Function: UPDATE_COLLISION                                          !
!                                                                      !
!  Purpose: Update the integrated (t0->t) tangential displacement.     !
!......................................................................!
      SUBROUTINE UPDATE_COLLISION(PFT, LLL, FACET_ID,                  &
         WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)

      implicit none

      DOUBLE PRECISION, INTENT(IN) :: PFT(DIMN)
      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(IN) :: WALL_COLLISION_FACET_ID(:,:)
      DOUBLE PRECISION, INTENT(INOUT) :: WALL_COLLISION_PFT(:,:,:)
      INTEGER :: CC

      do cc = 1, COLLISION_ARRAY_MAX
         if (facet_id == wall_collision_facet_id(cc,LLL)) then
            wall_collision_PFT(:,cc,LLL) = PFT(:)
            return
         endif
      enddo

      WRITE(ERR_MSG, 1100)
      CALL LOG_ERROR()

 1100 FORMAT('Error: COLLISION_ARRAY_MAX too small. ')

      END SUBROUTINE UPDATE_COLLISION

!......................................................................!
!  Function: REMOVE_COLLISION                                          !
!                                                                      !
!  Purpose: Clear the integrated (t0->t) tangential displacement once  !
!  the collision is over (contact ended).                              !
!......................................................................!
      SUBROUTINE REMOVE_COLLISION(LLL,FACET_ID,WALL_COLLISION_FACET_ID)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(INOUT) :: WALL_COLLISION_FACET_ID(:,:)
      INTEGER :: CC

      DO CC = 1, COLLISION_ARRAY_MAX
         IF (FACET_ID == WALL_COLLISION_FACET_ID(CC,LLL)) THEN
            WALL_COLLISION_FACET_ID(CC,LLL) = -1
            WALL_COLLISION_PFT(:,CC,LLL) = ZERO
            WALL_COLLISION_PFR(:,CC,LLL) = ZERO
            JKR_WALL_CONTACT_STATE(CC,LLL) = .FALSE.
            RETURN
         ENDIF
      ENDDO

      END SUBROUTINE REMOVE_COLLISION

      END SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CFRELVEL_WALL                                           !
!                                                                      !
!  Purpose: Calculate the normal and tangential components of the      !
!  relative velocity between a particle and wall contact.              !
!                                                                      !
!  Comments: Only the magnitude of the normal component is returned    !
!  whereas the full tangential vector is returned.                     !
!----------------------------------------------------------------------!
      SUBROUTINE CFRELVEL_WALL(BC_ID, N, LL, CLOSEST_PT, NORM, V_ROT, VRN, VRT)
      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Boundary Condition ID
      INTEGER, INTENT(IN) :: BC_ID
! STL Facet ID
      INTEGER, INTENT(IN) :: N
! Particle index.
      INTEGER, INTENT(IN) :: LL
! Point closest to the wall
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN)   :: CLOSEST_PT
! Unit normal from particle center to closest point on stl (wall)
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: NORM
! Additional relative translational motion due to rotation
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: V_ROT
! Magnitude of the total relative translational velocity.
      DOUBLE PRECISION, INTENT(OUT):: VRN
! Total relative translational velocity (vector).
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(OUT):: VRT

! Local variables
!---------------------------------------------------------------------//
! Wall tangential velocity at contact point
      DOUBLE PRECISION, DIMENSION(DIMN) :: WALL_TANG_VEL
! Total relative velocity at contact point
      DOUBLE PRECISION, DIMENSION(DIMN) :: VRELTRANS


      VRELTRANS(:) =  DES_VEL_NEW(LL,:) + V_ROT(:)
      VRELTRANS(:) =  VRELTRANS(:) - (STL_TANG_VEL(:,N) + STL_NORMAL_VEL(:,N))

! magnitude of normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! total relative translational slip velocity at the contact point
! Equation (8) in Tsuji et al. 1992
      VRT(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL_WALL


!----------------------------------------------------------------//
! SUBROUTINE: CALC_DEM_THERMO_WITH_WALL_STL
! By: Aaron M.
! Purpose: Compute heat transfer to particles due to boundaries
!----------------------------------------------------------------//

      SUBROUTINE CALC_DEM_THERMO_WITH_WALL_STL

      IMPLICIT NONE

      INTEGER :: LL             ! Loop index for particle
      INTEGER :: CELL_ID        ! Desgrid cell index
      INTEGER :: CELL_COUNT     ! Loop index for facets
      INTEGER :: IJK_FLUID      ! IJK index for fluid cell
      INTEGER :: I_FLUID, J_FLUID, K_FLUID
      INTEGER :: I_Facet, J_Facet, K_Facet, IJK_Facet
      INTEGER :: I1,J1,K1
      INTEGER :: phase_LL       ! Phase index for particle

      INTEGER, PARAMETER :: MAX_CONTACTS = 12
      INTEGER :: iFacet         ! loop index for facets
      INTEGER :: count_Facets   ! counter for number of facets particle contacts
      DOUBLE PRECISION, DIMENSION(3,MAX_CONTACTS) :: NORM_FAC_CONTACT
      LOGICAL :: USE_FACET
      LOGICAL :: DOMAIN_BDRY

      INTEGER :: NF             ! Facet ID
      INTEGER :: AXIS           ! Facet direction
      DOUBLE PRECISION :: RLENS_SQ ! lens radius squared
      DOUBLE PRECISION :: RLENS ! lens radius
      DOUBLE PRECISION :: RPART ! Particle radius

      ! vectors for minimum / maximum possible particle contacts
      DOUBLE PRECISION, DIMENSION(3) :: PARTPOS_MIN, PARTPOS_MAX
      DOUBLE PRECISION, DIMENSION(3) :: POS_TMP ! Temp vector
      DOUBLE PRECISION, DIMENSION(3) :: DIST, NORMAL

      DOUBLE PRECISION, DIMENSION(3) :: CLOSEST_PT

      DOUBLE PRECISION :: line_t  ! Normal projection for part/wall
      DOUBLE PRECISION :: DISTSQ ! Separation from particle and facet (squared)
      DOUBLE PRECISION :: PROJ

      INTEGER :: BC_ID ! BC ID
      INTEGER :: IBC   ! BC Loop Index

      DOUBLE PRECISION :: TWALL
      DOUBLE PRECISION :: K_gas
      DOUBLE PRECISION :: TPART
      DOUBLE PRECISION :: OVERLAP
      DOUBLE PRECISION :: QSWALL, AREA

      LOGICAL, SAVE :: OUTPUT_WARNING = .TRUE.

! Flag to distinguish point or edge intersection
      logical :: point_or_edge_int
      integer :: edge

      STL_Q(:) = ZERO


    !  IF(.NOT.DES_CONTINUUM_COUPLED.or.DES_EXPLICITLY_COUPLED)THEN
    !     CALL PARTICLES_IN_CELL
    !  ENDIF
      DO LL = 1, MAX_PIP
         ! Skip non-existent particles or ghost particles
         IF (.NOT.IS_NORMAL(LL)) CYCLE
         PHASE_LL = PIJK(LL,5)
         IF(.NOT.CALC_COND_DES(PHASE_LL))CYCLE
         ! Get desgrid cell index
         CELL_ID = DG_PIJK(LL)

         ! Skip cells that do not have neighboring facet
         IF (facets_at_dg(CELL_ID)%COUNT <1) CYCLE

         ! Store lens radius of particle
         RLENS = (ONE+FLPC)*DES_RADIUS(LL)
         RLENS_SQ = RLENS*RLENS

         RPART = DES_RADIUS(LL)

         ! Compute max/min particle locations
         PARTPOS_MAX(:) = des_pos_new(LL,:) + RLENS
         PARTPOS_MIN(:) = des_pos_new(LL,:) - RLENS

         ! Get fluid cell
         I_FLUID = PIJK(LL,1)
         J_FLUID = PIJK(LL,2)
         K_FLUID = PIJK(LL,3)
         IJK_FLUID = PIJK(LL,4)
         TPART = DES_T_S(LL)

         ! Sometimes PIJK is not updated every solids timestep and
         ! therefore doesn't give the correct fluid cell that
         ! contains the particles.  If so, determine actual fluid
         ! cell that contains the particle
         IF(.NOT.DES_CONTINUUM_COUPLED.or.DES_EXPLICITLY_COUPLED)THEN
            IF(I_FLUID <= ISTART3 .OR. I_FLUID >= IEND3) THEN
               CALL PIC_SEARCH(I_FLUID, DES_POS_NEW(LL,1), XE,       &
               DIMENSION_I, IMIN2, IMAX2)
            ELSE
               IF((DES_POS_NEW(LL,1) >= XE(I_FLUID-1)) .AND.         &
               (DES_POS_NEW(LL,1) <  XE(I_FLUID))) THEN
               I_FLUID = I_FLUID
               ELSEIF((DES_POS_NEW(LL,1) >= XE(I_FLUID)) .AND.       &
                  (DES_POS_NEW(LL,1) < XE(I_FLUID+1))) THEN
                  I_FLUID = I_FLUID+1
               ELSEIF((DES_POS_NEW(LL,1) >= XE(I_FLUID-2)) .AND.     &
                  (DES_POS_NEW(LL,1) < XE(I_FLUID-1))) THEN
                  I_FLUID = I_FLUID-1
               ELSE
                  CALL PIC_SEARCH(I_FLUID, DES_POS_NEW(LL,1), XE,    &
                  DIMENSION_I, IMIN2, IMAX2)
               ENDIF
            ENDIF !(I)
            IF(J_FLUID <= JSTART3 .OR. J_FLUID >= JEND3) THEN
               CALL PIC_SEARCH(J_FLUID, DES_POS_NEW(LL,2), YN,       &
               DIMENSION_J, JMIN2, JMAX2)
            ELSE
               IF((DES_POS_NEW(LL,2) >= YN(J_FLUID-1)) .AND.         &
                  (DES_POS_NEW(LL,2) <  YN(J_FLUID))) THEN
                  J_FLUID = J_FLUID
               ELSEIF((DES_POS_NEW(LL,2) >= YN(J_FLUID)) .AND.       &
                  (DES_POS_NEW(LL,2) < YN(J_FLUID+1))) THEN
                  J_FLUID = J_FLUID+1
               ELSEIF((DES_POS_NEW(LL,2) >= YN(J_FLUID-2)) .AND.     &
                  (DES_POS_NEW(LL,2) < YN(J_FLUID-1))) THEN
                  J_FLUID = J_FLUID-1
               ELSE
                  CALL PIC_SEARCH(J_FLUID, DES_POS_NEW(LL,2), YN,    &
                  DIMENSION_J, JMIN2, JMAX2)
               ENDIF
            ENDIF !(J)

            IF(NO_K) THEN
               K_FLUID = 1
            ELSE
               IF(K_FLUID <= KSTART3 .OR. K_FLUID >= KEND3) THEN
               CALL PIC_SEARCH(K_FLUID, DES_POS_NEW(LL,3), ZT,         &
                  DIMENSION_K, KMIN2, KMAX2)
               ELSE
                  IF((DES_POS_NEW(LL,3) >= ZT(K_FLUID-1)) .AND.        &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID))) THEN
                     K_FLUID = K_FLUID
                  ELSEIF((DES_POS_NEW(LL,3) >= ZT(K_FLUID)) .AND.      &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID+1))) THEN
                     K_FLUID = K_FLUID+1
                  ELSEIF((DES_POS_NEW(LL,3) >= ZT(K_FLUID-2)) .AND.    &
                     (DES_POS_NEW(LL,3) < ZT(K_FLUID-1))) THEN
                     K_FLUID = K_FLUID-1
                  ELSE
                     CALL PIC_SEARCH(K_FLUID, DES_POS_NEW(LL,3), ZT,   &
                     DIMENSION_K, KMIN2, KMAX2)
                  ENDIF
               ENDIF
            ENDIF !(K)
            IJK_FLUID = PIJK(LL,4)
         ENDIF ! (NOT CONTINUUM_COUPLED OR EXPLICIT)


! Initialize counter for number of facets
         COUNT_FACETS = 0

         ! Loop over potential facets
         DO CELL_COUNT = 1, facets_at_dg(CELL_ID)%count
            ! Get direction (axis) and facet id (NF)
            axis = facets_at_dg(cell_id)%dir(cell_count)
            NF = facets_at_dg(cell_id)%id(cell_count)

            ! Check to see if facet is out of reach of particle
            if (facets_at_dg(cell_id)%min(cell_count) >    &
               partpos_max(axis)) then
               cycle
            endif
            if (facets_at_dg(cell_id)%max(cell_count) <    &
               partpos_min(axis)) then
               cycle
            endif

            ! Compute projection (normal distance) from wall to particle
            line_t = DOT_PRODUCT(VERTEX(1,:,NF) - des_pos_new(LL,:),&
            &        NORM_FACE(:,NF))

            ! If normal exceeds particle radius, particle is not in contact
            if((line_t.lt.-RLENS))CYCLE

            ! Compute closest point on facet
            POS_TMP(:) = DES_POS_NEW(LL,:)
            CALL ClosestPtPointTriangle(NF, POS_TMP, VERTEX(:,:,NF), &
            &    CLOSEST_PT(:),point_or_edge_int, edge)
            ! Compute position vector from particle to closest point on facet
            DIST(:) = CLOSEST_PT(:)-POS_TMP(:)
            DISTSQ = DOT_PRODUCT(DIST,DIST)

            ! Skip particles that are more than lens radius from facet
            IF(DISTSQ .GE. (RLENS_SQ-SMALL_NUMBER))CYCLE

            ! Only do heat transfer for particles that are directly above facet
            ! Heat transfer routines not generalized for edge contacts, but
            ! those should normally yield negligible heat transfer contributions

            ! Normalize distance vector and compare to facet normal
            NORMAL(:)=-DIST(:)/sqrt(DISTSQ)
            PROJ = sqrt(abs(DOT_PRODUCT(NORMAL, NORM_FACE(:,NF))))
            IF(ABS(ONE-PROJ).gt.1.0D-6)CYCLE

            ! Get overlap
            OVERLAP = RPART - SQRT(DISTSQ)

            ! Initialize area for facet (for post-proc. flux)
            AREA = ZERO

            ! Initialize BDRY FLAG
            DOMAIN_BDRY = .FALSE.
            ! Get BC_ID
            BC_ID = BC_ID_STL_FACE(NF)

!JFD: This test (BC_ID.eq.0) may not be needed anynore because BC_ID_STL_FACE
! is set to -1 for default walls (see CONVERT_DEFAULT_WALLS_TO_STL
! in stl_preproc_des_mod.f)

            ! BC_ID is set to 0 in preproc if stl is a domain boundary
            IF(BC_ID.eq.0)then
               I1=I_FLUID
               J1=J_FLUID
               K1=K_FLUID
               DOMAIN_BDRY = .TRUE.

               ! Domain boundary, figure out real boundary ID
               IF(NORM_FACE(1,NF).ge.0.9999)THEN
                  ! WEST face
                  I1=IMIN2
                  IF(DO_K)THEN
                     AREA = DY(J_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DY(J_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(1,NF).le.-0.9999)THEN
                  ! EAST FACE
                  I1=IMAX2
                  IF(DO_K)THEN
                     AREA = DY(J_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DY(J_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(2,NF).ge.0.9999)THEN
                  ! SOUTH FACE
                  J1=JMIN2
                  IF(DO_K)THEN
                     AREA = DX(I_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DX(I_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(2,NF).le.-0.9999)THEN
                  ! NORTH FACE
                  J1=JMAX2
                  IF(DO_K)THEN
                     AREA = DX(I_FLUID)*DZ(K_FLUID)
                  ELSE
                     AREA = DX(I_FLUID)*ZLENGTH
                  ENDIF
               ELSEIF(NORM_FACE(3,NF).ge.0.9999)THEN
                  ! BOTTOM FACE
                  K1=KMIN2
                  AREA = DX(I_FLUID)*DY(J_FLUID)

               ELSEIF(NORM_FACE(3,NF).le.-0.9999)THEN
                  ! TOP FACE
                  K1=KMAX2
                  AREA = DX(I_FLUID)*DY(J_FLUID)
               ELSE
                  WRITE( *,*)'PROBLEM, COULD NOT FIND DOMAIN BOUNDARY'
                  WRITE(*,*)' In calc_thermo_des_wall_stl'
                  call mfix_exit(1)

               ENDIF


               ! Loop through defined BCs to see which one particle neighbors
               DO IBC = 1, DIMENSION_BC
                  IF(.NOT.BC_DEFINED(IBC))CYCLE
                  IF (I1.ge.BC_I_W(IBC).and.I1.le.BC_I_E(IBC).and.&
                      J1.ge.BC_J_S(IBC).and.J1.le.BC_J_N(IBC).and.&
                      K1.ge.BC_K_B(IBC).and.K1.le.BC_K_T(IBC))THEN
                      BC_ID = IBC
                      exit
                  ENDIF
               ENDDO
               IF(BC_ID.eq.0)then
                  IF(OUTPUT_WARNING)THEN
1111                  FORMAT("Warning: Could not find BC."/"Check input file to make sure domain boundaries are defined.",/ &
                          "DES_POS_NEW: ", 3(F12.5, 3X),/ &
                          "I,J,K: ", I4,I4,I4,          / &
                          "CLOSEST_PT: ", 3(F12.5, 3X), / &
                          "NORM_FACE: ",  3(F12.5, 3X), / &
                          'Suppressing further warnings.')
                      WRITE(ERR_MSG, 1111) (DES_POS_NEW(LL,IBC),IBC=1,3), &
                          I1, J1, K1,                                     &
                          (CLOSEST_PT(IBC),IBC=1,3),                      &
                          (NORM_FACE(IBC,NF),IBC=1,3)
                      OUTPUT_WARNING = .FALSE.
                      CALL LOG_WARNING()
                  ENDIF
                  CYCLE  !
               ENDIF

            ENDIF !Domain Boundary (facet ID was 0)

! JFD: At this point, if BC_ID==-1, this means the particle is
! contacting a default wall. Default wall are adiabatic, so we can cycle
! and move to the next facet.
            IF(BC_ID==-1) CYCLE ! Adiabatic default wall

            IF (BC_TYPE_ENUM(BC_ID) == NO_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == FREE_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == PAR_SLIP_WALL .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_NSW .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_FSW .OR. &
               BC_TYPE_ENUM(BC_ID) == CG_PSW) THEN

               ! CHECK TO MAKE SURE FACET IS UNIQUE
               USE_FACET=.TRUE.
               DO IFACET=1,count_facets
                  ! DO CHECK BY ENSURING NORMAL VECTOR IS NEARLY PARALLEL
                  PROJ = sqrt(abs(DOT_PRODUCT(NORMAL, NORM_FAC_CONTACT(:,IFACET))))
                  IF(ABS(ONE-PROJ).lt.1.0D-6)THEN
                     USE_FACET=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF(.NOT.USE_FACET)CYCLE

               ! FACET IS UNIQUE
               count_facets=count_facets+1
               NORM_FAC_CONTACT(:,count_facets)=NORMAL(:)

! Do heat transfer
               ! GET WALL TEMPERATURE
               TWALL = BC_TW_S(BC_ID,phase_LL)

               ! GET GAS THERMAL CONDUCTIVITY
               if(k_g0.eq.UNDEFINED)then
                  ! Compute gas conductivity as is done in calc_k_g
                  ! But use average of particle and wall temperature to be gas temperature

                  K_Gas = 6.02D-5*SQRT(HALF*(TWALL+TPART)/300.D0) ! cal/(s.cm.K)
                  ! 1 cal = 4.183925D0 J
                  IF (UNITS == 'SI') K_Gas = 418.3925D0*K_Gas !J/s.m.K
               else
                  K_Gas=k_g0
               endif
               IF(TWALL.eq.UNDEFINED)CYCLE
               QSWALL = DES_CONDUCTION_WALL(LL, OVERLAP,K_s0(phase_LL), &
               &        K_s0(phase_LL),K_Gas,TWALL, TPART, RPART, &
               &        RLENS, phase_LL)


               Q_Source(LL) = Q_Source(LL)+QSWALL

               IF(STL_AREA(NF)>ZERO) STL_Q(NF) = STL_Q(NF) + QSWALL/STL_AREA(NF)

               ! BELOW CODE IS ONLY NECESSARY FOR OUTPUTTING
               ! DATA.  Need to know fluid cell that contact
               ! point resides in so that wall flux can be
               ! output correctly.

               I_FACET = I_FLUID
               J_FACET = J_FLUID
               K_FACET = K_FLUID

               ! This checks to see if the contact point was NOT on a domain boundary

               IF(.NOT.DOMAIN_BDRY)THEN
                  IF(CLOSEST_PT(1) >= XE(I_FLUID))THEN
                     I_FACET = MIN(IMAX1, I_FLUID+1)
                  ELSEIF(CLOSEST_PT(1) < XE(I_FLUID-1))THEN
                     I_FACET = MAX(IMIN1, I_FLUID-1)
                  ENDIF

                  IF(CLOSEST_PT(2) >= YN(J_FLUID))THEN
                     J_FACET = MIN(JMAX1, J_FLUID+1)
                  ELSEIF(CLOSEST_PT(2) < YN(J_FLUID-1))THEN
                     J_FACET = MAX(JMIN1, J_FLUID-1)
                  ENDIF
                  IF(DO_K)THEN
                     IF(CLOSEST_PT(3) >= ZT(K_FLUID))THEN
                        K_FACET = MIN(KMAX1, K_FLUID+1)
                     ELSEIF(CLOSEST_PT(3) < ZT(K_FLUID-1))THEN
                        K_FACET = MAX(KMIN1, K_FLUID-1)
                     ENDIF
                  ENDIF
                  IJK_FACET=funijk(I_facet,J_facet,K_facet)
                  if (cartesian_grid) AREA=AREA_CUT(IJK_FACET)

               ! AREA is left undefined if contact was with cut-cell surface
                  IF (AREA.eq.ZERO)then
                     I_FACET = I_FLUID
                     J_FACET = J_FLUID
                     K_FACET = K_FLUID
                     IJK_FACET=funijk(I_facet,J_facet,K_facet)
                     IF(NORM_FACE(1,NF).ge.0.9999)THEN
                     ! WEST face
                        IF(DO_K)THEN
                           AREA = DY(J_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DY(J_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(1,NF).le.-0.9999)THEN
                     ! EAST FACE
                        IF(DO_K)THEN
                           AREA = DY(J_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DY(J_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(2,NF).ge.0.9999)THEN
                     ! SOUTH FACE
                        IF(DO_K)THEN
                           AREA = DX(I_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DX(I_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(2,NF).le.-0.9999)THEN
                     ! NORTH FACE
                        IF(DO_K)THEN
                           AREA = DX(I_FACET)*DZ(K_FACET)
                        ELSE
                           AREA = DX(I_FACET)*ZLENGTH
                        ENDIF
                     ELSEIF(NORM_FACE(3,NF).ge.0.9999)THEN
                     ! BOTTOM FACE
                        AREA = DX(I_FACET)*DY(J_FACET)
                     ELSEIF(NORM_FACE(3,NF).le.-0.9999)THEN
                     ! TOP FACE
                        AREA = DX(I_FACET)*DY(J_FACET)
                     ENDIF
                  ENDIF ! Area==0 (because cut-cell facet exists in non cut-cell)
               ENDIF ! NOT DOMAIN_BDRY
               IJK_FACET = FUNIJK(I_FACET,J_FACET,K_FACET)
               ! AM error check
               ! todo use proper error handling
               IF(cartesian_grid) THEN
                  IF(.NOT.FLUID_AT(IJK_FACET).AND. &
                     .NOT.BLOCKED_CELL_AT(IJK_FACET)) THEN
                     write(*,*)'ERROR: Cell containing facet is not a fluid', &
                     'or blocked cell'
                     write(*,*)FLUID_AT(IJK_FACET), BLOCKED_CELL_AT(IJK_FACET)
                     write(*,*)'PART POS',(DES_POS_NEW(LL,IBC),IBC=1,3)
                     write(*,*)'FACET NORM',(NORM_FACE(IBC,NF),IBC=1,3)
                     write(*,*)'BC_ID', BC_ID
                     write(*,*)'I,J,K (Facet)', I_FACET,J_FACET,K_FACET
                     call mfix_exit(1)
                  ENDIF
               ENDIF
               IF(FLUID_AT(IJK_FACET).AND.AREA>ZERO)THEN
                  DES_QW_Cond(IJK_FACET,phase_LL) = &
                     DES_QW_Cond(IJK_FACET, phase_LL) + QSWALL/AREA
               ENDIF

            ENDIF ! WALL BDRY
         ENDDO                  ! CELL_COUNT (facets)
      ENDDO  ! LL


! Sum across ranks
      call global_all_sum(stl_q)

      RETURN

   END SUBROUTINE CALC_DEM_THERMO_WITH_WALL_STL



!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: get_wall_tang_velocity                            !
!  Author: Jeff Dietiker                            Date: 18-Apr-2024  !
!                                                                      !
!  Purpose: Calculate the tangential velocity of a moving wall.        !
!           The velocity is computed at the point where a particle     !
!           and a wall are in contact (closest_pt).                    !
!           The wall can:                                              !
!             - Translate at velocity bc_wall_vel(:)                   !
!             - Rotate wrt center point bc_wall_rot_center(:)          !
!               with angular velocity bc_wall_omega(:)                 !
!                                                                      !
!----------------------------------------------------------------------!

      subroutine get_wall_tang_velocity(bc_id, closest_pt, wall_tang_vel)
      use keyframe_mod, only: omega_conv_fact


      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Boundary condition ID
        integer, intent(in)                             :: bc_id
! Point closest to the wall
        double precision, dimension(dimn), intent(in)   :: closest_pt
! Wall tangential velocity at contact point
        double precision, dimension(dimn), intent(out)  :: wall_tang_vel

! Local variables
!---------------------------------------------------------------------//
! Line joining center of rotation to particle center
        double precision, dimension(dimn)               :: OP
! Get BCID to index since BCIDs are not contiguous
        integer                                         :: ind
! Angular velocity (rad/sec)
        double precision, dimension(3) :: omega

        wall_tang_vel(:) = ZERO

        if(bc_id<=0) return

        omega = bc_wall_omega(bc_id, :) * omega_conv_fact
        OP = closest_pt - bc_wall_rot_center(bc_id, :)
        wall_tang_vel(:) = bc_wall_vel(bc_id, :) + cross(omega, OP)

      end subroutine get_wall_tang_velocity

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: gsp_wall_contact_implicit                               !
!  Author: Hang Zhou                                Date: April-2024   !
!  Revised: Renjie Ke                               Date: NOV-1-2024   !
!                                                                      !
!  Purpose: Calculate the contact point, normal direction and overlap  !
!           of a particle to wall. Note that those variables are       !
!           summed up here and will be averaged based on number of     !
!           contact.                                                   !
!----------------------------------------------------------------------!

      subroutine gsp_wall_contact_implicit(ll, NF, closest_pt, overlap_n, num_contact, normal, contact_global)

         use discretelement, only: cspart_data, glueEX, glueEY, glueEZ, pijk
         use physprop, only: d_p0
         use stl, only: norm_face

         integer, intent(in) :: ll
         integer, intent(in) :: NF
         double precision, intent(in) :: closest_pt(:)
         double precision, intent(out) :: overlap_n
         double precision, intent(out) :: normal(:)
         integer, intent(out) :: num_contact
         double precision, intent(out) :: contact_global(:)

         integer :: index_ll, pid, row_n
         double precision :: contact_ll(3), norm_w(3), dist(3), sphere_pos(3), distsq
         double precision :: sphere_radius, dist_n, overlap_ll, norm_w_sq, xloc, yloc, zloc
         double precision :: real_b_ll
         ! Find the closest point based on the location of the component sphere
         logical :: point_or_edge_int, moving_wall
         integer :: edge
         double precision :: closest_pt_sphere(3)
         double precision :: norm_w_tmp(3), norm_w_sq_tmp

         ! phase id
         pid = pijk(ll,5)

         ! locate components in cspart_data(:,:,pid)
         row_n = size(cspart_data(:,:,pid),1)

         ! initialization
         overlap_n = 0.0d0
         num_contact = 0
         normal(:) = (/0.d0,0.d0,0.d0/)
         contact_global(:) = (/0.d0,0.d0,0.d0/)
         closest_pt_sphere(:) = (/0.d0,0.d0,0.d0/)

         ! facet wall normal point to the fluid domain side
         norm_w_tmp(:) = norm_face(:,nf)
         norm_w_sq_tmp = dot_product(norm_w_tmp,norm_w_tmp)
         norm_w_tmp(:) = norm_w_tmp(:)/sqrt(norm_w_sq_tmp)

         real_b_ll = 2.0 * des_radius(ll)

         do index_ll = 1, row_n
            ! use component diameter to determine if this location represents a component sphere
            if(cspart_data(index_ll,4,pid) == 0.0) exit
            ! cycle if this component is not a surface component
            if(cspart_data(index_ll,5,pid) == 0.0) cycle

            ! cspart_data is filled with relative values, scale back to abs
            xloc = cspart_data(index_ll,1,pid) * real_b_ll
            yloc = cspart_data(index_ll,2,pid) * real_b_ll
            zloc = cspart_data(index_ll,3,pid) * real_b_ll

            ! calculate cmponent locations
            sphere_pos(:) = des_pos_new(ll,:) + xloc*glueEX(ll,:) &
                  + yloc*glueEY(ll,:) + zloc*glueEZ(ll,:)

            CALL ClosestPtPointTriangle(NF, sphere_pos, &
               VERTEX(:,:,NF), closest_pt_sphere(:),point_or_edge_int, edge)

            dist(:) =  sphere_pos(:) - closest_pt_sphere(:)
            distsq = dot_product(dist,dist)
            norm_w(:) = dist(:)/sqrt(distsq)
            ! the direction of norm_w now depends on whether the sphere is outside or inside
            ! norm_w is obtained from dist, so dist_n is always > 0
            dist_n = dot_product(dist, norm_w)
            dist_tmp = dot_product(dist, norm_w_tmp)

            sphere_radius = cspart_data(index_ll,4,pid) * real_b_ll / 2.0

            ! @renjieke, only allow the contact when component sphere is still within the domain
            ! so a large enough contact stiffness might have to use for that
            if (dist_tmp > 0.0 .and. dist_n < sphere_radius) then
               num_contact = num_contact + 1
               overlap_ll = sphere_radius - dist_n
               ! from sphere center --> perpendicular line to wall --> insecting point
               contact_ll(:) = sphere_pos(:) - (sphere_radius-overlap_ll)*norm_w
               overlap_n = overlap_n + overlap_ll
               normal = normal + norm_w
               contact_global = contact_global + contact_ll
            endif
            ! elseif(dist_tmp < 0.0) then
            !    num_contact = num_contact + 1
            !    overlap_ll = sphere_radius - dist_n
            !    norm_w = -norm_w
            !    ! from sphere center --> perpendicular line to wall --> insecting point
            !    contact_ll(:) = sphere_pos(:) - (sphere_radius-overlap_ll/2.0)*norm_w
            !    overlap_n = overlap_n + overlap_ll
            !    normal = normal + norm_w
            !    contact_global = contact_global + contact_ll
            ! elseif(dist_tmp == 0.0) then
            ! ! @renjiek in gsp implicit, actually the component sphere can be right at the triangle,
            ! ! so sphere_pos is exactly the cloestest_point
            !    num_contact = num_contact + 1
            !    overlap_n = overlap_n + sphere_radius
            !    normal = normal + norm_w_tmp
            !    contact_global = contact_global + sphere_pos(:)

         enddo
         return
      end subroutine gsp_wall_contact_implicit

END MODULE CALC_COLLISION_WALL
