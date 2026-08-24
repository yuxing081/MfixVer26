#include "error.inc"

MODULE CALC_FORCE_DEM_MOD

   USE CFRELVEL_MOD, ONLY: CFRELVEL
   USE error_manager

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_FORCE_DEM                                          !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Calculate contact force and torque on particle from        !
!           particle-particle and particle-wall collisions. Treats     !
!           wall interaction also as a two-particle interaction but    !
!           accounting for the wall properties                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CALC_FORCE_DEM

! Modules
!---------------------------------------------------------------------//
      USE calc_collision_wall
      use compar, only: pe_io, mype, numPEs
      USE constant, ONLY: Pi
      use desmpi_wrapper, only: DES_MPI_STOP,des_mpi_barrier
      USE des_thermo
      USE des_thermo_cond
      USE discretelement
      USE geometry, only: cyclic_x, cyclic_y, cyclic_z
      USE geometry, only: xlength, ylength, zlength
      USE mpi_unpack_des, ONLY: JKR_DIAG_PAIR_EVENT
      USE resize
      USE run
      use param1, only: half, one, small_number, zero, undefined
      use mpi_utility, only: global_all_sum

      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0 ! what is excess overlap?
! particle no. indices
      INTEGER :: I, LL, cc, CC_START, CC_END
! the overlap occurring between particle-particle or particle-wall
! collision in the normal direction
      DOUBLE PRECISION :: OVERLAP_N, OVERLAP_T(3)
! square root of the overlap
      DOUBLE PRECISION :: SQRT_OVERLAP
! distance vector between two particle centers or between a particle
! center and wall when the two surfaces are just at contact (i.e. no
! overlap)
      DOUBLE PRECISION :: R_LM,DIST_CI,DIST_CL
! the normal and tangential components of the translational relative
! velocity
      DOUBLE PRECISION :: V_REL_TRANS_NORM, rad
! distance vector between two particle centers or between a particle
! center and wall at current and previous time steps
      DOUBLE PRECISION :: DIST(3), NORMAL(3), DIST_MAG, POS(3)
! tangent to the plane of contact at current time step
      DOUBLE PRECISION :: VREL_T(3)
! normal and tangential forces
      DOUBLE PRECISION :: FN(3), FT(3)
! temporary storage of force
      DOUBLE PRECISION :: FC_TMP(3)
! temporary storage of force for torque
      DOUBLE PRECISION :: TOW_FORCE(3)
! temporary storage of torque
      DOUBLE PRECISION :: TOW_TMP(3,2)
! temporary storage of conduction/radiation
      DOUBLE PRECISION :: QQ_TMP

! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER :: PHASEI, PHASELL
! local values used spring constants and damping coefficients
      DOUBLE PRECISION :: ETAN_DES, ETAT_DES
      DOUBLE PRECISION :: KN_DES, KT_DES
! local values used for calculating cohesive forces
      DOUBLE PRECISION :: FORCE_COH, EQ_RADIUS, DistApart

      LOGICAL, PARAMETER :: report_excess_overlap = .FALSE.

      DOUBLE PRECISION :: FNMD, FTMD, MAG_OVERLAP_T, TANGENT(3)
      integer :: start_index, end_index, do_index

!{SBAN
! Effective particle radius and mass to recalculate collision parameters
      DOUBLE PRECISION :: rad_eff, m_eff
!}

! Glued particle DEM model
! the glued particle id which sphere LL and sphere I belong to, LL-I is a cntact sphere pair
      INTEGER :: gidI, gidLL, gid, io_err, ppdim1, ppdim2
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ppforce !(NGluedParticles, 3)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pptorque !(NGluedParticles, 3)
! contact center, contact center to glued particle LL belonging to, contact center to glued particle I belonging to
! translational relative velocity
      DOUBLE PRECISION :: Contact_PT(3),dist2centerVECLL(3),dist2centerVECI(3),VRELTRANS(3)

! Force chain data
      INTEGER :: FCHAINC_OLD           ! old value of FCHAINC
      INTEGER :: old_size,new_size ! New force chain array size when resizing (growing) arrays
      DOUBLE PRECISION :: FC_MAG_MAX, MAG_ORIENT
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Coh_force

      INTEGER :: FC_MAX_INDEX         ! FC index where max occurs (within LL-I pair)

! Rolling friction
      DOUBLE PRECISION :: RFTOW(3), omega_mag

! Spherical JKR contact and Chen rolling-resistance temporaries
      DOUBLE PRECISION :: A_OVER_A0, A_JKR, JKR_NORM(DIMN)
      DOUBLE PRECISION :: JKR_A0_LOCAL, JKR_DELTAC_LOCAL
      DOUBLE PRECISION :: JKR_FC_LOCAL, JKR_GAMMA_LOCAL
      DOUBLE PRECISION :: JKR_E_EFF_LOCAL
      DOUBLE PRECISION :: V_ROLL(3), XI_MAG, XI_CRIT, KR_JKR
      DOUBLE PRECISION :: MR_TORQUE(3)

! BONDED test parameters here
      INTEGER :: nbid
      DOUBLE PRECISION :: F1(3), F2(3), T1(3), T2(3)
!------------------calculations start here!--------------------------------------------------->

! Initialize cohesive forces
      IF(USE_COHESION) PostCohesive(:) = ZERO

! DES_USR_VAR(1,:) is a production particle property for JKR surface
! energy. Non-positive values select and are normalized to the phase default.
      IF(USE_JKR .AND. DES_USR_VAR_SIZE > 0) THEN
         DO LL=1,MAX_PIP
            IF(.NOT.IS_NORMAL(LL)) CYCLE
            IF(DES_USR_VAR(1,LL) <= ZERO) &
               DES_USR_VAR(1,LL) = DES_SURFACE_ENERGY(PIJK(LL,5))
         ENDDO
      ENDIF

! calculate particle wall contact
      CALL CALC_DEM_FORCE_WITH_WALL_STL

! at the beginning each time step, reset those arrays to zero
      IF (GSP_EXPLICIT) THEN
         IF(RIGID_GSP) THEN
            allocate(ppforce(NGluedParticles,3)); ppforce(:,:) = ZERO
            allocate(pptorque(NGluedParticles,3)); pptorque(:,:) = ZERO
         ENDIF
         des_coll_force(:,:) = ZERO
      ENDIF

      if (optflag1.eq.1) then
        start_index = 1
        end_index = MAX_PIP_EXIST
      else
        start_index = 1
        end_index = MAX_PIP
      endif


      IF(USE_COHESION) THEN
         allocate(Coh_force(MAX_PIP,3))
         Coh_force(:,:) = ZERO
      ENDIF

! Check particle LL neighbor contacts
!---------------------------------------------------------------------//

!$omp parallel default(none) private(pos,rad,cc,cc_start,cc_end,ll,i,  &
!$omp    overlap_n,vrel_t,v_rel_trans_norm,sqrt_overlap,dist,r_lm,     &
!$omp    kn_des,kt_des,phasell,phasei,etan_des,rad_eff,m_eff,          &
!$omp    gidLL,gidi,contact_pt,dist2centervecll,dist2centerveci,       &
!$omp    VRELTRANS, nbid, F1, F2, T1, T2,                              &
!$omp    etat_des,fn,ft,overlap_t,tangent,mag_overlap_t,               &
!$omp    eq_radius,distapart,force_coh,dist_mag,NORMAL,ftmd,fnmd,      &
!$omp    dist_cl, dist_ci, fc_tmp, tow_tmp, tow_force, qq_tmp,         &
!$omp    a_over_a0,a_jkr,jkr_norm,jkr_a0_local,jkr_deltac_local,      &
!$omp    jkr_fc_local,jkr_gamma_local,jkr_e_eff_local,v_roll,         &
!$omp    xi_mag,xi_crit,kr_jkr,mr_torque,                             &
!$omp    rftow,omega_mag )                                             &
!$omp shared(max_pip,neighbors,neighbor_index,des_pos_new,des_radius,  &
!$omp    omega_new, mew_r, numpes,                                     &
!$omp    des_coll_model_enum,kn,kt,pft_neighbor,pijk,hert_kn,hert_kt,  &
!$omp    des_etan,des_etat,mew,use_cohesion, calc_cond_des, dtsolid,   &
!$omp    van_der_waals,vdw_outer_cutoff,vdw_inner_cutoff,              &
!$omp    hamaker_constant,asperities,surface_energy, optflag1,         &
!$omp    start_index, end_index,rk_global_cnt_debug,                   &
!$omp    des_usr_var,GluedSphereDEM,gluepos,des_vel_new,glueomg,       &
!$omp    gsp_implicit, gsp_explicit, gp_neighsid, global2loc_id,       &
!$omp    flexible_gsp, rigid_gsp,                                      &
!$omp    ppforce,pptorque,gid_list, mype,                              &
!$omp    tow, fc, energy_eq, grav_mag, postcohesive, pmass, q_source,  &
!$omp    write_force_chain, des_coll_force, iglobal_id,particle_state, &
!$omp    old_size,new_size,io_err,                                     &
!$omp    fchain_fcmax,fc_mag_max,fc_max_index, mag_orient,coh_force,   &
!$omp    FCHAINC,FCHAINC_OLD,                                          &
!$omp    FCHAIN_MIDPOINT,FCHAIN_CONTACT_POINT,                         &
!$omp    FCHAIN_ORIENT, FCHAIN_FN,FCHAIN_FT,FCHAIN_FC,                 &
!$omp    FCHAIN_LENGTH,FCHAIN_FN_MAG,FCHAIN_LOCAL_ID1,FCHAIN_LOCAL_ID2,&
!$omp    FCHAIN_GLOBAL_ID1,FCHAIN_GLOBAL_ID2, FCHAIN_OVERLAP) &
!$omp shared(use_jkr,pfr_neighbor,jkr_contact_state,jkr_alpha_safety, &
!$omp    hert_kn_jkr,hert_kt_jkr,des_etan_jkr,des_etat_jkr,           &
!$omp    jkr_thetar_crit)                                             &
!$omp shared(max_pip_exist)

! Force chain initialization
      IF(WRITE_FORCE_CHAIN) THEN
         FCHAINC_OLD = min(FCHAINC, size(FCHAIN_MIDPOINT,1))
         FCHAINC = 0
         FCHAIN_MIDPOINT(:,:) = UNDEFINED
         FCHAIN_CONTACT_POINT(:,:) = UNDEFINED
         FCHAIN_ORIENT(:,:) = UNDEFINED
         FCHAIN_FCMAX(:) = .FALSE.
      ENDIF

!$omp do reduction (+: FCHAINC)
      DO do_index = start_index, end_index
         if (optflag1.eq.1) then
            ll = do_index
         else
            ll = do_index
            IF(IS_NONEXISTENT(LL)) CYCLE
         endif

         IF (gsp_explicit .and. energy_eq .and. numpes > 1) THEN
            global2loc_id(iglobal_id(LL)) = LL
         ENDIF

         IF(WRITE_FORCE_CHAIN) THEN
            FC_MAG_MAX = ZERO
            FC_MAX_INDEX = -1
         ENDIF

         POS = DES_POS_NEW(LL,:)
         rad = DES_RADIUS(LL)

! read the gid of the current sphere LL!
         IF (GSP_EXPLICIT) gidLL = gid_list(LL)
!-------------------------------------------------------------------------------------
! NEIGHBORS: store the neighbor component spheres in order
! NEIGHBOR_INDEX: the range index of neighbor component sphere
!-------------------------------------------------------------------------------------
         CC_START = 1
         IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
         CC_END   = NEIGHBOR_INDEX(LL)
         cc_loop: DO CC = CC_START, CC_END-1
! read the gid of the current sphere I!

            I  = NEIGHBORS(CC)
            ! only when those two are from the same GSP, then need to check if calc the bond force
            IF(FLEXIBLE_GSP) THEN
               IF(gid_list(LL) .eq. gid_list(I)) THEN
                  DO nbid = 1, 6
                     ! neighbors now should include the component in the same GSP
                     IF(iglobal_id(I) == gp_neighsid(LL,nbid)) THEN
                        CALL CALC_VIRTUAL_BOND_FORCE(LL, I, F1, F2, T1, T2, nbid)
                        ! FC may contain forces comeing from wall contact
!$omp atomic
                        FC(LL,1) = FC(LL,1) + F1(1)
!$omp atomic
                        FC(LL,2) = FC(LL,2) + F1(2)
!$omp atomic
                        FC(LL,3) = FC(LL,3) + F1(3)
!$omp atomic
                        FC(I,1) = FC(I,1) + F2(1)
!$omp atomic
                        FC(I,2) = FC(I,2) + F2(2)
!$omp atomic
                        FC(I,3) = FC(I,3) + F2(3)
!$omp atomic
                        TOW(LL,1) = TOW(LL,1) + T1(1)
!$omp atomic
                        TOW(LL,2) = TOW(LL,2) + T1(2)
!$omp atomic
                        TOW(LL,3) = TOW(LL,3) + T1(3)
!$omp atomic
                        TOW(I,1) = TOW(I,1) + T2(1)
!$omp atomic
                        TOW(I,2) = TOW(I,2) + T2(2)
!$omp atomic
                        TOW(I,3) = TOW(I,3) + T2(3)
                        CYCLE cc_loop
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            IF(RIGID_GSP) THEN
               gidI = gid_list(I)
               if(gidI == gidLL) CYCLE
            ENDIF

            IF(IS_NONEXISTENT(I)) CYCLE

! distance between LL and I center just about to contact, no overlap
            R_LM = rad + DES_RADIUS(I)

! distance between current and previous time steps
            DIST(:) = DES_POS_NEW(I,:) - POS(:)
            DIST_MAG = dot_product(DIST,DIST)

! reset temporary force of collision
            FC_TMP(:) = ZERO

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
               IF(DIST_MAG < (R_LM+VDW_OUTER_CUTOFF)**2) THEN
                  EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) /     &
                    (DES_RADIUS(LL)+DES_RADIUS(I))
                  IF(DIST_MAG > (VDW_INNER_CUTOFF+R_LM)**2) THEN
                     DistApart = (SQRT(DIST_MAG)-R_LM)
                     FORCE_COH = HAMAKER_CONSTANT * EQ_RADIUS /           &
                        (12d0*DistApart**2) * (Asperities/(Asperities+    &
                        EQ_RADIUS) + ONE/(ONE+Asperities/DistApart)**2 )
                  ELSE
                     FORCE_COH = 2d0 * PI * SURFACE_ENERGY * EQ_RADIUS *  &
                       (Asperities/(Asperities+EQ_RADIUS) + ONE/          &
                       (ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                  ENDIF
                  FC_TMP(:) = DIST(:)*FORCE_COH/SQRT(DIST_MAG)
                  Coh_Force(LL,:) = Coh_Force(LL,:) + FC_TMP(:)
                  Coh_Force(I,:)  = Coh_Force(I ,:) - FC_TMP(:)
                  TOW_TMP(:,:) = ZERO

               ENDIF
            ENDIF

            IF (ENERGY_EQ) THEN
               IF(gsp_implicit) THEN
                  ! No action needed, project2 method does not support for now
                  CONTINUE
               ELSE
! Calculate conduction and radiation for thermodynamic neighbors
                  IF(CALC_COND_DES(PIJK(LL,5))) THEN
                     QQ_TMP = DES_CONDUCTION(LL, I, sqrt(DIST_MAG), PIJK(LL,5), PIJK(LL,4))
!$omp atomic
                     Q_Source(LL) = Q_Source(LL) + QQ_TMP
!$omp atomic
                     Q_Source(I) = Q_Source(I) - QQ_TMP
                  ENDIF
               ENDIF
            ENDIF

! Compute actual-particle JKR parameters before applying the contact cutoff.
            phaseLL = PIJK(LL,5)
            phaseI = PIJK(I,5)
            rad_eff = (rad*DES_RADIUS(I))/(rad+DES_RADIUS(I))
            JKR_FC_LOCAL = ZERO
            JKR_A0_LOCAL = ZERO
            JKR_DELTAC_LOCAL = ZERO
            IF(USE_JKR) THEN
               JKR_GAMMA_LOCAL = SQRT(GET_DES_SURFACE_ENERGY(LL)* &
                  GET_DES_SURFACE_ENERGY(I))
               JKR_E_EFF_LOCAL = 0.75D0*HERT_KN_JKR(phaseLL,phaseI)
               JKR_A0_LOCAL = (9.D0*PI*JKR_GAMMA_LOCAL*rad_eff**2/ &
                  JKR_E_EFF_LOCAL)**(ONE/3.D0)
               JKR_DELTAC_LOCAL = JKR_A0_LOCAL**2 / &
                  (2.D0*6.D0**(ONE/3.D0)*rad_eff)
               JKR_FC_LOCAL = 3.D0*PI*rad_eff*JKR_GAMMA_LOCAL
            ENDIF

            IF(USE_JKR .AND. JKR_FC_LOCAL > ZERO) THEN
               IF(DIST_MAG > (R_LM + JKR_ALPHA_SAFETY* &
                  JKR_DELTAC_LOCAL)**2) THEN
                  PFT_NEIGHBOR(:,CC) = ZERO
                  PFR_NEIGHBOR(:,CC) = ZERO
                  JKR_CONTACT_STATE(CC) = .FALSE.
                  CYCLE
               ENDIF
            ELSE

            IF(DIST_MAG > (R_LM - SMALL_NUMBER)**2) THEN
               PFT_NEIGHBOR(:,CC) = 0.0
               CYCLE
            ENDIF
            ENDIF

            IF(DIST_MAG == 0) THEN
               WRITE(*,8550) LL, I
               ERROR STOP "division by zero"
 8550 FORMAT('distance between particles is zero:',2(2x,I10))
            ENDIF

            DIST_MAG = SQRT(DIST_MAG)
            NORMAL(:)= DIST(:)/DIST_MAG
! Calculate the normal overlap
            OVERLAP_N = R_LM-DIST_MAG

            IF(REPORT_EXCESS_OVERLAP) CALL PRINT_EXCESS_OVERLAP

            IF (RIGID_GSP) THEN
! translational relative velocity for glued particles
               Contact_PT = POS(:) +  NORMAL(:)*(rad-OVERLAP_N/2)
               dist2centerVECLL(:) = 0.0
               dist2centerVECI(:) = 0.0
               CALL GSP_MIN_IMAGE_VECTOR(Contact_PT, gluePos(gidLL,:), dist2centerVECLL)
               CALL GSP_MIN_IMAGE_VECTOR(Contact_PT, gluePos(gidI,:),  dist2centerVECI)
               VRELTRANS(:) = (DES_VEL_NEW(LL,:) - DES_VEL_NEW(I,:))
               VRELTRANS(:) =  VRELTRANS(:) + CROSS(glueOmg(gidLL,:),dist2centerVECLL)
               VRELTRANS(:) =  VRELTRANS(:) - CROSS(glueOmg(gidI,:),dist2centerVECI)
               V_REL_TRANS_NORM = DOT_PRODUCT(VRELTRANS,NORMAL)
               VREL_T(:) =  VRELTRANS(:) - V_REL_TRANS_NORM*NORMAL(:)
            ELSEIF (GSP_IMPLICIT) THEN
               CALL gsp_gsp_collison(LL, I, CC, V_REL_TRANS_NORM, VREL_T, NORMAL, OVERLAP_N, &
                                    dist2centerVECI, dist2centerVECLL)
               IF(OVERLAP_N == 0.0) CYCLE
            ELSE
! Calculate the components of translational relative velocity for a
! contacting particle pair and the tangent to the plane of contact
! flexible_gsp is also go into here
               CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, VREL_T, NORMAL(:), DIST_MAG)
            ENDIF

! A tensile JKR neck remains active only if this contact was established
! previously at nonnegative overlap.
            IF(USE_JKR .AND. JKR_FC_LOCAL > ZERO) THEN
               IF(OVERLAP_N >= ZERO) THEN
                  JKR_CONTACT_STATE(CC) = .TRUE.
               ELSEIF(.NOT.JKR_CONTACT_STATE(CC)) THEN
                  PFT_NEIGHBOR(:,CC) = ZERO
                  PFR_NEIGHBOR(:,CC) = ZERO
                  JKR_CONTACT_STATE(CC) = .FALSE.
                  CYCLE
               ENDIF
               A_OVER_A0 = NEWTON_RAPHSON_METHOD(OVERLAP_N / &
                  (6.D0**(ONE/3.D0)*JKR_DELTAC_LOCAL))
            ENDIF
            phaseLL = PIJK(LL,5)
            phaseI = PIJK(I,5)
!{SBAN
            rad_eff = (rad * DES_RADIUS(I)) / (rad + DES_RADIUS(I))
            m_eff = (PMASS(LL) * PMASS(I)) / (PMASS(LL) + PMASS(I))
! Hertz spring-dashpot contact model
            IF(USE_JKR .AND. JKR_FC_LOCAL > ZERO) THEN
               A_JKR = A_OVER_A0*JKR_A0_LOCAL
               KN_DES = HERT_KN_JKR(phaseLL,phaseI)*A_JKR
               KT_DES = HERT_KT_JKR(phaseLL,phaseI)*A_JKR
               ETAN_DES = DES_ETAN_JKR(phaseLL,phaseI)*SQRT(A_JKR)
               ETAT_DES = DES_ETAT_JKR(phaseLL,phaseI)*SQRT(A_JKR)
               JKR_NORM(:) = -4.D0*JKR_FC_LOCAL* &
                  (A_OVER_A0**3-A_OVER_A0**1.5D0)*NORMAL(:)
            ELSE
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)
               KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap * sqrt(rad_eff)
               KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap * sqrt(rad_eff)
               sqrt_overlap = SQRT(sqrt_overlap)
               ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap * sqrt(hert_kn(phaseLL,phaseI)*sqrt(rad_eff)*m_eff)
               ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap * sqrt(hert_kt(phaseLL,phaseI)*sqrt(rad_eff)*m_eff)
! Linear spring-dashpot contact model
            ELSE
               KN_DES = KN
               KT_DES = KT
               ETAN_DES = DES_ETAN(phaseLL,phaseI) * sqrt(m_eff)
               ETAT_DES = DES_ETAT(phaseLL,phaseI) * sqrt(m_eff)
!}
            ENDIF
            ENDIF

! Calculate the normal contact force
            IF(USE_JKR .AND. JKR_FC_LOCAL > ZERO) THEN
               FN(:) = JKR_NORM(:) - ETAN_DES*V_REL_TRANS_NORM*NORMAL(:)
            ELSE
               JKR_CONTACT_STATE(CC) = .FALSE.
               PFR_NEIGHBOR(:,CC) = ZERO
            FN(:) =  -(KN_DES * OVERLAP_N * NORMAL(:) + &
               ETAN_DES * V_REL_TRANS_NORM * NORMAL(:))
            ENDIF

! Calculate the tangential overlap
            OVERLAP_T(:) = DTSOLID*VREL_T(:) + PFT_NEIGHBOR(:,CC)
            MAG_OVERLAP_T = sqrt(dot_product(OVERLAP_T,OVERLAP_T))

! Calculate the tangential contact force.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Calculate the tangential contact force.
               FT = -KT_DES*OVERLAP_T - ETAT_DES*VREL_T
               FTMD = sqrt(dot_product(FT,FT))
! Max force before the on set of frictional slip.
               IF(USE_JKR .AND. JKR_FC_LOCAL > ZERO) THEN
                  FNMD = MEW*ABS(-DOT_PRODUCT(JKR_NORM,NORMAL) + &
                     2.D0*JKR_FC_LOCAL)
               ELSE
               FNMD = MEW*sqrt(dot_product(FN,FN))
               ENDIF
! Frictional slip
               IF(FTMD > FNMD) THEN
! Direction of tangential force.
                  IF(USE_JKR .AND. JKR_FC_LOCAL > ZERO) THEN
                     FT = (FNMD/FTMD)*FT
                     OVERLAP_T = -(FT + ETAT_DES*VREL_T)/KT_DES
                  ELSE
                  TANGENT = OVERLAP_T/MAG_OVERLAP_T
                  FT = -FNMD * TANGENT
                  OVERLAP_T = (FNMD/KT_DES) * TANGENT
                  ENDIF
               ENDIF
            ELSE
               FT = 0.0
            ENDIF

            IF(USE_JKR .AND. JKR_FC_LOCAL > ZERO) &
               CALL JKR_DIAG_PAIR_EVENT('pair_force_before', &
                  IGLOBAL_ID(LL),IGLOBAL_ID(I),-1,CC,OVERLAP_N, &
                  -DOT_PRODUCT(FN,NORMAL),FT,PFT_NEIGHBOR(:,CC), &
                  PFR_NEIGHBOR(:,CC),JKR_CONTACT_STATE(CC))

! Save tangential displacement history
            PFT_NEIGHBOR(:,CC) = OVERLAP_T(:)
            IF (RIGID_GSP .OR. GSP_IMPLICIT) THEN
               TOW_TMP(:,1) = CROSS(dist2centerVECLL,FN) + CROSS(dist2centerVECLL,FT)
               TOW_TMP(:,2) = CROSS(dist2centerVECI ,-FN) + CROSS(dist2centerVECI ,-FT)
            ELSE
! calculate the distance from the particles' centers to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2
               DIST_CL = DIST_MAG/2.d0 + (DES_RADIUS(LL)**2 - &
                  DES_RADIUS(I)**2)/(2.d0*DIST_MAG)

               DIST_CI = DIST_MAG - DIST_CL

               TOW_force(:) = CROSS(NORMAL(:), FT(:))
               TOW_TMP(:,1) = DIST_CL*TOW_force(:)
               TOW_TMP(:,2) = DIST_CI*TOW_force(:)
            ENDIF
! Calculate the total force FC of a collision pair
! total contact force ( FC_TMP may already include cohesive force)
            FC_TMP(:) = FC_TMP(:) + FN(:) + FT(:)

! force and torque summation, glued sphere can use the same format of summation

!$omp atomic
            FC(LL,1) = FC(LL,1) + FC_TMP(1)
!$omp atomic
            FC(LL,2) = FC(LL,2) + FC_TMP(2)
!$omp atomic
            FC(LL,3) = FC(LL,3) + FC_TMP(3)
!$omp atomic
            FC(I,1) = FC(I,1) - FC_TMP(1)
!$omp atomic
            FC(I,2) = FC(I,2) - FC_TMP(2)
!$omp atomic
            FC(I,3) = FC(I,3) - FC_TMP(3)

            IF(RIGID_GSP .AND. (particle_state(LL) .EQ. 1)) THEN
!$omp atomic
               ppforce(gidLL,1) = ppforce(gidLL,1) + FC_TMP(1)
!$omp atomic
               ppforce(gidLL,2) = ppforce(gidLL,2) + FC_TMP(2)
!$omp atomic
               ppforce(gidLL,3) = ppforce(gidLL,3) + FC_TMP(3)
!$omp atomic
               ppforce(gidI,1) = ppforce(gidI,1) - FC_TMP(1)
!$omp atomic
               ppforce(gidI,2) = ppforce(gidI,2) - FC_TMP(2)
!$omp atomic
               ppforce(gidI,3) = ppforce(gidI,3) - FC_TMP(3)
            ENDIF

! for each particle the signs of norm and ft both flip, so add the same torque

!$omp atomic
            TOW(LL,1) = TOW(LL,1) + TOW_TMP(1,1)
!$omp atomic
            TOW(LL,2) = TOW(LL,2) + TOW_TMP(2,1)
!$omp atomic
            TOW(LL,3) = TOW(LL,3) + TOW_TMP(3,1)
!$omp atomic
            TOW(I,1)  = TOW(I,1)  + TOW_TMP(1,2)
!$omp atomic
            TOW(I,2)  = TOW(I,2)  + TOW_TMP(2,2)
!$omp atomic
            TOW(I,3)  = TOW(I,3)  + TOW_TMP(3,2)

            IF(RIGID_GSP .AND. (particle_state(LL) .EQ. 1)) THEN
!$omp atomic
               pptorque(gidLL,1) = pptorque(gidLL,1) + TOW_TMP(1,1)
!$omp atomic
               pptorque(gidLL,2) = pptorque(gidLL,2) + TOW_TMP(2,1)
!$omp atomic
               pptorque(gidLL,3) = pptorque(gidLL,3) + TOW_TMP(3,1)
!$omp atomic
               pptorque(gidI,1) = pptorque(gidI,1) + TOW_TMP(1,2)
!$omp atomic
               pptorque(gidI,2) = pptorque(gidI,2) + TOW_TMP(2,2)
!$omp atomic
               pptorque(gidI,3) = pptorque(gidI,3) + TOW_TMP(3,2)

            ENDIF

            DES_COLL_FORCE(I,:) = FC(I,:)

!######################################################################
! Force chain data - Implemented by Eric Breard (U. of Oregon)
!                    Reviewed by Jeff Dietiker
!######################################################################
! Note: Force Chain Data does not work in SMP
!######################################################################

            IF(WRITE_FORCE_CHAIN) THEN
! Resize arrays if needed
               IF(FCHAINC>=FCHAINC_OLD.AND.MOD(FCHAINC,100)==0) THEN
                  old_size = size(FCHAIN_MIDPOINT,1)
                  IF(FCHAINC>=INT(0.9*old_size)) THEN
                     new_size = INT(1.5*old_size)
                     call real_grow2_reverse(FCHAIN_MIDPOINT,new_size)
                     FCHAIN_MIDPOINT(old_size+1:new_size,:) = UNDEFINED
                     call real_grow2_reverse(FCHAIN_CONTACT_POINT,new_size)
                     FCHAIN_CONTACT_POINT(old_size+1:new_size,:) = UNDEFINED
                     call real_grow2_reverse(FCHAIN_ORIENT,new_size)
                     call real_grow2_reverse(FCHAIN_FN,new_size)
                     call real_grow2_reverse(FCHAIN_FT,new_size)
                     call real_grow2_reverse(FCHAIN_FC,new_size)
                     call real_grow(FCHAIN_LENGTH,new_size)
                     call real_grow(FCHAIN_FN_MAG,new_size)
                     call logical_grow(FCHAIN_FCMAX,new_size)
                     call integer_grow(FCHAIN_LOCAL_ID1,new_size)
                     call integer_grow(FCHAIN_LOCAL_ID2,new_size)
                     call integer_grow(FCHAIN_GLOBAL_ID1,new_size)
                     call integer_grow(FCHAIN_GLOBAL_ID2,new_size)
                     call real_grow(FCHAIN_OVERLAP,new_size)
                  ENDIF
               ENDIF
! Save data
               FCHAINC                     = FCHAINC + 1  ! Force chain counter
               FCHAIN_MIDPOINT(FCHAINC,:)  = 0.5D0 * (DES_POS_NEW(LL,:) + DES_POS_NEW(I,:)) ! Mid point coordinates
               FCHAIN_CONTACT_POINT(FCHAINC,:)  = DES_POS_NEW(LL,:) + NORMAL(:) * (DES_RADIUS(LL) + 0.5*OVERLAP_N) ! Contact point
               MAG_ORIENT                  = dsqrt(dot_product(FN,FN))
               if(MAG_ORIENT/=ZERO) FCHAIN_ORIENT(FCHAINC,:)    = FN(:)/MAG_ORIENT ! Orientation (same as Normal force)
               FCHAIN_FN(FCHAINC,:)        = FN(:)    ! Normal force
               FCHAIN_FT(FCHAINC,:)        = FT(:)    ! Tangential force
               FCHAIN_FC(FCHAINC,:)        = FC_TMP(:)! Total force (including cohesion)
               FCHAIN_LENGTH(FCHAINC)      = DIST_MAG !distance between LL and I
               FCHAIN_FN_MAG(FCHAINC)      = DSQRT(FN(1)**2 + FN(2)**2 + FN(3)**2) ! Normal force magnitude
               FCHAIN_LOCAL_ID1(FCHAINC)   = LL
               FCHAIN_LOCAL_ID2(FCHAINC)   = I
               FCHAIN_GLOBAL_ID1(FCHAINC)  = iglobal_id(LL)
               FCHAIN_GLOBAL_ID2(FCHAINC)  = iglobal_id(I)
               FCHAIN_OVERLAP(FCHAINC)     = OVERLAP_N
               IF(FCHAIN_FN_MAG(FCHAINC)>FC_MAG_MAX) THEN
                  FC_MAX_INDEX = FCHAINC  ! Keep track of max magnitude for pair LL-I
                  FC_MAG_MAX = FCHAIN_FN_MAG(FCHAINC)
               ENDIF

            ENDIF
!Rolling friction
! JFD: Here the rolling friction coefficient is non-dimensional. This is the equivalent to
! the Zhou paper (Physica A 269 (1999) 536-553) 's definition (Model A) divided by the
! particle diameter. Therefore here mew_r is multiplied by the diameter.

! TBD: implementation of rolling friction for glued particles @renjieke 1-1-2024
            IF(USE_JKR .AND. JKR_FC_LOCAL > ZERO) THEN
! Chen rolling displacement is stored in a canonical global-ID frame.
               V_ROLL(:) = -rad_eff*CROSS( &
                  OMEGA_NEW(LL,:)-OMEGA_NEW(I,:),NORMAL(:)) - &
                  0.5D0*(DES_RADIUS(I)-DES_RADIUS(LL))/ &
                  (DES_RADIUS(LL)+DES_RADIUS(I))*VREL_T(:)
               PFR_NEIGHBOR(:,CC) = PFR_NEIGHBOR(:,CC) + DTSOLID* &
                  MERGE(ONE,-ONE,IGLOBAL_ID(LL)<IGLOBAL_ID(I))*V_ROLL(:)
               XI_MAG = SQRT(DOT_PRODUCT(PFR_NEIGHBOR(:,CC), &
                  PFR_NEIGHBOR(:,CC)))
               XI_CRIT = JKR_THETAR_CRIT*rad_eff
               IF(XI_MAG > XI_CRIT .AND. XI_MAG > ZERO) &
                  PFR_NEIGHBOR(:,CC) = XI_CRIT*PFR_NEIGHBOR(:,CC)/XI_MAG
               KR_JKR = 4.D0*JKR_FC_LOCAL*A_OVER_A0**1.5D0
               MR_TORQUE(:) = -KR_JKR*PFR_NEIGHBOR(:,CC)
!$omp atomic
               TOW(LL,1) = TOW(LL,1) + MR_TORQUE(1)
!$omp atomic
               TOW(LL,2) = TOW(LL,2) + MR_TORQUE(2)
!$omp atomic
               TOW(LL,3) = TOW(LL,3) + MR_TORQUE(3)
!$omp atomic
               TOW(I,1) = TOW(I,1) - MR_TORQUE(1)
!$omp atomic
               TOW(I,2) = TOW(I,2) - MR_TORQUE(2)
!$omp atomic
               TOW(I,3) = TOW(I,3) - MR_TORQUE(3)
               CALL JKR_DIAG_PAIR_EVENT('pair_force_after', &
                  IGLOBAL_ID(LL),IGLOBAL_ID(I),-1,CC,OVERLAP_N, &
                  -DOT_PRODUCT(FN,NORMAL),FT,PFT_NEIGHBOR(:,CC), &
                  PFR_NEIGHBOR(:,CC),JKR_CONTACT_STATE(CC))
            ELSE
            IF (.not. GluedSphereDEM) THEN
            if(MEW_R>ZERO) then
               ! Particle LL
               OMEGA_MAG = sqrt(dot_product(OMEGA_NEW(LL,:),OMEGA_NEW(LL,:)))
               if(OMEGA_MAG > ZERO) then
                  RFTOW = -MEW_R*(2.0*DES_RADIUS(LL))*KN_DES * OVERLAP_N*OMEGA_NEW(LL,:)/OMEGA_MAG
                  TOW(LL,:) = TOW(LL,:) + RFTOW
               endif
               !particle i
               OMEGA_MAG = sqrt(dot_product(OMEGA_NEW(I,:),OMEGA_NEW(I,:)))
               if(OMEGA_MAG > ZERO) then
                  RFTOW = -MEW_R*(2.0*DES_RADIUS(I))*KN_DES * OVERLAP_N*OMEGA_NEW(I,:)/OMEGA_MAG
                  TOW(I,:) = TOW(I,:) + RFTOW
               endif
            endif
            ENDIF
            ENDIF
         ENDDO cc_loop ! Neighbor CC i.e. particle I loop
         DES_COLL_FORCE(LL,:) = FC(LL,:)

         IF(WRITE_FORCE_CHAIN.AND.FC_MAX_INDEX>0) FCHAIN_FCMAX(FC_MAX_INDEX) = .TRUE. ! Record where max magnitude occurs (for a given particle LL
      ENDDO ! Particle LL loop
!$omp end do

!$omp end parallel

! just for post-processing mag. of cohesive forces on each particle
! Divide by weight or mass if gravity is turned off
      IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
         IF(GRAV_MAG > ZERO) THEN
            DO LL = 1, MAX_PIP
               PostCohesive(LL) = dsqrt(Coh_Force(LL,1)**2+Coh_Force(LL,2)**2+Coh_Force(LL,3)**2)/(PMASS(LL)*GRAV_MAG)
            ENDDO
         ELSE
            DO LL = 1, MAX_PIP
               PostCohesive(LL) = dsqrt(Coh_Force(LL,1)**2+Coh_Force(LL,2)**2+Coh_Force(LL,3)**2)/PMASS(LL)
            ENDDO
         ENDIF
         deallocate(Coh_force)
      ENDIF

      IF(RIGID_GSP) THEN
         DO gid = 1, NGluedParticles
            ! ppforce(gid) should be zero if this gsp does not contact, so remove if(glueColNum(gid)==0) cycle should be safe
            if(all(ppforce(gid,:) == 0.0)) CYCLE
            glueForce(gid,:)  = glueForce(gid,:)  + ppforce(gid,:)
            if(all(pptorque(gid,:) == 0.0)) CYCLE
            glueTorque(gid,:) = glueTorque(gid,:) + pptorque(gid,:)
         ENDDO
         IF(allocated(ppforce)) deallocate(ppforce)
         IF(allocated(pptorque)) deallocate(pptorque)
      ENDIF

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

        include 'functions.inc'

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
!$omp critical (jkr_nr_warn_pp)
               WRITE(*,'(A,ES12.5,A,ES12.5,A,I0)') &
                  'WARNING: JKR pair Newton derivative small; y=',Y, &
                  ' x=',X_NEW,' iter=',ITER
!$omp end critical (jkr_nr_warn_pp)
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
!$omp critical (jkr_nr_warn_pp)
         WRITE(*,'(A,ES12.5,A,ES12.5)') &
            'WARNING: JKR pair Newton did not converge; y=',Y,' x=',X_NEW
!$omp end critical (jkr_nr_warn_pp)
         NEWTON_RAPHSON_METHOD = X_NEW
      END FUNCTION NEWTON_RAPHSON_METHOD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: print_excess_overlap                                    !
!                                                                      !
!  Purpose: Print overlap warning messages.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE PRINT_EXCESS_OVERLAP

      use error_manager

      IF(OVERLAP_N > flag_overlap*DES_RADIUS(LL) .OR.                  &
         OVERLAP_N > flag_overlap*DES_RADIUS(I)) THEN

         WRITE(ERR_MSG,1000) trim(iVAL(LL)), trim(iVAL(I)), S_TIME,    &
            DES_RADIUS(LL), DES_RADIUS(I), OVERLAP_N

         CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1000 FORMAT('WARNING: Excessive overlap detected between ',          &
         'particles ',A,' and ',/A,' at time ',g11.4,'.',/             &
         'RADII:  ',g11.4,' and ',g11.4,4x,'OVERLAP: ',g11.4)

      END SUBROUTINE PRINT_EXCESS_OVERLAP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: gsp_gsp_collison                                        !
!                                                                      !
!  Purpose: Calculate the collision between glued particles            !
!           when internal property distributions are not considered.   !
!  Author: Hang Zhou                                  Date: April 2024 !
!  Revised: Renjie Ke                                 Date: 11/01/2024 !
!  Revised: Renjie Ke                                 Date: 03/15/2025 !
!           Force the contact normal in the implicit method the same   !
!           as in the explicit method; fix the MPI issue.              !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE gsp_gsp_collison(index_ll, index_i, index_CC, VRN, VSLIP, NORM, overlapn, dist2I, dist2LL)


      use GSP_MATH_MOD, only: gsp_quat_to_exyz
      use error_manager, only: err_msg

      IMPLICIT NONE

      INTEGER, INTENT(IN)::index_i, index_ll, index_CC !particle id
! slip velocity at point of contact
      DOUBLE PRECISION, INTENT(OUT) :: VSLIP(:)
! normal component of relative contact velocity (scalar)
      DOUBLE PRECISION, INTENT(OUT) :: VRN
! unit normal vector along the line of contact pointing from
! particle L to particle II
      DOUBLE PRECISION, INTENT(OUT) :: NORM(3)
! overlap between particles
      DOUBLE PRECISIOn, INTENT(OUT) :: OVERLAPN
! temporary storage of torque
      DOUBLE PRECISIOn, INTENT(OUT) :: dist2I(3), dist2LL(3)

      INTEGER :: ii, ll, mi, mll, ss
      INTEGER :: num_contact, n_rows
      INTEGER, DIMENSION (:), ALLOCATABLE :: si, sll, si_tmp, sll_tmp ! list of spheres
      DOUBLE PRECISION :: rel_omgrad(3)
      double precision :: startt

      double precision :: sphere_radius_i, sphere_radius_ll
      double precision :: sphere_pos_i(3), sphere_pos_ll(3)
! bounding diameter of particle I and LL, they are not always d_p0(phase)
      double precision :: real_b_ii, real_b_ll, xloc, yloc, zloc, r_sum
      double precision :: exone_i(3),eyone_i(3),ezone_i(3)
      double precision :: exone_ll(3),eyone_ll(3),ezone_ll(3)
      double precision :: distll2i(3),distll2i_mag,dist_tmp(3)
      double precision :: lli_contact_point(3)
      double precision :: rvel(3)

      mi = PIJK(index_i, 5) ! phase id of particle i
      mll = PIJK(index_ll, 5) ! phase id of partice j

! n_rows should be the same for all phase
      n_rows = size(cspart_data,1)


      NORM(:) = (/0.D0,0.D0,0.D0/)
      lli_contact_point(:) = (/0.D0,0.D0,0.D0/)
      OVERLAPN = 0.d0
      num_contact = 0

      real_b_ii = 2 * des_radius(index_i)
      real_b_ll = 2 * des_radius(index_ll)

      exone_i=(/1.0D0, 0.0D0, 0.0D0/)
      eyone_i=(/0.0D0, 1.0D0, 0.0D0/)
      ezone_i=(/0.0D0, 0.0D0, 1.0D0/)
      CALL gsp_quat_to_exyz(glueQuat(index_i,:),exone_i,eyone_i,ezone_i)

      exone_ll=(/1.0D0, 0.0D0, 0.0D0/)
      eyone_ll=(/0.0D0, 1.0D0, 0.0D0/)
      ezone_ll=(/0.0D0, 0.0D0, 1.0D0/)
      CALL gsp_quat_to_exyz(glueQuat(index_ll,:),exone_ll,eyone_ll,ezone_ll)
! loop over spheres in gsp index_ll
   ! loop over spheres in gsp index_i
      do ll=1, n_rows
         if(cspart_data(ll,4,mll) == 0.0) exit
         if(cspart_data(ll,5,mll) == 0.0) cycle

         sphere_radius_ll = cspart_data(ll,4,mll) * real_b_ll / 2.0
         xloc = cspart_data(ll,1,mll) * real_b_ll
         yloc = cspart_data(ll,2,mll) * real_b_ll
         zloc = cspart_data(ll,3,mll) * real_b_ll

         sphere_pos_ll(:) = des_pos_new(INDEX_LL,:) + xloc * exone_ll &
                  + yloc * eyone_ll + zloc * ezone_ll

         do ii = 1, n_rows
            if(cspart_data(ii,4,mi) == 0.0) exit
            if(cspart_data(ii,5,mi) == 0.0) cycle

            sphere_radius_i = cspart_data(ii,4,mi) * real_b_ii / 2.0
            xloc = cspart_data(ii,1,mi) * real_b_ii
            yloc = cspart_data(ii,2,mi) * real_b_ii
            zloc = cspart_data(ii,3,mi) * real_b_ii

            sphere_pos_i(:) = des_pos_new(INDEX_I,:) + xloc * exone_i &
                     + yloc * eyone_i + zloc * ezone_i

            r_sum = sphere_radius_i + sphere_radius_ll
            ! distll2i vector points to i
            distll2i(:) =  sphere_pos_i(:) - sphere_pos_ll(:)
            distll2i_mag = dot_product(distll2i,distll2i)
            ! if two components do not contact
            IF(distll2i_mag > (r_sum - SMALL_NUMBER)**2) CYCLE

            IF (ENERGY_EQ) THEN
            ! Calculate conduction and radiation for thermodynamic neighbors
               IF(CALC_COND_DES(mll)) THEN
                  QQ_TMP = DES_CONDUCTION(index_ll, index_i, sqrt(distll2i_mag), PIJK(index_ll,5), PIJK(index_ll,4), &
                     sphere_radius_i, sphere_radius_ll)
!$omp atomic
                  Q_Source(index_ll) = Q_Source(index_ll) + QQ_TMP
!$omp atomic
                  Q_Source(index_i) = Q_Source(index_i) - QQ_TMP
               ENDIF
            ENDIF

            distll2i_mag = sqrt(distll2i_mag)
            NORM(:) = NORM(:) + distll2i(:)/distll2i_mag
            OVERLAPN = OVERLAPN + r_sum - distll2i_mag
            ! the location of contact point should also be examed
            lli_contact_point(:)  = lli_contact_point(:) + sphere_pos_ll &
               + (sphere_radius_ll - (r_sum - distll2i_mag)/2)*distll2i(:)/distll2i_mag
            num_contact = num_contact + 1
         enddo ! end of ii
      enddo ! end of ll

      if(num_contact==0) then
         PFT_NEIGHBOR(:,index_CC) = 0.0
         return

      endif
      OVERLAPN = OVERLAPN/num_contact
      ! norm direction needs more carefully consideration!
      norm(:) = norm(:)/num_contact
      lli_contact_point(:) = lli_contact_point(:)/num_contact
      r_sum = DES_RADIUS(index_i) + DES_RADIUS(index_ll)

      dist2LL = lli_contact_point- des_pos_new(index_ll,:)
      dist2I = lli_contact_point- des_pos_new(index_i,:)
      rel_omgrad = CROSS(OMEGA_NEW(index_ll,:),dist2ll) - CROSS(OMEGA_NEW(index_i,:),dist2i)
      rvel(:) = des_vel_new(index_ll,:) - des_vel_new(index_i,:) + rel_omgrad

      VRN = DOT_PRODUCT(rvel,norm)
      VSLIP(:) = rvel(:) - VRN*norm

      END SUBROUTINE gsp_gsp_collison

    END SUBROUTINE CALC_FORCE_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_COORDINATION                                       C
!  Purpose: Compute coordination number and overlap statistics         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 17-SEP-24  C
!  Reviewer:                                                           C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CALC_COORDINATION

! Modules
!---------------------------------------------------------------------//
      USE discretelement
      USE param
      USE param1
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I, LID1, LID2
      double precision :: OVERLAP
!---------------------------------------------------------------------//

! Coordination number can be extracted from force chain data
! A force chain links two particles (Local IDs LID1 and LID2), so the
! coordination in incremented for each force chain. We can keep track of
! the min, max and mean overlap for each particle.

      IF(.NOT.WRITE_FORCE_CHAIN) RETURN

      COORDINATION(:) = ZERO
      MIN_OVERLAP(:) = UNDEFINED
      MAX_OVERLAP(:) = -UNDEFINED
      MEAN_OVERLAP(:) = ZERO


      DO I = 1, FCHAINC

         LID1 = FCHAIN_LOCAL_ID1(I)
         LID2 = FCHAIN_LOCAL_ID2(I)

         OVERLAP = FCHAIN_OVERLAP(I)

         COORDINATION(LID1) = COORDINATION(LID1) + 1
         COORDINATION(LID2) = COORDINATION(LID2) + 1

         MIN_OVERLAP(LID1) = DMIN1(MIN_OVERLAP(LID1),OVERLAP)
         MAX_OVERLAP(LID1) = DMAX1(MAX_OVERLAP(LID1),OVERLAP)
         MEAN_OVERLAP(LID1) = MEAN_OVERLAP(LID1) + OVERLAP

         MIN_OVERLAP(LID2) = DMIN1(MIN_OVERLAP(LID2),OVERLAP)
         MAX_OVERLAP(LID2) = DMAX1(MAX_OVERLAP(LID2),OVERLAP)
         MEAN_OVERLAP(LID2) = MEAN_OVERLAP(LID2) + OVERLAP

      ENDDO

      WHERE(COORDINATION>ZERO)
         MEAN_OVERLAP(:) = MEAN_OVERLAP(:) / COORDINATION(:)
      ELSEWHERE
         MIN_OVERLAP(:) = ZERO
         MAX_OVERLAP(:) = ZERO
      ENDWHERE

! SuperDEM force chain is made up of two segments, so we need to divide
! coordination by two (min, max, mean are the same)
      IF(SUPERDEM) COORDINATION(:) = COORDINATION(:) / 2.0D0

      RETURN

      END SUBROUTINE CALC_COORDINATION

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_VIRTUAL_BOND_FORCE                                 !
!  Author: Renjie Ke                                  Date: 01-Aug-25  !
!                                                                      !
!  Purpose: Calculate virtual bond force and torque between two        !
!           component spheres in the same GSP particles.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
   SUBROUTINE CALC_VIRTUAL_BOND_FORCE(pLL, pII, dF1, dF2, dT1, dT2, nbid)

      USE discretelement
      use constant, only: pi
      use gsp_math_mod, only: rotmat_from_normals
      use functions, only: is_nonexistent
      use compar, only: pe_io, mype, numPEs
      use run, only: energy_eq
      use param1, only: one, small_number, zero, undefined
      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: pLL, pII, nbid
      DOUBLE PRECISION, INTENT(INOUT) :: dF1(3), dF2(3), dT1(3), dT2(3)

      ! Local variables
      ! normal force and shear force
      DOUBLE PRECISION :: dFbond_elastic_norm(3), Fdamp_norm(3)
      DOUBLE PRECISION :: dFbond_elastic_tan(3),Fdamp_tan(3)
      DOUBLE PRECISION :: dM_shear_bending(3)

      ! torque
      DOUBLE PRECISION :: dM_elastic_norm(3), Mdamp_norm(3)
      DOUBLE PRECISION :: dM_elastic_tan(3), Mdamp_tan(3)

      ! shear directed torque
      DOUBLE PRECISION :: M_shearLL(3), M_shearII(3)
      DOUBLE PRECISION :: Ftcurrent(3)

      ! translational relative velocities
      DOUBLE PRECISION :: v_rel(3), v_rel_n(3),v_rel_t(3)
      ! angular relative velocities
      DOUBLE PRECISION :: omega_rel(3), omega_rel_n(3), omega_rel_t(3)
      ! sphere center to contact point
      ! if distance between 2 spheres are larger than (r1+r2), this contact point is actually the the bond center
      double precision :: r_pLL(3), r_pII(3)
      double precision :: vcp_LL(3), vcp_II(3), v_rel_total(3), v_pLL(3), v_pII(3)

      ! sphere positions and relative positions
      DOUBLE PRECISION :: x_rel(3), dist_mag
      DOUBLE PRECISION :: normal(3)
      DOUBLE PRECISION :: length_vector(3)

      ! Bond properties
      DOUBLE PRECISION :: initial_bond_length, area_bond
      DOUBLE PRECISION :: rad_eff, mass_eff, moi_eff
      DOUBLE PRECISION :: moi_bond, pmoi_bond
      DOUBLE PRECISION :: stiff_bond_norm, stiff_bond_tan, stiff_bond_bending, stiff_bond_twisting
      DOUBLE PRECISION :: eta_bond_norm, eta_bond_tan, eta_bond_bending, eta_bond_twisting
      ! Phase id
      INTEGER :: PID, PID_II
      ! Shear modulus
      DOUBLE PRECISION :: G_bond ! This should be calculate from E_YOUNG_BOND or it should be set by User
      ! Relative normal velocity and force/torque magnitude
      DOUBLE PRECISION :: dvn_rel, dFN, dwn_rel, dMN
      ! Rotational velocity
      DOUBLE PRECISION :: omegaLL(3),omegaII(3)
      ! Normal force/torque vector
      DOUBLE PRECISION :: FNvec(3), MNvec(3), FTvec(3), MTvec(3), Rmat(3,3)
      DOUBLE PRECISION :: FT_hist(3), MT_hist(3)
      ! Current normal and pre time step norma3l
      DOUBLE PRECISION :: n_old(3)
      ! Sum of two sphere radius
      DOUBLE PRECISION :: radius_sum
      ! The order of neighbors in gp_neighsid, this value should always be nbid-1 or nbid+1
      INTEGER :: nbid_ii
      ! Flexible GSP breakage model parameters
      DOUBLE PRECISION :: sigma_n, sigma_tau
      DOUBLE PRECISION :: tensile_strength, shear_strength
      logical :: isbreak

      ! Debug
      LOGICAL :: rk_debug = .False.
      ! Debug used values
      DOUBLE PRECISION :: POS1(3)
      DOUBLE PRECISION :: POS2(3)
      DOUBLE PRECISION :: POS3(3)
      DOUBLE PRECISION :: POS4(3)
      DOUBLE PRECISION :: centeroid(3)
      DOUBLE PRECISION :: damping_ratio
      DOUBLE PRECISION :: Eb
      double precision :: delta_vc, ell1, ell2
      double precision :: temp_hist_force

      ! Fdamp_norm(:) = 0.0
      ! Fdamp_tan(:) = 0.0

!----------------------------------data check----------------------------------
      ! a check for pijk(pII,5) == pijk(pLL,5) should be performed in case!
      PID = pijk(pLL, 5)
      PID_II = pijk(pII, 5)
      IF (PID /= PID_II) THEN
         WRITE(*,*) "Error: These two component spheres do not belong to same parent GSP"
         WRITE(*,*) "glueID of pLL: ", PID
         WRITE(*,*) "glueID of pII: ", PID_II
         ERROR STOP
      ENDIF
      IF(IS_NONEXISTENT(pII)) ERROR STOP
!----------------------------------data check----------------------------------

!----------------------------------relative velocity----------------------------------
      dF1(:) = 0.0D0
      dF2(:) = 0.0D0
      dT1(:) = 0.0D0
      dT2(:) = 0.0D0

      if(nbid == 1 .or. nbid == 3 .or. nbid == 5) then
         nbid_ii = nbid + 1
      else
         nbid_ii = nbid - 1
      endif

      ! relative position of two spheres
      x_rel(:) = des_pos_new(pII,:) - des_pos_new(pLL,:)
      ! Relative translational velocity of two spheres
      v_rel = des_vel_new(pII,:) - des_vel_new(pLL,:)
      ! Relative angular velocity of two spheres
      omega_rel = omega_new(pII,:) - omega_new(pLL,:)

      dist_mag = dot_product(x_rel,x_rel)
      dist_mag = SQRT(dist_mag)

      IF (dist_mag > 1.0e-16) THEN
         normal(:) = x_rel(:)/dist_mag
      ELSE
         normal(:) = 0.0D0
      ENDIF

      ! n_old(1) = n_histx(pLL, pII)
      ! n_old(2) = n_histy(pLL, pII)
      ! n_old(3) = n_histz(pLL, pII)

      n_old(1) = bnx(pLL, nbid)
      n_old(2) = bny(pLL, nbid)
      n_old(3) = bnz(pLL, nbid)

      ! first step history is empty
      IF (ALL(n_old == 0.0D0)) n_old(:) = normal(:)
      Rmat(:,:) = 0.0D0
      ! obtain the rotational matrix to rotate the contact plane
      CALL rotmat_from_normals(n_old, normal, Rmat)

      ! save the normal history for next step
      ! n_histx(pLL, pII) = normal(1)
      ! n_histy(pLL, pII) = normal(2)
      ! n_histz(pLL, pII) = normal(3)

      bnx(pLL, nbid) = normal(1)
      bny(pLL, nbid) = normal(2)
      bnz(pLL, nbid) = normal(3)

      bnx(pII, nbid_ii) = -normal(1)
      bny(pII, nbid_ii) = -normal(2)
      bnz(pII, nbid_ii) = -normal(3)


      ! method 1 --------------------------------------------------------
      radius_sum = des_radius(pLL) + des_radius(pII)
      ! this is to evaluate the location of the linear velocity
      ! but in whatever case, we still think that the contact plane is located at the bond center
      if (radius_sum .le. dist_mag) then
      !! if no overlap, the linear velocity is evaluate on the sphere surface
          r_pLL(:) =  des_radius(pLL) * normal
          r_pII(:) = -des_radius(pII) * normal
      else
      !! if has overlap, the linear velocity is evaluated at the bond center
         r_pLL(:) = 0.5D0 * dist_mag * normal
         r_pII(:) = -0.5D0 * dist_mag * normal
      endif

      ! delta_vc = radius_sum - dist_mag
      ! ell1 = des_radius(pLL) - 0.5D0 * delta_vc
      ! ell1 = des_radius(pII) - 0.5D0 * delta_vc
      ! r_pLL = ell1 * normal
      ! r_pII = -ell2 * normal
      ! method 2 --------------------------------------------------------
      ! this assume the contact plane is always in the center of the bond
      ! even if the bond is stretched, this might introduce extra linear velocities
      ! so make the shear force and shear induced moment bigger
      ! however, I need to test this @renjieke
      ! r_pLL(:) = 0.5D0 * dist_mag * normal
      ! r_pII(:) = -0.5D0 * dist_mag * normal
      omegaLL = omega_new(pLL,:)
      omegaII = omega_new(pII,:)
      v_pLL(:) = CROSS(omegaLL, r_pLL)
      v_pII(:) = CROSS(omegaII, r_pII)

      ! overall translational velocities at the bond center, some might also call this the contact point
      vcp_LL = des_vel_new(pLL,:) + v_pLL(:)
      vcp_II = des_vel_new(pII,:) + v_pII(:)

      ! calculate reltaive velocity in tangential direction

      ! Normal relative velocity
      v_rel_n(:) = dot_product(vcp_II(:)-vcp_LL(:), normal(:)) * normal(:) ! dot_product(v_rel, normal)
      dvn_rel = dot_product(vcp_II(:)-vcp_LL(:), normal(:))

      ! Tangential relative velocity
      v_rel_t(:) = vcp_II(:) - vcp_LL(:) - v_rel_n(:)

      ! Normal relative angular velocity, they affect the twisting
      omega_rel_n(:) = dot_product(omega_rel,normal) * normal(:)
      dwn_rel = dot_product(omega_rel,normal)

      ! Tangential relatve angular velocity, they affect the bending
      omega_rel_t(:) = omega_rel(:) - dot_product(omega_rel, normal) * normal(:)
!----------------------------------relative velocity----------------------------------


!----------------------------------bond properties----------------------------------

      ! @renjieke, I might want to change the sc2gpc_vec in later, so this might not be good enough
      ! the question is, should the stiffness using the initial_bond_length only or the current length
      length_vector(:) = sc2gpc_vec(pII,:) - sc2gpc_vec(pLL,:)
      initial_bond_length = dot_product(length_vector,length_vector)
      initial_bond_length = sqrt(initial_bond_length)

      ! The component should not overlap with each other, this will cause zero bond length and this is unacceptable
      IF(initial_bond_length == 0.0) THEN
         WRITE(*,*) "Error: Initially overlapping components"
         ERROR STOP
      ENDIF

      ! calculate the effective radius
      rad_eff = 2.0D0 * des_radius(pII)*des_radius(pLL)/(des_radius(pII) + des_radius(pLL))
      ! cross section of the bond, the bond is assumed to be a cylinder
      area_bond = pi*rad_eff**(2.0)
      ! calculate the effective mass
      mass_eff = 2.0D0 * pmass(pII)*pmass(pLL)/(pmass(pII) + pmass(pLL))

      ! calculate the bond normal stiffness, this is based on provided young's modulus
      stiff_bond_norm = E_YOUNG_BOND(PID)*area_bond/initial_bond_length
      ! shear modulus is derived from provided young's modulus
      G_bond = E_YOUNG_BOND(PID)/(2.0d0 * (1.0d0 + V_POISSON_BOND(PID)))
      ! calculate the bond shear stiffness based on the shear modulus
      stiff_bond_tan = G_bond*area_bond/initial_bond_length

      ! area moment of inertia, bending module
      moi_bond = pi*rad_eff**(4.0d0)/4.0d0
      stiff_bond_bending = E_YOUNG_BOND(PID)*moi_bond/initial_bond_length
      ! polar moment of inertia, twisting module
      pmoi_bond = pi*rad_eff**(4.0d0)/2.0d0
      stiff_bond_twisting = G_bond*pmoi_bond/initial_bond_length

      ! damping ratio == 1: critical damped
      ! damping ratio < 1: under damped
      ! damping ratio > 1: over damped
      damping_ratio = C_ratio_bond(PID)

      ! damping coefficient for normal and tangential displacement
      eta_bond_norm = 2.D0 * damping_ratio * SQRT(stiff_bond_norm * mass_eff)
      eta_bond_tan = 2.D0 * damping_ratio * SQRT(stiff_bond_tan * mass_eff)

      ! sphere moment of inertia
      ! the effective mass and radius are used here to calculate the sphere moment of inertia
      ! however, if two spheres have huge difference in size, this might cause issue
      ! I need to go back and check the reference paper
      moi_eff = 2.0d0/5.0d0 * mass_eff * rad_eff ** (2.0D0)
      ! damping coefficient for twisting and bending
      eta_bond_twisting = 2.D0 * damping_ratio * SQRT(stiff_bond_twisting * moi_eff)
      eta_bond_bending = 2.D0 * damping_ratio * SQRT(stiff_bond_bending * moi_eff)

!----------------------------------bond properties----------------------------------



!--------------------incremental bond force calculation starts here-----------------

      dFbond_elastic_norm(:) = 0.0D0
      dFbond_elastic_tan(:) = 0.0D0
      dM_elastic_norm(:) = 0.0D0
      dM_elastic_tan(:) = 0.0D0

      temp_hist_force = bfn(pLL, nbid)

      dFN = stiff_bond_norm * dvn_rel * dtsolid ! magnitude (local)
      bfn(pLL, nbid) = bfn(pLL, nbid) + dFN ! accumulated force magnitude (global array)
      bfn(pII, nbid_ii) = bfn(pLL, nbid)
      FNvec(:) = bfn(pLL, nbid) * normal ! accumulated force vector (local)

      ! Normal damping
      Fdamp_norm(:) = -eta_bond_norm * v_rel_n(:)

      ! Rotate the shear force before accumulating
      FT_hist(1) = bftx(pLL, nbid)
      FT_hist(2) = bfty(pLL, nbid)
      FT_hist(3) = bftz(pLL, nbid)
      FTvec(:) = MATMUL(Rmat, FT_hist)

      ! Tangential force
      dFbond_elastic_tan(:) = stiff_bond_tan * v_rel_t(:) * dtsolid

      bftx(pLL, nbid) = FTvec(1) + dFbond_elastic_tan(1)
      bfty(pLL, nbid) = FTvec(2) + dFbond_elastic_tan(2)
      bftz(pLL, nbid) = FTvec(3) + dFbond_elastic_tan(3)

      bftx(pII, nbid_ii) = - bftx(pLL, nbid)
      bfty(pII, nbid_ii) = - bfty(pLL, nbid)
      bftz(pII, nbid_ii) = - bftz(pLL, nbid)

      FTvec(1) = bftx(pLL, nbid)
      FTvec(2) = bfty(pLL, nbid)
      FTvec(3) = bftz(pLL, nbid)

      ! Tangential damping
      Fdamp_tan(:) = -eta_bond_tan * v_rel_t(:)

      dMN = stiff_bond_twisting * dwn_rel * dtsolid ! magnitude (local)
      bmn(pLL, nbid) = bmn(pLL, nbid) + dMN ! accumulated normal torque magnitude (global array)
      bmn(pII, nbid_ii) = bmn(pLL, nbid)
      MNvec(:) = bmn(pLL, nbid) * normal ! accumulated normal torque vector (local)

      ! damping for omega in normal direction
      Mdamp_norm(:) = -eta_bond_twisting * omega_rel_n(:)

      ! Rotate the shear moment before accumulating
      MT_hist(1) = bmtx(pLL, nbid)
      MT_hist(2) = bmty(pLL, nbid)
      MT_hist(3) = bmtz(pLL, nbid)
      MTvec(:) = MATMUL(Rmat, MT_hist)

      ! incremental moment in the bond tangential they affect the bending
      dM_elastic_tan = stiff_bond_bending * omega_rel_t(:) * dtsolid

      bmtx(pLL, nbid) = MTvec(1) + dM_elastic_tan(1)
      bmty(pLL, nbid) = MTvec(2) + dM_elastic_tan(2)
      bmtz(pLL, nbid) = MTvec(3) + dM_elastic_tan(3)

      bmtx(pII, nbid_ii) = - bmtx(pLL, nbid)
      bmty(pII, nbid_ii) = - bmty(pLL, nbid)
      bmtz(pII, nbid_ii) = - bmtz(pLL, nbid)

      MTvec(1) = bmtx(pLL, nbid)
      MTvec(2) = bmty(pLL, nbid)
      MTvec(3) = bmtz(pLL, nbid)

      ! damping for omege in tangential direction
      Mdamp_tan(:) = -eta_bond_bending * omega_rel_t(:)

      ! Total force for sphere LL and II
      dF1(:) = FNvec(:) + FTvec(:) - Fdamp_norm(:) - Fdamp_tan(:)
      dF2(:) = - dF1(:)

      ! Total Torque for sphere LL and II
      dT1(:) = MNvec(:) + MTvec(:) - Mdamp_norm(:) - Mdamp_tan(:)
      dT2(:) = -dT1(:)

      ! shear directed moment (r x F_tan)
      M_shearLL(:) =  0.0D0
      M_shearII(:) =  0.0D0

      ! use the total or just the shear component, theoretically they are the same
      M_shearLL(:) = CROSS(r_pLL, dF1)
      M_shearII(:) = CROSS(r_pII, dF2)
      ! M_shearLL(:) = CROSS(r_pLL, FTvec - Fdamp_tan)
      ! M_shearII(:) = CROSS(r_pII, -FTvec + Fdamp_tan)

      dT1(:) = dT1(:) + M_shearLL(:)
      dT2(:) = dT2(:) + M_shearII(:)

      sigma_n = sqrt(dot_product(FNvec, FNvec)) / area_bond + sqrt(dot_product(MTvec, MTvec)) * rad_eff / moi_bond
      sigma_tau = sqrt(dot_product(FTvec, FTvec)) / area_bond + sqrt(dot_product(MNvec, MNvec)) * rad_eff / pmoi_bond

      tensile_strength = tensile_strength_bond(PID)
      shear_strength = shear_strength_bond(PID)

      ! if virtual bond breakage is turned on, check if tensile/shear strength are provided.
      if(flexible_gsp_bond_break) then
         if(tensile_strength == undefined .and. shear_strength == undefined) then
            write(err_msg, 1495)
            call log_error()
         endif
         if(tensile_strength == undefined) then
            write(err_msg, 1496)
            call log_info()
         endif
         if(shear_strength == undefined) then
            write(err_msg, 1497)
            call log_info()
         endif
      endif

 1495  FORMAT('Both virtual bond strengths are undefined.', &
            /'Please verify the flexible gsp settings.')
 1496  FORMAT('Warning: only shear strength of the virtual bond is defined.')
 1497  FORMAT('Warning: only tensile strength of the virtual bond is defined.')

      isbreak = ((sigma_n .ge. tensile_strength) .or. (sigma_tau .ge. shear_strength)) .and. flexible_gsp_bond_break

      if(isbreak) then
      ! reset the force, torque and normal history
         bfn(pLL, nbid) = 0.0D0
         bfn(pII, nbid_ii) = 0.0D0

         bmn(pLL, nbid) = 0.0D0
         bmn(pII, nbid_ii) = 0.0D0

         bftx(pLL, nbid) = 0.0D0
         bfty(pLL, nbid) = 0.0D0
         bftz(pLL, nbid) = 0.0D0

         bftx(pII, nbid_ii) = 0.0D0
         bfty(pII, nbid_ii) = 0.0D0
         bftz(pII, nbid_ii) = 0.0D0

         bmtx(pLL, nbid) = 0.0D0
         bmty(pLL, nbid) = 0.0D0
         bmtz(pLL, nbid) = 0.0D0

         bmtx(pII, nbid_ii) = 0.0D0
         bmty(pII, nbid_ii) = 0.0D0
         bmtz(pII, nbid_ii) = 0.0D0

         bnx(pLL, nbid) = 0.0D0
         bny(pLL, nbid) = 0.0D0
         bnz(pLL, nbid) = 0.0D0

         bnx(pII, nbid_ii) = 0.0D0
         bny(pII, nbid_ii) = 0.0D0
         bnz(pII, nbid_ii) = 0.0D0

         dF1(:) = 0.0D0
         dF2(:) = 0.0D0

         dT1(:) = 0.0D0
         dT2(:) = 0.0D0
      ! remove the connection in gp_neighsid
         gp_neighsid(pLL, nbid) = -1
         gp_neighsid(pII, nbid_ii) = -1
         if(energy_eq) then
         ! the contact area now exposes to the surrounding
            gp_sa(pLL) = gp_sa(pLL) + gp_neighsa(pLL, nbid)
            gp_sa(pII) = gp_sa(pII) + gp_neighsa(pII, nbid_ii)
         ! reset the gp_neighsa
            gp_neighsa(pLL,nbid) = 0.0D0
            gp_neighsa(pII, nbid_ii) = 0.0D0
         endif
      endif

      RETURN
   END SUBROUTINE CALC_VIRTUAL_BOND_FORCE

 END MODULE CALC_FORCE_DEM_MOD
