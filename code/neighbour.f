MODULE NEIGHBOUR_MOD

   USE NSQUARE_MOD, ONLY: NSQUARE

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NEIGHBOUR                                              C
!  Purpose: DES - Neighbors search;
!           N-Square,
!           Quadtree(2D)/Octree(3D)  (use at own risk)
!           Cell linked
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

   SUBROUTINE NEIGHBOUR

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE discretelement
      use desgrid
      USE MPI_UNPACK_DES, ONLY: RESTORE_PENDING_PAIR_HISTORY, &
         JKR_DIAG_PAIR_EVENT, JKR_ROTATION_HISTORY_EVENT, &
         JKR_DIAGNOSTICS_ACTIVE
      Use des_thermo
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

INTEGER :: cc,ll,cc_start,cc_end,cc_start_old,cc_end_old,cc_old
INTEGER :: reverse_owner, reverse_start_old, reverse_end_old
LOGICAL :: found

!-----------------------------------------------
! Reset PPOS and NEIGHBOURS back to initialized values
      PPOS(:,:) = DES_POS_NEW(:,:)
      neighbor_index_old(:) = neighbor_index(:)

!$omp parallel do default(none) private(cc)                                   &
!$omp shared(neighbors, neighbors_old, pft_neighbor, pft_neighbor_old,        &
!$omp      pfr_neighbor,pfr_neighbor_old,pftw_neighbor,pftw_neighbor_old,     &
!$omp      jkr_frame_neighbor,jkr_frame_neighbor_old,                         &
!$omp      jkr_contact_state,                                                  &
!$omp      jkr_contact_state_old,                                             &
!$omp contact_point_a,contact_point_a_old,contact_point_b,contact_point_b_old,&
!$omp      contact_lambda_a,contact_lambda_a_old,                             &
!$omp      contact_lambda_b,contact_lambda_b_old, SuperDEM)
      do cc=1, size(neighbors)
         neighbors_old(cc) = neighbors(cc)
         pft_neighbor_old(:,cc) = pft_neighbor(:,cc)
         pfr_neighbor_old(:,cc) = pfr_neighbor(:,cc)
         jkr_frame_neighbor_old(:,cc) = jkr_frame_neighbor(:,cc)
         pftw_neighbor_old(cc) = pftw_neighbor(cc)
         jkr_contact_state_old(cc) = jkr_contact_state(cc)
         if (SuperDEM) then
            contact_point_a_old(:,cc)=contact_point_a(:,cc)
            contact_point_b_old(:,cc)=contact_point_b(:,cc)
            contact_lambda_a_old(cc)=contact_lambda_a(cc)
            contact_lambda_b_old(cc)=contact_lambda_b(cc)
         endif
      enddo
!$omp end parallel do
      IF(USE_JKR) THEN
         IF(JKR_DIAGNOSTICS_ACTIVE()) THEN
            DO LL=1,MAX_PIP
               CC_START = 1
               IF(LL > 1) CC_START = NEIGHBOR_INDEX_OLD(LL-1)
               CC_END = NEIGHBOR_INDEX_OLD(LL)
               DO CC=CC_START,CC_END-1
                  IF(NEIGHBORS_OLD(CC) <= 0) CYCLE
                  CALL JKR_DIAG_PAIR_EVENT('neighbor_before_rebuild', &
                     IGLOBAL_ID(LL),IGLOBAL_ID(NEIGHBORS_OLD(CC)), &
                     -1,CC,0D0,0D0,(/0D0,0D0,0D0/), &
                     PFT_NEIGHBOR_OLD(:,CC),PFR_NEIGHBOR_OLD(:,CC), &
                     JKR_CONTACT_STATE_OLD(CC))
               ENDDO
            ENDDO
         ENDIF
      ENDIF


      NEIGHBOR_INDEX(:) = 0

      IF (DES_NEIGHBOR_SEARCH.EQ.1) THEN
         CALL NSQUARE
      ELSEIF (DES_NEIGHBOR_SEARCH.EQ.4) THEN
          CALL DESGRID_NEIGH_BUILD
      ENDIF

!$omp parallel do default(none)                                         &
!$omp private(cc,cc_old,ll,found,cc_start,cc_end,cc_start_old,cc_end_old,&
!$omp    reverse_owner,reverse_start_old,reverse_end_old)                &
!$omp shared(max_pip,neighbors,neighbor_index,neighbor_index_old,       &
!$omp    neighbors_old,pft_neighbor,pft_neighbor_old,                   &
!$omp    pfr_neighbor,pfr_neighbor_old,pftw_neighbor,pftw_neighbor_old,&
!$omp    jkr_frame_neighbor,jkr_frame_neighbor_old,                     &
!$omp    jkr_contact_state,                                              &
!$omp    jkr_contact_state_old,use_jkr,                                 &
!$omp    contact_point_a,contact_point_a_old,                           &
!$omp    contact_point_b,contact_point_b_old,                           &
!$omp    contact_lambda_a,contact_lambda_a_old,                         &
!$omp    contact_lambda_b,contact_lambda_b_old,                         &
!$omp    SuperDEM)
      do ll = 1, max_pip

         CC_START = 1
         IF (LL.gt.1) CC_START = NEIGHBOR_INDEX(LL-1)
         CC_END   = NEIGHBOR_INDEX(LL)

         CC_START_OLD = 1
         IF (LL.gt.1) CC_START_OLD = NEIGHBOR_INDEX_OLD(LL-1)
         CC_END_OLD   = NEIGHBOR_INDEX_OLD(LL)

         DO CC = CC_START, CC_END-1
            found = .false.
            DO CC_OLD = CC_START_OLD, CC_END_OLD-1
               if (neighbors(cc) .eq. neighbors_old(cc_old)) then
                  pft_neighbor(:,cc) = pft_neighbor_old(:,cc_old)
                  pfr_neighbor(:,cc) = pfr_neighbor_old(:,cc_old)
                  jkr_frame_neighbor(:,cc) = jkr_frame_neighbor_old(:,cc_old)
                  pftw_neighbor(cc) = pftw_neighbor_old(cc_old)
                  jkr_contact_state(cc) = jkr_contact_state_old(cc_old)
!SuperDEM
                     if (SuperDEM) then
                         contact_point_a(:,cc)=contact_point_a_old(:,cc_old)
                         contact_point_b(:,cc)=contact_point_b_old(:,cc_old)
                         contact_lambda_a(cc)=contact_lambda_a_old(cc_old)
                         contact_lambda_b(cc)=contact_lambda_b_old(cc_old)
                     endif
                  found = .true.
                  exit
               endif
            enddo

! JKR pair ownership may reverse after sorting or owner migration.
! PFT is owner-directed; PFR and the neck state are pair-invariant.
            IF(.NOT.found .AND. USE_JKR) THEN
               reverse_owner = neighbors(cc)
               IF(reverse_owner >= 1 .AND. reverse_owner <= max_pip) THEN
                  reverse_start_old = 1
                  IF(reverse_owner > 1) &
                     reverse_start_old = neighbor_index_old(reverse_owner-1)
                  reverse_end_old = neighbor_index_old(reverse_owner)
                  DO cc_old = reverse_start_old, reverse_end_old-1
                     IF(neighbors_old(cc_old) == ll) THEN
                        pft_neighbor(:,cc) = -pft_neighbor_old(:,cc_old)
                        pfr_neighbor(:,cc) = pfr_neighbor_old(:,cc_old)
                        jkr_frame_neighbor(1:3,cc) = &
                           -jkr_frame_neighbor_old(1:3,cc_old)
                        jkr_frame_neighbor(4,cc) = &
                           jkr_frame_neighbor_old(4,cc_old)
                        pftw_neighbor(cc) = pftw_neighbor_old(cc_old)
                        jkr_contact_state(cc) = &
                           jkr_contact_state_old(cc_old)
                        found = .TRUE.
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF

            if (.not.found) then
               pft_neighbor(:,cc) = 0.0
               pfr_neighbor(:,cc) = 0.0
               jkr_frame_neighbor(:,cc) = 0.0
               pftw_neighbor(cc) = 0.0
               jkr_contact_state(cc) = .FALSE.
            endif
!SuperDEM
            if (.not.found) then
               if (SuperDEM) then
                    contact_point_a(:,cc)=0.0   !not contact
                    contact_point_b(:,cc)=0.0
                    contact_lambda_a(cc)=0.0
                    contact_lambda_b(cc)=0.0
               endif
            endif


         enddo
      enddo
!$omp end parallel do

! Apply DMP owner-migration histories only after valid CC slots exist.
      IF(USE_JKR) THEN
         CALL RESTORE_PENDING_PAIR_HISTORY
         IF(JKR_DIAGNOSTICS_ACTIVE()) THEN
            DO LL=1,MAX_PIP
               CC_START = 1
               IF(LL > 1) CC_START = NEIGHBOR_INDEX(LL-1)
               CC_END = NEIGHBOR_INDEX(LL)
               DO CC=CC_START,CC_END-1
                  IF(NEIGHBORS(CC) <= 0) CYCLE
                  CALL JKR_DIAG_PAIR_EVENT('neighbor_after_restore', &
                     IGLOBAL_ID(LL),IGLOBAL_ID(NEIGHBORS(CC)),-1,CC, &
                     0D0,0D0,(/0D0,0D0,0D0/),PFT_NEIGHBOR(:,CC), &
                     PFR_NEIGHBOR(:,CC),JKR_CONTACT_STATE(CC))
                  CALL JKR_ROTATION_HISTORY_EVENT('neighbor_after_restore', &
                     IGLOBAL_ID(LL),IGLOBAL_ID(NEIGHBORS(CC)),-1,CC,PFTW_NEIGHBOR(CC))
               ENDDO
            ENDDO
         ENDIF
      ENDIF


! resetting do_nsearch to false here since neighbor search will have
! just been invoked
      DO_NSEARCH = .FALSE.

      RETURN

   END SUBROUTINE NEIGHBOUR

END MODULE NEIGHBOUR_MOD
