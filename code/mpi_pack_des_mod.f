!----------------------------------------------------------------------!
!  Module: MPI_PACK_DES                                                !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for packing real and ghost particles     !
!     into the MPI send buffers.                                       !
!----------------------------------------------------------------------!
      MODULE MPI_PACK_DES

      PRIVATE
      PUBLIC :: DESMPI_PACK_PARCROSS, DESMPI_PACK_GHOSTPAR
      public :: pack_dbuf

      interface pack_dbuf
        module procedure pack_db0
        module procedure pack_db1
        module procedure pack_i0
        module procedure pack_i1
        module procedure pack_l0
      end interface pack_dbuf

      CONTAINS

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_PACK_GHOSTPAR                                    !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Packs ghost particle in the send buffer.                    !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_PACK_GHOSTPAR(PFACE)

! Global Variables:
!---------------------------------------------------------------------//
! Size of ghost particle packets
      use desmpi, only: iGhostPacketSize
! Particle information
      use discretelement
! Particle temperatures. current/previous
      use des_thermo, only: DES_T_s
! The MPI send buffer
      use desmpi, only: dSENDBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! The neighbor processor's rank
      use desmpi, only: iNEIGHPROC
! Function to convert DES grid IJK to new proc value.
      use desgrid, only: dg_ijkconv
! Size of the send buffer
      use desmpi, only: isendcnt
! Offset for particles with cycle BCs (otherwise zero)
      use desmpi, only: dcycl_offset
! Map of particles to DES grid
      use derived_types, only: DG_PIC
! Cell number of ghost particles
      use desmpi, only: iSENDINDICES
! Receive buf size
      use desmpi, only: iMAXBUF

      use functions, only: is_exiting
      use functions, only: is_ghost, is_entering_ghost, is_exiting_ghost
      use functions, only: is_any_ghost

      use vtk
      use compar, only: MyPE,istart1, iend1, jstart1, jend1, kstart1, kend1

      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: lijk,lindx,ltot_ind,lpicloc,lpar_cnt,lcurpar
      integer :: lbuf
      double precision :: cutoff_value, cutoff_sign
      integer ::  cutoff_coord
!......................................................................!

! Jeff Dietiker, Nov. 1st, 2021
! When the fluid grid spacing is less than half the des grid spacing
! at the process boundary, we can pack particles that will not
! belong to the neighbor process. The particles are deleted later on,
! but it is more efficient to not pack them in the first place.
! We use a generic test:
! cutoff_sign*des_pos_new(lcurpar,cutoff_coord) > cutoff_value
! so we do not need to test on pface inside the loop. This is equivalent to:
! North boundary: des_pos_new(lcurpar,2) < yn(jend1-2)
! South boundary: des_pos_new(lcurpar,2) > yn(jstart1+1)
! etc.

      SELECT CASE(PFACE)

         CASE(1) ! East
            cutoff_value = -xe(iend1-2)
            cutoff_sign  = -1.0D0
            cutoff_coord = 1

         CASE(2) ! West
            cutoff_value = xe(istart1+1)
            cutoff_sign  = 1.0D0
            cutoff_coord = 1

         CASE(3) ! North
            cutoff_value = -yn(jend1-2)
            cutoff_sign  = -1.0D0
            cutoff_coord = 2

         CASE(4) ! South
            cutoff_value = yn(jstart1+1)
            cutoff_sign  = 1.0D0
            cutoff_coord = 2

         CASE(5) ! Top
            cutoff_value = -zt(kend1-2)
            cutoff_sign  = -1.0D0
            cutoff_coord = 3

         CASE(6) ! Bottom
            cutoff_value = zt(kstart1+1)
            cutoff_sign  = 1.0D0
            cutoff_coord = 3

      END SELECT

      lpar_cnt = 0
      ltot_ind = isendindices(1,pface)
      do lindx = 2,ltot_ind+1
         lijk = isendindices(lindx,pface)
         do lpicloc =1,dg_pic(lijk)%isize
            lbuf = lpar_cnt*iGhostPacketSize + ibufoffset
            lcurpar = dg_pic(lijk)%p(lpicloc)

! Do not send particle data for a ghost particle whose owner has not yet
! updated the particle's data on this processor.
            if(is_any_ghost(lcurpar) .and. &
                (.not.ighost_updated(lcurpar))) cycle

! Do not pack particles that are not in the neighbor's ghost layers
          if(cutoff_sign*des_pos_new(lcurpar,cutoff_coord) > cutoff_value) cycle

! 1) Global ID
            call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
! 2) DES grid IJK
            call pack_dbuf(lbuf,dg_ijkconv(lijk,pface,                 &
               ineighproc(pface)),pface)
! 3) DES grid IJK - previous
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,  &
               ineighproc(pface)),pface)
! 4) Radius
            call pack_dbuf(lbuf,des_radius(lcurpar),pface)
! 5) Phase index
            call pack_dbuf(lbuf,pijk(lcurpar,5),pface)
! 6) Pmass for the use of new damping
            call pack_dbuf(lbuf,pmass(lcurpar),pface)
! 7) Position
            call pack_dbuf(lbuf,des_pos_new(lcurpar,:)+                &
               dcycl_offset(pface,:),pface)
! 8) Translational velocity
            call pack_dbuf(lbuf,des_vel_new(lcurpar,:),pface)
! 9) Rotational velocity
            call pack_dbuf(lbuf,omega_new(lcurpar,:),pface)
! 10) Exiting particle flag
            call pack_dbuf(lbuf,merge(1,0,is_exiting(lcurpar).or.is_exiting_ghost(lcurpar)),pface)
! 11) Temperature
            IF(ENERGY_EQ) &
               call pack_dbuf(lbuf,des_t_s(lcurpar),pface)
! 12) User variables
            IF(DES_USR_VAR_SIZE > 0) &
               call pack_dbuf(lbuf,des_usr_var(:,lcurpar),pface)
! 13) Parcel statwt
            IF(CGDEM) &
               call pack_dbuf(lbuf,des_cgp_stw(lcurpar),pface)
! 14) Real particle radius in parcel
            IF(CGDEM) &
               call pack_dbuf(lbuf,des_cgp_rpr(lcurpar),pface)
! 15) Reaction rates
            IF(SAVE_PART_RRATES)  call pack_dbuf(lbuf,Part_RRates_out(lcurpar,:),pface)

! 16) Residence time
            call pack_dbuf(lbuf,residence_time(lcurpar),pface)
! 17) SuperDEM
            IF(SuperDEM) THEN
               call pack_dbuf(lbuf,super_r(lcurpar,:),pface)
               call pack_dbuf(lbuf,super_mn(lcurpar,:),pface)
               call pack_dbuf(lbuf,super_q(lcurpar,:),pface)
            ENDIF
! 18) GluedSphereDEM
            IF(GSP_EXPLICIT) THEN
               call pack_dbuf(lbuf,gid_list(lcurpar),pface)
               call pack_dbuf(lbuf,gp_neighsid(lcurpar,:),pface)
               call pack_dbuf(lbuf,sc2gpc_vec(lcurpar,:),pface)
               IF(flexible_gsp) THEN
                  call pack_dbuf(lbuf, bfn(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bmn(lcurpar,:),pface)

                  call pack_dbuf(lbuf, bftx(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bfty(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bftz(lcurpar,:),pface)

                  call pack_dbuf(lbuf, bmtx(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bmty(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bmtz(lcurpar,:),pface)

                  call pack_dbuf(lbuf, bnx(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bny(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bnz(lcurpar,:),pface)
               ENDIF
               IF(ENERGY_EQ) THEN
                  call pack_dbuf(lbuf,gp_kps(lcurpar),pface)
               ENDIF
            ELSEIF(GSP_IMPLICIT) THEN
               call pack_dbuf(lbuf,glueSurface(lcurpar),pface)
               call pack_dbuf(lbuf,glueQuat(lcurpar,:),pface)
               call pack_dbuf(lbuf,glueDiameter(lcurpar),pface)
               call pack_dbuf(lbuf,glueBounding(lcurpar),pface)
            ENDIF

            lpar_cnt = lpar_cnt + 1
         end do
      end do
      dsendbuf(1+mod(pface,2))%facebuf(1)=lpar_cnt
      isendcnt(pface) = lpar_cnt*iGhostPacketSize + ibufoffset

! @renjieke, not all PE entering this subroutine, do not use global communication here
! remind user to increase resizing factor in GUI.
      if(gsp_explicit) then
         if(isendcnt(pface) > imaxbuf) then
            write(*,*) "Send buf size is larger than recv buf size MPI truncated errors will happen."
            write(*,*) "Please increase Buffer resizing factor."
         endif
      endif

      end subroutine desmpi_pack_ghostpar

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_PACK_PARCROSS                                    !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Packs real particle in the send buffer.                     !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_PACK_PARCROSS(PFACE)

! Global Variables:
!---------------------------------------------------------------------//
! The MPI send buffer
      use desmpi, only: dSENDBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! Runtime flag for solving species equations
      use run, only: ANY_SOLIDS_SPECIES_EQ
! Particle information
      use discretelement
! The neighbor processor's rank
      use desmpi, only: iNEIGHPROC
! The statistical weight of each particle.
      use mfix_pic, only: DES_STAT_WT
! Particle species composition
      use des_rxns, only: DES_X_s
! Particle temperatures. current/previous
      use des_thermo, only: DES_T_s
! Explicit convection and HOR
      use des_thermo, only: CONV_Qs, RXNS_Qs

      use desgrid, only: dg_ijkconv, icycoffset
      use desmpi, only: dcycl_offset, isendcnt
      use desmpi, only: irecvindices

      use desmpi, only: iParticlePacketSize
      use desmpi, only: iPairPacketSize
      use derived_types, only: dg_pic
      use functions
      use mpi_unpack_des, only: GET_PENDING_PAIR_COUNT, &
         GET_PENDING_PAIR_HISTORY, HAS_PENDING_PAIR_HISTORY

      use vtk

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: li, lj, lk
      integer :: ltot_ind,lindx,cc,cc_start,cc_end
      integer :: lneigh,lijk,&
                 lpicloc,lparcnt,lcurpar
      integer :: lbuf,num_neighborlists_to_send

      logical, allocatable, dimension(:) :: going_to_send

! Location in the buffer where the number of pair data is specified.
      integer :: num_neighborlists_send_buf_loc
      integer :: tmp_max_pip
      integer :: pair_owner, pair_neigh, pair_index, owner_gid, neigh_gid
      double precision :: pair_pft(3), pair_pfr(3), pair_twist, pair_sign
      double precision :: pair_frame(4)
      logical :: pair_state
!......................................................................!

! pack the particle crossing the boundary
      ltot_ind = irecvindices(1,pface)
      lparcnt = 0

      allocate(going_to_send(max_pip))
      going_to_send(:) = .false.

      do lindx = 2,ltot_ind+1
         lijk = irecvindices(lindx,pface)
         do lpicloc = 1,dg_pic(lijk)%isize
            lcurpar = dg_pic(lijk)%p(lpicloc)

! if ghost particle then cycle
            if(is_any_ghost(lcurpar)) cycle

            going_to_send(lcurpar) = .true.
            lbuf = lparcnt*iParticlePacketSize + ibufoffset

! 1) Global ID
            call pack_dbuf(lbuf,iglobal_id(lcurpar),pface)
! 2) DES Grid IJK
            call pack_dbuf(lbuf,dg_ijkconv(lijk,pface,                 &
               ineighproc(pface)),pface)
! 3) DES grid IJK - previous
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(lcurpar),pface,  &
               ineighproc(pface)),pface)
! 4) Radius
            call pack_dbuf(lbuf,des_radius(lcurpar),pface)
! 5) Fluid cell I index with cycle offset
            li = pijk(lcurpar,1) + icycoffset(pface,1)
            call pack_dbuf(lbuf,li,pface)
! 6) Fluid cell J index with cycle offset
            lj = pijk(lcurpar,2) + icycoffset(pface,2)
            call pack_dbuf(lbuf,lj,pface)
! 7) Fluid cell K index with cycle offset
            lk = pijk(lcurpar,3) + icycoffset(pface,3)
            call pack_dbuf(lbuf,lk,pface)
! 8) Fluid cell IJK on destination process
            call pack_dbuf(lbuf,funijk_proc(li,lj,lk,                  &
               ineighproc(pface)),pface)
! 9) Particle solids phase index
            call pack_dbuf(lbuf,pijk(lcurpar,5),pface)
! 10) Entering particle flag.
            call pack_dbuf(lbuf, is_entering(lcurpar).or.is_entering_ghost(lcurpar), pface)
! 11) Exiting particle flag.
            call pack_dbuf(lbuf, is_exiting(lcurpar).or.is_exiting_ghost(lcurpar), pface)
! 12) Density
            call pack_dbuf(lbuf,ro_sol(lcurpar),pface)
! 13) Volume
            call pack_dbuf(lbuf,pvol(lcurpar),pface)
! 14) Mass
            call pack_dbuf(lbuf,pmass(lcurpar),pface)
! 15) Moment of Inertia
            call pack_dbuf(lbuf,omoi(lcurpar),pface)
! 16) SuperDEM
            IF(SuperDEM .or. GSP_IMPLICIT) call pack_dbuf(lbuf,omoi3(lcurpar,:),pface)
! 17) Position with cyclic shift
            call pack_dbuf(lbuf,des_pos_new(lcurpar,:) +               &
               dcycl_offset(pface,:),pface)
! 18) Translational velocity
            call pack_dbuf(lbuf,des_vel_new(lcurpar,:),pface)
! 19) Rotational velocity
            call pack_dbuf(lbuf,omega_new(lcurpar,:),pface)
! 20) Accumulated translational forces
            call pack_dbuf(lbuf,fc(lcurpar,:),pface)
! 21) Accumulated torque forces
            call pack_dbuf(lbuf,tow(lcurpar,:),pface)

            IF(ENERGY_EQ) THEN
! 22) Temperature
               call pack_dbuf(lbuf,des_t_s(lcurpar),pface)
! 23) Species composition
               call pack_dbuf(lbuf,des_x_s(lcurpar,:),pface)
           ENDIF

! 24) User defined variable
            IF(DES_USR_VAR_SIZE > 0) call pack_dbuf(lbuf, des_usr_var(:,lcurpar),pface)
! 25) Particle orientation
            IF(PARTICLE_ORIENTATION) call pack_dbuf(lbuf,orientation(lcurpar,:),pface)
! -- Higher order integration variables
            IF (DO_OLD) THEN
! 26) Position (previous)
               call pack_dbuf(lbuf,des_pos_old(lcurpar,:) +            &
                  dcycl_offset(pface,:),pface)
! 27) Translational velocity (previous)
               call pack_dbuf(lbuf,des_vel_old(lcurpar,:),pface)
! 28) Rotational velocity (previous)
               call pack_dbuf(lbuf,omega_old(lcurpar,:),pface)
! 29) Translational acceleration (previous)
               call pack_dbuf(lbuf,des_acc_old(lcurpar,:),pface)
! 30) Rotational acceleration (previous)
               call pack_dbuf(lbuf,rot_acc_old(lcurpar,:),pface)
            ENDIF

            IF(DES_EXPLICITLY_COUPLED) THEN
! 31) Explicit drag force
               call pack_dbuf(lbuf, drag_fc(lcurpar,:),pface)
! 32) Explicit convective heat transfer
               IF(ENERGY_EQ) call pack_dbuf(lbuf,conv_qs(lcurpar),pface)
! 33) Explicit heat of reaction
               IF(ANY_SOLIDS_SPECIES_EQ) call pack_dbuf(lbuf, rxns_qs(lcurpar),pface)
           ENDIF
! 34) Parcel statwt
            IF(CGDEM) &
               call pack_dbuf(lbuf,des_cgp_stw(lcurpar),pface)
! 35) Real particle radius in parcel
            IF(CGDEM) &
               call pack_dbuf(lbuf,des_cgp_rpr(lcurpar),pface)
! 36) Reaction rates
            IF(SAVE_PART_RRATES)  call pack_dbuf(lbuf,Part_RRates_out(lcurpar,:),pface)
! 37) Residence time
            call pack_dbuf(lbuf,RESIDENCE_TIME(lcurpar),pface)
! 37a) Spherical JKR wall histories attached to this particle owner.
            IF(USE_JKR) THEN
               DO CC=1,COLLISION_ARRAY_MAX
                  CALL PACK_DBUF(lbuf,WALL_COLLISION_FACET_ID(CC,lcurpar),pface)
                  CALL PACK_DBUF(lbuf,WALL_COLLISION_PFT(:,CC,lcurpar),pface)
                  CALL PACK_DBUF(lbuf,WALL_COLLISION_PFR(:,CC,lcurpar),pface)
                  CALL PACK_DBUF(lbuf,WALL_COLLISION_PFTW(CC,lcurpar),pface)
                  CALL PACK_DBUF(lbuf,WALL_COLLISION_FRAME(:,CC,lcurpar),pface)
                  CALL PACK_DBUF(lbuf,JKR_WALL_CONTACT_STATE(CC,lcurpar),pface)
               ENDDO
            ENDIF
! 38) SuperDEM
            IF(SuperDEM) THEN
               call pack_dbuf(lbuf,super_r(lcurpar,:),pface)
               call pack_dbuf(lbuf,super_mn(lcurpar,:),pface)
               call pack_dbuf(lbuf,super_q(lcurpar,:),pface)
            ENDIF
! 38a) GluedSphereDEM
! renjieke 3-7-2024!
! when used those variables, particle_state is checked for only normal case, so do not need to pack them in ghost
! but for safety reason, test whether packing them in desmpi_pack_ghostpar will change the result, logically it should not change
! if really need to packing them in ghost, remember to unpacking ghost as well.
            IF(GSP_EXPLICIT) THEN
               call pack_dbuf(lbuf, gp_sa(lcurpar),pface)
               call pack_dbuf(lbuf, gid_list(lcurpar),pface)
               call pack_dbuf(lbuf, gp_neighsid(lcurpar,:),pface)
               call pack_dbuf(lbuf, sc2gpc_vec(lcurpar,:),pface)
               IF(flexible_gsp) THEN
                  call pack_dbuf(lbuf, bfn(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bmn(lcurpar,:),pface)

                  call pack_dbuf(lbuf, bftx(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bfty(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bftz(lcurpar,:),pface)

                  call pack_dbuf(lbuf, bmtx(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bmty(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bmtz(lcurpar,:),pface)

                  call pack_dbuf(lbuf, bnx(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bny(lcurpar,:),pface)
                  call pack_dbuf(lbuf, bnz(lcurpar,:),pface)
               ENDIF
               IF(ENERGY_EQ) THEN
                  call pack_dbuf(lbuf, gp_kps(lcurpar),pface)
                  call pack_dbuf(lbuf, gp_neighsa(lcurpar,:),pface)
               ENDIF
            ELSEIF(GSP_IMPLICIT) THEN
               call pack_dbuf(lbuf, glueSurface(lcurpar),pface)
               call pack_dbuf(lbuf, glueQuat(lcurpar,:),pface)
               call pack_dbuf(lbuf, glueDiameter(lcurpar),pface)
               call pack_dbuf(lbuf, glueBounding(lcurpar),pface)
            ENDIF
! DEM particles are converted to ghost particles. This action does not
! change the number of particles.
               if (is_entering(lcurpar)) then
                  call set_entering_ghost(lcurpar)
               elseif (is_exiting(lcurpar)) then
                  call set_exiting_ghost(lcurpar)
               elseif (is_rigid_motion(lcurpar)) then
                  call set_rigid_motion_ghost(lcurpar)
               else
                  call set_ghost(lcurpar)
               endif
               ighost_cnt = ighost_cnt + 1
! Clear out the force array.
            fc(lcurpar,:) = 0.
            lparcnt = lparcnt + 1
         end do
      end do

! 39a) Number of pair datasets
! Calculate the location in buffer where the number of pair data is
! stored and skip specifying the entry. After all the pair data is
! packed, then this value is set.

! @renjieke, if a particle is being packed and sended,
! we also packed all pairs of data associated with this particle and send them all
      lbuf = lparcnt*iParticlePacketSize + ibufoffset
      num_neighborlists_send_buf_loc = lbuf
      lbuf = lbuf+1
      num_neighborlists_to_send = 0

      tmp_max_pip = 0
      IF(USE_JKR .AND. lparcnt > 0) tmp_max_pip = max_pip
      ! if(flexible_gsp) tmp_max_pip = max_pip

      do lcurpar = 1, tmp_max_pip
! A stored pair can be incident on either migrating endpoint.
         if (.not.USE_JKR .and. .not.going_to_send(lcurpar)) cycle
         cc_start = 1
         if(lcurpar .gt. 1) cc_start = neighbor_index(lcurpar-1)
         cc_end = neighbor_index(lcurpar)
         do cc = cc_start, cc_end-1
            IF (0 .eq. NEIGHBORS(cc)) EXIT
! Do not send pairing data if the pair no longer exists or if the
! particle is exiting as it may be locatable during unpacking.
            lneigh = neighbors(cc)
            if(is_nonexistent(lneigh)) cycle
            if(is_exiting(lneigh)) cycle
            IF(USE_JKR) THEN
               IF(.NOT.(going_to_send(lcurpar) .OR. &
                  going_to_send(lneigh))) CYCLE
! Received history is newer than any pre-exchange CC at a corner rank.
               IF(HAS_PENDING_PAIR_HISTORY(iglobal_id(lcurpar), &
                  iglobal_id(lneigh))) CYCLE
            ENDIF
            pair_owner = lcurpar
            pair_neigh = lneigh
            pair_sign = 1D0
            IF(USE_JKR .AND. .NOT.going_to_send(lcurpar)) THEN
               pair_owner = lneigh
               pair_neigh = lcurpar
               pair_sign = -1D0
            ENDIF

            ! @renjieke, one thing could happens, that lcurpar and lneigh might not be
            ! within the same gsp. For those pair, we still send the virtual bond info,
            ! however, we should expect zero values there
! 40) Global ID of particle being packed.
            call pack_dbuf(lbuf,iglobal_id(pair_owner),pface)
! 41) DES grid IJK of cell receiving the particle.
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(pair_owner),pface,  &
               ineighproc(pface)),pface)
! 42) Global ID of neighbor particle.
            call pack_dbuf(lbuf,iglobal_id(pair_neigh),pface)
! 43) DES grid IJK of cell containing the neighbor particle.
            call pack_dbuf(lbuf,dg_ijkconv(dg_pijkprv(pair_neigh),pface,  &
               ineighproc(pface)),pface)
! 44) Spherical JKR histories: directed PFT, pair-invariant PFR/state, and twist.
            IF(USE_JKR) THEN
               CALL PACK_DBUF(lbuf,pair_sign*PFT_NEIGHBOR(:,CC),pface)
               CALL PACK_DBUF(lbuf,PFR_NEIGHBOR(:,CC),pface)
               CALL PACK_DBUF(lbuf,JKR_CONTACT_STATE(CC),pface)
               CALL PACK_DBUF(lbuf,PFTW_NEIGHBOR(CC),pface)
               CALL PACK_DBUF(lbuf,pair_sign*JKR_FRAME_NEIGHBOR(1:3,CC),pface)
               CALL PACK_DBUF(lbuf,JKR_FRAME_NEIGHBOR(4,CC),pface)
            ENDIF
! 45) Superquadric contact points
            IF(SuperDEM) THEN
               call  pack_dbuf(lbuf,contact_point_A(:,CC),pface)
               call  pack_dbuf(lbuf,contact_point_B(:,CC),pface)

               call  pack_dbuf(lbuf,contact_lambda_A(CC),pface)
               call  pack_dbuf(lbuf,contact_lambda_B(CC),pface)
            ENDIF

! Increment the number of pairs being sent.
            num_neighborlists_to_send = num_neighborlists_to_send + 1
         enddo
      enddo

! Forward histories received on an earlier migration face. Their CC
! slots do not exist until NEIGHBOUR, so use the same pending queue and
! native pair packet instead of creating a parallel contact table.
      IF(USE_JKR .AND. lparcnt > 0) THEN
         DO pair_index=1,GET_PENDING_PAIR_COUNT()
            CALL GET_PENDING_PAIR_HISTORY(pair_index,owner_gid,neigh_gid, &
               pair_pft,pair_pfr,pair_state,pair_twist,pair_frame)
            pair_owner = 0
            DO lcurpar=1,max_pip
               IF(.NOT.going_to_send(lcurpar)) CYCLE
               IF(iglobal_id(lcurpar) == owner_gid .OR. &
                  iglobal_id(lcurpar) == neigh_gid) THEN
                  pair_owner = lcurpar
                  EXIT
               ENDIF
            ENDDO
            IF(pair_owner == 0) CYCLE
! JKR unpack resolves IDs after all faces; no intermediate cell lookup.
            CALL PACK_DBUF(lbuf,owner_gid,pface)
            CALL PACK_DBUF(lbuf,0,pface)
            CALL PACK_DBUF(lbuf,neigh_gid,pface)
            CALL PACK_DBUF(lbuf,0,pface)
            CALL PACK_DBUF(lbuf,pair_pft,pface)
            CALL PACK_DBUF(lbuf,pair_pfr,pface)
            CALL PACK_DBUF(lbuf,pair_state,pface)
            CALL PACK_DBUF(lbuf,pair_twist,pface)
            CALL PACK_DBUF(lbuf,pair_frame,pface)
            num_neighborlists_to_send = num_neighborlists_to_send + 1
         ENDDO
      ENDIF

! Store the number of pair datasets being sent. This information is
! stored before the pairing data so the receiving process knows the
! amount of data to unpack.
      lbuf = num_neighborlists_send_buf_loc
! 39b) Number of pair datasets.
      call pack_dbuf(lbuf,num_neighborlists_to_send,pface)

      dsendbuf(1+mod(pface,2))%facebuf(1) = lparcnt
      isendcnt(pface) = lparcnt*iParticlePacketSize + &
         num_neighborlists_to_send*iPairPacketSize + ibufoffset + 3
      deallocate(going_to_send)

      RETURN
      END SUBROUTINE DESMPI_PACK_PARCROSS

!----------------------------------------------------------------------!
! PACK SUBROUTINE FOR SINGLE REAL VARIABLES                            !
!----------------------------------------------------------------------!
      subroutine pack_db0(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata

      dsendbuf(1+mod(pface,2))%facebuf(lbuf) = idata
      lbuf = lbuf + 1

      return
      end subroutine pack_db0

!----------------------------------------------------------------------!
! Pack subroutine for real arrays                                      !
!----------------------------------------------------------------------!
      subroutine pack_db1(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      dsendbuf(1+mod(pface,2))%facebuf(lbuf:lbuf+lsize-1) = idata
      lbuf = lbuf + lsize

      return
      end subroutine pack_db1

!----------------------------------------------------------------------!
! Pack subroutine for single integer variables                         !
!----------------------------------------------------------------------!
      subroutine pack_i0(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(in) :: idata

      dsendbuf(1+mod(pface,2))%facebuf(lbuf) = idata
      lbuf = lbuf + 1

      return
      end subroutine pack_i0

!----------------------------------------------------------------------!
! Pack subroutine for integer arrays                                   !
!----------------------------------------------------------------------!
      subroutine pack_i1(lbuf,idata,pface)
      use desmpi, only: dSENDBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(in) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      dsendbuf(1+mod(pface,2))%facebuf(lbuf:lbuf+lsize-1) = idata
      lbuf = lbuf + lsize

      return
      end subroutine pack_i1

!----------------------------------------------------------------------!
! Pack subroutine for logical scalars                                  !
!----------------------------------------------------------------------!
      subroutine pack_l0(lbuf,ldata,pface)
      use desmpi, only: dSENDBUF

      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      logical, intent(in) :: ldata

      dsendbuf(1+mod(pface,2))%facebuf(lbuf) = merge(1.0, 0.0, ldata)
      lbuf = lbuf + 1

      return
      end subroutine pack_l0

      end module mpi_pack_des
