!----------------------------------------------------------------------!
!  Module: MPI_UNPACK_DES                                              !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
!  Purpose: Contains routines for unpacking real and ghost particles   !
!     from the MPI recv buffers.                                       !
!----------------------------------------------------------------------!
      MODULE MPI_UNPACK_DES

      PRIVATE

      PUBLIC :: QUEUE_PAIR_HISTORY
      PUBLIC :: DESMPI_UNPACK_PARCROSS, DESMPI_UNPACK_GHOSTPAR
      PUBLIC :: unpack_dbuf
      PUBLIC :: RESTORE_PENDING_PAIR_HISTORY
      PUBLIC :: JKR_DIAG_PAIR_EVENT

      INTEGER :: PENDING_PAIR_COUNT = 0
      INTEGER :: PENDING_PAIR_CAPACITY = 0
      INTEGER, ALLOCATABLE :: PENDING_PAIR_OWNER_GID(:)
      INTEGER, ALLOCATABLE :: PENDING_PAIR_NEIGH_GID(:)
      DOUBLE PRECISION, ALLOCATABLE :: PENDING_PAIR_PFT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PENDING_PAIR_PFR(:,:)
      LOGICAL, ALLOCATABLE :: PENDING_PAIR_STATE(:)
      LOGICAL :: JKR_DIAG_CHECKED = .FALSE.
      LOGICAL :: JKR_DIAG_ENABLED = .FALSE.
      INTEGER :: JKR_DIAG_UNIT = -1

      interface unpack_dbuf
         module procedure unpack_db0 ! real scalars
         module procedure unpack_db1 ! real arrays
         module procedure unpack_i0  ! integer scalars
         module procedure unpack_i1  ! integer arrays
         module procedure unpack_l0  ! logical scalars
      end interface unpack_dbuf

      CONTAINS

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_UNPACK_GHOSTPAR                                  !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Unpacks ghost particle from the recv buffer.                !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_UNPACK_GHOSTPAR(pface)

! Global Variables:
!---------------------------------------------------------------------//
! Size of Particle data packet
      use desmpi, only: iGhostPacketSize
! Index of last particle added to this process.
      use desmpi, only: iSPOT
! Particle information
      use discretelement
! The MPI receive buffer
      use desmpi, only: dRECVBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! Dimensions of DES grid
      use desgrid, only: DG_IJKSIZE2
! Particle temperatures
      use des_thermo, only: DES_T_s

      use des_allocate
! Reactions
      use des_rxns, only:  Part_RRates_out, save_part_rrates
! Global Constants:
!---------------------------------------------------------------------//
      use constant, only: PI
      use compar, only: pe_io, mype, numPEs

      use functions, only: is_nonexistent
      use functions, only: is_normal,  set_normal
      use functions, only: is_exiting, set_exiting, set_exiting_ghost
      use functions, only: set_ghost, set_rigid_motion_ghost
      use sq_properties_mod
      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: lcurpar,lparid,lprvijk,lparijk,lparcnt
      integer :: lbuf,llocpar,lnewcnt
      logical,dimension(:),allocatable :: lfound
      integer,dimension(:),allocatable :: lnewspot,lnewpic
      logical :: tmp
      DOUBLE PRECISION :: axi(3), m,n, SVOLUME
!......................................................................!

! unpack the particles:
! if it already exists update the position
! if not and do_nsearch is true then add to the particle array

      lparcnt = drecvbuf(1+mod(pface,2))%facebuf(1)
      lnewcnt = lparcnt
      allocate (lfound(lparcnt),lnewspot(lparcnt),lnewpic(dg_ijksize2))
      lfound(:) = .false.
      lnewspot(:) =0
      lnewpic = 0

      do lcurpar = 1,lparcnt
         lbuf = (lcurpar-1)*iGhostPacketSize+ibufoffset

! 1) Global ID
         call unpack_dbuf(lbuf,lparid,pface)
! 2) DES grid IJK
         call unpack_dbuf(lbuf,lparijk,pface)
! 3) DES grid IJK - Previous
         call unpack_dbuf(lbuf,lprvijk,pface)

! Determine if this particle already exists on this process as a
! ghost particle. If so, (lfound), then the current information is
! updated on the current process. Otherwise (.NOT.lfound) a new
! ghost particle is created on this process.
         lfound(lcurpar) = locate_par(lparid,lprvijk,llocpar)
         if (lparijk .ne. lprvijk .and. .not.lfound(lcurpar)) then
            lfound(lcurpar) = locate_par(lparid,lparijk,llocpar)
         endif

         if(lfound(lcurpar)) then
! Store the local variables
            dg_pijk(llocpar) = lparijk
            dg_pijkprv(llocpar) = lprvijk

! 4) Radius
            call unpack_dbuf(lbuf,des_radius(llocpar),pface)
! 5) Phase index
            call unpack_dbuf(lbuf,pijk(llocpar,5),pface)
            if(des_rigid_motion(pijk(llocpar,5))) call set_rigid_motion_ghost(llocpar)
! 6) Particle mass for the use of the new damping method
            call unpack_dbuf(lbuf,pmass(llocpar),pface)
! 7) Position
            call unpack_dbuf(lbuf,des_pos_new(llocpar,1:dimn),pface)
! 8) Translational velocity
            call unpack_dbuf(lbuf,des_vel_new(llocpar,1:dimn),pface)
! 9) Rotational velocity
            call unpack_dbuf(lbuf,omega_new(llocpar,1:3),pface)
! 10) Exiting particle flag
            call unpack_dbuf(lbuf,tmp,pface)
            if (tmp) call set_exiting_ghost(llocpar)
! 11) Temperature
            IF(ENERGY_EQ) &
               call unpack_dbuf(lbuf,des_t_s(llocpar),pface)
! 12) User variables
            IF(DES_USR_VAR_SIZE > 0) &
               call unpack_dbuf(lbuf,des_usr_var(:,llocpar),pface)
! 13) Parcel statwt
            IF(CGDEM) &
               call unpack_dbuf(lbuf,des_cgp_stw(llocpar),pface)
! 14) Real particle radius in parcel
            IF(CGDEM) &
               call unpack_dbuf(lbuf,des_cgp_rpr(llocpar),pface)
! 15) Reaction rates
            IF(SAVE_PART_RRATES) call unpack_dbuf(lbuf,Part_RRates_out(llocpar,:),pface)
! 16) Residence time
            call unpack_dbuf(lbuf,residence_time(llocpar),pface)

! 17) SuperDEM
            IF(SuperDEM) THEN
               call unpack_dbuf(lbuf,super_r(llocpar,1:3),pface)
               call unpack_dbuf(lbuf,super_mn(llocpar,1:2),pface)
               call unpack_dbuf(lbuf,super_q(llocpar,1:4),pface)
            ENDIF
! Calculate the volume of the ghost particle.
            IF (SuperDEM) then
               axi(1)=super_r(llocpar,1)
               axi(2)=super_r(llocpar,2)
               axi(3)=super_r(llocpar,3)
               m=super_mn(llocpar,1) ! m
               n=super_mn(llocpar,2) ! n
               CALL SQ_VOLUME (axi,m,n,svolume)
               PVOL(llocpar)=svolume
            else
               PVOL(llocpar) = (4.0D0/3.0D0)*PI*DES_RADIUS(llocpar)**3
            endif
! GluedSpherDEM
            IF(GSP_EXPLICIT) THEN
               call unpack_dbuf(lbuf,gid_list(llocpar),pface)
               ! @renjieke, in flexible GSP, this value is never changed
               call unpack_dbuf(lbuf,gp_neighsid(llocpar,:),pface)
               call unpack_dbuf(lbuf,sc2gpc_vec(llocpar,:),pface)
               if(flexible_gsp) then
                  call unpack_dbuf(lbuf, bfn(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bmn(llocpar,:),pface)

                  call unpack_dbuf(lbuf, bftx(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bfty(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bftz(llocpar,:),pface)

                  call unpack_dbuf(lbuf, bmtx(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bmty(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bmtz(llocpar,:),pface)

                  call unpack_dbuf(lbuf, bnx(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bny(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bnz(llocpar,:),pface)
               endif
               IF(ENERGY_EQ) THEN
                  call unpack_dbuf(lbuf,gp_kps(llocpar),pface)
               ENDIF
            ELSEIF(GSP_IMPLICIT) THEN
               call unpack_dbuf(lbuf,glueSurface(llocpar),pface)
               call unpack_dbuf(lbuf,glueQuat(llocpar,:),pface)
               call unpack_dbuf(lbuf,glueDiameter(llocpar),pface)
               call unpack_dbuf(lbuf,glueBounding(llocpar),pface)
            ENDIF
! Flag that the ghost particle was updated.
            ighost_updated(llocpar) = .true.
            lnewcnt = lnewcnt-1

! Copy the current value to the previous value if needed.
            IF (DO_OLD) THEN
               des_pos_old(llocpar,:)= des_pos_new(llocpar,:)
               des_vel_old(llocpar,:)= des_vel_new(llocpar,:)
               omega_old(llocpar,:)= omega_new(llocpar,:)
            ENDIF

         else
            lnewpic(lparijk) = lnewpic(lparijk) + 1
         endif
      enddo

! iAdd new ghost particles
      if(lnewcnt > 0) then
         call PARTICLE_GROW(pip+lnewcnt)
         ighost_cnt = ighost_cnt + lnewcnt
         pip = pip + lnewcnt
         max_pip = max(pip,max_pip)
         do lcurpar = 1,lparcnt
            if(lfound(lcurpar)) cycle
            lbuf = (lcurpar-1)*iGhostPacketSize+ibufoffset

!  1) Global particle ID
            call unpack_dbuf(lbuf,lparid,pface)
!  2) DES grid IJK
            call unpack_dbuf(lbuf,lparijk,pface)
!  3) DES grid IJK - Previous
            call unpack_dbuf(lbuf,lprvijk,pface)
! Locate the first open space in the particle array.
            do while(.not.is_nonexistent(ispot))
               ispot = ispot + 1
            enddo

! Set the flags for the ghost particle and store the local variables.
            call set_ghost(ispot)
            iglobal_id(ispot)  = lparid
            dg_pijk(ispot) = lparijk
            dg_pijkprv(ispot) = lprvijk
!  4) Radius
            call unpack_dbuf(lbuf,des_radius(ispot),pface)
!  5) Phase index
            call unpack_dbuf(lbuf,pijk(ispot,5),pface)
            if(des_rigid_motion(pijk(ispot,5))) call set_rigid_motion_ghost(ispot)
!  6) Particle mass
            call unpack_dbuf(lbuf,pmass(ispot),pface)
!  7) Position
            call unpack_dbuf(lbuf,des_pos_new(ispot,1:dimn),pface)
!  8) Translational velocity
            call unpack_dbuf(lbuf,des_vel_new(ispot,1:dimn),pface)
!  9) Rotational velocity
            call unpack_dbuf(lbuf,omega_new(ispot,1:dimn),pface)
!  10) Exiting particle flag
            call unpack_dbuf(lbuf,tmp,pface)
            if (tmp) call set_exiting_ghost(ispot)
! 11) Temperature.
            IF(ENERGY_EQ) &
               call unpack_dbuf(lbuf,des_t_s(ispot),pface)
! 12) User variable
            IF(DES_USR_VAR_SIZE > 0)&
               call unpack_dbuf(lbuf,des_usr_var(:,ispot),pface)
! 13) Parcel statwt
            IF(CGDEM) &
               call unpack_dbuf(lbuf,des_cgp_stw(ispot),pface)
! 14) Real particle radius in parcel
            IF(CGDEM) &
               call unpack_dbuf(lbuf,des_cgp_rpr(ispot),pface)
! 15) Reaction rates
            IF(SAVE_PART_RRATES)  call unpack_dbuf(lbuf,Part_RRates_out(ispot,:),pface)
! 16) Residence time
            call unpack_dbuf(lbuf,residence_time(ispot),pface)

! 17) SuperDEM
            IF(SuperDEM) THEN
               call unpack_dbuf(lbuf,super_r(ispot,1:3),pface)
               call unpack_dbuf(lbuf,super_mn(ispot,1:2),pface)
               call unpack_dbuf(lbuf,super_q(ispot,1:4),pface)
            ENDIF
! GluedSphereDEM
            IF(GSP_EXPLICIT) THEN
               call unpack_dbuf(lbuf,gid_list(ispot),pface)
               call unpack_dbuf(lbuf,gp_neighsid(ispot,:),pface)
               call unpack_dbuf(lbuf,sc2gpc_vec(ispot,:),pface)
               if(flexible_gsp) then
                  call unpack_dbuf(lbuf, bfn(ispot,:),pface)
                  call unpack_dbuf(lbuf, bmn(ispot,:),pface)

                  call unpack_dbuf(lbuf, bftx(ispot,:),pface)
                  call unpack_dbuf(lbuf, bfty(ispot,:),pface)
                  call unpack_dbuf(lbuf, bftz(ispot,:),pface)

                  call unpack_dbuf(lbuf, bmtx(ispot,:),pface)
                  call unpack_dbuf(lbuf, bmty(ispot,:),pface)
                  call unpack_dbuf(lbuf, bmtz(ispot,:),pface)

                  call unpack_dbuf(lbuf, bnx(ispot,:),pface)
                  call unpack_dbuf(lbuf, bny(ispot,:),pface)
                  call unpack_dbuf(lbuf, bnz(ispot,:),pface)
               endif
               IF(ENERGY_EQ) THEN
                  call unpack_dbuf(lbuf,gp_kps(ispot),pface)
               ENDIF
            ELSEIF(GSP_IMPLICIT) THEN
               call unpack_dbuf(lbuf,glueSurface(ispot),pface)
               call unpack_dbuf(lbuf,glueQuat(ispot,:),pface)
               call unpack_dbuf(lbuf,glueDiameter(ispot),pface)
               call unpack_dbuf(lbuf,glueBounding(ispot),pface)
            ENDIF

            ighost_updated(ispot) = .true.
            lnewspot(lcurpar) = ispot

            IF(SuperDEM)   then
               axi(1)=super_r(ispot,1)
               axi(2)=super_r(ispot,2)
               axi(3)=super_r(ispot,3)
               m=super_mn(ispot,1)
               n=super_mn(ispot,2)
               CALL SQ_VOLUME (axi,m,n,svolume)
               PVOL(ispot)=svolume
            else
               PVOL(ispot) = (4.0D0/3.0D0)*PI*DES_RADIUS(ispot)**3
            endif
            IF (DO_OLD) THEN
               des_pos_old(ispot,1:dimn) = des_pos_new(ispot,1:dimn)
               des_vel_old(ispot,1:dimn) = des_vel_new(ispot,1:dimn)
               omega_old(ispot,1:3) = omega_new(ispot,1:3)
            ENDIF
         enddo
      endif

!deallocate temporary variablies
      deallocate (lfound,lnewspot,lnewpic)

      end subroutine desmpi_unpack_ghostpar

!----------------------------------------------------------------------!
!  Subroutine: DESMPI_UNPACK_PARCROSS                                  !
!  Author: Pradeep Gopalakrishnan                                      !
!                                                                      !
! Purpose: Unpacks real particle from the recv buffer.                 !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_UNPACK_PARCROSS(pface)

! Global Variables:
!---------------------------------------------------------------------//
! Size of ghost particle data packet
      use desmpi, only: iParticlePacketSize
! Index of last particle added to this process.
      use desmpi, only: iSPOT
! The MPI receive buffer
      use desmpi, only: dRECVBUF
! Buffer offset
      use desmpi, only: iBUFOFFSET
! Runtime flag for solving the energy equations
      use run, only: ENERGY_EQ
! Runtime flag for solving species equations
      use run, only: ANY_SOLIDS_SPECIES_EQ
! Particle information
      use discretelement
! Particle species composition and reaction rates
      use des_rxns, only: DES_X_s, Save_Part_RRates, Part_RRates_out
! Particle temperatures.
      use des_thermo, only: DES_T_s
! Explicit convection and HOR
      use des_thermo, only: CONV_Qs, RXNS_Qs
! The neighbor processor's rank
      use desmpi, only: iNEIGHPROC
! The ID of the current process
      use compar, only: myPE

! Module Procedures:
!---------------------------------------------------------------------//
      use des_allocate
      use desmpi_wrapper, only: DES_MPI_STOP
      use functions, only: IS_NORMAL, IS_NONEXISTENT
      use functions, only: SET_ENTERING, SET_EXITING, SET_NORMAL
      use functions, only: set_rigid_motion

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Processor boundary being packed (Top/Bottom/North/South/East/West)
      INTEGER, INTENT(IN) :: PFACE

! Local variables
!---------------------------------------------------------------------//
      integer :: lcurpar,lparcnt,llocpar,lparid,lparijk,lprvijk
      integer :: lneigh,&
                 lneighid,lneighijk
      logical :: lfound
      integer :: lbuf
      integer :: cc,num_neighborlists_sent,nn,ll

      logical :: tmp
      DOUBLE PRECISION :: PAIR_PFT_IN(3), PAIR_PFR_IN(3)
      LOGICAL :: PAIR_STATE_IN
!......................................................................!

! loop through particles and locate them and make changes
      lparcnt = drecvbuf(1+mod(pface,2))%facebuf(1)
! make sure enough space available
      call PARTICLE_GROW(pip + lparcnt)
      max_pip = max(pip+lparcnt,max_pip)

      do lcurpar =1,lparcnt

         lfound = .false.
         lbuf = (lcurpar-1)*iParticlePacketSize + ibufoffset
! 1) Global ID
         call unpack_dbuf(lbuf,lparid,pface)
! 2) DES Grid IJK
         call unpack_dbuf(lbuf,lparijk,pface)
! 3) DES grid IJK - previous
         call unpack_dbuf(lbuf,lprvijk,pface)

! A DEM particle should already exist on the current processor as a
! ghost particle. Match the sent particle to the local ghost particle
! by matching the global IDs. Decrement the iGHOST_CNT counter to
! account for the switch from ghost to real particle.
         lFOUND  = LOCATE_PAR(lPARID,lPRVIJK,lLOCPAR)
         IF (.NOT. lFOUND) THEN
            lFOUND = exten_locate_par(lPARID, lPARIJK, lLOCPAR)
            IF(.NOT.lFOUND) THEN
               WRITE(*,1000) iNEIGHPROC(PFACE), MYPE, lPARID
               CALL DES_MPI_STOP
            ENDIF
         ENDIF
         iGHOST_CNT = iGHOST_CNT - 1

 1000 FORMAT(2/1X,72('*'),/1x,'From: DESMPI_UNPACK_PARCROSS: ',/       &
         ' Error 1000: Unable to match particles crossing processor ', &
         'boundaries.',/3x,'Source Proc: ',I9,' ---> Destination ',    &
         'Proc: ', I9,/3x,'Global Particle ID: ',I12,/1x,72('*'))

! convert the local particle from ghost to existing and update its position
         call set_normal(llocpar)
         dg_pijk(llocpar) = lparijk
         dg_pijkprv(llocpar) = lprvijk
! 4) Radius
         call unpack_dbuf(lbuf,des_radius(llocpar),pface)
! 5-9) Fluid cell I,J,K,IJK, and solids phase index
         call unpack_dbuf(lbuf,pijk(llocpar,:),pface)
         if(des_rigid_motion(pijk(llocpar,5))) call set_rigid_motion(llocpar)
! 10) Entering particle flag.
         call unpack_dbuf(lbuf,tmp,pface)
         if (tmp) CALL SET_ENTERING(llocpar)
! 11) Exiting particle flag.
         call unpack_dbuf(lbuf,tmp,pface)
         if (tmp) CALL SET_EXITING(llocpar)
! 12) Density
         call unpack_dbuf(lbuf,ro_sol(llocpar),pface)
! 13) Volume
         call unpack_dbuf(lbuf,pvol(llocpar),pface)
! 14) Mass
         call unpack_dbuf(lbuf,pmass(llocpar),pface)
! 15) Moment of Inertia
         call unpack_dbuf(lbuf,omoi(llocpar),pface)
! 16) SuperDEM or GSP_IMPLICIT
         IF(SuperDEM .or. GSP_IMPLICIT) call unpack_dbuf(lbuf,omoi3(llocpar,1:3),pface)
! 17) Position with cyclic shift
         call unpack_dbuf(lbuf,des_pos_new(llocpar,:),pface)
! 18 Translational velocity
         call unpack_dbuf(lbuf,des_vel_new(llocpar,:),pface)
! 19) Rotational velocity
         call unpack_dbuf(lbuf,omega_new(llocpar,:),pface)
! 20) Accumulated translational forces
         call unpack_dbuf(lbuf,fc(llocpar,:),pface)
! 21) Accumulated torque forces
         call unpack_dbuf(lbuf,tow(llocpar,:),pface)

         IF(ENERGY_EQ) THEN
! 22) Temperature
            call unpack_dbuf(lbuf,des_t_s(llocpar),pface)
! 23) Species composition
            call unpack_dbuf(lbuf,des_x_s(llocpar,:),pface)
        ENDIF

! 24) User defined variable
         IF(DES_USR_VAR_SIZE > 0) call unpack_dbuf(lbuf,des_usr_var(:,llocpar),pface)
! 25) Particle orientation
         IF(PARTICLE_ORIENTATION) call unpack_dbuf(lbuf,orientation(llocpar,:),pface)

! -- Higher order integration variables
         IF (DO_OLD) THEN
! 26) Position (previous)
            call unpack_dbuf(lbuf,des_pos_old(llocpar,:),pface)
! 27) Translational velocity (previous)
            call unpack_dbuf(lbuf,des_vel_old(llocpar,:),pface)
! 28) Rotational velocity (previous)
            call unpack_dbuf(lbuf,omega_old(llocpar,:),pface)
! 29) Translational acceleration (previous)
            call unpack_dbuf(lbuf,des_acc_old(llocpar,:),pface)
! 30) Rotational acceleration (previous)
            call unpack_dbuf(lbuf,rot_acc_old(llocpar,:),pface)
         ENDIF

         IF(DES_EXPLICITLY_COUPLED) THEN
! 31) Explicit drag force
            call unpack_dbuf(lbuf,drag_fc(llocpar,:),pface)
! 32) Explicit convective heat transfer
            IF(ENERGY_EQ) call unpack_dbuf(lbuf,conv_qs(llocpar),pface)
! 33) Explicit heat of reaction
            IF(ANY_SOLIDS_SPECIES_EQ) call unpack_dbuf(lbuf,rxns_qs(llocpar),pface)
         ENDIF
! 34) Parcel statwt
         IF(CGDEM) call unpack_dbuf(lbuf,des_cgp_stw(llocpar),pface)
! 35) Real particle radius in parcel
         IF(CGDEM) call unpack_dbuf(lbuf,des_cgp_rpr(llocpar),pface)
! 36) Reaction rates
         IF(SAVE_PART_RRATES)  call unpack_dbuf(lbuf,Part_RRates_out(llocpar,:),pface)
! 37) Residence time
         call unpack_dbuf(lbuf,residence_time(llocpar),pface)
! 37a) Spherical JKR wall histories attached to this particle owner.
         IF(USE_JKR) THEN
            DO CC=1,COLLISION_ARRAY_MAX
               CALL UNPACK_DBUF(lbuf,WALL_COLLISION_FACET_ID(CC,llocpar),pface)
               CALL UNPACK_DBUF(lbuf,WALL_COLLISION_PFT(:,CC,llocpar),pface)
               CALL UNPACK_DBUF(lbuf,WALL_COLLISION_PFR(:,CC,llocpar),pface)
               CALL UNPACK_DBUF(lbuf,JKR_WALL_CONTACT_STATE(CC,llocpar),pface)
            ENDDO
         ENDIF
! 38) SuperDEM
         IF(SuperDEM) THEN
            call unpack_dbuf(lbuf,super_r(llocpar,1:3),pface)
            call unpack_dbuf(lbuf,super_mn(llocpar,1:2),pface)
            call unpack_dbuf(lbuf,super_q(llocpar,1:4),pface)
         ENDIF
! 38a) GluedSphereDEM
            IF(GSP_EXPLICIT) THEN
               call unpack_dbuf(lbuf, gp_sa(llocpar),pface)
               call unpack_dbuf(lbuf, gid_list(llocpar),pface)
               call unpack_dbuf(lbuf, gp_neighsid(llocpar,:),pface)
               call unpack_dbuf(lbuf, sc2gpc_vec(llocpar,:),pface)
               if(flexible_gsp) then
                  call unpack_dbuf(lbuf, bfn(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bmn(llocpar,:),pface)

                  call unpack_dbuf(lbuf, bftx(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bfty(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bftz(llocpar,:),pface)

                  call unpack_dbuf(lbuf, bmtx(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bmty(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bmtz(llocpar,:),pface)

                  call unpack_dbuf(lbuf, bnx(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bny(llocpar,:),pface)
                  call unpack_dbuf(lbuf, bnz(llocpar,:),pface)
               endif
               IF(ENERGY_EQ) THEN
                  call unpack_dbuf(lbuf, gp_kps(llocpar),pface)
                  call unpack_dbuf(lbuf, gp_neighsa(llocpar,:),pface)
               ENDIF
            ELSEIF(GSP_IMPLICIT) THEN
               call unpack_dbuf(lbuf, glueSurface(llocpar),pface)
               call unpack_dbuf(lbuf, glueQuat(llocpar,:),pface)
               call unpack_dbuf(lbuf, glueDiameter(llocpar),pface)
               call unpack_dbuf(lbuf, glueBounding(llocpar),pface)
            ENDIF
      end do

! 39) Number of pair datasets
      lbuf = lparcnt*iParticlePacketSize + ibufoffset
      call unpack_dbuf(lbuf,num_neighborlists_sent,pface)

      do nn = 1, num_neighborlists_sent
! 40) Global ID of packed particle.
         call unpack_dbuf(lbuf,lparid,pface)
! 41) DES grid IJK of cell receiving the particle.
         call unpack_dbuf(lbuf,lparijk,pface)

! Locate the particle on the current process.
         if (.not. locate_par(lparid,lparijk,llocpar)) then
             if (.not. exten_locate_par(lparid,lparijk,llocpar)) then
                 print *,"at buffer location",lbuf," pface = ",pface
                 print *,"Could not find particle ",lparid," IN IJK ",lparijk
                 call des_mpi_stop
             endif
         endif

! 42) Global ID of neighbor particle.
         call unpack_dbuf(lbuf,lneighid,pface)
! 43) DES grid IJK of cell containing the neighbor particle.
         call unpack_dbuf(lbuf,lneighijk,pface)

! Locate the neighbor particle on the current process.
         if (.not. locate_par(lneighid,lneighijk,lneigh)) then
            if (.not. exten_locate_par(lneighid,lparijk,lneigh)) then
               print *,"  "
               print *,"  "
               print *," fail on  ", myPE
               print *,"at buffer location",lbuf," pface = ",pface
               print *,"Could not find neighbor ",lneighid," IN IJK ",lneighijk
               call des_mpi_stop
            endif
         endif

! Queue pair history by global IDs. It is restored only after the next
! official neighbor rebuild has created valid CC slots on this owner.
! 44) Spherical JKR pair histories.
         IF(USE_JKR) THEN
            CALL UNPACK_DBUF(lbuf,PAIR_PFT_IN,pface)
            CALL UNPACK_DBUF(lbuf,PAIR_PFR_IN,pface)
            CALL UNPACK_DBUF(lbuf,PAIR_STATE_IN,pface)
            CALL QUEUE_PAIR_HISTORY(lparid,lneighid,PAIR_PFT_IN, &
               PAIR_PFR_IN,PAIR_STATE_IN)
         ENDIF
! 45) Superquadric contact points
         IF(SuperDEM) THEN
            call unpack_dbuf(lbuf,contact_point_a(:,cc),pface)
            call unpack_dbuf(lbuf,contact_point_b(:,cc),pface)
            call unpack_dbuf(lbuf,contact_lambda_a(cc),pface)
            call unpack_dbuf(lbuf,contact_lambda_b(cc),pface)
         ENDIF
      enddo

      END SUBROUTINE desmpi_unpack_parcross

!----------------------------------------------------------------------!
! Function: LOCATE_PAR                                                 !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: Return the local index of the particle matching the passed  !
!    global ID. The function returns TRUE if the particle is matched,  !
!    otherwise it returns FALSE.                                       !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION LOCATE_PAR(pGLOBALID, pIJK, pLOCALNO)

      use discretelement, only: iGLOBAL_ID
      use desgrid, only: DG_IJKStart2, DG_IJKEnd2
      use derived_types, only: dg_pic

      implicit none

! Dummy arguments:
!---------------------------------------------------------------------//
! Global ID of the particle
      INTEGER, INTENT(IN) :: pGlobalID
! IJK of DES grid cell containing the particle
      INTEGER, INTENT(IN) :: pIJK
! Local ID for the particle.
      INTEGER, INTENT(OUT) :: pLocalNO

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: lpicloc, lcurpar

! Initialize the result.
      locate_par = .false.

! Verify that the passied IJK value is within a valid range.
      if(pIJK < dg_ijkstart2 .or. pIJK > dg_ijkend2)  RETURN

! Loop the the particles in DES grid cell pIJK. Return to the calling
! routine if the passed global ID matches the global ID of one of
! the local particles.
      DO lpicloc = 1,dg_pic(pijk)%isize
         lcurpar = dg_pic(pijk)%p(lpicloc)
         IF(iGLOBAL_ID(lcurpar) == pGlobalID) THEN
            plocalno = lcurpar
            locate_par = .true.
            RETURN
         ENDIF
      ENDDO

      RETURN
      end function locate_par

!----------------------------------------------------------------------!
! Function: EXTEN_LOCATE_PAR                                           !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: Return the local index of the particle matching the passed  !
!    global ID. The function returns TRUE if the particle is matched,  !
!    otherwise it returns FALSE.                                       !
!------------------------------------------------------------------------
      LOGICAL FUNCTION EXTEN_LOCATE_PAR(pGlobalID, pIJK, pLocalNO)

      use derived_types, only: dg_pic
      use discretelement, only: iGLOBAL_ID
      use desgrid, only: DG_IJKStart2, DG_IJKEnd2
      use desgrid, only: dg_Iof_LO, DG_Jof_LO, DG_Kof_LO
      use geometry, only: NO_K

      use desgrid, only: dg_funijk

      implicit none

! Dummy variables:
!---------------------------------------------------------------------//
! The global ID of the particle to be matched locally
      INTEGER, INTENT(IN) :: pGlobalId
! The DES grid cell index expected to contain the particle.
      INTEGER, INTENT(IN) :: pIJK
! The local ID of the matching particle.
      INTEGER, INTENT(OUT) :: pLocalNo

! Local variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER :: lijk, li, lj, lk, lic, ljc, lkc, lkoffset
      INTEGER :: lpicloc,lcurpar

      exten_locate_par = .false.

      lic = dg_iof_lo(pijk)
      ljc = dg_jof_lo(pijk)
      lkc = dg_kof_lo(pijk)
      lkoffset = merge(0, 1, NO_K)
      DO  lk = lkc-lkoffset,lkc+lkoffset
      DO  lj = ljc-1,ljc+1
      DO  li = lic-1,lic+1
         lijk = dg_funijk(li,lj,lk)
         IF (lijk .lt. dg_ijkstart2 .or. lijk .gt. dg_ijkend2) CYCLE
         DO lpicloc = 1, dg_pic(lijk)%isize
            lcurpar = dg_pic(lijk)%p(lpicloc)
            IF (iglobal_id(lcurpar) .eq. pglobalid) THEN
               plocalno = lcurpar
               exten_locate_par = .true.
               RETURN
            END IF
         END DO
      END DO
      END DO
      END DO

      RETURN
      END FUNCTION EXTEN_LOCATE_PAR

!----------------------------------------------------------------------!
! Unpack subroutine for single real variables                          !
!----------------------------------------------------------------------!
      subroutine unpack_db0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata

      idata = drecvbuf(1+mod(pface,2))%facebuf(lbuf)
      lbuf = lbuf + 1

      return
      end subroutine unpack_db0

!----------------------------------------------------------------------!
! Unpack subroutine for real arrays                                    !
!----------------------------------------------------------------------!
      subroutine unpack_db1(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      idata = drecvbuf(1+mod(pface,2))%facebuf(lbuf:lbuf+lsize-1)
      lbuf = lbuf + lsize

      return
      end subroutine unpack_db1

!----------------------------------------------------------------------!
! Unpack subroutine for single integer variables                       !
!----------------------------------------------------------------------!
      subroutine unpack_i0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(inout) :: idata

      idata = drecvbuf(1+mod(pface,2))%facebuf(lbuf)
      lbuf = lbuf + 1

      return
      end subroutine unpack_i0

!----------------------------------------------------------------------!
! Unpack subroutine for integer arrays                                 !
!----------------------------------------------------------------------!
      subroutine unpack_i1(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      integer, intent(inout) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      idata = drecvbuf(1+mod(pface,2))%facebuf(lbuf:lbuf+lsize-1)
      lbuf = lbuf + lsize

      return
      end subroutine unpack_i1

!----------------------------------------------------------------------!
! Unpack subroutine for logical variables                              !
!----------------------------------------------------------------------!
      subroutine unpack_l0(lbuf,idata,pface)
      use desmpi, only: dRECVBUF
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      logical, intent(inout) :: idata

      idata = merge(.true.,.false.,0.5<drecvbuf(1+mod(pface,2))%facebuf(lbuf))
      lbuf = lbuf + 1

      return
      end subroutine unpack_l0



! Queue directed pair history until NEIGHBOUR has rebuilt valid CC slots.
      SUBROUTINE QUEUE_PAIR_HISTORY(OWNER_GID,NEIGH_GID,PFT,PFR,STATE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: OWNER_GID,NEIGH_GID
      DOUBLE PRECISION, INTENT(IN) :: PFT(3),PFR(3)
      LOGICAL, INTENT(IN) :: STATE

      IF(PENDING_PAIR_COUNT >= PENDING_PAIR_CAPACITY) &
         CALL GROW_PENDING_PAIR_HISTORY
      PENDING_PAIR_COUNT = PENDING_PAIR_COUNT + 1
      PENDING_PAIR_OWNER_GID(PENDING_PAIR_COUNT) = OWNER_GID
      PENDING_PAIR_NEIGH_GID(PENDING_PAIR_COUNT) = NEIGH_GID
      PENDING_PAIR_PFT(:,PENDING_PAIR_COUNT) = PFT
      PENDING_PAIR_PFR(:,PENDING_PAIR_COUNT) = PFR
      PENDING_PAIR_STATE(PENDING_PAIR_COUNT) = STATE
      CALL JKR_DIAG_PAIR_EVENT('unpack_queue',OWNER_GID,NEIGH_GID, &
         -1,-1,0D0,0D0,(/0D0,0D0,0D0/),PFT,PFR,STATE)
      END SUBROUTINE QUEUE_PAIR_HISTORY


      SUBROUTINE GROW_PENDING_PAIR_HISTORY
      IMPLICIT NONE
      INTEGER :: NEW_CAPACITY
      INTEGER, ALLOCATABLE :: TMP_I(:)
      DOUBLE PRECISION, ALLOCATABLE :: TMP_D(:,:)
      LOGICAL, ALLOCATABLE :: TMP_L(:)

      NEW_CAPACITY = MAX(16,2*PENDING_PAIR_CAPACITY)

      ALLOCATE(TMP_I(NEW_CAPACITY))
      IF(PENDING_PAIR_COUNT > 0) TMP_I(1:PENDING_PAIR_COUNT) = &
         PENDING_PAIR_OWNER_GID(1:PENDING_PAIR_COUNT)
      CALL MOVE_ALLOC(TMP_I,PENDING_PAIR_OWNER_GID)

      ALLOCATE(TMP_I(NEW_CAPACITY))
      IF(PENDING_PAIR_COUNT > 0) TMP_I(1:PENDING_PAIR_COUNT) = &
         PENDING_PAIR_NEIGH_GID(1:PENDING_PAIR_COUNT)
      CALL MOVE_ALLOC(TMP_I,PENDING_PAIR_NEIGH_GID)

      ALLOCATE(TMP_D(3,NEW_CAPACITY))
      IF(PENDING_PAIR_COUNT > 0) TMP_D(:,1:PENDING_PAIR_COUNT) = &
         PENDING_PAIR_PFT(:,1:PENDING_PAIR_COUNT)
      CALL MOVE_ALLOC(TMP_D,PENDING_PAIR_PFT)

      ALLOCATE(TMP_D(3,NEW_CAPACITY))
      IF(PENDING_PAIR_COUNT > 0) TMP_D(:,1:PENDING_PAIR_COUNT) = &
         PENDING_PAIR_PFR(:,1:PENDING_PAIR_COUNT)
      CALL MOVE_ALLOC(TMP_D,PENDING_PAIR_PFR)

      ALLOCATE(TMP_L(NEW_CAPACITY))
      IF(PENDING_PAIR_COUNT > 0) TMP_L(1:PENDING_PAIR_COUNT) = &
         PENDING_PAIR_STATE(1:PENDING_PAIR_COUNT)
      CALL MOVE_ALLOC(TMP_L,PENDING_PAIR_STATE)

      PENDING_PAIR_CAPACITY = NEW_CAPACITY
      END SUBROUTINE GROW_PENDING_PAIR_HISTORY


! Restore queued histories onto the official neighbor layout by global ID.
      SUBROUTINE RESTORE_PENDING_PAIR_HISTORY
      USE DISCRETELEMENT
      IMPLICIT NONE
      INTEGER :: N,LL,OWNER,NEIGH,CC,CC_START,CC_END
      LOGICAL :: FOUND

      DO N=1,PENDING_PAIR_COUNT
         OWNER = 0
         NEIGH = 0
         DO LL=1,MAX_PIP
            IF(IGLOBAL_ID(LL) == PENDING_PAIR_OWNER_GID(N)) OWNER = LL
            IF(IGLOBAL_ID(LL) == PENDING_PAIR_NEIGH_GID(N)) NEIGH = LL
         ENDDO
         IF(OWNER == 0 .OR. NEIGH == 0) CYCLE

         FOUND = .FALSE.
         CC_START = 1
         IF(OWNER > 1) CC_START = NEIGHBOR_INDEX(OWNER-1)
         CC_END = NEIGHBOR_INDEX(OWNER)
         DO CC=CC_START,CC_END-1
            IF(NEIGHBORS(CC) == NEIGH) THEN
               PFT_NEIGHBOR(:,CC) = PENDING_PAIR_PFT(:,N)
               PFR_NEIGHBOR(:,CC) = PENDING_PAIR_PFR(:,N)
               JKR_CONTACT_STATE(CC) = PENDING_PAIR_STATE(N)
               CALL JKR_DIAG_PAIR_EVENT('pending_restore_direct', &
                  PENDING_PAIR_OWNER_GID(N),PENDING_PAIR_NEIGH_GID(N), &
                  -1,CC,0D0,0D0,(/0D0,0D0,0D0/), &
                  PFT_NEIGHBOR(:,CC),PFR_NEIGHBOR(:,CC), &
                  JKR_CONTACT_STATE(CC))
               FOUND = .TRUE.
               EXIT
            ENDIF
         ENDDO

         IF(.NOT.FOUND) THEN
            CC_START = 1
            IF(NEIGH > 1) CC_START = NEIGHBOR_INDEX(NEIGH-1)
            CC_END = NEIGHBOR_INDEX(NEIGH)
            DO CC=CC_START,CC_END-1
               IF(NEIGHBORS(CC) == OWNER) THEN
                  PFT_NEIGHBOR(:,CC) = -PENDING_PAIR_PFT(:,N)
                  PFR_NEIGHBOR(:,CC) = -PENDING_PAIR_PFR(:,N)
                  JKR_CONTACT_STATE(CC) = PENDING_PAIR_STATE(N)
                  CALL JKR_DIAG_PAIR_EVENT('pending_restore_reverse', &
                     PENDING_PAIR_OWNER_GID(N),PENDING_PAIR_NEIGH_GID(N), &
                     -1,CC,0D0,0D0,(/0D0,0D0,0D0/), &
                     PFT_NEIGHBOR(:,CC),PFR_NEIGHBOR(:,CC), &
                     JKR_CONTACT_STATE(CC))
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      PENDING_PAIR_COUNT = 0
      END SUBROUTINE RESTORE_PENDING_PAIR_HISTORY

! Default-off, rank-local double-precision JKR diagnostic stream.
      SUBROUTINE JKR_DIAG_PAIR_EVENT(PHASE_TAG,OWNER_GID,NEIGH_GID, &
         WALL_ID,LOCAL_SLOT,OVERLAP,NORMAL_FORCE,TANGENTIAL_FORCE, &
         PFT,PFR,STATE)
      USE COMPAR, ONLY: MYPE
      USE RUN, ONLY: TIME
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: PHASE_TAG
      INTEGER, INTENT(IN) :: OWNER_GID,NEIGH_GID,WALL_ID,LOCAL_SLOT
      DOUBLE PRECISION, INTENT(IN) :: OVERLAP,NORMAL_FORCE
      DOUBLE PRECISION, INTENT(IN) :: TANGENTIAL_FORCE(3),PFT(3),PFR(3)
      LOGICAL, INTENT(IN) :: STATE
      CHARACTER(LEN=64) :: FILE_NAME
      CHARACTER(LEN=32) :: ENV_VALUE
      INTEGER :: ENV_LENGTH,ENV_STATUS

      IF(.NOT.JKR_DIAG_CHECKED) THEN
         ENV_VALUE = ''
         CALL GET_ENVIRONMENT_VARIABLE('MFIX_JKR_DIAGNOSTICS', &
            ENV_VALUE,ENV_LENGTH,ENV_STATUS)
         JKR_DIAG_ENABLED = ENV_STATUS == 0 .AND. ENV_LENGTH > 0 .AND. &
            TRIM(ENV_VALUE) /= '0'
         JKR_DIAG_CHECKED = .TRUE.
      ENDIF
      IF(.NOT.JKR_DIAG_ENABLED) RETURN

      IF(JKR_DIAG_UNIT == -1) THEN
         WRITE(FILE_NAME,'("JKR_DIAGNOSTIC_rank",I4.4,".dat")') MYPE
         OPEN(NEWUNIT=JKR_DIAG_UNIT,FILE=TRIM(FILE_NAME), &
            STATUS='REPLACE',ACTION='WRITE')
         WRITE(JKR_DIAG_UNIT,'(A)') &
            'time phase rank owner_gid neigh_gid wall_id local_slot ' // &
            'overlap normal_force tangential_force pft pfr jkr_state'
      ENDIF
      WRITE(JKR_DIAG_UNIT,*) TIME,TRIM(PHASE_TAG),MYPE,OWNER_GID, &
         NEIGH_GID,WALL_ID,LOCAL_SLOT,OVERLAP,NORMAL_FORCE, &
         TANGENTIAL_FORCE,PFT,PFR,STATE
      FLUSH(JKR_DIAG_UNIT)
      END SUBROUTINE JKR_DIAG_PAIR_EVENT

      end module mpi_unpack_des
