!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DES_ALLOCATE                                           C
!                                                                      C
!  Purpose: subroutines to allocate all DEM arrays                     C
!                                                                      C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE DES_ALLOCATE

   USE des_init_arrays_mod, only: des_init_arrays, des_init_particle_arrays

   PUBLIC:: DES_ALLOCATE_ARRAYS, ADD_PAIR, PARTICLE_GROW, ALLOCATE_DEM_MI, PARTICLE_GROW_GluedSphereDEM

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_ALLOCATE_ARRAYS                                     C
!  Purpose: Original allocate arrays subroutines for DES               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_ALLOCATE_ARRAYS

      USE compar
      USE constant
      USE cutcell
      USE des_psd
      USE derived_types, only: pic
      USE des_bc
      USE des_rxns
      USE des_thermo
      USE des_thermo_cond, only: DES_Qw_cond
      USE discretelement
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mfix_pic, only: mppic, des_stat_wt, ps_grad, epg_p
      USE mfix_pic, only: pic_bcmi
      USE param
      USE param1
      USE particle_filter, only: DES_INTERP_DPVM
      USE particle_filter, only: DES_INTERP_GARG
      USE particle_filter, only: DES_INTERP_GAUSS
      USE particle_filter, only: DES_INTERP_LHAT
      USE particle_filter, only: DES_INTERP_SATELLITE

      USE particle_filter, only: DES_INTERP_SCHEME_ENUM
      USE particle_filter, only: FILTER_CELL
      USE particle_filter, only: FILTER_SIZE
      USE particle_filter, only: FILTER_WEIGHT
      USE physprop
      USE run, only: ENERGY_EQ, ANY_SOLIDS_SPECIES_EQ

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      USE error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
      INTEGER :: GSP_ALLOC_SIZE
!-----------------------------------------------
! indices
      INTEGER :: IJK
!-----------------------------------------------

! For parallel processing the array size required should be either
! specified by the user or could be determined from total particles
! with some factor.

!It would be good if we could call read_part_input before this point
! so that PARTICLES would get set correctly.  Maybe we need a routine
! to get the number of particles from the file without fully loading it.

      MAX_PIP = merge(0, PARTICLES/numPEs, PARTICLES==UNDEFINED_I)
      MAX_PIP = MAX(MAX_PIP,4)

      WRITE(ERR_MSG,1000) trim(iVal(MAX_PIP))
      CALL LOG_MESSAGE(__FILE__, __LINE__, LOGLEVEL_INFO, HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT('Initial DES particle array size: ',A)

! DES Allocatable arrays
!-----------------------------------------------
! Dynamic particle info including another index for parallel
! processing for ghost
      ALLOCATE( PARTICLE_STATE (MAX_PIP) )
      ALLOCATE (iglobal_id(max_pip))

! R.Garg: Allocate necessary arrays for PIC mass inlet/outlet BCs
      IF(PIC_BCMI /= 0) CALL ALLOCATE_PIC_MI

! Particle attributes
! Radius, density, mass, moment of inertia
      Allocate(  DES_RADIUS (MAX_PIP) )
      Allocate(  RO_Sol (MAX_PIP) )
      Allocate(  PVOL (MAX_PIP) )
      Allocate(  PMASS (MAX_PIP) )
      Allocate(  OMOI (MAX_PIP) )
!SuperDEM
      IF(SuperDEM)  THEN
         Allocate(  OMOI3 (MAX_PIP,3) )
         Allocate(  super_r (MAX_PIP,3) )
         Allocate(  super_mn (MAX_PIP,2) )
         Allocate(  super_q (MAX_PIP,4) )
      ENDIF
! Glued Sphere DEM
      IF(GSP_EXPLICIT) THEN
! allocate size = number of glued spheres
! for flexible gsp, still need to calculate the number of glued particles
         IF(.NOT. GENER_PART_CONFIG) THEN
             IF(NGluedParticles == UNDEFINED_I) THEN
                 IF (particles == UNDEFINED_I) THEN
                     gsp_alloc_size = 1 ! will grow this in read_part_input
                 ELSE
                     gsp_alloc_size = MAX(particles,1)
                 ENDIF
             ENDIF
         ELSE
            NGluedParticles = 1
            ! set gsp_alloc_size to avoid memory leak
            gsp_alloc_size = 1
         ENDIF
         Allocate( particle_state_gsp(Gsp_Alloc_Size) ); particle_state_gsp(:) = 0
         Allocate( glueMass(Gsp_Alloc_Size) ); glueMass(:) = 0.0D0
         Allocate( glueVolume(Gsp_Alloc_Size) ); glueVolume(:) = 0.0D0
         Allocate( glueDiameter(Gsp_Alloc_Size) ); glueDiameter(:) = 0.0D0
         Allocate( glueBounding(Gsp_Alloc_Size) ); glueBounding(:) = 0.0D0
         Allocate( glueForce(Gsp_Alloc_Size, 3) ); glueForce(:,:) = 0.0D0
         Allocate( glueTorque(Gsp_Alloc_Size, 3) ); glueTorque(:,:) = 0.0D0

         Allocate( glueAcc(Gsp_Alloc_Size, 3) ); glueAcc(:,:) = 0.0D0
         Allocate( glueAngMom(Gsp_Alloc_Size, 3) ); glueAngMom(:,:) = 0.0D0
         Allocate( glueQuat(Gsp_Alloc_Size, 4) ); glueQuat(:,:) = 0.0D0
         Allocate( glueInertia(Gsp_Alloc_Size, 3) ); glueInertia(:,:) = 0.0D0
         Allocate( glueOmg(Gsp_Alloc_Size, 3) ); glueOmg(:,:) = 0.0D0
         Allocate( glueVel(Gsp_Alloc_Size, 3) ); glueVel(:,:) = 0.0D0
         Allocate( gluePos(Gsp_Alloc_Size, 3) ); gluePos(:,:) = 0.0D0

         Allocate( glueEX(Gsp_Alloc_Size, 3) ); glueEX(:,:) = ZERO
         Allocate( glueEY(Gsp_Alloc_Size, 3) ); glueEY(:,:) = ZERO
         Allocate( glueEZ(Gsp_Alloc_Size, 3) ); glueEZ(:,:) = ZERO

! allocate size = number of component spheres
         Allocate( gid_list (MAX_PIP) ); gid_list(:) = 0
         Allocate( gp_component_IC (MAX_PIP) ); gp_component_IC(:) = 0
         Allocate( cglobal_id (MAX_PIP) ); cglobal_id(:) = 0
         Allocate( gp_neighsid (MAX_PIP,6) ); gp_neighsid(:,:) = 0
         Allocate( gp_neighsa (MAX_PIP,6) ); gp_neighsa(:,:) = 0.0D0
         Allocate( gp_kps (MAX_PIP) ); gp_kps(:) = 0.0D0
         Allocate( gp_sa (MAX_PIP) ); gp_sa(:) = 0.0D0
         ! Allocate( gp_sdrag_mass (MAX_PIP) )
         Allocate( gp_squat (MAX_PIP,4) ); gp_squat(:,:) = 0.0D0
         Allocate( sc2gpc_vec (MAX_PIP,3) ); sc2gpc_vec(:,:) = 0.0D0

         if(flexible_gsp) then
            Allocate(bfn(MAX_PIP,6)); bfn = ZERO
            Allocate(bmn(MAX_PIP,6)); bmn = ZERO

            Allocate(bftx(MAX_PIP,6)); bftx = ZERO
            Allocate(bfty(MAX_PIP,6)); bfty = ZERO
            Allocate(bftz(MAX_PIP,6)); bftz = ZERO

            Allocate(bmtx(MAX_PIP,6)); bmtx = ZERO
            Allocate(bmty(MAX_PIP,6)); bmty = ZERO
            Allocate(bmtz(MAX_PIP,6)); bmtz = ZERO

            ! history of normal
            Allocate(bnx(MAX_PIP,6)); bnx = ZERO
            Allocate(bny(MAX_PIP,6)); bny = ZERO
            Allocate(bnz(MAX_PIP,6)); bnz = ZERO
         endif

         ! If using auto_seeing function, need to reset NGluedParticles = 0
         IF(GENER_PART_CONFIG) NGluedParticles = 0
      ELSEIF(GSP_IMPLICIT) THEN
         Allocate( OMOI3 (MAX_PIP,3) )
         Allocate( glueAngMom(MAX_PIP, 3) ); glueAngMom(:,:) = ZERO
         Allocate( glueQuat(MAX_PIP, 4) ); glueQuat(:,:) = ZERO
         Allocate( glueEX(MAX_PIP, 3) ); glueEX(:,:) = ZERO
         Allocate( glueEY(MAX_PIP, 3) ); glueEY(:,:) = ZERO
         Allocate( glueEZ(MAX_PIP, 3) ); glueEZ(:,:) = ZERO
         Allocate( glueSurface(MAX_PIP) ); glueSurface(:) = ZERO
         Allocate( glueBounding(MAX_PIP) ); glueBounding(:) = 0.0D0
         Allocate( glueDiameter(MAX_PIP) ); glueDiameter(:) = 0.0D0
         IF(GENER_PART_CONFIG) NGluedParticles = 0
         IF(GENER_PART_CONFIG) NSphereGSP = 0
      ENDIF
! Old and new particle positions, velocities (translational and
! rotational)
      Allocate(  DES_POS_NEW (MAX_PIP,DIMN) )
      Allocate(  DES_VEL_NEW (MAX_PIP,DIMN) )
      Allocate(  OMEGA_NEW (MAX_PIP,DIMN) )
      call create_custom_datatype

! Force chain arrays
! Arrays are initially allocated with 1000 elements.
! Arrays are grown when needed in calc_force_dem.f
      Allocate(FCHAIN_MIDPOINT (1000,DIMN) )
      FCHAIN_MIDPOINT(:,:) = UNDEFINED
      Allocate(FCHAIN_CONTACT_POINT (1000,DIMN) )
      FCHAIN_CONTACT_POINT(:,:) = UNDEFINED
      Allocate(FCHAIN_ORIENT (1000,DIMN) )
      Allocate(FCHAIN_FN (1000,DIMN) )
      Allocate(FCHAIN_FT (1000,DIMN) )
      Allocate(FCHAIN_FC (1000,DIMN) )
      Allocate(FCHAIN_LENGTH (1000))
      Allocate(FCHAIN_FN_MAG (1000))
      Allocate(FCHAIN_FCMAX (1000))
      Allocate(FCHAIN_LOCAL_ID1 (1000))
      Allocate(FCHAIN_LOCAL_ID2 (1000))
      Allocate(FCHAIN_GLOBAL_ID1 (1000))
      Allocate(FCHAIN_GLOBAL_ID2 (1000))
      Allocate(FCHAIN_OVERLAP (1000))
      Allocate(FCHAIN_JKR_STATE (1000))
      Allocate(FCHAIN_JKR_RADIUS (1000))
      Allocate(FCHAIN_JKR_GAMMA (1000))
      Allocate(FCHAIN_TANGENTIAL_DISPLACEMENT (1000,DIMN))
      Allocate(FCHAIN_ROLLING_DISPLACEMENT (1000,DIMN))
      Allocate(FCHAIN_TWIST_ANGLE (1000))
      FCHAIN_TWIST_ANGLE = ZERO

      Allocate(MIN_OVERLAP (MAX_PIP))
      Allocate(MAX_OVERLAP (MAX_PIP))
      Allocate(MEAN_OVERLAP (MAX_PIP))
      Allocate(COORDINATION (MAX_PIP))
      MIN_OVERLAP(:) = ZERO
      MAX_OVERLAP(:) = ZERO
      MEAN_OVERLAP(:) = ZERO
      COORDINATION(:) = ZERO

! Collision force
      Allocate(  DES_COLL_FORCE (MAX_PIP,DIMN) )

! Residence time
      Allocate(RESIDENCE_TIME (MAX_PIP))
      RESIDENCE_TIME(:) = ZERO

      IF(PARTICLE_ORIENTATION) Allocate(  ORIENTATION (MAX_PIP,DIMN) )

      IF (DO_OLD) THEN
         Allocate(  DES_POS_OLD (MAX_PIP,DIMN) )
         Allocate(  DES_VEL_OLD (MAX_PIP,DIMN) )
         Allocate(  DES_ACC_OLD (MAX_PIP,DIMN) )
         Allocate(  OMEGA_OLD (MAX_PIP,DIMN) )
         Allocate(  ROT_ACC_OLD (MAX_PIP,DIMN))
      ENDIF

! Allocating user defined array
      IF(DES_USR_VAR_SIZE > 0) &
         Allocate( DES_USR_VAR(DES_USR_VAR_SIZE,MAX_PIP) )

! Particle positions at the last call neighbor search algorithm call
      Allocate(  PPOS (MAX_PIP,DIMN) )

! Total, normal and tangential forces
      Allocate(  FC (MAX_PIP,DIMN) )

! Torque
      Allocate(  TOW (MAX_PIP,DIMN) )


! allocate variable for des grid binning
      allocate(dg_pijk(max_pip)); dg_pijk=0
      allocate(dg_pijkprv(max_pip)); dg_pijkprv=0

! allocate variables related to ghost particles
      allocate(ighost_updated(max_pip))



      Allocate(  wall_collision_facet_id (COLLISION_ARRAY_MAX, MAX_PIP) )
      wall_collision_facet_id(:,:) = -1
      Allocate(  wall_collision_PFT (DIMN, COLLISION_ARRAY_MAX, MAX_PIP) )
      Allocate(  wall_collision_PFTW(COLLISION_ARRAY_MAX, MAX_PIP) )
      wall_collision_PFTW = ZERO
      Allocate(  wall_collision_PFR (DIMN, COLLISION_ARRAY_MAX, MAX_PIP) )
      wall_collision_PFR = ZERO
      Allocate(  WALL_COLLISION_FRAME (4, COLLISION_ARRAY_MAX, MAX_PIP) )
      WALL_COLLISION_FRAME = ZERO
      Allocate(  JKR_WALL_CONTACT_STATE (COLLISION_ARRAY_MAX, MAX_PIP) )
      JKR_WALL_CONTACT_STATE = .FALSE.

! Temporary variables to store wall position, velocity and normal vector
      Allocate(  WALL_NORMAL  (NWALLS,DIMN) )

      Allocate(  NEIGHBOR_INDEX (MAX_PIP) )
      Allocate(  NEIGHBOR_INDEX_OLD (MAX_PIP) )
      Allocate(  NEIGHBORS (MAX_PIP) )
      NEIGHBORS(:) = 0

      Allocate(  NEIGHBORS_OLD (MAX_PIP) )
      Allocate(  PFT_NEIGHBOR (3,MAX_PIP) )
      Allocate(  PFT_NEIGHBOR_OLD (3,MAX_PIP) )
      Allocate(  PFR_NEIGHBOR (3,MAX_PIP) )
      PFR_NEIGHBOR = ZERO
      Allocate(  PFR_NEIGHBOR_OLD (3,MAX_PIP) )
      PFR_NEIGHBOR_OLD = ZERO
      Allocate(  JKR_FRAME_NEIGHBOR (4,MAX_PIP) )
      JKR_FRAME_NEIGHBOR = ZERO
      Allocate(  JKR_FRAME_NEIGHBOR_OLD (4,MAX_PIP) )
      JKR_FRAME_NEIGHBOR_OLD = ZERO
      Allocate(  PFTW_NEIGHBOR (MAX_PIP) )
      PFTW_NEIGHBOR = ZERO
      Allocate(  PFTW_NEIGHBOR_OLD (MAX_PIP) )
      PFTW_NEIGHBOR_OLD = ZERO
      Allocate(  JKR_CONTACT_STATE (MAX_PIP) )
      JKR_CONTACT_STATE = .FALSE.
      Allocate(  JKR_CONTACT_STATE_OLD (MAX_PIP) )
      JKR_CONTACT_STATE_OLD = .FALSE.

! superquadric model, to store contact points in global coordinates
      IF(SuperDEM)  THEN
         Allocate(  CONTACT_POINT_A (3,MAX_PIP) )
         Allocate(  CONTACT_POINT_A_OLD (3,MAX_PIP) )
         Allocate(  CONTACT_POINT_B (3,MAX_PIP) )
         Allocate(  CONTACT_POINT_B_OLD (3,MAX_PIP) )

         Allocate(  CONTACT_LAMBDA_A (MAX_PIP) )
         Allocate(  CONTACT_LAMBDA_A_OLD (MAX_PIP) )
         Allocate(  CONTACT_LAMBDA_B (MAX_PIP) )
         Allocate(  CONTACT_LAMBDA_B_OLD (MAX_PIP) )

      ENDIF
! Variable that stores the particle in cell information (ID) on the
! computational fluid grid defined by keywords IMAX, JMAX and KMAX.
      ALLOCATE(PIC(DIMENSION_3))
      DO IJK=1,DIMENSION_3
        NULLIFY(pic(ijk)%p)
      ENDDO

! Particles in a computational fluid cell (for volume fraction)
      Allocate(  PINC (DIMENSION_3) )

! Variables for masking and sorting particles
!Hari Sitaraman===============
      if(allocated(valid_fluid_indices)) deallocate(valid_fluid_indices)
      if(allocated(vol_surr_inv)) deallocate(vol_surr_inv)
      if(allocated(fluid_at_mask1)) deallocate(fluid_at_mask1)
      if(allocated(fluid_at_mask2)) deallocate(fluid_at_mask2)
      if(allocated(fluid_at_mask3)) deallocate(fluid_at_mask3)

      Allocate(valid_fluid_indices(DIMENSION_3))
      Allocate(vol_surr_inv(DIMENSION_3))
      Allocate(fluid_at_mask1(DIMENSION_3))
      Allocate(fluid_at_mask2(IEND2-ISTART2+1,JEND2-JSTART2+1,KEND2-KSTART2+1))
      Allocate(fluid_at_mask3(DIMENSION_3))
!============================

! Ghost Particles in a computational fluid cell (for volume fraction)
      Allocate(  GPINC (DIMENSION_3) )

! For each particle track its i,j,k location on computational fluid grid
! defined by keywords IMAX, JMAX and KMAX.
      Allocate(  PIJK (MAX_PIP,5) )

      ALLOCATE(DRAG_AM(DIMENSION_3))
      ALLOCATE(DRAG_BM(DIMENSION_3, DIMN))
      ALLOCATE(F_gp(MAX_PIP ))
      F_gp(1:MAX_PIP)  = ZERO

! Explicit drag force acting on a particle.
      Allocate(DRAG_FC (MAX_PIP,DIMN) )

! force due to gas-pressure gradient
      ALLOCATE(P_FORCE(DIMN, DIMENSION_3))

! Volume of nodes
      ALLOCATE(DES_VOL_NODE(DIMENSION_3))

      ALLOCATE(F_GDS(DIMENSION_3))
      ALLOCATE(VXF_GDS(DIMENSION_3))

      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_DPVM, DES_INTERP_GAUSS, DES_INTERP_LHAT, DES_INTERP_SATELLITE)
         ALLOCATE(FILTER_CELL(FILTER_SIZE, MAX_PIP))
         ALLOCATE(FILTER_WEIGHT(FILTER_SIZE, MAX_PIP))
      CASE(DES_INTERP_GARG)
         ALLOCATE(DES_ROPS_NODE(DIMENSION_3, DIMENSION_M))
         ALLOCATE(DES_VEL_NODE(DIMENSION_3, DIMN, DIMENSION_M))
      END SELECT

! Variables for hybrid model
      IF (DES_CONTINUUM_HYBRID) THEN
         ALLOCATE(SDRAG_AM(DIMENSION_3,DIMENSION_M))
         ALLOCATE(SDRAG_BM(DIMENSION_3, DIMN,DIMENSION_M))

         ALLOCATE(F_SDS(DIMENSION_3,DIMENSION_M))
         ALLOCATE(VXF_SDS(DIMENSION_3,DIMENSION_M))
      ENDIF

! coarse grained DEM
      IF (CGDEM) Allocate(  DES_CGP_STW (MAX_PIP) )
      IF (CGDEM) Allocate(  DES_CGP_RPR (MAX_PIP) )

! MP-PIC related
      IF(MPPIC) THEN
         ALLOCATE(DES_STAT_WT(MAX_PIP))
         ALLOCATE(DES_VEL_MAX(DIMN))
         ALLOCATE(PS_GRAD(3,MAX_PIP))
         ALLOCATE(EPG_P(MAX_PIP))
      ENDIF

! Averaged velocity obtained by averaging over all the particles
      ALLOCATE(DES_VEL_AVG(DIMN) )

! Global Granular Energy
      ALLOCATE(GLOBAL_GRAN_ENERGY(DIMN) )
      ALLOCATE(GLOBAL_GRAN_TEMP(DIMN) )

! variable for bed height of solids phase M
      ALLOCATE(BED_HEIGHT(DIMENSION_M))

! ---------------------------------------------------------------->>>
! BEGIN COHESION
      IF(USE_COHESION) THEN
! Matrix location of particle  (should be allocated in case user wishes
! to invoke routines in /cohesion subdirectory
         Allocate(  PostCohesive (MAX_PIP) )
      ENDIF
! END COHESION
! ----------------------------------------------------------------<<<

! ---------------------------------------------------------------->>>
! BEGIN Thermodynamic Allocation
      IF(ENERGY_EQ)THEN
! Particle temperature
         Allocate( DES_T_s( MAX_PIP ) )
! Species mass fractions comprising a particle. This array may not be
! needed for all thermo problems.
         Allocate( DES_X_s( MAX_PIP, DIMENSION_N_S))
! Specific heat
         Allocate( DES_C_PS( MAX_PIP ) )
! Total rate of heat transfer to individual particles.
         Allocate( Q_Source( MAX_PIP ) )
! Average solids temperature in fluid cell
         Allocate(avgDES_T_s(DIMENSION_3) )
! Gas/Solids convective heat transfer coupling
         IF(CALC_CONV_DES) THEN
! Fluid phase energy equation source terms
            Allocate(CONV_Sc(DIMENSION_3) )
            Allocate(CONV_Sp(DIMENSION_3) )
! Particle convection source term (explicit coupled)
            Allocate(CONV_Qs(MAX_PIP))
! Gas-particle heat transfer coefficient TIMES surface area
            Allocate(GAMMAxSA(MAX_PIP))
         ENDIF

! Allocate the history variables for Adams-Bashforth integration
         IF (INTG_ADAMS_BASHFORTH) &
            Allocate( Q_Source0( MAX_PIP ) )
! Allocate the array for storing particle-wall heat transfer per unit area
         IF (ANY(CALC_COND_DES)) &
            Allocate( DES_Qw_cond( DIMENSION_3, DIMENSION_M))
      ENDIF
! End Thermodynamic Allocation
! ----------------------------------------------------------------<<<


! ---------------------------------------------------------------->>>
! BEGIN Species Allocation
      IF(ANY_SOLIDS_SPECIES_EQ)THEN
! Species mass fractions comprising a particle. This array may not be
! needed for all thermo problems.
         if(.not.allocated(DES_X_s)) Allocate( DES_X_s( MAX_PIP, DIMENSION_N_S))
! Rate of solids phase production/consumption for each species
         Allocate( DES_R_s( MAX_PIP, DIMENSION_N_s) )
! Minimum particle mass of each phase.
         Allocate( DES_MIN_PMASS( DIMENSION_M))

         Allocate( DES_R_gp( DIMENSION_3, DIMENSION_N_g ) )
         Allocate( DES_R_gc( DIMENSION_3, DIMENSION_N_g ) )
         Allocate( DES_SUM_R_g( DIMENSION_3 ) )
         Allocate( DES_R_PHASE( DIMENSION_3, DIMENSION_LM+DIMENSION_M-1 ) )
         Allocate( DES_HOR_g( DIMENSION_3 ) )


! Allocate the history variables for Adams-Bashforth integration
         IF (INTG_ADAMS_BASHFORTH) THEN
! Rate of change of particle mass
            Allocate( dMdt_OLD( MAX_PIP ) )
! Rate of change of particle mass percent species
            Allocate( dXdt_OLD( MAX_PIP, DIMENSION_N_s) )
         ENDIF

! Energy generation from reaction (cal/sec)
         Allocate( RXNS_Qs( MAX_PIP ) )
      ENDIF
! End Species Allocation
! ----------------------------------------------------------------<<<


! Allocate DES reaction rates arrays used to write reaction rates in vtp files or monitors
      if (SAVE_PART_RRATES) then
         if(allocated(Part_RRates_out)) deallocate(Part_RRates_out)
         Allocate( Part_RRates_out( MAX_PIP, NO_OF_DES_RXNS))
         Part_RRates_out(:,:) = ZERO
      endif


      RETURN
      END SUBROUTINE DES_ALLOCATE_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_DEM_MIO                                        !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_DEM_MI

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: undefined
      USE des_bc, only: dem_bcmi
      USE des_bc, only: pi_factor, pi_count
      use des_bc, only: numfrac_limit
      use des_bc, only: dem_mi_time, dem_bc_poly_layout
      use des_bc, only: dem_mi
      use des_bc, only: dem_bcmi_ijkstart, dem_bcmi_ijkend
      IMPLICIT NONE
!-----------------------------------------------

! Particle injection factor
      Allocate( PI_FACTOR (DEM_BCMI) )
! Particle injection count (injection number)
      Allocate( PI_COUNT (DEM_BCMI) )
! Particle injection time scale
      Allocate( DEM_MI_TIME (DEM_BCMI) )
! Array used for polydisperse inlets: stores the particle number
! distribution of an inlet scaled with numfrac_limit
      Allocate( DEM_BC_POLY_LAYOUT( DEM_BCMI, NUMFRAC_LIMIT ) )
! Data structure for storing BC data.
      Allocate( DEM_MI(DEM_BCMI) )

! Initialization
! Integer arrays
      PI_FACTOR(:) = -1
      PI_COUNT(:) = -1
      DEM_BC_POLY_LAYOUT(:,:) = -1
! Double precision arrays
      DEM_MI_TIME(:) = UNDEFINED

      allocate( DEM_BCMI_IJKSTART(DEM_BCMI) )
      allocate( DEM_BCMI_IJKEND(DEM_BCMI) )

      DEM_BCMI_IJKSTART = -1
      DEM_BCMI_IJKEND   = -1

! Boundary classification
!         Allocate( PARTICLE_PLCMNT (DES_BCMI) )
! Character precision arrays
!         PARTICLE_PLCMNT(:) = UNDEFINED_C

      RETURN
      END SUBROUTINE ALLOCATE_DEM_MI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_PIC_MI                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: R. Garg                                    Date: 11-Jun-14  !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_PIC_MI

! Modules
!-----------------------------------------------
      USE mfix_pic, only: pic_bcmi

      use mfix_pic, only: pic_mi

      IMPLICIT NONE
!-----------------------------------------------

! Allocate/Initialize for inlets
      if(pic_bcmi /= 0) allocate( pic_mi(pic_bcmi))


      RETURN
      END SUBROUTINE ALLOCATE_PIC_MI



!``````````````````````````````````````````````````````````````````````!
! Subroutine: ADD_PAIR                                                 !
!                                                                      !
! Purpose: Adds a neighbor pair to the pairs array.                    !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION add_pair(ii,jj)
      USE discretelement
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ii,jj

      CALL NEIGHBOR_GROW(NEIGHBOR_INDEX(ii))

      NEIGHBORS(NEIGHBOR_INDEX(ii)) = jj
      NEIGHBOR_INDEX(ii) = NEIGHBOR_INDEX(ii) + 1
      add_pair = NEIGHBOR_INDEX(ii)

      RETURN
      END FUNCTION add_pair

!``````````````````````````````````````````````````````````````````````!
! Subroutine: NEIGHBOR_GROW                                            !
!                                                                      !
! Purpose: Grow neighbors arrays to new_neigh_max. Note that neighbor      !
! max should be increased before calling this routine. Also, no        !
! assumption to the previous array size is made as needed for restarts.!
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE NEIGHBOR_GROW(new_neigh_max)
        USE discretelement
        USE geometry
        IMPLICIT NONE

        integer, intent(in) :: new_neigh_max

        INTEGER :: lSIZE1
        INTEGER, DIMENSION(:), ALLOCATABLE :: neigh_tmp
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pf_tmp
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: frame_tmp
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pfw_tmp
        LOGICAL, DIMENSION(:), ALLOCATABLE :: jkr_tmp
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: contact_point_a_tmp
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: contact_point_b_tmp
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: contact_lambda_a_tmp
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: contact_lambda_b_tmp

        INTEGER new_size

        lSIZE1 = size(neighbors,1)

        IF ( new_neigh_max .le. lSIZE1 ) RETURN

        new_size = lSIZE1

        DO WHILE(new_size < new_neigh_max)
           new_size = 2*new_size
        ENDDO

        allocate(neigh_tmp(new_size))
        neigh_tmp(1:lSIZE1) = neighbors(1:lSIZE1)
        neigh_tmp(lSIZE1+1:) = 0
        call move_alloc(neigh_tmp,neighbors)

        allocate(neigh_tmp(new_size))
        neigh_tmp(1:lSIZE1) = neighbors_old(1:lSIZE1)
        neigh_tmp(lSIZE1+1:) = 0
        call move_alloc(neigh_tmp,neighbors_old)

        allocate(pf_tmp(3,new_size))
        pf_tmp(:,1:lSIZE1) = pft_neighbor(:,1:lSIZE1)
        pf_tmp(:,lSIZE1+1:) = 0
        call move_alloc(pf_tmp,pft_neighbor)

        allocate(pf_tmp(3,new_size))
        pf_tmp(:,1:lSIZE1) = pft_neighbor_old(:,1:lSIZE1)
        pf_tmp(:,lSIZE1+1:) = 0
        call move_alloc(pf_tmp,pft_neighbor_old)

        allocate(pf_tmp(3,new_size))
        pf_tmp(:,1:lSIZE1) = pfr_neighbor(:,1:lSIZE1)
        pf_tmp(:,lSIZE1+1:) = 0
        call move_alloc(pf_tmp,pfr_neighbor)

        allocate(pf_tmp(3,new_size))
        pf_tmp(:,1:lSIZE1) = pfr_neighbor_old(:,1:lSIZE1)
        pf_tmp(:,lSIZE1+1:) = 0
        call move_alloc(pf_tmp,pfr_neighbor_old)

        allocate(frame_tmp(4,new_size))
        frame_tmp(:,1:lSIZE1) = jkr_frame_neighbor(:,1:lSIZE1)
        frame_tmp(:,lSIZE1+1:) = 0
        call move_alloc(frame_tmp,jkr_frame_neighbor)

        allocate(frame_tmp(4,new_size))
        frame_tmp(:,1:lSIZE1) = jkr_frame_neighbor_old(:,1:lSIZE1)
        frame_tmp(:,lSIZE1+1:) = 0
        call move_alloc(frame_tmp,jkr_frame_neighbor_old)

        allocate(pfw_tmp(new_size))
        pfw_tmp(1:lSIZE1) = pftw_neighbor(1:lSIZE1)
        pfw_tmp(lSIZE1+1:) = 0
        call move_alloc(pfw_tmp,pftw_neighbor)

        allocate(pfw_tmp(new_size))
        pfw_tmp(1:lSIZE1) = pftw_neighbor_old(1:lSIZE1)
        pfw_tmp(lSIZE1+1:) = 0
        call move_alloc(pfw_tmp,pftw_neighbor_old)
        allocate(jkr_tmp(new_size))
        jkr_tmp(1:lSIZE1) = jkr_contact_state(1:lSIZE1)
        jkr_tmp(lSIZE1+1:) = .FALSE.
        call move_alloc(jkr_tmp,jkr_contact_state)

        allocate(jkr_tmp(new_size))
        jkr_tmp(1:lSIZE1) = jkr_contact_state_old(1:lSIZE1)
        jkr_tmp(lSIZE1+1:) = .FALSE.
        call move_alloc(jkr_tmp,jkr_contact_state_old)
! SupderDEM
        if (SuperDEM) then
           allocate(contact_point_A_tmp(3,new_size))
           contact_point_A_tmp(:,1:lSIZE1) = contact_point_A(:,1:lSIZE1)
           contact_point_A_tmp(:,lSIZE1+1:) = 0
           call move_alloc(contact_point_A_tmp,contact_point_A)

           allocate(contact_point_A_tmp(3,new_size))
           contact_point_A_tmp(:,1:lSIZE1) = contact_point_A_old(:,1:lSIZE1)
           contact_point_A_tmp(:,lSIZE1+1:) = 0
           call move_alloc(contact_point_A_tmp,contact_point_A_old)

           allocate(contact_point_B_tmp(3,new_size))
           contact_point_B_tmp(:,1:lSIZE1) = contact_point_B(:,1:lSIZE1)
           contact_point_B_tmp(:,lSIZE1+1:) = 0
           call move_alloc(contact_point_B_tmp,contact_point_B)

           allocate(contact_point_B_tmp(3,new_size))
           contact_point_B_tmp(:,1:lSIZE1) = contact_point_B_old(:,1:lSIZE1)
           contact_point_B_tmp(:,lSIZE1+1:) = 0
           call move_alloc(contact_point_B_tmp,contact_point_B_old)


           allocate(contact_lambda_a_tmp(new_size))
           contact_lambda_a_tmp(1:lSIZE1) = contact_lambda_a(1:lSIZE1)
           contact_lambda_a_tmp(lSIZE1+1:) = 0
           call move_alloc(contact_lambda_A_tmp,contact_lambda_A)

           allocate(contact_lambda_a_tmp(new_size))
           contact_lambda_a_tmp(1:lSIZE1) = contact_lambda_a_old(1:lSIZE1)
           contact_lambda_a_tmp(lSIZE1+1:) = 0
           call move_alloc(contact_lambda_A_tmp,contact_lambda_A_old)

           allocate(contact_lambda_b_tmp(new_size))
           contact_lambda_b_tmp(1:lSIZE1) = contact_lambda_b(1:lSIZE1)
           contact_lambda_b_tmp(lSIZE1+1:) = 0
           call move_alloc(contact_lambda_B_tmp,contact_lambda_B)

           allocate(contact_lambda_b_tmp(new_size))
           contact_lambda_b_tmp(1:lSIZE1) = contact_lambda_b_old(1:lSIZE1)
           contact_lambda_b_tmp(lSIZE1+1:) = 0
           call move_alloc(contact_lambda_B_tmp,contact_lambda_B_old)

        endif

      END SUBROUTINE NEIGHBOR_GROW

      SUBROUTINE clear_block(A, lo, hi)
         USE param1, only: undefined, zero
         IMPLICIT NONE
         DOUBLE PRECISION, INTENT(INOUT) :: A(:,:)
         INTEGER , INTENT(IN)    :: lo, hi

         A(:, lo:hi) = zero
         A(lo:hi, :) = zero
      END SUBROUTINE clear_block
!``````````````````````````````````````````````````````````````````````!
! Subroutine: PARTICLE_GROW                                            !
!                                                                      !
! Purpose: Grow particle arrays to new_max_pip. Note that pair         !
! max should be increased before calling this routine. Also, no        !
! assumption to the previous array size is made as needed for restarts.!
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE PARTICLE_GROW(new_max_pip)

        USE des_rxns
        USE des_thermo
        USE discretelement
        USE mfix_pic
        USE param1, only: undefined, zero
        USE particle_filter
        USE resize
        USE run
        USE VTK, ONLY: SAVE_PART_RRATES, PART_RRATES_OUT
        use mpi_utility

        IMPLICIT NONE

        integer, intent(in) :: new_max_pip
        integer :: old_size, new_size
        integer :: redundant
        integer :: lsz1, lsz2
        logical, allocatable :: jkrw_tmp(:,:)

        max_pip = max(max_pip, new_max_pip)
        IF (new_max_pip .le. size(des_radius)) RETURN

        old_size = size(des_radius)

        new_size = old_size

        DO WHILE (new_size < new_max_pip)
           new_size = 2*new_size
        ENDDO

        call real_grow(des_radius,new_size)
        call real_grow(RO_Sol,new_size)
        call real_grow(PVOL,new_size)
        call real_grow(PMASS,new_size)
        call real_grow(OMOI,new_size)
! SuperDEM
        IF(SuperDEM)  THEN
            call real_grow2_reverse(OMOI3,new_size)
            call real_grow2_reverse(super_r,new_size)
            call real_grow2_reverse(super_mn,new_size)
            call real_grow2_reverse(super_q,new_size)
        ENDIF
! GluedSphereDEM
! Case 1: grow due to auto-seeding
! Case 2: grow due to glued sphere mass inflow
      IF(GSP_EXPLICIT) THEN
      ! in mpi, those are locally arrays and each pe has its own copies
      ! now gid_list has local copies as well
         redundant = old_size + 1

         call integer_grow(gid_list, new_size)
         call real_grow(gp_sa, new_size)
         call real_grow(gp_kps, new_size)
         ! call real_grow(gp_sdrag_mass, new_size)
         call integer_grow2_reverse(gp_neighsid, new_size)
         call real_grow2_reverse(gp_neighsa, new_size)
         call real_grow2_reverse(gp_squat, new_size)
         call real_grow2_reverse(sc2gpc_vec, new_size)
         if(flexible_gsp) then
            call real_grow2_reverse(bfn, new_size)
            call real_grow2_reverse(bmn, new_size)

            call real_grow2_reverse(bftx, new_size)
            call real_grow2_reverse(bfty, new_size)
            call real_grow2_reverse(bftz, new_size)

            call real_grow2_reverse(bmtx, new_size)
            call real_grow2_reverse(bmty, new_size)
            call real_grow2_reverse(bmtz, new_size)

            call real_grow2_reverse(bnx, new_size)
            call real_grow2_reverse(bny, new_size)
            call real_grow2_reverse(bnz, new_size)
         endif

         ! begin reset
         gid_list(redundant:new_size) = ZERO
         gp_sa(redundant:new_size) = ZERO
         gp_kps(redundant:new_size) = ZERO
         gp_neighsid(redundant:new_size,:) = ZERO
         gp_neighsa(redundant:new_size,:) = ZERO
         gp_squat(redundant:new_size,:) = ZERO
         sc2gpc_vec(redundant:new_size,:) = ZERO
         if(flexible_gsp) then
            bfn(redundant:new_size,:) = ZERO
            bmn(redundant:new_size,:) = ZERO

            bftx(redundant:new_size,:) = ZERO
            bfty(redundant:new_size,:) = ZERO
            bftz(redundant:new_size,:) = ZERO

            bmtx(redundant:new_size,:) = ZERO
            bmty(redundant:new_size,:) = ZERO
            bmtz(redundant:new_size,:) = ZERO

            bnx(redundant:new_size,:) = ZERO
            bny(redundant:new_size,:) = ZERO
            bnz(redundant:new_size,:) = ZERO
         endif
         ! special case for cglobal_id, after initialization, cglobal was already deallocated (see make_arrays_des.f)
         ! so after initialization, cglobal_id should not be changed anymore
         if(allocated(cglobal_id)) then
            call integer_grow(cglobal_id, new_size)
            cglobal_id(redundant:new_size) = ZERO
         endif

         if(allocated(gp_component_IC)) then
            call integer_grow(gp_component_IC, new_size)
            gp_component_IC(redundant:new_size) = ZERO
         endif
      ELSEIF(GSP_IMPLICIT) THEN
         redundant = old_size + 1

         call real_grow2_reverse(OMOI3, new_size)
         call real_grow2_reverse(glueAngMom, new_size)
         call real_grow2_reverse(glueQuat, new_size)
         call real_grow2_reverse(glueEX, new_size)
         call real_grow2_reverse(glueEY, new_size)
         call real_grow2_reverse(glueEZ, new_size)
         call real_grow(glueSurface, new_size)
         call real_grow(glueDiameter, new_size)
         call real_grow(glueBounding, new_size)

         glueQuat(redundant:new_size,:)=ZERO
         glueAngMom(redundant:new_size,:)=ZERO
         glueEX(redundant:new_size,:)=ZERO
         glueEY(redundant:new_size,:)=ZERO
         glueEZ(redundant:new_size,:)=ZERO
         glueSurface(redundant:new_size)=ZERO
         glueDiameter(redundant:new_size)=ZERO
         glueBounding(redundant:new_size)=ZERO
      ENDIF

        call real_grow2_reverse(DES_POS_NEW,new_size)
        call real_grow2_reverse(DES_VEL_NEW,new_size)
        call real_grow2_reverse(OMEGA_NEW,new_size)
        call real_grow2_reverse(PPOS,new_size)

! Force chain arrays are grown when needed in calc_force_dem.f
! Note that force chain arrays do not have the same size as particle arrays
! The number of active elements is the number of contacts between particle pairs
        ! call real_grow2_reverse(FCHAIN_MIDPOINT,new_size)
        ! FCHAIN_MIDPOINT(:,:) = UNDEFINED
        ! call real_grow2_reverse(FCHAIN_FN,new_size)
        ! call real_grow(FCHAIN_LENGTH,new_size)
        ! call real_grow(FCHAIN_FN_MAG,new_size)

! Residence time
        call real_grow(RESIDENCE_TIME,new_size)

        call byte_grow(PARTICLE_STATE,new_size)
        call integer_grow(iglobal_id,new_size)
        call integer_grow2_reverse(pijk,new_size)
        call integer_grow(dg_pijk,new_size)
        call integer_grow(dg_pijkprv,new_size)
        call logical_grow(ighost_updated,new_size)
        call real_grow2_reverse(FC,new_size)
        call real_grow2_reverse(TOW,new_size)
        call real_grow2_reverse(DES_COLL_FORCE,new_size)
        call real_grow(F_GP,new_size)
        call integer_grow2(WALL_COLLISION_FACET_ID,new_size)
        call real_grow3(WALL_COLLISION_PFT,new_size)
        call real_grow3(WALL_COLLISION_PFR,new_size)
        call real_grow3(WALL_COLLISION_FRAME,new_size)
        call real_grow2(WALL_COLLISION_PFTW,new_size)
        WALL_COLLISION_PFTW(:,old_size+1:new_size) = ZERO
        WALL_COLLISION_PFR(:,:,old_size+1:new_size) = ZERO
        WALL_COLLISION_FRAME(:,:,old_size+1:new_size) = ZERO
        lsz1 = size(JKR_WALL_CONTACT_STATE,1)
        lsz2 = size(JKR_WALL_CONTACT_STATE,2)
        allocate(jkrw_tmp(lsz1,new_size))
        jkrw_tmp(:,1:lsz2) = JKR_WALL_CONTACT_STATE(:,1:lsz2)
        jkrw_tmp(:,lsz2+1:) = .FALSE.
        call move_alloc(jkrw_tmp,JKR_WALL_CONTACT_STATE)
        call real_grow2_reverse(DRAG_FC,new_size)

        call integer_grow(NEIGHBOR_INDEX,new_size)
        call integer_grow(NEIGHBOR_INDEX_OLD,new_size)

        IF(PARTICLE_ORIENTATION) call real_grow2_reverse(ORIENTATION,new_size)

        IF(FILTER_SIZE > 0) THEN
           call integer_grow2(FILTER_CELL,new_size)
           call real_grow2(FILTER_WEIGHT,new_size)
        ENDIF

        IF(CGDEM) call real_grow(des_cgp_stw,new_size)
        IF(CGDEM) call real_grow(des_cgp_rpr,new_size)

        IF(MPPIC) THEN
           call real_grow(DES_STAT_WT,new_size)
           call real_grow2(PS_GRAD,new_size)
           call real_grow(EPG_P,new_size)
        ENDIF

        IF(USE_COHESION) THEN
           call real_grow(PostCohesive,new_size)
        ENDIF

        IF (DO_OLD) THEN
           call real_grow2_reverse(DES_POS_OLD,new_size)
           call real_grow2_reverse(DES_VEL_OLD,new_size)
           call real_grow2_reverse(DES_ACC_OLD,new_size)
           call real_grow2_reverse(OMEGA_OLD,new_size)
           call real_grow2_reverse(ROT_ACC_OLD,new_size)
        ENDIF

        IF(ENERGY_EQ)THEN
           call real_grow(DES_T_s,new_size)
           call real_grow(DES_C_PS,new_size)
           call real_grow2_reverse(DES_X_s,new_size)
           call real_grow(Q_Source,new_size)
           IF(CALC_CONV_DES) THEN
              call real_grow(CONV_Qs, new_size)
              call real_grow(GAMMAxSA, new_size)
           ENDIF
           IF(INTG_ADAMS_BASHFORTH) &
              call real_grow(Q_Source0,new_size)
        ENDIF

        IF(ANY_SOLIDS_SPECIES_EQ)THEN

           call real_grow2_reverse(DES_X_s,new_size)
           call real_grow2_reverse( DES_R_s, new_size )

           IF (INTG_ADAMS_BASHFORTH) THEN
              call real_grow( dMdt_OLD, new_size )
              call real_grow2_reverse( dXdt_OLD, new_size )
           ENDIF

           call real_grow( RXNS_Qs, new_size )
        ENDIF

        IF(SAVE_PART_RRATES) call real_grow2_reverse(PART_RRATES_OUT,new_size)

        IF(DES_USR_VAR_SIZE > 0) &
        call real_grow2(DES_USR_VAR,new_size)

        CALL DES_INIT_PARTICLE_ARRAYS(old_size+1,new_size)

        call real_grow( COORDINATION, new_size )
        call real_grow( MIN_OVERLAP, new_size )
        call real_grow( MAX_OVERLAP, new_size )
        call real_grow( MEAN_OVERLAP, new_size )

        RETURN

      END SUBROUTINE PARTICLE_GROW

!``````````````````````````````````````````````````````````````````````!
! Subroutine: PARTICLE_GROW_GluedSphereDEM                             !
!                                                                      !
! Author: Hang Zhou                                     Date: Dec. 2023!
! Purpose: Grow glued particle arrays to new_ngluedparticle.           !
! Note that pair max should be increased before calling this routine.  !
! Also, no assumption to the previous array size is made as needed for !
! restarts.                                                            !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE PARTICLE_GROW_GluedSphereDEM(new_ngluedparticle)

        USE discretelement
        USE param1, only: undefined, zero
        USE resize
        USE run

        IMPLICIT NONE

        integer, intent(in) :: new_ngluedparticle
        integer :: old_size, new_size
        integer :: redundant
        IF (new_ngluedparticle .le. size(glueMass)) RETURN

        old_size = size(glueMass)

        new_size = old_size

        DO WHILE (new_size < new_ngluedparticle)
           new_size = 2 * new_size
        ENDDO

        call integer_grow(particle_state_gsp,new_size)
        call real_grow(glueMass,new_size)
	  call real_grow(glueVolume,new_size)
        call real_grow(glueDiameter,new_size)
        call real_grow(glueBounding,new_size)
        call real_grow2_reverse(glueForce,new_size)
        call real_grow2_reverse(glueTorque,new_size)
        !call real_grow2_reverse(glueForcePP,new_size)
        !call real_grow2_reverse(glueTorquePP,new_size)
        call real_grow2_reverse(glueAcc,new_size)
        call real_grow2_reverse(glueAngMom,new_size)
        call real_grow2_reverse(glueQuat,new_size)
        call real_grow2_reverse(glueInertia,new_size)
        call real_grow2_reverse(glueOmg,new_size)
        call real_grow2_reverse(glueVel,new_size)
        call real_grow2_reverse(gluePos,new_size)
        !call real_grow2_reverse(glueContactPoint,new_size)
        call real_grow2_reverse(glueEX,new_size)
        call real_grow2_reverse(glueEY,new_size)
        call real_grow2_reverse(glueEZ,new_size)

        redundant = old_size + 1
        ! begin reset
        glueMass(redundant:new_size)=ZERO
	  glueVolume(redundant:new_size)=ZERO
        glueDiameter(redundant:new_size)=ZERO
        glueBounding(redundant:new_size)=ZERO
        glueForce(redundant:new_size,:)=ZERO
        glueTorque(redundant:new_size,:)=ZERO
        particle_state_gsp(redundant:new_size)=0
        !glueForcePP(redundant:new_size,:)=ZERO
        !glueTorquePP(redundant:new_size,:)=ZERO
        glueAcc(redundant:new_size,:)=ZERO
        glueAngMom(redundant:new_size,:)=ZERO
        glueQuat(redundant:new_size,:)=ZERO
        glueInertia(redundant:new_size,:)=ZERO
        glueOmg(redundant:new_size,:)=ZERO
        glueVel(redundant:new_size,:)=ZERO
        gluePos(redundant:new_size,:)=ZERO
        !glueContactPoint(redundant:new_size,:)=ZERO
        glueEX(redundant:new_size,:)=ZERO
        glueEY(redundant:new_size,:)=ZERO
        glueEZ(redundant:new_size,:)=ZERO

      END SUBROUTINE PARTICLE_GROW_GluedSphereDEM

    END MODULE DES_ALLOCATE
