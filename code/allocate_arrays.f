#include "error.inc"

MODULE ALLOCATE_ARRAYS_MOD

   use error_manager

CONTAINS

! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ALLOCATE_ARRAYS                                         C
!                                                                      C
!  Author: M. Syamlal                                Date: 17-DEC-98   C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

   SUBROUTINE ALLOCATE_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------

      use ambm
      use cdist
      use cont, only: do_cont
      use des_rxns
      use drag
      use energy
      use fldvar
      use generate_particles, only: particle_count
      use geometry
      use ghdtheory
      use indices
      use kintheory
      use mflux
      use param
      use param1
      use pgcor
      use physprop
      use pscor
      use residual
      use run
      use rxns
      use scalars
      use iterate, only: errorpercent
      use tau_g
      use tau_s
      use trace
      use turb
      use visc_g
      use visc_s
      use vshear

      IMPLICIT NONE

!-----------------------------------------------
! Variables
!-----------------------------------------------

!ambm
      if (allocated(A_m)) deallocate(A_m); allocate(A_m(DIMENSION_3, -3:3, 0:DIMENSION_M))
      if (allocated(B_m)) deallocate(B_m); allocate(B_m(DIMENSION_3, 0:DIMENSION_M))

!cont
      if (allocated(DO_CONT)) deallocate(DO_CONT); allocate(DO_CONT(0:DIMENSION_M))

!drag
      if (allocated(F_gs)) deallocate(F_gs); allocate(F_gs(DIMENSION_3, DIMENSION_M))
      if (allocated(F_ss)) deallocate(F_ss); allocate(F_ss(DIMENSION_3, 0:DIMENSION_LM))

!Off diagonal friction coefficient in HYS drag relation
      IF(DRAG_TYPE_ENUM.EQ.HYS) THEN
         if (allocated(beta_ij)) deallocate(beta_ij); allocate(beta_ij(DIMENSION_3, 0:DIMENSION_M, 0:DIMENSION_M))
      ENDIF

!energy
      if (allocated(HOR_g)) deallocate(HOR_g); allocate(HOR_g(DIMENSION_3))
      if (allocated(HOR_s)) deallocate(HOR_s); allocate(HOR_s(DIMENSION_3, DIMENSION_M))
      if (allocated(GAMA_gs)) deallocate(GAMA_gs); allocate(GAMA_gs(DIMENSION_3, DIMENSION_M))
      if (allocated(GAMA_Rg)) deallocate(GAMA_Rg); allocate(GAMA_Rg(DIMENSION_3))
      if (allocated(GAMA_Rs)) deallocate(GAMA_Rs); allocate(GAMA_Rs(DIMENSION_3, DIMENSION_M))
      if (allocated(T_Rg)) deallocate(T_Rg); allocate(T_Rg(DIMENSION_3))
      if (allocated(T_Rs)) deallocate(T_Rs); allocate(T_Rs(DIMENSION_3, DIMENSION_M))

!fldvar
      if (allocated(EP_g)) deallocate(EP_g); allocate(EP_g(DIMENSION_3))
      if (allocated(epg_jfac)) deallocate(epg_jfac); allocate(epg_jfac(DIMENSION_3p))
      if (allocated(epg_ifac)) deallocate(epg_ifac); allocate(epg_ifac(DIMENSION_3p))
      if (allocated(eps_ifac)) deallocate(eps_ifac); allocate(eps_ifac(DIMENSION_3p, DIMENSION_M))
      if (allocated(EP_go)) deallocate(EP_go); allocate(EP_go(DIMENSION_3p))
      if (allocated(P_g)) deallocate(P_g); allocate(P_g(DIMENSION_3))
      if (allocated(P_go)) deallocate(P_go); allocate(P_go(DIMENSION_3p))
      if (allocated(RO_g)) deallocate(RO_g); allocate(RO_g(DIMENSION_3))
      if (allocated(RO_go)) deallocate(RO_go); allocate(RO_go(DIMENSION_3p))
      if (allocated(ROP_g)) deallocate(ROP_g); allocate(ROP_g(DIMENSION_3))
      if (allocated(ROP_go)) deallocate(ROP_go); allocate(ROP_go(DIMENSION_3p))
      if (allocated(RO_S)) deallocate(RO_S); allocate(RO_S(DIMENSION_3, DIMENSION_M))
      if (allocated(RO_So)) deallocate(RO_So); allocate(RO_So(DIMENSION_3p, DIMENSION_M))
      if (allocated(ROP_s)) deallocate(ROP_s); allocate(ROP_s(DIMENSION_3, DIMENSION_M))
      if (allocated(ROP_so)) deallocate(ROP_so); allocate(ROP_so(DIMENSION_3p, DIMENSION_M))

      if (allocated(EP_SS)) deallocate(EP_SS); allocate(EP_SS(DIMENSION_3,DIMENSION_M,DIMENSION_N_S))
      if (allocated(ERR_ARRAY)) deallocate(ERR_ARRAY); allocate(ERR_ARRAY(DIMENSION_3,DIMENSION_M))

      if (allocated(T_g)) deallocate(T_g); allocate(T_g(DIMENSION_3))
      if (allocated(T_s)) deallocate(T_s); allocate(T_s(DIMENSION_3, DIMENSION_M))
      if (allocated(T_go)) deallocate(T_go); allocate(T_go(DIMENSION_3p))
      if (allocated(T_so)) deallocate(T_so); allocate(T_so(DIMENSION_3p, DIMENSION_M))
      if (allocated(X_g)) deallocate(X_g); allocate(X_g(DIMENSION_3, DIMENSION_N_g))
      if (allocated(X_s)) deallocate(X_s); allocate(X_s(DIMENSION_3, DIMENSION_M, DIMENSION_N_s))
      if (allocated(X_go)) deallocate(X_go); allocate(X_go(DIMENSION_3p, DIMENSION_N_g))
      if (allocated(X_so)) deallocate(X_so); allocate(X_so(DIMENSION_3p, DIMENSION_M, DIMENSION_N_s))
      if (allocated(U_g)) deallocate(U_g); allocate(U_g(DIMENSION_3))
      if (allocated(U_go)) deallocate(U_go); allocate(U_go(DIMENSION_3p))
      if (allocated(U_s)) deallocate(U_s); allocate(U_s(DIMENSION_3, DIMENSION_M))
      if (allocated(U_so)) deallocate(U_so); allocate(U_so(DIMENSION_3p, DIMENSION_M))
      if (allocated(V_g)) deallocate(V_g); allocate(V_g(DIMENSION_3))
      if (allocated(V_go)) deallocate(V_go); allocate(V_go(DIMENSION_3p))
      if (allocated(V_s)) deallocate(V_s); allocate(V_s(DIMENSION_3, DIMENSION_M))
      if (allocated(V_so)) deallocate(V_so); allocate(V_so(DIMENSION_3p, DIMENSION_M))
      if (allocated(W_g)) deallocate(W_g); allocate(W_g(DIMENSION_3))
      if (allocated(W_go)) deallocate(W_go); allocate(W_go   (DIMENSION_3p))
      if (allocated(W_s)) deallocate(W_s); allocate(W_s   (DIMENSION_3, DIMENSION_M))
      if (allocated(W_so)) deallocate(W_so); allocate(W_so   (DIMENSION_3p, DIMENSION_M))
      if (allocated(P_s)) deallocate(P_s); allocate(P_s   (DIMENSION_3, DIMENSION_M))
      if (allocated(P_s_c)) deallocate(P_s_c); allocate(P_s_c(DIMENSION_3, DIMENSION_M))
      if (allocated(P_s_v)) deallocate(P_s_v); allocate(P_s_v(DIMENSION_3))
      if (allocated(P_s_f)) deallocate(P_s_f); allocate(P_s_f(DIMENSION_3))
      if (allocated(P_s_p)) deallocate(P_s_p); allocate(P_s_p(DIMENSION_3))
      if (allocated(P_star)) deallocate(P_star); allocate(P_star(DIMENSION_3))
      if (allocated(P_staro)) deallocate(P_staro); allocate(P_staro(DIMENSION_3p))
      if (allocated(THETA_m)) deallocate(THETA_m); allocate(THETA_m(DIMENSION_3, DIMENSION_M))
      if (allocated(THETA_mo)) deallocate(THETA_mo); allocate(THETA_mo(DIMENSION_3p, DIMENSION_M))

      IF(K_Epsilon)THEN
        if (allocated(K_Turb_G)) deallocate(K_Turb_G); allocate(K_Turb_G(DIMENSION_3))
        if (allocated(K_Turb_Go)) deallocate(K_Turb_Go); allocate(K_Turb_Go(DIMENSION_3p))
        if (allocated(E_Turb_G)) deallocate(E_Turb_G); allocate(E_Turb_G(DIMENSION_3))
        if (allocated(E_Turb_Go)) deallocate(E_Turb_Go); allocate(E_Turb_Go(DIMENSION_3p))
      ENDIF

      IF(DIMENSION_Scalar /= 0) THEN
        if (allocated(Scalar)) deallocate(Scalar); allocate(Scalar(DIMENSION_3,  DIMENSION_Scalar))
        if (allocated(Scalaro)) deallocate(Scalaro); allocate(Scalaro(DIMENSION_3p, DIMENSION_Scalar))
      ENDIF


!pgcor
      if (allocated(d_e)) deallocate(d_e); allocate(d_e(DIMENSION_3p, 0:DIMENSION_M))
      if (allocated(d_n)) deallocate(d_n); allocate(d_n(DIMENSION_3p, 0:DIMENSION_M))
      if (allocated(d_t)) deallocate(d_t); allocate(d_t(DIMENSION_3p, 0:DIMENSION_M))
      if (allocated(Pp_g)) deallocate(Pp_g); allocate(Pp_g(DIMENSION_3p))
      if (allocated(PHASE_4_P_g)) deallocate(PHASE_4_P_g); allocate(PHASE_4_P_g(DIMENSION_3p))

!physprop
      if (allocated(MU_g)) deallocate(MU_g); allocate(MU_g(DIMENSION_3))
      if (allocated(C_pg)) deallocate(C_pg); allocate(C_pg(DIMENSION_3))
      if (allocated(C_ps)) deallocate(C_ps); allocate(C_ps(DIMENSION_3, DIMENSION_M))
      if (allocated(K_g)) deallocate(K_g); allocate(K_g(DIMENSION_3))
      if (allocated(K_s)) deallocate(K_s); allocate(K_s(DIMENSION_3, DIMENSION_M))
      if (allocated(Kth_s)) deallocate(Kth_s); allocate(Kth_s(DIMENSION_3, DIMENSION_M))
      if (allocated(Kphi_s)) deallocate(Kphi_s); allocate(Kphi_s(DIMENSION_3, DIMENSION_M))
      if (allocated(DIF_g)) deallocate(DIF_g); allocate(DIF_g(DIMENSION_3p, DIMENSION_N_g))
      if (allocated(divJi)) deallocate(divJi); allocate(divJi(DIMENSION_3p, (DIMENSION_N_g-1)))
      if (allocated(DIF_s)) deallocate(DIF_s); allocate(DIF_s(DIMENSION_3p, DIMENSION_M, DIMENSION_N_s))
      if (allocated(MW_MIX_g)) deallocate(MW_MIX_g); allocate(MW_MIX_g(DIMENSION_3))

!pscor
      if (allocated(e_e)) deallocate(e_e); allocate(e_e(DIMENSION_3p))
      if (allocated(e_n)) deallocate(e_n); allocate(e_n(DIMENSION_3p))
      if (allocated(e_t)) deallocate(e_t); allocate(e_t(DIMENSION_3p))
      if (allocated(K_cp)) deallocate(K_cp); allocate(K_cp(DIMENSION_3p))
      if (allocated(EPp)) deallocate(EPp); allocate(EPp(DIMENSION_3p))
      if (allocated(PHASE_4_P_s)) deallocate(PHASE_4_P_s); allocate(PHASE_4_P_s(DIMENSION_3p))

!residual
      if (allocated(RESID)) deallocate(RESID); allocate(RESID  (NRESID, 0:DIMENSION_M))
      if (allocated(MAX_RESID)) deallocate(MAX_RESID); allocate(MAX_RESID  (NRESID, 0:DIMENSION_M))
      if (allocated(IJK_RESID)) deallocate(IJK_RESID); allocate(IJK_RESID  (NRESID, 0:DIMENSION_M))
      if (allocated(NUM_RESID)) deallocate(NUM_RESID); allocate(NUM_RESID  (NRESID, 0:DIMENSION_M))
      if (allocated(DEN_RESID)) deallocate(DEN_RESID); allocate(DEN_RESID  (NRESID, 0:DIMENSION_M))
      if (allocated(RESID_PACK)) deallocate(RESID_PACK); allocate(RESID_PACK  (NRESID*2*(DIMENSION_M+1)))

!rxns
      if (nRR .gt. 0) then
         if (allocated(ReactionRates)) deallocate(ReactionRates); allocate(ReactionRates(DIMENSION_3,nRR))
      endif
      if (allocated(R_gp)) deallocate(R_gp); allocate(R_gp(DIMENSION_3p, DIMENSION_N_g))
      if (allocated(R_sp)) deallocate(R_sp); allocate(R_sp(DIMENSION_3p, DIMENSION_M, DIMENSION_N_s))
      if (allocated(RoX_gc)) deallocate(RoX_gc); allocate(RoX_gc(DIMENSION_3p, DIMENSION_N_g))
      if (allocated(RoX_sc)) deallocate(RoX_sc); allocate(RoX_sc(DIMENSION_3p, DIMENSION_M, DIMENSION_N_s))
      if (allocated(SUM_R_g)) deallocate(SUM_R_g); allocate(SUM_R_g(DIMENSION_3p))
      if (allocated(SUM_R_s)) deallocate(SUM_R_s); allocate(SUM_R_s(DIMENSION_3p, DIMENSION_M))
      if (allocated(R_phase)) deallocate(R_phase); allocate(R_phase(DIMENSION_3, DIMENSION_LM+DIMENSION_M-1))


! Allocate Fluid and DES reaction rates arrays used to write reaction rates in vtk files or monitors
      if (SAVE_FLUID_RRATES) then
         if(allocated(Fluid_RRates_out)) deallocate(Fluid_RRates_out)
         allocate(Fluid_RRates_out(DIMENSION_3,NO_OF_RXNS))
         Fluid_RRates_out(:,:) = ZERO
      endif
      if (SAVE_DES_RRATES) then
         if(allocated(Des_RRates_out)) deallocate(Des_RRates_out)
         allocate(Des_RRates_out(DIMENSION_3,NO_OF_DES_RXNS))
         Des_RRates_out(:,:) = ZERO
      endif

!scalars
      IF(DIMENSION_Scalar /= 0) then
        if (allocated(Scalar_c)) deallocate(Scalar_c); allocate(Scalar_c(DIMENSION_3p,  DIMENSION_Scalar))
        if (allocated(Scalar_p)) deallocate(Scalar_p); allocate(Scalar_p(DIMENSION_3p,  DIMENSION_Scalar))
        if (allocated(Dif_Scalar)) deallocate(Dif_Scalar); allocate(Dif_Scalar(DIMENSION_3p, DIMENSION_Scalar))
      ENDIF

! add by rong for dqmom
      if (allocated(D_p)) deallocate(D_p); allocate(D_p(DIMENSION_3, DIMENSION_M))
      if (allocated(D_po)) deallocate(D_po); allocate(D_po(DIMENSION_3, DIMENSION_M))
      if (allocated(Source_a)) deallocate(Source_a); allocate(Source_a(DIMENSION_3, DIMENSION_M))
      if (allocated(S_bar)) deallocate(S_bar); allocate(S_bar(0:DIM_Scalar2-1))
      if (allocated(Matrix_a)) deallocate(Matrix_a); allocate(Matrix_a(DIM_Scalar2,DIM_scalar2))
      if (allocated(Matrix_b)) deallocate(Matrix_b); allocate(Matrix_b(DIM_Scalar2,DIM_scalar2))
      if (allocated(Matrix_c)) deallocate(Matrix_c); allocate(Matrix_c(DIM_Scalar2,DIM_scalar2))
      if (allocated(Inv_a)) deallocate(Inv_a); allocate(Inv_a(DIM_Scalar2,DIM_scalar2))
      if (allocated(A)) deallocate(A); allocate(A  (1:DIMENSION_Scalar))
      if (allocated(omega)) deallocate(omega); allocate(omega  (1:DIMENSION_m))
      if (allocated(beta_a)) deallocate(beta_a); allocate(beta_a(DIM_Scalar,DIM_Scalar))
      if (allocated(ystart)) deallocate(ystart); allocate(ystart(1:DIM_Scalar))

! K-Epsilon Turbulence model
      IF(K_Epsilon) THEN
        if (allocated(K_Turb_G_c)) deallocate(K_Turb_G_c); allocate(K_Turb_G_c(DIMENSION_3p))
        if (allocated(K_Turb_G_p)) deallocate(K_Turb_G_p); allocate(K_Turb_G_p(DIMENSION_3p))
        if (allocated(Dif_K_Turb_G)) deallocate(Dif_K_Turb_G); allocate(Dif_K_Turb_G(DIMENSION_3p))
        if (allocated(E_Turb_G_c)) deallocate(E_Turb_G_c); allocate(E_Turb_G_c(DIMENSION_3p))
        if (allocated(E_Turb_G_p)) deallocate(E_Turb_G_p); allocate(E_Turb_G_p(DIMENSION_3p))
        if (allocated(Dif_E_Turb_G)) deallocate(Dif_E_Turb_G); allocate(Dif_E_Turb_G(DIMENSION_3p))
      ENDIF

! Simonin or Ahmadi model
      IF(KT_TYPE_ENUM==SIMONIN_1996 .OR.&
         KT_TYPE_ENUM==AHMADI_1995) THEN
        if (allocated(K_12)) deallocate(K_12); allocate(K_12(DIMENSION_3))
        if (allocated(Tau_12)) deallocate(Tau_12); allocate(Tau_12(DIMENSION_3))
        if (allocated(Tau_1)) deallocate(Tau_1); allocate(Tau_1(DIMENSION_3))
      ENDIF

!tau_g
      if (allocated(TAU_U_g)) deallocate(TAU_U_g); allocate(TAU_U_g(DIMENSION_3p))
      if (allocated(TAU_V_g)) deallocate(TAU_V_g); allocate(TAU_V_g(DIMENSION_3p))
      if (allocated(TAU_W_g)) deallocate(TAU_W_g); allocate(TAU_W_g(DIMENSION_3p))
      if (allocated(DF_gu)) deallocate(DF_gu); allocate(DF_gu(DIMENSION_3p, -3:3))
      if (allocated(DF_gv)) deallocate(DF_gv); allocate(DF_gv(DIMENSION_3p, -3:3))
      if (allocated(DF_gw)) deallocate(DF_gw); allocate(DF_gw(DIMENSION_3p, -3:3))

! JFD TRY
      DF_GU(:,:) = 0.0D0
      DF_GV(:,:) = 0.0D0
      DF_GW(:,:) = 0.0D0

      if (allocated(CTAU_U_G)) deallocate(CTAU_U_G); allocate(CTAU_U_G(DIMENSION_3P))
      if (allocated(CTAU_V_G)) deallocate(CTAU_V_G); allocate(CTAU_V_G(DIMENSION_3P))
      if (allocated(CTAU_W_G)) deallocate(CTAU_W_G); allocate(CTAU_W_G(DIMENSION_3P))

!tau_s
      if (allocated(TAU_U_s)) deallocate(TAU_U_s); allocate(TAU_U_s(DIMENSION_3p, DIMENSION_M))
      if (allocated(TAU_V_s)) deallocate(TAU_V_s); allocate(TAU_V_s(DIMENSION_3p, DIMENSION_M))
      if (allocated(TAU_W_s)) deallocate(TAU_W_s); allocate(TAU_W_s(DIMENSION_3p, DIMENSION_M))

! generate_particles / particle_count
      if (allocated(PARTICLE_COUNT)) deallocate(PARTICLE_COUNT); allocate(PARTICLE_COUNT(DIMENSION_3))

!trace
      if (allocated(trD_s_C)) deallocate(trD_s_C); allocate(trD_s_C(DIMENSION_3, DIMENSION_M))
      if (allocated(trD_s2)) deallocate(trD_s2); allocate(trD_s2(DIMENSION_3, DIMENSION_M))
      if (allocated(trD_s_Co)) deallocate(trD_s_Co); allocate(trD_s_Co(DIMENSION_3, DIMENSION_M))
      if (allocated(trD_s_Co2)) deallocate(trD_s_Co2); allocate(trD_s_Co2(DIMENSION_3, DIMENSION_M))
!visc_g
      if (allocated(trD_g)) deallocate(trD_g); allocate(trD_g(DIMENSION_3))
      if (allocated(MU_gt)) deallocate(MU_gt); allocate(MU_gt(DIMENSION_3))
      if (allocated(EPMU_gt)) deallocate(EPMU_gt); allocate(EPMU_gt(DIMENSION_3p))
      if (allocated(LAMBDA_gt)) deallocate(LAMBDA_gt); allocate(LAMBDA_gt(DIMENSION_3p))
      if (allocated(EPLAMBDA_gt)) deallocate(EPLAMBDA_gt); allocate(EPLAMBDA_gt(DIMENSION_3))
      if (allocated(L_scale)) deallocate(L_scale); allocate(L_scale(DIMENSION_3))

!visc_s
      if (allocated(MU_s)) deallocate(MU_s); allocate(MU_s(DIMENSION_3, DIMENSION_M))
      if (allocated(EPMU_s)) deallocate(EPMU_s); allocate(EPMU_s(DIMENSION_3p, DIMENSION_M))
      if (allocated(LAMBDA_s)) deallocate(LAMBDA_s); allocate(LAMBDA_s(DIMENSION_3, DIMENSION_M))
      if (allocated(EPLAMBDA_s)) deallocate(EPLAMBDA_s); allocate(EPLAMBDA_s(DIMENSION_3p, DIMENSION_M))
      if (allocated(ALPHA_s)) deallocate(ALPHA_s); allocate(ALPHA_s(DIMENSION_3, DIMENSION_M))
      if (allocated(MU_s_c)) deallocate(MU_s_c); allocate(MU_s_c(DIMENSION_3, DIMENSION_M))
      if (allocated(LAMBDA_s_c)) deallocate(LAMBDA_s_c); allocate(LAMBDA_s_c(DIMENSION_3, DIMENSION_M))
      if (allocated(LAMBDA_s_v)) deallocate(LAMBDA_s_v); allocate(LAMBDA_s_v(DIMENSION_3))
      if (allocated(LAMBDA_s_f)) deallocate(LAMBDA_s_f); allocate(LAMBDA_s_f(DIMENSION_3))
      if (allocated(LAMBDA_s_p)) deallocate(LAMBDA_s_p); allocate(LAMBDA_s_p(DIMENSION_3))
      if (allocated(MU_s_v)) deallocate(MU_s_v); allocate(MU_s_v(DIMENSION_3))
      if (allocated(MU_s_f)) deallocate(MU_s_f); allocate(MU_s_f(DIMENSION_3))
      if (allocated(MU_s_p)) deallocate(MU_s_p); allocate(MU_s_p(DIMENSION_3))
      if (allocated(MU_b_v)) deallocate(MU_b_v); allocate(MU_b_v(DIMENSION_3))
      if (allocated(EP_star_array)) deallocate(EP_star_array); allocate(EP_star_array(DIMENSION_3))
      if (allocated(EP_g_blend_start)) deallocate(EP_g_blend_start); allocate(EP_g_blend_start(DIMENSION_3))
      if (allocated(EP_g_blend_end)) deallocate(EP_g_blend_end); allocate(EP_g_blend_end(DIMENSION_3))
      if (allocated(trD_s)) deallocate(trD_s); allocate(trD_s(DIMENSION_3, DIMENSION_M))
      if (allocated(I2_devD_s)) deallocate(I2_devD_s); allocate(I2_devD_s(DIMENSION_3))
      if (allocated(TrM_s)) deallocate(TrM_s); allocate(TrM_s(DIMENSION_3))
      if (allocated(TrDM_s)) deallocate(TrDM_s); allocate(TrDM_s(DIMENSION_3))

!shear quantities
      if (allocated(VSH)) deallocate(VSH); allocate(VSH(DIMENSION_3))
      if (allocated(VSHE)) deallocate(VSHE); allocate(VSHE(DIMENSION_3))

!mflux
      if (allocated(Flux_gE)) deallocate(Flux_gE); allocate(Flux_gE(DIMENSION_3p))
      if (allocated(Flux_sE)) deallocate(Flux_sE); allocate(Flux_sE(DIMENSION_3p, DIMENSION_M))
      if (allocated(Flux_gN)) deallocate(Flux_gN); allocate(Flux_gN(DIMENSION_3p))
      if (allocated(Flux_sN)) deallocate(Flux_sN); allocate(Flux_sN(DIMENSION_3p, DIMENSION_M))
      if (allocated(Flux_gT)) deallocate(Flux_gT); allocate(Flux_gT(DIMENSION_3p))
      if (allocated(Flux_sT)) deallocate(Flux_sT); allocate(Flux_sT(DIMENSION_3p, DIMENSION_M))
      IF(ADDED_MASS) THEN ! Fluxes calculated for just one 'bubble' species (M=M_AM)
         if (allocated(Flux_gSE)) deallocate(Flux_gSE); allocate(Flux_gSE(DIMENSION_3p))
         if (allocated(Flux_sSE)) deallocate(Flux_sSE); allocate(Flux_sSE(DIMENSION_3p))
         if (allocated(Flux_gSN)) deallocate(Flux_gSN); allocate(Flux_gSN(DIMENSION_3p))
         if (allocated(Flux_sSN)) deallocate(Flux_sSN); allocate(Flux_sSN(DIMENSION_3p))
         if (allocated(Flux_gST)) deallocate(Flux_gST); allocate(Flux_gST(DIMENSION_3p))
         if (allocated(Flux_sST)) deallocate(Flux_sST); allocate(Flux_sST(DIMENSION_3p))
      ENDIF
      if (allocated(ROP_gE)) deallocate(ROP_gE); allocate(ROP_gE(DIMENSION_3p))
      if (allocated(ROP_sE)) deallocate(ROP_sE); allocate(ROP_sE(DIMENSION_3p, DIMENSION_M))
      if (allocated(ROP_gN)) deallocate(ROP_gN); allocate(ROP_gN(DIMENSION_3p))
      if (allocated(ROP_sN)) deallocate(ROP_sN); allocate(ROP_sN(DIMENSION_3p, DIMENSION_M))
      if (allocated(ROP_gT)) deallocate(ROP_gT); allocate(ROP_gT(DIMENSION_3p))
      if (allocated(ROP_sT)) deallocate(ROP_sT); allocate(ROP_sT(DIMENSION_3p, DIMENSION_M))

! allocate variables for GHD Theory
      IF (KT_TYPE_ENUM == GHD_2007) THEN
        if (allocated(Flux_nE)) deallocate(Flux_nE); allocate(Flux_nE(DIMENSION_3p))
        if (allocated(Flux_nN)) deallocate(Flux_nN); allocate(Flux_nN(DIMENSION_3p))
        if (allocated(Flux_nT)) deallocate(Flux_nT); allocate(Flux_nT(DIMENSION_3p))
        if (allocated(Zeta0)) deallocate(Zeta0)
        allocate(Zeta0(DIMENSION_3p))   ! zeroth rate of cooling
        if (allocated(ZetaU)) deallocate(ZetaU)
        allocate(ZetaU(DIMENSION_3p))   ! 1st order cooling rate transport coefficient
        if (allocated(DiT)) deallocate(DiT)
        allocate(DiT(DIMENSION_3p, DIMENSION_M))   ! thermal diffusivity
        if (allocated(DijF)) deallocate(DijF)
        allocate(DijF(DIMENSION_3p, DIMENSION_M, DIMENSION_M))   ! mass mobility
        if (allocated(Lij)) deallocate(Lij)
        allocate(Lij(DIMENSION_3p, DIMENSION_M, DIMENSION_M))   ! thermal mobility
        if (allocated(Dij)) deallocate(Dij)
        allocate(Dij(DIMENSION_3p, DIMENSION_M, DIMENSION_M))   ! ordinary diffusion
        if (allocated(DijQ)) deallocate(DijQ)
        allocate(DijQ(DIMENSION_3p, DIMENSION_M, DIMENSION_M))   ! Dufour coeff.
        if (allocated(JoiX)) deallocate(JoiX)
        allocate(JoiX(DIMENSION_3p, DIMENSION_M))   ! X- species mass flux
        if (allocated(JoiY)) deallocate(JoiY)
        allocate(JoiY(DIMENSION_3p, DIMENSION_M))   ! Y- species mass flux
        if (allocated(JoiZ)) deallocate(JoiZ)
        allocate(JoiZ(DIMENSION_3p, DIMENSION_M))   ! Z- species mass flux
        if (allocated(FiX)) deallocate(FiX)
        allocate(FiX(DIMENSION_3p, DIMENSION_M))   ! X- external force
        if (allocated(FiY)) deallocate(FiY)
        allocate(FiY(DIMENSION_3p, DIMENSION_M))   ! Y- external force
        if (allocated(FiZ)) deallocate(FiZ)
        allocate(FiZ(DIMENSION_3p, DIMENSION_M))   ! Z- external force
        if (allocated(FiXvel)) deallocate(FiXvel)
        allocate(FiXvel(DIMENSION_3p, DIMENSION_M))   ! X- external force
        if (allocated(FiYvel)) deallocate(FiYvel)
        allocate(FiYvel(DIMENSION_3p, DIMENSION_M))   ! Y- external force
        if (allocated(FiZvel)) deallocate(FiZvel)
        allocate(FiZvel(DIMENSION_3p, DIMENSION_M))   ! Z- external force
        if (allocated(DELTAU)) deallocate(DELTAU)
        allocate(DELTAU(DIMENSION_3p, DIMENSION_M))
        if (allocated(DELTAV)) deallocate(DELTAV)
        allocate(DELTAV(DIMENSION_3p, DIMENSION_M))
        if (allocated(DELTAW)) deallocate(DELTAW)
        allocate(DELTAW(DIMENSION_3p, DIMENSION_M))
        if (allocated(dragFx)) deallocate(dragFx)
        allocate(dragFx(DIMENSION_3p, DIMENSION_M))   ! X- drag force
        if (allocated(dragFy)) deallocate(dragFy)
        allocate(dragFy(DIMENSION_3p, DIMENSION_M))   ! Y- drag force
        if (allocated(dragFz)) deallocate(dragFz)
        allocate(dragFz(DIMENSION_3p, DIMENSION_M))   ! Z- drag force
        if (allocated(dragFxflux)) deallocate(dragFxflux)
        allocate(dragFxflux(DIMENSION_3p, DIMENSION_M))   ! X- drag force
        if (allocated(dragFyflux)) deallocate(dragFyflux)
        allocate(dragFyflux(DIMENSION_3p, DIMENSION_M))   ! Y- drag force
        if (allocated(dragFzflux)) deallocate(dragFzflux)
        allocate(dragFzflux(DIMENSION_3p, DIMENSION_M))   ! Z- drag force
        if (allocated(FiMinusDragX)) deallocate(FiMinusDragX)
        allocate(FiMinusDragX(DIMENSION_3p, DIMENSION_M))   ! X- drag force
        if (allocated(JoiMinusDragX)) deallocate(JoiMinusDragX)
        allocate(JoiMinusDragX(DIMENSION_3p, DIMENSION_M))   ! X- drag force
        if (allocated(FiMinusDragY)) deallocate(FiMinusDragY)
        allocate(FiMinusDragY(DIMENSION_3p, DIMENSION_M))   ! Y- drag force
        if (allocated(JoiMinusDragY)) deallocate(JoiMinusDragY)
        allocate(JoiMinusDragY(DIMENSION_3p, DIMENSION_M))   ! Y- drag force
        if (allocated(FiMinusDragZ)) deallocate(FiMinusDragZ)
        allocate(FiMinusDragZ(DIMENSION_3p, DIMENSION_M))   ! Z- drag force
        if (allocated(JoiMinusDragZ)) deallocate(JoiMinusDragZ)
        allocate(JoiMinusDragZ(DIMENSION_3p, DIMENSION_M))   ! Z- drag force
        if (allocated(beta_cell_X)) deallocate(beta_cell_X)
        allocate(beta_cell_X(DIMENSION_3p, DIMENSION_M))   ! X- drag force
        if (allocated(beta_cell_Y)) deallocate(beta_cell_Y)
        allocate(beta_cell_Y(DIMENSION_3p, DIMENSION_M))   ! Y- drag force
        if (allocated(beta_cell_Z)) deallocate(beta_cell_Z)
        allocate(beta_cell_Z(DIMENSION_3p, DIMENSION_M))   ! Y- drag force
        if (allocated(beta_ij_cell_X)) deallocate(beta_ij_cell_X)
        allocate(beta_ij_cell_X(DIMENSION_3p, DIMENSION_M,DIMENSION_M))   ! X- drag force
        if (allocated(beta_ij_cell_Y)) deallocate(beta_ij_cell_Y)
        allocate(beta_ij_cell_Y(DIMENSION_3p, DIMENSION_M,DIMENSION_M))   ! Y- drag force
        if (allocated(beta_ij_cell_Z)) deallocate(beta_ij_cell_Z)
        allocate(beta_ij_cell_Z(DIMENSION_3p, DIMENSION_M,DIMENSION_M))   ! Y- drag force
        if (allocated(DEL_DOT_J)) deallocate(DEL_DOT_J)
        allocate(DEL_DOT_J(DIMENSION_3p, DIMENSION_M))
        if (allocated(DiT_HarmE)) deallocate(DiT_HarmE)
        allocate(DiT_HarmE(DIMENSION_3p))
        if (allocated(DiT_HarmN)) deallocate(DiT_HarmN)
        allocate(DiT_HarmN(DIMENSION_3p))
        if (allocated(DiT_HarmT)) deallocate(DiT_HarmT)
        allocate(DiT_HarmT(DIMENSION_3p))
        if (allocated(Dij_HarmE)) deallocate(Dij_HarmE)
        allocate(Dij_HarmE(DIMENSION_3p, DIMENSION_M))
        if (allocated(Dij_HarmN)) deallocate(Dij_HarmN)
        allocate(Dij_HarmN(DIMENSION_3p, DIMENSION_M))
        if (allocated(Dij_HarmT)) deallocate(Dij_HarmT)
        allocate(Dij_HarmT(DIMENSION_3p, DIMENSION_M))
        if (allocated(DijF_HarmE)) deallocate(DijF_HarmE)
        allocate(DijF_HarmE(DIMENSION_3p, DIMENSION_M))
        if (allocated(DijF_HarmN)) deallocate(DijF_HarmN)
        allocate(DijF_HarmN(DIMENSION_3p, DIMENSION_M))
        if (allocated(DijF_HarmT)) deallocate(DijF_HarmT)
        allocate(DijF_HarmT(DIMENSION_3p, DIMENSION_M))
      ENDIF


! We need to set this even when KT_TYPE is not set to IA_NONEP - at
! least in the current version of the code and needs to be revisited
      if (allocated(KTMOM_U_s)) deallocate(KTMOM_U_s); allocate(KTMOM_U_s(DIMENSION_3p, DIMENSION_M))
      if (allocated(KTMOM_V_s)) deallocate(KTMOM_V_s); allocate(KTMOM_V_s(DIMENSION_3p, DIMENSION_M))
      if (allocated(KTMOM_W_s)) deallocate(KTMOM_W_s); allocate(KTMOM_W_s(DIMENSION_3p, DIMENSION_M))

! allocate variables for Iddir & Arastoopour (2005) kinetic theory
! EDvel_sM_ip & EDT_s_ip are also used for Garzy & Dufty (1999) kinetic theory
      IF (KT_TYPE_ENUM == IA_2005) THEN
         if (allocated(trD_s2_ip)) deallocate(trD_s2_ip); allocate(trD_s2_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(MU_sM_ip)) deallocate(MU_sM_ip); allocate(MU_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(MU_sL_ip)) deallocate(MU_sL_ip); allocate(MU_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(XI_sM_ip)) deallocate(XI_sM_ip); allocate(XI_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(XI_sL_ip)) deallocate(XI_sL_ip); allocate(XI_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(Fnu_s_ip)) deallocate(Fnu_s_ip); allocate(Fnu_s_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M))
         if (allocated(FT_sM_ip)) deallocate(FT_sM_ip); allocate(FT_sM_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M))
         if (allocated(FT_sL_ip)) deallocate(FT_sL_ip); allocate(FT_sL_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M))
         if (allocated(Kth_sL_ip)) deallocate(Kth_sL_ip); allocate(Kth_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(Knu_sM_ip)) deallocate(Knu_sM_ip); allocate(Knu_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(Knu_sL_ip)) deallocate(Knu_sL_ip); allocate(Knu_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(Kvel_s_ip)) deallocate(Kvel_s_ip); allocate(Kvel_s_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M))
         if (allocated(EDvel_sL_ip)) deallocate(EDvel_sL_ip); allocate(EDvel_sL_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M))
         if (allocated(ED_ss_ip)) deallocate(ED_ss_ip); allocate(ED_ss_ip(DIMENSION_3p, 0:DIMENSION_LM))
      ENDIF
      IF (KT_TYPE_ENUM == GTSH_2012) THEN
         if (allocated(A2_gtsh)) deallocate(A2_gtsh); allocate(A2_gtsh(DIMENSION_3))
         if (allocated(xsi_gtsh)) deallocate(xsi_gtsh); allocate(xsi_gtsh(DIMENSION_3))
      ENDIF
      IF (KT_TYPE_ENUM == IA_2005 .OR. &
          KT_TYPE_ENUM == GD_1999 .OR. &
          KT_TYPE_ENUM == GTSH_2012) THEN
         if (allocated(EDT_s_ip)) deallocate(EDT_s_ip); allocate(EDT_s_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M))
         if (allocated(EDvel_sM_ip)) deallocate(EDvel_sM_ip); allocate(EDvel_sM_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M))
      ENDIF

      if (allocated(errorpercent)) deallocate(errorpercent); allocate(errorpercent(0:MMAX))

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_ARRAYS_GEOMETRY                               !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Domain decomposition and dimensions
      use geometry, only: oDX, oDX_E
      use geometry, only: oDZ, oDZ_T
      use geometry, only: oDY, oDY_N
      use geometry, only: X, X_E, oX, oX_E, cyl_X, cyl_X_E
      use geometry, only: Z, Z_T
! Averaging factors.
      use geometry, only: FX_E, FX_E_bar, FX, FX_bar
      use geometry, only: FY_N, FY_N_bar
      use geometry, only: FZ_T, FZ_T_bar
! Domain flags.
      use geometry, only: ICBC_FLAG
      use geometry, only: FLAG
      use geometry, only: FLAG_E, FLAG_N, FLAG_T
! Domain volumes and areas.
      use geometry, only: VOL, VOL_SURR, AYZ, AXZ, AXY, CENTER_S! Scalar grid
      use geometry, only: VOL_U, AYZ_U, AXZ_U, AXY_U  ! X-Momentum
      use geometry, only: VOL_V, AYZ_V, AXZ_V, AXY_V  ! Y-Momentum
      use geometry, only: VOL_W, AYZ_W, AXZ_W, AXY_W  ! Z-Momentum
! Cell Aspect Ratio
      use geometry, only: Aspect_Ratio, Aspect_Ratio_U, Aspect_Ratio_V, Aspect_Ratio_W
! Axis decomposition
      USE param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      USE param, only: DIMENSION_3
      USE param, only: DIMENSION_3L, DIMENSION_3P
! Flag for POSTMFIX
      use cdist, only: bDoing_postmfix

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use compar, only:ADJUST_PARTITION

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Error Flag
      INTEGER :: IER = 0
! Flag indicating that the arrays were previously allocated.
      INTEGER, SAVE :: CALLED = -1
!......................................................................!

      IF(ADJUST_PARTITION) CALLED = -1

      CALLED = CALLED + 1

      IF(CALLED > 0) THEN
         IF(.NOT.bDoing_postmfix) THEN
            RETURN
         ELSEIF(mod(CALLED,2) /= 0) THEN
            RETURN
         ENDIF
      ENDIF

! ALLOC geometry components related to the mesh. Check the
! allocation error status and abort if any failure is detected.
      ALLOCATE(X(0:DIMENSION_I))
      ALLOCATE(cyl_X(0:DIMENSION_I))
      ALLOCATE(X_E(0:DIMENSION_I))
      ALLOCATE(cyl_X_E(0:DIMENSION_I))
      ALLOCATE(oX(0:DIMENSION_I))
      ALLOCATE(oX_E(0:DIMENSION_I))
      ALLOCATE(oDX(0:DIMENSION_I))
      ALLOCATE(oDX_E(0:DIMENSION_I))
      IF(IER /= 0) goto 500

      ALLOCATE(oDY(0:DIMENSION_J))
      ALLOCATE(oDY_N(0:DIMENSION_J))
      IF(IER /= 0) goto 500

      ALLOCATE(Z(0:DIMENSION_K))
      ALLOCATE(Z_T(0:DIMENSION_K))
      ALLOCATE(oDZ(0:DIMENSION_K))
      ALLOCATE(oDZ_T(0:DIMENSION_K))
      IF(IER /= 0) goto 500

      ALLOCATE(FX(0:DIMENSION_I))
      ALLOCATE(FX_bar(0:DIMENSION_I))
      IF(IER /= 0) goto 500

      ALLOCATE(FX_E(0:DIMENSION_I))
      ALLOCATE(FX_E_bar(0:DIMENSION_I))
      IF(IER /= 0) goto 500

      ALLOCATE(FY_N(0:DIMENSION_J))
      ALLOCATE(FY_N_bar(0:DIMENSION_J))
      IF(IER /= 0) goto 500

      ALLOCATE(FZ_T(0:DIMENSION_K))
      ALLOCATE(FZ_T_bar(0:DIMENSION_K))
      IF(IER /= 0) goto 500

! Flags for the scalar grid.
      if (allocated(FLAG)) deallocate(FLAG); allocate(FLAG(DIMENSION_3))
      IF(IER /= 0) goto 500

! Flags for the momentum grids.
      if (allocated(FLAG_E)) deallocate(FLAG_E); allocate(FLAG_E(DIMENSION_3))
      if (allocated(FLAG_N)) deallocate(FLAG_N); allocate(FLAG_N(DIMENSION_3))
      if (allocated(FLAG_T)) deallocate(FLAG_T); allocate(FLAG_T(DIMENSION_3))
      IF(IER /= 0) goto 500

! Text flags for scalar grid.
      Allocate(ICBC_FLAG(DIMENSION_3L))
      ! if (allocated(ICBC_FLAG)) deallocate(ICBC_FLAG); allocate(ICBC_FLAG(DIMENSION_3L))
      IF(IER /= 0) goto 500

! Volume and face-areas of scalar grid.
      if (allocated(VOL)) deallocate(VOL); allocate(VOL(DIMENSION_3))
      if (allocated(AYZ)) deallocate(AYZ); allocate(AYZ(DIMENSION_3P))
      if (allocated(AXZ)) deallocate(AXZ); allocate(AXZ(DIMENSION_3P))
      if (allocated(AXY)) deallocate(AXY); allocate(AXY(DIMENSION_3P))

      if (allocated(CENTER_S)) deallocate(CENTER_S); allocate(CENTER_S(DIMENSION_3P,3))
      IF(IER /= 0) goto 500

      ! total volume of each cell's surrounding stencil cells
      if (allocated(VOL_SURR)) deallocate(VOL_SURR); allocate(VOL_SURR(DIMENSION_3))

! Volume and face-areas of X-Momentumn grid.
      if (allocated(VOL_U)) deallocate(VOL_U); allocate(VOL_U(DIMENSION_3))
      if (allocated(AYZ_U)) deallocate(AYZ_U); allocate(AYZ_U(DIMENSION_3P))
      if (allocated(AXZ_U)) deallocate(AXZ_U); allocate(AXZ_U(DIMENSION_3P))
      if (allocated(AXY_U)) deallocate(AXY_U); allocate(AXY_U(DIMENSION_3P))
      IF(IER /= 0) goto 500

! Volume and face-areas of Y-Momentum grid.
      if (allocated(VOL_V)) deallocate(VOL_V); allocate(VOL_V(DIMENSION_3))
      if (allocated(AYZ_V)) deallocate(AYZ_V); allocate(AYZ_V(DIMENSION_3P))
      if (allocated(AXZ_V)) deallocate(AXZ_V); allocate(AXZ_V(DIMENSION_3P))
      if (allocated(AXY_V)) deallocate(AXY_V); allocate(AXY_V(DIMENSION_3P))
      IF(IER /= 0) goto 500

! Volume and face-areas of Z-Momentum grid.
      if (allocated(VOL_W)) deallocate(VOL_W); allocate(VOL_W(DIMENSION_3))
      if (allocated(AYZ_W)) deallocate(AYZ_W); allocate(AYZ_W(DIMENSION_3P))
      if (allocated(AXZ_W)) deallocate(AXZ_W); allocate(AXZ_W(DIMENSION_3P))
      if (allocated(AXY_W)) deallocate(AXY_W); allocate(AXY_W(DIMENSION_3P))
      IF(IER /= 0) goto 500

! Cell Aspect Ratio
      if (allocated(Aspect_Ratio))   deallocate(Aspect_Ratio);   allocate(Aspect_Ratio(DIMENSION_3))
      if (allocated(Aspect_Ratio_U)) deallocate(Aspect_Ratio_U); allocate(Aspect_Ratio_U(DIMENSION_3))
      if (allocated(Aspect_Ratio_V)) deallocate(Aspect_Ratio_V); allocate(Aspect_Ratio_V(DIMENSION_3))
      if (allocated(Aspect_Ratio_W)) deallocate(Aspect_Ratio_W); allocate(Aspect_Ratio_W(DIMENSION_3))

! Collect the error flags from all ranks. If all allocations were
! successful, do nothing. Otherwise, flag the error and abort.
! Note that the allocation status is checked in groups. This can
! be increase if tracking the source of an allocation failure.
  500 CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL LOG_ERROR()
      ENDIF

 1100 FORMAT('Error 1100: Failure during array allocation.')

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutinee: DEALLOCATE_ARRAYS_                                     C
!  Purpose: Deallocate arrays                                          C
!                                                                      C
!  Author: Jeff Dietiker                             Date: 03-MAR-2016 C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DEALLOCATE_ARRAYS  !@@

         USE set_increments_mod, only: DEALLOCATE_ARRAYS_INCREMENTS

      CALL DEALLOCATE_ARRAYS_MAIN
      CALL DEALLOCATE_ARRAYS_GEOMETRY
      CALL DEALLOCATE_ARRAYS_INCREMENTS
      CALL DEALLOCATE_ARRAYS_PARALLEL
      CALL DEALLOCATE_CUT_CELL_ARRAYS
      CALL DEALLOCATE_VTU_ARRAYS
      CALL DEALLOCATE_DEM_MI
      CALL DEALLOCATE_PIC_MI
      CALL DES_DEALLOCATE_ARRAYS


      END SUBROUTINE DEALLOCATE_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutinee: DEALLOCATE_ARRAYS_MAIN                                 C
!  Purpose: Deallocate arrays                                          C
!                                                                      C
!  Author: Jeff Dietiker                             Date: 03-MAR-2016 C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DEALLOCATE_ARRAYS_MAIN

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use ambm
      use cdist
      use cont, only: do_cont
      use des_rxns
      use drag
      use energy
      use fldvar
      use generate_particles, only: particle_count
      use geometry
      use ghdtheory
      use indices
      use kintheory
      use mflux
      use param
      use param1
      use pgcor
      use physprop
      use pscor
      use residual
      use run
      use rxns
      use scalars
      use iterate, only: errorpercent
      use tau_g
      use tau_s
      use trace
      use turb
      use visc_g
      use visc_s
      use vshear

      IMPLICIT NONE

!-----------------------------------------------
! Variables
!-----------------------------------------------

!ambm
      if(allocated(A_m)) deallocate(A_m)
      if(allocated(B_m)) deallocate(B_m)

!cont
      if(allocated(DO_CONT)) deallocate(DO_CONT)

!drag
      if(allocated(F_gs)) deallocate(F_gs)
      if(allocated(F_ss)) deallocate(F_ss)


!Off diagonal friction coefficient in HYS drag relation
      if(allocated(beta_ij)) deallocate(beta_ij)


!energy
      if(allocated(HOR_g)) deallocate(HOR_g)
      if(allocated(HOR_s)) deallocate(HOR_s)
      if(allocated(GAMA_gs)) deallocate(GAMA_gs)
      if(allocated(GAMA_Rg)) deallocate(GAMA_Rg)
      if(allocated(GAMA_Rs)) deallocate(GAMA_Rs)
      if(allocated(T_Rg)) deallocate(T_Rg)
      if(allocated(T_Rs)) deallocate(T_Rs)


!fldvar

      if(allocated(EP_g)) deallocate(EP_g)
      if(allocated(epg_jfac)) deallocate(epg_jfac)
      if(allocated(epg_ifac)) deallocate(epg_ifac)
      if(allocated(eps_ifac)) deallocate(eps_ifac)
      if(allocated(EP_go)) deallocate(EP_go)
      if(allocated(P_g)) deallocate(P_g)
      if(allocated(P_go)) deallocate(P_go)
      if(allocated(RO_g)) deallocate(RO_g)
      if(allocated(RO_go)) deallocate(RO_go)
      if(allocated(ROP_g)) deallocate(ROP_g)
      if(allocated(ROP_go)) deallocate(ROP_go)
      if(allocated(RO_S)) deallocate(RO_S)
      if(allocated(RO_So)) deallocate(RO_So)
      if(allocated(ROP_s)) deallocate(ROP_s)
      if(allocated(ROP_so)) deallocate(ROP_so)

      if(allocated(EP_SS)) deallocate(EP_SS)
      if(allocated(ERR_ARRAY)) deallocate(ERR_ARRAY)

      if(allocated(T_g)) deallocate(T_g)
      if(allocated(T_s)) deallocate(T_s)
      if(allocated(T_go)) deallocate(T_go)
      if(allocated(T_so)) deallocate(T_so)
      if(allocated(X_g)) deallocate(X_g)
      if(allocated(X_s)) deallocate(X_s)
      if(allocated(X_go)) deallocate(X_go)
      if(allocated(X_so)) deallocate(X_so)
      if(allocated(U_g)) deallocate(U_g)
      if(allocated(U_go)) deallocate(U_go)
      if(allocated(U_s)) deallocate(U_s)
      if(allocated(U_so)) deallocate(U_so)
      if(allocated(V_g)) deallocate(V_g)
      if(allocated(V_go)) deallocate(V_go)
      if(allocated(V_s)) deallocate(V_s)
      if(allocated(V_so)) deallocate(V_so)
      if(allocated(W_g)) deallocate(W_g)
      if(allocated(W_go)) deallocate(W_go)
      if(allocated(W_s)) deallocate(W_s)
      if(allocated(W_so)) deallocate(W_so)
      if(allocated(P_s)) deallocate(P_s)
      if(allocated(P_s_c)) deallocate(P_s_c)
      if(allocated(P_s_v)) deallocate(P_s_v)
      if(allocated(P_s_f)) deallocate(P_s_f)
      if(allocated(P_s_p)) deallocate(P_s_p)
      if(allocated(P_star)) deallocate(P_star)
      if(allocated(P_staro)) deallocate(P_staro)
      if(allocated(THETA_m)) deallocate(THETA_m)
      if(allocated(THETA_mo)) deallocate(THETA_mo)


      IF(K_Epsilon)THEN
        if(allocated(K_Turb_G)) deallocate(K_Turb_G)
        if(allocated(K_Turb_Go)) deallocate(K_Turb_Go)
        if(allocated(E_Turb_G)) deallocate(E_Turb_G)
        if(allocated(E_Turb_Go)) deallocate(E_Turb_Go)
      ENDIF

      if(allocated(Scalar)) deallocate(Scalar)
      if(allocated(Scalaro)) deallocate(Scalaro)



!pgcor
      if(allocated(d_e)) deallocate(d_e)
      if(allocated(d_n)) deallocate(d_n)
      if(allocated(d_t)) deallocate(d_t)
      if(allocated(Pp_g)) deallocate(Pp_g)
      if(allocated(PHASE_4_P_g)) deallocate(PHASE_4_P_g)

!physprop
      if(allocated(MU_g)) deallocate(MU_g)
      if(allocated(C_pg)) deallocate(C_pg)
      if(allocated(C_ps)) deallocate(C_ps)
      if(allocated(K_g)) deallocate(K_g)
      if(allocated(K_s)) deallocate(K_s)
      if(allocated(Kth_s)) deallocate(Kth_s)
      if(allocated(Kphi_s)) deallocate(Kphi_s)
      if(allocated(DIF_g)) deallocate(DIF_g)
      if(allocated(DIF_s)) deallocate(DIF_s)
      if(allocated(MW_MIX_g)) deallocate(MW_MIX_g)

!pscor
      if(allocated(e_e)) deallocate(e_e)
      if(allocated(e_n)) deallocate(e_n)
      if(allocated(e_t)) deallocate(e_t)
      if(allocated(K_cp)) deallocate(K_cp)
      if(allocated(EPp)) deallocate(EPp)
      if(allocated(PHASE_4_P_s)) deallocate(PHASE_4_P_s)


!residual
      if(allocated(RESID)) deallocate(RESID)
      if(allocated(MAX_RESID)) deallocate(MAX_RESID)
      if(allocated(IJK_RESID)) deallocate(IJK_RESID)
      if(allocated(NUM_RESID)) deallocate(NUM_RESID)
      if(allocated(DEN_RESID)) deallocate(DEN_RESID)
      if(allocated(RESID_PACK)) deallocate(RESID_PACK)

!rxns
      if(allocated(R_gp)) deallocate(R_gp)
      if(allocated(R_sp)) deallocate(R_sp)
      if(allocated(RoX_gc)) deallocate(RoX_gc)
      if(allocated(RoX_sc)) deallocate(RoX_sc)
      if(allocated(SUM_R_g)) deallocate(SUM_R_g)
      if(allocated(SUM_R_s)) deallocate(SUM_R_s)
      if(allocated(R_phase)) deallocate(R_phase)

!scalars
      if(allocated(Scalar_c)) deallocate(Scalar_c)
      if(allocated(Scalar_p)) deallocate(Scalar_p)
      if(allocated(Dif_Scalar)) deallocate(Dif_Scalar)

! add by rong for dqmom
      if(allocated(D_p)) deallocate(D_p)
      if(allocated(D_po)) deallocate(D_po)
      if(allocated(Source_a)) deallocate(Source_a)
      if(allocated(S_bar)) deallocate(S_bar)
      if(allocated(Matrix_a)) deallocate(Matrix_a)
      if(allocated(Matrix_b)) deallocate(Matrix_b)
      if(allocated(Matrix_c)) deallocate(Matrix_c)
      if(allocated(Inv_a)) deallocate(Inv_a)
      if(allocated(A)) deallocate(A)
      if(allocated(omega)) deallocate(omega)
      if(allocated(beta_a)) deallocate(beta_a)
      if(allocated(ystart)) deallocate(ystart)

! K-Epsilon Turbulence model
      if(allocated(K_Turb_G_c)) deallocate(K_Turb_G_c)
      if(allocated(K_Turb_G_p)) deallocate(K_Turb_G_p)
      if(allocated(Dif_K_Turb_G)) deallocate(Dif_K_Turb_G)
      if(allocated(E_Turb_G_c)) deallocate(E_Turb_G_c)
      if(allocated(E_Turb_G_p)) deallocate(E_Turb_G_p)
      if(allocated(Dif_E_Turb_G)) deallocate(Dif_E_Turb_G)

! Simonin or Ahmadi model
      if(allocated(K_12)) deallocate(K_12)
      if(allocated(Tau_12)) deallocate(Tau_12)
      if(allocated(Tau_1)) deallocate(Tau_1)


!tau_g
      if(allocated(TAU_U_g)) deallocate(TAU_U_g)
      if(allocated(TAU_V_g)) deallocate(TAU_V_g)
      if(allocated(TAU_W_g)) deallocate(TAU_W_g)
      if(allocated(DF_gu)) deallocate(DF_gu)
      if(allocated(DF_gv)) deallocate(DF_gv)
      if(allocated(DF_gw)) deallocate(DF_gw)
      if(allocated(CTAU_U_G)) deallocate(CTAU_U_G)
      if(allocated(CTAU_V_G)) deallocate(CTAU_V_G)
      if(allocated(CTAU_W_G)) deallocate(CTAU_W_G)

!tau_s
      if(allocated(TAU_U_s)) deallocate(TAU_U_s)
      if(allocated(TAU_V_s)) deallocate(TAU_V_s)
      if(allocated(TAU_W_s)) deallocate(TAU_W_s)


! generate_particles / particle_count
      if(allocated(PARTICLE_COUNT)) deallocate(PARTICLE_COUNT)


!trace
      if(allocated(trD_s_C)) deallocate(trD_s_C)
      if(allocated(trD_s2)) deallocate(trD_s2)
      if(allocated(trD_s_Co)) deallocate(trD_s_Co)
      if(allocated(trD_s_Co2)) deallocate(trD_s_Co2)

!visc_g
      if(allocated(trD_g)) deallocate(trD_g)
      if(allocated(MU_gt)) deallocate(MU_gt)
      if(allocated(EPMU_gt)) deallocate(EPMU_gt)
      if(allocated(LAMBDA_gt)) deallocate(LAMBDA_gt)
      if(allocated(EPLAMBDA_gt)) deallocate(EPLAMBDA_gt)
      if(allocated(L_scale)) deallocate(L_scale)


!visc_s
      if(allocated(MU_s)) deallocate(MU_s)
      if(allocated(EPMU_s)) deallocate(EPMU_s)
      if(allocated(LAMBDA_s)) deallocate(LAMBDA_s)
      if(allocated(EPLAMBDA_s)) deallocate(EPLAMBDA_s)
      if(allocated(ALPHA_s)) deallocate(ALPHA_s)
      if(allocated(MU_s_c)) deallocate(MU_s_c)
      if(allocated(LAMBDA_s_c)) deallocate(LAMBDA_s_c)
      if(allocated(LAMBDA_s_v)) deallocate(LAMBDA_s_v)
      if(allocated(LAMBDA_s_f)) deallocate(LAMBDA_s_f)
      if(allocated(LAMBDA_s_p)) deallocate(LAMBDA_s_p)
      if(allocated(MU_s_v)) deallocate(MU_s_v)
      if(allocated(MU_s_f)) deallocate(MU_s_f)
      if(allocated(MU_s_p)) deallocate(MU_s_p)
      if(allocated(MU_b_v)) deallocate(MU_b_v)
      if(allocated(EP_star_array)) deallocate(EP_star_array)
      if(allocated(EP_g_blend_start)) deallocate(EP_g_blend_start)
      if(allocated(EP_g_blend_end)) deallocate(EP_g_blend_end)
      if(allocated(trD_s)) deallocate(trD_s)
      if(allocated(I2_devD_s)) deallocate(I2_devD_s)
      if(allocated(TrM_s)) deallocate(TrM_s)
      if(allocated(TrDM_s)) deallocate(TrDM_s)


!shear quantities
      if(allocated(VSH)) deallocate(VSH)
      if(allocated(VSHE)) deallocate(VSHE)

!mflux
      if(allocated(Flux_gE)) deallocate(Flux_gE)
      if(allocated(Flux_sE)) deallocate(Flux_sE)
      if(allocated(Flux_gN)) deallocate(Flux_gN)
      if(allocated(Flux_sN)) deallocate(Flux_sN)
      if(allocated(Flux_gT)) deallocate(Flux_gT)
      if(allocated(Flux_sT)) deallocate(Flux_sT)

      if(allocated(Flux_gSE)) deallocate(Flux_gSE)
      if(allocated(Flux_sSE)) deallocate(Flux_sSE)
      if(allocated(Flux_gSN)) deallocate(Flux_gSN)
      if(allocated(Flux_sSN)) deallocate(Flux_sSN)
      if(allocated(Flux_gST)) deallocate(Flux_gST)
      if(allocated(Flux_sST)) deallocate(Flux_sST)

      if(allocated(ROP_gE)) deallocate(ROP_gE)
      if(allocated(ROP_sE)) deallocate(ROP_sE)
      if(allocated(ROP_gN)) deallocate(ROP_gN)
      if(allocated(ROP_sN)) deallocate(ROP_sN)
      if(allocated(ROP_gT)) deallocate(ROP_gT)
      if(allocated(ROP_sT)) deallocate(ROP_sT)

! allocate variables for GHD Theory
      if(allocated(Flux_nE)) deallocate(Flux_nE)
      if(allocated(Flux_nN)) deallocate(Flux_nN)
      if(allocated(Flux_nT)) deallocate(Flux_nT)
      if(allocated(Zeta0)) deallocate(Zeta0)
      if(allocated(ZetaU)) deallocate(ZetaU)
      if(allocated(DiT)) deallocate(DiT)
      if(allocated(DijF)) deallocate(DijF)
      if(allocated(Lij)) deallocate(Lij)
      if(allocated(Dij)) deallocate(Dij)
      if(allocated(DijQ)) deallocate(DijQ)
      if(allocated(JoiX)) deallocate(JoiX)
      if(allocated(JoiY)) deallocate(JoiY)
      if(allocated(JoiZ)) deallocate(JoiZ)
      if(allocated(FiX)) deallocate(FiX)
      if(allocated(FiY)) deallocate(FiY)
      if(allocated(FiZ)) deallocate(FiZ)
      if(allocated(FiXvel)) deallocate(FiXvel)
      if(allocated(FiYvel)) deallocate(FiYvel)
      if(allocated(FiZvel)) deallocate(FiZvel)
      if(allocated(DELTAU)) deallocate(DELTAU)
      if(allocated(DELTAV)) deallocate(DELTAV)
      if(allocated(DELTAW)) deallocate(DELTAW)
      if(allocated(dragFx)) deallocate(dragFx)
      if(allocated(dragFy)) deallocate(dragFy)
      if(allocated(dragFz)) deallocate(dragFz)
      if(allocated(dragFxflux)) deallocate(dragFxflux)
      if(allocated(dragFyflux)) deallocate(dragFyflux)
      if(allocated(dragFzflux)) deallocate(dragFzflux)
      if(allocated(FiMinusDragX)) deallocate(FiMinusDragX)
      if(allocated(JoiMinusDragX)) deallocate(JoiMinusDragX)
      if(allocated(FiMinusDragY)) deallocate(FiMinusDragY)
      if(allocated(JoiMinusDragY)) deallocate(JoiMinusDragY)
      if(allocated(FiMinusDragZ)) deallocate(FiMinusDragZ)
      if(allocated(JoiMinusDragZ)) deallocate(JoiMinusDragZ)
      if(allocated(beta_cell_X)) deallocate(beta_cell_X)
      if(allocated(beta_cell_Y)) deallocate(beta_cell_Y)
      if(allocated(beta_cell_Z)) deallocate(beta_cell_Z)
      if(allocated(beta_ij_cell_X)) deallocate(beta_ij_cell_X)
      if(allocated(beta_ij_cell_Y)) deallocate(beta_ij_cell_Y)
      if(allocated(beta_ij_cell_Z)) deallocate(beta_ij_cell_Z)
      if(allocated(DEL_DOT_J)) deallocate(DEL_DOT_J)
      if(allocated(DiT_HarmE)) deallocate(DiT_HarmE)
      if(allocated(DiT_HarmN)) deallocate(DiT_HarmN)
      if(allocated(DiT_HarmT)) deallocate(DiT_HarmT)
      if(allocated(Dij_HarmE)) deallocate(Dij_HarmE)
      if(allocated(Dij_HarmN)) deallocate(Dij_HarmN)
      if(allocated(Dij_HarmT)) deallocate(Dij_HarmT)
      if(allocated(DijF_HarmE)) deallocate(DijF_HarmE)
      if(allocated(DijF_HarmN)) deallocate(DijF_HarmN)
      if(allocated(DijF_HarmT)) deallocate(DijF_HarmT)


! We need to set this even when KT_TYPE is not set to IA_NONEP - at
! least in the current version of the code and needs to be revisited
      if(allocated(KTMOM_U_s)) deallocate(KTMOM_U_s)
      if(allocated(KTMOM_V_s)) deallocate(KTMOM_V_s)
      if(allocated(KTMOM_W_s)) deallocate(KTMOM_W_s)

! allocate variables for Iddir & Arastoopour (2005) kinetic theory
! EDvel_sM_ip & EDT_s_ip are also used for Garzy & Dufty (1999) kinetic theory
      if(allocated(trD_s2_ip)) deallocate(trD_s2_ip)
      if(allocated(MU_sM_ip)) deallocate(MU_sM_ip)
      if(allocated(MU_sL_ip)) deallocate(MU_sL_ip)
      if(allocated(XI_sM_ip)) deallocate(XI_sM_ip)
      if(allocated(XI_sL_ip)) deallocate(XI_sL_ip)
      if(allocated(Fnu_s_ip)) deallocate(Fnu_s_ip)
      if(allocated(FT_sM_ip)) deallocate(FT_sM_ip)
      if(allocated(FT_sL_ip)) deallocate(FT_sL_ip)
      if(allocated(Kth_sL_ip)) deallocate(Kth_sL_ip)
      if(allocated(Knu_sM_ip)) deallocate(Knu_sM_ip)
      if(allocated(Knu_sL_ip)) deallocate(Knu_sL_ip)
      if(allocated(Kvel_s_ip)) deallocate(Kvel_s_ip)
      if(allocated(EDvel_sL_ip)) deallocate(EDvel_sL_ip)
      if(allocated(ED_ss_ip)) deallocate(ED_ss_ip)

      if(allocated(A2_gtsh)) deallocate(A2_gtsh)
      if(allocated(xsi_gtsh)) deallocate(xsi_gtsh)

      if(allocated(EDT_s_ip)) deallocate(EDT_s_ip)
      if(allocated(EDvel_sM_ip)) deallocate(EDvel_sM_ip)

      if(allocated(errorpercent)) deallocate(errorpercent)

      RETURN
      END SUBROUTINE DEALLOCATE_ARRAYS_MAIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DEALLOCATE_ARRAYS_GEOMETRY                             !
!  Author: Jeff DIetiker                              Date: 16-MAR-2016!
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_ARRAYS_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Domain decomposition and dimensions
      use geometry, only: oDX, oDX_E
      use geometry, only: oDZ, oDZ_T
      use geometry, only: oDY, oDY_N
      use geometry, only: X, X_E, oX, oX_E, cyl_X, cyl_X_E
      use geometry, only: Z, Z_T
! Averaging factors.
      use geometry, only: FX_E, FX_E_bar, FX, FX_bar
      use geometry, only: FY_N, FY_N_bar
      use geometry, only: FZ_T, FZ_T_bar
! Domain flags.
      use geometry, only: FLAG
      use geometry, only: FLAG_E, FLAG_N, FLAG_T
! Domain volumes and areas.
      use geometry, only: VOL, VOL_SURR, AYZ, AXZ, AXY, CENTER_S! Scalar grid
      use geometry, only: VOL_U, AYZ_U, AXZ_U, AXY_U  ! X-Momentum
      use geometry, only: VOL_V, AYZ_V, AXZ_V, AXY_V  ! Y-Momentum
      use geometry, only: VOL_W, AYZ_W, AXZ_W, AXY_W  ! Z-Momentum

      use geometry, only: icbc_flag

      use gridmap, only : IJK_ARRAY_OF, FUNIJK_MAP_C, DEAD_CELL_AT

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Allocate geometry components related to the mesh. Check the
! allocation error status and abort if any failure is detected.
      if(allocated(X)) deallocate(X)
      if(allocated(cyl_X)) deallocate(cyl_X)
      if(allocated(X_E)) deallocate(X_E)
      if(allocated(cyl_X_E)) deallocate(cyl_X_E)
      if(allocated(oX)) deallocate(oX)
      if(allocated(oX_E)) deallocate(oX_E)
      if(allocated(oDX)) deallocate(oDX)
      if(allocated(oDX_E)) deallocate(oDX_E)


      if(allocated(oDY)) deallocate(oDY)
      if(allocated(oDY_N)) deallocate(oDY_N)


      if(allocated(Z)) deallocate(Z)
      if(allocated(Z_T)) deallocate(Z_T)
      if(allocated(oDZ)) deallocate(oDZ)
      if(allocated(oDZ_T)) deallocate(oDZ_T)


      if(allocated(FX)) deallocate(FX)
      if(allocated(FX_bar)) deallocate(FX_bar)


      if(allocated(FX_E)) deallocate(FX_E)
      if(allocated(FX_E_bar)) deallocate(FX_E_bar)


      if(allocated(FY_N)) deallocate(FY_N)
      if(allocated(FY_N_bar)) deallocate(FY_N_bar)

      if(allocated(FZ_T)) deallocate(FZ_T)
      if(allocated(FZ_T_bar)) deallocate(FZ_T_bar)


! Flags for the scalar grid.
      if(allocated(FLAG)) deallocate(FLAG)

! Flags for the momentum grids.
      if(allocated(FLAG_E)) deallocate(FLAG_E)
      if(allocated(FLAG_N)) deallocate(FLAG_N)
      if(allocated(FLAG_T)) deallocate(FLAG_T)

! Text flags for scalar grid.
     if(associated(ICBC_FLAG)) nullify(ICBC_FLAG)
!     deallocate(ICBC_FLAG)

     if(associated(ICBC_FLAG)) then
        deallocate(ICBC_FLAG)
        nullify(ICBC_FLAG)
     endif

!gridmap
     if(allocated(IJK_ARRAY_OF)) deallocate(IJK_ARRAY_OF)
     if(allocated(FUNIJK_MAP_C)) deallocate(FUNIJK_MAP_C)
     if(allocated(DEAD_CELL_AT)) deallocate(DEAD_CELL_AT)

! Volume and face-areas of scalar grid.
      if(allocated(VOL)) deallocate(VOL)
      if(allocated(AYZ)) deallocate(AYZ)
      if(allocated(AXZ)) deallocate(AXZ)
      if(allocated(AXY)) deallocate(AXY)

      if(allocated(CENTER_S)) deallocate(CENTER_S)

! total volume of each cell's surrounding stencil cells
      if(allocated(VOL_SURR)) deallocate(VOL_SURR)

! Volume and face-areas of X-Momentumn grid.
      if(allocated(VOL_U)) deallocate(VOL_U)
      if(allocated(AYZ_U)) deallocate(AYZ_U)
      if(allocated(AXZ_U)) deallocate(AXZ_U)
      if(allocated(AXY_U)) deallocate(AXY_U)

! Volume and face-areas of Y-Momentum grid.
      if(allocated(VOL_V)) deallocate(VOL_V)
      if(allocated(AYZ_V)) deallocate(AYZ_V)
      if(allocated(AXZ_V)) deallocate(AXZ_V)
      if(allocated(AXY_V)) deallocate(AXY_V)

! Volume and face-areas of Z-Momentum grid.
      if(allocated(VOL_W)) deallocate(VOL_W)
      if(allocated(AYZ_W)) deallocate(AYZ_W)
      if(allocated(AXZ_W)) deallocate(AXZ_W)
      if(allocated(AXY_W)) deallocate(AXY_W)


      RETURN
      END SUBROUTINE DEALLOCATE_ARRAYS_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DEALLOCATE_ARRAYS_PARALLEL                             !
!  Author: Jeff Dietiker                              Date: 16-MAR-2016!
!                                                                      !
!  Purpose: The purpose of this module is to create increments to be   !
!           stored in the array STORE_INCREMENT which will be added    !
!           to cell index ijk to find the effective indices of its     !
!           neighbors. These increments are found using the 'class'    !
!           of cell ijk. The class is determined based on the          !
!           neighboring cell type, i.e. wall or fluid.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_ARRAYS_PARALLEL

      USE param
      USE param1
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use desgrid
      use derived_types, only: dg_pic
      use desmpi
      use sendrecv

      IMPLICIT NONE

      INTEGER ::lijk
      INTEGER :: Dstatus

! gridmap_mod

      if(allocated(ISIZE_ALL)) deallocate(ISIZE_ALL)
      if(allocated(JSIZE_ALL)) deallocate(JSIZE_ALL)
      if(allocated(KSIZE_ALL)) deallocate(KSIZE_ALL)

      if(allocated(imap)) deallocate(imap)
      if(allocated(jmap)) deallocate(jmap)
      if(allocated(kmap)) deallocate(kmap)

      if(allocated(imap_c)) deallocate(imap_c)
      if(allocated(jmap_c)) deallocate(jmap_c)
      if(allocated(kmap_c)) deallocate(kmap_c)

! desgrid_mod
      if(allocated(dg_istart1_all)) deallocate(dg_istart1_all)
      if(allocated(dg_iend1_all)) deallocate(dg_iend1_all)
      if(allocated(dg_istart2_all)) deallocate(dg_istart2_all)
      if(allocated(dg_iend2_all)) deallocate(dg_iend2_all)
      if(allocated(dg_isize_all)) deallocate(dg_isize_all)

      if(allocated(dg_jstart1_all)) deallocate(dg_jstart1_all)
      if(allocated(dg_jend1_all)) deallocate(dg_jend1_all)
      if(allocated(dg_jstart2_all)) deallocate(dg_jstart2_all)
      if(allocated(dg_jend2_all)) deallocate(dg_jend2_all)
      if(allocated(dg_jsize_all)) deallocate(dg_jsize_all)

      if(allocated(dg_kstart1_all)) deallocate(dg_kstart1_all)
      if(allocated(dg_kend1_all)) deallocate(dg_kend1_all)
      if(allocated(dg_kstart2_all)) deallocate(dg_kstart2_all)
      if(allocated(dg_kend2_all)) deallocate(dg_kend2_all)
      if(allocated(dg_ksize_all)) deallocate(dg_ksize_all)

      if(allocated(dg_dx_all)) deallocate(dg_dx_all)
      if(allocated(dg_dy_all)) deallocate(dg_dy_all)
      if(allocated(dg_dz_all)) deallocate(dg_dz_all)

      if(allocated(dg_cycoffset)) deallocate(dg_cycoffset)
      if(allocated(icycoffset)) deallocate(icycoffset)

      if(allocated(dg_c1_all)) deallocate(dg_c1_all)
      if(allocated(dg_c2_all)) deallocate(dg_c2_all)
      if(allocated(dg_c3_all)) deallocate(dg_c3_all)

     if(allocated(dg_pic)) then

         do lijk = 1,size(dg_pic)
            if(associated(dg_pic(lijk)%p)) then
               deallocate(dg_pic(lijk)%p,STAT=Dstatus)
               nullify(dg_pic(lijk)%p)
            endif
         end do
         deallocate(dg_pic,STAT=Dstatus)
      endif


! des/mpi_init_des_mod.f
      if(allocated(dsendbuf)) then
         do lijk=1,size(dsendbuf)
            if(allocated(dsendbuf(lijk)%facebuf)) deallocate(dsendbuf(lijk)%facebuf)
         enddo
         deallocate(dsendbuf,STAT=Dstatus)
      endif

      if(allocated(drecvbuf)) then
         do lijk=1,size(drecvbuf)
            if(allocated(drecvbuf(lijk)%facebuf)) deallocate(drecvbuf(lijk)%facebuf)
         enddo
         deallocate(drecvbuf,STAT=Dstatus)
      endif

      if(allocated(isendindices)) deallocate(isendindices)
      if(allocated(irecvindices)) deallocate(irecvindices)
      if(allocated(isendreq)) deallocate(isendreq)
      if(allocated(irecvreq)) deallocate(irecvreq)
      if(allocated(isendcnt)) deallocate(isendcnt)
      if(allocated(dcycl_offset)) deallocate(dcycl_offset)
      if(allocated(ineighproc)) deallocate(ineighproc)
      if(allocated(iexchflag)) deallocate(iexchflag)
      if(allocated(iscattercnts)) deallocate(iscattercnts)
      if(allocated(igathercnts)) deallocate(igathercnts)
      if(allocated(idispls)) deallocate(idispls)


! JFD: commented line below because deallocate statement fails
      ! deallocate(sendproc,STAT=Dstatus)
       deallocate(xsend,STAT=Dstatus)
       deallocate(sendijk,STAT=Dstatus)

      ! deallocate(recvproc,STAT=Dstatus)
       deallocate(xrecv,STAT=Dstatus)
       deallocate(recvijk,STAT=Dstatus)

       deallocate(sendproc1,STAT=Dstatus)
       deallocate(xsend1,STAT=Dstatus)
       deallocate(sendijk1,STAT=Dstatus)

       deallocate(recvproc1,STAT=Dstatus)
       deallocate(xrecv1,STAT=Dstatus)
       deallocate(recvijk1,STAT=Dstatus)

       deallocate(sendproc2,STAT=Dstatus)
      !  deallocate(xsend2,STAT=Dstatus)
      !  deallocate(sendijk2,STAT=Dstatus)

       deallocate(recvproc2,STAT=Dstatus)
!       deallocate(xrecv2,STAT=Dstatus)
!       deallocate(recvijk2,STAT=Dstatus)

       deallocate(send_persistent_request,STAT=Dstatus)
       deallocate(recv_persistent_request,STAT=Dstatus)
       ! deallocate(send_persistent_request1,STAT=Dstatus)
       ! deallocate(recv_persistent_request1,STAT=Dstatus)
       ! deallocate(send_persistent_request2,STAT=Dstatus)
       ! deallocate(recv_persistent_request2,STAT=Dstatus)

       deallocate(dsendbuffer,STAT=Dstatus)
       deallocate(isendbuffer,STAT=Dstatus)
       deallocate(csendbuffer,STAT=Dstatus)

       deallocate(drecvbuffer,STAT=Dstatus)
       deallocate(irecvbuffer,STAT=Dstatus)
       deallocate(crecvbuffer,STAT=Dstatus)

      deallocate(sendrequest,STAT=Dstatus)
      deallocate(recvrequest,STAT=Dstatus)



      RETURN
      END SUBROUTINE DEALLOCATE_ARRAYS_PARALLEL

      SUBROUTINE DEALLOCATE_CUT_CELL_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: ALLOCATE_ARRAYS
!  Purpose: allocate arrays
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      Use indices

      USE cutcell
      USE stl
      USE discretelement

      IMPLICIT NONE

      if(allocated(INTERIOR_CELL_AT)) deallocate(INTERIOR_CELL_AT)

      if(allocated(XG_E)) deallocate(XG_E)
      if(allocated(YG_N)) deallocate(YG_N)
      if(allocated(ZG_T)) deallocate(ZG_T)

      if(allocated(X_U)) deallocate(X_U)
      if(allocated(Y_U)) deallocate(Y_U)
      if(allocated(Z_U)) deallocate(Z_U)

      if(allocated(X_V)) deallocate(X_V)
      if(allocated(Y_V)) deallocate(Y_V)
      if(allocated(Z_V)) deallocate(Z_V)

      if(allocated(X_W)) deallocate(X_W)
      if(allocated(Y_W)) deallocate(Y_W)
      if(allocated(Z_W)) deallocate(Z_W)

      if(allocated(INTERSECT_X)) deallocate(INTERSECT_X)
      if(allocated(INTERSECT_Y)) deallocate(INTERSECT_Y)
      if(allocated(INTERSECT_Z)) deallocate(INTERSECT_Z)

      if(allocated(X_int)) deallocate(X_int)
      if(allocated(Y_int)) deallocate(Y_int)
      if(allocated(Z_int)) deallocate(Z_int)

      if(allocated(X_NEW_POINT)) deallocate(X_NEW_POINT)
      if(allocated(Y_NEW_POINT)) deallocate(Y_NEW_POINT)
      if(allocated(Z_NEW_POINT)) deallocate(Z_NEW_POINT)

      if(allocated(X_NEW_U_POINT)) deallocate(X_NEW_U_POINT)
      if(allocated(Y_NEW_U_POINT)) deallocate(Y_NEW_U_POINT)
      if(allocated(Z_NEW_U_POINT)) deallocate(Z_NEW_U_POINT)

      if(allocated(X_NEW_V_POINT)) deallocate(X_NEW_V_POINT)
      if(allocated(Y_NEW_V_POINT)) deallocate(Y_NEW_V_POINT)
      if(allocated(Z_NEW_V_POINT)) deallocate(Z_NEW_V_POINT)

      if(allocated(X_NEW_W_POINT)) deallocate(X_NEW_W_POINT)
      if(allocated(Y_NEW_W_POINT)) deallocate(Y_NEW_W_POINT)
      if(allocated(Z_NEW_W_POINT)) deallocate(Z_NEW_W_POINT)

      if(allocated(NUMBER_OF_NODES)) deallocate(NUMBER_OF_NODES)
      if(allocated(NUMBER_OF_U_NODES)) deallocate(NUMBER_OF_U_NODES)
      if(allocated(NUMBER_OF_V_NODES)) deallocate(NUMBER_OF_V_NODES)
      if(allocated(NUMBER_OF_W_NODES)) deallocate(NUMBER_OF_W_NODES)

      if(allocated(CONNECTIVITY)) deallocate(CONNECTIVITY)
      if(allocated(CONNECTIVITY_U)) deallocate(CONNECTIVITY_U)
      if(allocated(CONNECTIVITY_V)) deallocate(CONNECTIVITY_V)
      if(allocated(CONNECTIVITY_W)) deallocate(CONNECTIVITY_W)

      if(allocated(PARTITION)) deallocate(PARTITION)

      if(allocated(WALL_U_AT)) deallocate(WALL_U_AT)
      if(allocated(WALL_V_AT)) deallocate(WALL_V_AT)
      if(allocated(WALL_W_AT)) deallocate(WALL_W_AT)

      if(allocated(Area_CUT)) deallocate(Area_CUT)
      if(allocated(Area_U_CUT)) deallocate(Area_U_CUT)
      if(allocated(Area_V_CUT)) deallocate(Area_V_CUT)
      if(allocated(Area_W_CUT)) deallocate(Area_W_CUT)


      if(allocated(DELX_Ue)) deallocate(DELX_Ue)
      if(allocated(DELX_Uw)) deallocate(DELX_Uw)
      if(allocated(DELY_Un)) deallocate(DELY_Un)
      if(allocated(DELY_Us)) deallocate(DELY_Us)
      if(allocated(DELZ_Ut)) deallocate(DELZ_Ut)
      if(allocated(DELZ_Ub)) deallocate(DELZ_Ub)

      if(allocated(DELX_Ve)) deallocate(DELX_Ve)
      if(allocated(DELX_Vw)) deallocate(DELX_Vw)
      if(allocated(DELY_Vn)) deallocate(DELY_Vn)
      if(allocated(DELY_Vs)) deallocate(DELY_Vs)
      if(allocated(DELZ_Vt)) deallocate(DELZ_Vt)
      if(allocated(DELZ_Vb)) deallocate(DELZ_Vb)

      if(allocated(DELX_We)) deallocate(DELX_We)
      if(allocated(DELX_Ww)) deallocate(DELX_Ww)
      if(allocated(DELY_Wn)) deallocate(DELY_Wn)
      if(allocated(DELY_Ws)) deallocate(DELY_Ws)
      if(allocated(DELZ_Wt)) deallocate(DELZ_Wt)
      if(allocated(DELZ_Wb)) deallocate(DELZ_Wb)

      if(allocated(X_U_ec)) deallocate(X_U_ec)
      if(allocated(Y_U_ec)) deallocate(Y_U_ec)
      if(allocated(Z_U_ec)) deallocate(Z_U_ec)
      if(allocated(X_U_nc)) deallocate(X_U_nc)
      if(allocated(Y_U_nc)) deallocate(Y_U_nc)
      if(allocated(Z_U_nc)) deallocate(Z_U_nc)
      if(allocated(X_U_tc)) deallocate(X_U_tc)
      if(allocated(Y_U_tc)) deallocate(Y_U_tc)
      if(allocated(Z_U_tc)) deallocate(Z_U_tc)

      if(allocated(X_V_ec)) deallocate(X_V_ec)
      if(allocated(Y_V_ec)) deallocate(Y_V_ec)
      if(allocated(Z_V_ec)) deallocate(Z_V_ec)
      if(allocated(X_V_nc)) deallocate(X_V_nc)
      if(allocated(Y_V_nc)) deallocate(Y_V_nc)
      if(allocated(Z_V_nc)) deallocate(Z_V_nc)
      if(allocated(X_V_tc)) deallocate(X_V_tc)
      if(allocated(Y_V_tc)) deallocate(Y_V_tc)
      if(allocated(Z_V_tc)) deallocate(Z_V_tc)

      if(allocated(X_W_ec)) deallocate(X_W_ec)
      if(allocated(Y_W_ec)) deallocate(Y_W_ec)
      if(allocated(Z_W_ec)) deallocate(Z_W_ec)
      if(allocated(X_W_nc)) deallocate(X_W_nc)
      if(allocated(Y_W_nc)) deallocate(Y_W_nc)
      if(allocated(Z_W_nc)) deallocate(Z_W_nc)
      if(allocated(X_W_tc)) deallocate(X_W_tc)
      if(allocated(Y_W_tc)) deallocate(Y_W_tc)
      if(allocated(Z_W_tc)) deallocate(Z_W_tc)

      if(allocated(DELH_Scalar)) deallocate(DELH_Scalar)
      DELH_Scalar(:) = ZERO

      if(allocated(ALPHA_c)) deallocate(ALPHA_c)
      ALPHA_C(:) = ZERO

      if(allocated(CC_AE)) deallocate(CC_AE)
      CC_AE(:) = ZERO

      if(allocated(CC_AN)) deallocate(CC_AN)
      CC_AN(:) = ZERO

      if(allocated(CC_AT)) deallocate(CC_AT)
      CC_AT(:) = ZERO

      if(allocated(DELH_U)) deallocate(DELH_U)
      if(allocated(Theta_Ue)) deallocate(Theta_Ue)
      if(allocated(Theta_Ue_bar)) deallocate(Theta_Ue_bar)
      if(allocated(Theta_U_ne)) deallocate(Theta_U_ne)
      if(allocated(Theta_U_nw)) deallocate(Theta_U_nw)
      if(allocated(Theta_U_te)) deallocate(Theta_U_te)
      if(allocated(Theta_U_tw)) deallocate(Theta_U_tw)
      if(allocated(ALPHA_Ue_c)) deallocate(ALPHA_Ue_c)
      if(allocated(NOC_U_E)) deallocate(NOC_U_E)
      if(allocated(Theta_Un)) deallocate(Theta_Un)
      if(allocated(Theta_Un_bar)) deallocate(Theta_Un_bar)
      if(allocated(ALPHA_Un_c)) deallocate(ALPHA_Un_c)
      if(allocated(NOC_U_N)) deallocate(NOC_U_N)
      if(allocated(Theta_Ut)) deallocate(Theta_Ut)
      if(allocated(Theta_Ut_bar)) deallocate(Theta_Ut_bar)
      if(allocated(ALPHA_Ut_c)) deallocate(ALPHA_Ut_c)
      if(allocated(NOC_U_T)) deallocate(NOC_U_T)
      if(allocated(A_UPG_E)) deallocate(A_UPG_E)
      if(allocated(A_UPG_W)) deallocate(A_UPG_W)

      if(allocated(DELH_V)) deallocate(DELH_V)
      if(allocated(Theta_V_ne)) deallocate(Theta_V_ne)
      if(allocated(Theta_V_se)) deallocate(Theta_V_se)
      if(allocated(Theta_Vn)) deallocate(Theta_Vn)
      if(allocated(Theta_Vn_bar)) deallocate(Theta_Vn_bar)
      if(allocated(Theta_V_nt)) deallocate(Theta_V_nt)
      if(allocated(Theta_V_st)) deallocate(Theta_V_st)
      if(allocated(Theta_Ve)) deallocate(Theta_Ve)
      if(allocated(Theta_Ve_bar)) deallocate(Theta_Ve_bar)
      if(allocated(ALPHA_Ve_c)) deallocate(ALPHA_Ve_c)
      if(allocated(NOC_V_E)) deallocate(NOC_V_E)
      if(allocated(ALPHA_Vn_c)) deallocate(ALPHA_Vn_c)
      if(allocated(NOC_V_N)) deallocate(NOC_V_N)
      if(allocated(Theta_Vt)) deallocate(Theta_Vt)
      if(allocated(Theta_Vt_bar)) deallocate(Theta_Vt_bar)
      if(allocated(ALPHA_Vt_c)) deallocate(ALPHA_Vt_c)
      if(allocated(NOC_V_T)) deallocate(NOC_V_T)
      if(allocated(A_VPG_N)) deallocate(A_VPG_N)
      if(allocated(A_VPG_S)) deallocate(A_VPG_S)

      if(allocated(DELH_W)) deallocate(DELH_W)
      if(allocated(Theta_W_te)) deallocate(Theta_W_te)
      if(allocated(Theta_W_be)) deallocate(Theta_W_be)
      if(allocated(Theta_W_tn)) deallocate(Theta_W_tn)
      if(allocated(Theta_W_bn)) deallocate(Theta_W_bn)
      if(allocated(Theta_Wt)) deallocate(Theta_Wt)
      if(allocated(Theta_Wt_bar)) deallocate(Theta_Wt_bar)
      if(allocated(Theta_We)) deallocate(Theta_We)
      if(allocated(Theta_We_bar)) deallocate(Theta_We_bar)
      if(allocated(ALPHA_We_c)) deallocate(ALPHA_We_c)
      if(allocated(NOC_W_E)) deallocate(NOC_W_E)
      if(allocated(Theta_Wn)) deallocate(Theta_Wn)
      if(allocated(Theta_Wn_bar)) deallocate(Theta_Wn_bar)
      if(allocated(ALPHA_Wn_c)) deallocate(ALPHA_Wn_c)
      if(allocated(NOC_W_N)) deallocate(NOC_W_N)
      if(allocated(ALPHA_Wt_c)) deallocate(ALPHA_Wt_c)
      if(allocated(NOC_W_T)) deallocate(NOC_W_T)
      if(allocated(A_WPG_T)) deallocate(A_WPG_T)
      if(allocated(A_WPG_B)) deallocate(A_WPG_B)

      if(allocated(NORMAL_S)) deallocate(NORMAL_S)
      if(allocated(NORMAL_U)) deallocate(NORMAL_U)
      if(allocated(NORMAL_V)) deallocate(NORMAL_V)
      if(allocated(NORMAL_W)) deallocate(NORMAL_W)

      if(allocated(REFP_S)) deallocate(REFP_S)
      if(allocated(REFP_U)) deallocate(REFP_U)
      if(allocated(REFP_V)) deallocate(REFP_V)
      if(allocated(REFP_W)) deallocate(REFP_W)

      if(allocated(ONEoDX_E_U)) deallocate(ONEoDX_E_U)
      if(allocated(ONEoDY_N_U)) deallocate(ONEoDY_N_U)
      if(allocated(ONEoDZ_T_U)) deallocate(ONEoDZ_T_U)

      if(allocated(ONEoDX_E_V)) deallocate(ONEoDX_E_V)
      if(allocated(ONEoDY_N_V)) deallocate(ONEoDY_N_V)
      if(allocated(ONEoDZ_T_V)) deallocate(ONEoDZ_T_V)

      if(allocated(ONEoDX_E_W)) deallocate(ONEoDX_E_W)
      if(allocated(ONEoDY_N_W)) deallocate(ONEoDY_N_W)
      if(allocated(ONEoDZ_T_W)) deallocate(ONEoDZ_T_W)

      if(allocated(Xn_int)) deallocate(Xn_int)
      if(allocated(Xn_U_int)) deallocate(Xn_U_int)
      if(allocated(Xn_V_int)) deallocate(Xn_V_int)
      if(allocated(Xn_W_int)) deallocate(Xn_W_int)

      if(allocated(Ye_int)) deallocate(Ye_int)
      if(allocated(Ye_U_int)) deallocate(Ye_U_int)
      if(allocated(Ye_V_int)) deallocate(Ye_V_int)
      if(allocated(Ye_W_int)) deallocate(Ye_W_int)

      if(allocated(Zt_int)) deallocate(Zt_int)
      if(allocated(Zt_U_int)) deallocate(Zt_U_int)
      if(allocated(Zt_V_int)) deallocate(Zt_V_int)
      if(allocated(Zt_W_int)) deallocate(Zt_W_int)

      if(allocated(SNAP)) deallocate(SNAP)
      if(allocated(SNAP_SCALAR)) deallocate(SNAP_SCALAR)

      if(allocated(CUT_TREATMENT_AT)) deallocate(CUT_TREATMENT_AT)
      if(allocated(CUT_U_TREATMENT_AT)) deallocate(CUT_U_TREATMENT_AT)
      if(allocated(CUT_V_TREATMENT_AT)) deallocate(CUT_V_TREATMENT_AT)
      if(allocated(CUT_W_TREATMENT_AT)) deallocate(CUT_W_TREATMENT_AT)

      if(allocated(CUT_CELL_AT)) deallocate(CUT_CELL_AT)
      if(allocated(CUT_U_CELL_AT)) deallocate(CUT_U_CELL_AT)
      if(allocated(CUT_V_CELL_AT)) deallocate(CUT_V_CELL_AT)
      if(allocated(CUT_W_CELL_AT)) deallocate(CUT_W_CELL_AT)

      if(allocated(SMALL_CELL_AT)) deallocate(SMALL_CELL_AT)

      if(allocated(SMALL_CELL_FLAG)) deallocate(SMALL_CELL_FLAG)

      if(allocated(BLOCKED_CELL_AT)) deallocate(BLOCKED_CELL_AT)
      if(allocated(BLOCKED_U_CELL_AT)) deallocate(BLOCKED_U_CELL_AT)
      if(allocated(BLOCKED_V_CELL_AT)) deallocate(BLOCKED_V_CELL_AT)
      if(allocated(BLOCKED_W_CELL_AT)) deallocate(BLOCKED_W_CELL_AT)

      if(allocated(STANDARD_CELL_AT)) deallocate(STANDARD_CELL_AT)
      if(allocated(STANDARD_U_CELL_AT)) deallocate(STANDARD_U_CELL_AT)
      if(allocated(STANDARD_V_CELL_AT)) deallocate(STANDARD_V_CELL_AT)
      if(allocated(STANDARD_W_CELL_AT)) deallocate(STANDARD_W_CELL_AT)


      if(allocated(VORTICITY)) deallocate(VORTICITY)
      if(allocated(LAMBDA2)) deallocate(LAMBDA2)

      if(allocated(TRD_G_OUT)) deallocate(TRD_G_OUT)
      if(allocated(PP_G_OUT)) deallocate(PP_G_OUT)
      if(allocated(EPP_OUT)) deallocate(EPP_OUT)

      if(allocated(dudx_OUT)) deallocate(dudx_OUT)
      if(allocated(dvdy_OUT)) deallocate(dvdy_OUT)
      if(allocated(delv_OUT)) deallocate(delv_OUT)

      if(allocated(U_MASTER_OF)) deallocate(U_MASTER_OF)
      if(allocated(V_MASTER_OF)) deallocate(V_MASTER_OF)
      if(allocated(W_MASTER_OF)) deallocate(W_MASTER_OF)

      if(allocated(BC_ID)) deallocate(BC_ID)
      if(allocated(BC_U_ID)) deallocate(BC_U_ID)
      if(allocated(BC_V_ID)) deallocate(BC_V_ID)
      if(allocated(BC_W_ID)) deallocate(BC_W_ID)

      if(allocated(DEBUG_CG)) deallocate(DEBUG_CG)

      if(allocated(U_g_CC)) deallocate(U_g_CC)
      if(allocated(V_g_CC)) deallocate(V_g_CC)
      if(allocated(W_g_CC)) deallocate(W_g_CC)

      if(allocated(U_s_CC)) deallocate(U_s_CC)
      if(allocated(V_s_CC)) deallocate(V_s_CC)
      if(allocated(W_s_CC)) deallocate(W_s_CC)

      if(allocated(N_FACET_AT)) Deallocate(N_FACET_AT)


      if(allocated(LIST_FACET_AT)) Deallocate(LIST_FACET_AT)

      if(allocated(POTENTIAL_CUT_CELL_AT)) Deallocate(POTENTIAL_CUT_CELL_AT)


      if(allocated(F_AT)) deallocate(F_AT)
      if(allocated(F_AT_SCALAR)) deallocate(F_AT_SCALAR)

      if(allocated(DWALL)) deallocate(DWALL)

      IF(allocated(SCALAR_NODE_XYZ)) deallocate(SCALAR_NODE_XYZ)
      IF(allocated(Ovol_around_node)) deallocate(Ovol_around_node)
      IF(allocated(SCALAR_NODE_ATWALL)) deallocate(SCALAR_NODE_ATWALL)

      IF(allocated(MESH_MASK)) deallocate(MESH_MASK)
      RETURN
      END SUBROUTINE DEALLOCATE_CUT_CELL_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEALLOCATE_VTU_ARRAYS                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DEALLOCATE_VTU_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE vtk
      IMPLICIT NONE
!-----------------------------------------------


      if(allocated(CLEANED_CONNECTIVITY)) Deallocate(CLEANED_CONNECTIVITY)
      if(allocated(COORDS_OF_POINTS)) Deallocate(COORDS_OF_POINTS)
      if(allocated(GLOBAL_I_OF)) Deallocate(GLOBAL_I_OF)
      if(allocated(GLOBAL_J_OF)) Deallocate(GLOBAL_J_OF)
      if(allocated(GLOBAL_K_OF)) Deallocate(GLOBAL_K_OF)
      if(allocated(GLOBAL_CONNECTIVITY)) Deallocate(GLOBAL_CONNECTIVITY)
      if(allocated(GLOBAL_CLEANED_CONNECTIVITY)) Deallocate(GLOBAL_CLEANED_CONNECTIVITY)
      if(allocated(GLOBAL_NUMBER_OF_NODES)) Deallocate(GLOBAL_NUMBER_OF_NODES)
      if(allocated(GLOBAL_COORDS_OF_POINTS)) Deallocate(GLOBAL_COORDS_OF_POINTS)
      if(allocated(GLOBAL_INTERIOR_CELL_AT)) Deallocate(GLOBAL_INTERIOR_CELL_AT)
      if(allocated(GLOBAL_BLOCKED_CELL_AT)) Deallocate(GLOBAL_BLOCKED_CELL_AT)
      if(allocated(GLOBAL_STANDARD_CELL_AT)) Deallocate(GLOBAL_STANDARD_CELL_AT)
      if(allocated(GLOBAL_CUT_CELL_AT)) Deallocate(GLOBAL_CUT_CELL_AT)
      if(allocated(GLOBAL_SMALL_CELL_AT)) Deallocate(GLOBAL_SMALL_CELL_AT)
      if(allocated(GLOBAL_SNAP)) Deallocate(GLOBAL_SNAP)
      if(allocated(GLOBAL_F_AT)) Deallocate(GLOBAL_F_AT)
      if(allocated(GLOBAL_BC_ID)) Deallocate(GLOBAL_BC_ID)
      if(allocated(GLOBAL_FLAG)) Deallocate(GLOBAL_FLAG)
      if(allocated(GLOBAL_X_NEW_POINT)) Deallocate(GLOBAL_X_NEW_POINT)
      if(allocated(GLOBAL_Y_NEW_POINT)) Deallocate(GLOBAL_Y_NEW_POINT)
      if(allocated(GLOBAL_Z_NEW_POINT)) Deallocate(GLOBAL_Z_NEW_POINT)

      if(allocated(BELONGS_TO_VTK_SUBDOMAIN)) Deallocate(BELONGS_TO_VTK_SUBDOMAIN)



      RETURN
      END SUBROUTINE DEALLOCATE_VTU_ARRAYS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEALLOCATE_DEM_MIO                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DEALLOCATE_DEM_MI

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: undefined
      USE des_bc, only: pi_factor, pi_count
      use des_bc, only: numfrac_limit
      use des_bc, only: dem_mi_time, dem_bc_poly_layout
      use des_bc, only: dem_mi
      use des_bc, only: dem_bcmi_ijkstart, dem_bcmi_ijkend
      IMPLICIT NONE
!-----------------------------------------------

! Particle injection factor
      if(allocated(PI_FACTOR)) deallocate(PI_FACTOR)
! Particle injection count (injection number)
      if(allocated(PI_COUNT)) deallocate(PI_COUNT)
! Particle injection time scale
      if(allocated(DEM_MI_TIME)) deallocate(DEM_MI_TIME)
! Array used for polydisperse inlets: stores the particle number
! distribution of an inlet scaled with numfrac_limit
      if(allocated(DEM_BC_POLY_LAYOUT)) deallocate(DEM_BC_POLY_LAYOUT)
! Data structure for storing BC data.
      if(allocated(DEM_MI)) deallocate(DEM_MI)

      if(allocated(DEM_BCMI_IJKSTART)) deallocate(DEM_BCMI_IJKSTART)
      if(allocated(DEM_BCMI_IJKEND)) deallocate(DEM_BCMI_IJKEND)




      RETURN
      END SUBROUTINE DEALLOCATE_DEM_MI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEALLOCATE_PIC_MI                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_PIC_MI

! Modules
!-----------------------------------------------
      use mfix_pic, only: pic_mi

      implicit none
!-----------------------------------------------
      integer :: lc1
! Allocate/Initialize for inlets

      if(allocated(pic_mi)) then
         do lc1=lbound(pic_mi,1),ubound(pic_mi,1)
            if(allocated(pic_mi(lc1)%blocked)) &
                 deallocate(pic_mi(lc1)%blocked)
         enddo
         deallocate(pic_mi)
      endif

      RETURN
      END SUBROUTINE DEALLOCATE_PIC_MI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_DEALLOCATE_ARRAYS                                   C
!  Purpose: Deallocate arrays subroutines for DES                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DES_DEALLOCATE_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE compar
      USE constant
      USE cutcell
      USE derived_types, only: pic
      USE des_bc, only: dem_bcmo_ijk
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE functions
      USE funits
      USE geometry
      USE indices
      USE mfix_pic, only: des_stat_wt, des_stat_wt, ps_grad, epg_p
      USE param
      USE param1
      USE physprop
      USE mfix_pic, only: pic_bcmi

      USE particle_filter, only: DES_INTERP_GARG
      USE particle_filter, only: DES_INTERP_DPVM
      USE particle_filter, only: DES_INTERP_GAUSS
      USE particle_filter, only: DES_INTERP_LHAT
      USE particle_filter, only: FILTER_CELL
      USE particle_filter, only: FILTER_WEIGHT
      use des_bc, only: DEM_BCMO_IJKSTART, DEM_BCMO_IJKEND
      use stl, only: FACETS_AT_DG
      use vtk, only: Part_RRates_out

      use desgrid
      use sendrecvnode

      IMPLICIT NONE
      INTEGER:: lijk
      INTEGER:: Dstatus

! DES Allocatable arrays
!-----------------------------------------------
! Dynamic particle info including another index for parallel
! processing for ghost
      if(allocated(PARTICLE_STATE)) deallocate(PARTICLE_STATE)
      if(allocated(iglobal_id)) deallocate(iglobal_id)

! R.Garg: Allocate necessary arrays for PIC mass inlet/outlet BCs
      if(PIC_BCMI /= 0) CALL DEALLOCATE_PIC_MI

! Particle attributes
! Radius, density, mass, moment of inertia
      if(allocated(DES_RADIUS)) deallocate(DES_RADIUS)
      if(allocated(RO_Sol)) deallocate(RO_Sol)
      if(allocated(PVOL)) deallocate(PVOL)
      if(allocated(PMASS)) deallocate(PMASS)
      if(allocated(OMOI)) deallocate(OMOI)
! SuperDEM-Moment of inertia, semi-axi,roundness parameters,quaternions
      if(allocated(OMOI3)) deallocate(OMOI3)
      if(allocated(super_r)) deallocate(super_r)
      if(allocated(super_mn)) deallocate(super_mn)
      if(allocated(super_q)) deallocate(super_q)
! Glued particles
      if(GluedSphereDEM) then
         if(allocated(particle_state_gsp)) deallocate(particle_state_gsp)
         if(allocated(sphereNumber)) deallocate(sphereNumber)
         if(allocated(sphereRelPos)) deallocate(sphereRelPos)
         if(allocated(sphereRelDia)) deallocate(sphereRelDia)
         if(allocated(tmp_gp_sa)) deallocate(tmp_gp_sa)
         if(allocated(tmp_gp_neighsid)) deallocate(tmp_gp_neighsid)
         if(allocated(tmp_gp_neighsa)) deallocate(tmp_gp_neighsa)

         if(allocated(glueMass)) deallocate(glueMass)
	   if(allocated(glueVolume)) deallocate(glueVolume)
         if(allocated(glueDiameter)) deallocate(glueDiameter)
         if(allocated(glueBounding)) deallocate(glueBounding)
         if(allocated(glueForce)) deallocate(glueForce)
         if(allocated(glueTorque)) deallocate(glueTorque)
         !if(allocated(glueForcePP)) deallocate(glueForcePP)
         !if(allocated(glueTorquePP)) deallocate(glueTorquePP)
         if(allocated(glueAcc)) deallocate(glueAcc)
         if(allocated(glueAngMom)) deallocate(glueAngMom)
         if(allocated(glueQuat)) deallocate(glueQuat)
         if(allocated(glueInertia)) deallocate(glueInertia)
         if(allocated(glueOmg)) deallocate(glueOmg)
         if(allocated(glueVel)) deallocate(glueVel)
         if(allocated(gluePos)) deallocate(gluePos)
         !if(allocated(glueContactPoint)) deallocate(glueContactPoint)
         if(allocated(glueEX)) deallocate(glueEX)
         if(allocated(glueEY)) deallocate(glueEY)
         if(allocated(glueEZ)) deallocate(glueEZ)

         if(allocated(gid_list)) deallocate(gid_list)
         if(allocated(gp_component_IC)) deallocate(gp_component_IC)
         if(allocated(gp_neighsid)) deallocate(gp_neighsid)
         if(allocated(gp_neighsa)) deallocate(gp_neighsa)
         if(allocated(gp_kps)) deallocate(gp_kps)
         if(allocated(gp_sa)) deallocate(gp_sa)
         !if(allocated(gp_sdrag_mass)) deallocate(gp_sdrag_mass)
         if(allocated(gp_squat)) deallocate(gp_squat)
         if(allocated(sc2gpc_vec)) deallocate(sc2gpc_vec)
         if(allocated(cglobal_id)) deallocate(cglobal_id)
         if(allocated(marked_gsp)) deallocate(marked_gsp)
         if(allocated(marked_state)) deallocate(marked_state)
! for gsp_implicit
         if(allocated(cspart_data)) deallocate(cspart_data)
         if(allocated(glueSurface)) deallocate(glueSurface)
! for gsp output function in gsp implicit
         if(allocated(sphere_radius_out)) deallocate(sphere_radius_out)
         if(allocated(sphere_vel_out)) deallocate(sphere_vel_out)
         if(allocated(sphere_rot_out)) deallocate(sphere_rot_out)
         if(allocated(sphere_orient_out)) deallocate(sphere_orient_out)
         if(allocated(sphere_usr_var_out)) deallocate(sphere_usr_var_out)
         if(allocated(sphere_temp_out)) deallocate(sphere_temp_out)
         if(allocated(sphere_cps_out)) deallocate(sphere_cps_out)
         if(allocated(sphere_xs_out)) deallocate(sphere_xs_out)
         if(allocated(sphere_prates_out)) deallocate(sphere_prates_out)
         if(allocated(sphere_rho_out)) deallocate(sphere_rho_out)
         if(allocated(sphere_mass_out)) deallocate(sphere_mass_out)
         if(allocated(sphere_vol_out)) deallocate(sphere_vol_out)
         if(allocated(sphere_rank_out)) deallocate(sphere_rank_out)
         if(allocated(sphere_pid_out)) deallocate(sphere_pid_out)
         if(allocated(sphere_id_out)) deallocate(sphere_id_out)
         if(allocated(sphere_res_time_out)) deallocate(sphere_res_time_out)
         if(allocated(sphere_coordination_out)) deallocate(sphere_coordination_out)

         if(allocated(sphere_MIN_OVERLAP_out)) deallocate(sphere_MIN_OVERLAP_out)
         if(allocated(sphere_MAX_OVERLAP_out)) deallocate(sphere_MAX_OVERLAP_out)
         if(allocated(sphere_MEAN_OVERLAP_out)) deallocate(sphere_MEAN_OVERLAP_out)

         if(allocated(sphere_coll_force_out)) deallocate(sphere_coll_force_out)
         if(allocated(sphere_drag_out)) deallocate(sphere_drag_out)

         if(allocated(bftx)) deallocate(bftx)
         if(allocated(bfty)) deallocate(bfty)
         if(allocated(bftz)) deallocate(bftz)

         if(allocated(bmtx)) deallocate(bmtx)
         if(allocated(bmty)) deallocate(bmty)
         if(allocated(bmtz)) deallocate(bmtz)

         if(allocated(bfn)) deallocate(bfn)
         if(allocated(bmn)) deallocate(bmn)

         if(allocated(bnx)) deallocate(bnx)
         if(allocated(bny)) deallocate(bny)
         if(allocated(bnz)) deallocate(bnz)

      !    if(allocated(PFT_NEIGHBOR_X)) deallocate(PFT_NEIGHBOR_X)
      !    if(allocated(PFT_NEIGHBOR_Y)) deallocate(PFT_NEIGHBOR_Y)
      !    if(allocated(PFT_NEIGHBOR_Z)) deallocate(PFT_NEIGHBOR_Z)
         if(allocated(global2loc_id)) deallocate(global2loc_id)
      endif
! Old and new particle positions, velocities (translational and
! rotational)
      if(allocated(DES_POS_NEW)) deallocate(DES_POS_NEW)
      if(allocated(DES_VEL_NEW)) deallocate(DES_VEL_NEW)
      if(allocated(OMEGA_NEW)) deallocate(OMEGA_NEW)

      if(allocated(ORIENTATION)) deallocate(ORIENTATION)

      if(allocated(DES_POS_OLD)) deallocate(DES_POS_OLD)
      if(allocated(DES_VEL_OLD)) deallocate(DES_VEL_OLD)
      if(allocated(DES_ACC_OLD)) deallocate(DES_ACC_OLD)
      if(allocated(OMEGA_OLD)) deallocate(OMEGA_OLD)
      if(allocated(ROT_ACC_OLD)) deallocate(ROT_ACC_OLD)

! Force chain data
      if(allocated(FCHAIN_MIDPOINT)) deallocate(FCHAIN_MIDPOINT)
      if(allocated(FCHAIN_ORIENT)) deallocate(FCHAIN_ORIENT)
      if(allocated(FCHAIN_FN)) deallocate(FCHAIN_FN)
      if(allocated(FCHAIN_FT)) deallocate(FCHAIN_FT)
      if(allocated(FCHAIN_LENGTH)) deallocate(FCHAIN_LENGTH)
      if(allocated(FCHAIN_FN_MAG)) deallocate(FCHAIN_FN_MAG)
      if(allocated(FCHAIN_FCMAX)) deallocate(FCHAIN_FCMAX)
      if(allocated(FCHAIN_LOCAL_ID1)) deallocate(FCHAIN_LOCAL_ID1)
      if(allocated(FCHAIN_LOCAL_ID2)) deallocate(FCHAIN_LOCAL_ID2)
      if(allocated(FCHAIN_GLOBAL_ID1)) deallocate(FCHAIN_GLOBAL_ID1)
      if(allocated(FCHAIN_GLOBAL_ID2)) deallocate(FCHAIN_GLOBAL_ID2)
      if(allocated(FCHAIN_OVERLAP)) deallocate(FCHAIN_OVERLAP)

      if(allocated(COORDINATION)) deallocate(COORDINATION)
      if(allocated(MIN_OVERLAP)) deallocate(MIN_OVERLAP)
      if(allocated(MAX_OVERLAP)) deallocate(MAX_OVERLAP)
      if(allocated(MEAN_OVERLAP)) deallocate(MEAN_OVERLAP)

! Collision force
      if(allocated(DES_COLL_FORCE)) deallocate(DES_COLL_FORCE)

! Residence time
      if(allocated(RESIDENCE_TIME)) deallocate(RESIDENCE_TIME)

! Allocating user defined array
      if(allocated(DES_USR_VAR)) deallocate(DES_USR_VAR)

! Particle positions at the last call neighbor search algorithm call
      if(allocated(PPOS)) deallocate(PPOS)

! Total, normal and tangential forces
      if(allocated(FC)) deallocate(FC)

! Torque
      if(allocated(TOW)) deallocate(TOW)


! allocate variable for des grid binning
      if(allocated(dg_pijk)) deallocate(dg_pijk)
      if(allocated(dg_pijkprv)) deallocate(dg_pijkprv)

! allocate variables related to ghost particles
      if(allocated(ighost_updated)) deallocate(ighost_updated)


      if(allocated(wall_collision_facet_id)) deallocate(wall_collision_facet_id)

      if(allocated(wall_collision_PFT)) deallocate(wall_collision_PFT)

! Temporary variables to store wall position, velocity and normal vector
      if(allocated(WALL_NORMAL)) deallocate(WALL_NORMAL)

      if(allocated(NEIGHBOR_INDEX)) deallocate(NEIGHBOR_INDEX)
      if(allocated(NEIGHBOR_INDEX_OLD)) deallocate(NEIGHBOR_INDEX_OLD)
      if(allocated(NEIGHBORS)) deallocate(NEIGHBORS)

      if(allocated(NEIGHBORS_OLD)) deallocate(NEIGHBORS_OLD)
      if(allocated(PFT_NEIGHBOR)) deallocate(PFT_NEIGHBOR)
      if(allocated(PFT_NEIGHBOR_OLD)) deallocate(PFT_NEIGHBOR_OLD)

! SuperDEM, contact points on a pair of particles
      if(SuperDEM)  then
         if(allocated(CONTACT_POINT_A)) deallocate(CONTACT_POINT_A)
         if(allocated(CONTACT_POINT_A_OLD)) deallocate(CONTACT_POINT_A_OLD)
         if(allocated(CONTACT_POINT_B)) deallocate(CONTACT_POINT_B)
         if(allocated(CONTACT_POINT_B_OLD)) deallocate(CONTACT_POINT_B_OLD)

         if(allocated(CONTACT_lambda_A)) deallocate(CONTACT_LAMBDA_A)
         if(allocated(CONTACT_lambda_A_OLD)) deallocate(CONTACT_LAMBDA_A_OLD)
         if(allocated(CONTACT_lambda_B)) deallocate(CONTACT_LAMBDA_B)
         if(allocated(CONTACT_lambda_B_OLD)) deallocate(CONTACT_LAMBDA_B_OLD)
      endif


! Variable that stores the particle in cell information (ID) on the
! computational fluid grid defined by keywords IMAX, JMAX and KMAX.
      if(allocated(PIC)) then
         do lijk = 1,size(pic)
            if(associated(pic(lijk)%p)) then
               deallocate(pic(lijk)%p,STAT=Dstatus)
               nullify(pic(lijk)%p)
            endif
         enddo
         deallocate(PIC)
      endif

! Particles in a computational fluid cell (for volume fraction)
      if(allocated(PINC)) deallocate(PINC)

! Ghost Particles in a computational fluid cell (for volume fraction)
      if(allocated(GPINC)) deallocate(GPINC)

! For each particle track its i,j,k location on computational fluid grid
! defined by keywords IMAX, JMAX and KMAX.
      if(allocated(PIJK)) deallocate(PIJK)

      if(allocated(DRAG_AM)) deallocate(DRAG_AM)
      if(allocated(DRAG_BM)) deallocate(DRAG_BM)
      if(allocated(F_gp)) deallocate(F_gp)


! Explicit drag force acting on a particle.
      if(allocated(DRAG_FC)) deallocate(DRAG_FC)

! force due to gas-pressure gradient
      if(allocated(P_FORCE)) deallocate(P_FORCE)

! Volume of nodes
      if(allocated(DES_VOL_NODE)) deallocate(DES_VOL_NODE)

      if(allocated(F_GDS)) deallocate(F_GDS)
      if(allocated(VXF_GDS)) deallocate(VXF_GDS)

      if(allocated(FILTER_CELL)) deallocate(FILTER_CELL)
      if(allocated(FILTER_WEIGHT)) deallocate(FILTER_WEIGHT)
      if(allocated(DES_ROPS_NODE)) deallocate(DES_ROPS_NODE)
      if(allocated(DES_VEL_NODE)) deallocate(DES_VEL_NODE)

! Variables for hybrid model

      if(allocated(SDRAG_AM)) deallocate(SDRAG_AM)
      if(allocated(SDRAG_BM)) deallocate(SDRAG_BM)

      if(allocated(F_SDS)) deallocate(F_SDS)
      if(allocated(VXF_SDS)) deallocate(VXF_SDS)


! MP-PIC related
      if(allocated(DES_STAT_WT)) deallocate(DES_STAT_WT)
      if(allocated(DES_VEL_MAX)) deallocate(DES_VEL_MAX)
      if(allocated(PS_GRAD)) deallocate(PS_GRAD)
      if(allocated(EPG_P)) deallocate(EPG_P)

! Averaged velocity obtained by averaging over all the particles
      if(allocated(DES_VEL_AVG)) deallocate(DES_VEL_AVG)

! Global Granular Energy
      if(allocated(GLOBAL_GRAN_ENERGY)) deallocate(GLOBAL_GRAN_ENERGY)
      if(allocated(GLOBAL_GRAN_TEMP)) deallocate(GLOBAL_GRAN_TEMP)

! variable for bed height of solids phase M
      if(allocated(BED_HEIGHT)) deallocate(BED_HEIGHT)

! ---------------------------------------------------------------->>>
! BEGIN COHESION
! Matrix location of particle  (should be allocated in case user wishes
! to invoke routines in /cohesion subdirectory
      if(allocated(PostCohesive)) deallocate(PostCohesive)
! END COHESION
! ----------------------------------------------------------------<<<

! ---------------------------------------------------------------->>>
! BEGIN Thermodynamic Allocation
! Particle temperature
      if(allocated(DES_T_s)) deallocate(DES_T_s)
! Spec      ific heat
      if(allocated(DES_C_PS)) deallocate(DES_C_PS)
! Species mass fractions comprising a particle. This array may not be
! needed for all thermo problems.
      if(allocated(DES_X_s)) deallocate(DES_X_s)
! Minimum particle mass of each solid phase
      if(allocated(DES_MIN_PMASS)) deallocate(DES_MIN_PMASS)
! Total rate of heat transfer to individual particles.
      if(allocated(Q_Source)) deallocate(Q_Source)
! Average solids temperature in fluid cell
      if(allocated(avgDES_T_s)) deallocate(avgDES_T_s)
! Gas/Solids convective heat transfer coupling

! Fluid phase energy equation source terms
      if(allocated(CONV_Sc)) deallocate(CONV_Sc)
      if(allocated(CONV_Sp)) deallocate(CONV_Sp)
! Particle convection source term (explicit coupled)
      if(allocated(CONV_Qs)) deallocate(CONV_Qs)
! Gas-particle heat transfer coefficient TIMES surface area
      if(allocated(GAMMAxSA)) deallocate(GAMMAxSA)

! Allocate the history variables for Adams-Bashforth integration
      if(allocated(Q_Source0)) deallocate(Q_Source0)

! End Thermodynamic Allocation
! ----------------------------------------------------------------<<<


! ---------------------------------------------------------------->>>
! BEGIN Species Allocation
! Rate of solids phase production/consumption for each species
      if(allocated(DES_R_s)) deallocate(DES_R_s)

      if(allocated(DES_R_gp)) deallocate(DES_R_gp)
      if(allocated(DES_R_gc)) deallocate(DES_R_gc)
      if(allocated(DES_SUM_R_g)) deallocate(DES_SUM_R_g)
      if(allocated(DES_R_PHASE)) deallocate(DES_R_PHASE)
      if(allocated(DES_HOR_g)) deallocate(DES_HOR_g)

! Particle reaction rate for vtp file output
      if(allocated(Part_RRates_out)) deallocate(Part_RRates_out)

! Allocate the history variables for Adams-Bashforth integration

! Rate of change of particle mass
      if(allocated(dMdt_OLD)) deallocate(dMdt_OLD)
! Rate of change of particle mass percent species
      if(allocated(dXdt_OLD)) deallocate(dXdt_OLD)


! Energy generation from reaction (cal/sec)
      if(allocated(RXNS_Qs)) deallocate(RXNS_Qs)

! End Species Allocation
! ----------------------------------------------------------------<<<

      if(allocated(PARTICLE_STATE)) deallocate(PARTICLE_STATE)
      if(allocated(iglobal_id)) deallocate(iglobal_id)

      if(allocated(DEM_BCMO_IJK)) deallocate(DEM_BCMO_IJK)
      if(allocated(DEM_BCMO_IJKSTART)) deallocate(DEM_BCMO_IJKSTART)
      if(allocated(DEM_BCMO_IJKEND)) deallocate(DEM_BCMO_IJKEND)

! stl
      if(allocated(FACETS_AT_DG)) THEN
         do lijk = 1,size(FACETS_AT_DG)
            if (allocated(FACETS_AT_DG(lijk)%ID)) deallocate(FACETS_AT_DG(lijk)%ID)
            if (allocated(FACETS_AT_DG(lijk)%DIR)) deallocate(FACETS_AT_DG(lijk)%DIR)
            if (allocated(FACETS_AT_DG(lijk)%MIN)) deallocate(FACETS_AT_DG(lijk)%MIN)
            if (allocated(FACETS_AT_DG(lijk)%MAX)) deallocate(FACETS_AT_DG(lijk)%MAX)
         enddo
         deallocate(FACETS_AT_DG)
      endif


      if(allocated(GSTENCIL)) deallocate(GSTENCIL)
      if(allocated(VSTENCIL)) deallocate(VSTENCIL)
      if(allocated(PGRADSTENCIL)) deallocate(PGRADSTENCIL)
      if(allocated(PSGRADSTENCIL)) deallocate(PSGRADSTENCIL)
      if(allocated(VEL_SOL_STENCIL)) deallocate(VEL_SOL_STENCIL)
      if(allocated(SSTENCIL)) deallocate(SSTENCIL)

      if(allocated(dg_xe)) deallocate(dg_xe)
      if(allocated(dg_yn)) deallocate(dg_yn)
      if(allocated(dg_zt)) deallocate(dg_zt)

      if(allocated(xe)) deallocate(xe)
      if(allocated(yn)) deallocate(yn)
      if(allocated(zt)) deallocate(zt)

       call deallocate_des_nodes_pointers()

! Coarse Grain DEM
      if(allocated(des_cgp_stw)) deallocate(des_cgp_stw)
      if(allocated(des_cgp_rpr)) deallocate(des_cgp_rpr)

      RETURN

   END SUBROUTINE DES_DEALLOCATE_ARRAYS

END MODULE ALLOCATE_ARRAYS_MOD
