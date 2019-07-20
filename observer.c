#include "ACMSim.h"

#if MACHINE_TYPE == INDUCTION_MACHINE
struct InductionMachine im;
struct Observer ob;
void acm_init(){
    int i;
    for(i=0; i<2; ++i){
        im.us[i] = 0;
        im.is[i] = 0;
        im.us_curr[i] = 0;
        im.is_curr[i] = 0;
        im.us_prev[i] = 0;
        im.is_prev[i] = 0;        
    }

    im.Js = ACM.Js;
    im.Js_inv = 1./im.Js;
    im.rs = ACM.rs;
    im.rreq = ACM.rreq;
    im.alpha = ACM.alpha;
    im.Lsigma = ACM.Lsigma;
    im.Lsigma_inv = 1/im.Lsigma;
    im.Lmu = ACM.Lmu;
    im.Lmu_inv = ACM.Lmu_inv;
    im.npp = ACM.npp;
    im.omg = 0;
    im.theta_r = 0.0;
}
void ob_init(){

    ob.Js = im.Js;
    ob.Js_inv = im.Js_inv;
    ob.rs = im.rs;
    ob.rreq = im.rreq;
    ob.alpha = im.alpha;
    ob.Lsigma = im.Lsigma;
    ob.Lsigma_inv = im.Lsigma_inv;
    ob.Lmu = im.Lmu;
    ob.Lmu_inv = im.Lmu_inv;
    ob.npp = im.npp;
    ob.omg = im.omg;
    im.theta_r = 0.0;
    im.theta_d = 0.0;

    ob.psi_mu_al = 0.0;
    ob.psi_mu_be = 0.0;
}
void observation(){
    double rotor_flux_cmd, iMs, iTs, uMs_cmd, uTs_cmd;

    rotor_flux_cmd = CTRL.rotor_flux_cmd;
    iMs = CTRL.iMs;
    iTs = CTRL.iTs;
    uMs_cmd  = CTRL.uMs_cmd;
    uTs_cmd  = CTRL.uTs_cmd;

    // Speed estimation: Tajima1996
        ob.tajima.K_PEM = 2; // 0.1 ~ 2
        ob.tajima.omega_syn  = (uTs_cmd - ob.rs*iTs) / (rotor_flux_cmd + ob.Lsigma*iMs);
        // int i;
        // for(i=0;i<100;++i)
        {
            ob.tajima.e_M        = uMs_cmd - ob.rs*iMs + ob.tajima.omega_syn*ob.Lsigma*iTs;
            ob.tajima.omega_syn -= ob.tajima.K_PEM * ob.tajima.e_M;
        }
        ob.tajima.omega_sl   = ob.rreq*iTs / rotor_flux_cmd;
        ob.tajima.omg = ob.tajima.omega_syn - ob.tajima.omega_sl; // Instantaneous Velocity Computation

    // Flux estimation 1: Voltage model (implemented by shitty Euler method for now)
        double deriv_psi_s_al;
        double deriv_psi_s_be;
        static double psi_s_al = 0.0;
        static double psi_s_be = 0.0;
        deriv_psi_s_al = US_C(0) - ob.rs*IS_C(0);
        deriv_psi_s_be = US_C(1) - ob.rs*IS_C(1);
        psi_s_al += TS*deriv_psi_s_al;
        psi_s_be += TS*deriv_psi_s_be;
        ob.psi_mu_al = psi_s_al - ob.Lsigma*IS(0);
        ob.psi_mu_be = psi_s_be - ob.Lsigma*IS(1);

    // Flux estimation 2: Current model (this is a bad flux estimator)    
}

#elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
struct SynchronousMachine sm;
struct Observer ob;
void acm_init(){
    int i;
    for(i=0; i<2; ++i){
        sm.us[i] = 0;
        sm.is[i] = 0;
        sm.us_curr[i] = 0;
        sm.is_curr[i] = 0;
        sm.us_prev[i] = 0;
        sm.is_prev[i] = 0;        
    }

    sm.Js = ACM.Js;
    sm.Js_inv = 1./sm.Js;

    sm.R = ACM.R;
    sm.KE = ACM.KE;
    sm.Ld = ACM.Ld;
    sm.Lq = ACM.Lq;

    sm.npp = ACM.npp;
    sm.omg = 0;
}
void ob_init(){

    ob.Js = ACM.Js;
    ob.Js_inv = 1.0/ob.Js;

    ob.R = ACM.R;
    ob.KE = ACM.KE;
    ob.Ld = ACM.Ld;
    ob.Lq = ACM.Lq;

    ob.npp = ACM.npp;
    ob.omg = 0.0;

    ob.eemf_al = 0.0;
    ob.eemf_be = 0.0;
}
void observation(){
}
#endif

