#include "ACMSim.h"

struct InductionMachine im;
struct Observer ob;


void im_init(){
    int i;
    for(i=0; i<2; ++i){
        im.us[i] = 0;
        im.is[i] = 0;
        im.us_curr[i] = 0;
        im.is_curr[i] = 0;
        im.us_prev[i] = 0;
        im.is_prev[i] = 0;        
    }

    im.Js = IM.Js;
    im.Js_inv = 1./im.Js;
    im.rs = IM.rs;
    im.rreq = IM.rreq;
    im.alpha = IM.alpha;
    im.Lsigma = IM.Lsigma;
    im.Lsigma_inv = 1/im.Lsigma;
    im.Lmu = IM.Lmu;
    im.Lmu_inv = IM.Lmu_inv;
    im.npp = IM.npp;
    im.omg = 0;
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
    ob.tajima.e_M        = uMs_cmd - ob.rs*iMs + ob.tajima.omega_syn*ob.Lsigma*iTs;
    ob.tajima.omega_syn -= ob.tajima.K_PEM * ob.tajima.e_M;
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

