#include "ACMSim.h"


struct SynchronousMachine sm;
struct Observer ob;

void sm_init(){
    int i;
    for(i=0; i<2; ++i){
        sm.us[i] = 0;
        sm.is[i] = 0;
        sm.us_curr[i] = 0;
        sm.is_curr[i] = 0;
        sm.us_prev[i] = 0;
        sm.is_prev[i] = 0;
        sm.is_lpf[i]  = 0;
        sm.is_hpf[i]  = 0;
        sm.is_bpf[i]  = 0;

        sm.current_lpf_register[i] = 0;
        sm.current_hpf_register[i] = 0;
        sm.current_bpf_register1[i] = 0;
        sm.current_bpf_register2[i] = 0;
    }

    sm.npp     = PMSM_NUMBER_OF_POLE_PAIRS;
    sm.npp_inv = 1.0/sm.npp;

    sm.R      = PMSM_RESISTANCE;
    sm.Ld     = PMSM_D_AXIS_INDUCTANCE;
    sm.Ld_inv = 1/sm.Ld;
    sm.Lq     = PMSM_Q_AXIS_INDUCTANCE;
    sm.KE     = PMSM_PERMANENT_MAGNET_FLUX_LINKAGE; // Vs/rad

    sm.Js     = PMSM_SHAFT_INERTIA;
    sm.Js_inv = 1./sm.Js;

    sm.omg_elec = 0.0;
    sm.omg_mech = sm.omg_elec * sm.npp_inv;
    sm.theta_d = 0.0;
}



