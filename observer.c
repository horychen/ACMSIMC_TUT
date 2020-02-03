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
    }


    sm.npp = PMSM_NUMBER_OF_POLE_PAIRS;
    sm.npp_inv = 1.0/ACM.npp;

    sm.R  = PMSM_RESISTANCE;
    sm.Ld = PMSM_D_AXIS_INDUCTANCE;
    sm.Lq = PMSM_Q_AXIS_INDUCTANCE;
    sm.KE = PMSM_PERMANENT_MAGNET_FLUX_LINKAGE; // Vs/rad

    sm.Js = PMSM_SHAFT_INERTIA;
    sm.Js_inv = 1./sm.Js;

    sm.omg_elec = 0.0;
    sm.omg_mech = sm.omg_elec * sm.npp_inv;
    sm.theta_d = 0.0;
}

void ob_init(){

    ob.R  = PMSM_RESISTANCE;
    ob.KE = PMSM_D_AXIS_INDUCTANCE;
    ob.Ld = PMSM_Q_AXIS_INDUCTANCE;
    ob.Lq = PMSM_PERMANENT_MAGNET_FLUX_LINKAGE; // Vs/rad

    ob.Js = PMSM_SHAFT_INERTIA;
    ob.Js_inv = 1.0/ob.Js;

    ob.omg_elec = 0.0;
    ob.omg_mech = ob.omg_elec * sm.npp_inv;
    ob.theta_d = 0.0;

    ob.eemf_al = 0.0;
    ob.eemf_be = 0.0;
    ob.eemf_q = 0.0;
}


void observation(){

}
