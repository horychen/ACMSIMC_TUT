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

    sm.npp = ACM.npp;
    sm.npp_inv = 1.0/ACM.npp;

    sm.Js = ACM.Js;
    sm.Js_inv = 1./sm.Js;

    sm.R = ACM.R;
    sm.KE = ACM.KE;
    sm.Ld = ACM.Ld;
    sm.Lq = ACM.Lq;

    sm.omg_elec = 0;
    sm.omg_mech = sm.omg_elec * sm.npp_inv;
}

void ob_init(){

    ob.Js = ACM.Js;
    ob.Js_inv = 1.0/ob.Js;

    ob.R = ACM.R;
    ob.KE = ACM.KE;
    ob.Ld = ACM.Ld;
    ob.Lq = ACM.Lq;

    ob.omg_elec = 0;
    ob.omg_mech = ob.omg_elec * sm.npp_inv;

    ob.eemf_al = 0.0;
    ob.eemf_be = 0.0;
}


void observation(){
}
