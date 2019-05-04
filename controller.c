#include "ACMSim.h"

struct ControllerForExperiment CTRL;

/* Initialization */
void CTRL_init(){
    int i=0,j=0;

    CTRL.timebase = 0.0;

        /* Parameter (including speed) Adaptation */ 
        CTRL.rs     = IM.rs;
        CTRL.rreq   = IM.rreq;
        CTRL.Lsigma = IM.Lsigma;
        CTRL.alpha  = IM.alpha;
        CTRL.Lmu    = IM.Lmu;
        CTRL.Lmu_inv = 1.0/IM.Lmu;
        CTRL.Js     = IM.Js;
        CTRL.Js_inv = 1.0/IM.Js;    

    CTRL.ual = 0.0;
    CTRL.ube = 0.0;

    CTRL.rpm_cmd = 0.0;
}
