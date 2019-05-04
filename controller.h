
#ifndef CONTROLLER_H
#define CONTROLLER_H


struct ControllerForExperiment{

    double timebase;

    double ual;
    double ube;

    double rs;
    double rreq;
    double Lsigma;
    double alpha;
    double Lmu;
    double Lmu_inv;

    double Tload;
    double rpm_cmd;

    double Js;
    double Js_inv;
};
extern struct ControllerForExperiment CTRL;



void CTRL_init();


#endif
