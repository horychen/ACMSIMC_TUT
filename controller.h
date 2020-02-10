
#ifndef CONTROLLER_H
#define CONTROLLER_H

struct PID_Reg{
   double   Kp; // Proportional gain
   double   Ti; // Integral time constant
   double   Ki; // Ki = Kp/Ti*TS (note TS is included here for ease of discrete implementaion)
   double   i_state; // Integral internal state
   double   i_limit; // Output limit
};
double PID(struct PID_Reg *r, double err);

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

    double omg_fb;
    double ial_fb;
    double ibe_fb;
    double psi_mu_al_fb;
    double psi_mu_be_fb;

    double rotor_flux_cmd;

    double omg_ctrl_err;
    double speed_ctrl_err;

    double iMs;
    double iTs;

    double theta_M;
    double cosT;
    double sinT;

    double omega_syn;
    double omega_sl;

    double uMs_cmd;
    double uTs_cmd;
    double iMs_cmd;
    double iTs_cmd;

    struct PID_Reg PID_speed;
    struct PID_Reg PID_iMs;
    struct PID_Reg PID_iTs;
};
extern struct ControllerForExperiment CTRL;

void CTRL_init();
void control(double speed_cmd, double speed_cmd_dot);


void cmd_fast_speed_reversal(double timebase, double instant, double interval, double rpm_cmd);
void cmd_slow_speed_reversal(double timebase, double instant, double interval, double rpm_cmd);
#endif
