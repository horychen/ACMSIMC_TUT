
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

    double R;
    double KE;
    double Ld;
    double Lq;

    // double Tload;
    // double rpm_cmd;

    double npp;
    double Js;
    double Js_inv;

    double omg__fb;
    double ial__fb;
    double ibe__fb;
    double psi_mu_al__fb;
    double psi_mu_be__fb;

    double rotor_flux_cmd;

    double omg_ctrl_err;
    double speed_ctrl_err;

    double cosT;
    double sinT;

    double omega_syn;

    double theta_d__fb;
    double id__fb;
    double iq__fb;
    double ud_cmd;
    double uq_cmd;
    double id_cmd;
    double iq_cmd;

    double Tem;
    double Tem_cmd;

    // double theta_M;
    // double iMs;
    // double iTs;
    // double uMs_cmd;
    // double uTs_cmd;
    // double iMs_cmd;
    // double iTs_cmd;

    struct PID_Reg PID_speed;
    struct PID_Reg PID_id;
    struct PID_Reg PID_iq;
};
extern struct ControllerForExperiment CTRL;

void CTRL_init();
void control(double speed_cmd, double speed_cmd_dot);


void cmd_fast_speed_reversal(double timebase, double instant, double interval, double rpm_cmd);
void cmd_slow_speed_reversal(double timebase, double instant, double interval, double rpm_cmd);
#endif
