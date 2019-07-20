
#ifndef CONTROLLER_H
#define CONTROLLER_H

#define VC_LOOP_CEILING 40

struct PI_Reg{
   double   Kp;
   double   Ti;
   double   Ki; // Ki = Kp/Ti*TS
   double   i_state;
   double   i_limit;
};
double PI(struct PI_Reg *r, double err);

#if MACHINE_TYPE == INDUCTION_MACHINE
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

    struct PI_Reg pi_speed;
    struct PI_Reg pi_iMs;
    struct PI_Reg pi_iTs;
};
#elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
struct ControllerForExperiment{

    double timebase;

    double ual;
    double ube;

    double R;
    double KE;
    double Ld;
    double Lq;

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

    double uMs_cmd;
    double uTs_cmd;
    double iMs_cmd;
    double iTs_cmd;

    struct PI_Reg pi_speed;
    struct PI_Reg pi_iMs;
    struct PI_Reg pi_iTs;
};
#endif
extern struct ControllerForExperiment CTRL;

void CTRL_init();
void control(double speed_cmd, double speed_cmd_dot);


void cmd_fast_speed_reversal(double timebase, double instant, double interval, double rpm_cmd);
void cmd_slow_speed_reversal(double timebase, double instant, double interval, double rpm_cmd);
#endif
