#ifndef ADD_OBSERVER_H
#define ADD_OBSERVER_H

#define OB_COEF_K1    (1*1000)   //(10*1000) 
#define OB_COEF_K2    (4*20000)  //(10*20000)
#define OB_COEF_GAMMA (0.5*1e8)   //(160*1e8) 

/* Macro for External Access Interface */
#define US(X) sm.us[X]
#define IS(X) sm.is[X]
#define US_C(X) sm.us_curr[X]
#define IS_C(X) sm.is_curr[X]
#define US_P(X) sm.us_prev[X]
#define IS_P(X) sm.is_prev[X]

#define IS_LPF(X) sm.is_lpf[X]
#define IS_HPF(X) sm.is_hpf[X]
#define IS_BPF(X) sm.is_bpf[X]

#define OB_EEMF_AL ob.eemf_al
#define OB_EEMF_BE ob.eemf_be
#define OB_POS     ob.theta_d
#define OB_OMG     ob.xOmg
#define OB_LD      ob.Ld
#define OB_LQ      ob.Lq
#define OB_R       ob.R

struct SynchronousMachine{
    double us[2];
    double is[2];
    double us_curr[2];
    double is_curr[2];
    double us_prev[2];
    double is_prev[2];
    double is_lpf[2];
    double is_hpf[2];
    double is_bpf[2];

    double current_lpf_register[2];
    double current_hpf_register[2];
    double current_bpf_register1[2];
    double current_bpf_register2[2];

    double npp;
    double npp_inv;

    double R;
    double Ld;
    double Ld_inv;
    double Lq;
    double KE;

    double Js;
    double Js_inv;

    double omg_elec; // omg_elec = npp * omg_mech
    double omg_mech;
    double theta_d;
};
extern struct SynchronousMachine sm;

struct Observer{

    double R;
    double Ld;
    double Ld_inv;
    double Lq;
    double KE; // psi_PM;

    double DeltaL;

    double Js;
    double Js_inv;

    double omg_elec; // omg_elec = npp * omg_mech
    double omg_mech;
    double theta_d;

    double eemf_al;
    double eemf_be;



    double xPsi[2];
    double xChi[2];

    double xEta[2];
    double xVarSigma[2];

    double xUpsilon[2];
    double xZeta[2];

    double xOmg;


    double output_error[2];
    double output_error_eff[2];

    double k1;
    double k2;
    double gamma_omega;

};
extern struct Observer ob;


void sm_init();
void ob_init();
void observation();


// ----------------------------- HFSI


#define HFSI_ON
#define SENSORLESS_CONTROL_HFSI true

#define HFSI_VOLTAGE 5 // V
#define HFSI_CEILING 0 // Essential to be 1 because we assume h = TS, but in fact h = HFSI_CEILING * TS (?)

#define LPF_TIME_CONST_INVERSE (5*2*M_PI) // time constant is 1/400 <=> cutoff frequency is 400/(2*pi) ||| 换句话说，截止频率 * 2pi = 时间常数的倒数 |||| f_cutoff = 1 / (time constant * 2*pi)
#define LUENBERGER_GAIN_1 (45)   // 45       //45     // 30     // 30       // 30     // 30    // 20  // Large gain to position will cause steady state position error, but increase it close to limit
#define LUENBERGER_GAIN_2 (750)  // (500)    //(500)  // (900)  // (750)    // (300)  // (300) // 100 // If speed estimate has too much dynamics during reversal, you need to increase this gain actually...
#define LUENBERGER_GAIN_3 (6000) // (0*1800) //(2000) // (1500) // (0*6000) // (1500) // (790) // 500 // Tune reversal response to slight over-shoot

void hfsi_do();
void hfsi_init();
struct HFSI_Data{
    double test_signal_al;
    double test_signal_be;
    double test_signal_M;
    double test_signal_T;
    double M_lpf;
    double T_lpf;
    double M_hpf;
    double T_hpf;
    double theta_filter;
    double theta_d_raw;
    double theta_d;
    double omg_elec;
    double pseudo_load_torque;
    double mismatch;

    double tilde_theta_d;
    double square_wave_internal_register;
    double LAST_uM;

    double h;
};
extern struct HFSI_Data hfsi;
double difference_between_two_angles(double first, double second);
void dynamics_lpf(double input, double *state, double *derivative);
void RK4_111_general(void (*pointer_dynamics)(), double input, double *state, double hs);

#endif

