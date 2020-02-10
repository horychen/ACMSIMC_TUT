#ifndef ADD_OBSERVER_H
#define ADD_OBSERVER_H

#if OBSERVER_APPLIED == FIRST_APPROXIMATION_ANDERSON86 

/* Macro for External Access Interface */
#define US(X) im.us[X]
#define IS(X) im.is[X]
#define US_C(X) im.us_curr[X]
#define IS_C(X) im.is_curr[X]
#define US_P(X) im.us_prev[X]
#define IS_P(X) im.is_prev[X]

#define OB_RS ob.xTheta[0]
#define OB_RREQ (ob.xTheta[1]*im.Lmu)
#define OB_ALPHA ob.xTheta[1]
#define OB_OMG ob.xTheta[2]

#define OB_LMU (im.Lmu*1.0)
#define OB_LMU_INV (im.Lmu_inv/1.0)
#define OB_LSIGMA im.Lsigma

#define OB_FLUX(X) ob.xPsiMu[X]
// #define OB_TLOAD ob.taao_Tload
// #define OB_FILTERED_SPEED ob.taao_speed_lpf ?



// FLUX COMMAND
#define IM_FLUX_COMMAND_VALUE 1.2 // 1.7 for low speed rege 磁链太小发电带载会不稳 // 0.8 // 0.6 is worse than 1.0 in the sense of convergence arcuracy of alpha parameter
#define IM_FLUX_COMMAND_ON (true)
#define M1 (0*0.05) 
#define OMG1 (3.5*2*M_PI) 

#define K1_VALUE (1*1000) 
#define K2_VALUE (1*16000) // 2*16000 for simulated slow reversal
#define GAMMA_1 (0*5e6)
#define GAMMA_2 (0*5e6)
#define GAMMA_3 (0*5e7)

// #define K1_VALUE 10*5 
// #define K2_VALUE 300*3
// #define SCALE_I (0)
// #define GAMMA_1 (SCALE_I*1e5)
// #define GAMMA_2 (5e4)
// #define GAMMA_3 5e5

#define ADAPT_AT_30SEC false

struct InductionMachine{
    double us[2];
    double is[2];
    double us_curr[2];
    double is_curr[2];
    double us_prev[2];
    double is_prev[2];

    double Js;
    double Js_inv;
    // double Lm;
    // double Lm_inv;
    // double Lls;
    // double Llr;
    double Ls;
    // double Lr;
    // double sigma;
    double rs;
    // double rr;
    // double Tr;
    double alpha;
    double Lsigma;
    double Lsigma_inv;
    double Lmu;
    double Lmu_inv;
    double rreq;

    double npp;
    double omg;

    double theta_r;
    double theta_d;
};
extern struct InductionMachine im;

struct Observer{

    double k1; 
    double k2;

    double xPsiSigma[2];    // \psi_\sigma
    double xChi[2];         // \chi

    double xPsiMu[2];       // \psi_\mu
    double xUps_al[3];     // The alpha component of the filtered regressors
    double xUps_be[3];     // The beta component of the filtered regressors
    double xZeta_al[3];     // The alpha component of temporary states
    double xZeta_be[3];     // 

    double xEta[2];     // for effective mismtach
    double xVarsigma[2];     // The vector of temporary states

    double mismatch[3];
    double mismatch_eff[2];
    double error[2];

    double gamma[3];
    double xTheta[3];

    double taao_omg_integralPart;
    double taao_speed;

    double timebase;
    double Ts;

    double omega_e;
    double Tem;

    double taao_flux_cmd;
    int taao_flux_cmd_on;

    double cosT;
    double sinT;
    double theta_M;
};
extern struct Observer ob;

void acm_init();
void ob_init();
void observation();

double IM_FluxModulusCommand();

#endif

#endif
