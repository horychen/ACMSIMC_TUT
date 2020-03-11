#ifndef ADD_OBSERVER_H
#define ADD_OBSERVER_H

#define OB_COEF_G1    -1
#define OB_COEF_G2    0 //0
#define OB_COEF_ELL   1000 //1000
#define OB_COEF_GAMMA 50


/* Macro for External Access Interface */
#define US(X) sm.us[X]
#define IS(X) sm.is[X]
#define US_C(X) sm.us_curr[X]
#define IS_C(X) sm.is_curr[X]
#define US_P(X) sm.us_prev[X]
#define IS_P(X) sm.is_prev[X]

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
    double Lq;
    double KE; // psi_PM;

    double Js;
    double Js_inv;

    double omg_elec; // omg_elec = npp * omg_mech
    double omg_mech;
    double theta_d;

    double eemf_al;
    double eemf_be;
    double eemf_q;

    double xXi[2];
    double xEEMF_dummy[2];
    double xOmg;

    double pll_constructed_eemf_error[2];

    double DeltaL;
    double SigmaL;

    double g1;
    double g2;
    double ell;
    double gamma;
};
extern struct Observer ob;


void sm_init();
void ob_init();
void observation();

#endif
