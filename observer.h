#ifndef ADD_OBSERVER_H
#define ADD_OBSERVER_H




/* Macro for External Access Interface */
#define US(X) sm.us[X]
#define IS(X) sm.is[X]
#define US_C(X) sm.us_curr[X]
#define IS_C(X) sm.is_curr[X]
#define US_P(X) sm.us_prev[X]
#define IS_P(X) sm.is_prev[X]

struct SynchronousMachine{
    double us[2];
    double is[2];
    double us_curr[2];
    double is_curr[2];
    double us_prev[2];
    double is_prev[2];

    double npp;
    double npp_inv;

    double Js;
    double Js_inv;

    double R;
    double KE;
    double Ld;
    double Lq;

    double omg_elec; // omg_elec = npp * omg_mech
    double omg_mech;
    double theta_d;
};
extern struct SynchronousMachine sm;

struct Observer{
    double Js;
    double Js_inv;

    double R;
    double KE;
    double psi_PM;
    double Ld;
    double Lq;

    double omg_elec; // omg_elec = npp * omg_mech
    double omg_mech;
    double theta_d;

    double eemf_al;
    double eemf_be;
};
extern struct Observer ob;


void sm_init();
void ob_init();
void observation();

#endif
