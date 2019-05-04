#ifndef ADD_OBSERVER_H
#define ADD_OBSERVER_H

/* Macro for External Access Interface */
#define US(X) im.us[X]
#define IS(X) im.is[X]
#define US_C(X) im.us_curr[X]
#define IS_C(X) im.is_curr[X]
#define US_P(X) im.us_prev[X]
#define IS_P(X) im.is_prev[X]

struct InductionMachine{
    double us[2];
    double is[2];
    double us_curr[2];
    double is_curr[2];
    double us_prev[2];
    double is_prev[2];

    double Js;
    double Js_inv;

    double rs;
    double rreq;

    double alpha;
    double Lsigma;
    double Lsigma_inv;
    double Lmu;
    double Lmu_inv;

    double npp;
    double omg;
};
extern struct InductionMachine im;

#endif

