#ifndef ACMSIM_H
#define ACMSIM_H

/* standard lib */
// #include <stdbool.h> // bool for _Bool and true for 1
#include <stdio.h> // printf #include <stdbool.h> // bool for _Bool and true for 1
#include <process.h>//reqd. for system function prototype
#include <conio.h> // for clrscr, and getch()
#include "stdlib.h" // for rand()
#include "math.h"
#include "time.h"

/* Macro for Part transformation*/
#define AB2M(A, B, COS, SIN)  ( (A)*COS  + (B)*SIN )
#define AB2T(A, B, COS, SIN)  ( (A)*-SIN + (B)*COS )
#define MT2A(M, T, COS, SIN)  ( (M)*COS - (T)*SIN )
#define MT2B(M, T, COS, SIN)  ( (M)*SIN + (T)*COS )

// Macro for Power-invariant inverse Clarke transformation
#define AB2U(A, B) ( 0.816496580927726 * ( A ) )
#define AB2V(A, B) ( 0.816496580927726 * ( A*-0.5 + B*0.8660254037844387 ) )
#define AB2W(A, B) ( 0.816496580927726 * ( A*-0.5 + B*-0.8660254037844385 ) )

/* General Constants */
#define TRUE  True
#define FALSE False
#define True  true
#define False false
#define true  (1)
#define false (0)
#define ONE_OVER_2PI          0.15915494309189535 // 1/(2*pi)
#define TWO_PI_OVER_3         2.0943951023931953
#define SIN_2PI_SLASH_3       0.86602540378443871 // sin(2*pi/3)
#define SIN_DASH_2PI_SLASH_3 -0.86602540378443871 // sin(-2*pi/3)
#define SQRT_2_SLASH_3        0.81649658092772603 // sqrt(2.0/3.0)
#define abs                   use_fabs_instead_or_you_will_regret
// #define RPM_2_RAD_PER_SEC     0.20943951023931953 // (2/60*Tpi)
// #define RAD_PER_SEC_2_RPM     4.7746482927568605 // 1/(im.npp/60*Tpi)
#define RAD_PER_SEC_2_RPM ( 60.0/(2*M_PI*ACM.npp) )
#define RPM_2_RAD_PER_SEC ( (2*M_PI*ACM.npp)/60.0 )
#define M_PI_OVER_180   0.017453292519943295
// #define PI_D 3.1415926535897932384626433832795 /* double */

// 补充的宏，为了实现实验/仿真代码大一统
#define Uint32 unsigned long int
#define Uint16 unsigned int

#define NUMBER_OF_STATES 4 // valid for PMSM
#define PHASE_NUMBER 3 // 3 phase machine



// Everthing else is in here
#include "ACMConfig.h"



#define UAL_C_DIST ACM.ual_c_dist // alpha-component of the distorted phase voltage = sine voltage + distored component
#define UBE_C_DIST ACM.ube_c_dist
#define DIST_AL ACM.dist_al // alpha-component of the distorted component of the voltage
#define DIST_BE ACM.dist_be

struct SynchronousMachineSimulated{

    double Ts;

    double x[NUMBER_OF_STATES];
    double x_dot[NUMBER_OF_STATES];

    double omg_elec;
    double rpm;
    double rpm_cmd;
    double rpm_deriv_cmd;
    double Tload;
    double Tem;

    double npp;

    double R;
    double Ld;
    double Lq;
    double KE; // psi_PM
    double Js;

    double mu_m;
    double L0;
    double L1;

    double ual;
    double ube;
    double ial;
    double ibe;

    double ual_c_dist;
    double ube_c_dist;
    double dist_al;
    double dist_be;

    double theta_d;
    double ud;
    double uq;
    double id;
    double iq;

    double eemf_q;
    double eemf_al;
    double eemf_be;
    double theta_d__eemf;
};
extern struct SynchronousMachineSimulated ACM;

#include "controller.h"
#include "observer.h"
// #include "Add_TAAO.h"
// #include "utility.h"
#include "inverter.h"
#include "commissioning.h"

// Saturation
// #include "satlut.h"
// #define ACMSIMC_DEBUG false


/* Declaration of Utility Function */
void write_header_to_file(FILE *fw);
void write_data_to_file(FILE *fw);
int isNumber(double x);
double sign(double x);
double fabs(double x);
#endif