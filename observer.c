#include "ACMSim.h"


struct SynchronousMachine sm;
struct Observer ob;

void sm_init(){
    int i;
    for(i=0; i<2; ++i){
        sm.us[i] = 0;
        sm.is[i] = 0;
        sm.us_curr[i] = 0;
        sm.is_curr[i] = 0;
        sm.us_prev[i] = 0;
        sm.is_prev[i] = 0;        
    }

    sm.npp     = PMSM_NUMBER_OF_POLE_PAIRS;
    sm.npp_inv = 1.0/sm.npp;

    sm.R      = PMSM_RESISTANCE;
    sm.Ld     = PMSM_D_AXIS_INDUCTANCE;
    sm.Ld_inv = 1/sm.Ld;
    sm.Lq     = PMSM_Q_AXIS_INDUCTANCE;
    sm.KE     = PMSM_PERMANENT_MAGNET_FLUX_LINKAGE; // Vs/rad

    sm.Js     = PMSM_SHAFT_INERTIA;
    sm.Js_inv = 1./sm.Js;

    sm.omg_elec = 0.0;
    sm.omg_mech = sm.omg_elec * sm.npp_inv;
    sm.theta_d = 0.0;
}

void ob_init(){

    ob.R  = PMSM_RESISTANCE;
    ob.Ld = PMSM_D_AXIS_INDUCTANCE;
    ob.Lq = PMSM_Q_AXIS_INDUCTANCE;
    ob.KE = PMSM_PERMANENT_MAGNET_FLUX_LINKAGE; // Vs/rad

    ob.Js = PMSM_SHAFT_INERTIA;
    ob.Js_inv = 1.0/ob.Js;

    ob.omg_elec = 0.0;
    ob.omg_mech = ob.omg_elec * sm.npp_inv;
    ob.theta_d = 0.0;

    ob.eemf_al = 0.0;
    ob.eemf_be = 0.0;
    ob.eemf_q = 0.0;

    ob.xXi[0] = 0.0;
    ob.xXi[1] = 0.0;
    ob.xEEMF_dummy[0] = 0.0;
    ob.xEEMF_dummy[1] = 0.0;
    ob.xOmg = 0.0;

    ob.pll_constructed_eemf_error[0] = 0.0;
    ob.pll_constructed_eemf_error[1] = 0.0;

    // ob.DeltaL = (ob.Ld - ob.Lq) / 2;
    // ob.SigmaL = (ob.Ld + ob.Lq) / 2;

    ob.g1    = OB_COEF_G1;
    ob.g2    = OB_COEF_G2;
    ob.ell   = OB_COEF_ELL;
    ob.gamma = OB_COEF_GAMMA;
}

#define NUMBER_OF_EEMF_STATES 5
void rhs_func_eemf(double *xxn, double *xXi, double *xEEMF_dummy, double xOmg, double hs){

    #define R      ob.R
    #define LD     sm.Ld
    #define LD_INV sm.Ld_inv
    #define LQ     sm.Lq
    #define g1     ob.g1    // Diagonal gain for eemf observer
    #define g2     ob.g2    // Off-diagonal gain for eemf observer (=0)
    #define ell    ob.ell   // Gain for speed adaptive Luenberger observer
    #define gamma  ob.gamma // Gain for speed update rule
    // #define OMG_USED sm.omg_elec
    #define OMG_USED xOmg

    // EEMF (from reduced order observer) as reference model
    ob.eemf_al = xXi[0] + g1*IS(0);
    ob.eemf_be = xXi[1] + g1*IS(1);

    // Compute mismatch
    ob.pll_constructed_eemf_error[0] = ob.eemf_al - xEEMF_dummy[0];
    ob.pll_constructed_eemf_error[1] = ob.eemf_be - xEEMF_dummy[1];

    // f: 
    // f[0], f[1] are the derivatives of xXi.
    // f[2], f[3] are the derivatives of xEEMF_dummy.
    // f[4] is the derivative of xOmg.
    double f[NUMBER_OF_EEMF_STATES];
    double g1_deriv = 0;

    // xXi: Reduced-order eemf observer with intermediate state xXi
    f[0] = (g1*LD_INV)*xXi[0] + OMG_USED*-xXi[1] - g1*LD_INV*US(0) + ( R*LD_INV*g1 + g1*g1*LD_INV - g1_deriv)*IS(0) + ( LQ*LD_INV * OMG_USED*g1 ) *-IS(1);
    f[1] = (g1*LD_INV)*xXi[1] + OMG_USED* xXi[0] - g1*LD_INV*US(1) + ( R*LD_INV*g1 + g1*g1*LD_INV - g1_deriv)*IS(1) + ( LQ*LD_INV * OMG_USED*g1 ) * IS(0);

    // xEEMF_dummy: Speed adaptive Luenberger observer or PLL (extract frequency info from the reconstructed eemf/xXi)
    f[2] = xOmg*-xEEMF_dummy[1] + ell * ob.pll_constructed_eemf_error[0];
    f[3] = xOmg* xEEMF_dummy[0] + ell * ob.pll_constructed_eemf_error[1];

    // xOmg: Speed update rule
    f[4] = gamma * (-ob.pll_constructed_eemf_error[0]*xEEMF_dummy[1] + ob.pll_constructed_eemf_error[1]*xEEMF_dummy[0]);

    xxn[0] = ( f[0] )*hs;
    xxn[1] = ( f[1] )*hs;
    xxn[2] = ( f[2] )*hs;
    xxn[3] = ( f[3] )*hs;
    xxn[4] = ( f[4] )*hs;

    #undef R
    #undef LD    
    #undef LD_INV
    #undef DeltaL
    #undef SigmaL
    #undef g1    
    #undef g2    
    #undef gamma 
    #undef ell   
}
void rk4_eemf(double hs){
    static double xx1[NUMBER_OF_EEMF_STATES];
    static double xx2[NUMBER_OF_EEMF_STATES];
    static double xx3[NUMBER_OF_EEMF_STATES];
    static double xx4[NUMBER_OF_EEMF_STATES];
    static double x_temp[NUMBER_OF_EEMF_STATES];
    static double *p_x_temp=x_temp;

    /* Theoritically speaking, rhs_func should be time-varing like rhs_func(.,t).
       To apply codes in DSP, we do time-varing updating of IS(0) and IS(1) outside rhs_func(.) to save time. */

    /* 
     * Begin RK4 
     * */
    // time instant t
    US(0) = US_P(0);
    US(1) = US_P(1);
    IS(0) = IS_P(0);
    IS(1) = IS_P(1);
    rhs_func_eemf( xx1, ob.xXi, ob.xEEMF_dummy, ob.xOmg, hs ); 
    x_temp[0]  = ob.xXi[0]           + xx1[0]*0.5;
    x_temp[1]  = ob.xXi[1]           + xx1[1]*0.5;
    x_temp[2]  = ob.xEEMF_dummy[0]   + xx1[2]*0.5;
    x_temp[3]  = ob.xEEMF_dummy[1]   + xx1[3]*0.5;
    x_temp[4]  = ob.xOmg             + xx1[4]*0.5;

    // time instant t+hs/2
    IS(0) = 0.5*(IS_P(0)+IS_C(0));
    IS(1) = 0.5*(IS_P(1)+IS_C(1));
    rhs_func_eemf( xx2, p_x_temp, p_x_temp+2, *(p_x_temp+4), hs );
    x_temp[0]  = ob.xXi[0]           + xx2[0]*0.5;
    x_temp[1]  = ob.xXi[1]           + xx2[1]*0.5;
    x_temp[2]  = ob.xEEMF_dummy[0]   + xx2[2]*0.5;
    x_temp[3]  = ob.xEEMF_dummy[1]   + xx2[3]*0.5;
    x_temp[4]  = ob.xOmg             + xx2[4]*0.5;

    // time instant t+hs/2
    rhs_func_eemf( xx3, p_x_temp, p_x_temp+2, *(p_x_temp+4), hs );
    x_temp[0]  = ob.xXi[0]           + xx3[0];
    x_temp[1]  = ob.xXi[1]           + xx3[1];
    x_temp[2]  = ob.xEEMF_dummy[0]   + xx3[2];
    x_temp[3]  = ob.xEEMF_dummy[1]   + xx3[3];
    x_temp[4]  = ob.xOmg             + xx3[4];

    // time instant t+hs
    IS(0) = IS_C(0);
    IS(1) = IS_C(1);
    rhs_func_eemf( xx4, p_x_temp, p_x_temp+2, *(p_x_temp+4), hs );
    // \+=[^\n]*1\[(\d+)\][^\n]*2\[(\d+)\][^\n]*3\[(\d+)\][^\n]*4\[(\d+)\][^\n]*/ ([\d]+)
    // +=   (xx1[$5] + 2*(xx2[$5] + xx3[$5]) + xx4[$5])*0.166666666666667; // $5
    ob.xXi[0]         +=   (xx1[0] + 2*(xx2[0] + xx3[0]) + xx4[0])*0.166666666666667; // 0
    ob.xXi[1]         +=   (xx1[1] + 2*(xx2[1] + xx3[1]) + xx4[1])*0.166666666666667; // 1
    ob.xEEMF_dummy[0] +=   (xx1[2] + 2*(xx2[2] + xx3[2]) + xx4[2])*0.166666666666667; // 2
    ob.xEEMF_dummy[1] +=   (xx1[3] + 2*(xx2[3] + xx3[3]) + xx4[3])*0.166666666666667; // 3
    ob.xOmg           +=   (xx1[4] + 2*(xx2[4] + xx3[4]) + xx4[4])*0.166666666666667; // 4

    // EEMF
    ob.eemf_al  = ob.xXi[0] + ob.g1*IS(0);
    ob.eemf_be  = ob.xXi[1] + ob.g1*IS(1);
    ob.omg_elec = ob.xOmg;
    ob.omg_mech = ob.omg_elec * sm.npp_inv;
    // ob.theta_d  = atan2(-ob.eemf_al, 
    //                      ob.eemf_be); // 180 deg shift when speed is negative
    // ob.theta_d  = atan2(-ob.eemf_al*sign(ob.omg_elec), 
    //                      ob.eemf_be*sign(ob.omg_elec));
    ob.theta_d  = atan2(-ob.eemf_al*sign(ob.xOmg), 
                         ob.eemf_be*sign(ob.xOmg));
}

void observation(){

    ob.g1 = OB_COEF_G1;
    // if(CTRL.timebase>2){
    //     ob.g1 = -fabs(OB_OMG)*sm.Ld*2;
    // }
    OB_OMG = sm.omg_elec;

    /* OBSERVATION */
    rk4_eemf(TS); 

    /* 备份这个采样点的数据供下次使用。所以，观测的和实际的相比，是延迟一个采样周期的。 */
    //  2017年1月20日，将控制器放到了观测器的后面。
    // * 所以，上一步电压US_P的更新也要延后了。
    // US_P(0) = US_C(0); 
    // US_P(1) = US_C(1);
    IS_P(0) = IS_C(0);
    IS_P(1) = IS_C(1);
}

