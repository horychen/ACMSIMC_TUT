#include "ACMSim.h"

#if OBSERVER_APPLIED == FIRST_APPROXIMATION_ANDERSON86 

struct InductionMachine im;
struct Observer ob;
struct u_reg
{
    double u11,u12,u13,
            u21,u22,u23,
            u31,u32,u33;
} r;

/* In DSP, these shall be defines elsewhere. */
struct InductionMachine im;
struct Observer ob;

/* rk16 of 3ParamIdRevisited with Narendra1973
 * */
double deriv_theta[3];
void rhs_func_Zheng99(double *xxn, double *xPsiSigma, double *xChi,
        double *xUps_al, double *xUps_be, double *xZeta_al, double *xZeta_be, 
        double *xEta, double *xVarsigma, double *xTheta, double hs){

    // Compute mismatch
    ob.mismatch[0] = im.Lsigma*IS(0) - xPsiSigma[0];
    ob.mismatch[1] = im.Lsigma*IS(1) - xPsiSigma[1];

    /* The Filtered Regressor - Upsilon */
    #define xRs xTheta[0]
    #define xAlpha xTheta[1]
    #define xOmg xTheta[2]
    xxn[4]  = ( -ob.k1*xUps_al[0] + xZeta_al[0] - IS(0) )*hs;
    xxn[7]  = ( -ob.k1*xUps_be[0] + xZeta_be[0] - IS(1) )*hs;
    xxn[10] = ( -ob.k2*xUps_al[0] - xAlpha*IS(0) + xOmg*-IS(1) )*hs;
    xxn[13] = ( -ob.k2*xUps_be[0] - xAlpha*IS(1) + xOmg*IS(0) )*hs;

    xxn[5]  = ( -ob.k1*xUps_al[1] + xZeta_al[1] - im.Ls*IS(0) )*hs;
    xxn[8]  = ( -ob.k1*xUps_be[1] + xZeta_be[1] - im.Ls*IS(1) )*hs;
    xxn[11] = ( -ob.k2*xUps_al[1] + US(0) - xRs*IS(0) )*hs;
    xxn[14] = ( -ob.k2*xUps_be[1] + US(1) - xRs*IS(1) )*hs;

    xxn[6]  = ( -ob.k1*xUps_al[2] + xZeta_al[2] + im.Lsigma*IS(1) )*hs;
    xxn[9]  = ( -ob.k1*xUps_be[2] + xZeta_be[2] - im.Lsigma*IS(0) )*hs;
    xxn[12] = ( -ob.k2*xUps_al[2] + (US(1) - xRs*IS(1)) )*hs;
    xxn[15] = ( -ob.k2*xUps_be[2] - (US(0) - xRs*IS(0)) )*hs;

    /* The xEta */
    xxn[16] = ( -ob.k1*xEta[0] + xVarsigma[0] - IS(0)*xRs - im.Ls*IS(0)*xAlpha + im.Lsigma*IS(1)*xOmg )*hs;
    xxn[17] = ( -ob.k1*xEta[1] + xVarsigma[1] - IS(1)*xRs - im.Ls*IS(1)*xAlpha - im.Lsigma*IS(0)*xOmg )*hs;
    xxn[18] = ( -ob.k2*xEta[0] + (-xAlpha*IS(0) + xOmg*-IS(1))*xRs 
                               +         (US(0) - xRs*IS(0))  *xAlpha
                               +         (US(1) - xRs*IS(1))  *xOmg )*hs;
    xxn[19] = ( -ob.k2*xEta[1] + (-xAlpha*IS(1) + xOmg*IS(0))*xRs
                               +         (US(1) - xRs*IS(1)) *xAlpha
                               +       (-(US(0) - xRs*IS(0)))*xOmg )*hs;

    /* hat theta derivative */
    ob.mismatch_eff[0] = ob.mismatch[0] + xEta[0] - (xUps_al[0]*ob.xTheta[0]+xUps_al[1]*ob.xTheta[1]+xUps_al[2]*ob.xTheta[2]);
    ob.mismatch_eff[1] = ob.mismatch[1] + xEta[1] - (xUps_be[0]*ob.xTheta[0]+xUps_be[1]*ob.xTheta[1]+xUps_be[2]*ob.xTheta[2]);
    deriv_theta[0] = ob.gamma[0] * (ob.mismatch_eff[0]*xUps_al[0] + ob.mismatch_eff[1]*xUps_be[0]);
    deriv_theta[1] = ob.gamma[1] * (ob.mismatch_eff[0]*xUps_al[1] + ob.mismatch_eff[1]*xUps_be[1]);
    deriv_theta[2] = ob.gamma[2] * (ob.mismatch_eff[0]*xUps_al[2] + ob.mismatch_eff[1]*xUps_be[2]);

    xxn[20] = (deriv_theta[0])*hs;
    xxn[21] = (deriv_theta[1])*hs;
    xxn[22] = (deriv_theta[2])*hs;

    /* Full Order Observer */
    xxn[0] = ( xChi[0] + US(0) - IS(0)*xRs - im.Ls*IS(0)*xAlpha - im.Lsigma* IS(1) *xOmg + ob.k1*ob.mismatch[0] )*hs; 
    xxn[1] = ( xChi[1] + US(1) - IS(1)*xRs - im.Ls*IS(1)*xAlpha + im.Lsigma* IS(0) *xOmg + ob.k1*ob.mismatch[1] )*hs;
    xxn[2] = (                         (US(0)-xRs*IS(0))*xAlpha + (US(1)-xRs*IS(1))*xOmg + ob.k2*ob.mismatch[0] )*hs;
    xxn[3] = (                         (US(1)-xRs*IS(1))*xAlpha - (US(0)-xRs*IS(0))*xOmg + ob.k2*ob.mismatch[1] )*hs;

    #undef xRs 
    #undef xAlpha
    #undef xOmg
}
void rK4(double hs){
    static double xx1[23];
    static double xx2[23];
    static double xx3[23];
    static double xx4[23];
    static double x_temp[23];
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

    rhs_func_Zheng99( xx1, ob.xPsiSigma, ob.xChi, ob.xUps_al, ob.xUps_be, ob.xZeta_al, ob.xZeta_be, 
        ob.xEta, ob.xVarsigma, ob.xTheta, hs); 
    x_temp[0]  = ob.xPsiSigma[0]          + xx1[0]*0.5;
    x_temp[1]  = ob.xPsiSigma[1]          + xx1[1]*0.5;
    x_temp[2]  = ob.xChi[0]               + xx1[2]*0.5;
    x_temp[3]  = ob.xChi[1]               + xx1[3]*0.5;
    x_temp[4]  = ob.xUps_al[0]            + xx1[4]*0.5;
    x_temp[5]  = ob.xUps_al[1]            + xx1[5]*0.5;
    x_temp[6]  = ob.xUps_al[2]            + xx1[6]*0.5;
    x_temp[7]  = ob.xUps_be[0]            + xx1[7]*0.5;
    x_temp[8]  = ob.xUps_be[1]            + xx1[8]*0.5;
    x_temp[9]  = ob.xUps_be[2]            + xx1[9]*0.5;
    x_temp[10] = ob.xZeta_al[0]           + xx1[10]*0.5;
    x_temp[11] = ob.xZeta_al[1]           + xx1[11]*0.5;
    x_temp[12] = ob.xZeta_al[2]           + xx1[12]*0.5;
    x_temp[13] = ob.xZeta_be[0]           + xx1[13]*0.5;
    x_temp[14] = ob.xZeta_be[1]           + xx1[14]*0.5;
    x_temp[15] = ob.xZeta_be[2]           + xx1[15]*0.5;
    x_temp[16] = ob.xEta[0]               + xx1[16]*0.5;
    x_temp[17] = ob.xEta[1]               + xx1[17]*0.5;
    x_temp[18] = ob.xVarsigma[0]          + xx1[18]*0.5;
    x_temp[19] = ob.xVarsigma[1]          + xx1[19]*0.5;
    x_temp[20] = ob.xTheta[0]             + xx1[20]*0.5;
    x_temp[21] = ob.xTheta[1]             + xx1[21]*0.5;
    x_temp[22] = ob.xTheta[2]             + xx1[22]*0.5;

    // time instant t+hs/2
    IS(0) = 0.5*(IS_P(0)+IS_C(0));
    IS(1) = 0.5*(IS_P(1)+IS_C(1));
    rhs_func_Zheng99( xx2, p_x_temp, p_x_temp+2, p_x_temp+4, p_x_temp+7, p_x_temp+10, p_x_temp+13, 
        p_x_temp+16, p_x_temp+18, p_x_temp+20, hs );
    x_temp[0]  = ob.xPsiSigma[0]          + xx2[0]*0.5;
    x_temp[1]  = ob.xPsiSigma[1]          + xx2[1]*0.5;
    x_temp[2]  = ob.xChi[0]               + xx2[2]*0.5;
    x_temp[3]  = ob.xChi[1]               + xx2[3]*0.5;
    x_temp[4]  = ob.xUps_al[0]            + xx2[4]*0.5;
    x_temp[5]  = ob.xUps_al[1]            + xx2[5]*0.5;
    x_temp[6]  = ob.xUps_al[2]            + xx2[6]*0.5;
    x_temp[7]  = ob.xUps_be[0]            + xx2[7]*0.5;
    x_temp[8]  = ob.xUps_be[1]            + xx2[8]*0.5;
    x_temp[9]  = ob.xUps_be[2]            + xx2[9]*0.5;
    x_temp[10] = ob.xZeta_al[0]           + xx2[10]*0.5;
    x_temp[11] = ob.xZeta_al[1]           + xx2[11]*0.5;
    x_temp[12] = ob.xZeta_al[2]           + xx2[12]*0.5;
    x_temp[13] = ob.xZeta_be[0]           + xx2[13]*0.5;
    x_temp[14] = ob.xZeta_be[1]           + xx2[14]*0.5;
    x_temp[15] = ob.xZeta_be[2]           + xx2[15]*0.5;
    x_temp[16] = ob.xEta[0]               + xx2[16]*0.5;
    x_temp[17] = ob.xEta[1]               + xx2[17]*0.5;
    x_temp[18] = ob.xVarsigma[0]          + xx2[18]*0.5;
    x_temp[19] = ob.xVarsigma[1]          + xx2[19]*0.5;
    x_temp[20] = ob.xTheta[0]             + xx2[20]*0.5;
    x_temp[21] = ob.xTheta[1]             + xx2[21]*0.5;
    x_temp[22] = ob.xTheta[2]             + xx2[22]*0.5;

    // time instant t+hs/2
    rhs_func_Zheng99( xx3, p_x_temp, p_x_temp+2, p_x_temp+4, p_x_temp+7, p_x_temp+10, p_x_temp+13, 
        p_x_temp+16, p_x_temp+18, p_x_temp+20, hs );
    x_temp[0]  = ob.xPsiSigma[0]          + xx3[0];
    x_temp[1]  = ob.xPsiSigma[1]          + xx3[1];
    x_temp[2]  = ob.xChi[0]               + xx3[2];
    x_temp[3]  = ob.xChi[1]               + xx3[3];
    x_temp[4]  = ob.xUps_al[0]            + xx3[4];
    x_temp[5]  = ob.xUps_al[1]            + xx3[5];
    x_temp[6]  = ob.xUps_al[2]            + xx3[6];
    x_temp[7]  = ob.xUps_be[0]            + xx3[7];
    x_temp[8]  = ob.xUps_be[1]            + xx3[8];
    x_temp[9]  = ob.xUps_be[2]            + xx3[9];
    x_temp[10] = ob.xZeta_al[0]           + xx3[10];
    x_temp[11] = ob.xZeta_al[1]           + xx3[11];
    x_temp[12] = ob.xZeta_al[2]           + xx3[12];
    x_temp[13] = ob.xZeta_be[0]           + xx3[13];
    x_temp[14] = ob.xZeta_be[1]           + xx3[14];
    x_temp[15] = ob.xZeta_be[2]           + xx3[15];
    x_temp[16] = ob.xEta[0]               + xx3[16];
    x_temp[17] = ob.xEta[1]               + xx3[17];
    x_temp[18] = ob.xVarsigma[0]          + xx3[18];
    x_temp[19] = ob.xVarsigma[1]          + xx3[19];
    x_temp[20] = ob.xTheta[0]             + xx3[20];
    x_temp[21] = ob.xTheta[1]             + xx3[21];
    x_temp[22] = ob.xTheta[2]             + xx3[22];


    // time instant t+hs
    IS(0) = IS_C(0);
    IS(1) = IS_C(1);
    rhs_func_Zheng99( xx4, p_x_temp, p_x_temp+2, p_x_temp+4, p_x_temp+7, p_x_temp+10, p_x_temp+13, 
        p_x_temp+16, p_x_temp+18, p_x_temp+20, hs );
    // \+=[^\n]*1\[(\d+)\][^\n]*2\[(\d+)\][^\n]*3\[(\d+)\][^\n]*4\[(\d+)\][^\n]*/ ([\d]+)
    // +=   (xx1[$5] + 2*(xx2[$5] + xx3[$5]) + xx4[$5])*0.166666666666667; // $5
    ob.xPsiSigma[0]       +=   (xx1[0] + 2*(xx2[0] + xx3[0]) + xx4[0])*0.166666666666667; // 0
    ob.xPsiSigma[1]       +=   (xx1[1] + 2*(xx2[1] + xx3[1]) + xx4[1])*0.166666666666667; // 1
    ob.xChi[0]            +=   (xx1[2] + 2*(xx2[2] + xx3[2]) + xx4[2])*0.166666666666667; // 2
    ob.xChi[1]            +=   (xx1[3] + 2*(xx2[3] + xx3[3]) + xx4[3])*0.166666666666667; // 3
    ob.xUps_al[0]         +=   (xx1[4] + 2*(xx2[4] + xx3[4]) + xx4[4])*0.166666666666667; // 4
    ob.xUps_al[1]         +=   (xx1[5] + 2*(xx2[5] + xx3[5]) + xx4[5])*0.166666666666667; // 5
    ob.xUps_al[2]         +=   (xx1[6] + 2*(xx2[6] + xx3[6]) + xx4[6])*0.166666666666667; // 6
    ob.xUps_be[0]         +=   (xx1[7] + 2*(xx2[7] + xx3[7]) + xx4[7])*0.166666666666667; // 7
    ob.xUps_be[1]         +=   (xx1[8] + 2*(xx2[8] + xx3[8]) + xx4[8])*0.166666666666667; // 8
    ob.xUps_be[2]         +=   (xx1[9] + 2*(xx2[9] + xx3[9]) + xx4[9])*0.166666666666667; // 9
    ob.xZeta_al[0]        +=   (xx1[10] + 2*(xx2[10] + xx3[10]) + xx4[10])*0.166666666666667; // 10
    ob.xZeta_al[1]        +=   (xx1[11] + 2*(xx2[11] + xx3[11]) + xx4[11])*0.166666666666667; // 11
    ob.xZeta_al[2]        +=   (xx1[12] + 2*(xx2[12] + xx3[12]) + xx4[12])*0.166666666666667; // 12
    ob.xZeta_be[0]        +=   (xx1[13] + 2*(xx2[13] + xx3[13]) + xx4[13])*0.166666666666667; // 13
    ob.xZeta_be[1]        +=   (xx1[14] + 2*(xx2[14] + xx3[14]) + xx4[14])*0.166666666666667; // 14
    ob.xZeta_be[2]        +=   (xx1[15] + 2*(xx2[15] + xx3[15]) + xx4[15])*0.166666666666667; // 15
    ob.xEta[0]            +=   (xx1[16] + 2*(xx2[16] + xx3[16]) + xx4[16])*0.166666666666667; // 16
    ob.xEta[1]            +=   (xx1[17] + 2*(xx2[17] + xx3[17]) + xx4[17])*0.166666666666667; // 17
    ob.xVarsigma[0]       +=   (xx1[18] + 2*(xx2[18] + xx3[18]) + xx4[18])*0.166666666666667; // 18
    ob.xVarsigma[1]       +=   (xx1[19] + 2*(xx2[19] + xx3[19]) + xx4[19])*0.166666666666667; // 19
    ob.xTheta[0]          +=   (xx1[20] + 2*(xx2[20] + xx3[20]) + xx4[20])*0.166666666666667; // 20
    ob.xTheta[1]          +=   (xx1[21] + 2*(xx2[21] + xx3[21]) + xx4[21])*0.166666666666667; // 21
    ob.xTheta[2]          +=   (xx1[22] + 2*(xx2[22] + xx3[22]) + xx4[22])*0.166666666666667; // 22

    if(ob.xTheta[0]>10){
        ob.xTheta[0] = 10;
    }
    
    /* 其他电气参数 */
    // ob.xPsiMu[0] = IM.x[2]*IM.Lm_slash_Lr;
    // ob.xPsiMu[1] = IM.x[3]*IM.Lm_slash_Lr; 

    ob.xPsiMu[0] = (OB_ALPHA*ob.xChi[0]-OB_OMG*ob.xChi[1])/(OB_ALPHA*OB_ALPHA+OB_OMG*OB_OMG) - ob.xPsiSigma[0];
    ob.xPsiMu[1] = (OB_ALPHA*ob.xChi[1]+OB_OMG*ob.xChi[0])/(OB_ALPHA*OB_ALPHA+OB_OMG*OB_OMG) - ob.xPsiSigma[1]; 

    ob.taao_speed = OB_OMG * RAD_PER_SEC_2_RPM;

}


/* Main Flow of Observer
 * */
void observation(){

    #if ADAPT_AT_30SEC
        if(ob.timebase<20){
            ob.gamma[0] = 0;
            ob.gamma[1] = 0;
        }else if(fabs(ob.timebase*10000-200000)<TS){

            ob.gamma[0] = GAMMA_1;
            ob.gamma[1] = GAMMA_2;

            OB_RS = im.rs*0.5;
        }
    #endif

    /* OBSERVATION */
    if(GAMMA_3==0){
        OB_OMG = im.omg;
    }
    rK4(ob.Ts); 
}





/* Initialization Codes*/
void ob_init(){

    ob.k1 = K1_VALUE;
    ob.k2 = K2_VALUE;

    ob.xPsiSigma[0] = 0.0; 
    ob.xPsiSigma[1] = 0.0; 
    ob.xChi[0] = 0.0;
    ob.xChi[1] = 0.0;

    ob.xPsiMu[0] = 0.0;
    ob.xPsiMu[1] = 0.0;

    int i;
    for(i=0;i<3;++i){
        ob.xUps_al[i] = 0.0;
        ob.xUps_be[i] = 0.0;
        ob.xZeta_al[i] = 0.0;
        ob.xZeta_be[i] = 0.0;
    }

    ob.xEta[0] = 0.0;
    ob.xEta[1] = 0.0;
    ob.xVarsigma[0] = 0.0;
    ob.xVarsigma[1] = 0.0;

    ob.mismatch[0] = 0.0;
    ob.mismatch[1] = 0.0;
    ob.mismatch[2] = 0.0;
    ob.mismatch_eff[0] = 0.0;
    ob.mismatch_eff[1] = 0.0;
    ob.error[0] = 0.0;
    ob.error[1] = 0.0;

    ob.gamma[0] = GAMMA_1;
    ob.gamma[1] = GAMMA_2;
    ob.gamma[2] = GAMMA_3;

    ob.timebase = 0.0;
    ob.Ts = TS;

    ob.xTheta[0] = im.rs;
    ob.xTheta[1] = im.alpha;
    ob.xTheta[2] = im.omg;
    ob.taao_omg_integralPart = 0.0;
    ob.taao_speed = OB_OMG * RAD_PER_SEC_2_RPM; 

    ob.omega_e = 0.0;
    ob.Tem = 0.0;

    ob.taao_flux_cmd = IM_FLUX_COMMAND_VALUE;
    ob.taao_flux_cmd_on = IM_FLUX_COMMAND_ON;

    ob.cosT = 1.0;
    ob.sinT = 0.0;
    ob.theta_M = 0.0;
}

void acm_init(){
    int i;
    for(i=0; i<2; ++i){
        im.us[i] = 0;
        im.is[i] = 0;
        im.us_curr[i] = 0;
        im.is_curr[i] = 0;
        im.us_prev[i] = 0;
        im.is_prev[i] = 0;        
    }

    im.Js = ACM.Js;
    im.Js_inv = 1./im.Js;
    im.rs = ACM.rs;
    im.rreq = ACM.rreq;
    im.alpha = ACM.alpha;
    im.Lsigma = ACM.Lsigma;
    im.Lsigma_inv = 1/im.Lsigma;

    im.Lmu = ACM.Lmu;
    im.Ls = im.Lmu + im.Lsigma; //////////////////////// Dependence on this variable in dynamics is bad

    printf("%g, %g, %g", im.Ls, im.Lmu, ACM.Lmu);

    im.Lmu_inv = ACM.Lmu_inv;
    im.npp = ACM.npp;
    im.omg = 0;
    im.theta_r = 0.0;

    im.omg = 0.0;
    im.npp = ACM.npp; // no. of pole pair

    im.theta_r = 0.0;
    im.theta_d = 0.0;
}

double IM_FluxModulusCommand(){

    if(ob.taao_flux_cmd_on){
        return IM_FLUX_COMMAND_VALUE + M1*sin(OMG1*ob.timebase); // C + AC
    }else{
        return IM_FLUX_COMMAND_VALUE;   
    }
}

#endif

