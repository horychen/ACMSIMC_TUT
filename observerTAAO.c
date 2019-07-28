#include "ACMSim.h" 
 
#if OBSERVER_APPLIED == TOTALLY_ADAPTIVE_OBSERVER 
 
struct InductionMachine im; 
struct ObserverControl ob; 
 
/* Total Adaptation 
 * */ 
double deriv_rs; 
double deriv_Lsigma; 
double deriv_rreq; 
double deriv_Leq; 
double deriv_omg; 
double deriv_invJs; 
double deriv_TLslashJs; 
void rhs_func_14th(double *xxn, double *xVM, double *xCM, double *xB, 
    double xRs, double xLsigma, double xRreq, double xLeq, double xOmg,  
    double xVarPi, double xInvJs, double xTLslashJs, double hs){ 
 
    // Compute mismatch 
    ob.xVM[0] = (xVM[0] - xLsigma * IS(0)); // VM rotor flux = VM stator flux - Lsigma*i_s 
    ob.xVM[1] = (xVM[1] - xLsigma * IS(1)); 
    ob.mismatch[0] = ob.xVM[0] - xCM[0]; 
    ob.mismatch[1] = ob.xVM[1] - xCM[1]; 
    #ifdef PC_SIMULATION 
        ob.error[0] = IM.x[2]*IM.Lm_slash_Lr - xCM[0]; //(STEADY_FLUX_LINKAGE_MODULUS + ob.extraModulus)*ob.cosT - xCM[0]; 
        ob.error[1] = IM.x[3]*IM.Lm_slash_Lr - xCM[1]; //(STEADY_FLUX_LINKAGE_MODULUS + ob.extraModulus)*ob.sinT - xCM[1]; 
    #endif 
 
    /* Voltage model for stator flux */ 
    xxn[0] = ( US(0) - xRs*IS(0) - ob.k_VM*ob.mismatch[0] )*hs;
    xxn[1] = ( US(1) - xRs*IS(1) - ob.k_VM*ob.mismatch[1] )*hs; 
 
    /* Current model */ 
    double xLeq_inv = 1.0 / xLeq; 
    double xAlpha = xRreq*xLeq_inv; 
    xxn[2] = ( -xAlpha*xCM[0] + xRreq*IS(0) + xOmg*-xCM[1] + ob.k_CM*ob.mismatch[0] )*hs; 
    xxn[3] = ( -xAlpha*xCM[1] + xRreq*IS(1) + xOmg* xCM[0] + ob.k_CM*ob.mismatch[1] )*hs; 
 
    /* hat _xB derivative 
     * */ 
    xxn[4] = 0.0;
    xxn[5] = 0.0; 
 
    /* hat _xRs derivative */  
    deriv_rs = ob.gamma1 * (ob.mismatch[0]*IS(0) + ob.mismatch[1]*IS(1)); 
    xxn[6] = deriv_rs*hs; 
 
    /* hat Lsigma derivative*/ 
    xxn[7] = 0; 
 
    /* hat _xRreq derivative */ 
    deriv_rreq = ob.gamma3 * (ob.mismatch[0]*(IS(0)-xCM[0]*xLeq_inv) + ob.mismatch[1]*(IS(1)-xCM[1]*xLeq_inv)); 
    xxn[8] = deriv_rreq*hs; 
 
    /* hat _xLeq derivative */ 
    // deriv_Leq = ob.gamma4 * (ob.mismatch[0]*(xAlpha*xLeq_inv*xCM[0]) + ob.mismatch[1]*(xAlpha*xLeq_inv*xCM[1])); 
    deriv_Leq = ob.gamma4 * (ob.mismatch[0]*xCM[0] + ob.mismatch[1]*xCM[1]); // coefficient absorbsion 
    xxn[9] = deriv_Leq*hs; 
 
    /* hat _xOmg derivative */ 
    deriv_omg = ob.gamma5 * (ob.mismatch[0]*-xCM[1] + ob.mismatch[1]*xCM[0]);     
    xxn[10] = deriv_omg*hs; 
} 
void rK_fourteenth(double hs){ 
    double xx1[14]; 
    double xx2[14]; 
    double xx3[14]; 
    double xx4[14]; 
    double x_temp[14]; 
    double *p_x_temp=x_temp; 
 
    /* Theoritically speaking, rhs_func should be time-varing like rhs_func(.,t). 
       To apply codes in DSP, we do time-varing updating of IS(0) and IS(1) outside rhs_func(.) to save time. */ 
 
    // time instant t 
    US(0) = US_P(0); 
    US(1) = US_P(1); 
    IS(0) = IS_P(0); 
    IS(1) = IS_P(1); 
    rhs_func_14th( xx1, ob.statorFlux, ob.xCM, ob.taao_b, 
            ob.taao_rs, ob.taao_Lsigma, ob.taao_rreq, ob.taao_Leq, ob.taao_omg_integralPart,  
            ob.varPi, ob.taao_invJs, ob.taao_TLslashJs, hs); 
    x_temp[0]  = ob.statorFlux[0]  +   xx1[0]*0.5; 
    x_temp[1]  = ob.statorFlux[1]  +   xx1[1]*0.5; 
    x_temp[2]  = ob.xCM[0]         +   xx1[2]*0.5; 
    x_temp[3]  = ob.xCM[1]         +   xx1[3]*0.5; 
    x_temp[6]  = ob.taao_rs        +   xx1[6]*0.5; 
    x_temp[7]  = ob.taao_Lsigma    +   xx1[7]*0.5; 
    x_temp[8]  = ob.taao_rreq      +   xx1[8]*0.5; 
    x_temp[9]  = ob.taao_Leq       +   xx1[9]*0.5; 
    x_temp[10] = ob.taao_omg_integralPart + xx1[10]*0.5; 
    x_temp[11] = ob.varPi          +   xx1[11]*0.5; 
    x_temp[12] = ob.taao_invJs     +   xx1[12]*0.5; 
    x_temp[13] = ob.taao_TLslashJs +   xx1[13]*0.5; 
 
    // time instant t+hs/2 
    IS(0) = 0.5*(IS_P(0)+IS_C(0)); 
    IS(1) = 0.5*(IS_P(1)+IS_C(1)); 
    rhs_func_14th( xx2, p_x_temp, p_x_temp+2, p_x_temp+4, 
            *(p_x_temp+6), *(p_x_temp+7), *(p_x_temp+8), *(p_x_temp+9), *(p_x_temp+10),  
            *(p_x_temp+11), *(p_x_temp+12), *(p_x_temp+13), hs); 
    x_temp[0]  = ob.statorFlux[0]  +   xx2[0]*0.5; 
    x_temp[1]  = ob.statorFlux[1]  +   xx2[1]*0.5; 
    x_temp[2]  = ob.xCM[0]         +   xx2[2]*0.5; 
    x_temp[3]  = ob.xCM[1]         +   xx2[3]*0.5; 
    x_temp[6]  = ob.taao_rs        +   xx2[6]*0.5; 
    x_temp[7]  = ob.taao_Lsigma    +   xx2[7]*0.5; 
    x_temp[8]  = ob.taao_rreq      +   xx2[8]*0.5; 
    x_temp[9]  = ob.taao_Leq       +   xx2[9]*0.5; 
    x_temp[10] = ob.taao_omg_integralPart + xx2[10]*0.5; 
    x_temp[11] = ob.varPi          +   xx2[11]*0.5; 
    x_temp[12] = ob.taao_invJs     +   xx2[12]*0.5; 
    x_temp[13] = ob.taao_TLslashJs +   xx2[13]*0.5; 
 
    // time instant t+hs/2 
    rhs_func_14th( xx3, p_x_temp, p_x_temp+2, p_x_temp+4, 
            *(p_x_temp+6), *(p_x_temp+7), *(p_x_temp+8), *(p_x_temp+9), *(p_x_temp+10),  
            *(p_x_temp+11), *(p_x_temp+12), *(p_x_temp+13), hs); 
    x_temp[0]  = ob.statorFlux[0]  +   xx3[0]; 
    x_temp[1]  = ob.statorFlux[1]  +   xx3[1]; 
    x_temp[2]  = ob.xCM[0]         +   xx3[2]; 
    x_temp[3]  = ob.xCM[1]         +   xx3[3]; 
    x_temp[6]  = ob.taao_rs        +   xx3[6]; 
    x_temp[7]  = ob.taao_Lsigma    +   xx3[7]; 
    x_temp[8]  = ob.taao_rreq      +   xx3[8]; 
    x_temp[9]  = ob.taao_Leq       +   xx3[9]; 
    x_temp[10] = ob.taao_omg_integralPart + xx3[10]; 
    x_temp[11] = ob.varPi          +   xx3[11]; 
    x_temp[12] = ob.taao_invJs     +   xx3[12]; 
    x_temp[13] = ob.taao_TLslashJs +   xx3[13]; 
 
    // time instant t+hs 
    IS(0) = IS_C(0); 
    IS(1) = IS_C(1); 
    rhs_func_14th( xx4, p_x_temp, p_x_temp+2, p_x_temp+4, 
            *(p_x_temp+6), *(p_x_temp+7), *(p_x_temp+8), *(p_x_temp+9), *(p_x_temp+10),  
            *(p_x_temp+11), *(p_x_temp+12), *(p_x_temp+13), hs); 
    // \+=[^\n]*1\[(\d)\][^\n]*2\[(\d)\][^\n]*3\[(\d)\][^\n]*4\[(\d)\][^\n]*/ ([\d]+) 
    // +=   (xx1[$5] + 2*(xx2[$5] + xx3[$5]) + xx4[$5])*0.166666666666667; // $5 
    ob.statorFlux[0] +=   (xx1[0] + 2*(xx2[0] + xx3[0]) + xx4[0])*0.166666666666667; // 0 
    ob.statorFlux[1] +=   (xx1[1] + 2*(xx2[1] + xx3[1]) + xx4[1])*0.166666666666667; // 1 
    ob.xCM[0]        +=   (xx1[2] + 2*(xx2[2] + xx3[2]) + xx4[2])*0.166666666666667; // 2 
    ob.xCM[1]        +=   (xx1[3] + 2*(xx2[3] + xx3[3]) + xx4[3])*0.166666666666667; // 3 
    ob.taao_rs       +=   (xx1[6] + 2*(xx2[6] + xx3[6]) + xx4[6])*0.166666666666667; // 6 
    ob.taao_Lsigma   +=   (xx1[7] + 2*(xx2[7] + xx3[7]) + xx4[7])*0.166666666666667; // 7 
    ob.taao_rreq     +=   (xx1[8] + 2*(xx2[8] + xx3[8]) + xx4[8])*0.166666666666667; // 8 
    ob.taao_Leq      +=   (xx1[9] + 2*(xx2[9] + xx3[9]) + xx4[9])*0.166666666666667; // 9 
    ob.taao_omg_integralPart +=   (xx1[10] + 2*(xx2[10] + xx3[10]) + xx4[10])*0.166666666666667; // 10 
    ob.varPi                 +=   (xx1[11] + 2*(xx2[11] + xx3[11]) + xx4[11])*0.166666666666667; // 11 
    ob.taao_invJs            +=   (xx1[12] + 2*(xx2[12] + xx3[12]) + xx4[12])*0.166666666666667; // 12 
    ob.taao_TLslashJs        +=   (xx1[13] + 2*(xx2[13] + xx3[13]) + xx4[13])*0.166666666666667; // 13 
 
    // 限幅 
    if(ob.taao_Leq<0.2){ 
        ob.taao_Leq = 0.2; // Replace it with Proj() by Marino 
    } 
 
    ob.taao_alpha = ob.taao_rreq/ob.taao_Leq; 
 
    ob.taao_omg = ob.taao_omg_integralPart + ob.gamma5P * (ob.mismatch[0]*-ob.xCM[1] + ob.mismatch[1]*ob.xCM[0]); 
    ob.taao_speed = ob.taao_omg * RAD_PER_SEC_2_RPM;    // 60/3 
 
    ob.xVM[0] = ob.statorFlux[0] - ob.taao_Lsigma * IS(0); 
    ob.xVM[1] = ob.statorFlux[1] - ob.taao_Lsigma * IS(1); 
} 
 
/* Main Flow of Observer 
 * */ 
void observation(){ 

    // if(ob.timebase>10 && ob.timebase<12.5){
    //     // Regeneration mode
    //     ob.k_VM = 10;
    // }else{
    //     // Motoring operation
    //     ob.k_VM = KVM_VALUE;
    // }

    /* OBSERVATION */ 
    if(ob.gamma5==0){ 
        ob.taao_omg_integralPart = im.omg; 
    } 
    rK_fourteenth(ob.Ts); // 放在最后，保证im.omg已经更新*（NON_ADAPTIVE case） 
 
    IS_P(0) = IS_C(0); 
    IS_P(1) = IS_C(1); 
 
    ob.mismatch_mod = sqrt(ob.mismatch[0]*ob.mismatch[0]+ob.mismatch[1]*ob.mismatch[1]); 
} 
 
 
/* Initialization Codes*/ 
void acm_init(){ 
    im.us[0] = 0.0; 
    im.us[1] = 0.0; 
    im.is[0] = 0.0; 
    im.is[1] = 0.0; 
    im.us_curr[0] = 0.0; 
    im.us_curr[1] = 0.0; 
    im.is_curr[0] = 0.0; 
    im.is_curr[1] = 0.0; 
    im.us_prev[0] = 0.0; 
    im.us_prev[1] = 0.0; 
    im.is_prev[0] = 0.0; 
    im.is_prev[1] = 0.0; 
 
    im.Js = 0.032; 
    im.Js_inv = 1.0/im.Js; 
 
    im.Leq = 0.4482; 
    im.Lls = 0.0126; 
    #if NO_ROTOR_LEAKAGE 
        im.Llr = 0.0; 
    #else 
        im.Llr = 0.0126; 
    #endif 
    im.Lm = 0.5*(im.Leq+sqrtf(im.Leq*im.Leq+4*im.Llr*im.Leq)); 
    im.Lm_inv = 1.0/im.Lm; 
    im.Ls = im.Lm + im.Lls; 
    im.Lr = im.Lm + im.Llr; 
    im.sigma = 1.0 - im.Lm*im.Lm/im.Ls/im.Lr; 
 
    im.rs = 3.04; 
    im.rr = 1.69; 
    im.Lsigma = im.Ls*im.sigma; 
    im.alpha = im.rr/im.Lr; 
    im.rreq = im.Leq * im.alpha; 
    im.Leq_inv = 1.0/im.Leq; 
    im.Tr = im.Lr/im.rr;  
 
    im.omg = 0.0; 
    im.omg_curr = 0.0; 
    im.omg_prev = 0.0; 
 
    im.npp = 2.0; // no. of pole pair    
 
    im.theta_r = 0.0; 
} 
void ob_init(){ 
 
    ob.k_VM = KVM_VALUE; 
    ob.k_CM = 0.0; 
 
    ob.xVM[0] = 0;  
    ob.xVM[1] = 0;  
    ob.xCM[0] = 0.0; 
    ob.xCM[1] = 0.0; 
    ob.statorFlux[0] = 0.0; //0.01; // initial error exsits 
    ob.statorFlux[1] = 0.0; //0.01; 
    ob.varPi = 0.0; 
    ob.Tem = 0.0; 
 
    ob.mismatch[0] = 0.0; 
    ob.mismatch[1] = 0.0; 
    ob.mismatch[2] = 0.0; 
    ob.mismatch_mod = 0.0; 
    ob.error[0] = 0.0; 
    ob.error[1] = 0.0; 
 
    ob.gamma1 = GAMMA_RS; 
    ob.gamma2 = GAMMA_LSIGMA; 
    ob.gamma3 = GAMMA_RREQ; 
    ob.gamma4 = GAMMA_LEQ; 
    ob.gamma5 = GAMMA_OMEGAR; 
    ob.gamma5P = GAMMA_OMEGAR_P; 
 
    ob.timebase = 0.0; 
    ob.Ts = TS; 
 
    ob.taao_rs = im.rs*1.0; 
    ob.taao_Lsigma = im.Lsigma*1.0; 
    ob.taao_rreq = im.rreq*1.0; 
    ob.taao_Leq = im.Leq*1.0; 
 
    ob.taao_alpha = ob.taao_rreq/ob.taao_Leq; 
 
    ob.taao_omg = 0.0;              // Actual value im.omg is updated in eQEP 
    ob.taao_omg_integralPart = 0.0; 
    ob.taao_speed = ob.taao_omg * RAD_PER_SEC_2_RPM;  
 
 
    ob.omega_e = 0.0; 
 
    ob.taao_flux_cmd = TAAO_FLUX_COMMAND_VALUE; 
    ob.taao_flux_cmd_on = TAAO_FLUX_COMMAND_ON; 
 
 } 
 
double TAAO_FluxModulusCommand(){ 
    if(ob.taao_flux_cmd_on){ 
        // return TAAO_FLUX_COMMAND_VALUE + M1*sin(OMG1*ob.timebase); // C + AC 
        return TAAO_FLUX_COMMAND_VALUE + M1*sin(OMG1*ob.timebase) + M2*sin(OMG2*ob.timebase) + M3*sin(OMG3*ob.timebase); // C + AC 
    }else{ 
        return TAAO_FLUX_COMMAND_VALUE;    
    } 
} 
 
/* Consistent with DSP Programs Add_TDDA.c Terminate */ 
 
#endif 
