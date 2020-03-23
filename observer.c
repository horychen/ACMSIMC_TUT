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
        sm.is_lpf[i]  = 0;
        sm.is_hpf[i]  = 0;
        sm.is_bpf[i]  = 0;

        sm.current_lpf_register[i] = 0;
        sm.current_hpf_register[i] = 0;
        sm.current_bpf_register1[i] = 0;
        sm.current_bpf_register2[i] = 0;
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
    ob.Ld_inv = 1.0/ob.Ld;
    ob.Lq = PMSM_Q_AXIS_INDUCTANCE;
    ob.KE = PMSM_PERMANENT_MAGNET_FLUX_LINKAGE; // Vs/rad

    ob.DeltaL = ob.Ld - ob.Lq;

    ob.Js = PMSM_SHAFT_INERTIA;
    ob.Js_inv = 1.0/ob.Js;

    ob.omg_elec = 0.0;
    ob.omg_mech = ob.omg_elec * sm.npp_inv;
    ob.theta_d = 0.0;

    ob.eemf_al = 0.0;
    ob.eemf_be = 0.0;


    ob.xPsi[0] = 0.0;
    ob.xPsi[1] = 0.0;
    ob.xChi[0] = 0.0;
    ob.xChi[1] = 0.0;

    ob.xEta[0] = 0.0;
    ob.xEta[1] = 0.0;
    ob.xVarSigma[0] = 0.0;
    ob.xVarSigma[1] = 0.0;

    ob.xUpsilon[0] = 0.0;
    ob.xUpsilon[1] = 0.0;
    ob.xZeta[0] = 0.0;
    ob.xZeta[1] = 0.0;

    ob.xOmg = 0.0;



    ob.output_error[0] = 0.0;
    ob.output_error[1] = 0.0;
    ob.output_error_eff[0] = 0.0;
    ob.output_error_eff[1] = 0.0;

    ob.k1          = OB_COEF_K1;
    ob.k2          = OB_COEF_K2;
    ob.gamma_omega = OB_COEF_GAMMA;
}

#define NUMBER_OF_EDA_STATES 13
void rhs_func_eda(double *xxn, double *xPsi, double *xChi, 
                               double *xEta, double *xVarSigma,
                               double *xUpsilon, double *xZeta,
                               double xOmg, double hs){

    #define R           ob.R
    #define LD          ob.Ld
    #define LD_INV      ob.Ld_inv
    #define LQ          ob.Lq
    #define DeltaL      ob.DeltaL
    #define k1          ob.k1    // Diagonal gain for eemf observer
    #define k2          ob.k2    // Off-diagonal gain for eemf observer (=0)
    #define gamma_omega ob.gamma_omega // Gain for speed update rule
    // #define OMG_USED    sm.omg_elec
    #define OMG_USED xOmg



    // Compute mismatch
    ob.output_error[0] = LD*IS(0) - xPsi[0];
    ob.output_error[1] = LD*IS(1) - xPsi[1];
    ob.output_error_eff[0] = ob.output_error[0] + xEta[0] - xUpsilon[0]*OMG_USED;
    ob.output_error_eff[1] = ob.output_error[1] + xEta[1] - xUpsilon[1]*OMG_USED;

    double f[NUMBER_OF_EDA_STATES];

    // xPsi: d-axis leakage flux linkage
    f[0] = US(0) - R*IS(0) + xChi[0] + (DeltaL+LD)*OMG_USED*-IS(1) + k1*ob.output_error[0];
    f[1] = US(1) - R*IS(1) + xChi[1] + (DeltaL+LD)*OMG_USED* IS(0) + k1*ob.output_error[1];
    // xChi: New emf state
    f[2] = xOmg*xOmg*DeltaL*IS(0) - OMG_USED*-US(1) + OMG_USED*R*-IS(1) + k2*ob.output_error[0];
    f[3] = xOmg*xOmg*DeltaL*IS(1) - OMG_USED* US(0) + OMG_USED*R* IS(0) + k2*ob.output_error[1];

    // xEta: auxiliary state for effective output error
    f[4] = -k1*xEta[0] + xVarSigma[0] + OMG_USED*(DeltaL+LD)*-IS(1);
    f[5] = -k1*xEta[1] + xVarSigma[1] + OMG_USED*(DeltaL+LD)* IS(0);
    // xVarSigma: auxiliary state for xEta (output error like signal)
    f[6] = -k2*xEta[0]                + OMG_USED*(2*OMG_USED*DeltaL*IS(0) + US(1) + R*-IS(1));
    f[7] = -k2*xEta[1]                + OMG_USED*(2*OMG_USED*DeltaL*IS(1) - US(0) + R* IS(0));

    // xUpsilon: filtered regressor for effective output error
    f[8]  = -k1*xUpsilon[0] + xZeta[0] + (DeltaL+LD)*-IS(1);
    f[9]  = -k1*xUpsilon[1] + xZeta[1] + (DeltaL+LD)* IS(0);
    // xZeta: auxiliary state for xUpsilon (output error like signal)
    f[10] = -k2*xUpsilon[0]            + 2*OMG_USED*DeltaL*IS(0) + US(1) + R*-IS(1);
    f[11] = -k2*xUpsilon[1]            + 2*OMG_USED*DeltaL*IS(1) - US(0) + R* IS(0);

    // xOmg: Speed update rule
    f[12] = gamma_omega * (xUpsilon[0]*ob.output_error_eff[0] + xUpsilon[1]*ob.output_error_eff[1]);

    // convert derivative to incremental for current step
    xxn[0]  = ( f[0] )*hs;
    xxn[1]  = ( f[1] )*hs;
    xxn[2]  = ( f[2] )*hs;
    xxn[3]  = ( f[3] )*hs;
    xxn[4]  = ( f[4] )*hs;
    xxn[5]  = ( f[5] )*hs;
    xxn[6]  = ( f[6] )*hs;
    xxn[7]  = ( f[7] )*hs;
    xxn[8]  = ( f[8] )*hs;
    xxn[9]  = ( f[9] )*hs;
    xxn[10] = ( f[10] )*hs;
    xxn[11] = ( f[11] )*hs;
    xxn[12] = ( f[12] )*hs;

    #undef R
    #undef LD    
    #undef LD_INV
    #undef LQ
    #undef k1
    #undef k2
    #undef gamma_omega
    #undef OMG_USED
}
void rk4_eda(double hs){
    static double xx1[NUMBER_OF_EDA_STATES];
    static double xx2[NUMBER_OF_EDA_STATES];
    static double xx3[NUMBER_OF_EDA_STATES];
    static double xx4[NUMBER_OF_EDA_STATES];
    static double x_temp[NUMBER_OF_EDA_STATES];
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
    rhs_func_eda( xx1, ob.xPsi, ob.xChi, 
                        ob.xEta, ob.xVarSigma,
                        ob.xUpsilon, ob.xZeta,
                        ob.xOmg, hs ); 
    x_temp[0]   = ob.xPsi[0]      + xx1[0] *0.5;
    x_temp[1]   = ob.xPsi[1]      + xx1[1] *0.5;
    x_temp[2]   = ob.xChi[0]      + xx1[2] *0.5;
    x_temp[3]   = ob.xChi[1]      + xx1[3] *0.5;
    x_temp[4]   = ob.xEta[0]      + xx1[4] *0.5;
    x_temp[5]   = ob.xEta[1]      + xx1[5] *0.5;
    x_temp[6]   = ob.xVarSigma[0] + xx1[6] *0.5;
    x_temp[7]   = ob.xVarSigma[1] + xx1[7] *0.5;
    x_temp[8]   = ob.xUpsilon[0]  + xx1[8] *0.5;
    x_temp[9]   = ob.xUpsilon[1]  + xx1[9] *0.5;
    x_temp[10]  = ob.xZeta[0]     + xx1[10]*0.5;
    x_temp[11]  = ob.xZeta[1]     + xx1[11]*0.5;
    x_temp[12]  = ob.xOmg         + xx1[12]*0.5;

    // time instant t+hs/2
    IS(0) = 0.5*(IS_P(0)+IS_C(0));
    IS(1) = 0.5*(IS_P(1)+IS_C(1));
    rhs_func_eda( xx2, p_x_temp, p_x_temp+2, 
                        p_x_temp+4, p_x_temp+6,
                        p_x_temp+8, p_x_temp+10,
                        *(p_x_temp+12), hs );
    x_temp[0]   = ob.xPsi[0]      + xx2[0] *0.5;
    x_temp[1]   = ob.xPsi[1]      + xx2[1] *0.5;
    x_temp[2]   = ob.xChi[0]      + xx2[2] *0.5;
    x_temp[3]   = ob.xChi[1]      + xx2[3] *0.5;
    x_temp[4]   = ob.xEta[0]      + xx2[4] *0.5;
    x_temp[5]   = ob.xEta[1]      + xx2[5] *0.5;
    x_temp[6]   = ob.xVarSigma[0] + xx2[6] *0.5;
    x_temp[7]   = ob.xVarSigma[1] + xx2[7] *0.5;
    x_temp[8]   = ob.xUpsilon[0]  + xx2[8] *0.5;
    x_temp[9]   = ob.xUpsilon[1]  + xx2[9] *0.5;
    x_temp[10]  = ob.xZeta[0]     + xx2[10]*0.5;
    x_temp[11]  = ob.xZeta[1]     + xx2[11]*0.5;
    x_temp[12]  = ob.xOmg         + xx2[12]*0.5;

    // time instant t+hs/2
    rhs_func_eda( xx3, p_x_temp, p_x_temp+2, 
                        p_x_temp+4, p_x_temp+6,
                        p_x_temp+8, p_x_temp+10,
                        *(p_x_temp+12), hs );
    x_temp[0]   = ob.xPsi[0]      + xx3[0];
    x_temp[1]   = ob.xPsi[1]      + xx3[1];
    x_temp[2]   = ob.xChi[0]      + xx3[2];
    x_temp[3]   = ob.xChi[1]      + xx3[3];
    x_temp[4]   = ob.xEta[0]      + xx3[4];
    x_temp[5]   = ob.xEta[1]      + xx3[5];
    x_temp[6]   = ob.xVarSigma[0] + xx3[6];
    x_temp[7]   = ob.xVarSigma[1] + xx3[7];
    x_temp[8]   = ob.xUpsilon[0]  + xx3[8];
    x_temp[9]   = ob.xUpsilon[1]  + xx3[9];
    x_temp[10]  = ob.xZeta[0]     + xx3[10];
    x_temp[11]  = ob.xZeta[1]     + xx3[11];
    x_temp[12]  = ob.xOmg         + xx3[12];

    // time instant t+hs
    IS(0) = IS_C(0);
    IS(1) = IS_C(1);
    rhs_func_eda( xx4, p_x_temp, p_x_temp+2, 
                        p_x_temp+4, p_x_temp+6,
                        p_x_temp+8, p_x_temp+10,
                        *(p_x_temp+12), hs );
    // \+=[^\n]*1\[(\d+)\][^\n]*2\[(\d+)\][^\n]*3\[(\d+)\][^\n]*4\[(\d+)\][^\n]*/ ([\d]+)
    // +=   (xx1[$5] + 2*(xx2[$5] + xx3[$5]) + xx4[$5])*0.166666666666667; // $5
    ob.xPsi[0]      +=   (xx1[0] + 2*(xx2[0] + xx3[0]) + xx4[0])*0.166666666666667; // 0
    ob.xPsi[1]      +=   (xx1[1] + 2*(xx2[1] + xx3[1]) + xx4[1])*0.166666666666667; // 1
    ob.xChi[0]      +=   (xx1[2] + 2*(xx2[2] + xx3[2]) + xx4[2])*0.166666666666667; // 2
    ob.xChi[1]      +=   (xx1[3] + 2*(xx2[3] + xx3[3]) + xx4[3])*0.166666666666667; // 3
    ob.xEta[0]      +=   (xx1[4] + 2*(xx2[4] + xx3[4]) + xx4[4])*0.166666666666667; // 4
    ob.xEta[1]      +=   (xx1[5] + 2*(xx2[5] + xx3[5]) + xx4[5])*0.166666666666667; // 5
    ob.xVarSigma[0] +=   (xx1[6] + 2*(xx2[6] + xx3[6]) + xx4[6])*0.166666666666667; // 6
    ob.xVarSigma[1] +=   (xx1[7] + 2*(xx2[7] + xx3[7]) + xx4[7])*0.166666666666667; // 7
    ob.xUpsilon[0]  +=   (xx1[8] + 2*(xx2[8] + xx3[8]) + xx4[8])*0.166666666666667; // 8
    ob.xUpsilon[1]  +=   (xx1[9] + 2*(xx2[9] + xx3[9]) + xx4[9])*0.166666666666667; // 9
    ob.xZeta[0]     +=   (xx1[10] + 2*(xx2[10] + xx3[10]) + xx4[10])*0.166666666666667; // 10
    ob.xZeta[1]     +=   (xx1[11] + 2*(xx2[11] + xx3[11]) + xx4[11])*0.166666666666667; // 11
    ob.xOmg         +=   (xx1[12] + 2*(xx2[12] + xx3[12]) + xx4[12])*0.166666666666667; // 12

    // EEMF
    ob.eemf_al  = -ob.Ld * ob.xOmg *-IS(1) - ob.xChi[0];
    ob.eemf_be  = -ob.Ld * ob.xOmg * IS(0) - ob.xChi[1];
    // ob.theta_d  = atan2(-ob.eemf_al,  
    //                      ob.eemf_be); // 180 deg shift in rotor position when speed is negative.
    if(ob.xOmg == 0){ // 当xOmg恒为零时，用转速的指令的正负校正反转180度转子位置角度差
        ob.theta_d  = atan2(-ob.eemf_al*sign(ACM.rpm_cmd), 
                             ob.eemf_be*sign(ACM.rpm_cmd));
    }else{
        ob.theta_d  = atan2(-ob.eemf_al*sign(ob.xOmg), 
                             ob.eemf_be*sign(ob.xOmg)); // It seems very important to use sign(ob.xOmg) instead of sign(ACM.rpm_cmd) or else starting-up may be oscillating.
    }
    ob.omg_elec = ob.xOmg;
    ob.omg_mech = ob.omg_elec * sm.npp_inv;
}

void observation(){

    // ob.g1 = OB_COEF_K1;
    // if(CTRL.timebase>2){
    //     ob.g1 = -fabs(OB_OMG)*sm.Ld*2;
    // }
    // OB_OMG = sm.omg_elec;

    ob.k2 = (4*20000) * ACM.rpm_cmd / 1500;

    /* OBSERVATION */
    rk4_eda(TS); 

    /* 备份这个采样点的数据供下次使用。所以，观测的和实际的相比，是延迟一个采样周期的。 */
    //  2017年1月20日，将控制器放到了观测器的后面。
    // * 所以，上一步电压US_P的更新也要延后了。
    // US_P(0) = US_C(0); 
    // US_P(1) = US_C(1);
    IS_P(0) = IS_C(0);
    IS_P(1) = IS_C(1);
}


// ------------------------------------------------- HFSI


#ifdef HFSI_ON
    #define LPF_TIME_CONST_INVERSE (5*2*M_PI) // time constant is 1/400 <=> cutoff frequency is 400/(2*pi) ||| 换句话说，截止频率 * 2pi = 时间常数的倒数 |||| f_cutoff = 1 / (time constant * 2*pi)
    #define LUENBERGER_GAIN_1 30     // 30    // 30    // 20  // Large gain to position will cause steady state position error, but increase it close to limit
    #define LUENBERGER_GAIN_2 (750)  // (300) // (300) // 100 // If speed estimate has too much dynamics during reversal, you need to increase this gain actually...
    #define LUENBERGER_GAIN_3 (6000) // (1500) // (790) // 500 // Tune reversal response to slight over-shoot

    void dynamics_lpf(double input, double *state, double *derivative){
        derivative[0] = LPF_TIME_CONST_INVERSE * ( input - *state );
    }
    void RK4_111_general(void (*pointer_dynamics)(), double input, double *state, double hs){
        // 我把euler 改成rk4以后就一切正常了
        // 王彤:
        // 其实rk4也不是正解，因为运算量太大，一般在dsp里不用的
        // 按理应该用冲剂响应不变法离散
        // 或者用prewarp tustin离散
        // 冲击响应不变法类似与利用书里那个s域和z域的变换表变换，查查那个表就完了
        // 但是我想你手头就有现成的rk4
        // 所以应该更方便
        // 不过冲剂响应不变法在通带的衰减不是0db
        // 稍微差一点

        #define NS 1

        double k1[NS], k2[NS], k3[NS], k4[NS], intemediate_state[NS];
        double derivative[NS];
        int i;

        pointer_dynamics(input, state, derivative); // timer.t,
        for(i=0;i<NS;++i){        
            k1[i] = derivative[i] * hs;
            intemediate_state[i] = state[i] + k1[i]*0.5;
        }

        pointer_dynamics(input, intemediate_state, derivative); // timer.t+hs/2., 
        for(i=0;i<NS;++i){        
            k2[i] = derivative[i] * hs;
            intemediate_state[i] = state[i] + k2[i]*0.5;
        }
        
        pointer_dynamics(input, intemediate_state, derivative); // timer.t+hs/2., 
        for(i=0;i<NS;++i){        
            k3[i] = derivative[i] * hs;
            intemediate_state[i] = state[i] + k3[i];
        }
        
        pointer_dynamics(input, intemediate_state, derivative); // timer.t+hs, 
        for(i=0;i<NS;++i){        
            k4[i] = derivative[i] * hs;
            state[i] = state[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0;
        }
        #undef NS
    }
    void RK4_333_general(void (*pointer_dynamics)(), double input, double *state, double hs){
        #define NS 3

        double k1[NS], k2[NS], k3[NS], k4[NS], intemediate_state[NS];
        double derivative[NS];
        int i;

        pointer_dynamics(input, state, derivative); // timer.t,
        for(i=0;i<NS;++i){        
            k1[i] = derivative[i] * hs;
            intemediate_state[i] = state[i] + k1[i]*0.5;
        }

        pointer_dynamics(input, intemediate_state, derivative); // timer.t+hs/2., 
        for(i=0;i<NS;++i){        
            k2[i] = derivative[i] * hs;
            intemediate_state[i] = state[i] + k2[i]*0.5;
        }
        
        pointer_dynamics(input, intemediate_state, derivative); // timer.t+hs/2., 
        for(i=0;i<NS;++i){        
            k3[i] = derivative[i] * hs;
            intemediate_state[i] = state[i] + k3[i];
        }
        
        pointer_dynamics(input, intemediate_state, derivative); // timer.t+hs, 
        for(i=0;i<NS;++i){        
            k4[i] = derivative[i] * hs;
            state[i] = state[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0;
        }
        #undef NS
    }
    struct HFSI_Data hfsi;
    void hfsi_init(){
        hfsi.test_signal_al = 0.0;
        hfsi.test_signal_be = 0.0;
        hfsi.test_signal_M = 0.0;
        hfsi.test_signal_T = 0.0;
        hfsi.M_lpf = 0.0;
        hfsi.T_lpf = 0.0;
        hfsi.M_hpf = 0.0;
        hfsi.T_hpf = 0.0;
        hfsi.theta_filter = 0.0;
        hfsi.theta_d_raw = 0.0;
        hfsi.theta_d = 0.0;
        hfsi.omg_elec = 0.0;
        hfsi.pseudo_load_torque = 0.0;
        hfsi.mismatch = 0.0;
    }
    double difference_between_two_angles(double first, double second){
        while(first>2*M_PI){
            first-=2*M_PI;
        }
        while(second>2*M_PI){
            second-=2*M_PI;
        }

        while(first<0.0){
            first+=2*M_PI;
        }
        while(second<0.0){
            second+=2*M_PI;
        }

        if(fabs(first-second)<M_PI){
            return first-second;
        }else{
            if(first>second){
                return first-2*M_PI-second;
            }else{                
                return first+2*M_PI-second;
            }
        }
    }
    void dynamics_position_observer(double input, double *state, double *derivative){
        #define TEM_USED CTRL.Tem
        // #define TEM_USED ACM.Tem
        hfsi.mismatch = difference_between_two_angles(input, state[0]); // difference in angle
        derivative[0] = LUENBERGER_GAIN_1*hfsi.mismatch + state[1];
        derivative[1] = LUENBERGER_GAIN_2*hfsi.mismatch + TEM_USED*CTRL.npp/CTRL.Js - state[2];
        derivative[2] = -LUENBERGER_GAIN_3*hfsi.mismatch;
        // printf("%g, %g, %g\n", CTRL.timebase, state[2], derivative[2]);
    }
    void luenberger_filter(double theta_d_raw){
        static double state[3];

        RK4_333_general(dynamics_position_observer, theta_d_raw, state, TS);

        if(state[0]>M_PI){
            state[0] -= 2*M_PI;
        }else if(state[0]<-M_PI){
            state[0] += 2*M_PI;
        }
        hfsi.theta_d            = state[0] - 0.5*M_PI;
        hfsi.omg_elec           = state[1];
        hfsi.pseudo_load_torque = state[2];
    }
    void hfsi_do(){
        // luenberger_filter() needs CTRL.Tem so hfsi_do() should executed after contorl().
        // If filtered currents are used instead of IS_C, move everything before luenberger_filter() to measurement(). <- Not suggested: lpf ruins current loop control.

        hfsi.test_signal_al = ACM.ial;
        hfsi.test_signal_be = ACM.ibe;
        {
            // hfsi.theta_filter = sm.theta_d;
            hfsi.theta_filter = hfsi.theta_d;

            hfsi.test_signal_M = AB2M(hfsi.test_signal_al, hfsi.test_signal_be, cos(hfsi.theta_filter), sin(hfsi.theta_filter));
            hfsi.test_signal_T = AB2T(hfsi.test_signal_al, hfsi.test_signal_be, cos(hfsi.theta_filter), sin(hfsi.theta_filter));

            // LPF
            RK4_111_general(dynamics_lpf, hfsi.test_signal_M, &hfsi.M_lpf, TS);
            RK4_111_general(dynamics_lpf, hfsi.test_signal_T, &hfsi.T_lpf, TS);
            IS_LPF(0) = MT2A(hfsi.M_lpf, hfsi.T_lpf, cos(hfsi.theta_filter), sin(hfsi.theta_filter));
            IS_LPF(1) = MT2B(hfsi.M_lpf, hfsi.T_lpf, cos(hfsi.theta_filter), sin(hfsi.theta_filter));

            // HPF
            double LAST_IS_HPF[2];
            double DELTA_IS_HPF[2];
            LAST_IS_HPF[0] = IS_HPF(0);
            LAST_IS_HPF[1] = IS_HPF(1);
            hfsi.M_hpf = hfsi.test_signal_M - hfsi.M_lpf;
            hfsi.T_hpf = hfsi.test_signal_T - hfsi.T_lpf;
            IS_HPF(0) = MT2A(hfsi.M_hpf, hfsi.T_hpf, cos(hfsi.theta_filter), sin(hfsi.theta_filter));
            IS_HPF(1) = MT2B(hfsi.M_hpf, hfsi.T_hpf, cos(hfsi.theta_filter), sin(hfsi.theta_filter));

            DELTA_IS_HPF[0] = IS_HPF(0) - LAST_IS_HPF[0];
            DELTA_IS_HPF[1] = IS_HPF(1) - LAST_IS_HPF[1];
            hfsi.theta_d_raw = atan2(DELTA_IS_HPF[1], DELTA_IS_HPF[0]);
            if(hfsi.theta_d_raw>M_PI){
                printf("%g", hfsi.theta_d_raw/M_PI*180);
                hfsi.theta_d_raw -= 2*M_PI;
            }else if(hfsi.theta_d_raw<-M_PI){
                printf("%g", hfsi.theta_d_raw/M_PI*180);
                hfsi.theta_d_raw += 2*M_PI;
            }
        }
        // IS_C(0) = IS_LPF(0); // The waveform of filtered currents looks good but it is still delayed and that is detrimental to control system stability.
        // IS_C(1) = IS_LPF(1); 

        luenberger_filter(hfsi.theta_d_raw);
    }
#endif
