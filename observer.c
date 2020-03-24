#include "ACMSim.h"


struct SynchronousMachine sm;
// struct Observer ob;

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


// ------------------------------------------------- HFSI


#ifdef HFSI_ON
    #define LPF_TIME_CONST_INVERSE (5*2*M_PI) // time constant is 1/400 <=> cutoff frequency is 400/(2*pi) ||| 换句话说，截止频率 * 2pi = 时间常数的倒数 |||| f_cutoff = 1 / (time constant * 2*pi)
    #define LUENBERGER_GAIN_1 30       // 30       // 30     // 30    // 20  // Large gain to position will cause steady state position error, but increase it close to limit
    #define LUENBERGER_GAIN_2 (900)    // (750)    // (300)  // (300) // 100 // If speed estimate has too much dynamics during reversal, you need to increase this gain actually...
    #define LUENBERGER_GAIN_3 (1500) // (0*6000) // (1500) // (790) // 500 // Tune reversal response to slight over-shoot

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
