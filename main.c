#include "ACMSim.h"

double sign(double x){
    return (x > 0) - (x < 0);    
}
double fabs(double x){
    return (x >= 0) ? x : -x;
}

struct SynchronousMachineSimulated ACM;
void Machine_init(){

    ACM.Ts = MACHINE_TS;

    int i;
    for(i=0;i<NUMBER_OF_STATES;++i){
        ACM.x[i] = 0.0;
        ACM.x_dot[i] = 0.0;
    }

    ACM.omg_elec = 0.0;
    ACM.rpm = 0.0;
    ACM.rpm_cmd = 0.0;
    ACM.rpm_deriv_cmd = 0.0;
    ACM.Tload = 0.0;
    ACM.Tem = 0.0;

    ACM.npp = PMSM_NUMBER_OF_POLE_PAIRS;

    ACM.R  = PMSM_RESISTANCE;
    ACM.Ld = PMSM_D_AXIS_INDUCTANCE;
    ACM.Lq = PMSM_Q_AXIS_INDUCTANCE;
    ACM.KE = PMSM_PERMANENT_MAGNET_FLUX_LINKAGE; // Vs/rad
    ACM.Js = PMSM_SHAFT_INERTIA;

    ACM.mu_m = ACM.npp/ACM.Js;
    ACM.L0 = 0.5*(ACM.Ld + ACM.Lq);
    ACM.L1 = 0.5*(ACM.Ld - ACM.Lq);

    ACM.ual = 0.0;
    ACM.ube = 0.0;
    ACM.ial = 0.0;
    ACM.ibe = 0.0;

    ACM.theta_d = 0.0;
    ACM.ud = 0.0;
    ACM.uq = 0.0;
    ACM.id = 0.0;
    ACM.iq = 0.0;

    ACM.eemf_q  = 0.0;
    ACM.eemf_al = 0.0;
    ACM.eemf_be = 0.0;
    ACM.theta_d__eemf = 0.0;
}

/* Simple Model */
void RK_dynamics(double t, double *x, double *fx){
    // electromagnetic model
    fx[0] = (ACM.ud - ACM.R * x[0] + x[2]*ACM.Lq*x[1]) / ACM.Ld; // current-d
    fx[1] = (ACM.uq - ACM.R * x[1] - x[2]*ACM.Ld*x[0] - x[2]*ACM.KE) / ACM.Lq; // current-q

    // mechanical model
    ACM.Tem = ACM.npp*(x[1]*ACM.KE + (ACM.Ld - ACM.Lq)*x[0]*x[1]);
    fx[2] = (ACM.Tem - ACM.Tload)*ACM.mu_m; // elec. angular rotor speed
    fx[3] = x[2];                           // elec. angular rotor position
}
void RK_Linear(double t, double *x, double hs){
    #define NS NUMBER_OF_STATES

    double k1[NS], k2[NS], k3[NS], k4[NS], xk[NS];
    double fx[NS];
    int i;

    RK_dynamics(t, x, fx); // timer.t,
    for(i=0;i<NS;++i){        
        k1[i] = fx[i] * hs;
        xk[i] = x[i] + k1[i]*0.5;
    }
    
    RK_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<NS;++i){        
        k2[i] = fx[i] * hs;
        xk[i] = x[i] + k2[i]*0.5;
    }
    
    RK_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<NS;++i){        
        k3[i] = fx[i] * hs;
        xk[i] = x[i] + k3[i];
    }
    
    RK_dynamics(t, xk, fx); // timer.t+hs, 
    for(i=0;i<NS;++i){        
        k4[i] = fx[i] * hs;
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0;

        // derivatives
        ACM.x_dot[i] = (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0 / hs; 
    }
    #undef NS
}


int machine_simulation(){

    // solve for ACM.x with ACM.ud and ACM.uq as inputs
    RK_Linear(CTRL.timebase, ACM.x, ACM.Ts);

    // rotor position
    ACM.theta_d = ACM.x[3];
    if(ACM.theta_d > M_PI){
        ACM.theta_d -= 2*M_PI;
    }else if(ACM.theta_d < -M_PI){
        ACM.theta_d += 2*M_PI; // 反转！
    }
    ACM.x[3] = ACM.theta_d;

    // currents
    ACM.id  = ACM.x[0];
    ACM.iq  = ACM.x[1];
    ACM.ial = MT2A(ACM.id, ACM.iq, cos(ACM.theta_d), sin(ACM.theta_d));
    ACM.ibe = MT2B(ACM.id, ACM.iq, cos(ACM.theta_d), sin(ACM.theta_d));

    // speed
    ACM.omg_elec = ACM.x[2];
    ACM.rpm = ACM.x[2] * 60 / (2 * M_PI * ACM.npp);

    // extended emf
    ACM.eemf_q  = (ACM.Ld-ACM.Lq) * (ACM.omg_elec*ACM.id - ACM.x_dot[1]) + ACM.omg_elec*ACM.KE;
    ACM.eemf_al = ACM.eemf_q * -sin(ACM.theta_d);
    ACM.eemf_be = ACM.eemf_q *  cos(ACM.theta_d);
    // ACM.theta_d__eemf = atan2(-ACM.eemf_al, ACM.eemf_be);
    ACM.theta_d__eemf = atan2(-ACM.eemf_al*sign(ACM.omg_elec), ACM.eemf_be*sign(ACM.omg_elec));

    // detect bad simulation
    if(isNumber(ACM.rpm)){
        return false;
    }else{
        printf("ACM.rpm is %g\n", ACM.rpm);
        return true;        
    }
}

double _LPF_Euler(double input, double *LPF_y__last, double time_const_inv, double h){ // See HighFreqInjection_dq_frame__SQUAREWAVE.afx
    double y;
    y = *LPF_y__last + h * time_const_inv * (input - *LPF_y__last);
    *LPF_y__last = y;
    return y;
}
double _HPF_Euler(double input, double *HPF_y__last, double time_const_inv, double h){
    *HPF_y__last = _LPF_Euler(input, HPF_y__last, time_const_inv, h);
    return input - *HPF_y__last;
}
double _BPF_Euler(double input, double *LPF_y__last, double time_const_inv_lpf, 
                          double *HPF_y__last, double time_const_inv_hpf, double h){
    double y;
    y = _LPF_Euler(input, LPF_y__last, time_const_inv_lpf, h);
    y = _HPF_Euler(    y, HPF_y__last, time_const_inv_hpf, h);
    return y;
}

#define LPF_TIME_CONST_INVERSE (2*2*M_PI)
void dynamics_lpf(double input, double *state, double *fx){
    fx[0] = LPF_TIME_CONST_INVERSE * ( input - *state );
}
void RK4_general(void (*pointer_dynamics)(), double input, double *state, double hs){
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
double test_signal_al = 0;
double test_signal_be = 0;
double test_signal_M = 0;
double test_signal_T = 0;
double M_lpf = 0;
double T_lpf = 0;
double M_hpf = 0;
double T_hpf = 0;
double theta_filter = 0;
void measurement(){
    // Executed every TS

    // Voltage measurement
    US_C(0) = CTRL.ual;
    US_C(1) = CTRL.ube;
    US_P(0) = US_C(0);
    US_P(1) = US_C(1);


    // Current measurement
    // theta_filter = 5*2*M_PI*CTRL.timebase;
    theta_filter = sm.theta_d;

    // test_signal_al = 100*cos(theta_filter) + 0.0*cos(500*2*M_PI*CTRL.timebase);
    // test_signal_be = 100*sin(theta_filter);
    // static int square_wave_internal_register = 1;
    // static int dfe_counter = 0; 
    // #define HFSI_VOLTAGE 5 // V
    // if(dfe_counter++==2){
    //     dfe_counter = 0;
    //     square_wave_internal_register *= -1;
    // }
    // test_signal_al += HFSI_VOLTAGE*square_wave_internal_register;

    test_signal_al = ACM.ial;
    test_signal_be = ACM.ibe;

    {
        test_signal_M = AB2M(test_signal_al, test_signal_be, cos(theta_filter), sin(theta_filter));
        test_signal_T = AB2T(test_signal_al, test_signal_be, cos(theta_filter), sin(theta_filter));

        // LPF
        RK4_general(dynamics_lpf, test_signal_M, &M_lpf, TS);
        RK4_general(dynamics_lpf, test_signal_T, &T_lpf, TS);
        // M_lpf = _LPF_Euler(test_signal_M, (sm.current_lpf_register+0), 10*2*M_PI, TS);
        // T_lpf = _LPF_Euler(test_signal_T, (sm.current_lpf_register+1), 10*2*M_PI, TS);
        IS_LPF(0) = MT2A(M_lpf, T_lpf, cos(theta_filter), sin(theta_filter));
        IS_LPF(1) = MT2B(M_lpf, T_lpf, cos(theta_filter), sin(theta_filter));

        // HPF
        M_hpf = test_signal_M - M_lpf;
        T_hpf = test_signal_T - T_lpf;
        // M_hpf = _HPF_Euler(test_signal_M, (sm.current_hpf_register+0), 400*2*M_PI, TS);
        // T_hpf = _HPF_Euler(test_signal_T, (sm.current_hpf_register+1), 400*2*M_PI, TS);
        IS_HPF(0) = MT2A(M_hpf, T_hpf, cos(theta_filter), sin(theta_filter));
        IS_HPF(1) = MT2B(M_hpf, T_hpf, cos(theta_filter), sin(theta_filter));

        // BPF
        // IS_BPF(0) = _BPF_Euler(test_signal_al, (sm.current_bpf_register1+0), 500*2*M_PI,
        //                                  (sm.current_bpf_register2+0), 50*2*M_PI, TS); // time constant is 1/400 <=> cutoff frequency is 400/(2*pi)
    }
    IS_C(0) = IS_LPF(0);
    IS_C(1) = IS_LPF(1);
    // IS_C(0) = ACM.ial;
    // IS_C(1) = ACM.ibe;

    // Position and speed measurement
    sm.omg_elec = ACM.x[2];
    sm.omg_mech = sm.omg_elec * sm.npp_inv;
    sm.theta_d = ACM.x[3];
}
void inverter_model(){

    // 根据给定电压CTRL.ual和实际的电机电流ACM.ial，计算畸变的逆变器输出电压ACM.ual。
    #if INVERTER_NONLINEARITY

        InverterNonlinearity_SKSul96(CTRL.ual, CTRL.ube, ACM.ial, ACM.ibe);
        // InverterNonlinearity_Tsuji01
        ACM.ual = UAL_C_DIST;
        ACM.ube = UBE_C_DIST;

        // Distorted voltage (for visualization purpose)
        DIST_AL = ACM.ual - CTRL.ual;
        DIST_BE = ACM.ube - CTRL.ube;
    #else
        ACM.ual = CTRL.ual;
        ACM.ube = CTRL.ube;
        // 仿真是在永磁体磁场定向系下仿真的哦
        ACM.ud = AB2M(ACM.ual, ACM.ube, cos(ACM.theta_d), sin(ACM.theta_d));
        ACM.uq = AB2T(ACM.ual, ACM.ube, cos(ACM.theta_d), sin(ACM.theta_d));
    #endif
}

int main(){
    
    printf("NUMBER_OF_STEPS: %d\n\n", NUMBER_OF_STEPS);

    /* Initialization */
    Machine_init();
    CTRL_init();
    sm_init();
    ob_init();

    FILE *fw;
    fw = fopen(DATA_FILE_NAME, "w");
    printf("%s\n", DATA_FILE_NAME);
    write_header_to_file(fw);

    /* MAIN LOOP */
    clock_t begin, end;
    begin = clock();
    int _; // _ for the outer iteration
    int dfe_counter=0; // dfe_counter for down frequency execution
    for(_=0;_<NUMBER_OF_STEPS;++_){

        // printf("%d\n", _);

        /* Command (Speed or Position) */
        // cmd_fast_speed_reversal(CTRL.timebase, 5, 5, 1500); // timebase, instant, interval, rpm_cmd
        // cmd_fast_speed_reversal(CTRL.timebase, 5, 5, 200); // timebase, instant, interval, rpm_cmd
        ACM.rpm_cmd = 10;

        /* Load Torque */
        // ACM.Tload = 0 * sign(ACM.rpm); // No-load test
        // ACM.Tload = ACM.Tem; // Blocked-rotor test
        ACM.Tload = 2 * sign(ACM.rpm);

        /* Simulated ACM */
        if(machine_simulation()){ 
            printf("Break the loop.\n");
            break;
        }

        if(++dfe_counter == TS_UPSAMPLING_FREQ_EXE_INVERSE){
            dfe_counter = 0;

            /* Time in DSP */
            CTRL.timebase += TS;

            measurement();

            // observation();

            write_data_to_file(fw);

            control(ACM.rpm_cmd, 0);
        }

        inverter_model();
    }
    end = clock(); printf("The simulation in C costs %g sec.\n", (double)(end - begin)/CLOCKS_PER_SEC);
    fclose(fw);

    /* Fade out */
    system("python ./ACMPlot.py"); 
    // getch();
    // system("pause");
    // system("exit");
    return 0; 
}

/* Utility */
void write_header_to_file(FILE *fw){
    // no space is allowed!
    fprintf(fw, "x0(id)[A],x1(iq)[A],x2(speed)[rad/s],x3(position)[rad],ud_cmd[V],uq_cmd[V],id_cmd[A],id_err[A],iq_cmd[A],iq_err[A],|eemf|[V],eemf_be[V],theta_d[rad],theta_d__eemf[rad],mismatch[rad],sin(mismatch)[rad],OB_POS,sin(ER_POS),OB_EEMF_BE,error(OB_EEMF),OB_OMG,er_omg,test_signal_al,test_signal_be,IS_LPF_Euler(0),IS_LPF_Euler(1),IS_HPF_Euler(0),IS_HPF_Euler(1),M_hpf,T_hpf\n");
    {
        FILE *fw2;
        fw2 = fopen("info.dat", "w");
        fprintf(fw2, "TS,DOWN_SAMPLE,DATA_FILE_NAME\n");
        fprintf(fw2, "%g, %d, %s\n", TS, DOWN_SAMPLE, DATA_FILE_NAME);
        fclose(fw2);
    }
}
void write_data_to_file(FILE *fw){
    static int bool_animate_on = false;
    static int j=0,jj=0; // j,jj for down sampling

    // if(CTRL.timebase>20)
    {
        if(++j == DOWN_SAMPLE)
        {
            j=0;
            fprintf(fw, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
                    ACM.x[0], ACM.x[1], ACM.x[2], ACM.x[3], CTRL.ud_cmd, CTRL.uq_cmd, 
                    CTRL.id_cmd, CTRL.id__fb-CTRL.id_cmd, CTRL.iq_cmd, CTRL.iq__fb-CTRL.iq_cmd, ACM.eemf_q, ACM.eemf_be,
                    ACM.theta_d, ACM.theta_d__eemf,ACM.theta_d-ACM.theta_d__eemf,sin(ACM.theta_d-ACM.theta_d__eemf),
                    OB_POS, sin(ACM.theta_d-OB_POS), OB_EEMF_BE, ACM.eemf_be-OB_EEMF_BE, OB_OMG, ACM.omg_elec-OB_OMG,
                    test_signal_al, test_signal_be, IS_LPF(0), IS_LPF(1), IS_HPF(0), IS_HPF(1),
                    M_hpf, T_hpf
                    );
        }
    }

    // if(bool_animate_on==false){
    //     bool_animate_on = true;
    //     printf("Start ACMAnimate\n");
    //     system("start python ./ACMAnimate.py"); 
    // }
}

int isNumber(double x){
    // This looks like it should always be true, 
    // but it's false if x is an NaN (1.#QNAN0).
    return (x == x); 
    // see https://www.johndcook.com/blog/IEEE_exceptions_in_cpp/ cb: https://stackoverflow.com/questions/347920/what-do-1-inf00-1-ind00-and-1-ind-mean
}



