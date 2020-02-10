#include "ACMSim.h"

/* PI Control
 * */
double PID(struct PID_Reg *r, double err){
    #define I_STATE r->i_state
    #define I_LIMIT r->i_limit
    double output;
    I_STATE += err * r->Ki;    // 积分
    if( I_STATE > I_LIMIT)     // 添加积分饱和特性
        I_STATE = I_LIMIT; 
    else if( I_STATE < -I_LIMIT)
        I_STATE = -I_LIMIT;

    output = I_STATE + err * r->Kp;

    if(output > I_LIMIT)
        output = I_LIMIT;
    else if(output < -I_LIMIT)
        output = -I_LIMIT;
    return output;
    #undef I_STATE
    #undef I_LIMIT
}



/* Initialization */
struct ControllerForExperiment CTRL;
void CTRL_init(){
    int i=0,j=0;

    CTRL.timebase = 0.0;

        /* Parameter (including speed) Adaptation */ 
        CTRL.rs     = ACM.rs;
        CTRL.rreq   = ACM.rreq;
        CTRL.Lsigma = ACM.Lsigma;
        CTRL.alpha  = ACM.alpha;
        CTRL.Lmu    = ACM.Lmu;
        CTRL.Lmu_inv = 1.0/ACM.Lmu;
        CTRL.Js     = ACM.Js;
        CTRL.Js_inv = 1.0/ACM.Js;    

    CTRL.ual = 0.0;
    CTRL.ube = 0.0;

    CTRL.rpm_cmd = 0.0;
    CTRL.rotor_flux_cmd = 0.5;

    CTRL.omg_ctrl_err = 0.0;
    CTRL.speed_ctrl_err = 0.0;

    CTRL.omg_fb = 0.0;
    CTRL.ial_fb = 0.0;
    CTRL.ibe_fb = 0.0;
    CTRL.psi_mu_al_fb = 0.0;
    CTRL.psi_mu_be_fb = 0.0;

    CTRL.theta_M = 0.0;
    CTRL.cosT = 0.0;
    CTRL.sinT = 1;

    CTRL.uMs_cmd = 0.0;
    CTRL.uTs_cmd = 0.0;
    CTRL.iMs_cmd = 0.0;
    CTRL.iTs_cmd = 0.0;

    CTRL.omega_syn = 0.0;
    CTRL.omega_sl = 0.0;

    // ver. IEMDC
    CTRL.PID_speed.Kp = SPEED_LOOP_PID_PROPORTIONAL_GAIN;
    CTRL.PID_speed.Ti = SPEED_LOOP_PID_INTEGRAL_TIME_CONSTANT;
    CTRL.PID_speed.Ki = CTRL.PID_speed.Kp / CTRL.PID_speed.Ti * (TS*SPEED_LOOP_CEILING); // 4.77 = 1 / (npp*1/60*2*pi)
    CTRL.PID_speed.i_limit = SPEED_LOOP_LIMIT_NEWTON_METER;
    CTRL.PID_speed.i_state = 0.0;
    printf("Speed PID: Kp=%g, Ki=%g, limit=%g Nm\n", CTRL.PID_speed.Kp, CTRL.PID_speed.Ki/TS, CTRL.PID_speed.i_limit);

    CTRL.PID_iMs.Kp = CURRENT_LOOP_PID_PROPORTIONAL_GAIN; // cutoff frequency of 1530 rad/s
    CTRL.PID_iMs.Ti = CURRENT_LOOP_PID_INTEGRAL_TIME_CONSTANT;
    CTRL.PID_iMs.Ki = CTRL.PID_iMs.Kp/CTRL.PID_iMs.Ti*TS; // =0.025
    CTRL.PID_iMs.i_limit = CURRENT_LOOP_LIMIT_VOLTS; //350.0; // unit: Volt
    printf("Current PID: Kp=%g, Ki=%g, limit=%g V\n", CTRL.PID_iMs.Kp, CTRL.PID_iMs.Ki/TS, CTRL.PID_iMs.i_limit);
    CTRL.PID_iMs.i_state = 0.0;

    CTRL.PID_iTs.Kp = CURRENT_LOOP_PID_PROPORTIONAL_GAIN;
    CTRL.PID_iTs.Ti = CURRENT_LOOP_PID_INTEGRAL_TIME_CONSTANT;
    CTRL.PID_iTs.Ki = CTRL.PID_iTs.Kp/CTRL.PID_iTs.Ti*TS;
    CTRL.PID_iTs.i_limit = CURRENT_LOOP_LIMIT_VOLTS; // unit: Volt, 350V->max 1300rpm
    CTRL.PID_iTs.i_state = 0.0;
}
void control(double speed_cmd, double speed_cmd_dot){
    // OPEN LOOP CONTROL
    #if CONTROL_STRATEGY == VVVF_CONTROL
        #define VF_RATIO 18 //18.0 // 8 ~ 18 shows saturated phenomenon
        double freq = 2; // 0.15 ~ 0.5 ~ 2 （0.1时电压李萨茹就变成一个圆了）
        double volt = VF_RATIO*freq;
        CTRL.ual = volt*cos(2*M_PI*freq*CTRL.timebase);
        CTRL.ube = volt*sin(2*M_PI*freq*CTRL.timebase);
        return;
    #endif

    // Input 1 is feedback: estimated speed or measured speed
    #if SENSORLESS_CONTROL
        CTRL.omg_fb = OB_OMG;
        // #if OBSERVER_APPLIED == TAJIMA96
        //     CTRL.omega_syn = ob.tajima.omega_syn;
        //     CTRL.omega_sl  = ob.tajima.omega_sl;
        // #endif
    #else
        CTRL.omg_fb = im.omg;
    #endif
    // Input 2 is feedback: measured current 
    CTRL.ial_fb = IS_C(0);
    CTRL.ibe_fb = IS_C(1);
    // Input 3 differs for DFOC and IFOC
    #if CONTROL_STRATEGY == DFOC
        // DFOC: estimated flux components in alpha-beta frame
        CTRL.psi_mu_al_fb = OB_FLUX(0);
        CTRL.psi_mu_be_fb = OB_FLUX(1);
    #elif CONTROL_STRATEGY == IFOC
        // IFOC: estimated rotor resistance
        CTRL.rreq = OB_RREQ;
    #else
    #endif

    // Flux (linkage) command
    CTRL.rotor_flux_cmd = IM_FLUX_COMMAND_VALUE; // f(speed, dc bus voltage, last torque current command)
    // CTRL.rotor_flux_cmd = 3;
        // 1. speed is compared with the base speed to decide flux weakening or not
        // 2. dc bus voltage is required for certain application
        // 3. last torque current command is required for loss minimization

    // M-axis current command
    CTRL.iMs_cmd = CTRL.rotor_flux_cmd*CTRL.Lmu_inv + M1*OMG1*cos(OMG1*CTRL.timebase) / CTRL.rreq;
    // printf("%g, %g, %g\n", CTRL.Lmu_inv, CTRL.iMs_cmd, CTRL.iTs_cmd);

    // T-axis current command
    static int vc_count = 0;
    if(vc_count++==SPEED_LOOP_CEILING){ 
        vc_count = 0;
        CTRL.omg_ctrl_err = CTRL.omg_fb - speed_cmd*RPM_2_RAD_PER_SEC;
        CTRL.iTs_cmd = - PID(&CTRL.PID_speed, CTRL.omg_ctrl_err);

        CTRL.speed_ctrl_err = CTRL.omg_ctrl_err * RAD_PER_SEC_2_RPM;
    }


    #if CONTROL_STRATEGY == DFOC
        // feedback field orientation
        double modulus = sqrt(CTRL.psi_mu_al_fb*CTRL.psi_mu_al_fb + CTRL.psi_mu_be_fb*CTRL.psi_mu_be_fb);
        if(modulus<1e-3){
            CTRL.cosT = 1;
            CTRL.sinT = 0;
        }else{
            CTRL.cosT = CTRL.psi_mu_al_fb / modulus;
            CTRL.sinT = CTRL.psi_mu_be_fb / modulus;
        }
    #elif CONTROL_STRATEGY == IFOC
        // Feed-forward field orientation
        CTRL.theta_M += TS * CTRL.omega_syn;

        if(CTRL.theta_M > M_PI){
            CTRL.theta_M -= 2*M_PI;
        }else if(CTRL.theta_M < -M_PI){
            CTRL.theta_M += 2*M_PI; // 反转！
        }

        #if VOLTAGE_CURRENT_DECOUPLING_CIRCUIT
            CTRL.omega_sl = CTRL.rreq*CTRL.iTs / CTRL.rotor_flux_cmd;
        #else
            CTRL.omega_sl = CTRL.rreq*CTRL.iTs_cmd / CTRL.rotor_flux_cmd;
        #endif

        CTRL.omega_syn = CTRL.omg_fb + CTRL.omega_sl;

        CTRL.cosT = cos(CTRL.theta_M); 
        CTRL.sinT = sin(CTRL.theta_M);
    #endif

    // Measured current in M-T frame
    CTRL.iMs = AB2M(CTRL.ial_fb, CTRL.ibe_fb, CTRL.cosT, CTRL.sinT);
    CTRL.iTs = AB2T(CTRL.ial_fb, CTRL.ibe_fb, CTRL.cosT, CTRL.sinT);

    // Voltage command in M-T frame
    double vM, vT;
    vM = - PID(&CTRL.PID_iMs, CTRL.iMs-CTRL.iMs_cmd);
    vT = - PID(&CTRL.PID_iTs, CTRL.iTs-CTRL.iTs_cmd);

    // Current loop decoupling (skipped, see Chen.Huang-Stable)
    {   // Steady state dynamics based decoupling circuits for current regulation
        #if VOLTAGE_CURRENT_DECOUPLING_CIRCUIT == TRUE
            CTRL.uMs_cmd = vM + (CTRL.Lsigma)         *(-CTRL.omega_syn*CTRL.iTs); // Telford03/04
            CTRL.uTs_cmd = vT + CTRL.omega_syn*(CTRL.rotor_flux_cmd + CTRL.Lsigma*CTRL.iMs); // 这个行，但是无速度运行时，会导致M轴电流在转速暂态高频震荡。
            // CTRL.uMs_cmd = vM;
            // CTRL.uTs_cmd = vT;
        #else
            CTRL.uMs_cmd = vM;
            CTRL.uTs_cmd = vT;
        #endif
    }

    // Voltage command in alpha-beta frame
    CTRL.ual = MT2A(CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.cosT, CTRL.sinT);
    CTRL.ube = MT2B(CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.cosT, CTRL.sinT);
}



/* Command */
void cmd_fast_speed_reversal(double timebase, double instant, double interval, double rpm_cmd){
    if(timebase > instant+2*interval){
        ACM.rpm_cmd = rpm_cmd;
    }else if(timebase > instant+interval){
        ACM.rpm_cmd = -rpm_cmd;
    }else if(timebase > instant){
        ACM.rpm_cmd = rpm_cmd;
    }
}


void cmd_slow_speed_reversal(double timebase, double instant, double interval, double rpm_cmd){
    if(      timebase > instant+2*interval){
        ACM.rpm_cmd += TS*2.0;
        if(ACM.rpm_cmd<rpm_cmd){
            ACM.rpm_cmd = rpm_cmd;
        }
    }else if(timebase > instant+interval){
        ACM.rpm_cmd -= TS*2.0;
        if(ACM.rpm_cmd<-rpm_cmd){
            ACM.rpm_cmd = -rpm_cmd;
        }
    }else if(timebase > instant){
        ACM.rpm_cmd = rpm_cmd;
    }
}
