#include "ACMSim.h"

/* PI Control
 * */
double PI(struct PI_Reg *r, double err){
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


#if MACHINE_TYPE == INDUCTION_MACHINE
#elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
/* Initialization */
struct ControllerForExperiment CTRL;
void CTRL_init(){
    int i=0,j=0;

    CTRL.timebase = 0.0;

    CTRL.ual = 0.0;
    CTRL.ube = 0.0;

    CTRL.R = ACM.R;
    CTRL.KE = ACM.KE;
    CTRL.Ld = ACM.Ld;
    CTRL.Lq = ACM.Lq;

    CTRL.Tload = 0.0;
    CTRL.rpm_cmd = 0.0;

    CTRL.Js = ACM.Js;
    CTRL.Js_inv = 1.0 / CTRL.Js;

    CTRL.omg_fb = 0.0;
    CTRL.ial_fb = 0.0;
    CTRL.ibe_fb = 0.0;
    CTRL.psi_mu_al_fb = 0.0;
    CTRL.psi_mu_be_fb = 0.0;

    CTRL.rotor_flux_cmd = 0.0; // id=0 control

    CTRL.omg_ctrl_err = 0.0;
    CTRL.speed_ctrl_err = 0.0;

    CTRL.iMs = 0.0;
    CTRL.iTs = 0.0;

    CTRL.theta_M = 0.0;
    CTRL.cosT = 1.0;
    CTRL.sinT = 0.0;

    CTRL.omega_syn = 0.0;

    CTRL.uMs_cmd = 0.0;
    CTRL.uTs_cmd = 0.0;
    CTRL.iMs_cmd = 0.0;
    CTRL.iTs_cmd = 0.0;

    // ver. IEMDC
    CTRL.pi_speed.Kp = 0.5; 
    CTRL.pi_speed.Ti = 5;
    CTRL.pi_speed.Ki = (CTRL.pi_speed.Kp*4.77) / CTRL.pi_speed.Ti * (TS*VC_LOOP_CEILING*DOWN_FREQ_EXE_INVERSE);
    CTRL.pi_speed.i_state = 0.0;
    CTRL.pi_speed.i_limit = 8;

    printf("Kp_omg=%g, Ki_omg=%g\n", CTRL.pi_speed.Kp, CTRL.pi_speed.Ki);

    CTRL.pi_iD.Kp = 15; // cutoff frequency of 1530 rad/s
    CTRL.pi_iD.Ti = 0.08;
    CTRL.pi_iD.Ki = CTRL.pi_iD.Kp/CTRL.pi_iD.Ti*TS; // =0.025
    CTRL.pi_iD.i_state = 0.0;
    CTRL.pi_iD.i_limit = 650; //350.0; // unit: Volt

    CTRL.pi_iQ.Kp = 15;
    CTRL.pi_iQ.Ti = 0.08;
    CTRL.pi_iQ.Ki = CTRL.pi_iQ.Kp/CTRL.pi_iQ.Ti*TS;
    CTRL.pi_iQ.i_state = 0.0;
    CTRL.pi_iQ.i_limit = 650; // unit: Volt, 350V->max 1300rpm



    CTRL.pi_iD_PR_pose.Kp = 15/3; // cutoff frequency of 1530 rad/s
    CTRL.pi_iD_PR_pose.Ti = 0.08;
    CTRL.pi_iD_PR_pose.Ki = CTRL.pi_iD_PR_pose.Kp/CTRL.pi_iD_PR_pose.Ti*TS; // =0.025
    CTRL.pi_iD_PR_pose.i_state = 0.0;
    CTRL.pi_iD_PR_pose.i_limit = 650; //350.0; // unit: Volt

    CTRL.pi_iD_PR_nese.Kp = 15/3; // cutoff frequency of 1530 rad/s
    CTRL.pi_iD_PR_nese.Ti = 0.08;
    CTRL.pi_iD_PR_nese.Ki = CTRL.pi_iD_PR_nese.Kp/CTRL.pi_iD_PR_nese.Ti*TS; // =0.025
    CTRL.pi_iD_PR_nese.i_state = 0.0;
    CTRL.pi_iD_PR_nese.i_limit = 650; //350.0; // unit: Volt



    CTRL.pi_iQ_PR_pose.Kp = 15/3; // cutoff frequency of 1530 rad/s
    CTRL.pi_iQ_PR_pose.Ti = 0.08;
    CTRL.pi_iQ_PR_pose.Ki = CTRL.pi_iQ_PR_pose.Kp/CTRL.pi_iQ_PR_pose.Ti*TS; // =0.025
    CTRL.pi_iQ_PR_pose.i_state = 0.0;
    CTRL.pi_iQ_PR_pose.i_limit = 650; //350.0; // unit: Volt

    CTRL.pi_iQ_PR_nese.Kp = 15/3; // cutoff frequency of 1530 rad/s
    CTRL.pi_iQ_PR_nese.Ti = 0.08;
    CTRL.pi_iQ_PR_nese.Ki = CTRL.pi_iQ_PR_nese.Kp/CTRL.pi_iQ_PR_nese.Ti*TS; // =0.025
    CTRL.pi_iQ_PR_nese.i_state = 0.0;
    CTRL.pi_iQ_PR_nese.i_limit = 650; //350.0; // unit: Volt

    printf("Kp_cur=%g, Ki_cur=%g\n", CTRL.pi_iD.Kp, CTRL.pi_iD.Ki);
}
void control(double speed_cmd, double speed_cmd_dot){
    // Input 1 is feedback: estimated speed or measured speed
    #if SENSORLESS_CONTROL
        getch("Not Implemented");
        // CTRL.omg_fb    ;
        // CTRL.omega_syn ;
    #else
        CTRL.omg_fb = sm.omg;
    #endif
    // Input 2 is feedback: measured current 
    CTRL.ial_fb = IS_C(0);
    CTRL.ibe_fb = IS_C(1);
    // Input 3 is the rotor d-axis position
    #if SENSORLESS_CONTROL
        getch("Not Implemented");
    #else
        CTRL.theta_M = sm.theta_d;
    #endif

    #if CONTROL_STRATEGY == NULL_D_AXIS_CURRENT_CONTROL
        // Flux (linkage) command
        CTRL.rotor_flux_cmd = 0.0;
    #endif

    // M-axis current command
    CTRL.iMs_cmd = CTRL.rotor_flux_cmd / CTRL.Ld;

    // T-axis current command
    static int vc_count = 0;
    if(vc_count++==VC_LOOP_CEILING*DOWN_FREQ_EXE_INVERSE){ 
        vc_count = 0;
        CTRL.omg_ctrl_err = CTRL.omg_fb - speed_cmd*RPM_2_RAD_PER_SEC;
        CTRL.iTs_cmd = - PI(&CTRL.pi_speed, CTRL.omg_ctrl_err);

        CTRL.speed_ctrl_err = CTRL.omg_ctrl_err * RAD_PER_SEC_2_RPM;
    }

    #if CONTROL_STRATEGY == NULL_D_AXIS_CURRENT_CONTROL
        CTRL.cosT = cos(CTRL.theta_M); 
        CTRL.sinT = sin(CTRL.theta_M);
    #endif

    // Measured current in M-T frame
    CTRL.iMs = AB2M(CTRL.ial_fb, CTRL.ibe_fb, CTRL.cosT, CTRL.sinT);
    CTRL.iTs = AB2T(CTRL.ial_fb, CTRL.ibe_fb, CTRL.cosT, CTRL.sinT);

    // Voltage command in M-T frame
    double vM, vT;
    vM = - PI(&CTRL.pi_iD, CTRL.iMs-CTRL.iMs_cmd);
    vT = - PI(&CTRL.pi_iQ, CTRL.iTs-CTRL.iTs_cmd);

    // Current loop decoupling (skipped for now)
    CTRL.uMs_cmd = vM;
    CTRL.uTs_cmd = vT;

    // Voltage command in alpha-beta frame
    CTRL.ual = MT2A(CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.cosT, CTRL.sinT);
    CTRL.ube = MT2B(CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.cosT, CTRL.sinT);
}
double regulator(struct PI_Reg *reg, double error, int bool_turn_on){
    double temp=0.0;
    if(bool_turn_on){
        temp = - PI(reg, error);
    }else{
        temp = 0;
        reg->i_state = 0;
    }
    return temp;
}
void selfcommissioning(){

    #define RATED_VOLTAGE 200 // Vrms (line-to-line)
    #define RATED_CURRENT 10 // Arms (line-to-line)
    #define RATED_FREQUENCY 50 // Hz
    #define RATED_MECHANICAL_SPEED 1500 // rpm

    int bool_alpha_axis_turn_on=true;
    int bool_beta_axis_turn_on=true;
    int bool_inductance_identification=false;

    double d_axis_current_amplitude_used_in_inductance_identification=0.0, q_axis_current_amplitude_used_in_inductance_identification = 0.0; // A
    double vd_PR_pose=0.0, vq_PR_pose=0.0, vd_PR_nese=0.0, vq_PR_nese=0.0;

    if(CTRL.timebase < 5){
        CTRL.ial_cmd = 0;
        CTRL.ibe_cmd = 0.2*RATED_CURRENT;
        bool_alpha_axis_turn_on = false;
        bool_beta_axis_turn_on = true;
    }else if(CTRL.timebase < 10){
        CTRL.ial_cmd = 0.2*RATED_CURRENT;
        CTRL.ibe_cmd = 0;
        bool_alpha_axis_turn_on = true;
        bool_beta_axis_turn_on = false;
    }else if(CTRL.timebase < 30){
        CTRL.ial_cmd = 0.2*RATED_CURRENT;
        CTRL.ibe_cmd = 0;

        bool_alpha_axis_turn_on = true;
        bool_beta_axis_turn_on = true;

        bool_inductance_identification = true;
        d_axis_current_amplitude_used_in_inductance_identification = 0.01*RATED_CURRENT;
        q_axis_current_amplitude_used_in_inductance_identification = 0.00*RATED_CURRENT;
    }else{
        CTRL.ial_cmd = 0;
        CTRL.ibe_cmd = 0;
        bool_alpha_axis_turn_on = true;
        bool_beta_axis_turn_on = true;

        bool_inductance_identification = false;
        d_axis_current_amplitude_used_in_inductance_identification = 0.0;
        q_axis_current_amplitude_used_in_inductance_identification = 0.0;
    }

    double val, vbe;
    val    = regulator(&CTRL.pi_iD, IS_C(0) - CTRL.ial_cmd, bool_alpha_axis_turn_on);
    vbe    = regulator(&CTRL.pi_iQ, IS_C(1) - CTRL.ibe_cmd, bool_beta_axis_turn_on);

    // 正序测量电流
    CTRL.ids_pose = AB2M(IS_C(0), IS_C(1), cos(50*2*M_PI*CTRL.timebase), sin(50*2*M_PI*CTRL.timebase));
    CTRL.iqs_pose = AB2T(IS_C(0), IS_C(1), cos(50*2*M_PI*CTRL.timebase), sin(50*2*M_PI*CTRL.timebase));

    vd_PR_pose = regulator(&CTRL.pi_iD_PR_pose, CTRL.ids_pose - d_axis_current_amplitude_used_in_inductance_identification, bool_inductance_identification);
    vq_PR_pose = regulator(&CTRL.pi_iQ_PR_pose, CTRL.iqs_pose - q_axis_current_amplitude_used_in_inductance_identification, bool_inductance_identification);

    // 负序测量电流
    CTRL.ids_nese = AB2M(IS_C(0), IS_C(1), cos(-50*2*M_PI*CTRL.timebase), sin(-50*2*M_PI*CTRL.timebase));
    CTRL.iqs_nese = AB2T(IS_C(0), IS_C(1), cos(-50*2*M_PI*CTRL.timebase), sin(-50*2*M_PI*CTRL.timebase));

    vd_PR_nese = regulator(&CTRL.pi_iD_PR_nese, CTRL.ids_nese - d_axis_current_amplitude_used_in_inductance_identification, bool_inductance_identification);
    vq_PR_nese = regulator(&CTRL.pi_iQ_PR_nese, CTRL.iqs_nese - q_axis_current_amplitude_used_in_inductance_identification, bool_inductance_identification);

    double val_PR_pose, vbe_PR_pose, val_PR_nese, vbe_PR_nese;
    val_PR_pose = MT2A(vd_PR_pose, vq_PR_pose, cos(50*2*M_PI*CTRL.timebase), sin(50*2*M_PI*CTRL.timebase));
    vbe_PR_pose = MT2B(vd_PR_pose, vq_PR_pose, cos(50*2*M_PI*CTRL.timebase), sin(50*2*M_PI*CTRL.timebase));
    val_PR_nese = MT2A(vd_PR_nese, vq_PR_nese, cos(-50*2*M_PI*CTRL.timebase), sin(-50*2*M_PI*CTRL.timebase));
    vbe_PR_nese = MT2B(vd_PR_nese, vq_PR_nese, cos(-50*2*M_PI*CTRL.timebase), sin(-50*2*M_PI*CTRL.timebase));

    // Current loop decoupling (skipped for now)
    CTRL.ual_cmd = val + val_PR_pose + val_PR_nese;
    CTRL.ube_cmd = vbe + vbe_PR_pose + vbe_PR_nese;

    // CTRL.ual_cmd = val + cos(50*2*M_PI*CTRL.timebase);
    // CTRL.ube_cmd = vbe + 0*cos(50*2*M_PI*CTRL.timebase);

    // // Voltage command in alpha-beta frame
    // CTRL.ual = MT2A(CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.cosT, CTRL.sinT);
    // CTRL.ube = MT2B(CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.cosT, CTRL.sinT);
    CTRL.ual = CTRL.ual_cmd;
    CTRL.ube = CTRL.ube_cmd;
}

#endif



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


