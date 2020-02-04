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

    CTRL.omg__fb = 0.0;
    CTRL.ial__fb = 0.0;
    CTRL.ibe__fb = 0.0;
    CTRL.psi_mu_al__fb = 0.0;
    CTRL.psi_mu_be__fb = 0.0;

    CTRL.rotor_flux_cmd = 0.0; // id=0 control

    CTRL.omg_ctrl_err = 0.0;
    CTRL.speed_ctrl_err = 0.0;

    CTRL.cosT = 1.0;
    CTRL.sinT = 0.0;

    CTRL.omega_syn = 0.0;

    CTRL.theta_d__fb = 0.0;
    CTRL.id__fb = 0.0;
    CTRL.iq__fb = 0.0;
    CTRL.ud_cmd = 0.0;
    CTRL.uq_cmd = 0.0;
    CTRL.id_cmd = 0.0;
    CTRL.iq_cmd = 0.0;

    // ver. IEMDC
    CTRL.PID_speed.Kp = 0.5; 
    CTRL.PID_speed.Ti = 5;
    CTRL.PID_speed.Ki = (CTRL.PID_speed.Kp*4.77) / CTRL.PID_speed.Ti * (TS*VC_LOOP_CEILING*DOWN_FREQ_EXE_INVERSE);
    CTRL.PID_speed.i_state = 0.0;
    CTRL.PID_speed.i_limit = 8;

    printf("Kp_omg=%g, Ki_omg=%g\n", CTRL.PID_speed.Kp, CTRL.PID_speed.Ki);

    CTRL.PID_id.Kp = 15; // cutoff frequency of 1530 rad/s
    CTRL.PID_id.Ti = 0.08;
    CTRL.PID_id.Ki = CTRL.PID_id.Kp/CTRL.PID_id.Ti*TS; // =0.025
    CTRL.PID_id.i_state = 0.0;
    CTRL.PID_id.i_limit = 650; //350.0; // unit: Volt

    CTRL.PID_iq.Kp = 15;
    CTRL.PID_iq.Ti = 0.08;
    CTRL.PID_iq.Ki = CTRL.PID_iq.Kp/CTRL.PID_iq.Ti*TS;
    CTRL.PID_iq.i_state = 0.0;
    CTRL.PID_iq.i_limit = 650; // unit: Volt, 350V->max 1300rpm

    printf("Kp_cur=%g, Ki_cur=%g\n", CTRL.PID_id.Kp, CTRL.PID_id.Ki);
}
void control(double speed_cmd, double speed_cmd_dot){
    // Input 1 is feedback: estimated speed/position or measured speed/position
    #if SENSORLESS_CONTROL
        getch("Not Implemented");
        // CTRL.omg__fb    ;
        // CTRL.omega_syn ;
    #else
        // from measurement() in main.c
        CTRL.omg__fb = sm.omg_elec;
        CTRL.theta_d__fb = sm.theta_d;
    #endif

    // Input 2 is feedback: measured current 
    CTRL.ial__fb = IS_C(0);
    CTRL.ibe__fb = IS_C(1);

    // Input 3 is the flux linkage command 
    #if CONTROL_STRATEGY == NULL_D_AXIS_CURRENT_CONTROL
        CTRL.rotor_flux_cmd = 0.0;
        CTRL.cosT = cos(CTRL.theta_d__fb); 
        CTRL.sinT = sin(CTRL.theta_d__fb);
    #else
        getch("Not Implemented");        
    #endif

    // d-axis current command
    CTRL.id_cmd = CTRL.rotor_flux_cmd / CTRL.Ld;

    // q-axis current command
    static int vc_count = 0;
    if(vc_count++==VC_LOOP_CEILING*DOWN_FREQ_EXE_INVERSE){ 
        vc_count = 0;
        CTRL.omg_ctrl_err = CTRL.omg__fb - speed_cmd*RPM_2_RAD_PER_SEC;
        CTRL.iq_cmd = - PID(&CTRL.PID_speed, CTRL.omg_ctrl_err);

        // for plot
        CTRL.speed_ctrl_err = CTRL.omg_ctrl_err * RAD_PER_SEC_2_RPM;
    }

    // Measured current in d-q frame
    CTRL.id__fb = AB2M(CTRL.ial__fb, CTRL.ibe__fb, CTRL.cosT, CTRL.sinT);
    CTRL.iq__fb = AB2T(CTRL.ial__fb, CTRL.ibe__fb, CTRL.cosT, CTRL.sinT);

    // Voltage command in d-q frame
    double vd, vq;
    vd = - PID(&CTRL.PID_id, CTRL.id__fb-CTRL.id_cmd);
    vq = - PID(&CTRL.PID_iq, CTRL.iq__fb-CTRL.iq_cmd);

    // Current loop decoupling (skipped for now)
    CTRL.ud_cmd = vd;
    CTRL.uq_cmd = vq;

    // Voltage command in alpha-beta frame
    CTRL.ual = MT2A(CTRL.ud_cmd, CTRL.uq_cmd, CTRL.cosT, CTRL.sinT);
    CTRL.ube = MT2B(CTRL.ud_cmd, CTRL.uq_cmd, CTRL.cosT, CTRL.sinT);
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

