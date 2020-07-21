#include "ACMSim.h"

/* PI Control
 * */
#define INCREMENTAL_PID TRUE
#if INCREMENTAL_PID
double PID(struct PID_Reg *r, double err){

    #define O_STATE r->o_state
    #define E_STATE r->e_state
    #define I_LIMIT r->i_limit

    double delta_u_n;
    delta_u_n = r->Kp * ( err - E_STATE ) + r->Ki * err;

    double output;
    output = O_STATE + delta_u_n;

    if(output > I_LIMIT)
        output = I_LIMIT;
    else if(output < -I_LIMIT)
        output = -I_LIMIT;

    E_STATE = err; 
    O_STATE = output;

    return output;

    #undef O_STATE
    #undef E_STATE
    #undef I_LIMIT
}
#else
double PID(struct PID_Reg *r, double err){

    #define DYNAMIC_CLAPMING TRUE

    #define I_STATE r->i_state
    #define I_LIMIT r->i_limit
    double output;
    double P_output;

    P_output = err * r->Kp; // 比例

    I_STATE += err * r->Ki; // 积分

    // 添加积分饱和特性
    #if DYNAMIC_CLAPMING
        // dynamic clamping
        if( I_STATE > I_LIMIT - P_output)
            I_STATE = I_LIMIT - P_output;
        else if( I_STATE < -I_LIMIT - P_output)
            I_STATE = -I_LIMIT - P_output;
    #else
        // static clamping
        if( I_STATE > I_LIMIT)
            I_STATE = I_LIMIT; 
        else if( I_STATE < -I_LIMIT)
            I_STATE = -I_LIMIT;
    #endif

    output = I_STATE + P_output;

    if(output > I_LIMIT)
        output = I_LIMIT;
    else if(output < -I_LIMIT)
        output = -I_LIMIT;
    return output;
    #undef I_STATE
    #undef I_LIMIT
}
#endif

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

    // CTRL.Tload = 0.0;
    // CTRL.rpm_cmd = 0.0;

    CTRL.npp = ACM.npp;
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

    CTRL.Tem = 0.0;
    CTRL.Tem_cmd = 0.0;

    // ver. IEMDC
    CTRL.PID_speed.Kp = SPEED_LOOP_PID_PROPORTIONAL_GAIN;
    CTRL.PID_speed.Ti = SPEED_LOOP_PID_INTEGRAL_TIME_CONSTANT;
    CTRL.PID_speed.Ki = CTRL.PID_speed.Kp / CTRL.PID_speed.Ti * (TS*SPEED_LOOP_CEILING); // 4.77 = 1 / (npp*1/60*2*pi)
    CTRL.PID_speed.i_limit = SPEED_LOOP_LIMIT_NEWTON_METER;
    CTRL.PID_speed.i_state = 0.0;
    printf("Speed PID: Kp=%g, Ki=%g, limit=%g Nm\n", CTRL.PID_speed.Kp, CTRL.PID_speed.Ki/TS, CTRL.PID_speed.i_limit);

    CTRL.PID_id.Kp = CURRENT_LOOP_PID_PROPORTIONAL_GAIN; // cutoff frequency of 1530 rad/s
    CTRL.PID_id.Ti = CURRENT_LOOP_PID_INTEGRAL_TIME_CONSTANT;
    CTRL.PID_id.Ki = CTRL.PID_id.Kp/CTRL.PID_id.Ti*TS; // =0.025
    CTRL.PID_id.i_limit = CURRENT_LOOP_LIMIT_VOLTS; //350.0; // unit: Volt
    CTRL.PID_id.i_state = 0.0;
    printf("Current PID: Kp=%g, Ki=%g, limit=%g V\n", CTRL.PID_id.Kp, CTRL.PID_id.Ki/TS, CTRL.PID_id.i_limit);

    CTRL.PID_iq.Kp = CURRENT_LOOP_PID_PROPORTIONAL_GAIN;
    CTRL.PID_iq.Ti = CURRENT_LOOP_PID_INTEGRAL_TIME_CONSTANT;
    CTRL.PID_iq.Ki = CTRL.PID_iq.Kp/CTRL.PID_iq.Ti*TS;
    CTRL.PID_iq.i_limit = CURRENT_LOOP_LIMIT_VOLTS; // unit: Volt, 350V->max 1300rpm
    CTRL.PID_iq.i_state = 0.0;
}
double theta_d_harnefors = 0.0;
double omg_harnefors = 0.0;
void harnefors_scvm(){
    #define KE_MISMATCH 1.0 // 0.7
    double d_axis_emf;
    double q_axis_emf;
    #define LAMBDA 2 // 2
    #define CJH_TUNING_A 25 // 1
    #define CJH_TUNING_B 1 // 1
    double lambda_s = LAMBDA * sign(omg_harnefors);
    double alpha_bw_lpf = CJH_TUNING_A*0.1*(1500*RPM_2_RAD_PER_SEC) + CJH_TUNING_B*2*LAMBDA*fabs(omg_harnefors);
    // d_axis_emf = CTRL.ud_cmd - 1*CTRL.R*CTRL.id_cmd + omg_harnefors*1.0*CTRL.Lq*CTRL.iq_cmd; // If Ld=Lq.
    // q_axis_emf = CTRL.uq_cmd - 1*CTRL.R*CTRL.iq_cmd - omg_harnefors*1.0*CTRL.Ld*CTRL.id_cmd; // If Ld=Lq.
    d_axis_emf = CTRL.ud_cmd - 1*CTRL.R*CTRL.id_cmd + omg_harnefors*1.0*CTRL.Lq*CTRL.iq_cmd; // eemf
    q_axis_emf = CTRL.uq_cmd - 1*CTRL.R*CTRL.iq_cmd - omg_harnefors*1.0*CTRL.Lq*CTRL.id_cmd; // eemf
    // Note it is bad habit to write numerical integration explictly like this. The states on the right may be accencidentally modified on the run.
    theta_d_harnefors += TS * omg_harnefors;
    omg_harnefors += TS * alpha_bw_lpf * ( (q_axis_emf - lambda_s*d_axis_emf)/(CTRL.KE*KE_MISMATCH+(CTRL.Ld-CTRL.Lq)*CTRL.id_cmd) - omg_harnefors );
    while(theta_d_harnefors>M_PI) theta_d_harnefors-=2*M_PI;
    while(theta_d_harnefors<-M_PI) theta_d_harnefors+=2*M_PI;   
}
void control(double speed_cmd, double speed_cmd_dot){
    // Input 1 is feedback: estimated speed/position or measured speed/position
    #if SENSORLESS_CONTROL
        #if SENSORLESS_CONTROL_HFSI
            CTRL.omg__fb     = hfsi.omg_elec;
            CTRL.theta_d__fb = hfsi.theta_d;
        #else
            // getch("Not Implemented");
            // CTRL.omg__fb    ;
            // CTRL.omega_syn ;

            // CTRL.omg__fb     = OB_OMG;
            // CTRL.theta_d__fb = OB_POS;

            harnefors_scvm();
            CTRL.omg__fb     = omg_harnefors;
            CTRL.theta_d__fb = theta_d_harnefors;
        #endif
    #else

        // from measurement() in main.c
        CTRL.omg__fb     = sm.omg_elec;
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
    if(vc_count++ == SPEED_LOOP_CEILING){
        // velocity control loop execution frequency is 40 times slower than current control loop execution frequency
        vc_count = 0;
        CTRL.omg_ctrl_err = CTRL.omg__fb - speed_cmd*RPM_2_RAD_PER_SEC;
        CTRL.iq_cmd = - PID(&CTRL.PID_speed, CTRL.omg_ctrl_err);

        // for plot
        CTRL.speed_ctrl_err = CTRL.omg_ctrl_err * RAD_PER_SEC_2_RPM;
    }

    // Measured current in d-q frame
    CTRL.id__fb = AB2M(CTRL.ial__fb, CTRL.ibe__fb, CTRL.cosT, CTRL.sinT);
    CTRL.iq__fb = AB2T(CTRL.ial__fb, CTRL.ibe__fb, CTRL.cosT, CTRL.sinT);

    // For luenberger position observer for HFSI
    CTRL.Tem     = CTRL.npp * (CTRL.KE*CTRL.iq__fb + (CTRL.Ld-CTRL.Lq)*CTRL.id__fb*CTRL.iq__fb);
    CTRL.Tem_cmd = CTRL.npp * (CTRL.KE*CTRL.iq_cmd + (CTRL.Ld-CTRL.Lq)*CTRL.id_cmd*CTRL.iq_cmd);

    // Voltage command in d-q frame
    double vd, vq;
    vd = - PID(&CTRL.PID_id, CTRL.id__fb-CTRL.id_cmd);
    vq = - PID(&CTRL.PID_iq, CTRL.iq__fb-CTRL.iq_cmd);

    // Current loop decoupling (skipped for now)
    CTRL.ud_cmd = vd;
    CTRL.uq_cmd = vq;

    #ifdef HFSI_ON
        // Extra excitation for observation
        {
            static int dfe_counter = 0; 
            if(dfe_counter++==HFSI_CEILING){
                dfe_counter = 0;
                hfsi.square_wave_internal_register *= -1;
            }
            // hfsi.square_wave_internal_register *= -1;
            CTRL.ud_cmd += HFSI_VOLTAGE*hfsi.square_wave_internal_register;
        }
    #endif

    // Voltage command in alpha-beta frame
    CTRL.ual = MT2A(CTRL.ud_cmd, CTRL.uq_cmd, CTRL.cosT, CTRL.sinT);
    CTRL.ube = MT2B(CTRL.ud_cmd, CTRL.uq_cmd, CTRL.cosT, CTRL.sinT);
}



/* Command */
void cmd_fast_speed_reversal(double timebase, double instant, double interval, double rpm_cmd){
    if(timebase > instant+2*interval){
        ACM.rpm_cmd = 1*1500 + rpm_cmd;
    }else if(timebase > instant+interval){
        ACM.rpm_cmd = 1*1500 + -rpm_cmd;
    }else if(timebase > instant){
        ACM.rpm_cmd = 1*1500 + rpm_cmd;
    }else{
        ACM.rpm_cmd = 20; // default initial command
    }
}

