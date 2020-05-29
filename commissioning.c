#include "ACMSim.h"

struct CommissioningDataStruct COMM;

void COMM_init(){
    COMM.timebase = 0.0;

    COMM.PID_id.Kp = CURRENT_LOOP_PID_PROPORTIONAL_GAIN; // cutoff frequency of 1530 rad/s
    COMM.PID_id.Ti = CURRENT_LOOP_PID_INTEGRAL_TIME_CONSTANT;
    COMM.PID_id.Ki = COMM.PID_id.Kp/COMM.PID_id.Ti*TS; // =0.025
    COMM.PID_id.i_limit = CURRENT_LOOP_LIMIT_VOLTS; //350.0; // unit: Volt
    COMM.PID_id.i_state = 0.0;
    printf("Current PID: Kp=%g, Ki=%g, limit=%g V\n", COMM.PID_id.Kp, COMM.PID_id.Ki/TS, COMM.PID_id.i_limit);

    COMM.PID_iq.Kp = CURRENT_LOOP_PID_PROPORTIONAL_GAIN;
    COMM.PID_iq.Ti = CURRENT_LOOP_PID_INTEGRAL_TIME_CONSTANT;
    COMM.PID_iq.Ki = COMM.PID_iq.Kp/COMM.PID_iq.Ti*TS;
    COMM.PID_iq.i_limit = CURRENT_LOOP_LIMIT_VOLTS; // unit: Volt, 350V->max 1300rpm
    COMM.PID_iq.i_state = 0.0;

    // R
    COMM.current_sum = 0.0;
    COMM.voltage_sum = 0.0;
    COMM.counterSS = 0;
    COMM.bool_collecting = false;
    int i;
    for(i=0;i<100;++i){
        COMM.iv_data[i][0] = 0.0;
        COMM.iv_data[i][1] = 0.0;
    }

    // L
    COMM.L = 0.999;
}

/* standard lib */
#include <stdio.h> // printf #include <stdbool.h> // bool for _Bool and true for 1
#include <conio.h> // for clrscr, and getch()
#include "stdlib.h" // for rand()
#include "math.h"

#define SS_RATED_CURRENT_RATIO 1e-2
#define SS_CEILING ((long int)(0.2/TS))
int reachSteadyStateCurrent(double current_error, double rated_current){
    static long int counterSS = 0;

    if(fabs(current_error) < rated_current * SS_RATED_CURRENT_RATIO){
        counterSS += 1;

        // Avoid to collect over-shoot data
        if(counterSS > SS_CEILING){
            counterSS = 0;
            return true; // 目前是一旦判断为稳态，就永远返回true。当然，也可以设计成回差的形式，达到稳态后还要判断是否脱离稳态。
        }
    }
    return false;
}

struct CommissioningDataStruct COMM;

/* event array and enum below must be in sync! */
enum state_codes { _nameplateData=0, _currentSensor=1, _resistanceId=2, _inductanceId=3, _initialPosId=4,  _PMFluxLinkageId=5, _inertiaId=6, _unexpected=7, _end=8 };
enum ret_codes   { ok=0, repeat=1, quit=2 };

int event_nameplateData(void){
    printf("Start self-commissioning!\n");
    COMM.npp = 2;
    COMM.IN = 8;
    printf("Name plate data: npp=%d, IN=%g.\n", COMM.npp, COMM.IN);

    CTRL.ual = 0.0;
    CTRL.ube = 0.0;

    return ok;
}
int event_currentSensor(void){
    static long int counterEntered = 0;
    if(counterEntered++==0)
        printf("Current Sensor Calibration (pass).\n");

    CTRL.ual = 0.0;
    CTRL.ube = 0.0;

    return ok;
}
int event_resistanceId(void){
    static long int counterEntered = 0;
    if(counterEntered++==0)
        printf("Resistance Identification.\n Current [A] | Voltage [V]\n");

    #define RS_ID_NUMBER_OF_STAIRS 10
    #define RS_ID_MAXIMUM_CURRENT COMM.IN
    double current_increase_per_step = RS_ID_MAXIMUM_CURRENT / RS_ID_NUMBER_OF_STAIRS;
    static int i = 0;
    COMM.current_command = - RS_ID_MAXIMUM_CURRENT;
    COMM.current_command += current_increase_per_step*i;

    double error = IS_C(0) - COMM.current_command;
    CTRL.ual = - PID(&COMM.PID_id, error);
    CTRL.ube = 0.0;

    // collect steady state data
    if(COMM.bool_collecting){
        COMM.current_sum += IS_C(0);
        COMM.voltage_sum += CTRL.ual;
        COMM.counterSS += 1;

        if(COMM.counterSS>800){
            COMM.iv_data[i][0] = COMM.current_sum/(double)COMM.counterSS;
            COMM.iv_data[i][1] = COMM.voltage_sum/(double)COMM.counterSS;
            printf("%f, %f\n", COMM.iv_data[i][0], COMM.iv_data[i][1]);
            COMM.bool_collecting = false;
            ++i;
        }
    }else{
        // reset
        COMM.current_sum = 0.0;
        COMM.voltage_sum = 0.0;
        COMM.counterSS = 0;

        // check steady state and assign boolean variable
        if(reachSteadyStateCurrent(error, RS_ID_MAXIMUM_CURRENT)){
            if(COMM.current_command > RS_ID_MAXIMUM_CURRENT){

                // Get resistance value
                COMM.R = (COMM.iv_data[0][1] - COMM.iv_data[1][1]) / (COMM.iv_data[0][0] - COMM.iv_data[1][0]);
                // COMM.R = (COMM.iv_data[18+0][1] - COMM.iv_data[18+1][1]) / (COMM.iv_data[18+0][0] - COMM.iv_data[18+1][0]);
                printf("R=%g Ohm (including inverter's on-resistance)\n", COMM.R);
                return ok;
            }

            COMM.bool_collecting = true;
        }
    }

    return repeat;
}
int event_inductanceId(void){
    static double last_voltage_command;
    static long int counterEntered = 0;
    if(counterEntered++==0){
        printf("Inductance Identification.\n");
        last_voltage_command = CTRL.ual;
    }

    double Delta_IS_al;
    double Delta_US_al;
    double IS_slope;

    Delta_US_al = 0.5*last_voltage_command;
    Delta_IS_al = IS_C(0) - IS_P(0);
    IS_slope = Delta_IS_al / TS;

    if(counterEntered<20){
        printf("L=%g\n", Delta_US_al/IS_slope );
        if(Delta_US_al/IS_slope < COMM.L){
            COMM.L = Delta_US_al/IS_slope;
        }
    }else{
        if(counterEntered==20){
            printf("COMM.L = %g\n", COMM.L);
        }
    }

    CTRL.ual = last_voltage_command + Delta_US_al;
    CTRL.ube = 0.0;

    if(counterEntered>800){
        return ok;
    }else{
        return repeat;
    }
}
int event_initialPosId(void){
    printf("Initial Position Identification (pass).\n");
    return ok;
}
int event_PMFluxLinkageId(void){
    static long int counterEntered = 0;
    if(counterEntered++==0)
        printf("PM Flux Linkage Identification.\n");

    double rpm_speed_command = 100; 
    double ud_cmd;
    double uq_cmd;
    ud_cmd = CTRL.ud_cmd;
    uq_cmd = CTRL.uq_cmd;
    control(rpm_speed_command, 0);

    static double last_KE;
    last_KE = COMM.KE;

    COMM.KE = (uq_cmd - COMM.R*CTRL.iq__fb) / (rpm_speed_command*RPM_2_RAD_PER_SEC) - COMM.L*CTRL.id__fb;

    if(fabs(COMM.KE - last_KE)>=1e-3){
        counterEntered = 1; // 复用一下这个计数器用于判断KE是否收敛，可以少声明一个计数器。
    }

    if(counterEntered>16000){
        printf("KE = %g\n", COMM.KE);
        return ok;
    }else{
        return repeat;
    }
}
#define TEST_SIGNAL_FREQUENCY 0.5
#define TEST_SIGNAL_PERIOD (2)
#define SPEED_COMMAND_BIAS 200
#define SPEED_COMMAND_RANGE 10
int event_inertiaId(void){

    static long int counterEntered = 0;
    if(counterEntered++==0)
        printf("Inertia Identification.\n");

    #define STATOR_CURRENT_ALPHA IS_C(0)
    #define STATOR_CURRENT_BETA IS_C(1)
    #define ROTOR_FLUX_ALPHA MT2A(COMM.KE, 0.0, CTRL.cosT, CTRL.sinT)
    #define ROTOR_FLUX_BETA  MT2B(COMM.KE, 0.0, CTRL.cosT, CTRL.sinT)
    #define NORMINAL_INERTIA (0.10) // 随便给个数量级差不多的数，不同数值只对画图有影响
    #define SPEED_SIGNAL (CTRL.omg__fb)

    #define AWAYA_LAMBDA (31.4)

    // Inertia Observer 2 Awaya1992 
    static double t = 0.0;
    t += TS;

    static double q0 = 0.0;
    static double q1 = 0.0;
    static double q1_dot = 0.0;
    static double q2 = 0.0;

    static double tau_est = 0.0;

    // double Tem = CTRL.Tem;
    double Tem = COMM.npp*( STATOR_CURRENT_ALPHA * -ROTOR_FLUX_BETA + STATOR_CURRENT_BETA * ROTOR_FLUX_ALPHA);
    q0 += TS * AWAYA_LAMBDA*( -q0 + Tem);

    q1_dot = AWAYA_LAMBDA*( -q1 + SPEED_SIGNAL/COMM.npp );
    q1 += TS * q1_dot;

    tau_est = -(NORMINAL_INERTIA)*q1_dot + q0;

    static double sum_A = 0.0;
    static double sum_B = 0.0;
    static double est_Js_variation = 0.0;
    static double est_Js = NORMINAL_INERTIA;

    sum_A += TS * (tau_est*q1_dot);
    sum_B += TS * (q1_dot*q1_dot);
    if(t >= TEST_SIGNAL_PERIOD){      
        if(sum_B<0.001){
            sum_B = 0.001; // avoid zero denominator
        }
        est_Js_variation = + sum_A / sum_B;
        est_Js = NORMINAL_INERTIA + est_Js_variation;
        t = 0.0;
        sum_A = 0.0;
        sum_B = 0.0;
    }

    // 限幅，只影响画图
    if(est_Js>0.1){
        est_Js = 0.1;
    }else if(est_Js<0.001){
        est_Js = 0.001;
    }
    COMM.Js = est_Js;

    // double *rpm_cmd;
    // double *rpm_cmd_dot;
    // periodicRampSpeedCommand(rpm_cmd, rpm_cmd_dot);
    // printf("%g %g\n", *rpm_cmd, *rpm_cmd_dot);
    // getch();
    // control(*rpm_cmd, *rpm_cmd_dot);

    control(SPEED_COMMAND_BIAS+SPEED_COMMAND_RANGE*cos(2*M_PI*TEST_SIGNAL_FREQUENCY*CTRL.timebase), 0.0);

    return repeat;
}
int event_unexpected(void){
    printf("Debug!");
    return quit;
}
int event_end(void){
    printf("Exit!");
    printf("END.");
    return quit;
}

// Reference: https://www.geeksforgeeks.org/enumeration-enum-c/
int (* event[])(void) = { event_nameplateData, event_currentSensor, event_resistanceId, event_inductanceId, event_initialPosId, 
                          event_PMFluxLinkageId, event_inertiaId, event_unexpected, event_end };
int lookup_transitions[][3] = {
                // return codes:
                       //      ok            repeat           quit
    [_nameplateData]   = {_currentSensor,    _nameplateData , _end},
    [_currentSensor]   = {_resistanceId ,    _currentSensor , _end},
    [_resistanceId]    = {_inductanceId ,    _resistanceId , _end},
    [_inductanceId]    = {_initialPosId ,    _inductanceId , _end},
    [_initialPosId]    = {_PMFluxLinkageId , _initialPosId , _end},
    [_PMFluxLinkageId] = {_inertiaId ,       _PMFluxLinkageId , _end},
    [_inertiaId]       = {_end ,             _inertiaId , _end},
    [_unexpected]      = {_end ,             _unexpected , _end},
    /* transitions from end state aren't needed */
};

#define ENTRY_STATE _nameplateData
#define END_STATE   _end

void commissioning(){

    static enum state_codes cur_state = ENTRY_STATE;
    static enum ret_codes rc;
    static int (* state_func)(void);

    if(cur_state!=END_STATE){
        state_func = event[cur_state];
        rc = state_func();
        cur_state = lookup_transitions[cur_state][rc];
    }

    IS_P(0) = IS_C(0);
    IS_P(1) = IS_C(1);
}
