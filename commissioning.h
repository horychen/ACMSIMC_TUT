
#ifndef COMMISSIONING_H
#define COMMISSIONING_H

struct CommissioningDataStruct{

    double timebase;

    int npp; // number of pole pairs
    double IN; // rated line current in Ampere RMS

    double R;
    double L;
    // double Ld;
    // double Lq;
    double KE;
    double Js; // shaft inertia

    struct PID_Reg PID_speed;
    struct PID_Reg PID_id;
    struct PID_Reg PID_iq;

    double current_command;

    double current_sum;
    double voltage_sum;
    long int counterSS;
    int bool_collecting;
    double iv_data[100][2];
};
extern struct CommissioningDataStruct COMM;

void COMM_init();
void commissioning();


#endif
