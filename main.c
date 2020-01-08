#include "ACMSim.h"

double sign(double x){
    return (x > 0) - (x < 0);    
}
double fabs(double x){
    return (x >= 0) ? x : -x;
}

#if MACHINE_TYPE == INDUCTION_MACHINE
    struct InductionMachineSimulated ACM;
    void Machine_init(){
        int i;
        for(i=0;i<5;++i){
            ACM.x[i] = 0.0;
        }
        ACM.rpm = 0.0;
        ACM.rpm_cmd = 0.0;
        ACM.rpm_deriv_cmd = 0.0;
        ACM.Tload = 0.0;
        ACM.Tem = 0.0;

        ACM.Lmu    = 0.4482;
            // Those parameters are introduced since branch saturation
            ACM.Lls = 0.0126;
            ACM.Llr = 0.0126;
            ACM.Lm  = 0.5*(ACM.Lmu + sqrt(ACM.Lmu*ACM.Lmu + 4*ACM.Llr*ACM.Lmu));
            ACM.Lm_slash_Lr = ACM.Lm/(ACM.Lm+ACM.Llr);
            ACM.Lr_slash_Lm = (ACM.Lm+ACM.Llr)/ACM.Lm; // i_dreq = ACM.idr*ACM.Lr_slash_Lm
            ACM.LSigmal = 1.0 / (1.0 / ACM.Lls + 1.0 / ACM.Llr);
        ACM.Lsigma = ACM.Lm + ACM.Lls - ACM.Lmu; // = Ls * (1.0 - Lm*Lm/Ls/Lr);
        printf("Validate: %g = %g?\n", ACM.Lsigma, (ACM.Lm+ACM.Lls) * (1.0 - ACM.Lm*ACM.Lm/(ACM.Lm+ACM.Lls)/(ACM.Lm+ACM.Llr)) );

        ACM.rreq   = 1.69;
        ACM.rs     = 3.04;
            ACM.rr = ACM.rreq * ACM.Lr_slash_Lm*ACM.Lr_slash_Lm;

        ACM.alpha  = ACM.rreq / (ACM.Lmu);
        ACM.Lmu_inv= 1.0/ACM.Lmu;

        ACM.Js = 0.0636; // Awaya92 using im.omg
        ACM.npp = 2;
        ACM.mu_m = ACM.npp/ACM.Js;

        ACM.Ts  = MACHINE_TS;

        ACM.ial = 0.0;
        ACM.ibe = 0.0;

        ACM.ual = 0.0;
        ACM.ube = 0.0;

            // Those variables are introduced since branch saturation
            ACM.iqs = 0.0;
            ACM.ids = 0.0;
            ACM.iqr = 0.0;
            ACM.idr = 0.0;
            ACM.psimq = 0.0;
            ACM.psimd = 0.0;
    }

    /* Saturation Model (also considers inverter model) | Could be adapted for PMSM*/
    void collectCurrents(double *x){
        // Generalised Current by Therrien2013
        ACM.izq = x[1]/ACM.Lls + x[3]/ACM.Llr;
        ACM.izd = x[0]/ACM.Lls + x[2]/ACM.Llr;
        ACM.iz = sqrt(ACM.izd*ACM.izd + ACM.izq*ACM.izq);

        if(ACM.iz>1e-8){
            #if SATURATED_MAGNETIC_CIRCUIT
                ACM.psim = sat_lookup(ACM.iz, satLUT);
                ACM.im = ACM.iz - ACM.psim/ACM.LSigmal;
                {
                    ACM.Lm = ACM.psim/ACM.im;
                    ACM.Lmu = ACM.Lm*ACM.Lm/(ACM.Lm+ACM.Llr);
                    ACM.Lmu_inv = 1.0/ACM.Lmu;
                    ACM.alpha = ACM.rr/(ACM.Lm+ACM.Llr);
                    ACM.rreq = ACM.Lmu*ACM.alpha;
                    ACM.Lsigma = (ACM.Lls+ACM.Lm) - ACM.Lmu;
                    ACM.Lm_slash_Lr = ACM.Lm/(ACM.Lm+ACM.Llr);
                    ACM.Lr_slash_Lm = (ACM.Lm+ACM.Llr)/ACM.Lm;
                }
            #else
                ACM.psim = 1.0/(1.0/ACM.Lm+1.0/ACM.Lls+1.0/ACM.Llr)*ACM.iz;
            #endif

            ACM.psimq = ACM.psim/ACM.iz*ACM.izq;
            ACM.psimd = ACM.psim/ACM.iz*ACM.izd;
        }else{
            #if ACMSIMC_DEBUG
                printf("how to handle zero iz?\n");
            #endif
            ACM.psimq = 0;
            ACM.psimd = 0;
        }

        ACM.iqs = (x[1] - ACM.psimq) / ACM.Lls;
        ACM.ids = (x[0] - ACM.psimd) / ACM.Lls;
        ACM.iqr = (x[3] - ACM.psimq) / ACM.Llr;
        ACM.idr = (x[2] - ACM.psimd) / ACM.Llr;    

        // /* Direct compute is ir from psis psir */
        // ACM.iqs = (x[1] - (ACM.Lm/(ACM.Lm+ACM.Llr))*x[3])/ACM.Lsigma;
        // ACM.ids = (x[0] - (ACM.Lm/(ACM.Lm+ACM.Llr))*x[2])/ACM.Lsigma;
        // ACM.iqr = (x[3] - (ACM.Lm/(ACM.Lm+ACM.Lls))*x[1])/(ACM.Lm+ACM.Llr-ACM.Lm*ACM.Lm/(ACM.Lm+ACM.Lls));
        // ACM.idr = (x[2] - (ACM.Lm/(ACM.Lm+ACM.Lls))*x[0])/(ACM.Lm+ACM.Llr-ACM.Lm*ACM.Lm/(ACM.Lm+ACM.Lls));
    }
    void rK5_satDynamics(double t, double *x, double *fx){
        /* argument t is omitted*/

        /* STEP ZERO: collect all the currents: is, ir, iz */
        collectCurrents(x);

        /* STEP ONE: Inverter Nonlinearity Considered */
        ACM.ual; // not CTRL.ual
        ACM.ube; // not CTRL.ube

        /* STEP TWO: electromagnetic model with full flux states in alpha-beta frame */
        fx[1] = ACM.ube - ACM.rs*ACM.iqs; // 这里是反过来的！Q轴在前面！
        fx[0] = ACM.ual - ACM.rs*ACM.ids;
        fx[3] =         - ACM.rr*ACM.iqr + x[4]*x[2]; // 这里是反过来的！Q轴在前面！
        fx[2] =         - ACM.rr*ACM.idr - x[4]*x[3];

        /* STEP THREE: mechanical model */
        // ACM.Tem = ACM.npp*(ACM.Lm/(ACM.Lm+ACM.Llr))*(ACM.iqs*x[2]-ACM.ids*x[3]); // this is not better 
        ACM.Tem = ACM.npp*(ACM.iqs*x[0]-ACM.ids*x[1]);
        fx[4] = (ACM.Tem - ACM.Tload)*ACM.mu_m;
    }
    double one_over_six = 1.0/6.0;
    void rK555_Sat(double t, double *x, double hs){
        double k1[5], k2[5], k3[5], k4[5], xk[5];
        double fx[5];
        int i;

        rK5_satDynamics(t, x, fx); // timer.t,
        for(i=0;i<5;++i){        
            k1[i] = fx[i] * hs;
            xk[i] = x[i] + k1[i]*0.5;
        }
        
        rK5_satDynamics(t, xk, fx); // timer.t+hs/2., 
        for(i=0;i<5;++i){        
            k2[i] = fx[i] * hs;
            xk[i] = x[i] + k2[i]*0.5;
        }
        
        rK5_satDynamics(t, xk, fx); // timer.t+hs/2., 
        for(i=0;i<5;++i){        
            k3[i] = fx[i] * hs;
            xk[i] = x[i] + k3[i];
        }
        
        rK5_satDynamics(t, xk, fx); // timer.t+hs, 
        for(i=0;i<5;++i){        
            k4[i] = fx[i] * hs;
            x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])*one_over_six;
        }
        
        collectCurrents(x);
    }
#elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
    struct SynchronousMachineSimulated ACM;
    void Machine_init(){
        int i;
        for(i=0;i<5;++i){
            ACM.x[i] = 0.0;
        }
        ACM.rpm = 0.0;
        ACM.rpm_cmd = 0.0;
        ACM.rpm_deriv_cmd = 0.0;
        ACM.Tload = 0.0;
        ACM.Tem = 0.0;

        ACM.R  = 0.45;
        ACM.Ld = 4.15*1e-3;
        ACM.Lq = 16.74*1e-3;
        ACM.KE = 0.504; // Vs/rad
        ACM.L0 = 0.5*(ACM.Ld + ACM.Lq);
        ACM.L1 = 0.5*(ACM.Ld - ACM.Lq);

        ACM.Js = 0.06; // Awaya92 using ACM.omg
        ACM.npp = 2;
        ACM.mu_m = ACM.npp/ACM.Js;

        ACM.Ts  = MACHINE_TS;

        ACM.id = 0.0;
        ACM.iq = 0.0;

        ACM.ial = 0.0;
        ACM.ibe = 0.0;

        ACM.ud = 0.0;
        ACM.uq = 0.0;

        ACM.ual = 0.0;
        ACM.ube = 0.0;

        ACM.theta_d = 0.0;
    }
#endif

/* Simple Model */
void rK5_dynamics(double t, double *x, double *fx){
    #if MACHINE_TYPE == INDUCTION_MACHINE
        // electromagnetic model
        fx[2] = ACM.rreq*x[0] - ACM.alpha*x[2] - x[4]*x[3]; // flux-alpha
        fx[3] = ACM.rreq*x[1] - ACM.alpha*x[3] + x[4]*x[2]; // flux-beta
        fx[0] = (ACM.ual - ACM.rs*x[0] - fx[2])/ACM.Lsigma; // current-alpha
        fx[1] = (ACM.ube - ACM.rs*x[1] - fx[3])/ACM.Lsigma; // current-beta

        // mechanical model
        ACM.Tem = ACM.npp*(x[1]*x[2]-x[0]*x[3]);
        fx[4] = (ACM.Tem - ACM.Tload)*ACM.mu_m; // elec. angular rotor speed
        fx[5] = x[4];                           // elec. angular rotor position
    #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
        // electromagnetic model
        fx[0] = (ACM.ud - ACM.R * x[0] + x[2]*ACM.Lq*x[1]) / ACM.Ld; // current-d
        fx[1] = (ACM.uq - ACM.R * x[1] - x[2]*ACM.Ld*x[0] - x[2]*ACM.KE) / ACM.Lq; // current-q

        // mechanical model
        ACM.Tem = ACM.npp*(x[1]*ACM.KE + (ACM.Ld - ACM.Lq)*x[0]*x[1]);
        fx[2] = (ACM.Tem - ACM.Tload)*ACM.mu_m; // elec. angular rotor speed
        fx[3] = x[2];                           // elec. angular rotor position
    #endif
}
void rK555_Lin(double t, double *x, double hs){
    #if MACHINE_TYPE == INDUCTION_MACHINE
        #define NUMBER_OF_STATES 6
    #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
        #define NUMBER_OF_STATES 4
    #endif
    #define NS NUMBER_OF_STATES

    double k1[NS], k2[NS], k3[NS], k4[NS], xk[NS];
    double fx[NS];
    int i;

    rK5_dynamics(t, x, fx); // timer.t,
    for(i=0;i<NS;++i){        
        k1[i] = fx[i] * hs;
        xk[i] = x[i] + k1[i]*0.5;
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<NS;++i){        
        k2[i] = fx[i] * hs;
        xk[i] = x[i] + k2[i]*0.5;
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs/2., 
    for(i=0;i<NS;++i){        
        k3[i] = fx[i] * hs;
        xk[i] = x[i] + k3[i];
    }
    
    rK5_dynamics(t, xk, fx); // timer.t+hs, 
    for(i=0;i<NS;++i){        
        k4[i] = fx[i] * hs;
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6.0;
    }
}


int machine_simulation(){

    // API for explicit access
    #if MACHINE_TYPE == INDUCTION_MACHINE
        // rK555_Lin(CTRL.timebase, ACM.x, ACM.Ts);
        // ACM.ial    = ACM.x[0]; // rK555_Lin
        // ACM.ibe    = ACM.x[1]; // rK555_Lin
        // ACM.psi_al = ACM.x[2]; // rK555_Lin
        // ACM.psi_be = ACM.x[3]; // rK555_Lin

        rK555_Sat(CTRL.timebase, ACM.x, ACM.Ts);
        ACM.ial    = ACM.ids; // rK555_Sat
        ACM.ibe    = ACM.iqs; // rK555_Sat
        ACM.psi_al = ACM.x[2]*ACM.Lm_slash_Lr; // rK555_Sat
        ACM.psi_be = ACM.x[3]*ACM.Lm_slash_Lr; // rK555_Sat

        ACM.rpm    = ACM.x[4] * 60 / (2 * M_PI * ACM.npp);

    #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
        rK555_Lin(CTRL.timebase, ACM.x, ACM.Ts);

        ACM.theta_d = ACM.x[3];
        if(ACM.theta_d > M_PI){
            ACM.theta_d -= 2*M_PI;
        }else if(ACM.theta_d < -M_PI){
            ACM.theta_d += 2*M_PI; // 反转！
        }
        ACM.x[3] = ACM.theta_d;

        ACM.id  = ACM.x[0];
        ACM.iq  = ACM.x[1];
        ACM.ial = MT2A(ACM.id, ACM.iq, cos(ACM.theta_d), sin(ACM.theta_d));
        ACM.ibe = MT2B(ACM.id, ACM.iq, cos(ACM.theta_d), sin(ACM.theta_d));
        ACM.rpm = ACM.x[2] * 60 / (2 * M_PI * ACM.npp);
    #endif

    if(isNumber(ACM.rpm)){
        return false;
    }else{
        printf("ACM.rpm is %g\n", ACM.rpm);
        return true;        
    }
}
void measurement(){
    US_C(0) = CTRL.ual;
    US_C(1) = CTRL.ube;
    US_P(0) = US_C(0);
    US_P(1) = US_C(1);

    #if MACHINE_TYPE == INDUCTION_MACHINE
        IS_C(0) = ACM.ial;
        IS_C(1) = ACM.ibe;
        im.omg = ACM.x[4];
        im.theta_r = ACM.x[5];
    #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
        IS_C(0) = ACM.ial;
        IS_C(1) = ACM.ibe;
        sm.omg = ACM.x[2];
        sm.theta_d = ACM.x[3];
        sm.theta_r = sm.theta_d;
    #endif
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
    #endif

    // 冗余变量赋值
    #if MACHINE_TYPE == INDUCTION_MACHINE
        ;
    #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
        ACM.ud = AB2M(ACM.ual, ACM.ube, cos(ACM.theta_d), sin(ACM.theta_d));
        ACM.uq = AB2T(ACM.ual, ACM.ube, cos(ACM.theta_d), sin(ACM.theta_d));
    #endif
}

int main(){
    
    printf("NUMBER_OF_LINES: %d\n\n", NUMBER_OF_LINES);

    /* Initialization */
    Machine_init();
    CTRL_init();
    acm_init();
    ob_init();

    FILE *fw;
    fw = fopen("algorithm.dat", "w");
    write_header_to_file(fw);

    /* MAIN LOOP */
    clock_t  begin, end;
    begin = clock();
    int _; // _ for the outer iteration
    int dfe=0; // dfe for down frequency execution
    for(_=0;_<NUMBER_OF_LINES;++_){

        /* Command and Load Torque */
        // cmd_fast_speed_reversal(CTRL.timebase, 5, 5, 1500); // timebase, instant, interval, rpm_cmd
        cmd_fast_speed_reversal(CTRL.timebase, 5, 5, 100); // timebase, instant, interval, rpm_cmd
        // ACM.Tload = 5 * sign(ACM.rpm); 

        ACM.Tload = 0 * sign(ACM.rpm); // No-load test
        // ACM.Tload = ACM.Tem; // Blocked-rotor test

        /* Simulated ACM */
        if(machine_simulation()){ 
            printf("Break the loop.\n");
            break;
        }

        if(++dfe==DOWN_FREQ_EXE){
            dfe = 0;

            /* Time */
            CTRL.timebase += TS;

            measurement();

            observation();

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
    #if MACHINE_TYPE == INDUCTION_MACHINE
        // no space is allowed!
        // fprintf(fw, "x0,x1,x2,x3,rpm,uMs_cmd,uTs_cmd,iMs_cmd,iMs,iTs_cmd,iTs,psi_mu_al,tajima_rpm\n");
        // fprintf(fw, "$x_0$,$x_1$,$x_2$,$x_3$,Speed [rpm],$u_{Ms}^*$,$u_{Ts}^*$,$i_{Ms}^*$,$i_{Ms}$,$i_{Ts}^*$,$i_{Ts}$,$\\psi_{\\alpha\\mu}$,tajima_rpm\n");
        // fprintf(fw, "ACM.x[0],ACM.x[1],ACM.x[2],ACM.x[3],ACM.x[4],ACM.Tem,CTRL.uMs_cmd,CTRL.uTs_cmd,CTRL.iMs_cmd,CTRL.iMs,CTRL.iTs_cmd,CTRL.iTs,ob.psi_mu_al,ob.tajima.omg*RAD_PER_SEC_2_RPM,ACM.rpm\n");
        fprintf(fw, "ACM.x[0],ACM.x[1],ACM.x[2],ACM.x[3],ACM.x[4],ACM.Tem,ACM.ual,ACM.ube,ACM.rpm_cmd,e_omega\n");
    #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
        // no space is allowed!
        fprintf(fw, "x0,x1,x2,x3,uMs_cmd,uTs_cmd,iMs_cmd,iMs,iTs_cmd,iTs\n");
    #endif

    {
        FILE *fw2;
        fw2 = fopen("info.dat", "w");
        fprintf(fw2, "TS,DOWN_SAMPLE\n");
        fprintf(fw2, "%g, %d\n", TS, DOWN_SAMPLE);
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
            #if MACHINE_TYPE == INDUCTION_MACHINE
                // 数目必须对上，否则ACMAnimate会失效，但是不会影响ACMPlot
                fprintf(fw, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
                        ACM.x[0],ACM.x[1],ACM.x[2],ACM.x[3],ACM.x[4],ACM.Tem,
                        ACM.ual,ACM.ube,ACM.rpm_cmd,ACM.rpm_cmd-ACM.x[4]*RAD_PER_SEC_2_RPM
                        );
            #elif MACHINE_TYPE == SYNCHRONOUS_MACHINE
                fprintf(fw, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
                        ACM.x[0], ACM.x[1], ACM.x[2], ACM.x[3],
                        CTRL.uMs_cmd, CTRL.uTs_cmd, CTRL.iMs_cmd, CTRL.iMs, CTRL.iTs_cmd, CTRL.iTs
                        );
            #endif
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
    // but it's false if x is a NaN (1.#QNAN0).
    return (x == x); 
    // see https://www.johndcook.com/blog/IEEE_exceptions_in_cpp/ cb: https://stackoverflow.com/questions/347920/what-do-1-inf00-1-ind00-and-1-ind-mean
}



