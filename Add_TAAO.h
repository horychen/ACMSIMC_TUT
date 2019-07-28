/* 
 *  Add_TAAO.h 
 * 
 *  by hory chen 
 *  Email: horychen@qq.com 
 *  Created on: Nov. 26th 2016 
 * 
 *  Warning: Do not declare any variables in h files, use extern instead. 
 */ 
 
/* Consistent with DSP Programs Begin 
 * Except for type double             */ 
#ifndef ADD_TAAO_H 
#define ADD_TAAO_H 
#if OBSERVER_APPLIED == TOTALLY_ADAPTIVE_OBSERVER 
 
 
#define TAAO_FLUX_COMMAND_VALUE (1.2) // 1.2 
#define TAAO_FLUX_COMMAND_ON (false) 
#define TIME_DIVISION_MULTIPLEXING (false) 
 
#define KVM_VALUE (300) // (100) 
#define KMEC_VALUE (0) 

/* 第三轮调试（无速度传感器运行、供稿IEMDC） 
 * */ 
#define M1 (0.01*0) //0.01//(0.05) // 转速辨识大震荡的罪魁祸首竟然是中频注入？ 
#define OMG1 (24*2*M_PI) //(0.5*2*M_PI) 
#define M3 (0.01*0) //(0.02) 
#define OMG3 (10*2*M_PI) /* 10Hz竟然可以？！低频注入仍然不行，只能猜测25Hz可能刚好被逆非补偿给削弱了？不！就是越高频稳态辨识误差越大！ */ 
#define M2 (0.0) // 0.1 
#define OMG2 (6*2*M_PI) //(1*2*M_PI) 
 
/* 四参数辨识（有偏） */ 
#define GROUP_M 0*0.7 
#define GROUP_T 10
#define GAMMA_LSIGMA 0.0 //0.1 //0.5 
#define GAMMA_RS     0*200 //(GROUP_M+0*GROUP_T)*(215) // rs小12.26~12.275大rs@TL=6Nm, //13.5@TL=20Nm //1.5 //2.5 // 5  
#define GAMMA_LEQ    0*GROUP_M*2.5*10*1e1 //13 for only Leq Id. // 10, 14.01~15就太快了，（加入饱和特性后）转速暂态时导致系统不稳定！是在150rpm下不稳定，用gamma4=14.01可以体现出来，14.00也不行。 
#define GAMMA_RREQ   0*1000 //0*GROUP_T*122*1e1 //*2.5 
#define GAMMA_OMEGAR GROUP_T*20520 // 6e5 for decoupled IFOC??? // 6e4 for IOFLC 500rpm //6e3 for IFOC 50,80rpm //2e3 //(7e2) // 1e3 for no ripple; 5e4 for fast response; 5e5 for good simulation. With kVM=4000 
#define GAMMA_OMEGAR_P 0// (100) 稳态转速自振荡 // 后面的都是IFOC加入解耦控制之前的// 600 for 0.8 flux command; 100 for 1.2 flux command 

 
 
#define AB2M(A, B, COS, SIN)  ( (A)*COS  + (B)*SIN ) 
#define AB2T(A, B, COS, SIN)  ( (A)*-SIN + (B)*COS ) 
#define MT2A(M, T, COS, SIN)  ( (M)*COS - (T)*SIN ) 
#define MT2B(M, T, COS, SIN)  ( (M)*SIN + (T)*COS ) 
 
 
/* Macro for External Access Interface */ 
#define US(X) im.us[X] 
#define IS(X) im.is[X] 
#define US_C(X) im.us_curr[X] 
#define IS_C(X) im.is_curr[X] 
#define US_P(X) im.us_prev[X] 
#define IS_P(X) im.is_prev[X] 
 
#define OB_FLUX(X) ob.xCM[X] //ob.xCM[0]*IM.Lr_slash_Lm; 
#define OB_OMG ob.taao_omg 
#define OB_ALPHA ob.taao_alpha 
#define OB_LMU ob.taao_Leq 
#define OB_LSIGMA ob.taao_Lsigma 
#define OB_RS ob.taao_rs 
#define OB_RREQ ob.taao_rreq 
 
#define OB_LMU_INV (1.0/OB_LMU) 
#define OB_TLOAD IM.Tload 
 
struct InductionMachine{ 
    double us[2]; 
    double is[2]; 
    double us_curr[2]; 
    double is_curr[2]; 
    double us_prev[2]; 
    double is_prev[2]; 
 
    double Js; 
    double Js_inv; 
    double Lm; 
    double Lm_inv; 
 
    double Lls; 
    double Llr; 
    double Ls; 
    double Lr; 
    double sigma; 
    double rs; 
    double rr; 
    double Tr; 
    double alpha; 
    double Lsigma; 
    double Leq; 
    double Leq_inv; 
    double rreq; 
 
    double npp; 
    double omg; 
    double omg_curr; 
    double omg_prev; 
 
    double theta_r; 
}; 
extern struct InductionMachine im; 
 
struct ObserverControl{ 
    double k_VM; 
    double k_CM; 
    double k_mec; 
 
    double xVM[2]; // rotor flux of VM 
    double xCM[2]; // rotor flux of CM 
    double statorFlux[2]; // stator flux of VM 
    double varPi; // cascaded mechanical model 
    double Tem; 
 
    double mismatch[3]; 
    // double pmismatch[2]; 
    double mismatch_mod; 
    // double mismatch_phase; 
    double error[2]; 
 
    double gamma1; // rs 
    double gamma2; // Lsigma 
    double gamma3; // rreq 
    double gamma4; // Leq 
    double gamma5; // omg 
    double gamma5P; // omg 
 
    double timebase; 
    double Ts; 
 
    double taao_b[2]; 
    double taao_rs; 
    double taao_Lsigma; 
    double taao_rreq; 
    double taao_Leq; 
 
    double taao_alpha; 
 
    double taao_omg; 
    double taao_omg_integralPart; 
    double taao_speed; 
    double taao_invJs; 
    double taao_TLslashJs; 
 
    double omega_e; 
 
    double taao_flux_cmd; 
    int taao_flux_cmd_on; 
 
    double taao_VMFlux_fed; 
    double taao_CMHpfFlux_fed; 
 
    double taao_thetaA; 
    double taao_cos_thetaA; 
    double taao_sin_thetaA; 
 
    // double uMs; // for Tsuji's CM whose magnitude will never fail, which is worth a try. 
    // double uTs; 
    // double iMs; 
    // double iTs; 
}; 
extern struct ObserverControl ob; 
 
void acm_init(); 
void ob_init(); 
void observation(); 
double TAAO_FluxModulusCommand(); 
 
#endif 
#endif 
/* Consistent with DSP Programs End */ 
 
