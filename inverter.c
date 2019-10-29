#include "ACMSim.h"

#if INVERTER_NONLINEARITY

void InverterNonlinearity_SKSul96(double ual, double ube, double ial, double ibe){
    double ua,ub,uc;
    double ia,ib,ic;
    double Udist;
    double TM;
    double Rce=0.04958, Rdiode=0.05618;

    TM = _Toff - _Ton - _Tdead + _Tcomp; // Sul1996
    Udist = (_Udc*TM*TS_INVERSE - _Vce0 - _Vd0) / 6.0; // Udist = (_Udc*TM/1e-4 - _Vce0 - _Vd0) / 6.0;
    // Udist = (_Udc*TM*TS_INVERSE) / 6.0;
    // Udist = 0.0;

    ia = SQRT_2_SLASH_3 * (       ial                              );
    ib = SQRT_2_SLASH_3 * (-0.5 * ial - SIN_DASH_2PI_SLASH_3 * ibe );
    ic = SQRT_2_SLASH_3 * (-0.5 * ial - SIN_2PI_SLASH_3      * ibe );

    /* compute in abc frame */
    // ua = SQRT_2_SLASH_3 * (       ual                              );
    // ub = SQRT_2_SLASH_3 * (-0.5 * ual - SIN_DASH_2PI_SLASH_3 * ube );
    // uc = SQRT_2_SLASH_3 * (-0.5 * ual - SIN_2PI_SLASH_3      * ube );
    // ua += Udist * (2*sign(ia) - sign(ib) - sign(ic)) - 0.5*(Rce+Rdiode)*ia;
    // ub += Udist * (2*sign(ib) - sign(ic) - sign(ia)) - 0.5*(Rce+Rdiode)*ib;
    // uc += Udist * (2*sign(ic) - sign(ia) - sign(ib)) - 0.5*(Rce+Rdiode)*ic;
    // UAL_C_DIST = SQRT_2_SLASH_3      * (ua - 0.5*ub - 0.5*uc); // sqrt(2/3.)
    // UBE_C_DIST = 0.70710678118654746 * (         ub -     uc); // sqrt(2/3.)*sin(2*pi/3) = sqrt(2/3.)*(sqrt(3)/2)

    /* directly compute in alpha-beta frame */
    // CHECK the sign of the distortion voltage!
    // Sul把Udist视为补偿的电压（假定上升下降时间都已经知道了而且是“补偿”上去的）
    UAL_C_DIST = ual + sqrtf(1.5)*Udist*(2*sign(ia) - sign(ib) - sign(ic)) - 0.5*(Rce+Rdiode)*ial;
    UBE_C_DIST = ube + 3/sqrtf(2)*Udist*(             sign(ib) - sign(ic)) - 0.5*(Rce+Rdiode)*ibe; 

}


#endif