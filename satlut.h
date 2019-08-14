#ifndef SATURATION_LOOPUPTABLE_H
#define SATURATION_LOOPUPTABLE_H


#define IZ_STEP 0.01
#define IZ_STEP_INV 100

float sat_lookup(double iz, float satLUT[]);

extern float satLUT[];




#endif