/*
 *  inverter.h
 *
 *  By: Hory Chen
 *  Email: horychen@qq.com
 *  Created on: June 27th 2017
 *
 *  Tips: Do not declare any variables in h files.
 */ 

#ifndef INVERTER_H
#define INVERTER_H

void InverterNonlinearity_SKSul96(double ual, double ube, double ial, double ibe);

// Inverter Model
#define _Vce0  1.8 // V
#define _Vd0   1.3 // V
#define _Udc   300 // V
#define _Toff  0.32e-6 // sec
#define _Ton   0.15e-6 // sec
#define _Tdead 3.30e-6 // sec
#define _Tcomp 0.0*(_Tdead+_Ton-_Toff) //3.13e-6; //8e-6; // 过补偿 //3.13e-6; // 只补偿死区

#endif
