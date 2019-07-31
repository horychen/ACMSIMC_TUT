# ACMSIMC_TUT
> AC Machine Simulation in C (Tutorial Version)

## Introduction
I have used C to simulate motor control and adptive observers for over 4 years now.
This is a tutorial for those who hate using Simulink to simulate ac motor control.

## Numerial Methods
The numerical integration method is currently RK4, which is quite enough. 
DoPri54 will be included in future version.

## Introduction to Current Branches:
- [IM] master: vvvf branch plus other utility features.
- [IM] vvvf: the skeleton with induction motor simulation and VVVF control.
- [IM] foc: field oriented control (direct/indirect) with basic sensorless control.
- [IM] animate: test the feature of waveform animation.
- [PMSM] pmsm: id=0 control for interior permanent magnet synchronous motor.
- [IM] mras: model reference adaptive system based sensorless control (my 2017-Chen.Huang-Online paper).

## Dependency under Windows
- Anaconda 3 (you should be able to call python from cmd.exe)
- MinGW (you should be able to call gcc from cmd.exe)
- FEMM (femm.info, you do not need this if you are not interested in Finite Element Analysis and motor design)
  - PyFEMM (use pip to install this)

## Dependency under Linux
- Not tested yet. (I believe that linux users should be able to figure it out...)

## Video Tutorials
For unfamiliated users, I have been creating video tutorials. However, they are currently in Chinese. 
In near future, tge Engligh version will be produced.

> If you speak Chinese, I have a dedicated tutorial video on how to make this thing work from the ground up.
> Please take a look at this link to [知乎](https://zhuanlan.zhihu.com/p/64445558).
> In fact, now you can check out [my personal page](https://horychen.github.io) for a list of tutorial videos.
