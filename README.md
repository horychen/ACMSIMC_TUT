# ACMSIMC_TUT
> AC Machine Simulation in C (Tutorial Version)

I have used C to simulate motor control and adptive observers for 4 years now.
This is a tutorial for those who hate using Simulink to simulate motor control.

The numerical integration method is currently RK4, which is quite enough. 
DoPri54 will be included in future revision.

**Introduction on current branches:**
- master: vvvf branch plus other utility features.
- vvvf: the skeleton with induction motor simulation and VVVF control.
- foc: field oriented control (direct/indirect) with basic sensorless control.
- animate: test the feature of waveform animation.

> If you speak Chinese, I have a dedicated tutorial video on how to make this thing work from the ground up.
> Please take a look at this link to 知乎: https://zhuanlan.zhihu.com/p/64445558
> In fact, now you can check out my personal page for a list of tutorial videos at [My Page](https://horychen.github.io)
