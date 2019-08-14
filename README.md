# ACMSIMC_TUT
> AC Machine Simulation in C (Tutorial Version)

I have used C to simulate motor control and adptive observers for over 4 years now.
This is a tutorial for those who hate using Simulink to simulate motor control.

The numerical integration method is currently RK4, which is quite enough. 
DoPri54 will be included in future version (including stiffness detection and variable step numerical integration).

**Introduction on current branches:**
- [IM & PMSM] master: the branch with full feature (for most of times).
- [IM] vvvf: the skeleton with induction motor simulation and VVVF control.
- [IM] foc: field oriented control (direct/indirect) with basic sensorless control.
- [IM] animate: test the feature of waveform (results) animation.
- [PMSM] pmsm: id=0 control for interior permanent magnet synchronous motor.
- [IM] vi_decouple: add voltage-current decoupling circuit for improved control performance during high speed reversal.
- [IM] mras: model reference adaptive system based sensorless control (my 2017-Chen.Huang-Online paper).
- [IM] \_femm: (This branch does not really belong here but I don't want to create a new repository for it...) It is about the design of the induction motor using free softwares as well as fitting the design to the equivalent circuit parameters for further control simulation.
- [IM] saturation: include iron core saturation effect into the induction motor model simulation.

> If you speak Chinese, I have a series of dedicated tutorial videos on how to make this thing work from the ground up.
> Please take a look at this link to [知乎](https://zhuanlan.zhihu.com/p/64445558).
> In fact, now you can check out [my personal page](https://horychen.github.io) for a list of tutorial videos.
