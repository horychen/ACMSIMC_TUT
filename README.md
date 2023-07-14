# ACMSIMC_TUT
> AC Machine Simulation in C (Tutorial Version)

## Introduction
I have used C to simulate motor control and adptive observers for over 4 years now.
This is a tutorial for those who hate using Simulink to simulate ac motor control.
The benefit is that you can directly reuse the codes in DSP based motor drive.

## Numerial Methods
The numerical integration method is currently RK4, which is quite enough. 
DoPri54 will be included in future version (including stiffness detection and variable step numerical integration).

## Introduction to Current Branches (in time order):
- [IM] vvvf: the skeleton with induction motor simulation and VVVF control.
- [IM] foc: field oriented control (direct/indirect) with basic sensorless control.
- [IM] animate: test the feature of waveform animation like an oscilloscope.
- [PMSM] pmsm: id=0 control for interior permanent magnet synchronous motor.
- [IM] vi_decouple: add voltage-current decoupling circuit for improved control performance during high speed reversal.
- [IM] mras: model reference adaptive system based sensorless control (my 2017-Chen.Huang-Online paper).
- [IM] \_femm: (This branch does not really belong here but I don't want to create a new repository for it...) It is about the design of the induction motor using free softwares as well as fitting the design to the equivalent circuit parameters for further control simulation.
- [IM] saturation: include iron core saturation effect into the induction motor model simulation.
- [Both] inverter_model: simple inverter modeling based on the paper 1996-Choi.Sul-Inverter.
- [Both] **master: contain all the features of the branches mentioned above. The master branch is not updated anymore from this point, because I realized that having both IM and PMSM codes in one place is a silly idea.**
- [PMSM] commissioning_pmsm: (under developing) self-commissioning procedure for permanent magnet motor.
- [PMSM] eemf: sensorless control based on extended emf method proposed by Zhiqian Chen et al. (2003). Sensorless open loop works but closed-loop has problems.
- [PMSM] hfsi-yoon: pulsating type square wave high frequency signal injection method based sensorless control for motor with saliency.
- [PMSM] commissioning_spmsm: self-commissioning procedure for surface mounted permanent magnet synchronous motor.
- [PMSM] spmsm_scvm_harnerfors06: statically compensated voltage model based sensorless control for SPMSM. This branch also includes new dynamic anti-windup and incremental PID regulator.

(I recommend to use GitKraken to view this repository, so you know clearly which branch is the newest one. Note that master is obsolete. Sorry about my poor code management skills T-T.)

## Visualization
- The plots are made using package matplotlib. 
    - In branch animate, I tested the feature of waveform animation with matplotlib.
- You can also try browser based libraries, e.g., plotly_express.

## Dependency under Windows
- Anaconda 3 (you should be able to call python from cmd.exe)
- MinGW (you should be able to call gcc from cmd.exe)
- FEMM ([femm.info](http://www.femm.info/wiki/HomePage), you do not need this if you are not interested in Finite Element Analysis and motor design)
  - PyFEMM (a wrapper for FEMM API; use pip to install this)
- Editor (optional): Sublime Text (preferred) or Visual Studio Code.

## Dependency under Linux
- FEMM is a Windows-only FEA software. 
    - Alternative is ElmerFEM for Linux, but it is poorly documented. Please DO NOT try it out unless you are a ドM.
- Others are not tested yet. (Linux users should be able to figure it out...)

## Compile

### Compile in Sublime Text
Create a file named "**!C_GCC.sublime-build**" in this folder "C:\Users\yourUserName\AppData\Roaming\Sublime Text 3\Packages\User" with the following content:
```json
    {
        "working_dir": "$file_path",
        "cmd": "gcc -Wall $file_name -o $file_base_name",
        "file_regex": "^(..[^:]*):([0-9]+):?([0-9]+)?:? (.*)$",
        "selector": "source.c",
        "variants": 
        [
            {   
            "name": "ACMSIMC_TUT",
                // "shell_cmd": "gcc $file controller.c observer.c -L. -o $file_base_name && start cmd /c \"${file_path}/${file_base_name}\""
                "shell_cmd": "gcc $file commissioning.c inverter.c controller.c observer.c -L. -o $file_base_name && start cmd /c \"${file_path}/${file_base_name}\""
            }
        ]
    }
```
Press down Ctrl + Shift + B to select "ACMSIMC_TUT" and next time you only need to press down Ctrl + B to build.

### Compile in cmd.exe
Change directory to where main.c is and open the cmd.exe there and type:
```batch
gcc main.c commissioning.c inverter.c controller.c observer.c -I. -L. -o main
```
Then, double click main.exe to run.
If you have python added to system path, plots will also pop up.


## Video Tutorials
For unfamiliar users, I have been creating video tutorials. However, they are currently in Chinese. 
> *See [my Bilibili space](https://space.bilibili.com/7132537) (ACTIVELY UPDATED).*
> 
> Also take a look at this link to [知乎](https://zhuanlan.zhihu.com/p/64445558) (not acvtively updated).
> 
> Also check out [my personal page](https://horychen.github.io) for a list of tutorial videos (not acvtively updated). 

In near future, the English version will be brought about. But I would like to do it in high quality so I procrastinate.

<!-- **[Important Update]** I made a video in English to explain how to use C codes to simulate sensorless drive using the method from one of my papers ("Resistances and Speed Estimation in Sensorless Induction Motor Drives Using a Model with Known Regressors"). Here is the [link]() to it.
Jiahao
2020/02/10
 -->
