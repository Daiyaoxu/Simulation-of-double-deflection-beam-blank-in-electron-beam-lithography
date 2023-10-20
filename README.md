# Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography

**工作的具体结果可以查看PPT 文件**

This project focuses on the description of Lens excitation effects and conjugate blanking optics to understand the effect of the electron beam blanking on the electron beam through the trajectory equation.We can ascertain the impact of the Beam Blank design on exposure processing through simulation calculations.

Electron optics, with its inherent optical properties, suggests that deflection might not necessarily impact the electron beam's spot. Under specific conditions, where an electron lies at a certain off-axis position with a particular lateral spatial velocity, it can still converge at a designated focus point, even disregarding potential lens aberrations like spherical and chromatic aberrations. However, any spatial deviation of the electron beam can lead to a beam spot dispersion, reminiscent of the imaging principles seen in optical lenses. Such behavior is typically articulated using concepts of phase space distribution and emittance.

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/c4cef2ad-4957-4ff0-b044-97791df76124" width="800px">


The primary function of the Beam Blank system is to deviate the electron beam by applying voltage to a plate, steering it outside of an aperture to cease the beam flow. During this deflection, inherent limitations of the circuitry and anomalies during the capacitive charging and discharging, such as overflow and latency effects, can lead to voltage surges and excessive feedback. Such dynamics can disturb the phase space distribution of the electron beam. Any resulting deviation in the beam spot on the sample could subsequently distort the exposure lines. For intricate micro and nanostructures, these distortions could be detrimental, even potentially causing unintended connections between two distinct structures. As scanning speeds amplify, these effects became worse. To facilitate high-frequency scanning speeds, it becomes imperative to employ faster response times and a more stable voltage control system to mitigate these influences. To tackle these challenges, we use a double-deflection blanker. This optical design helps reduce the errors that might happen. So, to increase the scanning speed without causing error. The structure of double-deflection blanker is shown below.

![image](https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/a661ccf4-6a19-4753-bb90-531377d8da54)

Gesley, M. (1993). An electron optical theory of beam blanking. Review of Scientific Instruments, 64(11), 3169–3190. https://doi.org/10.1063/1.1144326


### For the electronic offset calculation process

![image](https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/e2df3a5c-3162-4485-aec1-1028d22e2088)

$$\Delta s_i$$是电子在离开beam blank末端的位移，其中M为物镜的放大倍数，a(t)时间相关加速度$$-\frac{e_0 V(t)}{m d}$$，τ和λ为用于积分的临时变量，通过最终的横向位移和横向速度可以计算出





**Fundamental equation **

$$
\Delta s_i=\frac{1}{M}\left[\int_0^{t_0}\left(\int_0^\lambda a(\tau) d \tau\right) d \lambda-\frac{1}{2} t_0 \int_0^{t_0} a(\tau) d \tau\right]
$$

### Error analysis of a single deflector


<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/4e3daed9-fb45-44f1-afb4-4a8daca09341" width="500px">


<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/ba6ca019-8fae-453f-b09d-40f010494a9d" width="500px">


<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/8ca08e89-4c77-479c-8cf9-be8e0d77b94a" width="500px">


<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/c9779f57-c00f-4939-8459-ba331bdaa7c3" width="500px">

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/a9cdf4fc-1895-4374-a276-eae5af5987a8" width="500px">

