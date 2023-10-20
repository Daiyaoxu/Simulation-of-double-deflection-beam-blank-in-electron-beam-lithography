# Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography

**The exact results of the work can be viewed in a PowerPoint file**

This project focuses on the description of Lens excitation effects and conjugate blanking optics to understand the effect of the electron beam blanking on the electron beam through the trajectory equation.We can ascertain the impact of the Beam Blank design on exposure processing through simulation calculations.

Electron optics, with its inherent optical properties, suggests that deflection might not necessarily impact the electron beam's spot. Under specific conditions, where an electron lies at a certain off-axis position with a particular lateral spatial velocity, it can still converge at a designated focus point, even disregarding potential lens aberrations like spherical and chromatic aberrations. However, any spatial deviation of the electron beam can lead to a beam spot dispersion, reminiscent of the imaging principles seen in optical lenses. Such behavior is typically articulated using concepts of phase space distribution and emittance.

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/c4cef2ad-4957-4ff0-b044-97791df76124" width="800px">


The primary function of the Beam Blank system is to deviate the electron beam by applying voltage to a plate, steering it outside of an aperture to cease the beam flow. During this deflection, inherent limitations of the circuitry and anomalies during the capacitive charging and discharging, such as overflow and latency effects, can lead to voltage surges and excessive feedback. Such dynamics can disturb the phase space distribution of the electron beam. Any resulting deviation in the beam spot on the sample could subsequently distort the exposure lines. For intricate micro and nanostructures, these distortions could be detrimental, even potentially causing unintended connections between two distinct structures. As scanning speeds amplify, these effects became worse. To facilitate high-frequency scanning speeds, it becomes imperative to employ faster response times and a more stable voltage control system to mitigate these influences. To tackle these challenges, we use a double-deflection blanker. This optical design helps reduce the errors that might happen. So, to increase the scanning speed without causing error. The structure of double-deflection blanker is shown below.

![image](https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/a661ccf4-6a19-4753-bb90-531377d8da54)

Gesley, M. (1993). An electron optical theory of beam blanking. Review of Scientific Instruments, 64(11), 3169–3190. https://doi.org/10.1063/1.1144326

** For the electronic offset calculation process**

![image](https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/e2df3a5c-3162-4485-aec1-1028d22e2088)


**Fundamental equation **

$$
\Delta s_i=\frac{1}{M}\left[\int_0^{t_0}\left(\int_0^\lambda a(\tau) d \tau\right) d \lambda-\frac{1}{2} t_0 \int_0^{t_0} a(\tau) d \tau\right]
$$

s is the displacement of the electron at the end of leaving the beam blank, where M is the magnification of the objective lens, and a(t) the time-dependent acceleration, where

$$a(t)=-\frac{e_0 V(t)}{m d}$$

τ and λ are temporary variables used for integration, and the exposure error can be calculated from the final transverse displacement and transverse velocity

### Error analysis of a single deflector

With actual circuit measurements, we can obtain a specific representation of the step function, as shown in the figure on the left. The whole step response process is about 130ns. there will be 5% over response. By fitting we can get a specific expression for the step function. The transfer function of the system can be obtained by calculating the step response function. This process can help us for the futher beam blank system design process.


<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/4e3daed9-fb45-44f1-afb4-4a8daca09341" width="500px">

After obtaining the voltage data, we can follow the structure of the device and the following figure shows the basic structure of the beam blank we calculated.

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/ba6ca019-8fae-453f-b09d-40f010494a9d" width="500px">

Single deflector pole plate beam blank, end transverse deflection displacement magnitude, middle arrow site is beam open state.

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/8ca08e89-4c77-479c-8cf9-be8e0d77b94a" width="500px">

By calculating the inverse extension of the trajectory, we can calculate the offset error of the image plane, and by multiplying this error by the magnification of the objective lens, we can calculate the error under the actual electron-beam lithography, and we can see that the error mainly occurs when the deflection voltage is turned on and off, with an offset of about 20 nm.

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/c9779f57-c00f-4939-8459-ba331bdaa7c3" width="500px">

Subsequently, through the Gaussian distribution function, we can practically restore the results of the electron beam lithography, as shown in the following figure, the left side is the ideal target pattern at a scanning frequency of 5Mhz, and the right figure is the calculated results based on the computational simulation.

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/a9cdf4fc-1895-4374-a276-eae5af5987a8" width="500px">

This part of the detailed process can be viewed in the PowerPoint file.
