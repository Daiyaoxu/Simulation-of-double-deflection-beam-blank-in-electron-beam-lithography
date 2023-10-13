# Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography

This project focuses on the description of Lens excitation effects and conjugate blanking optics to understand the effect of the electron beam blanking on the electron beam through the trajectory equation.We can ascertain the impact of the Beam Blank design on exposure processing through simulation calculations.

由于电子光学产生的 光学特性，我们认为偏转本身并不一定会对电子束斑产生影响，就像下图第一张所展示的那样。在不考虑磁透镜的畸变（球差，慧差，像散等） 当电子所在的旁轴位置和横向空间速度符合某种关系时，通过透镜对汇聚，电子依然会汇聚在某一个点。但如果这时候这种如果电子束在空间中位置偏移，如下图第二张所示，像的位置发生偏移，导致束斑扩散。这种现象与光学透镜的成像原理相似。在这里我们一般用相空间分布和发射度这个概念描述。粒子的这种信息，也可以叫Twiss参数。

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/c4cef2ad-4957-4ff0-b044-97791df76124" width="800px">

**对于电子偏移量的计算流程**


<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/ba6ca019-8fae-453f-b09d-40f010494a9d" width="800px">

$$
\Delta s_i=\frac{1}{M}\left[\int_0^{t_0}\left(\int_0^\lambda a(\tau) d \tau\right) d \lambda-\frac{1}{2} t_0 \int_0^{t_0} a(\tau) d \tau\right]
$$

<img src="https://github.com/Daiyaoxu/Simulation-of-double-deflection-beam-blank-in-electron-beam-lithography/assets/130887176/c2202cfd-fdf1-4a20-9464-128c90f9e203" width="800px">


