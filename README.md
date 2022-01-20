# Numerical Simulation of 1-D Sod Shock Tube

**The MATLAB codes for the realization of the numerical simulation of 1-D Sod Shock Tube**

## Information

***Author:*** pkufzh (Small Shrimp)

***Course:*** Fundamentals of Computational Fluid Dynamics (CFD)

***Submit:*** 2021/12/28

## Description

All the developed *MATLAB* codes are saved under the ***Codes*** folder.

### Main Program

- <font color = 'red'>**Program_Sod_Shock_Tube_Main**</font>
  - **Main Program:** Numerical simulation of 1-D compressible flow (Sod Shock Tube)

### Attached Function Modules

***Important Note:*** Please ensure the following files are placed in the same folder with the main program!

- <font color = 'blue'>**Flux_Vect_Split_Common.m**</font>
  - **Flux Vector Splitting (FVS)** with different methods *(Optional)*
  - **Steger-Warming (S-W)**
    - **Lax-Friedrichs (L-F)** 
    - **van Leer**
    - **Liou-Steffen (Advection Upstream Splitting Method, AUSM)**
  
- <font color = 'blue'>**Flux_Diff_Split_Common.m**</font>
  - **Flux Difference Splitting (FDS)** with different methods *(Optional)*
  - ***Roe*** Scheme
  
- <font color = 'blue'>**Diff_Cons_Common.m**</font>
  - Calculate the  flux difference $ \frac{\partial \mathbf{F}}{\partial x} $ from positive flux $ \frac{\partial \mathbf{F}^{+}}{\partial x} $ and negative flux  $ \frac{\partial \mathbf{F}^{-}}{\partial x} $ with conservation form through **FVS** or **FDS**, i.e. 

$$
\mathbf{F}_{j + \frac{1}{2}} = \mathbf{F}_{j + \frac{1}{2} L}^{+} + \mathbf{F}_{j + \frac{1}{2} R}^{-}
$$

  - **Shock Capturing Methods** *(Optional)*
    - **(TVD) *Total Variation Diminishing*** Scheme with **van Leer Limiter**
    - **(NND, *H. X. Zhang*) *Non-oscillatory, Non-free-parameters Dissipative Difference*** Scheme
    - **(Original WENO, 5 order, *Jiang & Shu*) *Weighted Essentially Non-oscillatory Scheme*** Scheme
  - **First Level Upwind Schemes** *(Optional)*
    - **1 order** (2 points)
    - **2 order** (3 points)
    - **3 order** (4 points with *bias*)
    - **5 order** (6 points with *bias*)
  - *Note:* All the **upwind schemes** used in this program had been converted into the conservative form.
  
- <font color = 'blue'>**Cal_Minmod.m**</font>
  - **Cal. minmod(a, b)**
  - *Sign Definition*

$$
\operatorname{minmod}(a,b) = \frac{1}{2}\left[\operatorname{sgn}(a) + \operatorname{sgn}(b)\right]\cdot\operatorname{min}\left(\left|a\right|, \left|b\right|\right)
$$

​		where $ a $, $ b $ is 1-Dimensional array with same length.

- <font color = 'blue'>**Plot_Props.m**</font>
  - Plot the properties of fluid with preset and uniform axis coordinates

### Exact Riemann Solution: Referred Functions by Gogol (2021)

***Reference:*** *Gogol* (2021). Sod Shock Tube Problem Solver [Click to the Website](https://www.mathworks.com/matlabcentral/fileexchange/46311-sod-shock-tube-problem-solver), From MATLAB Central File Exchange. Retrieved December 28, 2021. The main codes were developed by the original author.

- <font color = 'blue'>**analytic_sod.m**</font>
  - Solve Sod's Shock Tube problem using exact Riemann solution
  - [Reference Page](http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html)

- <font color = 'blue'>**sod_func.m**</font>
  - Define functions to be used in *analytic_sod.m*
  - Initial conditions

- <font color = 'blue'>**sod_demo.m**</font>
  - A demo script file to show the use of *analytic_sod.m*

***Note:***

- The above *MATLAB* codes were tested and passed on *MATLAB R2021b*, Windows 64-bit system. 
- For the vector images saved in the ***Paper.pdf*** are large, loading may be slow. Thanks for your patient waiting!!!

## Reference

1. Gogol (2021). Sod Shock Tube Problem Solver ([[Click to the Website](https://www.mathworks.com/matlabcentral/fileexchange/46311-sod-shock-tube-problem-solver)), *MATLAB Central File Exchange.* Retrieved December 28, 2021.
2. Steger, J. L., \& Warming, R. F. (1981). Flux vector splitting of the inviscid gasdynamic equations with application to finite-difference methods. *Journal of computational physics*, 40(2), 263-293.
3. Van Leer, B. (1997). Flux-vector splitting for the Euler equation. In Upwind and high-resolution schemes (pp. 80-89). *Springer*, Berlin, Heidelberg.
4. Liou, M. S., \& Steffen Jr, C. J. (1993). A new flux splitting scheme. *Journal of Computational physics*, 107(1), 23-39.
5. Roe, P. L. (1981). Approximate Riemann solvers, parameter vectors, and difference schemes. *Journal of computational physics*, 43(2), 357-372.
6. Godunov, S., \& Bohachevsky, I. (1959). Finite difference method for numerical computation of discontinuous solutions of the equations of fluid dynamics. *Matematičeskij sbornik*, 47(3), 271-306.
7. Jennings, G. (1974). Discrete shocks. *Communications on pure and applied mathematics*, 27(1), 25-37.
8. Van Leer, B. (1979). Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method. *Journal of computational Physics*, 32(1), 101-136.
9. Yee, H. C., Warming, R. F., \& Harten, A. (1985). Implicit total variation diminishing (TVD) schemes for steady-state calculations. *Journal of Computational Physics*, 57(3), 327-360.
10. Sweby, P. K. (1984). High resolution schemes using flux limiters for hyperbolic conservation laws. *SIAM journal on numerical analysis*, 21(5), 995-1011.
11. Fu, D., \& Ma, Y. (1997). A high order accurate difference scheme for complex flow fields. *Journal of Computational physics*, 134(1), 1-15.
12. 马延文, \& 傅德薰. (1992). 计算空气动力学中一个新的激波捕捉法——耗散比拟法. *中国科学(A辑 数学 物理学 天文学 技术科学)*, 35(3), 263-271.
13. 张涵信. (1984). 差分计算中激波上、下游解出现波动的探讨. *空气动力学学报*(01), 14-21.
14. 张涵信. (1988). 无波动,无自由参数的耗散差分格式. *空气动力学学报*(2).
15. ZHUANG, F., \& ZHANG, H. (1987). Computational fluid dynamics in China. *In 8th Computational Fluid Dynamics Conference* (p. 1134).
16. Harten, A., Engquist, B., Osher, S., \& Chakravarthy, S. R. (1987). Uniformly high order accurate essentially non-oscillatory schemes, III. In Upwind and high-resolution schemes (pp. 218-290). *Springer*, Berlin, Heidelberg.
17. Shu, C. W., \& Osher, S. (1988). Efficient implementation of essentially non-oscillatory shock-capturing schemes. *Journal of computational physics*, 77(2), 439-471.
18. Chakravarthy, S. R. (1990). Some Aspects of Essentially Nonoscillatory (ENO) Formulations for the Euler Equations. *National Aeronautics and Space Administration, Office of Management, Scientific and Technical Information Division.*
19. Jiang, G. S., \& Shu, C. W. (1996). Efficient implementation of weighted ENO schemes. *Journal of computational physics*, 126(1), 202-228.
20. 刘儒勋, \& 舒其望. (2003). 计算流体力学的若干新方法. *科学出版社*.
21. 吴望一,蔡庆东.(2000).时间空间均为二阶的新型NND差分格式. *应用数学和力学* (06),561-572.
22. Liu, X. D., Osher, S., \& Chan, T. (1994). Weighted essentially non-oscillatory schemes. *Journal of computational physics*, 115(1), 200-212.

------

<font size = 2.5>^_^ This is the END of the page / document. Thank you for reading! </font>

<font size = 2.5>If you think this project is helpful to you, do not hesitate to light up the 'Star' or 'Fork'!</font>

<font size = 2.5>Developed or Finished by *pkufzh (Small Shrimp)* on 2022/01/20.</font>

**<u>Contact me:</u>**

Github Page: https://github.com/pkufzh

ResearchGate: https://www.researchgate.net/profile/Zhenghao-Feng

Bilibili Space: https://space.bilibili.com/167343763  




<center><font size = 2.5>This project is protected by the MIT license. Please obey the open source rules.</font></center>