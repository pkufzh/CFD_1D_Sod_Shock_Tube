# Numerical Simulation of 1-D Sod Shock Tube

**Realize the numerical simulation of 1-D Sod Shock Tube (MATLAB Codes)**

## Information

***Author:*** pkufzh (Small Shrimp)

***Course Name:*** Fundamentals of Computational Fluid Dynamics (CFD)

***Submit Date:*** 2021/12/28

## MATLAB Codes Description

### Main Program

- <font color = 'red'>**Program_Sod_Shock_Tube_Main**</font>
  - **Main Program:** Numerical simulation of 1-D compressible flow (Sod Shock Tube)

### Attached Function Modules

***Important Note:*** Please ensure the following files are placed in **the same folder** with the main program!

- <font color = 'blue'>**Flux_Vect_Split_Common.m**</font>
  - Flux Vector Splitting (FVS) with different methods

- <font color = 'blue'>**Flux_Diff_Split_Common.m**</font>
  - Flux Difference Splitting (FDS) with different methods

- <font color = 'blue'>**Diff_Cons_Common.m**</font>
  - The approach to cal. difference for F_x from F_p and F_n with conservation form
  - Note: the upwind schemes are converted into conservative form

- <font color = 'blue'>**Cal_Minmod.m**</font>
  - Cal. minmod(a, b)
  - Definition: *minmod*(a, b) = 0.5 * (sign(a) + siagn(b)) * min(abs(a), abs(b)), where a, b is 1-D array with same length

- <font color = 'blue'>**Plot_Props.m**</font>
  - Plot the properties with axis coordinate

### Exact Riemann Solution: Referred Functions by Gogol (2021)

***Reference:*** *Gogol* (2021). Sod Shock Tube Problem Solver [Click to the Website](https://www.mathworks.com/matlabcentral/fileexchange/46311-sod-shock-tube-problem-solver), From MATLAB Central File Exchange. Retrieved December 28, 2021.

- **analytic_sod.m**
  - Solve Sod's Shock Tube problem
  - [Reference Page](http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html)

- **sod_func.m**
  - Define functions to be used in analytic_sod.m
  - Initial conditions

- **sod_demo.m**
  - A demo script file to show the use of analytic_sod.m

***Note:***

- The above programs were tested on *MATLAB R2021b*, Windows 64-bit system. 
- For the vector images saved in the assignment PDF are large, loading may be slow.  Thanks for your patient waiting!!!

------

^_^ This is the END of the page / document. Thank you for reading!

*Finished by **pkufzh (Small Shrimp)** on* **2022/01/20**.

**<u>Contact me:</u>**

Github Homepage: https://github.com/pkufzh

Research Gate: https://www.researchgate.net/profile/Zhenghao-Feng

Bilibili Space: https://space.bilibili.com/167343763

-- Who am I? A happy shrimp from Peking University!