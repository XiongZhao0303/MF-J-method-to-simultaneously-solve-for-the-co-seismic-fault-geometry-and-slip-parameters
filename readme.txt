**Package Overview**

This package utilizes the **MF-J method** to simultaneously solve for the co-seismic fault geometry and slip parameters. The fault geometric parameters are obtained through **Markov Chain Monte Carlo (MCMC)** sampling, while the **hyperparameters** and **slip parameters** are determined using the **Variance Component Estimation (VCE)** method.

### To Use the Package

Modify the **main_ file** to input the following parameters:

- **Observation data information**
- **Initial values for the parameters to be estimated**
- **Upper and lower limits for the parameters**
- **Sampling step size**
- **Number of iterations**
- **Approximate maximum slip amount and slip angle information**

---

### Subroutine Functions

1. **burning_phase.m**  
   A MATLAB script for the burn-in phase, designed to stabilize model parameters through a series of MCMC iterations. The goal of the burn-in phase is to eliminate bias from the initial state by performing sufficient sampling, allowing the chain's state to converge to the target posterior distribution.

2. **calculateG.m**  
   Subroutine for calculating the Green's matrix.

3. **calculateL.m**  
   Subroutine for calculating the second-order Laplace smoothing matrix.

4. **compute_log_posterior.m**  
   Subroutine for calculating the posterior probability density of the parameters.

5. **fault_geometry.m**  
   Constructs fault geometry based on the number of faults. It calculates geometric properties (positions, strike, dip, etc.) based on input data like fault lengths, widths, and subdivisions.

6. **logdet.m**  
   Computes the logarithm of the determinant.

7. **okada_InSAR.m**  
   Computes the OKADA model for SAR data.

8. **plot_fault_slip.m**  
   Visualizes the distribution of fault slip.

9. **plot_pdf.m**  
   Plots the posterior probability density distributions of parameters during the sampling period.

10. **plot_scatter_data.m**  
    Used for plotting data and visualizing the fitting or misfit between observed and modeled values.

11. **posterior.m**  
    Determines whether to accept the current sampled parameters.

12. **rigorous_vce.m**  
    Provides a rigorous formulation of the **Variance Component Estimation (VCE)** method for solving hyperparameters.

13. **simply_vce.m**  
    A simplified (recommended) approach for the **Variance Component Estimation (VCE)** method for solving hyperparameters.

14. **sampling_phase.m**  
    In this phase, model parameters are considered to have stabilized. The goal is to obtain the optimal model parameters and uncertainty estimates.

15. **simulatefault.m**  
    Constructs a rectangular fault grid.

16. **slip_bound_constraints.m**  
    Applies boundary constraints on fault slip, including constraints on slip direction and magnitude.

17. **slip_results.m**  
    Displays the results of the fault slip distribution and calculates key quantities such as **seismic moment (Mo)** and **magnitude (Mw)**.

---

### Example Earthquakes Provided for Reference

The program includes four actual earthquake cases for reference:

1. **2017 Iran M7.3 earthquake**
2. **2022 Kashi earthquake (Mw 6.0)**
3. **2024 Wushi earthquake, Xinjiang (Mw 7.0)**
4. **2022 Mw 6.0 doublet-fault earthquake in Hormozgan, Iran**

Among these, the first three are **single-fault** events, while the last one is a **double-fault** earthquake.

### To Run the Program

Simply execute the following scripts:

- `main_iran_2017.m`
- `main_JS.m`
- `main_WS.m`
- `main_iran_2022_double_faults.m`

---

### Contact Information

For any issues regarding the program, please contact the first author, **Zhao Xiong**, at **xzhao2019@whu.edu.cn** (English/Chinese).

---

### Citation

Zhao, X., Zhou, L., Xu, C., et al. (2024). Modified Bayesian method for simultaneously imaging fault geometry and slip distribution with reduced uncertainty, applied to the 2017 Mw 7.3 Sarpol-e Zahab (Iran) earthquake. *Journal of Geodesy*, 98, 106. [https://link.springer.com/article/10.1007/s00190-024-01906-6)

---
