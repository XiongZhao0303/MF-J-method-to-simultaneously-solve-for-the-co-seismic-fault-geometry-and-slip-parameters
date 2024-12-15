This package utilizes the MF-J method to simultaneously solve for the co-seismic fault geometry and slip parameters. The fault geometric parameters are obtained through MCMC sampling, while the hyperparameters and slip parameters are determined using the VCE method.


Here, you only need to modify the main_MF-J file to input the observation data information, initial values for the parameters to be estimated, the upper and lower limits, the sampling step size, and the number of iterations.

Below is an introduction to the functions of each subroutine:


calculateG.m is the subroutine for calculating the Green's matrix.

calculateL.m is the subroutine for calculating the second-order Laplace smoothing matrix.


compute_log_posterior.m is the subroutine for calculating the posterior probability density of the parameters.

logdet.m computes the logarithm of the determinant.

okada_InSAR.m is used to compute the OKADA model for SAR data.

plot_fault_slip.m is used to visualize the distribution of fault slip.

plot_pdf.m is used to plot the posterior probability density distributions of parameters during the sampling period.


posterior.m is used to determine whether to accept the current sampled parameters.


rigorous_vce.m and simply_vce.m provide two equivalent forms of the Variance Component Estimation (VCE) method for solving hyperparameters: one using a rigorous formulation and the other a simplified approach.


simulatefault.m constructs a rectangular fault grid.


slip_bound_constraints.m applies boundary constraints on fault slip, as well as constraints on slip direction and magnitude.

For any issues regarding the program, please contact the first author, Zhao Xiong, at xzhao2019@whu.edu.cn (English/Chinese).
