% =========================================================================
% =================== Bayesian MCMC with VCE Framework ====================
% ========================== Author: Zhao Xiong ============================
% =========================================================================
% Cite this code:
% Zhao, X., Zhou, L., Xu, C., et al. 
% Modified Bayesian method for simultaneously imaging fault geometry and slip distribution 
% with reduced uncertainty, applied to the 2017 Mw 7.3 Sarpol-e Zahab (Iran) earthquake. 
% Journal of Geodesy, 98, 106 (2024). https://doi.org/10.1007/s00190-024-01906-6
% ========================================================================

% ========================================================================
% 
% Fault Parameters and Variance Components:
% - EX, EY, EZ:      Coordinates of the upper-left vertex of the fault in the 
%                    east-west, north-south, and vertical directions.
% - Strike, Dip:     Fault orientation parameters.
% - Sig_d, Sig_alpha: Variance components for observed and virtual data.

% Unknown Parameters:
% - Bayesian MCMC method:
%   * EX, EY, EZ, Strike, Dip
% - VCE method:
%   * Slip, Sig_d, Sig_alpha

clc; 
close all; 
clear all;

% Output directories
folderName1 = 'out/JS_2020/burn/';
folderName2 = 'out/JS_2020/sampling/';
mkdir(folderName1);
mkdir(folderName2);

% Load custom colormap
mycolormap = load('mycolormap.mat');
mycolormap = mycolormap.ans;

%% ======================== Input Fault Information ==========================
% Set the number of faults, the length and width of each fault, 
% and the number of subdivisions for each fault.
fault_lengths = 40;  % Fault lengths
fault_widths = 20;   % Fault widths
fault_params = {[21, 11]}; % Subdivisions along strike and dip

% Set initial values and parameter bounds for estimation (EX, EY, EZ, Strike, Dip)
x_initial = [0;   0;   0; 220;10]; % Fault 

% Parameter Bounds
lb = [ -30; -30;  -15;  200; 5];   % Lower bounds
ub = [  30;  30;   0;   360; 30];   % Upper bounds
% Step Sizes for Parameter Sampling
step_x = [0.5; 0.5; 0.5; 1; 1];  % Fault 1

% MCMC Parameters
N_step_burning = 5000;  % Burn-in steps
N_step_sampling = 5000;  % Sampling steps
sampling_interval = 25;  % Sampling interval

% Fault Slip Constraints(This information can be obtained from
% organizations such as GCMT, USGS, GFZ...)
rake = 100;     % Slip angle
MAXSLIP = 1; % Maximum slip

%% ===================== Load Observation Data ============================
% Observation Data File Path
data_file_path = './inp/JS_2020/insar_4tree_llh.inp';

% Load the observation data from the specified file
[lon, lat, X, Y, def, lose, losn, losu, var] = textread(data_file_path, '%f %f %f %f %f %f %f %f %f');

% Explanation of variables:
% - lon, lat: Longitude and latitude of the observation data.
% - X, Y (km): UTM (Universal Transverse Mercator) coordinates converted from 
%   longitude and latitude. To ensure the fault remains close to the deformation zone, 
%   the epicenter coordinates are typically used as the reference point for the UTM conversion.
% - def (m): SAR line-of-sight (LOS) deformation observed in the data.
% - lose, losn, losu: Conversion coefficients for LOS deformation in the 
%   east, north, and vertical directions, respectively (ASC/DES: about -0.6/0.6, -0.11, 0.75).

D = diag(var);                 % Variance matrix
sig_d_inv = inv(D);            % Weight matrix
sig_d_inv_sqrt = sqrt(sig_d_inv); % Square root of the weight matrix

% Prepare additional data
Num = size(def, 1);            % Number of observations
los = [lose, losn, losu];      % Line-of-sight (LOS) vector
XY = [X, Y];                   % Coordinate data

%% ===================== Model Initialization ==============================
num_faults = length(fault_lengths); % Number of faults
% Specify the number of parameters to be estimated for each fault
num_eachF_para = length(x_initial)/num_faults; 
% NOTE:The variable `num_eachF_para` represents the number of parameters to be estimated for each fault.
% By default, this number is five: EX, EY, EZ, STRIKE, and DIP. 
% If you need to change the number of parameters, ensure you carefully adjust 
% the corresponding positions of each parameter when calling the `simulatefault.m` function.
[Faults, subdivisions, M_total] = fault_geometry(num_faults, fault_lengths, fault_widths, x_initial, fault_params,num_eachF_para);


% Combine all faults into one array
Fault = vertcat(Faults{:});

% Initialize Slip Values
Tzeros = zeros(M_total, 1);

% Compute Constraints
L = calculateL(Faults, fault_params);
G = calculateG(Fault, XY, los);
[ls, us] = slip_bound_constraints(fault_params, M_total, rake, MAXSLIP);
[s, Sig_d, Sig_alpha] = simply_vce(sig_d_inv, def, G, L, Num, M_total, ls, us, Tzeros);

% Calculate Initial Log-Posterior
log_posterior_initial = compute_log_posterior(Num, Sig_d, Sig_alpha, D, L, G, sig_d_inv, def, s);

% Plot initial Fault Slip Results
plot_fault_slip(Fault, s, subdivisions,mycolormap);

%% ===================== Burn-In Phase ====================================
n_accept_burning = 0; % Count accepted samples
n_reject_burning = 0; % Count rejected samples
log_posterior_burning = zeros(N_step_burning, 1); % Log-posterior during burn-in

[x_burned, log_posterior_burned, x_acc_burn, acceptance_rate_burning] = burning_phase(...
    x_initial, lb, ub, step_x, N_step_burning, fault_lengths, fault_widths, fault_params, ...
    XY, los, D, sig_d_inv, def, Num, M_total, ls, us, Tzeros, log_posterior_initial,num_eachF_para);

fileName1 = fullfile(folderName1, 'x_acc_burn.txt');
save(fileName1, 'x_acc_burn', '-ascii'); 
%% ===================== Sampling Phase ===================================
[x_sampled, log_posterior_sampled, x_acc_sampling, mm_sampling, Slip, acceptance_rate_sampling] = sampling_phase(...
    x_burned, lb, ub, step_x, N_step_burning, N_step_sampling, fault_lengths, fault_widths, fault_params, ...
    XY, los, D, sig_d_inv, def, Num, M_total, ls, us, Tzeros, sampling_interval, log_posterior_burned,num_eachF_para);

% Save Sampling Results
fileName2 = fullfile(folderName2, 'x_acc_sampling.txt');
save(fileName2, 'x_acc_sampling', '-ascii'); 

fileName3 = fullfile(folderName2, 'mm_sampling.txt');
save(fileName3, 'mm_sampling', '-ascii'); 

%% Burn + Sampling Acceptance
para_acc = [x_acc_burn; x_acc_sampling];
plot_convergence_para(para_acc, num_faults,num_eachF_para);
plot_pdf(mm_sampling, num_faults,num_eachF_para);

%% Final Results and Visualization

% Compute mean and standard deviation of sampled parameters
X_o = mean(mm_sampling, 1);          % Mean of model parameters
X_std = std(mm_sampling, 0, 1);     % Standard deviation of model parameters

% Generate fault geometry based on estimated parameters
Faults = fault_geometry(num_faults, fault_lengths, fault_widths, X_o, fault_params, num_eachF_para);
Fault = vertcat(Faults{:});         % Combine fault geometries into a single matrix

% Calculate the Green's function matrix
G = calculateG(Fault, XY, los);

% Compute mean and standard deviation of slip
s_o = mean(Slip, 1);                % Mean slip values
s_std = std(Slip, 0, 1);            % Standard deviation of slip values

% Compute the modeled deformation
d_o = G * s_o';

% Calculate seismic moment and magnitude for each fault
slip_results(Faults, num_faults, subdivisions, s_o);

% Compute Root Mean Square (RMS) error
RMS = sqrt((def - d_o)' * (def - d_o) / length(def));

% Calculate Variance Reduction (VR)
VR = (1 - (def - d_o)' * (def - d_o) / (d_o' * d_o)) * 100;

% Visualization
% Plot fault slip (mean values)
plot_fault_slip(Fault, s_o, subdivisions, mycolormap);

% Plot fault slip (standard deviation)
plot_fault_slip(Fault, s_std, subdivisions, mycolormap);

% Plot scatter data to visualize misfit
plot_scatter_data(los, def, d_o, lon, lat); % Misfit visualization


