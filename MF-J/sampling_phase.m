function [x, log_posterior, x_acc_sampling, mm_sampling, Slip, acceptance_rate_sampling] = sampling_phase(x_burned, lb, ub, step_x, N_step_burning,N_step_sampling, fault_lengths, fault_widths, fault_params, XY, los, D,sig_d_inv, def, Num, M, ls, us, Tzeros, sampling_interval,log_posterior_burned,num_eachF_para)
    % Sampling Phase: Perform sampling after the burning phase
    %
    % Inputs:
    %   - x_initial: Initial parameters after burning
    %   - lb, ub: Lower and upper bounds for parameters
    %   - step_x: Step size for parameter updates
    %   - N_step_sampling: Number of sampling steps
    %   - fault_lengths, fault_widths, fault_params: Fault simulation parameters
    %   - XY, los: Data for fault simulation
    %   - sig_d_inv, def, Num, M, ls, us, Tzeros: Additional data for computing the posterior
    %   - sampling_interval: Interval at which to store intermediate samples
    %
    % Outputs:
    %   - x: Final parameters after the sampling phase
    %   - log_posterior: Final log-posterior value
    %   - x_acc_sampling: Array of accepted parameters during the sampling phase
    %   - mm_sampling: Parameters at specified intervals during sampling
    %   - Slip: Slip data at sampling intervals
    %   - acceptance_rate_sampling: Acceptance rate during the sampling phase

    x = x_burned; % Initialize parameters
    log_posterior = log_posterior_burned; % Initial log-posterior
    x_acc_sampling = []; % Initialize storage for accepted parameters
    mm_sampling = []; % Storage for parameters at sampling intervals
    Slip = []; % Storage for slip data

    n_accept = 0; % Counter for accepted samples
    n_reject = 0; % Counter for rejected samples

    for jj = 1:N_step_sampling
        % Generate candidate parameters
        u = 2 * rand(size(x)) - 1;
        x_cand = x + step_x .* u;

        % Check bounds
        if any(x_cand < lb | x_cand > ub)
            accept = 0; % Reject if out of bounds
        else
            % Simulate fault with candidate parameters
            num_faults = length(fault_lengths); % Number of faults                     
            Faults = fault_geometry(num_faults, fault_lengths, fault_widths, x_cand, fault_params,num_eachF_para);
            Fault = vertcat(Faults{:});


            % Calculate L and G for the candidate fault
            L = calculateL(Faults, fault_params);
            G = calculateG(Fault, XY, los);



            % Compute slip and posterior
            [s, Sig_d, Sig_alpha] = simply_vce(sig_d_inv, def, G, L, Num, M, ls, us, Tzeros);
            log_posterior_cand = compute_log_posterior(Num, Sig_d, Sig_alpha, D, L, G, sig_d_inv, def, s);

            % Compute posterior probability
            accept = posterior(log_posterior, log_posterior_cand);
        end

        % Update parameters if accepted
        if accept == 1
            x = x_cand
            log_posterior = log_posterior_cand;
            x_acc_sampling = [x_acc_sampling; x_cand', Sig_d, Sig_alpha, jj+N_step_burning, log_posterior]; % Save accepted sample
            n_accept = n_accept + 1; % Increment acceptance count
        else
            n_reject = n_reject + 1; % Increment rejection count
        end

        % Store results at specified sampling intervals
        if mod(jj, sampling_interval) == 0
            mm_sampling = [mm_sampling; x_cand', Sig_d, Sig_alpha];
            Slip = [Slip; s'];
        end
    end

    % Compute acceptance rate
    acceptance_rate_sampling = n_accept / N_step_sampling;
end
