function [x, log_posterior, x_acc_burn, acceptance_rate_burning] = burning_phase(x_initial, lb, ub, step_x, N_step_burning, fault_lengths, fault_widths, fault_params,  XY, los,D, sig_d_inv, def, Num, M, ls, us, Tzeros,log_posterior_initial,num_eachF_para)
    % Burning Phase: Adjusts the parameters to reach a stable distribution
    %
    % Inputs:
    %   - x_initial: Initial parameters
    %   - lb, ub: Lower and upper bounds for parameters
    %   - step_x: Step size for parameter updates
    %   - N_step_burning: Number of burning steps
    %   - Length1, Width1, f1, Length2, Width2, f2: Fault simulation parameters
    %   - XY, los: Data for fault simulation
    %   - sig_d_inv, def, Num, M, ls, us, Tzeros: Additional data for computing the posterior
    %
    % Outputs:
    %   - x: Final parameters after the burning phase
    %   - log_posterior: Final log-posterior value
    %   - x_acc_burn: Array of accepted parameters during the burning phase
    %   - acceptance_rate_burning: Acceptance rate during the burning phase
    
    x = x_initial; % Initialize parameters
    log_posterior = log_posterior_initial; % Initial log-posterior (implement this)
    x_acc_burn = []; % Initialize storage for accepted parameters

    n_accept_burning = 0; % Counter for accepted samples
    n_reject_burning = 0; % Counter for rejected samples

    for ii = 1:N_step_burning
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
            x = x_cand; % Update parameters
            log_posterior = log_posterior_cand; % Update log-posterior
            temp=[x_cand', Sig_d, Sig_alpha, ii, log_posterior]
            x_acc_burn = [x_acc_burn; temp]; % Save accepted sample
            n_accept_burning = n_accept_burning + 1; % Increment acceptance count
        else
            n_reject_burning = n_reject_burning + 1; % Increment rejection count
        end
        
        % Log posterior for each step
        log_posterior_burning(ii) = log_posterior;
    end

    % Compute acceptance rate
    acceptance_rate_burning = n_accept_burning / N_step_burning;
end
