function [s, Sig_d, Sig_alpha] = simply_vce(sig_d_inv, def, G, L, num, M,ls,us,Tzeros)
    % SIMPLY_VCE performs simplified variance component estimation (VCE) to determine
    % the relative weights of two variance components, adjusting them iteratively
    % until convergence. This method is particularly useful in geodetic adjustments.
    %
    % Inputs:
    %   sig_d_inv - Inverse of the data covariance matrix
    %   def       - Deformation vector
    %   Tzeros    - Zero vector for regularization
    %   G         - Design matrix for the problem
    %   L         - Regularization matrix
    %   ls        - Lower bounds for slip values
    %   us        - Upper bounds for slip values
    %   n         - Number of data points
    %   M         - Number of regularization constraints
    %
    % Outputs:
    %   s         - Estimated slip vector
    %   Sig_d     - Updated variance component for data
    %   Sig_alpha  - Updated variance component for regularization
    %  regu_para
    % Initialize -regularization parameters
    Sig_d = 1; 
    Sig_alpha = 1;
    max_iterations = 10; % Maximum iterations to prevent infinite loops
    tolerance = 1e-4; % Convergence tolerance for variance components
    iteration = 1;

    while true
        % Construct the u and Gall matrices for the least squares problem
        Tinv = L / sqrt(Sig_alpha);
        u = [sqrt(sig_d_inv) * def / sqrt(Sig_d); Tzeros];
        Gall = [sqrt(sig_d_inv) * G / sqrt(Sig_d); Tinv];

        % Handle NaN and Inf values
        Gall(isnan(Gall)) = 0; Gall(isinf(Gall)) = 0;
        u(isnan(u)) = 0; u(isinf(u)) = 0;

        % Set optimization options
        options = optimoptions('lsqlin', 'Display', 'none');

        % Solve the least squares problem using lsqlin
        s= lsqlin(Gall, u, [], [], [], [], ls, us, [], options);

        % Calculate residuals for data fitting and regularization terms
        v1 = G * s - def; % Residual for data fitting
        v2 = L * s;       % Residual for regularization
        N1 = G' * sig_d_inv * G / Sig_d;
        N2 = L' * L / Sig_alpha;
        N = N1 + N2;

        % Compute variance components based on current residuals
        c = [
            (v1' * sig_d_inv * v1 / (num - trace(inv(N) * N1))) / Sig_d, 
            (v2' * v2 / (M - trace(inv(N) * N2))) / Sig_alpha
        ];

        % Check for convergence: stop if components are close enough
        if abs(c(1) - c(2)) < tolerance || iteration >= max_iterations || Sig_alpha < 0
            break;
        end

        % Update variance components based on computed values
        Sig_d = c(1) * Sig_d;
        Sig_alpha = c(2) * Sig_alpha;
        iteration = iteration + 1; 
    end

    % Output slip ratio
   regu_para= Sig_alpha /Sig_d 
end
