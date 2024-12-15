function [Slip, Sig_d, Sig_alpha] = rigorous_vce(sig_d_inv, def, G, L, num, M,ls,us,Tzeros)
    % Initialize variables
    Sig_d = 1; 
    Sig_alpha = 1;
    max_iterations = 10; % Maximum iterations to prevent infinite loops
    tolerance = 1e-4; % Convergence tolerance for variance components
    iterations = 1;
    while true  
    Tinv = L / sqrt(Sig_alpha);
    u = [sqrt(sig_d_inv) * def / sqrt(Sig_d);  Tzeros];
    % Compute Gall matrix
    Gall = [sqrt(sig_d_inv) * G / sqrt(Sig_d); Tinv];
    % Handle any NaN or Inf values in Gall and u
    Gall(isnan(Gall)) = 0;
    Gall(isinf(Gall)) = 0;
    u(isnan(u)) = 0;
    u(isinf(u)) = 0;
    % Set optimization options
    options = optimoptions('lsqlin', 'Display', 'none');
    % Solve the least squares problem using lsqlin
    Slip = lsqlin(Gall, u, [], [], [], [], ls, us, [], options);
%     Slip=inv(G'*sig_d_inv*G/Sig_d+L'*L/Sig_alpha)*G'*sig_d_inv*def/Sig_d;
        % Calculate v1 and v2
        v1 = G * Slip - def;
        v2 = L * Slip;
        % Calculate w
        w = [v1'  * sig_d_inv * v1/Sig_d; v2' * v2/Sig_alpha];
        % Compute N1 and N2
        N1 = G'* sig_d_inv * G /Sig_d ;
        N2 = L' * L/ Sig_alpha;
        N = N1 + N2;

        % Compute s1 to s4
        s1 = 2 * trace(inv(N) * N1) - trace(inv(N) * N1 * inv(N) * N1);
        s2 = trace(inv(N) * N1 * inv(N) * N2);
        s3 = trace(inv(N) * N2 * inv(N) * N1);
        s4 = 2 * trace(inv(N) * N2) - trace(inv(N) * N2 * inv(N) * N2);
        s =[ num - s1, s2; s3, M - s4];
        % Compute c
        c = inv(s' * s) * s' * w;
        % Convergence check or other break conditions
        if (abs(c(1) - c(2)) < tolerance) || Sig_alpha< 0 || iterations >= max_iterations
            break;
        else
            % Update p1 and p2 based on c
            c1 = c(1);
            Sig_d =  c(1) * Sig_d;
            Sig_alpha = c(2)  *Sig_alpha;
            iterations = iterations + 1; 
        end
    end
p2=Sig_d/Sig_alpha
end