function log_posterior = compute_log_posterior(Num,Sig_d, Sig_alpha, D, L, G, sig_d_inv, def, s)
    % Step 1: Compute p1 and alpha
    % Step 2: Compute each term in the log-posterior function(Fukuda and
    % johnson,2010)
    LL1 = (-Num / 2) * log(Sig_d);
    LL2 = (-1/2) * logdet(D);
    LL3 = (1/2) * logdet(L' * L / Sig_alpha);
    LL4 = (-1/2) * logdet((G' * sig_d_inv * G / Sig_d) + (L' * L / Sig_alpha));
    % Step 3: Compute the negative log-likelihood part (f)
    f = (-1/2) * (((def - G * s)' * sig_d_inv * (def - G * s)) / Sig_d + s' * L' * L * s / Sig_alpha);
    % Step 4: Sum up all terms to get the log-posterior
    log_posterior = LL1 + LL2 + LL3 + LL4 + f;
end
