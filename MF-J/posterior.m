function [accept] = posterior(log_posterior,log_posterior_cand)
%
%   accept = 1: accept the candidate
%   accept = 0: reject the candidate
%

%
%   Compute ratio of posterior probability
%
log_r = log_posterior_cand - log_posterior;

%
%   Accept or reject the candidate
%
if log_r >= 0
    accept = 1;
else
    u1 = rand;
    if log(u1) < log_r
        accept = 1;
    else
        accept = 0;
    end
end
