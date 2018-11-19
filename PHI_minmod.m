function PHI_out=PHI_minmod(r)

% Flux limiting function for minmod limiter

PHI_out = max(0, min(1,r));