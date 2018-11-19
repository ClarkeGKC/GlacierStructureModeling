function [D_out, v_out] = D_vel_calc_new(H, grad_S, T, dS_dl, THERMAL)

% This revised approach puts all the suspect elements in the same cage so
% they can be closely watched. It also simplifies the Diffusion code (a lot)
% Note that Q_out is not used but output for possible analysis

global  A_0 A_GLEN n_GLEN R_gas Q_creep RHO g 
global  XI nx ny n_XI N_xy dx XI_term

% Note that h is dlstance above the bed B but h is ordered top-down
% as is XI

if THERMAL
  A = A_0*exp(-Q_creep./(R_gas*T));
else
  A = A_GLEN.*ones(n_XI, N_xy);
end 

% Calculated once only in SetParameters.m
% XI_term = repmat((1-XI).^n_GLEN, 1, N_xy);

F_1  = A.*XI_term;

I_1       = cumtrapz(XI, F_1);
I_2    = trapz(XI, I_1)';

H_term = H.^(n_GLEN+1).*grad_S.^(n_GLEN-1).*dS_dl;

v_out  = -2*(RHO*g)^n_GLEN*repmat(H_term', n_XI, 1).*I_1;
D      = 2*(RHO*g).^n_GLEN.*H.^(n_GLEN+2).*grad_S.^(n_GLEN-1).*I_2;

D_out  = D/dx^2;

% Q_out is calculated here and could be output to calling program but not used there
% It is only to allow comparison with Q calculated in Diffusion.m
% The following is correct but no longer computed nor output from this function

% Q_out = -2*(RHO*g)^n_GLEN*(H/A_XI).^(n_GLEN+2).*grad_S.^(n_GLEN-1).*J_2.*dS_dl;


