function [T, T_S, q_G] = InitializeThermalModel(B, S)

global T_KELVIN BETA k_0 k_1 ZETA n_ZETA d_ZETA YRSEC
global IC_JC

ZETA     = linspace(0, 1, n_ZETA)';    % NOTE NEW definition (Sept 2017) has ZETA = (z_B)/H
d_ZETA   = 1/(n_ZETA-1);

N         = IC_JC.N;
ic_jc     = IC_JC.ic_jc;

z_HI      = max(S);
z_LO      = min(S);
z_ref     = 0.5*(z_HI+z_LO);
T_ref     = -3.0 + T_KELVIN;
dT_atm_dz = -0.006;
q_G       = 0.07*ones(N, 1)*YRSEC;              % Note units are J/(m^2 yr)
T_S       = T_ref + dT_atm_dz*(S-z_ref);

ic_jc_I   = ic_jc(S>B);

T_S_I     = T_S(ic_jc_I);
T_S_I_mean = mean(T_S(ic_jc_I));
T.A        = repmat(T_S, 1, n_ZETA)';                             % default

if numel(ic_jc_I>0)
  H_I      = S(ic_jc_I)-B(ic_jc_I);
  q_G_I    = q_G(ic_jc_I);
  k_mean   = k_0*exp(-k_1*T_S_I)*YRSEC;    % Not units are J/(m K yr)
  T_B_I    = T_S_I + q_G_I.*H_I./k_mean;
  T_MP_B_I = T_KELVIN-BETA*H_I;
  T_B_min  = min(T_B_I, T_MP_B_I);
  T.A(:,ic_jc_I) = transpose(repmat(T_B_min, 1, n_ZETA)+kron(ZETA', T_S_I-T_B_min));
end

T.grid = 'ZETA';
