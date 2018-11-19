function par=SetParameters(MODEL)

global      n_GLEN A_GLEN RHO RHO_w g OMEGA ETA_0
global      nx_sub_default ny_sub_default nx_skirt ny_skirt nx_pad ny_pad N_xy
global      CFL A_0 Q_creep R_gas
global      d_XI XI XI_term
global      BETA k_0 k_1 c_0 c_1 YRSEC T_KELVIN
global      N_xy wwt_kc wwt_kp wwt_km wwt_kpp wwt_kmm 
global      KAPPA_ref_yr H_min dH_min LAMBDA_T_max dt_max_T

% Function to set the basic physical parameters of the model and make them
% visible to other parts of the program by making them global

% RGM properties and parameters

n_GLEN   = 3;                         % Glen's flow law exponent
A_GLEN   = 7.5738e-17;                % Cuffey & Paterson (4th ed) Glen's law parameter in Pa^{-n} yr^{-1} units (same as A_GLEN=2.4e-24 Pa^{-3} s^{-1})
m_SLIDE  = 2;                         % Sliding law exponent

RHO      = 900;                       % Ice density (kg/m^3)
RHO_w    = 1000;                      % Water density (kg/m^3)
g        = 9.80;                      % Gravity (m/s^2)

Q_creep  = 115e3;
R_gas    = 8.314462;
T_KELVIN = 273.15;
ETA_0    = 1.0e15;                    % Residual (low-stress limit) ice viscosity (Pa s)   

OMEGA    = 1.5;                       % Hindmarsh parameter
CFL      = 0.125;
nx_skirt = 4;
ny_skirt = 4;
nx_pad   = 2;
ny_pad   = 2;

par.n_GLEN    = n_GLEN;                % Glen's flow law exponent
par.A_GLEN    = A_GLEN;                % Cuffey & Paterson (4th ed) Glen's law parameter in Pa^{-n} yr^{-1} units (same as A_GLEN=2.4e-24 Pa^{-3} s^{-1})

par.RHO       = RHO;                   % Ice density (kg/m^3)
par.RHO_w     = RHO_w;                 % Water density (kg/m^3)
par.g         = g;                     % Gravity (m/s^2)
par.ETA_0     = ETA_0;                 % Residual ice viscosity (Pa s)
par.OMEGA     = OMEGA;                 % Hindmarsh parameter
par.CFL       = CFL;
par.nx_skirt  = nx_skirt;
par.ny_skirt  = ny_skirt;
par.nx_pad    = nx_pad;
par.ny_pad    = ny_pad;

A_0           = A_GLEN/exp(-Q_creep/(R_gas*T_KELVIN));   % Sets flow to isothermal law 

par.A_0      = A_0;
par.Q_creep  = Q_creep;
par.R_gas    = R_gas;

YRSEC        = 3600*24*365.25;
par.YRSEC    = YRSEC;

% Thermal parameters

k_0      = 9.828;    % W/(m K)
k_1      = 0.0057;
c_0      = 152.5;
c_1      = 7.122;

n_ZETA   = MODEL.n_ZETA;
n_XI     = MODEL.n_XI;

H_min    = 20;   % For thinner ice the Crank-Nicholson solver is not used (Try H_min = 200 for Cordilleran Ice Sheet)

KAPPA_ref_yr = YRSEC*1.0e-6;

LAMBDA_T_max = 0.25;   %% Increasing this while maintaining stability & accuracy will speed things up

dH_min       = H_min/(n_ZETA-1);
dt_max_T     = LAMBDA_T_max*dH_min^2/KAPPA_ref_yr;

BETA     = 8.7e-4;     % K/m
T_KELVIN = 273.15;

par.T_KELVIN  = T_KELVIN;
par.BETA      = BETA;

XI        = MODEL.XI';                                 % Follow Greve & Blatter not Marshall !!
d_XI      = MODEL.d_XI;

XI_term   = repmat((1-XI).^n_GLEN, 1, N_xy);

wt_kc   = [-3  zeros(1, n_XI-2)  3]/(2*d_XI);
wt_kp   = [+4  +ones(1, n_XI-2)  0]/(2*d_XI);
wt_km   = [0   -ones(1, n_XI-2) -4]/(2*d_XI);

wt_kpp  = [-1  zeros(1, n_XI-1)  ]/(2*d_XI);
wt_kmm  = [zeros(1, n_XI-1)    +1]/(2*d_XI);

wwt_kpp = repmat(wt_kpp, N_xy, 1)';
wwt_kmm = repmat(wt_kmm, N_xy, 1)';

wwt_kp  = repmat(wt_kp, N_xy, 1)';
wwt_kc  = repmat(wt_kc, N_xy, 1)';
wwt_km  = repmat(wt_km, N_xy, 1)';

par.nx_sub_default = 2*nx_skirt+1;
par.ny_sub_default = 2*ny_skirt+1;

nx_sub_default = par.nx_sub_default;
ny_sub_default = par.ny_sub_default;

par.k_0  = k_0;    % W/(m K)
par.k_1  = k_1;
par.c_0  = c_0;
par.c_1  = c_1;

H_min        = 20;   % For thinner ice the Crank-Nicholson solver is not used (just assume linear ramp)

par.H_min    = H_min;
