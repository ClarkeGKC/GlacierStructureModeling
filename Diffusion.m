function [D_IC_jc, D_IP_jc, D_ic_JC, D_ic_JP, V, L, PHI_ic_jc, E, H_struct] = Diffusion(S, B, C_t, T, b_dot, t)

% Limited MUSCL scheme diffusion calculation

% This version of diffusion is provisional and reinstates the four MUSCL
% cases in order to establish whether the treatment of the u_IC_jc_m and
% u_IC_jc_p fields, etc. leads to symmetry breaking

% Previous versions of diffusion.m did not treat the non-isothermal
% velocity calculations properly (although the corresponding flux
% calculations were correct)

% This version aims to be fully consistent in how H values are used to
% construct the velocity field components. The aim is to use same logic for
% finding H_IC_jc ... etc as for D_IC_jc ... etc.
% This change has now been tested and it has been confirmed that
% isothermal and non-isothermal code delivers the same results (for
% isothermal ice)

global  A_const n_GLEN RHO g dx m_SLIDE dt
global  A_0 R_gas Q_creep YRSEC
global  IC_JC
global  kc kp km kpp kmm wwt_kc wwt_kp wwt_km wwt_kpp wwt_kmm
global  THERMAL n_XI XI d_XI
global  B_IC_jc B_IP_jc B_ic_JC B_ic_JP  dB_dx_IC_jc dB_dx_IP_jc dB_dy_ic_JC dB_dy_ic_JP

ic_jc  = IC_JC.ic_jc;
ip_jc  = IC_JC.ip_jc;
ip_jp  = IC_JC.ip_jp;
ipp_jc = IC_JC.ipp_jc;
im_jc  = IC_JC.im_jc;
imm_jc = IC_JC.imm_jc;
ic_jp  = IC_JC.ic_jp;
ic_jpp = IC_JC.ic_jpp;
ic_jm  = IC_JC.ic_jm;
ic_jmm = IC_JC.ic_jmm;
im_jm  = IC_JC.im_jm;
ip_jm  = IC_JC.ip_jm;
im_jp  = IC_JC.im_jp;
dx     = IC_JC.dx;
dy     = IC_JC.dy;
N_xy   = IC_JC.N;
nx     = IC_JC.nx;  % only for debugging
ny     = IC_JC.ny;

dS_dx_IC_jc = (S(ic_jc)-S(im_jc))/dx;
dS_dx_IP_jc = (S(ip_jc)-S(ic_jc))/dx;

dS_dy_ic_JC = (S(ic_jc)-S(ic_jm))/dy;
dS_dy_ic_JP = (S(ic_jp)-S(ic_jc))/dy;

A_tilde = A_const/dx^2;
C_tilde = C_t*(RHO*g)^m_SLIDE/dx^2;    % Here is where surging is controlled. Note that C_tilde is now an array not a constant

% Defining C_slide on a staggered grid is essential for proper handling of surge activation zone

C_IC_jc = 0.5*(C_tilde(im_jc)+C_tilde(ic_jc));
C_IP_jc = 0.5*(C_tilde(ic_jc)+C_tilde(ip_jc));
C_ic_JC = 0.5*(C_tilde(ic_jm)+C_tilde(ic_jc));
C_ic_JP = 0.5*(C_tilde(ic_jc)+C_tilde(ic_jp));

H      = S-B;

rx_ic_jc  = (H(ic_jc)-H(im_jc))./(H(ip_jc)-H(ic_jc));
rx_ip_jc  = (H(ip_jc)-H(ic_jc))./(H(ipp_jc)-H(ip_jc));
rx_im_jc  = (H(im_jc)-H(imm_jc))./(H(ic_jc)-H(im_jc));

H_IP_jc_m = H(ic_jc) + 0.5*PHI_minmod(rx_ic_jc).*(H(ip_jc)-H(ic_jc));
H_IP_jc_p = H(ip_jc) - 0.5*PHI_minmod(rx_ip_jc).*(H(ipp_jc)-H(ip_jc));
H_IC_jc_m = H(im_jc) + 0.5*PHI_minmod(rx_im_jc).*(H(ic_jc)-H(im_jc));
H_IC_jc_p = H(ic_jc) - 0.5*PHI_minmod(rx_ic_jc).*(H(ip_jc)-H(ic_jc));

ry_ic_jc  = (H(ic_jc)-H(ic_jm))./(H(ic_jp)-H(ic_jc));
ry_ic_jp  = (H(ic_jp)-H(ic_jc))./(H(ic_jpp)-H(ic_jp));
ry_ic_jm  = (H(ic_jm)-H(ic_jmm))./(H(ic_jc)-H(ic_jm));

H_ic_JP_m = H(ic_jc) + 0.5*PHI_minmod(ry_ic_jc).*(H(ic_jp)-H(ic_jc));
H_ic_JP_p = H(ic_jp) - 0.5*PHI_minmod(ry_ic_jp).*(H(ic_jpp)-H(ic_jp));
H_ic_JC_m = H(ic_jm) + 0.5*PHI_minmod(ry_ic_jm).*(H(ic_jc)-H(ic_jm));
H_ic_JC_p = H(ic_jc) - 0.5*PHI_minmod(ry_ic_jc).*(H(ic_jp)-H(ic_jc));

% Use the same interpolation scheme for the temperature field as for the H field

% 1. MUSCL scheme UP ===================================================

L_H                       = false(N_xy, 1);
L_H(H_ic_JP_m<=H_ic_JP_p) = true;

L_S                       = false(N_xy, 1);
L_S(S(ic_jp)<=S(ic_jc))   = true;

grad_S2_ic_JP = ((S(ip_jc) - S(im_jc) + S(ip_jp) - S(im_jp))/(4*dx)).^2 ...
  + ((S(ic_jp) - S(ic_jc))/dy).^2;

grad_S_ic_JP  = sqrt(grad_S2_ic_JP);

if THERMAL
  T_ic_JP = 0.5*(T(:,ic_jc)+T(:,ic_jp));
else
  T_ic_JP = [];
end

[D_ic_JP_m_A, v_ic_JP_m_A] = D_vel_calc_new(H_ic_JP_m, grad_S_ic_JP, T_ic_JP, dS_dy_ic_JP, THERMAL);
[D_ic_JP_p_A, v_ic_JP_p_A] = D_vel_calc_new(H_ic_JP_p, grad_S_ic_JP, T_ic_JP, dS_dy_ic_JP, THERMAL);

D_ic_JP_m_C = C_ic_JP.*H_ic_JP_m.^(m_SLIDE+1).*grad_S2_ic_JP.^((m_SLIDE-1)/2);
D_ic_JP_m   = D_ic_JP_m_A + D_ic_JP_m_C;

D_ic_JP_p_C = C_ic_JP.*H_ic_JP_p.^(m_SLIDE+1).*grad_S2_ic_JP.^((m_SLIDE-1)/2);
D_ic_JP_p   = D_ic_JP_p_A + D_ic_JP_p_C;

L_min                = D_ic_JP_m<=D_ic_JP_p;

D_ic_JP_min          = D_ic_JP_p;
D_ic_JP_min(L_min)   = D_ic_JP_m(L_min);

D_ic_JP_min_A        = D_ic_JP_p_A;
D_ic_JP_min_A(L_min) = D_ic_JP_m_A(L_min);

D_ic_JP_min_C        = D_ic_JP_p_C;
D_ic_JP_min_C(L_min) = D_ic_JP_m_C(L_min);

H_ic_JP_min          = H_ic_JP_p;
H_ic_JP_min(L_min)   = H_ic_JP_m(L_min);

L_max                = D_ic_JP_m>=D_ic_JP_p;

D_ic_JP_max          = D_ic_JP_p;
D_ic_JP_max(L_max)   = D_ic_JP_m(L_max);

D_ic_JP_max_A        = D_ic_JP_p_A;
D_ic_JP_max_A(L_max) = D_ic_JP_m_A(L_max);

D_ic_JP_max_C        = D_ic_JP_p_C;
D_ic_JP_max_C(L_max) = D_ic_JP_m_C(L_max);

H_ic_JP_max          = H_ic_JP_p;
H_ic_JP_max(L_max)   = H_ic_JP_m(L_max);

v_ic_JP_min_A          = v_ic_JP_p_A;
v_ic_JP_min_A(:,L_min) = v_ic_JP_m_A(:,L_min);
v_ic_JP_max_A          = v_ic_JP_p_A;
v_ic_JP_max_A(:,L_max) = v_ic_JP_m_A(:,L_max);

D_ic_JP              = zeros(N_xy, 1);

% Use local slope to identify upstream direction

D_ic_JP( L_S &  L_H) = D_ic_JP_min( L_S &  L_H);
D_ic_JP( L_S & ~L_H) = D_ic_JP_max( L_S & ~L_H);
D_ic_JP(~L_S &  L_H) = D_ic_JP_max(~L_S &  L_H);
D_ic_JP(~L_S & ~L_H) = D_ic_JP_min(~L_S & ~L_H);

D_ic_JP_A              = zeros(N_xy, 1);

D_ic_JP_A( L_S &  L_H) = D_ic_JP_min_A( L_S &  L_H);
D_ic_JP_A( L_S & ~L_H) = D_ic_JP_max_A( L_S & ~L_H);
D_ic_JP_A(~L_S &  L_H) = D_ic_JP_max_A(~L_S &  L_H);
D_ic_JP_A(~L_S & ~L_H) = D_ic_JP_min_A(~L_S & ~L_H);

D_ic_JP_C              = zeros(N_xy, 1);

D_ic_JP_C( L_S &  L_H) = D_ic_JP_min_C( L_S &  L_H);
D_ic_JP_C( L_S & ~L_H) = D_ic_JP_max_C( L_S & ~L_H);
D_ic_JP_C(~L_S &  L_H) = D_ic_JP_max_C(~L_S &  L_H);
D_ic_JP_C(~L_S & ~L_H) = D_ic_JP_min_C(~L_S & ~L_H);

H_ic_JP                = zeros(N_xy, 1);

H_ic_JP( L_S &  L_H)   = H_ic_JP_min( L_S &  L_H);
H_ic_JP( L_S & ~L_H)   = H_ic_JP_max( L_S & ~L_H);
H_ic_JP(~L_S &  L_H)   = H_ic_JP_max(~L_S &  L_H);
H_ic_JP(~L_S & ~L_H)   = H_ic_JP_min(~L_S & ~L_H);

v_ic_JP_A                = zeros(n_XI, N_xy);
v_ic_JP_A(:, L_S &  L_H) = v_ic_JP_min_A(:, L_S &  L_H);
v_ic_JP_A(:, L_S & ~L_H) = v_ic_JP_max_A(:, L_S & ~L_H);
v_ic_JP_A(:,~L_S &  L_H) = v_ic_JP_max_A(:,~L_S &  L_H);
v_ic_JP_A(:,~L_S & ~L_H) = v_ic_JP_min_A(:,~L_S & ~L_H);


% 2. MUSCL scheme DOWN ==================================================

L_H                       = false(N_xy, 1);
L_H(H_ic_JC_m<=H_ic_JC_p) = true;

L_S                     = false(N_xy, 1);
L_S(S(ic_jc)<=S(ic_jm)) = true;

grad_S2_ic_JC = ((S(ip_jm) - S(im_jm) + S(ip_jc) - S(im_jc))/(4*dx)).^2 ...
              + ((S(ic_jc) - S(ic_jm))/dy).^2;
grad_S_ic_JC  = sqrt(grad_S2_ic_JC);            
                      
if THERMAL 
  T_ic_JC = 0.5*(T(:,ic_jm)+T(:,ic_jc));
else
  T_ic_JC = [];
end

[D_ic_JC_m_A, v_ic_JC_m_A] = D_vel_calc_new(H_ic_JC_m, grad_S_ic_JC, T_ic_JC, dS_dy_ic_JC, THERMAL);
[D_ic_JC_p_A, v_ic_JC_p_A] = D_vel_calc_new(H_ic_JC_p, grad_S_ic_JC, T_ic_JC, dS_dy_ic_JC, THERMAL);

D_ic_JC_m_C = C_ic_JC.*H_ic_JC_m.^(m_SLIDE+1).*grad_S2_ic_JC.^((m_SLIDE-1)/2);
D_ic_JC_m   = D_ic_JC_m_A + D_ic_JC_m_C;

D_ic_JC_p_C = C_ic_JC.*H_ic_JC_p.^(m_SLIDE+1).*grad_S2_ic_JC.^((m_SLIDE-1)/2);
D_ic_JC_p   = D_ic_JC_p_A + D_ic_JC_p_C;

L_min              = D_ic_JC_m<=D_ic_JC_p;

D_ic_JC_min        = D_ic_JC_p;
D_ic_JC_min(L_min) = D_ic_JC_m(L_min);

D_ic_JC_min_A        = D_ic_JC_p_A;
D_ic_JC_min_A(L_min) = D_ic_JC_m_A(L_min);

D_ic_JC_min_C        = D_ic_JC_p_C;
D_ic_JC_min_C(L_min) = D_ic_JC_m_C(L_min);

H_ic_JC_min        = H_ic_JC_p;
H_ic_JC_min(L_min) = H_ic_JC_m(L_min);

L_max              = D_ic_JC_m>=D_ic_JC_p;

D_ic_JC_max        = D_ic_JC_p;
D_ic_JC_max(L_max) = D_ic_JC_m(L_max);

D_ic_JC_max_A        = D_ic_JC_p_A;
D_ic_JC_max_A(L_max) = D_ic_JC_m_A(L_max);

D_ic_JC_max_C        = D_ic_JC_p_C;
D_ic_JC_max_C(L_max) = D_ic_JC_m_C(L_max);

H_ic_JC_max        = H_ic_JC_p;
H_ic_JC_max(L_max) = H_ic_JC_m(L_max);

v_ic_JC_min_A          = v_ic_JC_p_A;
v_ic_JC_min_A(:,L_min) = v_ic_JC_m_A(:,L_min);
v_ic_JC_max_A          = v_ic_JC_p_A;
v_ic_JC_max_A(:,L_max) = v_ic_JC_m_A(:,L_max);

D_ic_JC              = zeros(N_xy, 1);

% Use local slope to identify upstream direction

D_ic_JC( L_S &  L_H) = D_ic_JC_min( L_S &  L_H);
D_ic_JC( L_S & ~L_H) = D_ic_JC_max( L_S & ~L_H);
D_ic_JC(~L_S &  L_H) = D_ic_JC_max(~L_S &  L_H);
D_ic_JC(~L_S & ~L_H) = D_ic_JC_min(~L_S & ~L_H);

D_ic_JC_A              = zeros(N_xy, 1);

D_ic_JC_A( L_S &  L_H) = D_ic_JC_min_A( L_S &  L_H);
D_ic_JC_A( L_S & ~L_H) = D_ic_JC_max_A( L_S & ~L_H);
D_ic_JC_A(~L_S &  L_H) = D_ic_JC_max_A(~L_S &  L_H);
D_ic_JC_A(~L_S & ~L_H) = D_ic_JC_min_A(~L_S & ~L_H);

D_ic_JC_C              = zeros(N_xy, 1);

D_ic_JC_C( L_S &  L_H) = D_ic_JC_min_C( L_S &  L_H);
D_ic_JC_C( L_S & ~L_H) = D_ic_JC_max_C( L_S & ~L_H);
D_ic_JC_C(~L_S &  L_H) = D_ic_JC_max_C(~L_S &  L_H);
D_ic_JC_C(~L_S & ~L_H) = D_ic_JC_min_C(~L_S & ~L_H);

% The following is open for discussion. I use H_ic_JC etc for calculating
% ice velocity components from discharge. There is some argument for
% using the limited H_ic_JC values rather than simply computed ones but
% it probably makes little practical difference (certainly not to the
% geometric evolution -- only to the flow and tracer evolution

H_ic_JC              = zeros(N_xy, 1);

H_ic_JC( L_S &  L_H) = H_ic_JC_min( L_S &  L_H);
H_ic_JC( L_S & ~L_H) = H_ic_JC_max( L_S & ~L_H);
H_ic_JC(~L_S &  L_H) = H_ic_JC_max(~L_S &  L_H);
H_ic_JC(~L_S & ~L_H) = H_ic_JC_min(~L_S & ~L_H);
 
v_ic_JC_A                = zeros(n_XI, N_xy);
v_ic_JC_A(:, L_S &  L_H) = v_ic_JC_min_A(:, L_S &  L_H);
v_ic_JC_A(:, L_S & ~L_H) = v_ic_JC_max_A(:, L_S & ~L_H);
v_ic_JC_A(:,~L_S &  L_H) = v_ic_JC_max_A(:,~L_S &  L_H);
v_ic_JC_A(:,~L_S & ~L_H) = v_ic_JC_min_A(:,~L_S & ~L_H);


% 3. MUSCL scheme UP ====================================================

L_H                       = false(N_xy, 1);
L_H(H_IP_jc_m<=H_IP_jc_p) = true;

L_S                       = false(N_xy, 1);
L_S(S(ip_jc)<=S(ic_jc))   = true;

grad_S2_IP_jc = ((S(ic_jp) - S(ic_jm) + S(ip_jp) - S(ip_jm))/(4*dy)).^2 ...
              + ((S(ip_jc) - S(ic_jc))/dx).^2;
            
grad_S_IP_jc  = sqrt(grad_S2_IP_jc);

if THERMAL            
  T_IP_jc = 0.5*(T(:,ic_jc) + T(:,ip_jc));
else
  T_IP_jc = [];
end

[D_IP_jc_m_A, u_IP_jc_m_A] = D_vel_calc_new(H_IP_jc_m, grad_S_IP_jc, T_IP_jc, dS_dx_IP_jc, THERMAL);
[D_IP_jc_p_A, u_IP_jc_p_A] = D_vel_calc_new(H_IP_jc_p, grad_S_IP_jc, T_IP_jc, dS_dx_IP_jc, THERMAL);

D_IP_jc_m_C    = C_IP_jc.*H_IP_jc_m.^(m_SLIDE+1).*grad_S2_IP_jc.^((m_SLIDE-1)/2);
D_IP_jc_m      = D_IP_jc_m_A + D_IP_jc_m_C;

D_IP_jc_p_C    = C_IP_jc.*H_IP_jc_p.^(m_SLIDE+1).*grad_S2_IP_jc.^((m_SLIDE-1)/2);
D_IP_jc_p      = D_IP_jc_p_A + D_IP_jc_p_C;

L_min          = D_IP_jc_m<=D_IP_jc_p;

D_IP_jc_min          = D_IP_jc_p;
D_IP_jc_min(L_min)   = D_IP_jc_m(L_min);

D_IP_jc_min_A        = D_IP_jc_p_A;
D_IP_jc_min_A(L_min) = D_IP_jc_m_A(L_min);

D_IP_jc_min_C        = D_IP_jc_p_C;
D_IP_jc_min_C(L_min) = D_IP_jc_m_C(L_min);

H_IP_jc_min          = H_IP_jc_p;
H_IP_jc_min(L_min)   = H_IP_jc_m(L_min);

L_max                = D_IP_jc_m>=D_IP_jc_p;

D_IP_jc_max          = D_IP_jc_p;
D_IP_jc_max(L_max)   = D_IP_jc_m(L_max);

D_IP_jc_max_A        = D_IP_jc_p_A;
D_IP_jc_max_A(L_max) = D_IP_jc_m_A(L_max);

D_IP_jc_max_C        = D_IP_jc_p_C;
D_IP_jc_max_C(L_max) = D_IP_jc_m_C(L_max);

H_IP_jc_max          = H_IP_jc_p;
H_IP_jc_max(L_max)   = H_IP_jc_m(L_max);

u_IP_jc_min_A          = u_IP_jc_p_A;
u_IP_jc_min_A(:,L_min) = u_IP_jc_m_A(:,L_min);
u_IP_jc_max_A          = u_IP_jc_p_A;
u_IP_jc_max_A(:,L_max) = u_IP_jc_m_A(:,L_max);

D_IP_jc              = zeros(N_xy, 1);

D_IP_jc( L_S &  L_H) = D_IP_jc_min( L_S &  L_H);
D_IP_jc( L_S & ~L_H) = D_IP_jc_max( L_S & ~L_H);
D_IP_jc(~L_S &  L_H) = D_IP_jc_max(~L_S &  L_H);
D_IP_jc(~L_S & ~L_H) = D_IP_jc_min(~L_S & ~L_H);

D_IP_jc_A              = zeros(N_xy, 1);

D_IP_jc_A( L_S &  L_H) = D_IP_jc_min_A( L_S &  L_H);
D_IP_jc_A( L_S & ~L_H) = D_IP_jc_max_A( L_S & ~L_H);
D_IP_jc_A(~L_S &  L_H) = D_IP_jc_max_A(~L_S &  L_H);
D_IP_jc_A(~L_S & ~L_H) = D_IP_jc_min_A(~L_S & ~L_H);

D_IP_jc_C              = zeros(N_xy, 1);

D_IP_jc_C( L_S &  L_H) = D_IP_jc_min_C( L_S &  L_H);
D_IP_jc_C( L_S & ~L_H) = D_IP_jc_max_C( L_S & ~L_H);
D_IP_jc_C(~L_S &  L_H) = D_IP_jc_max_C(~L_S &  L_H);
D_IP_jc_C(~L_S & ~L_H) = D_IP_jc_min_C(~L_S & ~L_H);

H_IP_jc                = zeros(N_xy, 1);

H_IP_jc( L_S &  L_H)   = H_IP_jc_min( L_S &  L_H);
H_IP_jc( L_S & ~L_H)   = H_IP_jc_max( L_S & ~L_H);
H_IP_jc(~L_S &  L_H)   = H_IP_jc_max(~L_S &  L_H);
H_IP_jc(~L_S & ~L_H)   = H_IP_jc_min(~L_S & ~L_H);

u_IP_jc_A                = zeros(n_XI, N_xy);
u_IP_jc_A(:, L_S &  L_H) = u_IP_jc_min_A(:, L_S &  L_H);
u_IP_jc_A(:, L_S & ~L_H) = u_IP_jc_max_A(:, L_S & ~L_H);
u_IP_jc_A(:,~L_S &  L_H) = u_IP_jc_max_A(:,~L_S &  L_H);
u_IP_jc_A(:,~L_S & ~L_H) = u_IP_jc_min_A(:,~L_S & ~L_H);


% 4. MUSCL scheme DOWN ==================================================

L_H      = false(N_xy, 1);
L_H(H_IC_jc_m<=H_IC_jc_p) = true;

L_S      = false(N_xy, 1);
L_S(S(ic_jc)<=S(im_jc)) = true;

grad_S2_IC_jc = ((S(im_jp) - S(im_jm) + S(ic_jp) - S(ic_jm))/(4*dy)).^2 ...
              + ((S(ic_jc) - S(im_jc))/dx).^2;    
grad_S_IC_jc  = sqrt(grad_S2_IC_jc);            
           
if THERMAL  
  T_IC_jc = 0.5*(T(:,im_jc)+T(:,ic_jc));
else
  T_IC_jc = [];
end

[D_IC_jc_m_A, u_IC_jc_m_A] = D_vel_calc_new(H_IC_jc_m, grad_S_IC_jc, T_IC_jc, dS_dx_IC_jc, THERMAL);
[D_IC_jc_p_A, u_IC_jc_p_A] = D_vel_calc_new(H_IC_jc_p, grad_S_IC_jc, T_IC_jc, dS_dx_IC_jc, THERMAL);

D_IC_jc_m_C = C_IC_jc.*H_IC_jc_m.^(m_SLIDE+1).*grad_S2_IC_jc.^((m_SLIDE-1)/2);
D_IC_jc_m   = D_IC_jc_m_A + D_IC_jc_m_C;

D_IC_jc_p_C = C_IC_jc.*H_IC_jc_p.^(m_SLIDE+1).*grad_S2_IC_jc.^((m_SLIDE-1)/2);
D_IC_jc_p   = D_IC_jc_p_A + D_IC_jc_p_C;

L_min              = D_IC_jc_m<=D_IC_jc_p;

D_IC_jc_min        = D_IC_jc_p;
D_IC_jc_min(L_min) = D_IC_jc_m(L_min);

D_IC_jc_min_A        = D_IC_jc_p_A;
D_IC_jc_min_A(L_min) = D_IC_jc_m_A(L_min);

D_IC_jc_min_C        = D_IC_jc_p_C;
D_IC_jc_min_C(L_min) = D_IC_jc_m_C(L_min);

H_IC_jc_min          = H_IC_jc_p;
H_IC_jc_min(L_min)   = H_IC_jc_m(L_min);

L_max              = D_IC_jc_m>=D_IC_jc_p;

D_IC_jc_max        = D_IC_jc_p;
D_IC_jc_max(L_max) = D_IC_jc_m(L_max);

D_IC_jc_max_A        = D_IC_jc_p_A;
D_IC_jc_max_A(L_max) = D_IC_jc_m_A(L_max);

D_IC_jc_max_C        = D_IC_jc_p_C;
D_IC_jc_max_C(L_max) = D_IC_jc_m_C(L_max);

H_IC_jc_max          = H_IC_jc_p;
H_IC_jc_max(L_max)   = H_IC_jc_m(L_max);

u_IC_jc_min_A          = u_IC_jc_p_A;
u_IC_jc_min_A(:,L_min) = u_IC_jc_m_A(:,L_min);
u_IC_jc_max_A          = u_IC_jc_p_A;
u_IC_jc_max_A(:,L_max) = u_IC_jc_m_A(:,L_max);

D_IC_jc              = zeros(N_xy, 1);

D_IC_jc( L_S &  L_H) = D_IC_jc_min( L_S &  L_H);
D_IC_jc( L_S & ~L_H) = D_IC_jc_max( L_S & ~L_H);
D_IC_jc(~L_S &  L_H) = D_IC_jc_max(~L_S &  L_H);
D_IC_jc(~L_S & ~L_H) = D_IC_jc_min(~L_S & ~L_H);

D_IC_jc_A              = zeros(N_xy, 1);

D_IC_jc_A( L_S &  L_H) = D_IC_jc_min_A( L_S &  L_H);
D_IC_jc_A( L_S & ~L_H) = D_IC_jc_max_A( L_S & ~L_H);
D_IC_jc_A(~L_S &  L_H) = D_IC_jc_max_A(~L_S &  L_H);
D_IC_jc_A(~L_S & ~L_H) = D_IC_jc_min_A(~L_S & ~L_H);

D_IC_jc_C              = zeros(N_xy, 1);

D_IC_jc_C( L_S &  L_H) = D_IC_jc_min_C( L_S &  L_H);
D_IC_jc_C( L_S & ~L_H) = D_IC_jc_max_C( L_S & ~L_H);
D_IC_jc_C(~L_S &  L_H) = D_IC_jc_max_C(~L_S &  L_H);
D_IC_jc_C(~L_S & ~L_H) = D_IC_jc_min_C(~L_S & ~L_H);

H_IC_jc                = zeros(N_xy, 1);

H_IC_jc( L_S &  L_H)   = H_IC_jc_min( L_S &  L_H);
H_IC_jc( L_S & ~L_H)   = H_IC_jc_max( L_S & ~L_H);
H_IC_jc(~L_S &  L_H)   = H_IC_jc_max(~L_S &  L_H);
H_IC_jc(~L_S & ~L_H)   = H_IC_jc_min(~L_S & ~L_H);

u_IC_jc_A                = zeros(n_XI, N_xy);
u_IC_jc_A(:, L_S &  L_H) = u_IC_jc_min_A(:, L_S &  L_H);
u_IC_jc_A(:, L_S & ~L_H) = u_IC_jc_max_A(:, L_S & ~L_H);
u_IC_jc_A(:,~L_S &  L_H) = u_IC_jc_max_A(:,~L_S &  L_H);
u_IC_jc_A(:,~L_S & ~L_H) = u_IC_jc_min_A(:,~L_S & ~L_H);


% End ot MUSCL Scheme code sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(dB_dx_IC_jc)
  B_IC_jc     = 0.5*(B(ic_jc) + B(im_jc));
  B_IP_jc     = 0.5*(B(ic_jc) + B(ip_jc));
  B_ic_JC     = 0.5*(B(ic_jc) + B(ic_jm));
  B_ic_JP     = 0.5*(B(ic_jc) + B(ic_jp));
  
  dB_dx_IC_jc = (B(ic_jc)-B(im_jc))/dx;
  dB_dx_IP_jc = (B(ip_jc)-B(ic_jc))/dx;

  dB_dy_ic_JC = (B(ic_jc)-B(ic_jm))/dy;
  dB_dy_ic_JP = (B(ic_jp)-B(ic_jc))/dy;
end

% Calculate dS/dt and dH/dt for KBC calculations

Qx_IC_jc     = -dx^2*D_IC_jc.*dS_dx_IC_jc;
Qx_IP_jc     = -dx^2*D_IP_jc.*dS_dx_IP_jc;
Qy_ic_JC     = -dy^2*D_ic_JC.*dS_dy_ic_JC;
Qy_ic_JP     = -dy^2*D_ic_JP.*dS_dy_ic_JP;

Qx_IC_jc_C   = -dx^2*D_IC_jc_C.*dS_dx_IC_jc;
Qx_IP_jc_C   = -dx^2*D_IP_jc_C.*dS_dx_IP_jc;
Qy_ic_JC_C   = -dy^2*D_ic_JC_C.*dS_dy_ic_JC;
Qy_ic_JP_C   = -dy^2*D_ic_JP_C.*dS_dy_ic_JP;

u_IC_jc_C        = Qx_IC_jc_C./H_IC_jc;

L_H              = H_IC_jc<=0;
u_IC_jc_C(L_H)   = 0;
u_IC_jc_A(:,L_H) = 0;
u_IC_jc          = u_IC_jc_A + repmat(u_IC_jc_C', n_XI, 1);

u_IP_jc_C        = Qx_IP_jc_C./H_IP_jc;

L_H              = H_IP_jc<=0;
u_IP_jc_C(L_H)   = 0;
u_IP_jc_A(:,L_H) = 0;
u_IP_jc          = u_IP_jc_A + repmat(u_IP_jc_C', n_XI, 1);

v_ic_JC_C        = Qy_ic_JC_C./H_ic_JC;

L_H              = H_ic_JC<=0;
v_ic_JC_C(L_H)   = 0;
v_ic_JC_A(:,L_H) = 0;
v_ic_JC          = v_ic_JC_A + repmat(v_ic_JC_C', n_XI, 1);

v_ic_JP_C        = Qy_ic_JP_C./H_ic_JP;

L_H              = H_ic_JP<=0;
v_ic_JP_C(L_H)   = 0;
v_ic_JP_A(:,L_H) = 0;
v_ic_JP          = v_ic_JP_A + repmat(v_ic_JP_C', n_XI, 1);

u_ic_jc   = 0.5*(u_IC_jc + u_IP_jc);
v_ic_jc   = 0.5*(v_ic_JC + v_ic_JP);

L_H            = H<=0;
u_ic_jc(:,L_H) = 0;
v_ic_jc(:,L_H) = 0;

L_H = H<=0;

dH_dx_IC_jc = (H(ic_jc) - H(im_jc))/dx;
dH_dx_IP_jc = (H(ip_jc) - H(ic_jc))/dx;

dH_dy_ic_JC = (H(ic_jc) - H(ic_jm))/dy;
dH_dy_ic_JP = (H(ic_jp) - H(ic_jc))/dy;

dH_dx       = 0.5*(dH_dx_IC_jc + dH_dx_IP_jc);
dH_dy       = 0.5*(dH_dy_ic_JC + dH_dy_ic_JP);

dB_dx       = 0.5*(dB_dx_IC_jc + dB_dx_IP_jc);
dB_dy       = 0.5*(dB_dy_ic_JC + dB_dy_ic_JP);

dS_dx       = dB_dx + dH_dx;
dS_dy       = dB_dy + dH_dy;

H_use       = max(H, 1);  % special treatment for ice thinner than 1 m
H_use(L_H)  = 0;

d_XI_dx     = -repmat(dB_dx./H_use, 1, n_XI)' - kron((dH_dx./H_use)', XI);
d_XI_dy     = -repmat(dB_dy./H_use, 1, n_XI)' - kron((dH_dy./H_use)', XI);
d_XI_dz     = +repmat(1./H_use, 1, n_XI)';
  
XI_KP            = 0.5*(XI(2:n_XI) + XI(1:n_XI-1));

[w_ic_jc, du_dx_ic_jc_KP, dv_dy_ic_jc_KP, dw_dz_ic_jc_KP, div_v_ic_jc_KP] = w_elisa(u_IC_jc, u_IP_jc, v_ic_JC, v_ic_JP, dB_dx, dB_dy, dH_dx, dH_dy, H_use, H_IC_jc, H_IP_jc, H_ic_JC, H_ic_JP, XI);

% The above values for the velocity divergence terms are perfect -- nothing compares
% Thus I use these to compute the same components on the kc (non-staggered) grid
% in order to ensure that these points also obey the incompressibility condition
% and then bootstrap in values at B and S using other means

du_ic_jc_dx_B = (1./H_use).*((u_IP_jc(1,:)'.*H_IP_jc - u_IC_jc(1,:)'.*H_IC_jc)/dx ...
  -(-3*u_ic_jc(1,:)'.*dB_dx + 4*u_ic_jc(2,:)'.*(dB_dx + XI(2)*dH_dx) - u_ic_jc(3,:)'.*(dB_dx + XI(3)*dH_dx))/(2*d_XI));

dv_ic_jc_dy_B = (1./H_use).*((v_ic_JP(1,:)'.*H_ic_JP - v_ic_JC(1,:)'.*H_ic_JC)/dy ...
  -(-3*v_ic_jc(1,:)'.*dB_dy + 4*v_ic_jc(2,:)'.*(dB_dy + XI(2)*dH_dy) - v_ic_jc(3,:)'.*(dB_dy + XI(3)*dH_dy))/(2*d_XI));

dw_ic_jc_dz_B = (1./H_use).*(-3*w_ic_jc(1,:)' + 4*w_ic_jc(2,:)' - w_ic_jc(3,:)')/(2*d_XI);

du_ic_jc_dx_S = (1./H_use).*((u_IP_jc(n_XI,:)'.*H_IP_jc - u_IC_jc(n_XI,:)'.*H_IC_jc)/dx ...
    -(+u_ic_jc(n_XI-2,:)'.*(dB_dx + XI(n_XI-2)*dH_dx) - 4*u_ic_jc(n_XI-1,:)'.*(dB_dx + XI(n_XI-1)*dH_dx) + 3*u_ic_jc(n_XI,:)'.*(dB_dx + XI(n_XI)*dH_dx))/(2*d_XI));
  
dv_ic_jc_dy_S = (1./H_use).*((v_ic_JP(n_XI,:)'.*H_ic_JP - v_ic_JC(n_XI,:)'.*H_ic_JC)/dy ...
    -(+v_ic_jc(n_XI-2,:)'.*(dB_dy + XI(n_XI-2)*dH_dy) - 4*v_ic_jc(n_XI-1,:)'.*(dB_dy + XI(n_XI-1)*dH_dy) + 3*v_ic_jc(n_XI,:)'.*(dB_dy + XI(n_XI)*dH_dy))/(2*d_XI));

dw_ic_jc_dz_S = (1./H).*(w_ic_jc(n_XI-2,:)' - 4*w_ic_jc(n_XI-1,:)' + 3*w_ic_jc(n_XI,:)')/(2*d_XI);

du_ic_jc_dx = cat(1, du_ic_jc_dx_B', 0.5*(du_dx_ic_jc_KP(2:n_XI-1,:) + du_dx_ic_jc_KP(1:n_XI-2,:)), du_ic_jc_dx_S');
dv_ic_jc_dy = cat(1, dv_ic_jc_dy_B', 0.5*(dv_dy_ic_jc_KP(2:n_XI-1,:) + dv_dy_ic_jc_KP(1:n_XI-2,:)), dv_ic_jc_dy_S');
dw_ic_jc_dz = cat(1, dw_ic_jc_dz_B', 0.5*(dw_dz_ic_jc_KP(2:n_XI-1,:) + dw_dz_ic_jc_KP(1:n_XI-2,:)), dw_ic_jc_dz_S');

% Calculate the error tweak terms to enforce incompressibility at B and S surfaces

div_v_ic_jc = du_ic_jc_dx + dv_ic_jc_dy + dw_ic_jc_dz;
div_v2      = du_ic_jc_dx.^2 + dv_ic_jc_dy.^2 + dw_ic_jc_dz.^2;

LAMBDA_B    = 2*div_v_ic_jc(1,:)./div_v2(1,:);
LAMBDA_S    = 2*div_v_ic_jc(n_XI,:)./div_v2(n_XI,:);

du_ic_jc_dx(1,:)    = (1 - 0.5*LAMBDA_B.*du_ic_jc_dx(1,:)).*du_ic_jc_dx(1,:);
dv_ic_jc_dy(1,:)    = (1 - 0.5*LAMBDA_B.*dv_ic_jc_dy(1,:)).*dv_ic_jc_dy(1,:);
dw_ic_jc_dz(1,:)    = (1 - 0.5*LAMBDA_B.*dw_ic_jc_dz(1,:)).*dw_ic_jc_dz(1,:);

du_ic_jc_dx(n_XI,:) = (1 - 0.5*LAMBDA_S.*du_ic_jc_dx(n_XI,:)).*du_ic_jc_dx(n_XI,:);
dv_ic_jc_dy(n_XI,:) = (1 - 0.5*LAMBDA_S.*dv_ic_jc_dy(n_XI,:)).*dv_ic_jc_dy(n_XI,:);
dw_ic_jc_dz(n_XI,:) = (1 - 0.5*LAMBDA_S.*dw_ic_jc_dz(n_XI,:)).*dw_ic_jc_dz(n_XI,:);

du_ic_jc_dx(:,L_H)  = 0;
dv_ic_jc_dy(:,L_H)  = 0;
dw_ic_jc_dz(:,L_H)  = 0;

w_ic_jc(:,L_H)      = 0;

% The following XI derivatives are NOT used for computing components for flow divergence
% but are use to calculate non-diagonal component of velocity gradient tensor

du_ic_jc_d_XI_B  = (-3*u_ic_jc(1,:)' + 4*u_ic_jc(2,:)' - u_ic_jc(3,:)')/(2*d_XI);
dv_ic_jc_d_XI_B  = (-3*v_ic_jc(1,:)' + 4*v_ic_jc(2,:)' - v_ic_jc(3,:)')/(2*d_XI);
dw_ic_jc_d_XI_B  = (-3*w_ic_jc(1,:)' + 4*w_ic_jc(2,:)' - w_ic_jc(3,:)')/(2*d_XI);

du_ic_jc_d_XI_S  = (u_ic_jc(n_XI-2,:)' - 4*u_ic_jc(n_XI-1,:)' + 3*u_ic_jc(n_XI,:)')/(2*d_XI);
dv_ic_jc_d_XI_S  = (v_ic_jc(n_XI-2,:)' - 4*v_ic_jc(n_XI-1,:)' + 3*v_ic_jc(n_XI,:)')/(2*d_XI);
dw_ic_jc_d_XI_S  = (w_ic_jc(n_XI-2,:)' - 4*w_ic_jc(n_XI-1,:)' + 3*w_ic_jc(n_XI,:)')/(2*d_XI);

du_ic_jc_d_XI    = cat(1, du_ic_jc_d_XI_B', (u_ic_jc(3:n_XI,:)-u_ic_jc(1:n_XI-2,:))/(2*d_XI), du_ic_jc_d_XI_S');
dv_ic_jc_d_XI    = cat(1, dv_ic_jc_d_XI_B', (v_ic_jc(3:n_XI,:)-v_ic_jc(1:n_XI-2,:))/(2*d_XI), dv_ic_jc_d_XI_S');
dw_ic_jc_d_XI    = cat(1, dw_ic_jc_d_XI_B', (w_ic_jc(3:n_XI,:)-w_ic_jc(1:n_XI-2,:))/(2*d_XI), dw_ic_jc_d_XI_S');

% Fill out the remainining terms of velocity gradient tensor (in Cartesian system)

du_ic_jc_dy = (u_IP_jc(:,ic_jp)+u_IC_jc(:,ic_jp)-u_IP_jc(:,ic_jm)-u_IC_jc(:,ic_jm))/(4*dy) +du_ic_jc_d_XI.*d_XI_dy;
du_ic_jc_dz = du_ic_jc_d_XI.*d_XI_dz;

dv_ic_jc_dx = (v_ic_JP(:,ip_jc)+v_ic_JC(:,ip_jc)-v_ic_JP(:,im_jc)-v_ic_JC(:,im_jc))/(4*dx) + dv_ic_jc_d_XI.*d_XI_dx;
dv_ic_jc_dz = dv_ic_jc_d_XI.*d_XI_dz;

dw_ic_jc_dx = (w_ic_jc(:,ip_jc)-w_ic_jc(:,im_jc))/(2*dx) + dw_ic_jc_d_XI.*d_XI_dx;
dw_ic_jc_dy = (w_ic_jc(:,ic_jp)-w_ic_jc(:,ic_jm))/(2*dy) + dw_ic_jc_d_XI.*d_XI_dy;
  
du_ic_jc_dy(:,L_H) = 0;
du_ic_jc_dz(:,L_H) = 0;

dv_ic_jc_dx(:,L_H) = 0;
dv_ic_jc_dz(:,L_H) = 0;

dw_ic_jc_dx(:,L_H) = 0;
dw_ic_jc_dy(:,L_H) = 0;

V.u_ic_jc   = u_ic_jc;
V.v_ic_jc   = v_ic_jc;
V.w_ic_jc   = w_ic_jc;

V.u_IC_jc_A = u_IC_jc_A;
V.u_IC_jc_C = u_IC_jc_C;
V.u_IP_jc_A = u_IP_jc_A;
V.u_IP_jc_C = u_IP_jc_C;

V.v_ic_JC_A = v_ic_JC_A;
V.v_ic_JC_C = v_ic_JC_C;
V.v_ic_JP_A = v_ic_JP_A;
V.v_ic_JP_C = v_ic_JP_C;

V.grid      = 'XI';

L          = zeros(3,3, n_XI, N_xy);

L(1,1,:,:) = du_ic_jc_dx;
L(1,2,:,:) = du_ic_jc_dy;
L(1,3,:,:) = du_ic_jc_dz;

L(2,1,:,:) = dv_ic_jc_dx;
L(2,2,:,:) = dv_ic_jc_dy;
L(2,3,:,:) = dv_ic_jc_dz;

L(3,1,:,:) = dw_ic_jc_dx;
L(3,2,:,:) = dw_ic_jc_dy;
L(3,3,:,:) = dw_ic_jc_dz;

D          = 0.5*(L + permute(L, [2 1 3 4]));

% First invariant of strain-rate (1/yr)

E1       = squeeze(D(1,1,:,:) + D(2,2,:,:) + D(3,3,:,:));
   
% Terms of the antisymmetrical spin rate tensor

W       = 0.5*(L - permute(L, [2 1 3 4]));

% Second invariant of strain-rate (1/yr^2)

E2   = 0.5*(D(1,1,:,:).^2 + D(2,2,:,:).^2 + D(3,3,:,:).^2) + D(1,2,:,:).^2 + D(1,3,:,:).^2 + D(2,3,:,:).^2;
E2   = squeeze(E2);

if THERMAL
  A         = A_0*exp(-Q_creep./(R_gas*T));
  s2        = (E2./A.^2).^(1/n_GLEN);
  TAU_eff   = sqrt(s2);
  PHI_ic_jc = 2*A.*TAU_eff.^(n_GLEN+1)/YRSEC;        % units are W/m^3
else  
  PHI_ic_jc = 0;
end

E = Build_E_Structure(T, D, W, N_xy);

H_struct.IC_jc = H_IC_jc;
H_struct.IP_jc = H_IP_jc;
H_struct.ic_JC = H_ic_JC;
H_struct.ic_JP = H_ic_JP;

