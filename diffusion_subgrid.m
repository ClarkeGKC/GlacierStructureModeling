function [D_x_dn, D_x_up, D_y_dn, D_y_up] = diffusion_subgrid(S, B, sub)

% Limited MUSCL scheme diffusion calculation

global dx dy A_const C_const n_GLEN m_SLIDE N_sub
global IC_JC

dx      = IC_JC.dx;
dy      = IC_JC.dy;

A_tilde = A_const/dx^2;
C_tilde = C_const/dx^2;

H = S-B;

% MUSCL scheme UP

H_y_up_m = H(sub.ic_jc) + 0.5*minmod( H(sub.ic_jp)-H(sub.ic_jc), H(sub.ic_jc)-H(sub.ic_jm));
H_y_up_p = H(sub.ic_jp) - 0.5*minmod( H(sub.ic_jpp)-H(sub.ic_jp), H(sub.ic_jp)-H(sub.ic_jc));

s2_y_grad_up = ((S(sub.ip_jc) - S(sub.im_jc) + S(sub.ip_jp) - S(sub.im_jp))/(4*dx)).^2 ...
             + ((S(sub.ic_jp) - S(sub.ic_jc))/dy).^2;
          
D_y_up_m   = A_tilde*H_y_up_m.^(n_GLEN+2) .*s2_y_grad_up.^((n_GLEN-1)/2) ...
           + C_tilde*H_y_up_m.^(m_SLIDE+1).*s2_y_grad_up.^((m_SLIDE-1)/2);
D_y_up_p   = A_tilde*H_y_up_p.^(n_GLEN+2) .*s2_y_grad_up.^((n_GLEN-1)/2) ...
           + C_tilde*H_y_up_p.^(m_SLIDE+1).*s2_y_grad_up.^((m_SLIDE-1)/2);

D_y_up_min = min(D_y_up_m, D_y_up_p);
D_y_up_max = max(D_y_up_m, D_y_up_p);

D_y_up     = zeros(N_sub, 1);

% Use local slope to identify upstream direction

D_y_up(S(sub.ic_jp)<=S(sub.ic_jc) & H_y_up_m<=H_y_up_p) = D_y_up_min(S(sub.ic_jp)<=S(sub.ic_jc) & H_y_up_m<=H_y_up_p);
D_y_up(S(sub.ic_jp)<=S(sub.ic_jc) & H_y_up_m> H_y_up_p) = D_y_up_max(S(sub.ic_jp)<=S(sub.ic_jc) & H_y_up_m> H_y_up_p);
D_y_up(S(sub.ic_jp)> S(sub.ic_jc) & H_y_up_m<=H_y_up_p) = D_y_up_max(S(sub.ic_jp)> S(sub.ic_jc) & H_y_up_m<=H_y_up_p);
D_y_up(S(sub.ic_jp)> S(sub.ic_jc) & H_y_up_m> H_y_up_p) = D_y_up_min(S(sub.ic_jp)> S(sub.ic_jc) & H_y_up_m> H_y_up_p);

% MUSCL scheme DOWN

H_y_dn_m = H(sub.ic_jm) + 0.5*minmod( H(sub.ic_jc)-H(sub.ic_jm), H(sub.ic_jm)-H(sub.ic_jmm));
H_y_dn_p = H(sub.ic_jc) - 0.5*minmod( H(sub.ic_jp)-H(sub.ic_jc), H(sub.ic_jc)-H(sub.ic_jm));

s2_y_grad_dn = ((S(sub.ip_jm) - S(sub.im_jm) + S(sub.ip_jc) - S(sub.im_jc))/(4*dx)).^2 ...
             + ((S(sub.ic_jc) - S(sub.ic_jm))/dy).^2;
         
D_y_dn_m = A_tilde*H_y_dn_m.^(n_GLEN+2) .*s2_y_grad_dn.^((n_GLEN-1)/2) ...
         + C_tilde*H_y_dn_m.^(m_SLIDE+1).*s2_y_grad_dn.^((m_SLIDE-1)/2);
D_y_dn_p = A_tilde*H_y_dn_p.^(n_GLEN+2) .*s2_y_grad_dn.^((n_GLEN-1)/2);
         + C_tilde*H_y_dn_p.^(m_SLIDE+1).*s2_y_grad_dn.^((m_SLIDE-1)/2);

D_y_dn_min = min(D_y_dn_m, D_y_dn_p);
D_y_dn_max = max(D_y_dn_m, D_y_dn_p);

D_y_dn     = zeros(N_sub, 1);

% Use local slope to identify upstream direction

D_y_dn(S(sub.ic_jc)<=S(sub.ic_jm) & H_y_dn_m<=H_y_dn_p) = D_y_dn_min(S(sub.ic_jc)<=S(sub.ic_jm) & H_y_dn_m<=H_y_dn_p);
D_y_dn(S(sub.ic_jc)<=S(sub.ic_jm) & H_y_dn_m> H_y_dn_p) = D_y_dn_max(S(sub.ic_jc)<=S(sub.ic_jm) & H_y_dn_m> H_y_dn_p);
D_y_dn(S(sub.ic_jc)> S(sub.ic_jm) & H_y_dn_m<=H_y_dn_p) = D_y_dn_max(S(sub.ic_jc)> S(sub.ic_jm) & H_y_dn_m<=H_y_dn_p);
D_y_dn(S(sub.ic_jc)> S(sub.ic_jm) & H_y_dn_m> H_y_dn_p) = D_y_dn_min(S(sub.ic_jc)> S(sub.ic_jm) & H_y_dn_m> H_y_dn_p);

% MUSCL scheme UP

H_x_up_m = H(sub.ic_jc) + 0.5*minmod( H(sub.ip_jc) -H(sub.ic_jc), H(sub.ic_jc)-H(sub.im_jc));
H_x_up_p = H(sub.ip_jc) - 0.5*minmod( H(sub.ipp_jc)-H(sub.ip_jc), H(sub.ip_jc)-H(sub.ic_jc));

s2_x_grad_up = ((S(sub.ic_jp) - S(sub.ic_jm) + S(sub.ip_jp) - S(sub.ip_jm))/(4*dy)).^2 ...
             + ((S(sub.ip_jc) - S(sub.ic_jc))/dx).^2;
          
D_x_up_m   = A_tilde*H_x_up_m.^(n_GLEN+2) .*s2_x_grad_up.^((n_GLEN-1)/2) ...
           + C_tilde*H_x_up_m.^(m_SLIDE+1).*s2_x_grad_up.^((m_SLIDE-1)/2);
D_x_up_p   = A_tilde*H_x_up_p.^(n_GLEN+2) .*s2_x_grad_up.^((n_GLEN-1)/2) ...
           + C_tilde*H_x_up_p.^(m_SLIDE+1).*s2_x_grad_up.^((m_SLIDE-1)/2);

D_x_up_min = min(D_x_up_m, D_x_up_p);
D_x_up_max = max(D_x_up_m, D_x_up_p);

D_x_up     = zeros(N_sub, 1);

D_x_up(S(sub.ip_jc)<=S(sub.ic_jc) & H_x_up_m<=H_x_up_p) = D_x_up_min(S(sub.ip_jc)<=S(sub.ic_jc) & H_x_up_m<=H_x_up_p);
D_x_up(S(sub.ip_jc)<=S(sub.ic_jc) & H_x_up_m> H_x_up_p) = D_x_up_max(S(sub.ip_jc)<=S(sub.ic_jc) & H_x_up_m> H_x_up_p);
D_x_up(S(sub.ip_jc)> S(sub.ic_jc) & H_x_up_m<=H_x_up_p) = D_x_up_max(S(sub.ip_jc)> S(sub.ic_jc) & H_x_up_m<=H_x_up_p);
D_x_up(S(sub.ip_jc)> S(sub.ic_jc) & H_x_up_m> H_x_up_p) = D_x_up_min(S(sub.ip_jc)> S(sub.ic_jc) & H_x_up_m> H_x_up_p);

% MUSCL scheme DOWN

H_x_dn_m = H(sub.im_jc) + 0.5*minmod( H(sub.ic_jc)-H(sub.im_jc), H(sub.im_jc)-H(sub.imm_jc));
H_x_dn_p = H(sub.ic_jc) - 0.5*minmod( H(sub.ip_jc)-H(sub.ic_jc), H(sub.ic_jc)-H(sub.im_jc));

s2_x_grad_dn = ((S(sub.im_jp) - S(sub.im_jm) + S(sub.ic_jp) - S(sub.ic_jm))/(4*dy)).^2 ...
             + ((S(sub.ic_jc) - S(sub.im_jc))/dx).^2;
          
D_x_dn_m = A_tilde*H_x_dn_m.^(n_GLEN+2) .*s2_x_grad_dn.^((n_GLEN-1)/2) ...
         + C_tilde*H_x_dn_m.^(m_SLIDE+1).*s2_x_grad_dn.^((m_SLIDE-1)/2);
D_x_dn_p = A_tilde*H_x_dn_p.^(n_GLEN+2) .*s2_x_grad_dn.^((n_GLEN-1)/2) ...
         + C_tilde*H_x_dn_p.^(m_SLIDE+1).*s2_x_grad_dn.^((m_SLIDE-1)/2);

D_x_dn_min = min(D_x_dn_m, D_x_dn_p);
D_x_dn_max = max(D_x_dn_m, D_x_dn_p);

D_x_dn  = zeros(N_sub, 1);

D_x_dn(S(sub.ic_jc)<=S(sub.im_jc) & H_x_dn_m<=H_x_dn_p) = D_x_dn_min(S(sub.ic_jc)<=S(sub.im_jc) & H_x_dn_m<=H_x_dn_p);
D_x_dn(S(sub.ic_jc)<=S(sub.im_jc) & H_x_dn_m> H_x_dn_p) = D_x_dn_max(S(sub.ic_jc)<=S(sub.im_jc) & H_x_dn_m> H_x_dn_p);
D_x_dn(S(sub.ic_jc)> S(sub.im_jc) & H_x_dn_m<=H_x_dn_p) = D_x_dn_max(S(sub.ic_jc)> S(sub.im_jc) & H_x_dn_m<=H_x_dn_p);
D_x_dn(S(sub.ic_jc)> S(sub.im_jc) & H_x_dn_m> H_x_dn_p) = D_x_dn_min(S(sub.ic_jc)> S(sub.im_jc) & H_x_dn_m> H_x_dn_p);

% This stopped working when I switched to MATLAB 2015b

% D_x_dn = D_x_dn';
% D_x_up = D_x_up';
% D_y_dn = D_y_dn';
% D_y_up = D_y_up';
