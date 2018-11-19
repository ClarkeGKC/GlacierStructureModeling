function [w_ic_jc, du_dx_ic_jc_KP, dv_dy_ic_jc_KP, dw_dz_ic_jc_KP, div_v_ic_jc_KP]  = w_elisa(u_IC_jc, u_IP_jc, v_ic_JC, v_ic_JP, dB_dx, dB_dy, dH_dx, dH_dy, H, H_IC_jc, H_IP_jc, H_ic_JC, H_ic_JP, XI)

global nx ny dx dy d_XI n_XI

u_ic_jc_kc = 0.5*(u_IC_jc + u_IP_jc);
v_ic_jc_kc = 0.5*(v_ic_JC + v_ic_JP);

u_IC_jc_KP = 0.5*(u_IC_jc(1:n_XI-1,:) + u_IC_jc(2:n_XI,:));
u_IP_jc_KP = 0.5*(u_IP_jc(1:n_XI-1,:) + u_IP_jc(2:n_XI,:));

v_ic_JC_KP = 0.5*(v_ic_JC(1:n_XI-1,:) + v_ic_JC(2:n_XI,:));
v_ic_JP_KP = 0.5*(v_ic_JP(1:n_XI-1,:) + v_ic_JP(2:n_XI,:));

u_BH_XI_x    = u_ic_jc_kc.*(repmat(dB_dx, 1, n_XI)' + kron(XI, dH_dx'));
v_BH_XI_y    = v_ic_jc_kc.*(repmat(dB_dy, 1, n_XI)' + kron(XI, dH_dy'));

uH_KP_term    = u_IP_jc_KP.*repmat(H_IP_jc, 1, n_XI-1)' - u_IC_jc_KP.*repmat(H_IC_jc, 1, n_XI-1)';
vH_KP_term    = v_ic_JP_KP.*repmat(H_ic_JP, 1, n_XI-1)' - v_ic_JC_KP.*repmat(H_ic_JC, 1, n_XI-1)';

w_ic_jc       = zeros(size(u_IC_jc));
w_ic_jc(1,:) = u_ic_jc_kc(1,:)'.*dB_dx + v_ic_jc_kc(1,:)'.*dB_dy;

for k=1:n_XI-1
  w_ic_jc(k+1,:) = w_ic_jc(k,:) +  u_BH_XI_x(k+1,:) + v_BH_XI_y(k+1,:) - u_BH_XI_x(k,:) - v_BH_XI_y(k,:) ...
                 - d_XI*uH_KP_term(k, :)/dx - d_XI*vH_KP_term(k,:)/dy;
end  

du_dx_ic_jc_KP = repmat(1./H, 1, n_XI-1)'.*(uH_KP_term/dx - (u_BH_XI_x(2:n_XI,:) - u_BH_XI_x(1:n_XI-1,:))/d_XI);
dv_dy_ic_jc_KP = repmat(1./H, 1, n_XI-1)'.*(vH_KP_term/dy - (v_BH_XI_y(2:n_XI,:) - v_BH_XI_y(1:n_XI-1,:))/d_XI);
dw_dz_ic_jc_KP = repmat(1./H, 1, n_XI-1)'.*(w_ic_jc(2:n_XI,:)-w_ic_jc(1:n_XI-1,:))/d_XI;

div_v_ic_jc_KP = du_dx_ic_jc_KP + dv_dy_ic_jc_KP + dw_dz_ic_jc_KP;
