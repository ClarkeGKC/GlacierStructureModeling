function [T_zeta_out, t_end] = step_T(T_zeta, B, H, u, v, w, PHI, dt, T_S, q_G, H_t, t)

% Fully implicit version

global RHO BETA k_0 k_1 c_0 c_1 YRSEC d_ZETA ZETA n_ZETA T_KELVIN VERBOSITY
global IC_JC km_T kc_T kp_T
global H_min dt_max_T
global B_x_E B_x_W B_y_N B_y_S     % Not shared but saved between time steps

ic_jc = IC_JC.ic_jc;
ic_jp = IC_JC.ic_jp;
ic_jm = IC_JC.ic_jm;
im_jc = IC_JC.im_jc;
ip_jc = IC_JC.ip_jc;
dx    = IC_JC.dx;
dy    = IC_JC.dy;
nx    = IC_JC.nx;
ny    = IC_JC.ny;
N     = IC_JC.N;

% First treat the cells that have thin ice cover and set full column
% temperature to ice surface temperature (or better still to a linear
% function of depth)

if strcmp(T_zeta.grid, 'ZETA')==0
  error('step_T(): Input temperature structure should be ZETA grid not "%s"', T_zeta.grid)
end  

T = T_zeta;

ic_jc_0 = ic_jc(H<=0);
T.A(:,ic_jc_0) = repmat(T_S(ic_jc_0), 1 , n_ZETA)';  % Set ice-free temperatures to T_S

ic_jc_L = ic_jc(H>0 & H<=H_min);                      % index to select for thin-ice-covered cells
k_L     = k_0*exp(-k_1*T_S(ic_jc_L))*YRSEC;           % Note units are J/(m K yr)

T_B_L   = T_S(ic_jc_L)+H(ic_jc_L).*q_G(ic_jc_L)./k_L;

T.A(:,ic_jc_L) = transpose(repmat(T_B_L, 1, n_ZETA)) + kron(ZETA, (T_S(ic_jc_L)-T_B_L)');  % Set "thin" ice to linear gradient

% Next solve for temperature in the cells that have ice thickness greater than H_min_T_solver. 
% If necessary do this by multiple passes with time step dt_T<dt

ic_jc_I = ic_jc(H>H_min);              % Index to select for thick-ice-covered cells

if isempty(ic_jc_I) || numel(ic_jc_I)<10
  T_zeta_out = T_zeta;
  t_end = t+dt;
  return
else
  T_zeta_out = T;                           % For all thick-ice points will be overwritten by solver output. This just sets the default values
end  

nt_T = ceil(dt/dt_max_T);
dt_T = dt/nt_T;

ic_jp_I = ic_jp(ic_jc_I);
ic_jm_I = ic_jm(ic_jc_I);
im_jc_I = im_jc(ic_jc_I);
ip_jc_I = ip_jc(ic_jc_I);

H_I     = H(ic_jc_I);
H_t_I   = H_t(ic_jc_I);
T_S_I   = T_S(ic_jc_I);

T_S_I   = min(T_S_I, T_KELVIN);      % Restrict surface ice temperate to 0 C or lower

T_MP_I  = T_KELVIN-BETA*H_I;

N_xy_I  = numel(ic_jc_I);
N_I     = n_ZETA*N_xy_I;

T_I   = T.A(:,ic_jc_I);
u_I   = u.A(:,ic_jc_I);
v_I   = v.A(:,ic_jc_I);
w_I   = w.A(:,ic_jc_I);
PHI_I = PHI.A(:, ic_jc_I);

T_ic_jp_I = T.A(:,ic_jp_I);
T_ic_jm_I = T.A(:,ic_jm_I);
T_im_jc_I = T.A(:,im_jc_I);
T_ip_jc_I = T.A(:,ip_jc_I);

c_I     = c_0 + c_1*T_I;
k_I     = k_0*exp(-k_1*T_I)*YRSEC;
KAPPA_I = k_I./(RHO*c_I);
dk_dT_I = -k_1*k_I;

u_I     = reshape(u_I, N_I, 1);
v_I     = reshape(v_I, N_I, 1);
w_I     = reshape(w_I, N_I, 1);   

% UP_X, UP_Y, UP_Z are upwinding flags (UP==Upsilon in the TeX document T_notes_v200.tex) 
% UP_X is true when the flow is aligned with the positive X axis 
% similarly for the other flag
 
UP_x    = u_I>0;    % Eastward flow
UP_y    = v_I>0;    % Northward flow
UP_z    = w_I>0;    % Upward flow

u_I     = reshape(u_I, n_ZETA, N_xy_I);
v_I     = reshape(v_I, n_ZETA, N_xy_I);
w_I     = reshape(w_I, n_ZETA, N_xy_I);

UP_x    = reshape(UP_x, n_ZETA, N_xy_I);
UP_y    = reshape(UP_y, n_ZETA, N_xy_I);
UP_z    = reshape(UP_z, n_ZETA, N_xy_I);

HH_I    = repmat(H_I, 1, n_ZETA)';

% Note: Moving A, B, C, evaluations outside the for loop makes them constant 
% over the full loop time. Thus, for example, any change in bottom BC will
% not affect these tridiagonal coefficients

if isempty(B_x_E)
  B_x_E   = (B(ic_jc) - B(im_jc))/dx;   % Upwinded derivative for Eastward flow
  B_x_W   = (B(ip_jc) - B(ic_jc))/dx;   % Upwinded derivative for Westward flow

  B_y_N   = (B(ic_jc) - B(ic_jm))/dy;   % Upwinded derivative for Northward flow
  B_y_S   = (B(ic_jp) - B(ic_jc))/dy;   % Upwinded derivative for Southward flow
end

H_x_E   = (H(ic_jc) - H(im_jc))/dx;
H_x_W   = (H(ip_jc) - H(ic_jc))/dx;

H_y_N   = (H(ic_jc) - H(ic_jm))/dy;
H_y_S   = (H(ic_jp) - H(ic_jc))/dy;

B_x_E_I = B_x_E(ic_jc_I);
B_x_W_I = B_x_W(ic_jc_I);
B_y_N_I = B_y_N(ic_jc_I);
B_y_S_I = B_y_S(ic_jc_I);

H_x_E_I = H_x_E(ic_jc_I);
H_x_W_I = H_x_W(ic_jc_I);
H_y_N_I = H_y_N(ic_jc_I);
H_y_S_I = H_y_S(ic_jc_I);

% Use surface velocity to determinine upwinding of X and Y components

B_x_I   = UP_x(n_ZETA,:)'.*B_x_E_I + ~UP_x(n_ZETA,:)'.*B_x_W_I;
B_y_I   = UP_y(n_ZETA,:)'.*B_y_N_I + ~UP_y(n_ZETA,:)'.*B_y_S_I;

H_x_I   = UP_x(n_ZETA,:)'.*H_x_E_I + ~UP_x(n_ZETA,:)'.*H_x_W_I;
H_y_I   = UP_y(n_ZETA,:)'.*H_y_N_I + ~UP_y(n_ZETA,:)'.*H_y_S_I;

F          = u_I.*(repmat(B_x_I, 1, n_ZETA)' + kron(ZETA, H_x_I')) ...
           + v_I.*(repmat(B_y_I, 1, n_ZETA)' + kron(ZETA, H_y_I')) ...
           - w_I + kron(ZETA, H_t_I');

A          = (dt_T/d_ZETA)*UP_z.*F./HH_I - (dt_T/d_ZETA^2)*KAPPA_I./HH_I.^2;
B_part     = (dt_T/d_ZETA)*F./HH_I;
B          = 1 + ~UP_z.*B_part - UP_z.*B_part + 2*(dt_T/d_ZETA^2)*KAPPA_I./HH_I.^2;
C          = -(dt_T/d_ZETA)*~UP_z.*F./HH_I - (dt_T/d_ZETA^2)*KAPPA_I./HH_I.^2;

A(n_ZETA,:) = 0;
B(n_ZETA,:) = 1;
C(n_ZETA,:) = 0;

q_G_I  = q_G(ic_jc_I);

% Seems that d_ZETA.*H_I should be divided rather than multiplied in above expression

q_B_I = -k_I(1,:)'.*(T_I(2,:)-T_I(1,:))'./(H_I*d_ZETA);  % Used for bottom melting test (below)

L_melt = T_I(1,:)'>=T_MP_I & q_G_I>=q_B_I;

% Apply Dirichlet condition for melting bed

if any(L_melt)
  A(1, L_melt)  = 0;
  B(1, L_melt)  = 1;
  C(1, L_melt)  = 0;
end

if any(~L_melt)
  A(1, ~L_melt) = 0;
  B(1, ~L_melt) = 1 + 2*(dt_T/d_ZETA^2)*KAPPA_I(1,~L_melt)'./H_I(~L_melt).^2;
  C(1, ~L_melt) = -2*(dt_T/d_ZETA^2)*KAPPA_I(1,~L_melt)'./H_I(~L_melt).^2; 
end

A = reshape(A, N_I, 1);
B = reshape(B, N_I, 1);
C = reshape(C, N_I, 1);

ic_jc_kc_I = (1:N_I)';

row   = [ic_jc_kc_I; ic_jc_kc_I; ic_jc_kc_I];
col_A = max(1, ic_jc_kc_I-1);
col_B = ic_jc_kc_I;
col_C = min(N_I, ic_jc_kc_I+1);

col = [col_A; col_B; col_C];
val = [A; B; C];

M = sparse(row, col, val);

for it=1:nt_T
  
  dT_dX_I   = (UP_x.*(T_I - T_im_jc_I) + ~UP_x.*(T_ip_jc_I - T_I))/dx;
  dT_dY_I   = (UP_y.*(T_I - T_ic_jm_I) + ~UP_y.*(T_ic_jp_I - T_I))/dy;
  dT_dZ_I   = (UP_z.*(T_I(kc_T,:) - T_I(km_T,:)) + ~UP_z.*(T_I(kp_T,:) - T_I(kc_T,:)))/d_ZETA;
  
  D       = T_I - dt_T*u_I.*dT_dX_I ...
                - dt_T*v_I.*dT_dY_I ...
          + (dt_T/RHO)*dk_dT_I./(c_I.*HH_I.^2).*dT_dZ_I.^2 ...
          + (dt_T/RHO)*PHI_I./c_I;
      
  D(n_ZETA,:) = T_S_I;
  
  if any(L_melt)
    D(1, L_melt)  = T_MP_I(L_melt);     
  end
  
  % Apply Neumann condition for frozen bed
  
  % Note there are some subtle pitfalls relating to array shapes when
  % using the 2D arrays like H_I and
  % q_G_I in conjunction with 3D arrays this explains the extensive use of transpose
  % operator below
  
  if any(~L_melt)
    D(1, ~L_melt) = T_I(1,~L_melt) ...
                  + 2*(dt_T/d_ZETA)*KAPPA_I(1,~L_melt).*q_G_I(~L_melt)'./(k_I(1,~L_melt).*H_I(~L_melt)');
                  + (dt_T/(RHO*d_ZETA^2))./(c_I(1,~L_melt).*H_I(~L_melt)'.^2).*dk_dT_I(1,~L_melt).*(T_I(2,~L_melt)-T_I(1,~L_melt)).^2 ...
                  + (dt_T/RHO)*PHI_I(1,~L_melt)./c_I(1,~L_melt);
  end
  
  D = reshape(D, N_I, 1);
               
  T_out_I = M\D;
  T_lim_I = T_KELVIN-BETA*kron(ZETA, H_I);     % Limit ice temperature to pressure melting temperature
  T_lim_I = reshape(T_lim_I, N_I, 1);
  T_out_I(T_out_I>T_lim_I) = T_lim_I(T_out_I>T_lim_I);

  T_out_I = reshape(T_out_I, n_ZETA, N_xy_I);
  T_I     = T_out_I;
end

T_zeta_out.A(:,ic_jc_I) = T_out_I;
T_zeta_out.grid = 'ZETA';

[Br_max, i_max] = max(H_I'.^2.*PHI_I(1,:)./(k_I(1,:).*T_I(1,:)));
if Br_max>1
  fprintf(1,'**step_T(%.2f): High Brinkman number (Br=%.4e) with H_I(%d)=%.2f m\n', t+dt, Br_max, i_max, H_I(i_max));
  H_map = nan(N, 1);
  H_map(ic_jc_I) = H_I;
  H_map = reshape(H_map, ny, nx);
end  
Area   = numel(H_I)*dx*dy/10^6;
Vol    = sum(H_I)*dx*dy/10^9;

t_end = t+dt;

if VERBOSITY>=1
  fprintf(1,'  step_T(%.2f): min(T_I)=%.2f K; max(T_I)=%.2f K; max(T_S_I)=%.2f K; max(T_B_I)=%.2f K; max(H)=%.2f m; max(Br)=%9.3e; A(t)=%.2f km^2; V(t)=%.3f km^3\n', ...
    t_end, min(min(T_out_I)), max(max(T_out_I)), max(T_out_I(1,:)), max(T_out_I(n_ZETA, :)), max(H_I), Br_max, Area, Vol);
end  
