function [S_new, t_n, LAMBDA_max, T, V, L, DD, E] = step(S_old, B, b_dot, C_t, T, T_S, q_G, dt, t)

% Definition of staggered-grid indices (noting that dx=dy)

% IC = ic-0.5
% IP = ic+0.5
% JC = jc-0.5
% JP = jc+0.5

global  OMEGA
global  PatchStruct
global  THERMAL TRACER TRACER_METHOD MB_CHECK t_tracer_START VERBOSITY FLAG_tracer
global  IC_JC
global  X2 Y2 X3 Y3 Z3
global  ZETA n_XI XI_vec
global  B3
global  A_last V_last Bal_last t_err_plot H_err_plot N_edge_plot

global  H_m1 H_m2 V_m1
global  BOOT xd_m1 xd_m2 yd_m1 yd_m2 td_m1 td_m2
global  F_m1 L_m1
global  ji_cnt i_val j_val k_val

h_thin  = 0.10;

ic_jc   = IC_JC.ic_jc;
ip_jc   = IC_JC.ip_jc;
im_jc   = IC_JC.im_jc;
ic_jp   = IC_JC.ic_jp;
ic_jm   = IC_JC.ic_jm;
nx      = IC_JC.nx;
ny      = IC_JC.ny;
N_xy    = IC_JC.N;

dx      = IC_JC.dx;
dy      = IC_JC.dy;

if ~exist('B3', 'var') || isempty(B3)
  B3 = repmat(B', n_XI, 1);
  B3 = reshape(B3, n_XI, ny, nx);
end  

if THERMAL
  if isempty(XI_vec)
    XI_vec = linspace(0, 1, n_XI)';
  end  
  T_zeta = T;

  if strcmp(T_zeta.grid, 'ZETA')==0
    error('step(): Input temperature field must be on ZETA grid not "%s"', T.grid)
  end

  % now convert temperature to XI grid for calculation of mechanical diffusion terms

  T_xi.grid = 'XI';
  T_xi.A    = interp1q(ZETA, T_zeta.A, XI_vec);
end

timer_start_implicit = tic;

if THERMAL
  [D_IC_jc, D_IP_jc, D_ic_JC, D_ic_JP, V, L, PHI_ic_jc.A, E, H] = Diffusion(S_old, B, C_t, T_xi.A, b_dot, t);
else
  [D_IC_jc, D_IP_jc, D_ic_JC, D_ic_JP, V, L, PHI_ic_jc.A, E, H] = Diffusion(S_old, B, C_t, [], b_dot, t);
end

V.grid    = 'XI';

u_ic_jc.A = V.u_ic_jc;
v_ic_jc.A = V.v_ic_jc;
w_ic_jc.A = V.w_ic_jc;

u_ic_jc.grid = V.grid;
v_ic_jc.grid = V.grid;
w_ic_jc.grid = V.grid;

PHI_ic_jc.grid = 'XI';

D_sum    = D_IC_jc+D_IP_jc+D_ic_JC+D_ic_JP;

if any(isinf(D_sum))
  L_inf = false(numel(D_sum), 1);
  L_inf(isnan(D_sum)) = true;
  ic_jc_inf = ic_jc(L_inf);
end

D_max    = max([D_IC_jc' D_IP_jc' D_ic_JC' D_ic_JP']);

LAMBDA_max = dt*D_max;            % Monitor the diffusion parameter even when not used

row = [ic_jc;    ic_jc;    ic_jc;    ic_jc;    ic_jc];
col = [im_jc;    ip_jc;    ic_jm;    ic_jp;    ic_jc];
val = [-OMEGA*D_IC_jc; -OMEGA*D_IP_jc; -OMEGA*D_ic_JC; -OMEGA*D_ic_JP; ...
      1/dt+OMEGA*D_sum];

A   = sparse(row,col,val,N_xy,N_xy);         % matrix A is symmetric positive definite
C   = (1-OMEGA)*(D_IC_jc.*S_old(im_jc) + D_IP_jc.*S_old(ip_jc) + D_ic_JC.*S_old(ic_jm) + D_ic_JP.*S_old(ic_jp)) ...
    + (1/dt -(1-OMEGA)*D_sum).*S_old(ic_jc) + b_dot;
  
A_cond = condest(A);

if A_cond>10^6
  disp('Condition WARNING')
end

S_new  = A\C;  
H_star = S_new-B;

S_new  = max(S_new, B);
B_dot  = dt*b_dot;

% search for cells where H_star<0
   
L_star = false(N_xy,1);

H_test = H_star-B_dot;

L_star(H_star<0 & H_test<0) = true;

% If problem cell has no ice-covered neighbour cells at lower elevation then it cannot
% be a mass fountain and any mass balance error is innocent and thus ignored

% Additionally, I impose the requirement that the cell

ic_jc_err = ic_jc(L_star);

D_eps     = 0;

L_out     = D_IC_jc(L_star)>D_eps & S_new(im_jc(L_star))<S_new(ic_jc(L_star)) | ...
            D_IP_jc(L_star)>D_eps & S_new(ip_jc(L_star))<S_new(ic_jc(L_star)) | ...
            D_ic_JC(L_star)>D_eps & S_new(ic_jm(L_star))<S_new(ic_jc(L_star)) | ...
            D_ic_JP(L_star)>D_eps & S_new(ic_jp(L_star))<S_new(ic_jc(L_star));
        
L_inp    =  D_IC_jc(L_star)>D_eps & S_new(im_jc(L_star))>S_new(ic_jc(L_star)) | ...
            D_IP_jc(L_star)>D_eps & S_new(ip_jc(L_star))>S_new(ic_jc(L_star)) | ...
            D_ic_JC(L_star)>D_eps & S_new(ic_jm(L_star))>S_new(ic_jc(L_star)) | ...
            D_ic_JP(L_star)>D_eps & S_new(ic_jp(L_star))>S_new(ic_jc(L_star));
          
ic_jc_err = ic_jc_err(L_inp & L_out);

% special treatment of points that have been denuded, entirely or partially, by flow rather than by ablation

PatchStruct.t_implicit = toc(timer_start_implicit)+PatchStruct.t_implicit;

timer_start_explicit = tic;

if ~isempty(ic_jc_err)
  S_implicit = S_new;
  PatchStruct.e_cnt = numel(ic_jc_err) + PatchStruct.e_cnt;
  PatchStruct.calls = PatchStruct.calls+1;
  [S_new, F] = PatchSolverPackage(ic_jc_err, S_old, S_new, B, b_dot, t, dt);
  PatchStruct.t_explicit = toc(timer_start_explicit) + PatchStruct.t_explicit;
  
  dS = reshape(S_new-S_implicit, ny, nx);
  
  figure(901)
  
  imagesc(dS), colorbar
  hold on
  
  for k=1:numel(F)
    if F(k).cnt>0
      plot(F(k).ic, F(k).jc, 'w+')
      text(F(k).ic, F(k).jc, num2str(k));
    end
  end
  hold off
  title(sprintf('Patch grid at t=%.2f', t))
end

if THERMAL
  H_t = (S_new - S_old)/dt;
  [T_zeta, t_end] = step_T(T_zeta, B, S_new-B, XI_to_ZETA(u_ic_jc), XI_to_ZETA(v_ic_jc), XI_to_ZETA(w_ic_jc), XI_to_ZETA(PHI_ic_jc), dt, T_S, q_G, H_t, t);
  T = T_zeta;
end

t_n = t+dt;

if MB_CHECK && isempty(V_last)
  A_last   = NaN;
  V_last   = NaN;
  Bal_last = NaN;
end

if MB_CHECK
  H_new    = S_new-B;
  A_new    = dx*dy*sum(sum(H_new>0));
  V_new    = dx*dy*sum(sum(H_new));
  Bal_new  = dx*dy*sum(max(dt*b_dot, -H_new));
  A_mid    = 0.5*(A_last+A_new); 
  V_mid    = 0.5*(V_last+V_new);
  dV       = V_new-V_last;
  A_last   = A_new;
  V_last   = V_new;
  dV_bal   = 0.5*(Bal_last+Bal_new);
  Bal_last = Bal_new;
  dV_err   = dV-dV_bal;
  dH_err   = dV_err/A_mid;
  
  L_I      = reshape(H_new>0, ny, nx);
  N_edge   = sum(L_I(1,:)) + sum(L_I(end,:)) + sum(L_I(:,1))+sum(L_I(:,end));
  if ~isnan(dV_err)
    fprintf(1,'step(%d): V_mid=%.0f; V_dV=%.0f; dV_bal=%.0f; dV_err=%.0f; dH_err=%.6f\n', N_edge, V_mid, dV, dV_bal, dV_err, dH_err)
  end  
  t_err_plot(end+1) = t;
  H_err_plot(end+1) = dH_err;
  N_edge_plot(end+1) = N_edge;
end

% Bail out at this point if t_n less than tracer start time t_tracer_START

if t_n<t_tracer_START
  F      = [];
  DD.xd  = [];
  DD.yd  = [];
  DD.td  = [];
  return
end

% -------------------------------------------------------------------------------
% IF NO BAIL-OUT THEN TRACER TRANSPORT HAS BEEN ACTIVATED FOR ALL SUBSEQUENT CODE
%--------------------------------------------------------------------------------

if FLAG_tracer
  if isempty(BOOT) 
    fprintf(1,'\nLAUNCHING TRACER MARKERS AT %.2f <%s>\n\n', t_n, datestr(now));
  else
    fprintf(1,'\nRELAUNCHING TRACER MARKERS AT %.2f <%s>\n\n', t_n, datestr(now));
    switch TRACER_METHOD
      case '1S2T'
        xd_m1 = BOOT.xd;
        yd_m1 = BOOT.yd;
        td_m1 = BOOT.td;    
      case '2S3T'
        xd_m1 = BOOT.xd;
        yd_m1 = BOOT.yd;
        td_m1 = BOOT.td;
    
        xd_m2 = xd_m1;
        yd_m2 = yd_m1;
        td_m2 = td_m1;
      otherwise
        error('step(): Unprogrammed TRACER_METHOD "%s"', TRACER_METHOD)
    end  
  end
  FLAG_tracer = 0;
end

H = S_new-B;

u = V.u_ic_jc;
v = V.v_ic_jc;
w = V.w_ic_jc;

if isempty(H_m1)
  H_m1    = H;
  H_m2    = H;
  V_m1.x  = u;
  V_m1.y  = v;
  V_m1.z  = w;
end

L_H    = H>0;
L_acc  = b_dot>0;

switch TRACER_METHOD
  case '1S2T' 
    [X_m1, Y_m1, Z_m1, XI_m1] = SemiLagrangeBackwardFullStep_1S2T(B, H, H_m1, V_m1, t);
  case '2S3T'
    [X_m1, Y_m1, Z_m1, XI_m1, X_m2, Y_m2, Z_m2, XI_m2]  = SemiLagrangeBackwardFullStep_2S3T(B, H, H_m1, H_m2, V_m1, t); 
  otherwise
    error('step(): Unprogrammed TRACER_METHOD "%s"', TRACER_METHOD)
end      
  
% Save the following for the next call to step()

H_m2   = H_m1;
H_m1   = H;

V_m1.x = u;
V_m1.y = v;
V_m1.z = w;

% It is possible that this can be accelerated by only interpolating the
% glacier-covered points -- not proven

if TRACER
  switch TRACER_METHOD
    case '1S2T'
      if isempty(xd_m1)
        xd_m1 = X3;
        yd_m1 = Y3;
        td_m1 = ones(size(X3))*(t-dt);    
      end
  
      xd = interp3(Y3, Z3, X3, xd_m1, Y_m1, XI_m1, X_m1);
      yd = interp3(Y3, Z3, X3, yd_m1, Y_m1, XI_m1, X_m1);
      td = interp3(Y3, Z3, X3, td_m1, Y_m1, XI_m1, X_m1);
  
      % Accumulated surface mass has local provenance
  
      xd(n_XI,L_acc) = X3(n_XI,L_acc);
      yd(n_XI,L_acc) = Y3(n_XI,L_acc);
      td(n_XI,L_acc) = t;           
  
      DD.xd = xd;
      DD.yd = yd;
      DD.td = td;
  
      % save these arrays for next call to step()
  
      xd_m1 = xd;  
      yd_m1 = yd;
      td_m1 = td;    
      
    case '2S3T'
      if isempty(xd_m1)
        xd_m1 = X3;
        yd_m1 = Y3;
        td_m1 = ones(size(X3))*(t-dt);
    
        xd_m2 = X3;
        yd_m2 = Y3;
        td_m2 = ones(size(X3))*(t-2*dt);
      end
  
      xd = interp3(Y3, Z3, X3, xd_m2, Y_m2, XI_m2, X_m2);
      yd = interp3(Y3, Z3, X3, yd_m2, Y_m2, XI_m2, X_m2);
      td = interp3(Y3, Z3, X3, td_m2, Y_m2, XI_m2, X_m2);
  
      % Accumulated surface mass has local provenance
  
      xd(n_XI,L_acc) = X3(n_XI,L_acc);
      yd(n_XI,L_acc) = Y3(n_XI,L_acc);
      td(n_XI,L_acc) = t;           
  
      DD.xd = xd;
      DD.yd = yd;
      DD.td = td;
  
      % save these arrays for next call to step()
  
      xd_m2 = xd_m1;
      xd_m1 = xd;
  
      yd_m2 = yd_m1;
      yd_m1 = yd;
  
      td_m2 = td_m1;
      td_m1 = td;
    otherwise
      error('step(): Unprogrammed TRACER_METHOD "%s"', TRACER_METHOD)
  end    
end
