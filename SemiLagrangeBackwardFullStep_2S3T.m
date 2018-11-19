function [X_m1, Y_m1, Z_m1, XI_m1, X_m2, Y_m2, Z_m2, XI_m2] = SemiLagrangeBackwardFullStep_2S3T(B, H, H_m1, H_m2, V_m1, t)

% This script uses the 2S3T stepping and "corrects" for any off-grid values by moving points onto the grid

global dt nx ny N_xy n_XI X2 Y2 X3 Y3 Z3
global x_MIN x_MAX y_MIN y_MAX

if isempty(x_MIN)
  x_MIN = min(min(min(X3)));
  x_MAX = max(max(max(X3)));
  y_MIN = min(min(min(Y3)));
  y_MAX = max(max(max(Y3)));
end

MAX_CNT_FLAG = 5;
MAX_CNT_EXIT = 25;

d_max   = 0.0001;
h_thin  = 0.20;

XI_vec  = linspace(0, 1, n_XI)';

B     = reshape(B, ny, nx);
H     = reshape(H, ny, nx);
H_m1  = reshape(H_m1, ny, nx);
H_m2  = reshape(H_m2, ny, nx);
u_m1  = reshape(V_m1.x, n_XI, ny, nx);
v_m1  = reshape(V_m1.y, n_XI, ny, nx);
w_m1  = reshape(V_m1.z, n_XI, ny, nx);

X  = X3;
Y  = Y3;

Z  = repmat(reshape(B, 1, N_xy), n_XI,1) + reshape(Z3, n_XI, N_xy).*repmat(reshape(H, 1, N_xy), n_XI, 1);
Z  = reshape(Z, n_XI, ny, nx);

L_I = H>0; 

B3  = repmat(reshape(B, 1, N_xy), n_XI, 1);
B3  = reshape(B3, n_XI, ny, nx);
H3  = repmat(reshape(H, 1, N_xy), n_XI, 1);
H3  = reshape(H3, n_XI, ny, nx);

Z3_I = Z3(:,L_I);

L_thin = H3(:,L_I)<h_thin;

XI_0         = (Z(:,L_I)-B3(:,L_I))./H3(:,L_I);
XI_0(L_thin) = Z3_I(L_thin);

a     = 0.5*dt*interp3(Y3, Z3, X3, u_m1, Y(:,L_I), XI_0, X(:,L_I));
b     = 0.5*dt*interp3(Y3, Z3, X3, v_m1, Y(:,L_I), XI_0, X(:,L_I));
c     = 0.5*dt*interp3(Y3, Z3, X3, w_m1, Y(:,L_I), XI_0, X(:,L_I));

cnt   = 1;

H_mid = 0.5*(H+H_m1);

% Loop processing might be speeded up if only non-converged points interated -- Not proven

while 1
  a_last = a;
  b_last = b;
  c_last = c;
  
  X_m05  = X3(:,L_I) - a_last;
  Y_m05  = Y3(:,L_I) - b_last;
  Z_m05  = Z(:,L_I)  - c_last;
  
  L_x_MIN = X_m05<x_MIN;
  L_x_MAX = X_m05>x_MAX;
  L_y_MIN = Y_m05<y_MIN;
  L_y_MAX = Y_m05>y_MAX;
  
  X_m05(L_x_MIN) = x_MIN;
  X_m05(L_x_MAX) = x_MAX;
  Y_m05(L_y_MIN) = y_MIN;
  Y_m05(L_y_MAX) = y_MAX;
  
  B_m05  = interp2(X2, Y2, B, X_m05, Y_m05);
  h_m05  = interp2(X2, Y2, H_mid, X_m05, Y_m05);
  
  L_thin = h_m05<h_thin;     % special handling of thin ice
  
  XI_m05         = (Z_m05-B_m05)./h_m05;
  XI_m05(L_thin) = Z3_I(L_thin);
  
  L_S_out = XI_m05>1;
  L_B_out = XI_m05<0;
  
  XI_m05(L_S_out) = 1;
  XI_m05(L_B_out) = 0;
  
  a     = 0.5*dt*interp3(Y3, Z3, X3, u_m1, Y_m05, XI_m05, X_m05);
  b     = 0.5*dt*interp3(Y3, Z3, X3, v_m1, Y_m05, XI_m05, X_m05);
  c     = 0.5*dt*interp3(Y3, Z3, X3, w_m1, Y_m05, XI_m05, X_m05);
  
  if any(isnan([reshape(a, 1, numel(a)) reshape(b, 1, numel(b)) reshape(c, 1, numel(c)) ]))
    error('SemiLagrangeBackwardFullStep(#1): NaN values in one or more of (a,b,c) at t=%.2f yr\n', t);
  end  
  
  if cnt==MAX_CNT_EXIT
    fprintf(1,'SemiLagrangeBackwardFullStep_1=2S3T(1): Loop count equals MAX_CNT_EXIT at t=%.2f yr. Exiting loop\n', MAX_CNT_EXIT, t);
  end  

  cnt = cnt+1;
  
  L_conv = abs(a-a_last)<d_max & abs(b-b_last)<d_max & abs(c-c_last)<d_max;
  L_test = L_conv | L_thin;
  if all(reshape(L_test, numel(L_test), 1))
    break
  end
end  
  
if cnt>MAX_CNT_FLAG
  fprintf(1,'SemiLagrangeBackwardFullStep_2S3T(1): First iteration count=%d exceeds MAX_CNT=%d at t=%.2f for d_max=%e\n', cnt, MAX_CNT, t, d_max); 
end

X_m1  = X3(:,L_I) - 2*a;
Y_m1  = Y3(:,L_I) - 2*b;
Z_m1  = Z(:,L_I)  - 2*c;

L_x_MIN = X_m1<x_MIN;
L_x_MAX = X_m1>x_MAX;
L_y_MIN = Y_m1<y_MIN;
L_y_MAX = Y_m1>y_MAX;
  
X_m1(L_x_MIN) = x_MIN;
X_m1(L_x_MAX) = x_MAX;
Y_m1(L_y_MIN) = y_MIN;
Y_m1(L_y_MAX) = y_MAX;

B_m1  = interp2(X2, Y2, B, X_m1, Y_m1);
h_m1  = interp2(X2, Y2, H_m1, X_m1, Y_m1);

L_thin = h_m1<h_thin;     % special handling of thin ice
  
XI_m1         = (Z_m1-B_m1)./h_m1;
XI_m1(L_thin) = Z3_I(L_thin);

L_S_out = XI_m1>1;
L_B_out = XI_m1<0;

XI_m1(L_S_out) = 1;
XI_m1(L_B_out) = 0;

N_I = sum(sum(L_I));

IJ = repmat(1:N_I, n_XI, 1);

Z_m1(n_XI, IJ(L_S_out)) = interp2(X2, Y2, B, X_m1(IJ(L_S_out)), Y_m1(IJ(L_S_out))) + interp2(X2, Y2, H_m1, X_m1(IJ(L_S_out)), Y_m1(IJ(L_S_out)));
Z_m1(1, IJ(L_B_out))    = interp2(X2, Y2, B, X_m1(IJ(L_B_out)), Y_m1(IJ(L_B_out)));

% First backward step is completed. Now execute the second backward step
% starting from the point arrived at in step 1

a     = 0.5*dt*interp3(Y3, Z3, X3, u_m1, Y_m1, XI_m1, X_m1);
b     = 0.5*dt*interp3(Y3, Z3, X3, v_m1, Y_m1, XI_m1, X_m1);
c     = 0.5*dt*interp3(Y3, Z3, X3, w_m1, Y_m1, XI_m1, X_m1);

cnt   = 1;

H_mid = 0.5*(H_m1+H_m2); % Grid values of H at t-(3/2)*dt

while 1
  a_last = a;
  b_last = b;
  c_last = c;
  
  X_m15  = X_m1 - a_last;
  Y_m15  = Y_m1 - b_last;
  Z_m15  = Z_m1 - c_last;
  
  L_x_MIN = X_m15<x_MIN;
  L_x_MAX = X_m15>x_MAX;
  L_y_MIN = Y_m15<y_MIN;
  L_y_MAX = Y_m15>y_MAX;
  
  X_m15(L_x_MIN) = x_MIN;
  X_m15(L_x_MAX) = x_MAX;
  Y_m15(L_y_MIN) = y_MIN;
  Y_m15(L_y_MAX) = y_MAX;
  
  B_m15  = interp2(X2, Y2, B, X_m15, Y_m15);
  h_m15  = interp2(X2, Y2, H_mid, X_m15, Y_m15);

  L_thin = h_m15<h_thin;     % special handling of thin ice
  
  XI_m15         = (Z_m15-B_m15)./h_m15;
  XI_m15(L_thin) = Z3_I(L_thin);

  L_S_out = XI_m15>1;
  L_B_out = XI_m15<0;
  
  XI_m15(L_S_out) = 1;
  XI_m15(L_B_out) = 0;
  
  a     = 0.5*dt*interp3(Y3, Z3, X3, u_m1, Y_m15, XI_m15, X_m15);
  b     = 0.5*dt*interp3(Y3, Z3, X3, v_m1, Y_m15, XI_m15, X_m15);
  c     = 0.5*dt*interp3(Y3, Z3, X3, w_m1, Y_m15, XI_m15, X_m15);
  
  if any(isnan([reshape(a, 1, numel(a)) reshape(b, 1, numel(b)) reshape(c, 1, numel(c)) ]))
    error('SemiLagrangeBackwardFullStep(#2): NaN values in one or more of (a,b,c) at t=%.2f yr\n', t);
  end  
  
  if cnt==MAX_CNT_EXIT
    fprintf(1,'SemiLagrangeBackwardFullStep_2S3T(2): Loop count equals MAX_CNT_EXIT at t=%.2f yr. Exiting loop\n', MAX_CNT_EXIT, t);
  end  

  cnt = cnt+1;
  
  L_conv = abs(a-a_last)<d_max & abs(b-b_last)<d_max & abs(c-c_last)<d_max;
  L_test = L_conv | L_thin;
  if all(reshape(L_test, numel(L_test), 1))
    break
  end
end  

if cnt>MAX_CNT_FLAG
  fprintf(1,'SemiLagrangeBackwardFullStep_2S3T(2): Second iteration count=%d exceeds MAX_CNT=%d at t=%.2f for d_max=%e\n', cnt, MAX_CNT, t, d_max);   
end

X_m2  = X_m1 - 2*a;
Y_m2  = Y_m1 - 2*b;
Z_m2  = Z_m1 - 2*c;

L_x_MIN = X_m2<x_MIN;
L_x_MAX = X_m2>x_MAX;
L_y_MIN = Y_m2<y_MIN;
L_y_MAX = Y_m2>y_MAX;
  
X_m2(L_x_MIN) = x_MIN;
X_m2(L_x_MAX) = x_MAX;
Y_m2(L_y_MIN) = y_MIN;
Y_m2(L_y_MAX) = y_MAX;

B_m2  = interp2(X2, Y2, B, X_m2, Y_m2);
h_m2  = interp2(X2, Y2, H_m2, X_m2, Y_m2);

XI_m2 = (Z_m2-B_m2)./h_m2;

L_S_out = XI_m2>1;
L_B_out = XI_m2<0;

XI_m2(L_S_out) = 1;
XI_m2(L_B_out) = 0;

Z_m2(n_XI, IJ(L_S_out)) = interp2(X2, Y2, B, X_m2(IJ(L_S_out)), Y_m2(IJ(L_S_out))) + interp2(X2, Y2, H_m2, X_m2(IJ(L_S_out)), Y_m2(IJ(L_S_out)));
Z_m2(1, IJ(L_B_out))    = interp2(X2, Y2, B, X_m2(IJ(L_B_out)), Y_m2(IJ(L_B_out)));

XX_m2   = X3;
YY_m2   = Y3;
ZZ_m2   = Z;
XIXI_m2 = Z3;

XX_m2(:,L_I)   = X_m2;
YY_m2(:,L_I)   = Y_m2;
ZZ_m2(:,L_I)   = Z_m2;
XIXI_m2(:,L_I) = XI_m2;

% Now straighten out (simplify) the notation for output variables

X_m2  = XX_m2;
Y_m2  = YY_m2;
Z_m2  = ZZ_m2;
XI_m2 = XIXI_m2;

XX_m1   = X3;
YY_m1   = Y3;
ZZ_m1   = Z;
XIXI_m1 = Z3;

XX_m1(:,L_I)   = X_m1;
YY_m1(:,L_I)   = Y_m1;
ZZ_m1(:,L_I)   = Z_m1;
XIXI_m1(:,L_I) = XI_m1;

% Now straighten out (simplify) the notation for output variables

X_m1  = XX_m1;
Y_m1  = YY_m1;
Z_m1  = ZZ_m1;
XI_m1 = XIXI_m1;