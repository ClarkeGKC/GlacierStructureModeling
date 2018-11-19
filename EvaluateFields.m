function [B_pt, H_pt, S_pt, b_dot_pt, M_pt, u_pt, v_pt, w_pt, L_pt, s_pt, A_pt, ETA_pt, TAU_pt, T_pt] = ...
  EvaluateFields(L_ok, X_pt, Y_pt, XI_pt, B, b_dot, DD, V, LL, EE, Year)

global X2 Y2 X3 Y3 Z3 n_XI nx ny

H             = zeros(ny, nx);
H(DD.ic_jc_I) = DD.H;
L_I           = H>0;

M             = zeros(ny, nx);
M(L_I)        = DD.M;

B_pt(L_ok)  = interp2(X2, Y2, B, X_pt(L_ok), Y_pt(L_ok));
H_pt(L_ok)  = interp2(X2, Y2, H, X_pt(L_ok), Y_pt(L_ok));
S_pt(L_ok)  = B_pt(L_ok) + H_pt(L_ok);

b_dot_pt(L_ok) = interp2(X2, Y2, b_dot, X_pt(L_ok), Y_pt(L_ok));

M_pt(L_ok)   = interp2(X2, Y2, M, X_pt(L_ok), Y_pt(L_ok));

u            = zeros(n_XI, ny, nx);
v            = zeros(n_XI, ny, nx);
w            = zeros(n_XI, ny, nx);
   
u(:,L_I)     = V.x;
v(:,L_I)     = V.y;
w(:,L_I)     = V.z;
  
L            = zeros(3, 3, n_XI, ny, nx);
L(:,:,:,L_I) = LL.L;
   
%T           = NaN(n_XI, ny, nx);
T            = zeros(n_XI, ny, nx);
T(:,L_I)     = EE.T; 
  
%s            = NaN(3,3, n_XI, ny, nx);
s            = zeros(3 ,3, n_XI, ny, nx);
s(:,:,:,L_I) = EE.s;
   
%A            = NaN(n_XI, ny, nx);
%ETA          = NaN(n_XI, ny, nx);
%TAU          = NaN(n_XI, ny, nx);
   
A            = zeros(n_XI, ny, nx);
ETA          = zeros(n_XI, ny, nx);
TAU          = zeros(n_XI, ny, nx);

A_TEST       = NaN(n_XI, ny, nx);
A_TEST(:,L_I)= EE.A;

A(:,L_I)     = EE.A;
ETA(:,L_I)   = EE.ETA;
TAU(:,L_I)   = EE.TAU;

u_pt(L_ok)  = interp3(Y3, Z3, X3, u, Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
v_pt(L_ok)  = interp3(Y3, Z3, X3, v, Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
w_pt(L_ok)  = interp3(Y3, Z3, X3, w, Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

L_pt(1,1,L_ok) = interp3(Y3, Z3, X3, squeeze(L(1,1,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
L_pt(1,2,L_ok) = interp3(Y3, Z3, X3, squeeze(L(1,2,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
L_pt(1,3,L_ok) = interp3(Y3, Z3, X3, squeeze(L(1,3,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

L_pt(2,1,L_ok) = interp3(Y3, Z3, X3, squeeze(L(2,1,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
L_pt(2,2,L_ok) = interp3(Y3, Z3, X3, squeeze(L(2,2,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
L_pt(2,3,L_ok) = interp3(Y3, Z3, X3, squeeze(L(2,3,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

L_pt(3,1,L_ok) = interp3(Y3, Z3, X3, squeeze(L(3,1,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
L_pt(3,2,L_ok) = interp3(Y3, Z3, X3, squeeze(L(3,2,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
L_pt(3,3,L_ok) = interp3(Y3, Z3, X3, squeeze(L(3,3,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

s_pt(1,1,L_ok) = interp3(Y3, Z3, X3, squeeze(s(1,1,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
s_pt(1,2,L_ok) = interp3(Y3, Z3, X3, squeeze(s(1,2,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
s_pt(1,3,L_ok) = interp3(Y3, Z3, X3, squeeze(s(1,3,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

s_pt(2,1,L_ok) = interp3(Y3, Z3, X3, squeeze(s(2,1,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
s_pt(2,2,L_ok) = interp3(Y3, Z3, X3, squeeze(s(2,2,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
s_pt(2,3,L_ok) = interp3(Y3, Z3, X3, squeeze(s(2,3,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

s_pt(3,1,L_ok) = interp3(Y3, Z3, X3, squeeze(s(3,1,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
s_pt(3,2,L_ok) = interp3(Y3, Z3, X3, squeeze(s(3,2,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
s_pt(3,3,L_ok) = interp3(Y3, Z3, X3, squeeze(s(3,3,:,:,:)), Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

A_TEST_pt      = interp3(Y3, Z3, X3, A_TEST, Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

if any(isnan(A_TEST_pt))
  fprintf(1, 'EvaluateFields(): Warning at t=%.3f yr -- Interpolation of one or more near-boundary points uses zero values for off-glacier points to avoid NaN problems\n', Year);
end  

A_pt(L_ok)     = interp3(Y3, Z3, X3, A, Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
ETA_pt(L_ok)   = interp3(Y3, Z3, X3, ETA, Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
TAU_pt(L_ok)   = interp3(Y3, Z3, X3, TAU, Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));
T_pt(L_ok)     = interp3(Y3, Z3, X3, T, Y_pt(L_ok), XI_pt(L_ok), X_pt(L_ok));

