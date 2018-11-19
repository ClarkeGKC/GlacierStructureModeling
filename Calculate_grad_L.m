function grad_L = Calculate_grad_L(MODEL, D, V_struct)

% Note that the interchangeability of partial derivatives (e.g., u_xy = u_yx
% has been demonstrated in the TeX document Velocity2ndDerivatives.tex 
% and is assumed here to reduce computing time

global X2 Y2 X3 Y3 Z3 IC_JC XI
global B B_x B_y B_xx B_yy B_xy B_yx
global ic_jc ip_jc im_jc ic_jp ic_jm im_jm im_jp ip_jm ip_jp nx ny n_XI N_xy dx dy d_XI kc km kp

if isempty(B)
  load(MODEL.file_dem, 'B')
  nx      = MODEL.nx;
  ny      = MODEL.ny;
  n_XI    = MODEL.n_XI;
  N_xy    = MODEL.N_xy;
  XI      = MODEL.XI;
  dx      = MODEL.dx;
  dy      = MODEL.dy;
  
  d_XI    = 1/(n_XI-1);

  x_grid  = MODEL.x_grid;
  y_grid  = MODEL.y_grid;

  [Y3, Z3, X3] = meshgrid(MODEL.y_grid, MODEL.XI, MODEL.x_grid);
  [X2, Y2]     = meshgrid(MODEL.x_grid, MODEL.y_grid); 
    
  SetupIndexArrays(MODEL);
  
  ic_jc   = IC_JC.ic_jc;
  im_jc   = IC_JC.im_jc;
  ip_jc   = IC_JC.ip_jc;
  ic_jm   = IC_JC.ic_jm;
  ic_jp   = IC_JC.ic_jp;
  im_jm   = IC_JC.im_jm;
  im_jp   = IC_JC.im_jp;
  ip_jm   = IC_JC.ip_jm;
  ip_jp   = IC_JC.ip_jp;
  
  B_x     = (B(ip_jc)-B(im_jc))/(2*dx);
  B_y     = (B(ic_jp)-B(ic_jm))/(2*dy);
  B_xx    = (B(ip_jc) - 2*B(ic_jc) + B(im_jc))/dx^2;
  B_yy    = (B(ic_jp) - 2*B(ic_jc) + B(ic_jm))/dy^2;
  B_xy    = (B(ip_jp) - B(im_jp) - B(ip_jm) + B(im_jm))/(4*dx*dy);
  B_yx    = B_xy;
end  

ic_jc_I   = D.ic_jc_I;

H         = zeros(ny, nx);
H(ic_jc_I) = D.H;

H_x     = (H(ip_jc)-H(im_jc))/(2*dx);
H_y     = (H(ic_jp)-H(ic_jm))/(2*dy);
H_xx    = (H(ip_jc) - 2*H(ic_jc) + H(im_jc))/dx^2;
H_yy    = (H(ic_jp) - 2*H(ic_jc) + H(ic_jm))/dy^2;
H_xy    = (H(ip_jp) - H(im_jp) - H(ip_jm) + H(im_jm))/(4*dx*dy);
H_yx    = H_xy;

H       = reshape(H, N_xy, 1);

% Define velocity components on (x,y,XI) grid

U       = zeros(n_XI, N_xy);
V       = zeros(n_XI, N_xy);
W       = zeros(n_XI, N_xy);

U(:,ic_jc_I) = V_struct.x;
V(:,ic_jc_I) = V_struct.y;
W(:,ic_jc_I) = V_struct.z;

XI_x      = -repmat(B_x./H, 1, n_XI)' - kron(XI, H_x./H)';
XI_y      = -repmat(B_y./H, 1, n_XI)' - kron(XI, H_y./H)';
XI_z      = +repmat(1./H, 1, n_XI)';

XI_xx     = +repmat(H_x.*B_x./H.^2, 1, n_XI)' + kron(XI, H_x.*H_x./H.^2)' - repmat(B_xx./H, 1, n_XI)' - XI_x.*repmat(H_x./H, 1, n_XI)' - kron(XI, H_xx./H)';
XI_xy     = +repmat(H_x.*B_y./H.^2, 1, n_XI)' + kron(XI, H_x.*H_y./H.^2)' - repmat(B_xy./H, 1, n_XI)' - XI_x.*repmat(H_y./H, 1, n_XI)' - kron(XI, H_xy./H)';
XI_xz     = -repmat(H_x./H.^2, 1, n_XI)';

XI_yx     = +repmat(H_y.*B_x./H.^2, 1, n_XI)' + kron(XI, H_y.*H_x./H.^2)' - repmat(B_yx./H, 1, n_XI)' - XI_y.*repmat(H_x./H, 1, n_XI)' - kron(XI, H_yx./H)';
XI_yy     = +repmat(H_y.*B_y./H.^2, 1, n_XI)' + kron(XI, H_y.*H_y./H.^2)' - repmat(B_yy./H, 1, n_XI)' - XI_y.*repmat(H_y./H, 1, n_XI)' - kron(XI, H_yy./H)';
XI_yz     = -repmat(H_y./H.^2, 1, n_XI)';

dXI_XI_x  = -repmat(H_x./H, 1, n_XI)';
dXI_XI_y  = -repmat(H_y./H, 1, n_XI)';

% Use subscript "k" to denote XI and compute second derivatives in the sigma "grid" (as opposed to Cartesian space derivatives)

U_xx      = (U(:,ip_jc) - 2*U(:,ic_jc) + U(:,im_jc))/dx^2;
U_yx      = (U(:,ip_jp) - U(:,im_jp) - U(:,ip_jm) + U(:,im_jm))/(4*dx*dy);

U_kx         = (U(kp,ip_jc) - U(kp,im_jc) - U(km,ip_jc) + U(km, im_jc))/(4*dx*d_XI);
U_kx(1,:)    = (U( 2,ip_jc) - U( 2,im_jc) - U( 1,ip_jc) + U( 1, im_jc))/(2*dx*d_XI);
U_kx(n_XI,:) = (U(n_XI,ip_jc) - U(n_XI,im_jc) - U(n_XI-1, ip_jc) + U(n_XI-1,im_jc))/(2*dx*d_XI);

U_xy         = U_yx;
U_yy      = (U(:,ic_jp) - 2*U(:,ic_jc) + U(:,ic_jm))/dy^2;

U_ky         = (U(kp,ic_jp) - U(kp,ic_jm) - U(km,ic_jp) + U(km,ic_jm))/(4*dy*d_XI);
U_ky(1,:)    = (U(2,ic_jp) - U(2,ic_jm) - U(1,ic_jp) + U(1,ic_jm))/(2*dy*d_XI);
U_ky(n_XI,:) = (U(n_XI,ic_jp) - U(n_XI,ic_jm) - U(n_XI-1,ic_jp) + U(n_XI-1,ic_jm))/(2*dy*d_XI);

U_xk         = U_kx; 
U_yk         = U_ky;

U_kk         = (U(kp,:) - 2*U(kc,:) + U(km,:))/d_XI^2;
U_kk(1,:)    = (U(3,:) - 2*U(2,:) + U(1,:))/d_XI^2;
U_kk(n_XI,:) = (U(n_XI,:) - 2*U(n_XI-1,:) + U(n_XI-2,:))/d_XI^2;

U_k          = (U(kp,:) - U(km,:))/(2*d_XI);
U_k(1,:)     = (U(2,:)  - U(1,:))/d_XI;
U_k(n_XI,:)  = (U(n_XI,:) - U(n_XI-1,:))/d_XI;

u_xx         = U_xx + XI_xx.*U_k + 2*XI_x.*U_xk + XI_x.^2.*U_kk;
u_yx         = U_yx + XI_yx.*U_k + XI_x.*U_yk + XI_y.*U_kx + XI_y.*XI_x.*U_kk;
u_zx         = XI_z.*U_kx + XI_z.*dXI_XI_x.*U_k + XI_z.*XI_x.*U_kk;

u_xy         = u_yx;
u_yy         = U_yy + XI_yy.*U_k + 2*XI_y.*U_ky + XI_y.^2.*U_kk;
u_zy         = XI_z.*U_ky + XI_z.*dXI_XI_y.*U_k + XI_z.*XI_y.*U_kk;

u_xz         = u_zx;
u_yz         = u_zy;
u_zz         = XI_z.^2.*U_kk;

V_xx         = (V(:,ip_jc) - 2*V(:,ic_jc) + V(:,im_jc))/dx^2;
V_yx         = (V(:,ip_jp) - V(:,im_jp) - V(:,ip_jm) + V(:,im_jm))/(4*dx*dy);

V_kx         = (V(kp,ip_jc) - V(kp,im_jc) - V(km,ip_jc) + V(km, im_jc))/(4*dx*d_XI);
V_kx(1,:)    = (V( 2,ip_jc) - V( 2,im_jc) - V( 1,ip_jc) + V( 1, im_jc))/(2*dx*d_XI);
V_kx(n_XI,:) = (V(n_XI,ip_jc) - V(n_XI,im_jc) - V(n_XI-1, ip_jc) + V(n_XI-1,im_jc))/(2*dx*d_XI);

V_xy         = V_yx;
V_yy         = (V(:,ic_jp) - 2*V(:,ic_jc) + V(:,ic_jm))/dy^2;

V_ky         = (V(kp,ic_jp) - V(kp,ic_jm) - V(km,ic_jp) + V(km,ic_jm))/(4*dy*d_XI);
V_ky(1,:)    = (V(2,ic_jp) - V(2,ic_jm) - V(1,ic_jp) + V(1,ic_jm))/(2*dy*d_XI);
V_ky(n_XI,:) = (V(n_XI,ic_jp) - V(n_XI,ic_jm) - V(n_XI-1,ic_jp) + V(n_XI-1,ic_jm))/(2*dy*d_XI);

V_xk         = V_kx;
V_yk         = V_ky;

V_kk         = (V(kp,:) - 2*V(kc,:) + V(km,:))/d_XI^2;
V_kk(1,:)    = (V(3,:) - 2*V(2,:) + V(1,:))/d_XI^2;
V_kk(n_XI,:) = (V(n_XI,:) - 2*V(n_XI-1,:) + V(n_XI-2,:))/d_XI^2;

V_k          = (V(kp,:) - V(km,:))/(2*d_XI);
V_k(1,:)     = (V(2,:)  - V(1,:))/d_XI;
V_k(n_XI,:)  = (V(n_XI,:) - V(n_XI-1,:))/d_XI;

v_xx         = V_xx + XI_xx.*V_k + 2*XI_x.*V_xk + XI_x.^2.*V_kk;
v_yx         = V_yx + XI_yx.*V_k + XI_x.*V_yk + XI_y.*V_kx + XI_y.*XI_x.*V_kk;
v_zx         = XI_z.*V_kx + XI_z.*dXI_XI_x.*V_k + XI_z.*XI_x.*V_kk;

v_xy         = v_yx;
v_yy         = V_yy + XI_yy.*V_k + 2*XI_y.*V_ky + XI_y.^2.*V_kk;
v_zy         = XI_z.*V_ky + XI_z.*dXI_XI_y.*V_k + XI_z.*XI_y.*V_kk;

v_xz         = v_zx; 
v_yz         = v_zy;
v_zz         = XI_z.^2.*V_kk;

W_xx          = (W(:,ip_jc) - 2*W(:,ic_jc) + W(:,im_jc))/dx^2;
W_yx          = (W(:,ip_jp) - W(:,im_jp) - W(:,ip_jm) + W(:,im_jm))/(4*dx*dy);

W_kx         = (W(kp,ip_jc) - W(kp,im_jc) - W(km,ip_jc) + W(km, im_jc))/(4*dx*d_XI);
W_kx(1,:)    = (W( 2,ip_jc) - W( 2,im_jc) - W( 1,ip_jc) + W( 1, im_jc))/(2*dx*d_XI);
W_kx(n_XI,:) = (W(n_XI,ip_jc) - W(n_XI,im_jc) - W(n_XI-1, ip_jc) + W(n_XI-1,im_jc))/(2*dx*d_XI);

W_xy         = W_yx;  
W_yy         = (W(:,ic_jp) - 2*W(:,ic_jc) + W(:,ic_jm))/dy^2;

W_ky         = (W(kp,ic_jp) - W(kp,ic_jm) - W(km,ic_jp) + W(km,ic_jm))/(4*dy*d_XI);
W_ky(1,:)    = (W(2,ic_jp) - W(2,ic_jm) - W(1,ic_jp) + W(1,ic_jm))/(2*dy*d_XI);
W_ky(n_XI,:) = (W(n_XI,ic_jp) - W(n_XI,ic_jm) - W(n_XI-1,ic_jp) + W(n_XI-1,ic_jm))/(2*dy*d_XI);

W_xk         = W_kx;  
W_yk         = W_ky;

W_kk         = (W(kp,:) - 2*W(kc,:) + W(km,:))/d_XI^2;
W_kk(1,:)    = (W(3,:) - 2*W(2,:) + W(1,:))/d_XI^2;
W_kk(n_XI,:) = (W(n_XI,:) - 2*W(n_XI-1,:) + W(n_XI-2,:))/d_XI^2;

W_k          = (W(kp,:) - W(km,:))/(2*d_XI);
W_k(1,:)     = (W(2,:)  - W(1,:))/d_XI;
W_k(n_XI,:)  = (W(n_XI,:) - W(n_XI-1,:))/d_XI;

w_xx         = W_xx + XI_xx.*W_k + 2*XI_x.*W_xk + XI_x.^2.*W_kk;
w_yx         = W_yx + XI_yx.*W_k + XI_x.*W_yk + XI_y.*W_kx + XI_y.*XI_x.*W_kk;
w_zx         = XI_z.*W_kx + XI_z.*dXI_XI_x.*W_k + XI_z.*XI_x.*W_kk;

w_xy         = w_yx;
w_yy         = W_yy + XI_yy.*W_k + 2*XI_y.*W_ky + XI_y.^2.*W_kk;
w_zy         = XI_z.*W_ky + XI_z.*dXI_XI_y.*W_k + XI_z.*XI_y.*W_kk;

w_xz         = w_zx;
w_yz         = w_zy;
w_zz         = XI_z.^2.*W_kk;

grad_L       = zeros(3, 3, 3, n_XI, N_xy);

grad_L(1,1,1,:,:) = u_xx;
grad_L(1,1,2,:,:) = u_xy;
grad_L(1,1,3,:,:) = u_xz;

grad_L(1,2,1,:,:) = v_xx;
grad_L(1,2,2,:,:) = v_xy;
grad_L(1,2,3,:,:) = v_xz;

grad_L(1,3,1,:,:) = w_xx;
grad_L(1,3,2,:,:) = w_xy;
grad_L(1,3,3,:,:) = w_xz;
                   
grad_L(2,1,1,:,:) = u_yx;
grad_L(2,1,2,:,:) = u_yy;
grad_L(2,1,3,:,:) = u_yz;

grad_L(2,2,1,:,:) = v_yx;
grad_L(2,2,2,:,:) = v_yy;
grad_L(2,2,3,:,:) = v_yz;

grad_L(2,3,1,:,:) = w_yx;
grad_L(2,3,2,:,:) = w_yy;
grad_L(2,3,3,:,:) = w_yz;
                   
grad_L(3,1,1,:,:) = u_zx;
grad_L(3,1,2,:,:) = u_zy;
grad_L(3,1,3,:,:) = u_zz;

grad_L(3,2,1,:,:) = v_zx;
grad_L(3,2,2,:,:) = v_zy;
grad_L(3,2,3,:,:) = v_zz;

grad_L(3,3,1,:,:) = w_zx;
grad_L(3,3,2,:,:) = w_zy;
grad_L(3,3,3,:,:) = w_zz;
                   
grad_L = reshape(grad_L, 3, 3, 3, n_XI, ny, nx);                   














