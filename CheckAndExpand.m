function [DD, V, LL, EE] = CheckAndExpand(DD, V, LL, EE)

ic_jc_I = DD.ic_jc_I;

nx   = DD.nx;
ny   = DD.ny;
n_XI = DD.n_XI;
N_xy = DD.N;

H    = zeros(N_xy, 1);
H(ic_jc_I) = DD.H;
DD.H = reshape(H, ny, nx);

M = zeros(N_xy, 1);
M(ic_jc_I) = DD.M;
DD.M  = reshape(M, ny, nx);


u  = NaN(n_XI, N_xy);
v  = NaN(n_XI, N_xy);
w  = NaN(n_XI, N_xy);

u(:,ic_jc_I) = V.x;
v(:,ic_jc_I) = V.y;
w(:,ic_jc_I) = V.z;

V.x = reshape(u, n_XI, ny, nx);
V.y = reshape(v, n_XI, ny, nx);
V.z = reshape(w, n_XI, ny, nx);

if isfield(DD, 'xd')
  xd = NaN(n_XI, N_xy);
  yd = NaN(n_XI, N_xy);
  td = NaN(n_XI, N_xy);
  
  xd(:,ic_jc_I) = DD.xd;
  yd(:,ic_jc_I) = DD.yd;
  td(:,ic_jc_I) = DD.td;
  
  DD.xd = reshape(xd, n_XI, ny, nx);
  DD.yd = reshape(yd, n_XI, ny, nx);
  DD.td = reshape(td, n_XI, ny, nx);
end  

if nargin==2
  LL = [];
  EE = [];
  return
end

L                = zeros(3, 3, n_XI, ny, nx);
L(:,:,:,ic_jc_I) = LL.L;
LL.L             = L;

if nargin==3
  EE = [];
  return
end

T       = NaN(n_XI, N_xy);
A       = NaN(n_XI, N_xy);
ETA     = NaN(n_XI, N_xy);
TAU     = NaN(n_XI, N_xy);

if isempty(EE.T)
  EE = rmfield(EE, 'T');
else  
  T(:,ic_jc_I) = EE.T;
  EE.T          = reshape(T, n_XI, ny, nx);
end

A(:,ic_jc_I)   = EE.A;
ETA(:,ic_jc_I) = EE.ETA;
TAU(:,ic_jc_I) = EE.TAU;

EE.A      = reshape(A, n_XI, ny, nx);
EE.ETA    = reshape(ETA, n_XI, ny, nx);
EE.TAU    = reshape(TAU, n_XI, ny, nx);

e   = NaN(3, 3, n_XI, N_xy);
s   = NaN(3, 3, n_XI, N_xy);
W   = NaN(3, 3, n_XI, N_xy);

e(:,:,:,ic_jc_I) = EE.e;
s(:,:,:,ic_jc_I) = EE.s;
W(:,:,:,ic_jc_I) = EE.W;

EE.e = reshape(e, 3, 3, n_XI, ny, nx);
EE.s = reshape(s, 3, 3, n_XI, ny, nx);
EE.W = reshape(W, 3, 3, n_XI, ny, nx);
