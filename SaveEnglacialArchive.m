function SaveEnglacialArchive(file_arc, yr, H, Vel, L, DD, T, E, M)

global  ZETA XI n_XI IC_JC

DD.nx      = IC_JC.nx;
DD.ny      = IC_JC.ny;
DD.n_XI    = n_XI;
DD.N       = IC_JC.N;
DD.dx      = IC_JC.dx;
DD.dy      = IC_JC.dy;
DD.N_all   = DD.N*DD.n_XI;

L_I        = H>0;
DD.ic_jc_I = IC_JC.ic_jc(L_I);

DD.H       = H(L_I);
DD.M       = M(L_I);
DD.t       = yr;

if ~isempty(DD.xd)
  DD.xd      = DD.xd(:,L_I);
  DD.yd      = DD.yd(:,L_I);
  DD.td      = DD.td(:,L_I);
else
  DD.xd      = [];
  DD.yd      = [];
  DD.td      = [];
end

N_xy       = IC_JC.N;

V.x   = Vel.u_ic_jc(:,L_I);
V.y   = Vel.v_ic_jc(:,L_I);
V.z   = Vel.w_ic_jc(:,L_I);
V.t   = yr;

T_tmp = T;                 % T on ZETA grid
T.A   = T_tmp.A;

if isempty(T.A)
  EE.T = [];
else  
  T_xi.grid = 'XI';
  T_xi.A = interp1q(ZETA, T.A, XI);
  EE.T    = T_xi.A(:,L_I);
end

EE.t     = yr;
EE.A     = E.A(:,L_I);
EE.ETA   = E.ETA(:,L_I);
EE.TAU   = E.TAU(:,L_I);
EE.s     = E.s(:,:,:,L_I);

LL.t     = yr;
LL.L     = L(:,:,:,L_I);

save(file_arc, 'DD', 'EE', 'LL', 'T', 'V')