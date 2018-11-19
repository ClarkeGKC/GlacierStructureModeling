function SetupIndexArrays(MODEL)

global IC_JC  km kc kp kmm kpp km_T kc_T kp_T

nx = MODEL.nx;
ny = MODEL.ny;
dx = MODEL.dx;
dy = MODEL.dy;
N  = MODEL.N_xy;

THERMAL = MODEL.THERMAL;

% Note that (j,i) are (r,c) indices such that

%   (ic,jc) <=> (r, c)
%   (ip,jc) <=> (r,c+1)
%   (im_jc) <=> (r,c-1)
%   (ic,jp) <=> (r-1,c)
%   (ic,jm) <=> (r+1,c)

ic_jc  = 1:N;
ic_jc  = reshape(ic_jc, ny, nx);

ic     = 1:nx;
ip     = [2:nx nx];
ipp    = [3:nx nx nx];

im     = [1 1:nx-1];
imm    = [1 1 1:nx-2];

jc     = 1:ny;
jm     = [2:ny ny];
jmm    = [3:ny ny ny];

jp     = [1 1:ny-1];
jpp    = [1 1 1:ny-2];

im_jc  = ic_jc(jc,  im);
imm_jc = ic_jc(jc,  imm);
ip_jc  = ic_jc(jc,  ip);
ipp_jc = ic_jc(jc,  ipp);

ic_jm  = ic_jc(jm,  ic);
ic_jmm = ic_jc(jmm, ic);
ic_jp  = ic_jc(jp,  ic);
ic_jpp = ic_jc(jpp, ic);

im_jm  = ic_jc(jm, im);
ip_jm  = ic_jc(jm, ip);
im_jp  = ic_jc(jp, im);
ip_jp  = ic_jc(jp, ip);

ic_jc  = reshape(ic_jc,  N, 1);
im_jc  = reshape(im_jc,  N, 1);
imm_jc = reshape(imm_jc, N, 1);
ip_jc  = reshape(ip_jc,  N, 1);
ipp_jc = reshape(ipp_jc, N, 1);

ic_jm  = reshape(ic_jm,  N, 1);
ic_jmm = reshape(ic_jmm, N, 1);
ic_jp  = reshape(ic_jp,  N, 1);
ic_jpp = reshape(ic_jpp, N, 1);

im_jm  = reshape(im_jm, N, 1);
ip_jm  = reshape(ip_jm, N, 1);
im_jp  = reshape(im_jp, N, 1);
ip_jp  = reshape(ip_jp, N, 1);

IC_JC.ic     = ic;
IC_JC.jc     = jc;
IC_JC.ic_jc  = ic_jc;
IC_JC.ip_jc  = ip_jc;
IC_JC.ipp_jc = ipp_jc;
IC_JC.im_jc  = im_jc;
IC_JC.imm_jc = imm_jc;
IC_JC.ic_jp  = ic_jp;
IC_JC.ic_jpp = ic_jpp;
IC_JC.ic_jm  = ic_jm;
IC_JC.ic_jmm = ic_jmm;
IC_JC.im_jm  = im_jm;
IC_JC.ip_jm  = ip_jm;
IC_JC.im_jp  = im_jp;
IC_JC.ip_jp  = ip_jp;
IC_JC.nx     = nx;
IC_JC.ny     = ny;
IC_JC.N      = N;
IC_JC.dx     = dx;
IC_JC.dy     = dy;

IC_JC.ic_jc_2d = reshape(ic_jc, ny, nx);

n_XI = MODEL.n_XI;

kc   = 1:n_XI;
kp   = min(kc+1, n_XI);
km   = max(kc-1, 1);

kpp  = min(kc+2, n_XI);
kmm  = max(kc-2, 1);

if THERMAL
  n_ZETA = MODEL.n_ZETA;
  kc_T = 1:n_ZETA;
  kp_T = min(kc_T+1, n_ZETA);
  km_T = max(kc_T-1, 1);
end  






