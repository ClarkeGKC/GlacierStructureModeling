function [N, idx] = SetupSubgridIndexArrays(nx, ny)

N = nx*ny;

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

idx.ic_jc  = reshape(ic_jc,  N, 1);
idx.im_jc  = reshape(im_jc,  N, 1);
idx.imm_jc = reshape(imm_jc, N, 1);
idx.ip_jc  = reshape(ip_jc,  N, 1);
idx.ipp_jc = reshape(ipp_jc, N, 1);

idx.ic_jm  = reshape(ic_jm,  N, 1);
idx.ic_jmm = reshape(ic_jmm, N, 1);
idx.ic_jp  = reshape(ic_jp,  N, 1);
idx.ic_jpp = reshape(ic_jpp, N, 1);

idx.im_jm  = reshape(im_jm, N, 1);
idx.ip_jm  = reshape(ip_jm, N, 1);
idx.im_jp  = reshape(im_jp, N, 1);
idx.ip_jp  = reshape(ip_jp, N, 1);

