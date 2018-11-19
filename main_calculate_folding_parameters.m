function main_calculate_folding_parameters(ModelName)

% This script calculates and save the potential folding parameters
% (linear and quadratic scalar invariants of the grad_F tensor)

% Note that the linear and quadratic scalar invariants of rank 3 tensors
% are derived and listed in Ahmad (2001)

% Amhmad, F. 2001. Invariants of a Cartesian tensor of rank 3
%   Arch. Mech., 63(4), 383--392.

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
  ModelDirName = fullfile(pwd, ModelDirName);
elseif nargin==1
  clear global
  clearvars -except ModelName
  ModelDirName = fullfile(pwd, sprintf('Archive-%s', ModelName));
else
  error('main_calculate_and_save_folding_parameters(): Unprogrammed number of input arguments')
end  

t_start = tic;


setpref('Internet', 'E_mail', 'clarke@eos.ubc.ca');
setpref('Internet', 'SMTP_Server', 'smtp.eos.ubc.ca');

file_date   = datestr(now);
file_source = fullfile(pwd, mfilename);

dir_ARC     = ModelDirName;
dir_FOLD    = fullfile(dir_ARC, 'FoldingParameter');
file_DAT    = fullfile(dir_FOLD, 'FOLD.mat');

DIR         = dir(file_DAT);

fprintf(1,'\n1. Loading F-tracking data (%.2f GB) from "%s"\n', DIR.bytes/10^9, file_DAT);

load(file_DAT, 'MODEL', 'F_par')

ns     = numel(F_par);
s_LIST = 1:ns;

% Prune the structure of superfluous elements

dt     = MODEL.dt;
nx     = MODEL.nx;
ny     = MODEL.ny;
n_XI   = MODEL.n_XI;
N_xy   = MODEL.N_xy;
x_grid = MODEL.x_grid;
y_grid = MODEL.y_grid;
XI     = MODEL.XI;

[Y3, Z3, X3] = meshgrid(y_grid, XI, x_grid);

% Define terms of alternating tensor

e_ijk        = zeros(3,3,3);

e_ijk(1,2,3) = 1;
e_ijk(3,1,2) = 1;
e_ijk(2,3,1) = 1;
e_ijk(3,2,1) = -1;
e_ijk(2,1,3) = -1;
e_ijk(1,3,2) = -1;

n_JJ  = 11;  

for s=1:ns
  t      = F_par(s).t;
  nt     = numel(t);
  
  F      = F_par(s).F;
  grad_F = zeros(3,3,3,nt);
  J      = zeros(1, nt);
  JJ     = zeros(n_JJ, nt);
  L      = F_par(s).L;
  grad_L = F_par(s).grad_L;
    
  for t=2:nt 
    A    = eye(3)-0.5*dt*L(:,:,t-1);
    B    = eye(3)+0.5*dt*L(:,:,t-1);
 
    grad_F(1,:,:,t) = inv(A)*B*squeeze(grad_F(1,:,:,t-1)) + dt*inv(A)*squeeze(grad_L(1,:,:,t-1))*squeeze(F(:,:,t-1));
    grad_F(2,:,:,t) = inv(A)*B*squeeze(grad_F(2,:,:,t-1)) + dt*inv(A)*squeeze(grad_L(2,:,:,t-1))*squeeze(F(:,:,t-1));
    grad_F(3,:,:,t) = inv(A)*B*squeeze(grad_F(3,:,:,t-1)) + dt*inv(A)*squeeze(grad_L(3,:,:,t-1))*squeeze(F(:,:,t-1));

    J(t)      = sum(sum(sum(e_ijk.*squeeze(grad_F(:,:,:,t)) )));
    
    grad_F_iik = squeeze(grad_F(1,1,:,t) + grad_F(2,2,:,t) + grad_F(3,3,:,t));
    grad_F_ppk = grad_F_iik;
    grad_F_ppj = grad_F_iik;
   
    JJ(1,t)     = dot(grad_F_iik, grad_F_ppk);
   
    grad_F_iji = squeeze(grad_F(1,:,1,t) + grad_F(2,:,2,t) + grad_F(3,:,3,t));
    grad_F_pjp = grad_F_iji;
    grad_F_pkp = grad_F_iji;
    grad_F_pip = grad_F_iji;
   
    JJ(2,t)     = dot(grad_F_iji, grad_F_pjp);
   
    grad_F_ijj = squeeze(grad_F(:,1,1,t) + grad_F(:,2,2,t) + grad_F(:,3,3,t));
    grad_F_iqq = grad_F_ijj;
    grad_F_kpp = grad_F_ijj;
    grad_F_jqq = grad_F_ijj;
   
    JJ(3,t)    = dot(grad_F_ijj, grad_F_iqq);
   
    grad_F_ijk = squeeze(grad_F(:,:,:,t));
   
    JJ(4,t)    = sum(sum(sum(grad_F_ijk.*grad_F_ijk)));
   
    grad_F_ikj = permute(grad_F_ijk, [1 3 2]);
    grad_F_jik = permute(grad_F_ijk, [2 1 3]);
    grad_F_kji = permute(grad_F_ijk, [3 2 1]);
    grad_F_jki = permute(grad_F_ijk, [2 3 1]);
    grad_F_kij = permute(grad_F_ijk, [3 1 2]);
   
    JJ(5,t)    = sum(sum(sum(grad_F_ijk.*grad_F_ikj)));
    JJ(6,t)    = sum(sum(sum(grad_F_ijk.*grad_F_jik)));
    JJ(7,t)    = sum(sum(sum(grad_F_ijk.*grad_F_kji)));
    JJ(8,t)    = dot(grad_F_iik, grad_F_kpp);
    JJ(9,t)    = 0.5*(dot(grad_F_iik, grad_F_pkp) + dot(grad_F_iji, grad_F_ppj));
    JJ(10,t)   = 0.5*(dot(grad_F_iji, grad_F_jqq) + dot(grad_F_ijj, grad_F_pip));
    JJ(11,t)   = 0.5*(sum(sum(sum(grad_F_ijk.*grad_F_kij + grad_F_ijk.*grad_F_jki))));       
  end
  F_par(s).grad_F = grad_F;
  F_par(s).J      = J;
  F_par(s).JJ     = JJ;
end
 
save(fullfile(dir_FOLD, 'FOLD_PAR.mat'), 'MODEL', 'F_par', 'file_date', 'file_source', '-v7.3')

fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);

sendmail('clarke@eos.ubc.ca', sprintf('main_calculate_folding_parameters(): Structure evolution run "%s" completed', ModelName), ...
  sprintf('Model %s completed at %s. Elapsed time %.2f hr.', ModelName, datestr(now), toc(t_start)/3600));






