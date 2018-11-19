function main_calculate_and_save_F_par_structure(ModelName)

% This script plots the folding parameter

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
  error('main_calculate_and_save_grad_F_par_structure(): Unprogrammed number of input arguments')
end  

fprintf(1,'\nLAUNCHING :: main_calculate_and_save_F_par_structure(''%s'') - Calculate grad_F tensor and grad_F scalar invariants along particle tracks (for use in folding parameterization)\n', ModelName);

t_start = tic;

file_date   = datestr(now);
file_source = fullfile(pwd, mfilename);

setpref('Internet', 'E_mail', 'clarke@eos.ubc.ca');
setpref('Internet', 'SMTP_Server', 'smtp.eos.ubc.ca');

dir_ARC     = ModelDirName;
dir_GRAD    = fullfile(dir_ARC, 'grad_L');
file_DAT    = fullfile(dir_ARC, 'DOWN_TRACKS.mat');
dir_FOLD    = fullfile(dir_ARC, 'FoldingParameter');

if ~exist(dir_FOLD)
  mkdir(dir_FOLD)
end  

DIR         = dir(file_DAT);

fprintf(1,'\n1. Loading down-tracking data (%.2f GB) from "%s"\n', DIR.bytes/10^9, file_DAT);

load(file_DAT, 'MODEL', 'P_dn')

F_par  = P_dn;

ns     = numel(P_dn);
s_LIST = 1:ns;

% Prune the structure of superfluous elements

F_par  = rmfield(F_par, {'b_dot', 'M', 's', 'A', 'ETA', 'TAU', 'T', 'sigma_zz', 'p', 'sigma', 'D', 'W', 'E_D', 'vecs_D', ...
                         'E_s', 'E_sigma', 'E_R', 'R_pre', 'R_post', 'U', 'V', 'E_U', 'vecs_U', 'R_vecs_U'});

for s=s_LIST                       
  F_par(s).grad_L = [];                       
end

dt     = MODEL.dt;
nx     = MODEL.nx;
ny     = MODEL.ny;
n_XI   = MODEL.n_XI;
N_xy   = MODEL.N_xy;
x_grid = MODEL.x_grid;
y_grid = MODEL.y_grid;
XI     = MODEL.XI;

[Y3, Z3, X3] = meshgrid(y_grid, XI, x_grid);
[X2, Y2]     = meshgrid(x_grid, y_grid); 

t_MIN  = NaN(1,ns);
t_MAX  = NaN(1,ns);

for s=1:ns
  t_MIN(s) = P_dn(s).t(1);
  t_MAX(s) = P_dn(s).t(end);
end

t_LIST        = min(t_MIN):dt:max(t_MAX);
yr_last_surge = F_par(1).yr_last_surge;

DIR           = dir(fullfile(dir_GRAD, '*.mat'));
ModelFileList = {DIR.name};
ModelYearList = ExtractYearListFromFileList(ModelFileList);

for t=1:numel(t_LIST)  
  yr_TGT   = t_LIST(t); 
  
  if mod(yr_TGT, 100)==0
    fprintf(1,'%.0f\n', yr_TGT);
  elseif mod(yr_TGT,50)==0
    fprintf(1,'50');
  elseif mod(yr_TGT,10)==0
    fprintf(1,'|');    
  elseif mod(yr_TGT, 1)==0
    fprintf(1,'.');
  end
  if t==numel(t_LIST)
    fprintf(1,'\n');
  end 
  [FileSwap, YearSwap] = GetYearSwapFileName(MODEL, yr_last_surge, ModelFileList, ModelYearList, yr_TGT);
  load(fullfile(dir_GRAD, FileSwap), 'grad_L')
  
  L_alive  = t_MIN<=yr_TGT & t_MAX>=yr_TGT;
  s_alive  = s_LIST(L_alive);  
  F_use    = F_par(L_alive);
      
  np       = numel(F_use);
  
  i_use    = NaN(1,np);
  
  X_p      = NaN(1,np);
  Y_p      = NaN(1,np);
  XI_p     = NaN(1,np);
  
  for f=1:np
    i_use(f) = find(abs(F_use(f).t-yr_TGT)<0.5*dt);
    if isempty(i_use(f))
      error('TILT: yr_TGT=%.3f yr', yr_TGT)
    else
      X_p(1,f) = F_use(f).X(i_use(f));
      Y_p(1,f) = F_use(f).Y(i_use(f));
      XI_p(1,f) = F_use(f).XI(i_use(f));
    end
  end
  
  grad_L_p = zeros(3,3,3,np); 
  
  % Here we exploit the proven symmetry of the cross-derivatives (see GKCC
  % TeX document  Velocity2ndDerivatives.tex)
  
% 2nd-derivatives of u

  grad_L_p(1,1,1,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,1,1,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(1,1,2,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,1,2,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(1,1,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,1,3,:,:,:)), Y_p, XI_p, X_p);
  
  grad_L_p(2,1,1,:) = grad_L_p(1,1,2,:,:);
  grad_L_p(2,1,2,:) = interp3(Y3, Z3, X3, squeeze(grad_L(2,1,2,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(2,1,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(2,1,3,:,:,:)), Y_p, XI_p, X_p);
  
  grad_L_p(3,1,1,:) = grad_L_p(1,1,3,:,:);
  grad_L_p(3,1,2,:) = grad_L_p(2,1,3,:,:);
  grad_L_p(3,1,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(3,1,3,:,:,:)), Y_p, XI_p, X_p);
  
  % 2nd derivatives of v
  
  grad_L_p(1,2,1,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,2,1,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(1,2,2,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,2,2,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(1,2,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,2,3,:,:,:)), Y_p, XI_p, X_p);

  grad_L_p(2,2,1,:) = grad_L_p(1,2,2,:,:);
  grad_L_p(2,2,2,:) = interp3(Y3, Z3, X3, squeeze(grad_L(2,2,2,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(2,2,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(2,2,3,:,:,:)), Y_p, XI_p, X_p);

  grad_L_p(3,2,1,:) = grad_L_p(1,2,3,:,:);
  grad_L_p(3,2,2,:) = grad_L_p(2,2,3,:,:);
  grad_L_p(3,2,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(3,2,3,:,:,:)), Y_p, XI_p, X_p);
  
 % 2nd derivatives of w
 
  grad_L_p(1,3,1,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,3,1,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(1,3,2,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,3,2,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(1,3,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(1,3,3,:,:,:)), Y_p, XI_p, X_p);

  grad_L_p(2,3,1,:) = grad_L_p(1,3,2,:,:);
  grad_L_p(2,3,2,:) = interp3(Y3, Z3, X3, squeeze(grad_L(2,3,2,:,:,:)), Y_p, XI_p, X_p);
  grad_L_p(2,3,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(2,3,3,:,:,:)), Y_p, XI_p, X_p);
  
  grad_L_p(3,3,1,:) = grad_L_p(1,3,3,:,:);
  grad_L_p(3,3,2,:) = grad_L_p(2,3,3,:,:);
  grad_L_p(3,3,3,:) = interp3(Y3, Z3, X3, squeeze(grad_L(3,3,3,:,:,:)), Y_p, XI_p, X_p);
 
  for p=1:np
     F_par(s_alive(p)).grad_L(:,:,:,i_use(p)) = grad_L_p(:,:,:,p); 
  end 

end  
 
save(fullfile(dir_FOLD, 'FOLD.mat'), 'MODEL', 'F_par', 'file_date', 'file_source', '-v7.3')

fprintf(1,'\nALL DONE. Elapsed time %.2f hr\n', toc(t_start)/3600);

sendmail('clarke@eos.ubc.ca', sprintf('main_calculate_and_save_F_par_structure(): Structure evolution run "%s" completed', ModelName), ...
  sprintf('Model %s completed at %s. Elapsed time %.2f hr.', ModelName, datestr(now), toc(t_start)/3600));







