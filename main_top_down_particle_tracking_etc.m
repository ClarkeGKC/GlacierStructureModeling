function main_top_down_particle_tracking_etc(ModelName)

% This script is used to check the bottom-up tracking for consistency
% It uses non-reversed bottom-up-step output (i.e., from "main_bottom_up_nonreversed.m")
% as its starting point

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
elseif nargin==1
  clear global
  clearvars -except ModelName
  ModelDirName  = fullfile(pwd, sprintf('Archive-%s', ModelName));
end

global dt dx dy X2 Y2 X3 Y3 Z3 XI n_XI d_XI nx ny RHO g YRSEC

fprintf(1,'\nLAUNCHING :: main_top_down_particle_tracking_etc(''%s'') - Particle tracking from deposition points to downstream sampling sites plus additional post-processing\n', ModelName);

VERBOSE = 0;

t_start = tic;

file_date      = datestr(now);
file_source    = fullfile(pwd, mfilename);

setpref('Internet', 'E_mail', 'clarke@eos.ubc.ca');
setpref('Internet', 'SMTP_Server', 'smtp.eos.ubc.ca');

SourceFileName = fullfile(ModelDirName, 'UP_TRACKS.mat');

dir_arch     = fullfile(pwd, sprintf('Archive-%s', ModelName));
dir_SNAP     = fullfile(dir_arch, 'Snapshots');

load(SourceFileName, 'MODEL', 'SITE', 'P_up')

fprintf(1,'==================================================================================\n');
fprintf(1,'TOP-DOWN PARTICLE PATHS USING 1S2T SCHEME FOR MODEL %s\n\n', ModelName);
fprintf(1,'1. Loading data structures saved from "main_bottom_up_particle_tracking()"\n');

fprintf(1,'2. Initializing top-down tracking structure P_dn\n');

ns  = numel(P_up);

t_src  = zeros(1,ns);   t_tgt  = zeros(1, ns);
X_src  = zeros(1,ns);   X_tgt  = zeros(1, ns);
Y_src  = zeros(1,ns);   Y_tgt  = zeros(1, ns);
Z_src  = zeros(1,ns);   Z_tgt  = zeros(1, ns);
XI_src = zeros(1,ns);   XI_tgt = zeros(1,ns);

for s=1:ns
  
  t_src(s)  = P_up(s).t(end);
  X_src(s)  = P_up(s).X(end);
  Y_src(s)  = P_up(s).Y(end);
  Z_src(s)  = P_up(s).Z(end);
  XI_src(s) = P_up(s).XI(end);
  
  t_tgt(s)  = P_up(s).t(1);
  X_tgt(s)  = P_up(s).X(1);
  Y_tgt(s)  = P_up(s).Y(1);
  Z_tgt(s)  = P_up(s).Z(1);
  XI_tgt(s) = P_up(s).XI(1);
end

if any(X_tgt~=[SITE.X_site]) || any(Y_tgt~=[SITE.Y_site]) || any(Z_tgt~=[SITE.Z_init_site]) || any(XI_tgt~=[SITE.XI_init_site]) || any(t_tgt~=[SITE.t_site])
  error('main_top_down_particle_tracking_etc(): The archived P_up(s) site data do not match those in SITE structure') 
end

L_sub   = XI_tgt<1;

nx      = MODEL.nx;
ny      = MODEL.ny;
XI      = MODEL.XI;
n_XI    = MODEL.n_XI;
n_ZETA  = MODEL.n_ZETA;

d_XI    = 1/(n_XI-1);
dx      = MODEL.dx;
dy      = MODEL.dy;

dt      = MODEL.dt;

RHO     = MODEL.par.RHO;
g       = MODEL.par.g;
YRSEC   = MODEL.YRSEC;

IPREC   = GetTimeStepPrecision(dt);

load(MODEL.file_dem, 'B')

[X2, Y2] = meshgrid(MODEL.x_grid, MODEL.y_grid);
[Y3, Z3, X3] = meshgrid(MODEL.y_grid, XI, MODEL.x_grid);

load(MODEL.file_balance, 'b_dot')

if MODEL.GRID_shift
  B = B-MODEL.z_0;
end  

ns       = numel(X_src);
Site_cnt = 1:ns;
B_src    = interp2(X2, Y2, B, X_src, Y_src);

t_last_surge = SITE(1).yr_last_surge;

DIR          = dir(fullfile(dir_SNAP, '*.mat'));

ModelFileList = {DIR.name};
ModelYearList = ExtractYearListFromFileList(ModelFileList);

YR_start = min(t_src);

site_counter = 1:ns;
  
L_init       = t_src>YR_start-0.1*dt & t_src<YR_start+0.1*dt;   % Allow for slightly inaccurate timing; Defer action on L_init until end of loop
site_init    = site_counter(L_init);

fprintf(1,'3. Start top-down particle tracking\n\n');

for l=site_init
  fprintf(1,'   ... Launching site #%d at time=%.2f\n', site_counter(l), YR_start);
end  
   
[DD, V, LL, EE, ~, ~] = GetAliasYearData(YR_start, t_last_surge, dir_SNAP, ModelFileList, ModelYearList, MODEL, 0);

H            = zeros(ny,nx);
H(DD.ic_jc_I) = DD.H;
 
B_pt    = B_src(L_init);
H_pt    = interp2(X2, Y2, H, X_src(L_init), Y_src(L_init));
S_pt    = B_pt + H_pt;

L_ok    = L_init;
L_done  = false(1, ns);

X_pt    = NaN(1, ns);
Y_pt    = NaN(1, ns);
Z_pt    = NaN(1, ns);
XI_pt   = NaN(1, ns);

X_pt(L_init)  = X_src(L_init);
Y_pt(L_init)  = Y_src(L_init);
Z_pt(L_init)  = Z_src(L_init);
XI_pt(L_init) = XI_src(L_init);

% Evaluate required archived fields at points on flow trajectory

[B_pt, H_pt, S_pt, b_dot_pt, M_pt, u_pt, v_pt, w_pt, L_pt, s_pt, A_pt, ETA_pt, TAU_pt, T_pt] = EvaluateFields(L_init, X_pt, Y_pt, XI_pt, B, b_dot, DD, V, LL, EE, YR_start);

t       = 0;    

P_dn(ns)  = struct('t', [], 'X', [], 'Y', [], 'Z', [], 'XI', [], 'B', [], 'H', [], 'S', [], 'b_dot', [], 'M', [], 'u', [], 'v', [], 'w', [], ...
  'L', [], 's', [], 'A', [], 'ETA', [], 'TAU', [], 'T', []);

P_dn    = PlugValues_to_P_dn_Structure(P_dn, YR_start, X_pt, Y_pt, Z_pt, XI_pt, B_pt, H_pt, S_pt, b_dot_pt, M_pt, ...
  u_pt, v_pt, w_pt, L_pt, s_pt, A_pt, ETA_pt, TAU_pt, T_pt);

r_dist_last = Inf(1, ns);

while 1
  t       = t+1;
  YearGet = YR_start + t*dt;
  YearGet = round(YearGet,IPREC);
  
  [DD, V, LL, EE, ~, ~] = GetAliasYearData(YearGet, t_last_surge, dir_SNAP, ModelFileList, ModelYearList, MODEL, VERBOSE);
      
  H            = zeros(ny, nx);
  L_I          = false(ny, nx);
  L_I(DD.ic_jc_I) = true;
  H(L_I)       = DD.H;
    
  u            = zeros(n_XI, ny, nx);
  v            = zeros(n_XI, ny, nx);
  w            = zeros(n_XI, ny, nx);
   
  u(:,L_I)     = V.x;
  v(:,L_I)     = V.y;
  w(:,L_I)     = V.z;
   
  L_init       = t_src>YearGet-0.1*dt & t_src<YearGet+0.1*dt;   % Allow for slightly inaccurate timing; Defer action on L_init until end of loop
  
  if any(L_init)
   B_pt    = B_src(L_init);
   H_pt    = interp2(X2, Y2, H, X_src(L_init), Y_src(L_init));
   S_pt    = B_pt + H_pt;

   X_pt(L_init) = X_src(L_init);
   Y_pt(L_init) = Y_src(L_init);
   Z_pt(L_init) = S_pt;
   L_ok         = L_ok | L_init;
  end
  
  [X_p, Y_p, Z_p, XI_p, L_ok_p] = ForwardFullStep_1S2T(X_pt, Y_pt, Z_pt, B, H, b_dot, u, v, w, YearGet, L_ok);  
    
  % Introduce proximity bail-out to allow for sub-surface flow marker sites (if any)
  
  if any(L_sub)
    r_dist_p = sqrt((X_p-X_tgt).^2 + (Y_p-Y_tgt).^2 + (Z_p-Z_tgt).^2);
    L_test_p = r_dist_p<2*dx & (r_dist_p>r_dist_last) & L_sub;            % Limit the bail-out to when somewhat close to P_tgt
    if any(L_test_p)
      X_p(L_test_p)    = NaN;
      Y_p(L_test_p)    = NaN;
      Z_p(L_test_p)    = NaN;
      XI_p(L_test_p)   = NaN;
      L_ok_p(L_test_p) = 0;
    end
    r_dist_last = r_dist_p;
  end
  
  % Subtle point here: Use L_ok not L_ok_p so as to include trajectory endpoints
  
  [B_p, H_p, S_p, b_dot_p, M_p, u_p, v_p, w_p, L_p, s_p, A_p, ETA_p, TAU_p, T_p] = EvaluateFields(L_ok, X_p, Y_p, XI_p, B, b_dot, DD, V, LL, EE, YearGet);
  
  P_dn      = PlugValues_to_P_dn_Structure(P_dn, YearGet, X_p, Y_p, Z_p, XI_p, B_p, H_p, S_p, b_dot_p, M_p, ...
              u_p, v_p, w_p, L_p, s_p, A_p, ETA_p, TAU_p, T_p);
  
  % introduce logic to allow situations where no active markers but run is not completed

  L_changed = L_ok & ~L_ok_p;
  L_done    = L_done | L_changed;
  
  if ~any(L_ok_p) && all(L_done)
    break
  end
  
  X_pt = X_p;
  Y_pt = Y_p;
  Z_pt = Z_p;
  L_ok = L_ok_p;
  
  if mod(round(YearGet, IPREC), 10)==0
    fprintf(1,'   ... %d active flow markers at time=%.2f yr <%s>\n', sum(L_ok), YearGet, datestr(now));
  end
  
end

fprintf(1, '4. Correcting the shape of saved tensor time series for "L" and "s"\n');

for s=1:ns
  nt        = numel(P_dn(s).t);
  P_dn(s).L = reshape(P_dn(s).L, 3, 3, nt);
  P_dn(s).s = reshape(P_dn(s).s, 3, 3, nt);
end  

fprintf(1,'5. Top-down tracking completed\n');
fprintf(1,'6. Evaluating D tensor and its eigenproperties along tracks\n')

P_dn = Compute_D_tensor_etc(P_dn);

fprintf(1,'7. Evaluating F tensor and its polar decomposition along tracks\n')

P_dn = Compute_F_tensor_etc(P_dn);

file_mat = fullfile(dir_arch, 'DOWN_TRACKS.mat');

fprintf(1,'8. Merging SITE data with P_dn structure\n');

P_dn  = MergeSiteData(P_dn, SITE);

fprintf(1,'8. Saving results in %s\n', file_mat)

save(file_mat, 'MODEL', 'P_dn', 'file_date', 'file_source', '-v7.3')  % Use late version of save to avoid 2GB length limit

fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);
 
sendmail('clarke@eos.ubc.ca', sprintf('main_top_down_particle_tracking_etc(): Top down particle tracking run "%s" completed', ModelName), ...
  sprintf('Model %s completed at %s. Elapsed time %.2f hr.', ModelName, datestr(now), toc(t_start)/3600));
