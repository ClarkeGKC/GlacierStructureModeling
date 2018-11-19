function main(ModelName, BOOTSTRAP, BootYear)

% This script runs the thermo-mechanical tracer-tracking model with add-on
% features such as the orientation matrix

close all
clear global

if nargin==3
  clearvars -except ModelName BOOTSTRAP BootYear
else
  clearvars -except ModelName
  BOOTSTRAP = 0;
end

global  n_GLEN A_GLEN RHO g dx dy m_SLIDE C_SLIDE C_SURGE ETA_0 YRSEC OMEGA nx ny N_xy M
global  A_const C_const
global  THERMAL ZETA d_ZETA n_ZETA n_XI XI XI_min XI_max
global  t_tracer_START TRACER TRACER_METHOD MB_CHECK VERBOSITY
global  n_sub_steps 
global  X3 Y3 Z3 X2 Y2 xc_min xc_max yc_min yc_max dt
global  FLAG FLAG_tracer BOOT

t_start     = tic;

fprintf(1,'\nLAUNCHING :: main(''%s'') - Thermodynamics flow model\n', ModelName);

file_model  = strcat(ModelName, '.dat');
file_source = fullfile(pwd, mfilename);
file_date   = datestr(now);


if BOOTSTRAP
  MODEL_check = GetModelData(file_model);
  load(fullfile(pwd, sprintf('Archive-%s', ModelName), 'Model.mat'), 'MODEL')
  NAMES_check = fieldnames(MODEL_check);
  NAMES_all = fieldnames(MODEL);
  
  if strcmp(MODEL.file_balance, MODEL_check.file_balance)==0
    fprintf(1, 'main(): Swapping the mass balance file to "%s"\n', MODEL_check.file_balance)
    MODEL.file_balance = MODEL_check.file_balance;
  end
  
  if ~strcmp(MODEL.name, ModelName)
    fprintf(1,'WARNING: Model name in loaded MODEL structure differs from that of input argument. Renaming MODEL.name to ModelName\n');
    MODEL.name      = ModelName;
  end  
else
  MODEL      = GetModelData(file_model);
end

setpref('Internet', 'E_mail', 'clarke@eos.ubc.ca');
setpref('Internet', 'SMTP_Server', 'smtp.eos.ubc.ca');

file_dem     = MODEL.file_dem;

for f=1:numel(MODEL.file_slide)  
  file_slide{f} = MODEL.file_slide{f};
end

file_balance  = MODEL.file_balance;
file_exclude  = MODEL.file_exclude;  % balance exclusion mask (if used)
dx            = MODEL.dx;
dy            = MODEL.dy;
nx            = MODEL.nx;
ny            = MODEL.ny;
N_xy          = MODEL.N_xy;
dt            = MODEL.dt;
n_sub_steps   = MODEL.n_sub_steps;
FLAG_tracer   = 1;
FLAG_archive  = 1;

if BOOTSTRAP
  t_START     = dt*round(BootYear/dt);
else  
  t_START     = MODEL.t_START;
end  

t_STOP        = MODEL.t_STOP;
t_tracer_START  = MODEL.t_tracer_START;
t_archive_START = MODEL.t_archive_START;
m_SLIDE       = MODEL.m_SLIDE;
C_SLIDE       = MODEL.C_SLIDE;
C_SURGE       = MODEL.C_SURGE;
t_CYCLE       = MODEL.t_CYCLE;
THERMAL       = MODEL.THERMAL;
TRACER        = MODEL.TRACER;
TRACER_METHOD = MODEL.TRACER_METHOD;
PLOT          = MODEL.PLOT;
MB_CHECK      = MODEL.MB_CHECK;
VERBOSITY     = MODEL.VERBOSITY;

arc_subdir    = sprintf('Archive-%s', ModelName);

YRSEC      = 365.25*24*3600;
par        = SetParameters(MODEL);

n_XI       = MODEL.n_XI;
n_ZETA     = MODEL.n_ZETA;
ETA_0      = par.ETA_0;

b_dot_FRY  = -5;

disp('==================================================================================')
fprintf(1,'LAUNCHING TRAPRIDGE SURGE MODEL "%s" FOR STRUCTURAL ANALYSIS\n', ModelName);
fprintf(1,'MODEL - Ver 4.00\n\n');

SetupIndexArrays(MODEL);
N_xy = nx*ny;

if BOOTSTRAP
  [B, S, T, T_S, q_G, xd, yd, td] = LoadBootstrap(MODEL, t_START);
  x_grid = MODEL.x_grid;
  y_grid = MODEL.y_grid;
  if ~isempty(xd)
    BOOT.xd      = xd;
    BOOT.yd      = yd;
    BOOT.td      = td;
    BOOT.yr_done = BootYear;
  end  
else  
  load(file_dem, 'B', 'S', 'x_grid', 'y_grid')
  if MODEL.GRID_shift 
    MODEL.x_0 = x_grid(1);
    MODEL.y_0 = y_grid(end);
    MODEL.z_0 = min(min(B));
    S         = S - MODEL.z_0;
    B         = B - MODEL.z_0;
    x_grid    = x_grid - MODEL.x_0;
    y_grid    = y_grid - MODEL.y_0;
  else  
    MODEL.GRID_shift = 0;
    MODEL.x_0 = 0;
    MODEL.y_0 = 0;
    MODEL.z_0 = 0;
  end
  [Y3, Z3, X3] = meshgrid(y_grid, XI, x_grid);
  [X2, Y2]     = meshgrid(x_grid, y_grid); 
  MODEL.x_grid = x_grid;
  MODEL.y_grid = y_grid;
end

switch numel(file_slide)
  case 1
    load(char(file_slide), 'M', 'X_surge_W', 'Y_surge_W', 'X_surge_E', 'Y_surge_E')
    M       = reshape(M, N_xy, 1);  
  otherwise
    M_tmp = NaN(numel(file_slide), ny, nx);
    for p=1:numel(file_slide)
      if p==1
       load(file_slide{p}, 'M', 'X_surge_W', 'Y_surge_W', 'X_surge_E', 'Y_surge_E')
      else    
        load(file_slide{p}, 'M')
      end  
      M_tmp(p,:,:) = M;
    end
    M = reshape(M_tmp, numel(file_slide), N_xy, 1);
end

if exist('X_surge_W', 'var') && exist('X_surge_E', 'var') && exist('Y_surge_W', 'var') && exist('Y_surge_E')
  if isfield(MODEL, 'X_mask_min') && abs(X_surge_W-MODEL.X_mask_min)>1
      fprintf(1,'MODEL.X_mask_min in file "%s.dat" inconsistent with X_surge_W in file "%s" and is reset\n', ModelName, file_slide{1});
      fprintf(1,'  MODEL.X_mask_min = %.2f (discarded)\n', MODEL.X_mask_min);
      fprintf(1,'  X_surge_W        = %.2f (applied)\n', X_surge_W);
      MODEL.X_mask_min = X_surge_W;
  end  
  if isfield(MODEL, 'X_mask_max') && abs(X_surge_E-MODEL.X_mask_max)>1
      fprintf(1,'MODEL.X_mask_max in file "%s.dat" inconsistent with X_surge_E in file "%s" and is reset\n', ModelName, file_slide{1});
      fprintf(1,'  MODEL.X_mask_max = %.2f (discarded)\n', MODEL.X_mask_max);
      fprintf(1,'  X_surge_E        = %.2f (applied)\n', X_surge_E);
      MODEL.X_mask_max = X_surge_E;
  end  
  if isfield(MODEL, 'Y_mask_min') && abs(Y_surge_W-MODEL.Y_mask_min)>1
      fprintf(1,'MODEL.Y_mask_min in file "%s.dat" inconsistent with Y_surge_W in file "%s" and is reset\n', ModelName, file_slide{1});
      fprintf(1,'  MODEL.Y_mask_min = %.2f (discarded)\n', MODEL.Y_mask_min);
      fprintf(1,'  Y_surge_W        = %.2f (applied)\n', Y_surge_W);
      MODEL.Y_mask_min = Y_surge_W;
  end  
  if isfield(MODEL, 'Y_mask_max') && abs(Y_surge_E-MODEL.Y_mask_max)>1
      fprintf(1,'MODEL.Y_mask_max in file "%s.dat" inconsistent with Y_surge_E in file "%s" and is reset\n', ModelName, file_slide{1});
      fprintf(1,'  MODEL.Y_mask_max = %.2f (discarded)\n', MODEL.Y_mask_max);
      fprintf(1,'  Y_surge_E        = %.2f (applied)\n', Y_surge_E);
      MODEL.Y_mask_max = Y_surge_E;
  end      
end

xc_min = x_grid(1);
xc_max = x_grid(end);
yc_max = y_grid(1);
yc_min = y_grid(end);
 
XI_min = XI(1);
XI_max = XI(end);

if ~isempty(MODEL.descriptor)
  fprintf(1,'> %s\n\n', MODEL.descriptor);
end
fprintf(1,'  OMEGA             = %.2f\n', OMEGA);
fprintf(1,'  A_GLEN            = %e\n', A_GLEN);
fprintf(1,'  m_SLIDE           = %.2f\n', m_SLIDE);
fprintf(1,'  C_SLIDE           = %e\n', C_SLIDE);
fprintf(1,'  C_SURGE           = %e\n', C_SURGE);
fprintf(1,'  ETA_0             = %e\n', ETA_0);
fprintf(1,'  nx                = %d\n', nx);
fprintf(1,'  ny                = %d\n', ny);
fprintf(1,'  nx                = %.1f m\n', dx);
fprintf(1,'  dy                = %.1f m\n', dy);
fprintf(1,'  t_START           = %.1f yr\n', t_START);
fprintf(1,'  t_tracer_START    = %.1f yr\n', t_tracer_START);
fprintf(1,'  t_archive_START   = %.1f yr\n', t_archive_START);
fprintf(1,'  t_STOP            = %.1f yr\n', t_STOP);
fprintf(1,'  dt                = %.3f yr\n', dt);
fprintf(1,'  n_sub_steps       = %d\n', n_sub_steps);
fprintf(1,'  n_XI              = %d\n', n_XI);
fprintf(1,'  t_CYCLE           = %.2f yr\n', MODEL.t_CYCLE);
fprintf(1,'  dt_SURGE          = %.2f yr\n', MODEL.dt_SURGE);
fprintf(1,'  dt_PEAK           = %.2f yr\n', MODEL.dt_PEAK);
fprintf(1,'  dL_front          = %.1f m\n', MODEL.dL_front);
fprintf(1,'  THERMAL           = %d\n', MODEL.THERMAL);
if THERMAL
  fprintf(1,'  n_ZETA            = %d\n', n_ZETA);
end
fprintf(1,'  TRACER            = %d\n', MODEL.TRACER);

if TRACER
  fprintf(1,'  TRACER_METHOD     = %s\n', MODEL.TRACER_METHOD);
end  

fprintf(1,'  MB_CHECK          = %d\n', MODEL.MB_CHECK);
fprintf(1,'  VERBOSITY         = %d\n', MODEL.VERBOSITY);

disp('==================================================================================')

A_const = 2*A_GLEN*(RHO*g)^n_GLEN/(n_GLEN+2);
C_const = C_SLIDE*(RHO*g)^m_SLIDE;                 %% C_surge vs C_slide

load(file_balance, 'b_dot');

if ~isempty(file_exclude)
  load(file_exclude, 'I_exclude')
  b_dot(I_exclude==1) = b_dot_FRY;
end

[x_ELA, y_ELA] = GetELA(x_grid, y_grid, b_dot);

timestamp = datestr(now,31);

b_dot       = reshape(b_dot, N_xy, 1);

if THERMAL
  if ~BOOTSTRAP
    [T, T_S, q_G] = InitializeThermalModel(reshape(B, numel(B), 1), reshape(S, numel(S), 1));
  else
    ZETA   = linspace(0, 1, n_ZETA)';    % Note NEW (Sept 2017) definition of ZETA matches that of Greve & Blatter and is opposite to Jenssen and to Marshall. ZETA = (z-B)/H in the NEW system
    d_ZETA = 1/(n_ZETA-1);
  end
else
  T.A     = [];
  T.grid = {};
  T_S    = [];
  q_G    = [];
end

MODEL.q_G = q_G;
MODEL.T_S = T_S;

% Display some diagnostics before proceeding with model integration

t       = t_START;

% Assign times for command window output, plot generation and full matfile output of debris and ice DEMs

dt_PLOT = min(10, t_STOP);
dt_OUT  = min(10, t_STOP);

t_OUT   = t_START+dt_OUT;
t_PLOT  = t_START+dt_PLOT;

FLAG    = true;

B       = reshape(B, numel(B), 1);
S       = reshape(S, numel(S), 1);
b_dot   = reshape(b_dot, numel(b_dot), 1);

TEST         = true;
ARCHIVE_init = false;

t_test  = [40 100 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 5000];
t_test  = t_test(t_test>=t_START);
 
i_test = 1;

IPREC  = GetTimeStepPrecision(dt);

handle = 100;

while 1

  % Allow multiple sliding masks to enable intricate surge patterns

  C_t = zeros(N_xy, 1);
  
  for p=1:numel(t_CYCLE)
    [M_t, MODEL] = ActivateDeactivateSlidingMask(p, M, X2, Y2, MODEL, t);         
    C_t          = C_t + C_SLIDE + (C_SURGE-C_SLIDE)*reshape(M_t, N_xy, 1);
  end
    
  [S, H, yr_done, LAMBDA_max, T, V, L, DD, E] = step_multi_year(ModelName, S, B, b_dot, M, C_t,  T, T_S, q_G, dt, n_sub_steps, t);
  
  % In addition to archive save output every 100 years to facilitate bootstrapping
  
  if yr_done>=t_archive_START || mod(round(yr_done, IPREC), 100)==0

    if FLAG_archive && yr_done>=t_archive_START
      fprintf('\nSTARTING RUN ARCHIVE AT %.2f <%s>\n\n', yr_done, datestr(now));
      FLAG_archive = 0;
    end
        
    if ~ARCHIVE_init
      dir_arch   = fullfile(pwd, arc_subdir);
      dir_arch_out = fullfile(dir_arch, 'Snapshots');
      
      if ~exist(dir_arch, 'dir')
        mkdir(dir_arch)
      end
      if ~exist(dir_arch_out, 'dir')
        mkdir(dir_arch_out)
      end
 
      file_mat   = fullfile(dir_arch, 'Model.mat');
      MODEL.RunTime_hr = toc(t_start)/3600;
      SaveModelParameters(file_mat, MODEL, par, file_source, file_date)
      copyfile(file_model, fullfile(dir_arch, file_model))   % copy text file to run archive
      
      ARCHIVE_init=1;
    end
    
    if BOOTSTRAP
      file_arch_out = fullfile(dir_arch_out, sprintf('BootOut(%08.3f).mat', yr_done));     
    else
      file_arch_out = fullfile(dir_arch_out, sprintf('Out(%08.3f).mat', yr_done));     
    end  
    SaveEnglacialArchive(file_arch_out, yr_done, H, V, L, DD, T, E, M_t);
  end    
   
  t = yr_done;
  
  u_ic_jc = V.u_ic_jc;
  v_ic_jc = V.v_ic_jc;
  w_ic_jc = V.w_ic_jc;
  
  v_max   = max(max(sqrt(u_ic_jc.^2+v_ic_jc.^2)));
  
  if t>=t_OUT-0.5*dt
    t_OUT = t_OUT+dt_OUT;
    fprintf(1,'t=%8.2f yr: max(H)=%.2f m; ALPHA_I=%.2f%%; VOL_I=%.5f km^3; <%s>\n', ...
       t, max(S-B), 100*sum(S>B)/N_xy, dx^2*sum(S-B)/10^9, datestr(now));
  end
  
  if PLOT && t>=t_PLOT-0.5*dt
    t_PLOT = t_PLOT+dt_PLOT;      
    H_plot = S-B;
    
    if ~exist('X_site', 'var')
      file_site = fullfile(fullfile(pwd, arc_subdir), 'SiteList.dat');
      if exist(file_site, 'file')
        [Grid, t_site, t_last_surge, X_site, Y_site, XI_site, dZ_tweak_site, Name_site, Type_site, ~] = GetSiteList(file_site, MODEL);
        if strcmp(upper(Grid), 'WGS84')
          fprintf(1,'Converting site coordinates to NAD27 from WGS84 using jiffy GKCC shift\n');
          
          X_WGS84_site = X_site;
          Y_WGS84_site = Y_site;
          
          X_site = X_WGS84_site + 110.557;
          Y_site = Y_WGS84_site -165.460;
        end
        
        X_site  = X_site - MODEL.x_0;
        Y_site  = Y_site - MODEL.y_0;
      else
        X_site = [];
      end
    end
          
    if ~exist('handle_1', 'var')
      handle_1      = 1;
    else
      if handle_1==1
        handle_1 = 3;
      elseif handle_1==3
        handle_1 = 1;
      else
        error('main(): Unprogrammed value of handle_1')
      end  
    end  
    
    if max(H_plot)-min(H_plot)>0.01
      figure(handle_1)
     
      H_plot(H_plot==0) = NaN;
      H_plot = reshape(H_plot, ny, nx);
      imagescnan(x_grid, y_grid, H_plot), axis equal, axis tight, colorbar
      set(gca, 'YDir', 'normal')
      hold on
      plot(x_ELA, y_ELA, 'y')
      
      if ~isempty(X_site)
        plot(X_site, Y_site, 'wo', X_site, Y_site, 'k+')
      end  
      hold off
      
      title(sprintf('%s: Ice thickness (m asl) at %6.1f yr',ModelName, t), 'interpreter', 'none');
      drawnow
    
      figure(handle_1+1)
    
      v_surf = reshape(squeeze(sqrt(u_ic_jc(n_XI,:).^2 + v_ic_jc(n_XI,:).^2 + w_ic_jc(n_XI,:).^2)), ny, nx);  % XI(1) is surface value
    
      v_surf(isnan(H_plot)) = NaN;
      imagescnan(x_grid, y_grid, v_surf), axis equal, axis tight, colorbar
      set(gca, 'YDir', 'normal')
      hold on
      plot(x_ELA, y_ELA, 'y')
      if ~isempty(X_site)
        plot(X_site, Y_site, 'wo', X_site, Y_site, 'k+')
      end  
      
      % Draw the endpoints of the ice mask (if specified)
      
      if isfield(MODEL, 'X_mask_min')
        hold on
        plot([MODEL.X_mask_min-MODEL.x_0  MODEL.X_mask_max-MODEL.x_0], [MODEL.Y_mask_min-MODEL.y_0 MODEL.Y_mask_max-MODEL.y_0], 'y')
        plot([MODEL.X_mask_min-MODEL.x_0  MODEL.X_mask_max-MODEL.x_0], [MODEL.Y_mask_min-MODEL.y_0 MODEL.Y_mask_max-MODEL.y_0], 'kd')         
      end  
      hold off
      title(sprintf('%s: Surface flow speed (m/yr) at %6.1f yr',ModelName, t), 'interpreter', 'none');
      drawnow
    end
  end
  
  if t>=t_STOP
    break
  end

end

fprintf('\nALL DONE: Forward modelling of MODEL "%s". Elapsed time %.2f hr\n', ModelName, toc(t_start)/3600)

sendmail('clarke@eos.ubc.ca', sprintf('Stucture evolution run "%s" completed', ModelName), ...
  sprintf('Model %s completed at %s. Elapsed time %.2f hr.', ModelName, datestr(now), toc(t_start)/3600));

