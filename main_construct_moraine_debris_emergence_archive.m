function main_constuct_moraine_debris_emergence_archive(ModelName)

% This script is mainly useful for generating supraglacial moraines

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
  ModelDirName = fullfile(pwd, ModelDirName);
else
  clear global
  clearvars -except ModelName
  ModelDirName = fullfile(pwd, sprintf('Archive-%s', ModelName));
end 

global dt dx dy nx ny N_xy X2 Y2 X3 Y3 Z3 XI n_XI d_XI

fprintf(1,'\nLAUNCHING :: main_construct_moraine_debris_emergence_archive(''%s'') - Englacial particle tracking from source to emergence points\n', ModelName);

setpref('Internet', 'E_mail', 'clarke@eos.ubc.ca');
setpref('Internet', 'SMTP_Server', 'smtp.eos.ubc.ca');

VERBOSE    = 0;
PLOT       = 1;
t_start    = tic;

file_date      = datestr(now);
file_source    = fullfile(pwd, mfilename);

ModelFileName  = fullfile(ModelDirName, 'Model.mat');
SourceFileName = fullfile(ModelDirName, 'SourceList.dat');
SiteFileName   = fullfile(ModelDirName, 'SiteList.dat');

dir_arch     = fullfile(pwd, sprintf('Archive-%s', ModelName));
dir_SNAP     = fullfile(dir_arch, 'Snapshots');

dir_Moraine  = fullfile(dir_arch, 'MoraineTracks');

% Check existence of required files and directories

OK = CheckExistence({SourceFileName, ModelFileName, dir_SNAP});

if ~OK
  return
end 

fprintf(1,'\n01. Initializing SOURCE structure describing non-ice tracer input locations\n');

load(ModelFileName, 'MODEL')

nx      = MODEL.nx;
ny      = MODEL.ny;
N_xy    = MODEL.N_xy;
XI      = MODEL.XI;
n_XI    = MODEL.n_XI;
d_XI    = MODEL.d_XI;
N_xy    = MODEL.N_xy;

dx      = MODEL.dx;
dy      = MODEL.dy;
dt      = MODEL.dt;

t_CYCLE = MODEL.t_CYCLE;

[X3, Y3, Z3] = MakeGrids(MODEL);
X3           = X3.ic_jc_kc;
Y3           = Y3.ic_jc_kc;
Z3           = Z3.ic_jc_kc;
[X2, Y2]     = meshgrid(MODEL.x_grid, MODEL.y_grid);

load(MODEL.file_dem, 'B');

load(MODEL.file_balance, 'b_dot')

[SourceGrid, X_NAD27_source, Y_NAD27_source, X_WGS84_source, Y_WGX84_source, Name_source, YR_start, YR_stop, t_last_surge_alt] = GetSourceList(SourceFileName, MODEL);

X_source = X_NAD27_source;
Y_source = Y_NAD27_source;

if MODEL.GRID_shift   
  X_source = X_source-MODEL.x_0;
  Y_source = Y_source-MODEL.y_0;
end

% Get t_last_surge structure from SiteList file

[~, ~, t_last_surge, ~, ~, ~, ~, ~, ~, ~, ~] = GetSiteList(SiteFileName, MODEL);

if t_last_surge.start~=t_last_surge_alt.start || t_last_surge.stop~=t_last_surge_alt.stop
  error('main_construct_moraine_debris_emergence_archive(): t_last_surge values differ between "%s" and "%s"', SourceFileName, SiteFileName)
end

YR_view   = 2006.6;
YR_view   = round(YR_view, GetTimeStepPrecision(dt));    % Round to a computed model time

YR_start  = YR_view - t_CYCLE + dt;
YearList  = YR_start:dt:YR_view;

DIR          = dir(fullfile(dir_SNAP, '*.mat'));

ModelFileList = {DIR.name};
ModelYearList = ExtractYearListFromFileList(ModelFileList);

[DD, ~, ~, ~, ~, YearSwap] = GetAliasYearData(YR_view, t_last_surge, dir_SNAP, ModelFileList, ModelYearList, MODEL, 0);

DD_ref      = DD;
clear DD

SOURCE = InitializeSourceStructure(X_source, Y_source, Name_source, t_last_surge, MODEL);

clear X_source Y_source Name_source

% Prune the unused fields from structure SOURCE

SOURCE = rmfield(SOURCE, {'B_p', 'H_p', 'S_p', 'X_p', 'Y_p', 'Z_p', 'XI_p', 'u_p', 'v_p', 'w_p'});

% Shift the vertical grid to be consistent with X,Y shifting of source
% coordinates performed in the above function call

if MODEL.GRID_shift
  B = B-MODEL.z_0;
end  

ns    = numel(SOURCE);
nt    = numel(YearList);

B_inj = interp2(X2, Y2, B, [SOURCE.X_source], [SOURCE.Y_source]);

X_p    = [];
Y_p    = [];
Z_p    = [];
S_p    = [];    % Only for debug tracking information

Site_p = [];
t_p    = [];
t_model_p = [];

X_e    = [];
Y_e    = [];
Z_e    = [];
t_e    = [];
t_model_e  = [];

t_in       = [];
t_model_in = [];

Site_e = [];   % This is the SOURCE site for an emergent debris point

DONE   = 0;

t      = 0;

fprintf(1,'02. Inserting flow markers at %d injection points and %d time steps\n', ns, nt);

INJECT_cnt = max(round(0.1/dt),1);

while 1
  if mod(t,INJECT_cnt)==0 && t<=nt
    INJECT = 1;
 else
    INJECT = 0;
  end  
  t = t+1;
  
  if t>nt
    YearGet = YearList(nt)+(t-nt)*dt;
    YearGet = round(YearGet, GetTimeStepPrecision(dt));
  else  
    YearGet = YearList(t);
  end
  
  if VERBOSE>1 || mod(YearGet, 10)<0.1*dt
    fprintf(1,'At start of time step t=%d: At year = %.2f with %d englacial markers and %d emerged markers (total=%d)\n', t, YearGet, numel(X_p), numel(X_e), numel(X_p)+numel(X_e));
  end
  
  clear DD V
  
  [DD, V, ~, ~, ~, YearSwap] = GetAliasYearData(YearGet, t_last_surge, dir_SNAP, ModelFileList, ModelYearList, MODEL, VERBOSE);
  
  H            = zeros(ny, nx);
  L_I          = false(ny, nx);
  L_I(DD.ic_jc_I) = true;
  H(L_I)       = DD.H;
     
  v_x          = zeros(n_XI, ny, nx);
  v_y          = zeros(n_XI, ny, nx);
  v_z          = zeros(n_XI, ny, nx);
  
  v_x(:,L_I)   = V.x;
  v_y(:,L_I)   = V.y;
  v_z(:,L_I)   = V.z;
  
  % The following is inefficient but perhaps tolerable
  
  if INJECT
    Site_inj   = 1:ns;
    X_inj      = [SOURCE.X_source];
    Y_inj      = [SOURCE.Y_source];
    H_inj      = interp2(X2, Y2, H, X_inj, Y_inj);
    S_inj      = B_inj+H_inj;
    Z_inj      = S_inj;
    
    t_inj       = YearGet*ones(1, ns);
    t_model_inj = YearSwap*ones(1, ns);
  else
    clear Site_inj X_inj Y_inj H_inj S_inj Z_inj t_inj t_model_inj
  end
  
  if t>1
    if mod(t,100)==0
      fprintf(1,'%d\n', round(YearGet));
    elseif mod(t,10)==0
      fprintf(1,'.');
    end
        
    X_m = X_p;
    Y_m = Y_p;
    Z_m = Z_p;
    
    L_ok = true(size(X_m));
       
    [X_p, Y_p, Z_p, XI_p, L_ok_p] = ForwardFullStep_1S2T(X_m, Y_m, Z_m, B, H, b_dot, v_x, v_y, v_z, YearGet, L_ok);
    
    L_done = ~L_ok_p;
       
    if any(L_done)
      
      L_bad   = Site_e==4 & X_e<1000;
      
      if any(L_bad)
        fprintf(1,'ACHTUNG(1) at t=%d\n', t)
      end  

      X_exit = X_m(L_done);   % Set the exiting values to the previous "live" values since X_p(L_done) etc are set to NaN in ForwardFullStep routine
      Y_exit = Y_m(L_done);
      
      B_exit = interp2(X2, Y2, B, X_m(L_done), Y_m(L_done));  % Debug lint
      H_exit = interp2(X2, Y2, H, X_m(L_done), Y_m(L_done));  % Debug lint
      S_exit = B_exit+H_exit;                       % Debug lint

      Z_exit = S_exit;
      t_inp  = t_p(L_done);
      t_exit       = YearGet*ones(1, sum(L_done));
      t_model_exit = YearSwap*ones(1,sum(L_done));
      Site_exit    = Site_p(L_done);
      
      X_p = X_p(~L_done);
      Y_p = Y_p(~L_done);
      Z_p = Z_p(~L_done);
      S_p = S_p(~L_done);    % Debug lint
      t_p = t_p(~L_done);
      t_model_p = t_model_p(~L_done);
      Site_p    = Site_p(~L_done);
      
      X_e     = [X_e X_exit];
      Y_e     = [Y_e Y_exit];
      Z_e     = [Z_e Z_exit];
      t_e     = [t_e t_exit];
      t_model_e = [t_model_e t_model_exit];
      t_in    = [t_in t_inp];
      t_model_in = [t_model_in t_model_p];
      Site_e  = [Site_e Site_exit];
      
      L_bad   = Site_e==4 & X_e<1000;
      
      if any(L_bad)
        fprintf(1,'ACHTUNG(2) at t=%d\n', t)
      end  
      
    end  
  
    if all(L_done)
      fprintf(1,'04. All englacial flow markers have exited the glacier at loop count t=%d\n', t);
      DONE = 1;
    end  
  end
  
  % Add marker points at each source site and time step until a full surge
  % cycle has been completed. Then stop adding markers and simply continue
  % tracking the englacial markers until none remain
  
  if INJECT
    X_p    = [X_p X_inj];
    Y_p    = [Y_p Y_inj];
    Z_p    = [Z_p Z_inj];
    S_p    = [S_p Z_inj];    
    t_p    = [t_p t_inj];
    t_model_p = [t_model_p t_model_inj];
    Site_p = [Site_p Site_inj];
    
    if t==nt
      fprintf(1,'03. Insertion of flow markers completed\n');
    end  
  end
  
  if DONE
    break
  end  
end

fprintf(1,'05. Organizing and saving flow marker emergence results\n');

for s=1:ns
  
  % Extract results for each site
  
  L_s   = Site_e==s;
  Xe_s  = X_e(L_s);
  Ye_s  = Y_e(L_s);
  Ze_s  = Z_e(L_s);
  te_s  = t_e(L_s);
  tin_s = t_in(L_s);
    
  te_model_s = t_model_e(L_s);
  
  % Sort by emergence time
  
  [~, i_sort] = sort(te_s, 'ascend');
   
  SOURCE(s).X_e  = Xe_s(i_sort);
  SOURCE(s).Y_e  = Ye_s(i_sort);
  SOURCE(s).Z_e  = Ze_s(i_sort);
  SOURCE(s).t_e  = te_s(i_sort);
  SOURCE(s).t_model_e = te_model_s(i_sort);
  SOURCE(s).t_in = tin_s(i_sort);

  % For emergence times that exceed YR_view swap alias times 
  
  while SOURCE(s).t_e>YR_view
    SOURCE(s).t_e  = SOURCE(s).t_e-t_CYCLE;
    SOURCE(s).t_in = SOURCE(s).t_in-t_CYCLE;
  end
end

if ~exist(dir_Moraine, 'dir')
  mkdir(dir_Moraine)
end  

file_mat = fullfile(dir_Moraine, 'DEBRIS_EXIT.mat');
DD        = DD_ref;

save(file_mat, 'MODEL', 'SOURCE', 'DD', 'file_date', 'file_source')

if PLOT
  handle = 0;
  handle = PlotDebrisEmergencePoints(handle, MODEL, SOURCE, DD, dir_Moraine, 1);
end

fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);

sendmail('clarke@eos.ubc.ca', sprintf('main_construct_moraine_debris_emergence_archive(): Construction of moraine debris emergence archive "%s" completed', ModelName), ...
  sprintf('Model %s completed at %s. Elapsed time %.2f hr.', ModelName, datestr(now), toc(t_start)/3600));
