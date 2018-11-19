function main_bottom_up_particle_tracking(ModelName, SiteFileName)

% This script is to calculate the trajectory of flow markers from their downstream 
% endpoints to the their upstream surface sources

% All potential sites are listed and these are grouped by the Site_type vectors
% So they can be subsequently selected for particular analyses (e.g.
% surface maps, sites along specified flow-lines and sites on a flow
% cross-section. The purpose of the grouping is to avoid multiple runs to
% determine all possible site trajectories

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
  SiteFileName  = 'SiteList.dat';
  SiteShortName = 'SiteList.dat';
  SiteFileName  = fullfile(pwd, sprintf('Archive-%s', ModelName), SiteFileName);
elseif nargin==1
  clear global
  clearvars -except ModelName
  SiteFileName  = 'SiteList.dat';
  SiteShortName = 'SiteList.dat';
  SiteFileName  = fullfile(pwd, sprintf('Archive-%s', ModelName), SiteFileName);  
end  

global nx ny dx dy dt X2 Y2 X3 Y3 Z3 RHO g XI d_XI n_XI

fprintf(1,'\nLAUNCHING :: main_bottom_up_particle_tracking(''%s'') - Particle tracking from sampling sites to upstream deposition sites\n', ModelName);

t_start = tic;

VERBOSE      = 0;

setpref('Internet', 'E_mail', 'clarke@eos.ubc.ca');
setpref('Internet', 'SMTP_Server', 'smtp.eos.ubc.ca');

dir_ROOT   = fullfile(pwd, sprintf('Archive-%s', ModelName));
dir_SNAP   = fullfile(dir_ROOT, 'Snapshots');

dir_OUT  = dir_ROOT;

file_model = fullfile(pwd, sprintf('Archive-%s', ModelName), 'Model.mat');

load(file_model, 'MODEL')

SetupIndexArrays(MODEL)

RHO     = MODEL.par.RHO;
g       = MODEL.par.g;

XI      = MODEL.XI;
d_XI    = 1./(MODEL.n_XI-1);

fprintf(1,'==================================================================================\n');
fprintf(1,'BOTTOM-UP PARTICLE PATHS USING %s SCHEME FOR MODEL %s\n\n', MODEL.TRACER_METHOD, ModelName);
if isfield(MODEL, 'descriptor') && ~isempty(MODEL.descriptor)
  fprintf(1,'> %s\n\n', MODEL.descriptor);
end
if isfield(MODEL, 'OMEGA') && ~isempty(MODEL.OMEGA)
  fprintf(1,'  OMEGA       = %.2f\n', MODEL.OMEGA);
end
if isfield(MODEL, 'A_GLEN') && ~isempty(MODEL.A_GLEN)
  fprintf(1,'  A_GLEN      = %e\n', MODEL.A_GLEN);
end  
fprintf(1,'  C_SLIDE     = %e\n', MODEL.C_SLIDE);
fprintf(1,'  C_SURGE     = %e\n', MODEL.C_SURGE);
fprintf(1,'  nx          = %d\n', MODEL.nx);
fprintf(1,'  ny          = %d\n', MODEL.ny);
fprintf(1,'  nx          = %.1f m\n', MODEL.dx);
fprintf(1,'  dy          = %.1f m\n', MODEL.dy);
fprintf(1,'  t_START     = %.1f yr\n', MODEL.t_START);
fprintf(1,'  t_STOP      = %.1f yr\n', MODEL.t_STOP);
fprintf(1,'  dt          = %.2f yr\n', MODEL.dt);
fprintf(1,'  n_sub_steps = %d\n', MODEL.n_sub_steps);
fprintf(1,'  n_XI        = %d\n', MODEL.n_XI);
fprintf(1,'  t_CYCLE     = %.2f yr\n', MODEL.t_CYCLE);
fprintf(1,'  dt_SURGE    = %.2f yr\n', MODEL.dt_SURGE);
fprintf(1,'  dt_PEAK     = %.2f yr\n', MODEL.dt_PEAK);
fprintf(1,'  dL_front    = %.1f m\n', MODEL.dL_front);
fprintf(1,'  THERMAL     = %d\n', MODEL.THERMAL);

if MODEL.THERMAL
  fprintf(1,'  n_ZETA      = %d\n', MODEL.n_ZETA);
end

if exist('FOLD_INDEX', 'var')
  fprintf(1,'  FOLD_INDEX        = %d\n', MODEL.FOLD_INDEX);
end

fprintf(1,'==================================================================================\n');

fprintf(1,'01. Initializing SITE structure using "%s"\n', SiteFileName);

if ~isfield(MODEL, 'GRID_shift')
  GRID_shift = 0;
  x_0        = 0;
  y_0        = 0;
  z_0        = 0;
  MODEL.GRID_shift = GRID_shift;
  MODEL.x_0  = x_0;
  MODEL.y_0  = y_0;
  MODEL.z_0  = z_0;
else
  GRID_shift = MODEL.GRID_shift;
  x_0        = MODEL.x_0;
  y_0        = MODEL.y_0;
  z_0        = MODEL.z_0;
end  

nx        = MODEL.nx;
ny        = MODEL.ny;
n_XI      = MODEL.n_XI;

dx        = MODEL.dx;
dy        = MODEL.dy;
x_grid    = MODEL.x_grid;
y_grid    = MODEL.y_grid;

dt        = MODEL.dt;
t_CYCLE   = MODEL.t_CYCLE;

N_xy      = nx*ny;
N_all     = N_xy*n_XI;

[Y3, Z3, X3] = meshgrid(y_grid, XI, x_grid);
[X2, Y2]     = meshgrid(x_grid, y_grid);

load(MODEL.file_dem, 'B') 

load(MODEL.file_balance, 'b_dot')

[SiteGrid, t_site, t_last_surge, X_NAD27_site, Y_NAD27_site, X_WGS84_site, Y_WGS84_site, XI_site, dZ_tweak_site, Name_site, Type_site, Group] = GetSiteList(SiteFileName, MODEL);

X_site = X_NAD27_site;
Y_site = Y_NAD27_site;

if MODEL.GRID_shift
  X_site = X_site - x_0;
  Y_site = Y_site - y_0;
  B      = B - z_0;
end

DIR = dir(fullfile(dir_SNAP, '*.mat'));

ModelFileList = {DIR.name};
ModelYearList = ExtractYearListFromFileList(ModelFileList);

IPREC    = GetTimeStepPrecision(dt);

YR_stop = t_site;  
[FileSwap, YR_swap] = GetYearSwapFileName(MODEL, t_last_surge, ModelFileList, ModelYearList, YR_stop);   % Use alias years to match stop year with an output file
  
SITE = InitializeSiteStructure(t_site, t_last_surge, X_site, Y_site, XI_site, dZ_tweak_site, Name_site, Type_site, Group, ...
       b_dot, fullfile(dir_SNAP, FileSwap), MODEL);

YR_backtrack = 800;
YR_backstop  = t_site-YR_backtrack; 
  
% Round to nearest year that matches a viable file name

YR_start = round(SITE(1).t_site, IPREC);   % Assume all sampling sites have same year 
YR_min   = round(YR_backstop, IPREC);
    
YearList = YR_start:-dt:YR_min;      % This list runs backward in time (from site time to deposition time)  

ns   = numel(SITE);

% Check SITES to ensure that all are in the ablation area

if any([SITE.b_dot_site]>0)
  for s=1:ns
    if SITE(s).b_dot_site>0
      fprintf(1,'main_bottom_up_nonreversed(): BAD SITE data for SITE(%d) in file "%s". Site is in accumulation region\n', s, SiteShortName);
    end
  end  
  error('Program halted')
end

t_step   = tic;
t_int    = tic;

site_seq  = 1:ns;

fprintf(1,'02. Extracting and setting initial values for backtracking archive\n');

[DD, V, ~, ~, ~, YearSwap] = GetAliasYearData(YearList(1), t_last_surge, dir_SNAP, ModelFileList, ModelYearList, MODEL);

[DD, V, ~, ~] = CheckAndExpand(DD, V);
   
H            = DD.H;
S            = B+H;

handle = 1;
figure(handle)
H_plot = H;
H_plot(H==0) = NaN;
      
imagescnan(x_grid, y_grid, H_plot), axis equal, colorbar
set(gca, 'YDir', 'normal')
hold on
plot(X_site, Y_site, 'wo', X_site, Y_site, 'k+');
hold off

X_pt         = [SITE.X_site];
Y_pt         = [SITE.Y_site];
Z_pt         = [SITE.Z_init_site];

H_pt         = interp2(X2, Y2, H, X_pt, Y_pt);

L_bad        = H_pt<=0;

if any(L_bad)
  fprintf(1,'\nmain_bottom_up_particle_tracking(): %d site markers lie beyond the ice margin at YR_start=%.2f yr. These are deleted from input SITE structure\n\n', sum(L_bad), YR_start);
  X_pt      = X_pt(~L_bad);
  Y_pt      = Y_pt(~L_bad);
  Z_pt      = Z_pt(~L_bad);
  H_pt      = H_pt(~L_bad);
  SITE      = SITE(~L_bad);
  ns        = numel(SITE);
end
% %%%%%%%%%%%%%%%%%%%%%%%%
% i_use = [274 303];
% X_pt  = X_pt(i_use);
% Y_pt  = Y_pt(i_use);
% Z_pt  = Z_pt(i_use);
% H_pt  = H_pt(i_use);
% SITE  = SITE(i_use);
% ns    = numel(SITE);
%%%%%%%%%%%%%%%%%%%%
B_pt         = interp2(X2, Y2, B, X_pt, Y_pt);
S_pt         = B_pt + H_pt;

L_out        = Z_pt>S_pt;
Z_pt(L_out)  = S_pt(L_out);

% Initialize or reinitialize P_up structure

L_OK         = true(1, ns);
nt           = numel(YearList);

XI_pt        = [SITE.XI_init_site];

t_pt         = YR_start;   % Value of t_site rounded to correspond to a time step;

P_up(ns)     = struct('t', [], 'X', [], 'Y', [], 'Z', [], 'XI', [], 'S', []);
P_up         = PlugValues_to_P_up_Structure(P_up, t_pt, X_pt, Y_pt, Z_pt, XI_pt, S_pt);

X_m          = X_pt;
Y_m          = Y_pt;
Z_m          = Z_pt;
 
fprintf(1,'03. Backtrack initialization completed\n');
  
DONE  = 0;

for t=2:nt
   
  if all(~L_OK)
    DONE = 1;
    break
  end
  
  i_site = site_seq(L_OK);   % List of "active" sites for backtracking
    
  t_step = tic;
   
  [DD, V, ~, ~, ~, YearSwap] = GetAliasYearData(YearList(t), t_last_surge, dir_SNAP, ModelFileList, ModelYearList, MODEL, VERBOSE);
  [DD, V, ~, ~] = CheckAndExpand(DD, V);
       
  H          = DD.H;
  S          = B+H; 
  
  X_pt     = X_m;
  Y_pt     = Y_m;
  Z_pt     = Z_m;
    
  t_pt     = YearList(t);
  
  L_OK_save = L_OK;
  
  [X_m, Y_m, Z_m, XI_m, S_m, L_OK] = BackwardFullStep_1S2T(X_pt, Y_pt, Z_pt, B, H, b_dot, V.x, V.y, V.z, t_pt, L_OK); 
  
  % Test for an emerged marker and check that the ice column at the emergence point is not trivially thin
  % (in which case purge this point
     
  if ~all(L_OK==L_OK_save)
    L_pt   = L_OK_save==~L_OK;
    cnt    = 1:numel(L_OK);
    i_pt   = cnt(L_pt);
    B_m_pt = interp2(X2, Y2, B, X_m(i_pt), Y_m(i_pt));
    H_m_pt = S_m(i_pt)-B_m_pt;
    
    L_thin = H_m_pt<1;
    i_thin = i_pt(L_thin);
    L_OK_save(i_thin) = 0;
  end  
  
  % Test for an emerged marker and check the mass balance at the emergence
  % point (we require emergence to occur in accumulation zone)

  if ~all(L_OK==L_OK_save)
    L_test = ~(L_OK==L_OK_save);
    b_dot_test = interp2(X2, Y2, b_dot, X_m(L_test), Y_m(L_test));

    if any(b_dot_test<0)
      error('main_bottom_up_particle_tracking(): Emerged marker(s) #%d are in ablation zone. Program halted. Adjust dZ_tweak_site in "%s"\n', find(L_test), SiteFileName)
    end
  end    
  
  P_up         = PlugValues_to_P_up_Structure(P_up, t_pt, X_m, Y_m, Z_m, XI_m, S_m);

  if mod(YearList(t), 10)<1.0e-06
    fprintf(1,' Processing model year %.2f (swap year %.2f) with %d active site(s). Interval time=%.1f min; Total elapsed time=%.2f h <%s>\n', YearList(t), YearSwap, sum(L_OK), ...
      toc(t_int)/60, toc(t_start)/3600, datestr(now));
    t_int = tic;
  elseif mod(YearList(t), 1)<1.0e-06 && ~VERBOSE
    fprintf(1,'.');
  end  

end

if all(~L_OK)
  fprintf(1,'04. Backtracking completed\n');
else
  fprintf(1,'\n*** WARNING ***\n');
  fprintf(1,'  YearList(nt)=%.2f endpoint reached YR_backstop=%.2f with %d active backtracking sites.\n', YearList(nt), YR_backstop, sum(L_OK));
  fprintf(1,'  Parameter YR_backtrack=%.2f must be set to a larger value and script rerun\n\n', YR_backtrack);
  fprintf(1,'  The marker identifications and XI values for unfinished markers are as follows:\n')
  
  site_OK    = {SITE(L_OK).Name_site};
  P_up_OK    = P_up(L_OK);
  
   for p=1:numel(P_up_OK)
    fprintf(1,'    %-12s  %7.4f\n', char(site_OK{p}), P_up_OK(p).XI(end));
  end  
end

fprintf(1,'05. Checking backtrack trajectory to ensure marker emergence is in accumulation zone and not in very thin ice\n')

b_dot_test = NaN(1, ns);
B_test     = NaN(1, ns);
Z_test     = NaN(1, ns);

for s=1:ns
  b_dot_test(s) = interp2(X2, Y2, b_dot, P_up(s).X(end), P_up(s).Y(end));
  B_test(s)     = interp2(X2, Y2, B, P_up(s).X(end), P_up(s).Y(end));
  Z_test(s)     = P_up(s).Z(end);
end

% First purge any trajectories that exit because of zero ice thickness or which terminate in the ablation zone (e.g., because of stalling)

L_thin  = Z_test-B_test < 1;
L_abl   = b_dot_test<=0;    % Exiting in ablation zone with non-vanishing ice column is not permitted

L_remove = L_thin | L_abl;

P_up    = P_up(~L_remove);
SITE    = SITE(~L_remove); 

L_fail  = L_abl & ~L_thin;

if any(L_fail)
  fprintf(1, 'main_bottom_up_particle_tracking(): One or more tracking sites leads to erroneous emergence in ablation zone\n\n')
  is_fail = site_seq(L_fail);
  for is=is_fail
    fprintf(1,'    Site #%d : X=%.3f, Y=%.3f, b_dot=%.3f\n', is, P_up(is).X(end), P_up(is).Y(end), b_dot_test(is));
  end
end

file_mat = fullfile(dir_ROOT, 'UP_TRACKS.mat');

file_source = fullfile(pwd, mfilename);
file_date   = datestr(now);
 
save(file_mat, 'MODEL', 'SITE', 'P_up', 'file_source', 'file_date', '-v7.3')
 
fprintf(1,'\nALL DONE. Trajectories extracted and processed for Model %s. Elapsed time %.2f h\n', ModelName, toc(t_start)/3660);

sendmail('clarke@eos.ubc.ca', sprintf('main_bottom_up_particle_tracking(): Bottom up particle tracking run "%s" completed', ModelName), ...
  sprintf('Model %s completed at %s. Elapsed time %.2f hr.', ModelName, datestr(now), toc(t_start)/3600));
