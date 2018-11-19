function main_calculate_and_plot_strain_foliation_S1_map(ModelName)

% This script plot the strain ellipsoids using calculation performed in "main_evolve_displacement_gradient_tensor()"

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
  error('main_calculate_and_plot_strain_foliation_S1_map(): Unprogrammed number of input arguments')
end  

t_start = tic;

file_source = fullfile(pwd, mfilename);
file_date   = datestr(now);

dir_ROOT    = pwd;
dir_ARC     = ModelDirName;
dir_SNAP    = fullfile(dir_ARC, 'Snapshots');
file_DAT    = fullfile(dir_ARC, 'DOWN_TRACKS.mat');

DIR         = dir(file_DAT);

fprintf(1,'\n1. Loading down-tracking data (%.2f GB) from "%s"\n', DIR.bytes/10^9, file_DAT);

load(file_DAT, 'MODEL', 'P_dn')

P_save = P_dn;

dt     = MODEL.dt;
ns     = numel(P_dn);
nx     = MODEL.nx;
ny     = MODEL.ny;
x_grid = MODEL.x_grid;
y_grid = MODEL.y_grid;

PLOT    = InitializePlotStructure(fullfile(dir_ARC, 'PlotList.dat'), MODEL);
MAP     = PLOT.Map.StrainFoliation;

PlotSubDir = fullfile('PLOTS', 'Maps', 'StrainFoliation');

n_PLOTS = numel(MAP);

yr_site = [P_dn.t_site];

if any(yr_site~=yr_site(1))
  error('main_calculate_and_plot_strain_ellipsoid_axes(): All strain sites must have same observation year')
else
  yr_site       = yr_site(1);
  yr_last_surge = P_dn(1).yr_last_surge;
end

DIR = dir(fullfile(dir_SNAP, '*.mat'));

ModelFileList = {DIR.name};
ModelYearList = ExtractYearListFromFileList(ModelFileList);

[D, ~, ~, ~, ~, YearSwap] = GetAliasYearData(yr_site, yr_last_surge, dir_SNAP, ModelFileList, ModelYearList, MODEL, 0);

H            = zeros(ny,nx);
H(D.ic_jc_I) = D.H;

handle = 0;

for m=1:numel(MAP)
  
  P_dn     = P_save; 
  
  x_W      = min(x_grid);
  x_E      = max(x_grid);
  y_S      = min(y_grid);
  y_N      = max(y_grid);
  
  if isfield(MAP, 'x_W')
    x_W = max(x_W, MAP(m).x_W);
  end  
  if isfield(MAP, 'x_E')
    x_E = min(x_E, MAP(m).x_E);
  end
  if isfield(MAP, 'y_S')
    y_S = max(y_S, MAP(m).y_S);
  end
  if isfield(MAP, 'y_N')
    y_N = min(y_N, MAP(m).y_N);
  end

  L_x_use = x_grid>=x_W & x_grid<=x_E;
  L_y_use = y_grid>=y_S & y_grid<=y_N;  
  
  L_map = false(1, ns);
  
  if ~isfield(MAP(m), 'XI')
    MAP(m).XI = [];
  end
  
  if isfield(MAP(m), 'Filter') && ~isempty(MAP(m).Filter)
    L_map = false(1, ns);
    t_min = NaN(1, ns);
    
    for s=1:ns
      t_min(s) = P_dn(s).t(1);
      ng       = numel(P_dn(s).Group);
      for gg=1:ng
        if strcmp(P_dn(s).Group(gg).Text, MAP(m).Filter)
          L_map(s) = 1;
        end
      end  
    end
    if ~any(L_map)
      fprintf(1,'main_calculate_and_plot_strain_foliation_S1_map(): No group points found for n=%d with search text "%s"', m, MAP(m).Filter)
      continue
    end
        
    P_dn  = P_dn(L_map);
    ns    = numel(P_dn);
  end
  
  if ~isempty(MAP(m).XI)
    L_xi = abs([P_dn.XI_site] - MAP(m).XI)<0.0001;
    P_dn = P_dn(L_xi);
    ns   = numel(P_dn);
  end
  
  handle = handle+1;
  figure(handle)
  
  t_list = yr_site;

  R_factor     = MAP(m).R_factor;
  Ratio        = MAP(m).Ratio;
  TextSize     = MAP(m).TextSize;
  c_map        = colormap('jet');
  c_map(1,:)   = MAP(m).Gray;
  c_map(end,:) = [1 1 1];
    
  colormap(c_map)
   
  imagesc(x_grid(L_x_use), y_grid(L_y_use), H(L_y_use,L_x_use)>0), axis equal, axis tight
  set(gca, 'YDir', 'normal')
  
  if isfield(MAP(m), 'x_W') && ~isempty(MAP(m).x_W)
    xlim([MAP(m).x_W MAP(m).x_E])
  end
  
  if isfield(MAP(m), 'y_S') && ~isempty(MAP(m).y_S)
    ylim([MAP(m).y_S MAP(m).y_N])
  end
   
  hold on

  for s=1:ns
    x_IN = P_dn(s).X(1);
    y_IN = P_dn(s).Y(1);
    t_IN = P_dn(s).t(1);
  
    x_SITE = P_dn(s).X_site;
    y_SITE = P_dn(s).Y_site;
    
    if isfield(MAP, 'PlotSources') && MAP(m).PlotSources==1
      p1 = plot(x_IN, y_IN, 'ko', 'MarkerSize', 0.75*MAP(m).SiteMarkerSize);
      p2 = plot(x_IN, y_IN, 'k+', 'MarkerSize', MAP(m).SiteMarkerSize);
    end
    
    t_trk = P_dn(s).t;
    
    L_use = false(1, numel(t_trk));
    L_use(end) = 1;
    StrainFoliationPlane(P_dn(s), L_use, MAP)
    
    plot(x_SITE, y_SITE, 'k.')
  end
  hold off
  
  xlabel(MAP(m).xlabel)
  ylabel(MAP(m).ylabel)
  title(horzcat({sprintf('MODEL %s', MODEL.name), MAP(m).title}), 'interpreter', 'none')
  
  if MAP(m).HARDCOPY
    PlotDir = fullfile(dir_ARC, PlotSubDir);
    if ~exist(PlotDir, 'dir')
      mkdir(PlotDir)
    end
    print(handle, fullfile(PlotDir, sprintf('StrainFoliation(%d).pdf', handle)), '-dpdf')
  end  
end

fprintf(1,'\nALL DONE\n')




