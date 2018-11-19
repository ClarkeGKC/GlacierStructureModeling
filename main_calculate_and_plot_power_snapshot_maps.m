function main_calculate_and_plot_power_snapshot_maps(ModelName)

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
  fprintf(1,'\n');
else
  clear global
  clearvars -except ModelName
  ModelDirName = fullfile(pwd, 'Archive-%s', ModelName);
end 

global X2 Y2

t_start = tic;

VERBOSE      = 0;

dir_ROOT   = fullfile(pwd, sprintf('Archive-%s', ModelName));
dir_SNAP   = fullfile(dir_ROOT, 'Snapshots');
dir_PLOT   = fullfile(dir_ROOT, 'PLOTS', 'Maps', 'Power');

DIR = dir(fullfile(dir_SNAP, 'Out(*).mat'));

ModelFileList     = {DIR.name};
ModelFileYearList = ExtractYearListFromFileList(ModelFileList);

SiteFileName = fullfile(dir_ROOT, 'SiteList.dat');
load(fullfile(dir_ROOT, 'Model.mat'), 'MODEL') 

RHO_w = MODEL.par.RHO_w;
RHO   = MODEL.par.RHO;
g     = MODEL.par.g;

nx    = MODEL.nx;
ny    = MODEL.ny;
N_xy  = MODEL.N_xy;

dx    = MODEL.dx;
dy    = MODEL.dy;

n_XI  = MODEL.n_XI;
XI    = MODEL.XI;

YRSEC = MODEL.par.YRSEC;

load(MODEL.file_dem, 'B');

B     = B-MODEL.z_0;

PLOT      = InitializePlotStructure(fullfile(dir_ROOT, 'PlotList.dat'), MODEL);
PLOT      = PLOT.Map.Power;

if PLOT.HARDCOPY
  if ~exist(dir_PLOT, 'dir')
    mkdir(dir_PLOT)
  end
end

[~, ~, t_last_surge,~, ~, ~, ~, ~, ~, ~, ~, ~] = GetSiteList(SiteFileName, MODEL);

PlotYearList = PLOT.PlotYearList;
nt       = numel(PlotYearList);

fprintf(1,'1. Extracting snapshots and constructing vertically-integrated power maps for each time slice\n');

PHI_map = NaN(nt, ny, nx);
H_map   = NaN(nt, ny, nx);
L_H     = false(nt, ny, nx);

for t=1:nt
  
  YearGet   = PlotYearList(t);
  
  [DD, VV, LL, EE, TT, YearSwap] = GetAliasYearData(YearGet, t_last_surge, dir_SNAP, ModelFileList, ModelFileYearList, MODEL);
  
  H = zeros(MODEL.ny, MODEL.nx);
  ic_jc_I    = DD.ic_jc_I;
  H(ic_jc_I) = DD.H;
  L_H(t,:,:)     = H>0;
  
  H_map(t,:,:) = H;
  
  % A_glacier(t) = dx*dy*sum(sum(L_I));
  
  s   = zeros(3,3,n_XI, ny, nx);
  L   = zeros(3,3,n_XI, ny, nx);
  
  s_I = EE.s;
  s(:,:,:,ic_jc_I) = s_I;
  L_I = LL.L;
  L(:,:,:,ic_jc_I) = L_I;
  
  % convert L to 1/s 
  s_dot_L = s.*L/YRSEC;
  PHI     = reshape(s_dot_L, 9, n_XI, ny, nx);        % Rate of creep/frictional energy dissipation per unit volume W/m^3
  PHI     = squeeze(sum(PHI, 1));
  
  PHI        = reshape(PHI, n_XI, N_xy);
  PHI_column = reshape(H, 1, N_xy).*trapz(XI, PHI);   % Rate of energy dissipation per unit area of bed (W/m^2)
  
  PHI_cell       = dx*dy*reshape(PHI_column, ny, nx);     % Rate of frictional energy dissipation in a grid cell     [W]
  
  PHI_mW_per_m2  = 1000*reshape(PHI_column, ny, nx);
  
  PHI_map(t, L_H(t,:,:)) = PHI_mW_per_m2(L_H(t,:,:));

end

% Fiddle the color map to make gray background for no-data/NaN regions

PHI_map_max = PLOT.PhiMax_mW;

nc    = 256; 
d_PHI = PHI_map_max/(nc-1); 

c_map = colormap(jet(nc-1));
c_map = [PLOT.Gray; c_map];

CLIM = [-d_PHI PHI_map_max];

for t=1:nt
  if t==1
    x_grid  = MODEL.x_grid;
    y_grid  = MODEL.y_grid;
      
    if isfield(PLOT, 'x_W') && ~isempty(PLOT.x_W)
      x_MIN = PLOT.x_W;
    else
      x_MIN = x_grid(1);
    end
    if isfield(PLOT, 'x_E') && ~isempty(PLOT.x_E)
      x_MAX = PLOT.x_E;
    else
      x_MAX = x_grid(end);
    end
    if isfield(PLOT, 'y_S') && ~isempty(PLOT.y_S)
      y_MIN = PLOT.y_S;
    else
      y_MIN = y_grid(end);
    end
    if isfield(PLOT, 'y_N') && ~isempty(PLOT.y_N)
      y_MAX = PLOT.y_N;
    else
      y_MAX = y_grid(1);
    end  
    
    Lx_grid = x_grid>=x_MIN & x_grid<=x_MAX;
    Ly_grid = y_grid>=y_MIN & y_grid<=y_MAX;
    
    x_grid  = x_grid(Lx_grid);
    y_grid  = y_grid(Ly_grid);
   
    x_MIN = x_grid(1)-0.5*MODEL.dx;
    x_MAX = x_grid(end)+0.5*MODEL.dx;
    y_MIN = y_grid(end)-0.5*MODEL.dy;
    y_MAX = y_grid(1)+0.5*MODEL.dx;
  end    

  figure(t)
  
  plot([x_MIN x_MAX x_MAX x_MIN x_MIN], [y_MIN y_MIN y_MAX y_MAX y_MIN], 'k')
  
  hold on
  
  PHI_plot = squeeze(PHI_map(t,:,:));
  
  PHI_plot(isnan(PHI_plot)) = -d_PHI;
  
  colormap(c_map)
  
  PHI_plot = PHI_plot( Ly_grid, Lx_grid);

  imagesc(x_grid, y_grid, PHI_plot, CLIM), colorbar

  axis equal
  axis tight
  xlabel(PLOT.xlabel)
  ylabel(PLOT.ylabel);
  
  if isfield(PLOT, 'title') && ~isempty(PLOT.title)
    title({sprintf('%s "%s"',PLOT.title, ModelName), 'Vertically-integrated power density (mW/m^2)'}, 'interpreter', 'none')
  end
  
  if isfield(PLOT, 'YearPlot') && ~isempty(PLOT.YearPlot)
    text(x_grid(end)-450, y_grid(1)-75, sprintf('%.1f', PlotYearList(t)), 'FontSize', 11, 'FontWeight', 'bold')
  end
  hold off
  
  if PLOT.HARDCOPY
    print(t,fullfile(dir_PLOT, sprintf('Fig_P(%.2f).pdf', PlotYearList(t))), '-dpdf')
  end
    
end

fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);

