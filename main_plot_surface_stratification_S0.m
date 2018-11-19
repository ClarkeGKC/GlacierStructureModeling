function main_plot_surface_stratification_S0(ModelName)

% This script aims to plot the near-horizontal stratification exposed at the glacier surface
% Hambrey identifies this as the primary stratification S_0 and I assume it is an isochronal surface

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
elseif nargin==1
  clear global
  clearvars -except ModelName
  ModelDirName = fullfile(pwd, sprintf('Archive-%s', ModelName));
else
  error('main_plot_longitudinal_foliation(): Unprogrammed number of input arguments')
end  

t_start = tic;

VERBOSE    = 0;

dir_ARC    = fullfile(pwd, sprintf('Archive-%s', ModelName));
dir_SNAP   = fullfile(dir_ARC, 'Snapshots');

dir_OUT  = dir_ARC;
dir_PLOT = fullfile(dir_ARC, 'PLOTS', 'Maps', 'Primary Stratification_S0');

if ~exist(dir_PLOT, 'dir')
  mkdir(dir_PLOT);
end  

file_model = fullfile(pwd, sprintf('Archive-%s', ModelName), 'Model.mat');

load(file_model, 'MODEL')

SetupIndexArrays(MODEL)

PLOT    = InitializePlotStructure(fullfile(dir_ARC, 'PlotList.dat'), MODEL);
MAP     = PLOT.Map.t_d;

SiteFileName = fullfile(dir_ARC, 'SiteList.dat');

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
XI        = MODEL.XI;

dx        = MODEL.dx;
dy        = MODEL.dy;
x_grid    = MODEL.x_grid;
y_grid    = MODEL.y_grid;

dt        = MODEL.dt;
t_CYCLE   = MODEL.t_CYCLE;

% Minimum and maximum tracer dates (in model years)

yr_STOP   = MODEL.t_STOP;
yr_START  = MODEL.t_tracer_START;

N_xy      = nx*ny;
N_all     = N_xy*n_XI;

[Y3, Z3, X3] = meshgrid(y_grid, XI, x_grid);

load(MODEL.file_dem, 'B') 

[SiteGrid, t_site, t_last_surge, X_NAD27_site, Y_NAD27_site, X_WGS84_site, Y_WGS84_site, XI_site, Name_site, Type_site, Group] = GetSiteList(SiteFileName, MODEL);

X_site = X_NAD27_site;
Y_site = Y_NAD27_site;

% Perform grid shift

B      =  B-z_0;
X_site = X_site - x_0;
Y_site = Y_site - y_0;

DIR = dir(fullfile(dir_SNAP, '*.mat'));

ModelFileList = {DIR.name};
ModelYearList = ExtractYearListFromFileList(ModelFileList);

IPREC    = GetTimeStepPrecision(dt);

[FileSwap, YR_swap] = GetYearSwapFileName(MODEL, t_last_surge, ModelFileList, ModelYearList, t_site);   % Use alias years to match stop year with an output file

td_shift_yr = t_site-YR_swap;
  
load(fullfile(dir_SNAP, FileSwap), 'DD', 'V', 'LL', 'EE')

H       = zeros(ny, nx);
H(DD.ic_jc_I) = DD.H;
L_I     = H>0;

S       = B+H;

xd      = X3;
yd      = Y3;
td      = t_site*ones(size(X3));

xd(:,L_I) = DD.xd;
yd(:,L_I) = DD.yd;
td(:,L_I) = DD.td;

td      = td+td_shift_yr;

td_plot = squeeze(td(n_XI, :, :));
xd_plot = squeeze(xd(n_XI, :, :));
yd_plot = squeeze(yd(n_XI, :, :));

fprintf(1,'1. Using same "Lower Trapridge" framing as for "Map of lower part of Trapridge Glacier\n\n')

x_W = min(x_grid);
x_E = max(x_grid);
y_S = min(y_grid);
y_N = max(y_grid);

if isfield(MAP, 'x_W')
  x_W = max(x_W, MAP.x_W);
end  
if isfield(MAP, 'x_E')
  x_E = min(x_E, MAP.x_E);
end
if isfield(MAP, 'y_S')
  y_S = max(y_S, MAP.y_S);
end
if isfield(MAP, 'y_N')
  y_N = min(y_N, MAP.y_N);
end

L_x_use = x_grid>=x_W & x_grid<=x_E;
L_y_use = y_grid>=y_S & y_grid<=y_N;

handle = 1;
figure(handle)

if isfield(MAP, 'ShadedBackground')
  C_gray = MAP.ShadedBackground;
  C_map  = colormap('gray');
  C_map(end,:) = [1 1 1];
  C_map(1, :)  = C_gray;
  colormap(C_map)
  B_plot       = real(H(L_y_use, L_x_use)>0);
  imagesc(x_grid(L_x_use), y_grid(L_y_use), B_plot), axis equal, axis tight
  set(gca, 'YDir', 'normal')
end

C       = contourc(x_grid(L_x_use), y_grid(L_y_use), td_plot(L_y_use, L_x_use), MAP.ContourLevels);

if isfield(MAP, 'ThickLineLevels')
  C_thick = contourc(x_grid(L_x_use), y_grid(L_y_use), td_plot(L_y_use, L_x_use), MAP.ThickLineLevels);
else
  C_thick = [];
end

axis equal
axis tight

hold on

LineSpec = 'k';

if isfield(MAP, 'ContourColors')
  if MAP.ContourColors
    LineSpec = MAP.ContourColors';
  end  
end

if isfield(MAP, 'DefaultLineWidth')
  DefaultLineWidth = MAP.DefaultLineWidth;
else
  DefaultLineWidth = 1;
end  

i = 1;

while i<size(C,2)
  n = C(2,i);
  L = ones(1,n)*C(1,i);
  X = C(1, i+1:i+n);
  Y = C(2, i+1:i+n);
  plot(X, Y, LineSpec, 'LineWidth', DefaultLineWidth);
  i     = i+n+1;
end

if ~isempty(C_thick)
  if isfield(MAP, 'ThickLineWidth')
    ThickLineWidth = MAP.ThickLineWidth;
  else
    ThickLineWidth = 5;
  end
  
  i = 1;
  while i<size(C_thick,2)
    n = C_thick(2,i);
    L = ones(1,n)*C_thick(1,i);
    X = C_thick(1, i+1:i+n);
    Y = C_thick(2, i+1:i+n);  
    plot(X, Y, 'k', 'LineWidth', ThickLineWidth)
    i    = i+n+1;
  end
end  

if MAP.PlotSites 
  plot(X_site, Y_site, 'wo', X_site, Y_site, 'k+');
end

hold off

title({sprintf('Lower "%s" at %.1f CE', ModelName, t_site), MAP.title}, 'interpreter', 'none')
xlabel(MAP.xlabel)
ylabel(MAP.ylabel)

if MAP.HARDCOPY
  print(handle, fullfile(dir_PLOT, 'Fig_S0.pdf'), '-dpdf')
end

handle = handle+1;
figure(handle)

j_list = 1:ny;
j_use  = j_list(L_y_use);
j_mid  = ceil(mean(j_use));

td_profile = squeeze(td(:,j_mid, L_x_use));

contour(td_profile, MAP.ContourLevels(1:5:end), 'ShowText', 'on')

xlabel('Distance (m)')
title({sprintf('Lower "%s" at %.1f CE', ModelName, t_site), 'Profile of depositational year'}, 'interpreter', 'none')

if MAP.HARDCOPY
  print(handle, fullfile(dir_PLOT, 'Fig_td_profile.pdf'), '-dpdf')
end

handle = handle+1;
figure(handle)
contour(x_grid(L_x_use), y_grid(L_y_use), H(L_y_use, L_x_use), 51), axis equal, axis tight
xlabel('Easting (m)')
ylabel('Northing (m)')
title({sprintf('Lower "%s" at%.1f CE', ModelName, t_site), 'Ice thickness (m)'}, 'interpreter', 'none')

if MAP.HARDCOPY
  print(handle, fullfile(dir_PLOT, 'Fig_H_ice.pdf'), '-dpdf')
end

fprintf(1,'\nALL DONE\n');


