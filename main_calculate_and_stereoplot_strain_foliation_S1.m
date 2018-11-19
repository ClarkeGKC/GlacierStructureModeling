function main_calculate_and_stereoplot_strain_foliation_S1(ModelName)

% This script plot the strain ellipsoids using calculation performed in "main_evolve_displacement_gradient_tensor()"

% Note that the "DIRECT" calculation method is truly straightforward and
% avoids some subtle pitfalls associated with the "STRIKE-DIP" method

% Furthermore: The "STRIKE-DIP" method calculations have some subtle and uncorrected errors
% Thus       : Do not use "STRIKE-DIP" option

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
  error('main_calculate_and_plot_strain_ellipsoid_axes(): Unprogrammed number of input arguments')
end

METHOD = 'STRIKE-DIP';
METHOD = 'DIRECT';   

t_start = tic;

file_source = fullfile(pwd, mfilename);
file_date   = datestr(now);

dir_ARC     = ModelDirName;
dir_SNAP    = fullfile(dir_ARC, 'Snapshots');
dir_PLOT    = fullfile(dir_ARC, 'PLOTS', 'StrainFoliation_S1_Stereoplots');
file_DAT    = fullfile(dir_ARC, 'DOWN_TRACKS.mat');

DIR         = dir(file_DAT);

file_model = fullfile(pwd, sprintf('Archive-%s', ModelName), 'Model.mat');

load(file_model, 'MODEL')

x_grid = MODEL.x_grid;
y_grid = MODEL.y_grid;

SiteFileName  = fullfile(pwd, sprintf('Archive-%s', ModelName), 'SiteList.dat');
  
fprintf(1,'\n1. Loading site data (not necessarily up-to-date) in down-tracking data file\n');

[~, ~, t_last_surge, ~, ~, ~, ~, ~, ~, Name_site, Type_site, Group] = GetSiteList(SiteFileName, MODEL);

GroupText = {Group.Text};
GroupID   = {Group.ID};

OK = 0;

OPTIONS  = InitializePlotStructure(fullfile(dir_ARC, 'PlotList.dat'));
MAP      = OPTIONS.Map.StrainFoliationSites;

clear OPTIONS

fprintf(1,'\n2. Loading down-tracking data (%.2f GB) from "%s"\n', DIR.bytes/10^9, file_DAT);

load(file_DAT, 'MODEL', 'P_dn')

NS     = numel(Type_site);

if NS~=numel(P_dn)
  fprintf(1,'\nmain_calculate_and_stereoplot_strain_foliation(): WARNING >> Number of sites in "SiteList.dat" differs from size of P_dn structure\n\n');
  fprintf(1,'  Number of sites in "SiteList.dat"      = %d\n', NS);
  fprintf(1,'  Number of elements of structure P_dn   = %d\n\n', numel(P_dn)); 
end

yr_site = [P_dn.t_site];

if any(yr_site~=yr_site(1))
  error('main_calculate_and_stereoplot_strain_foliation_S1(): All strain sites must have same observation year')
end

yr_site  = yr_site(1);

% OPTIONS for Diptych_S1 plots

OPTIONS.PlotTitle       = 0;

OPTIONS.HardCopy        = 1;   % Stereoplots HARDCOPY only;
OPTIONS.AzGrid          = 0;
OPTIONS.AzTicks         = 1;
OPTIONS.AzTickInterval  = 90;
OPTIONS.AzTickRfract    = 0.07;
OPTIONS.AzLabelRfract   = 0.25;
OPTIONS.AzLabels        = 'alphabet';
OPTIONS.AzLabelFontScale = 1.0;
OPTIONS.VelArrow        = 1;
OPTIONS.TitleFontScale  = 1.20;
OPTIONS.TextFontScale   = 1.0;
OPTIONS.MarkerScale     = 1.0;
OPTIONS.MinDip          = 30;
OPTIONS.MinAxisRatio    = 5;

OPTIONS.Color_1         = MAP.Color_1;
OPTIONS.Color_2         = MAP.Color_2;
OPTIONS.Color_3         = MAP.Color_3;
OPTIONS.Color_4         = MAP.Color_4;

handle = 0;

if ~isfield(MAP, 'XI')
  MAP.XI = [];
end

if isfield(MAP, 'Filter') && ~isempty(MAP.Filter)
  L_zone = false(1, NS);

  for s=1:NS
    if ~isempty(findstr(GroupCode, Type_site{s}))
      L_zone(s) = 1;
    end
  end
  P_dn  = P_dn(L_zone);
  ns    = numel(P_dn);
end  
  
if ~isempty(MAP.XI)
  L_xi = abs([P_dn.XI_site] - MAP.XI)<0.0001;
  P_dn = P_dn(L_xi);
  ns   = numel(P_dn);
end
 

X_site = zeros(1, ns);
Y_site = zeros(1, ns);
    
X_c  = zeros(1, ns);
Y_c  = zeros(1, ns);
R    = zeros(3, 3, ns);
E_U  = zeros(3, ns);
W_U  = zeros(3, 3, ns);
  
for s=1:ns       
  t_trk = P_dn(s).t;
    
  L_use = false(1, numel(t_trk));
  L_use(end) = 1;
  
  Name_site{s} = P_dn(s).Name_site;
  
  X_site(s)   = P_dn(s).X_site;
  Y_site(s)   = P_dn(s).Y_site;
  
  X_c(s)      = P_dn(s).X(end);
  Y_c(s)      = P_dn(s).Y(end);
  R(:,:,s)    = P_dn(s).R_pre(:,:,end);
  E_U(:,s)    = P_dn(s).E_U(:,end);
  W_U(:,:,s)  = P_dn(s).vecs_U(:,:,end);
end
 
[x_p, y_p, z_p, phi, lambda, Ratio_23] = StrainFoliationStereoplotPoints(X_c, Y_c, R, E_U, W_U, METHOD);
 
N_pts = numel(x_p);

SITE.Name_site = 'Model';

Color_p        = ones(size(x_p));
L_use          = true(size(x_p));

L_angle        = phi<pi*OPTIONS.MinDip/180;
L_ratio        = Ratio_23<OPTIONS.MinAxisRatio;
L_both         = L_angle&L_ratio;

Color_p(L_angle) = 2;
Color_p(L_ratio) = 3;
Color_p(L_both)  = 4;
 
handle = Diptych_S1(handle, MODEL, SITE, N_pts, x_p, y_p, z_p, L_use, Color_p,  [], 'S_1 foliation', OPTIONS, 0);
SaveEigenvalueTable(dir_PLOT, SITE, 'Modelled S_1', x_p, y_p, z_p, [], [], 0);

if OPTIONS.HardCopy
  if ~exist(dir_PLOT, 'dir')
    mkdir(dir_PLOT)
  end
  print(handle, fullfile(dir_PLOT, sprintf('StrainFoliationStereoPlot(%s-ALL).pdf', METHOD)), '-dpdf')
end

L_use   = ~L_angle&~L_ratio;

N_use = sum(L_use);
SITE.Name_site = 'Edited';

handle = Diptych_S1(handle, MODEL, SITE, N_use, x_p, y_p, z_p, L_use, ones(size(x_p)),  [], 'S_1 foliation (edited)', OPTIONS, 0);
SaveEigenvalueTable(dir_PLOT, SITE, 'Edited S_1', x_p(L_use), y_p(L_use), z_p(L_use), [], [], 0);

if OPTIONS.HardCopy
  if ~exist(dir_PLOT, 'dir')
    mkdir(dir_PLOT)
  end
  print(handle, fullfile(dir_PLOT, sprintf('StrainFoliationStereoPlot(%s-EDITED).pdf', METHOD)), '-dpdf')
end

% Plot measurement site locations (for simulation)

DIR = dir(fullfile(dir_SNAP, 'Out(*).mat'));

ModelFileList     = {DIR.name};
ModelFileYearList = ExtractYearListFromFileList(ModelFileList);

YearGet   = yr_site;

[D, ~, ~, ~, ~, YearSwap] = GetAliasYearData(YearGet, t_last_surge, dir_SNAP, ModelFileList, ModelFileYearList, MODEL);
  
H = zeros(MODEL.ny, MODEL.nx);
H(D.ic_jc_I) = D.H;

if isfield(MAP, 'x_W') && ~isempty(MAP.x_W)
  x_MIN = MAP.x_W;
else
  x_MIN = x_grid(1);
end

if isfield(MAP, 'x_E') && ~isempty(MAP.x_E)
  x_MAX = MAP.x_E;
else
  x_MAX = x_grid(end);
end

if isfield(MAP, 'y_S') && ~isempty(MAP.y_S)
  y_MIN = MAP.y_S;
else
  y_MIN = y_grid(end);
end

if isfield(MAP, 'y_N') && ~isempty(MAP.y_N)
  y_MAX = MAP.y_N;
else
  y_MAX = y_grid(1);
end  
  
Lx_grid = x_grid>=x_MIN & x_grid<=x_MAX;
Ly_grid = y_grid>=y_MIN & y_grid<=y_MAX;

H  = H(Ly_grid, Lx_grid);

x_grid = x_grid(Lx_grid);
y_grid = y_grid(Ly_grid);

x_MIN = x_grid(1) - 0.5*MODEL.dx;
x_MAX = x_grid(end) + 0.5*MODEL.dx;
y_MIN = y_grid(end)-0.5*MODEL.dy;
y_MAX = y_grid(1)+0.5*MODEL.dy;

I  = H>0;

handle = handle+1;
  
figure(handle)
plot([x_MIN x_MAX x_MAX x_MIN x_MIN], [y_MIN y_MIN y_MAX y_MAX y_MIN], 'k')
hold on
  
colormap('jet')
c_map = colormap;
if isfield(MAP, 'Gray') && ~isempty(MAP.Gray)
  c_map(1,:)   = MAP.Gray;
  c_map(end,:) = [1 1 1];
else  
  c_map(1,:)   = [0.85 0.85 0.85];
  c_map(end,:) = [1 1 1];
end

colormap(c_map)
imagesc(x_grid, y_grid, double(I)), axis equal, axis tight
hold on

plot(X_site(Color_p==1), Y_site(Color_p==1), 'ko', 'MarkerFaceColor', OPTIONS.Color_1);
plot(X_site(Color_p==2), Y_site(Color_p==2), 'ko', 'MarkerFaceColor', OPTIONS.Color_2);
plot(X_site(Color_p==3), Y_site(Color_p==3), 'ko', 'MarkerFaceColor', OPTIONS.Color_3);
plot(X_site(Color_p==4), Y_site(Color_p==4), 'ko', 'MarkerFaceColor', OPTIONS.Color_4);

if isfield(MAP, 'x_label') && ~isempty(MAP.x_label)
  xlabel(MAP.x_label)
end
if isfield(MAP, 'y_label') && ~isempty(MAP.y_label)
  ylabel(MAP.y_label)
end  

if isfield(MAP, 'title') && ~isempty(MAP.title)
  title(MAP.title, 'Interpreter', 'none')
end  
    
if isfield(MAP, 'HARDCOPY') && ~isempty(MAP.HARDCOPY)
  print(handle, fullfile(dir_PLOT, sprintf('S_1 sites(%s).pdf', METHOD)), '-dpdf')
end   
 
str = {};
str{end+1} = ' ';
str{end+1} = sprintf('TABLE OF S1 POINTS : %s METHOD', METHOD);
str{end+1} = sprintf('MODEL "%s"', ModelName);
str{end+1} = ' ';
str{end+1} = 'Site         PHI    LAMBDA     x_p      y_p     z_p';
str{end+1} = '            (deg)   (deg)';
str{end+1} = ' ';

switch METHOD
  case 'STRIKE-DIP'
    for s=1:N_pts
      str{end+1} = sprintf('%-10s  %5.2f   %6.2f   %7.4f  %7.4f  %7.4f', char(Name_site{s}), 180*phi(s)/pi, 180*lambda(s)/pi, x_p(s), y_p(s), z_p(s));
    end
  case 'DIRECT'
    for s=1:N_pts
      str{end+1} = sprintf('%-10s  %5.2f  -------  %7.4f  %7.4f  %7.4f', char(Name_site{s}), 180*phi(s)/pi,  x_p(s), y_p(s), z_p(s));
    end
  otherwise
    error('main_calculate_and_stereoplot_strain_foliation_S1(): Unprogrammed case ')
end

str{end+1} = '';

f_out = fopen(fullfile(dir_PLOT, sprintf('SiteTable(%s).dat', METHOD)), 'w');

for s=1:numel(str)
  fprintf(1,'%s\n', str{s})
  fprintf(f_out, '%s\n', str{s});
end

fclose(f_out);

fprintf(1,'\nALL DONE\n')
