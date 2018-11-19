function main_calculate_and_plot_power_timeseries(ModelName)

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
dir_PLOT   = fullfile(dir_ROOT, 'PLOTS', 'PowerTimeSeries');

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
PLOT      = PLOT.TimeSeries.Power;

if PLOT.HARDCOPY
  if ~exist(dir_PLOT, 'dir')
    mkdir(dir_PLOT)
  end
end

[~, ~, t_last_surge,~, ~, ~, ~, ~, ~, ~, ~, ~] = GetSiteList(SiteFileName, MODEL);

% The non reframed grid is potentially required for some subsequent 2D interpolations
  
[X2, Y2] = meshgrid(MODEL.x_grid-MODEL.x_0, MODEL.y_grid-MODEL.y_0);

dt_PLOT  = PLOT.dt_PLOT;
dt_SEC   = YRSEC*dt_PLOT;
YearList = PLOT.t_MIN : dt_PLOT : PLOT.t_MAX;
nt       = numel(YearList);

A_glacier   = zeros(nt, 1);
V_glacier   = zeros(nt, 1);
PHI_glacier = zeros(nt, 1);
U_glacier   = zeros(nt, 1);

fprintf(1,'1. Extracting snapshots and constructing time series elements from them\n');

for t=1:nt
  
  YearGet   = YearList(t);
  
  [DD, VV, LL, EE, TT, YearSwap] = GetAliasYearData(YearGet, t_last_surge, dir_SNAP, ModelFileList, ModelFileYearList, MODEL);
  
  H = zeros(MODEL.ny, MODEL.nx);
  ic_jc_I    = DD.ic_jc_I;
  H(ic_jc_I) = DD.H;
  L_I        = H>0;
  
  A_glacier(t) = dx*dy*sum(sum(L_I));
  
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
  PHI_glacier(t) = sum(sum(PHI_cell));              % Rate of frictional energy dissipation for entire glacier [W]
  
  V_glacier(t)   = sum(sum(H))*dx*dy;
  m_glacier      = RHO*V_glacier;                   % glacier mass [kg]
  
  Z_ice          = repmat(reshape(B, 1, N_xy), n_XI, 1) + kron(XI',reshape(H, 1, N_xy));
  
  U_ice          = RHO*g*Z_ice;  % Gravitational potential energy of ice particle at elevation Z_ice   [kg/m^3 * m/s^2 * m] =  N/m^2 = Pa]
  
  U_col          = reshape(H, 1, N_xy).*trapz(XI, U_ice);  % The vertically-integrated gravitational potential  [Pa m]
  U_cell         = dx*dy*reshape(U_col, ny, nx);  
  U_glacier(t)   = sum(sum(U_cell));              % Gravitational potential energy of entire glacier   [Pa m^3 = J]
     
end

U_glacier = U_glacier';     
A_glacier = A_glacier';

dU_dt = [ U_glacier(2)-U_glacier(1)  0.5*(U_glacier(3:nt)-U_glacier(1:nt-2))   U_glacier(nt)-U_glacier(nt-1) ]/dt_SEC;   % W
dA_dt = [ A_glacier(2)-A_glacier(1)  0.5*(A_glacier(3:nt)-A_glacier(1:nt-2))   A_glacier(nt)-A_glacier(nt-1) ]/dt_PLOT;  % m^2/s

u_glacier   = U_glacier./V_glacier;    % [J/m^3]  
phi_glacier = PHI_glacier./V_glacier;  % [W/m^3]

du_dt = [ u_glacier(2)-u_glacier(1)  0.5*(u_glacier(3:nt)-u_glacier(1:nt-2))   u_glacier(nt)-u_glacier(nt-1) ]/dt_SEC;  % W/m^3

fprintf(1,'\n2. Done with integrations and time series construction\n')

YearSpan = YearList-YearList(1);
L_use    = YearSpan < MODEL.t_CYCLE-0.1*MODEL.dt;

str        = {};
str{end+1} = sprintf('Table 1. Results of energy analysis for Model "%s"', ModelName);
str{end+1} = ' ';
str{end+1} = sprintf('Time average rate of change of gravitational potential energy  = %12.4e (W/m^3)', mean(du_dt(L_use)));
str{end+1} = sprintf('Time average rate of change of frictional energy dissipation   = %12.4e (W/m^3)', mean(phi_glacier(L_use)));
str{end+1} = sprintf('Time average rate of change of ice area                        = %12.4e (m^2/yr)', mean(dA_dt(L_use)));

file_tbl = fullfile(dir_PLOT, 'Table.txt');
f_tbl    = fopen(file_tbl, 'w');

for l=1:numel(str)
  fprintf(1,'%s\n', str{l});
  fprintf(f_tbl, '%s\n', str{l});
end

figure(1)
h_line = [];

h_line(end+1) = line(YearList, dU_dt, 'color', 'k', 'linestyle', '--'); hold on
h_line(end+1) = line(YearList, PHI_glacier, 'color', 'k', 'linestyle', '-');
ax_L = gca;
set(ax_L, 'Xcolor', 'k', 'Ycolor', 'k')
ax_R = axes('Position', get(ax_L, 'Position'), 'XaxisLocation', 'bottom', 'YAxisLocation', 'right', 'Color', 'none', 'Xcolor', 'none', 'Ycolor', 'b');
h_line(end+1) = line(YearList, dA_dt, 'color', 'b', 'Parent', ax_R);
hold off
xlabel(PLOT.xlabel)
ylabel(ax_L, PLOT.ylabel_L1)
ylabel(ax_R, PLOT.ylabel_R1)
xlim(ax_R, [PLOT.t_MIN PLOT.t_MAX])
xlim(ax_L, [PLOT.t_MIN PLOT.t_MAX])

title(sprintf('%s "%s"', PLOT.title, ModelName), 'interpreter', 'none')

if PLOT.Legend
  legend(h_line, {'dU/dt (W)', '\Phi(t) (W)', 'dA/dt (m^2/yr)'}, 'location', 'northwest');
end 
if PLOT.HARDCOPY
  print(1,fullfile(dir_PLOT, 'Fig_01.pdf'), '-dpdf')
end

figure(2)
h_line = [];  

du_dt_mW_m3       = 1000*du_dt;
phi_glacier_mW_m3 = 1000*phi_glacier;

h_line(end+1) = line(YearList, du_dt_mW_m3, 'color', 'k', 'linestyle', '--'); hold on
h_line(end+1) = line(YearList, phi_glacier_mW_m3, 'color', 'k', 'linestyle', '-');
h_line(end+1) = line([YearList(1) YearList(nt)], [mean(du_dt_mW_m3(L_use)) mean(du_dt_mW_m3(L_use))], 'color', 'k', 'linestyle', '-.'); hold on
h_line(end+1) = line([YearList(1) YearList(nt)], [mean(phi_glacier_mW_m3(L_use)) mean(phi_glacier_mW_m3(L_use))], 'color', 'k', 'linestyle', ':'); hold on

ax_L = gca;
set(ax_L, 'Xcolor', 'k', 'Ycolor', 'k')
ax_R = axes('Position', get(ax_L, 'Position'), 'XaxisLocation', 'bottom', 'YAxisLocation', 'right', 'Color', 'none', 'Xcolor', 'none', 'Ycolor', 'b');
h_line(end+1) = line(YearList, dA_dt, 'color', 'b', 'Parent', ax_R);
h_line(end+1) = line([YearList(1) YearList(nt)], [mean(dA_dt(L_use)) mean(dA_dt(L_use))], 'color', 'b', 'linestyle', 'none', 'marker', '.'); 
hold off
xlabel(PLOT.xlabel)
ylabel(ax_L, PLOT.ylabel_L2)
ylabel(ax_R, PLOT.ylabel_R2)
xlim(ax_R, [PLOT.t_MIN PLOT.t_MAX])
xlim(ax_L, [PLOT.t_MIN PLOT.t_MAX])

title(sprintf('%s "%s"', PLOT.title, ModelName), 'interpreter', 'none')

if PLOT.Legend
  legend(h_line, {'du/dt (mW/m^3)', '\phi (mW/m^3)', '<du/dt>', '<\phi>', 'dA/dt (m^2/yr)', '<dA/dt>'}, 'location', 'northwest');
  %legend('boxoff')
end 
if PLOT.HARDCOPY
  print(2,fullfile(dir_PLOT, 'Fig_02.pdf'), '-dpdf')
end


fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);

