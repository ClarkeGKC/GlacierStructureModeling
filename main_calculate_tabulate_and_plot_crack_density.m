function main_calculate_tabulate_and_plot_crack_density(ModelName)

% Read enhanced track file, perform and even analysis and generate event stereoplots

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
elseif nargin==1
  clear global
  clearvars -except ModelName
  ModelDirName  = fullfile(pwd, sprintf('Archive-%s', ModelName));
end  

global  RHO RHO_w g VERBOSITY

VERBOSITY   = 1;

t_start     = tic;

file_source = fullfile(pwd, mfilename);
file_date   = datestr(now);

dir_ARC  = ModelDirName;
file_DAT = fullfile(dir_ARC, 'DOWN_TRACKS.mat');

DIR      = dir(file_DAT);

fprintf(1,'\n1. Loading down-tracking data (%.2f GB) from "%s"\n', DIR.bytes/10^9, file_DAT);

load(file_DAT, 'MODEL', 'P_dn')

RHO    = MODEL.par.RHO;
RHO_w  = MODEL.par.RHO_w;
g      = MODEL.par.g;

PLOT    = InitializePlotStructure(fullfile(dir_ARC, 'PlotList.dat'), MODEL);
DENSITY = PLOT.Density.Crack;

TargetGroup = DENSITY.Filter;  
L_match     = true(1, numel(P_dn));

for s=1:numel(P_dn)
  GroupText = P_dn(s).Group.Text;
  L_match(s) = strcmp(TargetGroup, GroupText);
end

P_dn  = P_dn(L_match);

if isfield(DENSITY, 'NoPlotList') && ~isempty(DENSITY.NoPlotList)
  P_dn = SitePurge(P_dn, DENSITY.NoPlotList);
end

if numel(DENSITY)~=1
  error('main_calcululate_tabulate_and_plot_crack_density(): Present version of this code allows only a single density plot per run')
end  

ns      = numel(P_dn);

CRACK_MODEL = GetCrackModelParameters(fullfile(ModelDirName, 'CrackModel.dat'));

RULE        = CRACK_MODEL.RULE;
RuleStr     = CRACK_MODEL.RuleStr;
ParamStr    = CRACK_MODEL.ParamStr;

dir_PLOT    = fullfile(dir_ARC, 'PLOTS', 'CrackDensity', sprintf('Rule_#%02d', RULE));

if ~exist(dir_PLOT, 'dir')
  mkdir(dir_PLOT);
end  

BLANKING_RULE  = CRACK_MODEL.BLANKING_RULE;
BLANKING_MODEL = GetBlankingModelParameters(BLANKING_RULE, fullfile(ModelDirName, 'BlankingModel.dat'));
BlankStr       = BLANKING_MODEL.ParamStr;

% Define azimuthal bins in terms of dot product values

switch CRACK_MODEL.BIN_DEF
  case 'Hambrey'
    theta = pi*[20 90-20]/180;  % Consistent with Hambrey's transverse and longitudinal bins 
  case 'Equal'
    theta = pi*[22.5 90-22.5]/180;  % All bins have equal angle span
  otherwise
    error('main_calculate_tabulate_and_plot_crack_density(): Unprogrammed CRACK_BIN_DEF')
end

cos_theta = cos(theta);

% Eigenvalue table file grows by "appending". It should be destroyed before starting a new run

file_tbl = fullfile(dir_PLOT, sprintf('EigenvalueTable_%s.dat', ModelName));

if exist(file_tbl, 'file')
  delete(file_tbl)
end  

L_density  = false(1, ns);

TextSearch = DENSITY.Filter;

for s=1:ns
  t_min(s) = P_dn(s).t(1);
  ng       = numel(P_dn(s).Group);
  for gg=1:ng
    if strcmp(P_dn(s).Group(gg).Text, TextSearch)
      L_density(s) = 1;
    end
  end  
end

if sum(L_density)==0
  fprintf(1,'main_calculate_tabulate_and_plot_crack_density(): No group points found for p=%d with search text "%s"', p, TextSearch)
  error('main_calculate_tabulate_and_plot_crack_density(): Run aborted')
end

PP_dn  = P_dn(L_density);
t_min  = t_min(L_density);
  
N_pts   = zeros(numel(PP_dn),1);
N_long  = zeros(numel(PP_dn),1);
N_trans = zeros(numel(PP_dn),1);
N_diag  = zeros(numel(PP_dn),1);
N_sum   = zeros(numel(PP_dn),1);

PathTime = zeros(numel(PP_dn),1);
PathDist = zeros(numel(PP_dn),1);

for s=1:numel(PP_dn)
  nt     = numel(PP_dn(s).t);  
  t_pt   = PP_dn(s).t;
  
  % Eigenvectors D,  s, sigmma and R are identical (but not their eigenvalues)
     
  X_vecs = PP_dn(s).vecs_D;
      
  [L_crack, ~, ~] = ApplyCrackRule(CRACK_MODEL, PP_dn(s));
       
  L_use = L_crack;  
  t_use = t_pt(L_use);
  
  Ext_pts = squeeze(X_vecs(:,1,L_use));
  
  R_Ext   = P_dn(s).R_post(:,:,L_use);
  
  for k=1:sum(L_use)
    Ext_rot(:,k) = squeeze(R_Ext(:,:,k))*squeeze(Ext_pts(:,k));
    Ext_rot(:,k) = Ext_rot(:,k)/norm(Ext_rot(:,k));
  end
  
  % The data from V_p(s) are organized with time flowing from past to
  % present. Thus the final point corresponds to the sampling site and
  % NOT as previously assumed the 1st point
  
  n_flow_x = PP_dn(s).u(end)/sqrt(PP_dn(s).u(end)^2 + PP_dn(s).v(end)^2 + PP_dn(s).w(end)^2);
  n_flow_y = PP_dn(s).v(end)/sqrt(PP_dn(s).u(end)^2 + PP_dn(s).v(end)^2 + PP_dn(s).w(end)^2);
  n_flow_z = PP_dn(s).w(end)/sqrt(PP_dn(s).u(end)^2 + PP_dn(s).v(end)^2 + PP_dn(s).w(end)^2);
  
  n_flow_2d_x = n_flow_x/sqrt(n_flow_x^2 + n_flow_y^2);
  n_flow_2d_y = n_flow_y/sqrt(n_flow_x^2 + n_flow_y^2);
  
  % Note that that above might be better in 2D since Hambrey was not
  % concerned with vertical component of ice flow at the observation site
  % (and this could be quite large)
      
  N_pts(s) = sum(L_use);
  
  X_pt     = PP_dn(s).X;
  Y_pt     = PP_dn(s).Y;
  Z_pt     = PP_dn(s).Z;
  
  dX_pt    = X_pt(2:end) - X_pt(1:end-1);
  dY_pt    = Y_pt(2:end) - Y_pt(1:end-1);
  dZ_pt    = Z_pt(2:end) - Z_pt(1:end-1);
  
  ds_pt    = sqrt(dX_pt.^2 + dY_pt.^2 + dZ_pt.^2);
  
  PathTime(s) = PP_dn(s).t(end) - PP_dn(s).t(1); 
  PathDist(s) = sum(ds_pt);

  ds_xy_pt = [0 sqrt(dX_pt.^2 + dY_pt.^2)];
  s_xy_pt  = cumsum(ds_xy_pt);

  if sum(L_use)>0           
    Ext_rot = ProjectToLowerHemisphere(Ext_rot);
    t_p    = t_use;
    x_p    = Ext_rot(1,:);           % Coordinates (x_p, y_p, z_p) are for a unit vector orthogonal to a rotated fault plane
    y_p    = Ext_rot(2,:);
    z_p    = Ext_rot(3,:);
   
    s_xy_p = s_xy_pt(L_use);
   
    % Note that L_use winnows the x_p etc list (see
    
    if BLANKING_RULE>0
      [x_p, y_p, z_p, t_p] = ApplyBlanking(BLANKING_MODEL, x_p, y_p, z_p, t_p, s_xy_p);
    end  

    nc_x = x_p;
    nc_y = y_p;
    nc_z = z_p;
    
    % Convert to 2D (for comparison with Hambrey work)
        
    nc_2d_x = x_p./sqrt(x_p.^2 + y_p.^2);
    nc_2d_y = y_p./sqrt(x_p.^2 + y_p.^2);
    
    % The following could require more thought -- I think this is the best
    % comparison to Hambrey methodology
    
%   dot_prod   = nc_x*n_flow_x + nc_y*n_flow_y + nc_z*n_flow_z;
    dot_prod   = nc_2d_x*n_flow_2d_x + nc_2d_y*n_flow_2d_y;
    
    N_trans(s) = sum(abs(dot_prod)>cos_theta(1));
    N_long(s)  = sum(abs(dot_prod)<cos_theta(2));
    N_diag(s)  = sum(abs(dot_prod)<=cos_theta(1) & abs(dot_prod)>=cos_theta(2));
    N_sum(s)   = N_trans(s)+N_long(s)+N_diag(s);
  end  

  clear Ext_rot n_flow_x n_flow_y n_flow_z  
end

dir_OUT  = fullfile(dir_ARC, 'CrackDensityTable', sprintf('CrackRule_#%02d', RULE));
if ~exist(dir_OUT, 'dir')
  mkdir(dir_OUT)
end

str        = {};

str{end+1} = ' ';
str{end+1} = sprintf('TABLE 1. Crack "density" for model %s', ModelName);
str{end+1} = sprintf('         Run date: %s', datestr(now));
str{end+1} = ' ';
str{end+1} = sprintf('Crack rule : %s', RuleStr);
  
for l=1:numel(ParamStr)
  str{end+1} = ParamStr{l};
end
  
str{end+1} = BlankStr;
str{end+1} = ' ';
str{end+1} = '   Site  Z_site   H_site  Travel time  Path length  Transverse  Longitudinal  Diagonal  Total';
str{end+1} = '         (masl)    (m)       (yr)         (m)         (cnt)       (cnt)        (cnt)   (cnt)';
str{end+1} = ' ';

for s=1:numel(PP_dn)
  if MODEL.GRID_shift==1
    Z_site_masl = P_dn(s).Z_site + MODEL.z_0;
  else
    Z_site_masl = P_dn(s).Z_site;
  end
  str{end+1} = sprintf('   %-4s   %4.0f    %4.0f     %6.2f       %6.1f      %5d       %5d        %5d    %5d ', ...
    char(PP_dn(s).Name_site), Z_site_masl, P_dn(s).H_site, PathTime(s), PathDist(s), N_trans(s), N_long(s), N_diag(s), N_sum(s));
end

file_tbl = fullfile(dir_OUT, 'Table.txt');

f_tbl    = fopen(file_tbl, 'w');

for s=1:numel(str)
  fprintf(1,'%s\n', char(str{s}));
  fprintf(f_tbl, '%s\n', char(str{s}));
end  

fclose(f_tbl);

SITE_list = [P_dn.Name_site];

if isfield(DENSITY,'Special') && strcmp(DENSITY.Special, 'Trapridge')
  FC_list   = {'FC1', 'FC2', 'FC3', 'FC4', 'FC5', 'FC6', 'FC7'};
  M_list    = {'M1', 'M2', 'M3'};
  
  % Test to determine whether P_dn includes all members on the FC and M_lists
  % Some sites could be beyond the glacier margins depending on the ice
  % dynamics model specification (e.g. FC7 for "Trapridge_1.00NoSurge")
  
  L_FC      = ismember(FC_list, SITE_list);
  L_M       = ismember(M_list, SITE_list);
  
  FC_list   = FC_list(L_FC);
  M_list    = M_list(L_M);

  Alt_list  = {};
  TRAPRIDGE = 1;
else
  FC_list   = {};
  M_list    = {};
  Alt_list  = SITE_list;
  TRAPRIDGE = 0;
end  

handle  = 0;

switch TRAPRIDGE
  case 0
    All_trans = N_trans;
    All_long  = N_long;
    All_diag  = N_diag;
    All_names = Alt_list;
    
    if ~isempty(All_names)
      for s=1:numel(All_names)
        CRACK(s).N_trans = All_trans(s);
        CRACK(s).N_long  = All_long(s);
        CRACK(s).N_diag  = All_diag(s);
        CRACK(s).N_sum   = All_trans(s)+All_long(s)+All_diag(s);
        CRACK(s).Site    = All_names{s};
      end
  
      Plot_trans = All_trans;
      Plot_long  = All_long;
      Plot_diag  = All_diag;
      Plot_bars  = [All_trans All_long All_diag];
      Plot_names = All_names;
  
      handle = handle+1;
      figure(handle)
  
      h_bar = bar(Plot_bars, 'stacked'); 
      % axis tight
  
      if isfield(DENSITY, 'TransBarColor')
        h_bar(1).FaceColor = DENSITY.TransBarColor;
      end
      if isfield(DENSITY, 'LongBarColor')
        h_bar(2).FaceColor = DENSITY.LongBarColor;
      end
      if isfield(DENSITY, 'DiagBarColor')
        h_bar(3).FaceColor = DENSITY.DiagBarColor;
      end  
  
      set(gca, 'XTick', 1:numel(Plot_names), 'XTickLabel', Plot_names);

      LegList = {'Transverse', 'Longitudinal', 'Diagonal'};
      legend(h_bar(end:-1:1), LegList{end:-1:1}, 'Location', 'northeast');
      legend boxoff
      title({sprintf('MODEL : %s', ModelName), RuleStr, ParamStr{1:end}, BlankStr}, 'interpreter', 'none')
  
      print(handle, fullfile(dir_PLOT, 'Fig_ALL_bars.pdf'), '-dpdf')

      save(fullfile(dir_OUT, 'CrackCount.mat'), 'MODEL', 'CRACK', 'CRACK_MODEL', 'file_source', 'file_date') 
    end   
  case 1
    if ~any(ismember(FC_list, SITE_list)) && ~any(ismember(M_list, SITE_list))
      fprintf(1,'\nNo plots generated because data are not for observed Trapridge sites\n')
      fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);
      return
    end

    if any(ismember(FC_list, SITE_list))
  
      [LOC_1a, LOC_1b] = ismember(SITE_list, FC_list);

      FC_trans = N_trans(LOC_1a);
      FC_long  = N_long(LOC_1a);
      FC_diag  = N_diag(LOC_1a);
  
      LOC_1b   = LOC_1b(LOC_1b~=0);
  
      FC_site  = FC_list(LOC_1b);

      FC_bars  = [FC_trans FC_long  FC_diag];
  
      handle   = handle+1;

      figure(handle)
  
      h_FC = bar(FC_bars, 'stacked');
  
      if isfield(DENSITY, 'TransBarColor')
        h_FC(1).FaceColor = DENSITY.TransBarColor;
      end
      if isfield(DENSITY, 'LongBarColor')
        h_FC(2).FaceColor = DENSITY.LongBarColor;
      end
      if isfield(DENSITY, 'DiagBarColor')
        h_FC(3).FaceColor = DENSITY.DiagBarColor;
      end  
    
      set(gca, 'XTick', 1:numel(FC_list), 'XTickLabel', FC_list);
      LegList_FC = {'Transverse', 'Longitudinal', 'Diagonal'};
      legend(h_FC(end:-1:1), LegList_FC{end:-1:1}, 'Location', 'northeast')
      legend boxoff
      title({sprintf('MODEL : %s', ModelName), RuleStr, ParamStr{1:end}, BlankStr}, 'interpreter', 'none')
  
      print(handle,fullfile(dir_PLOT, 'Fig_FC_bars.pdf'), '-dpdf')
  
      All_trans = FC_trans';
      All_long  = FC_long';
      All_diag  = FC_diag';
      All_names = {FC_site{1:end}};
    else
      All_trans = [];
      All_long  = [];
      All_diag  = [];
      All_names = {};
    end

    if any(ismember(M_list, SITE_list)) 
      [LOC_2a, LOC_2b] = ismember(SITE_list, M_list);
  
      M_trans = N_trans(LOC_2a);
      M_long  = N_long(LOC_2a);
      M_diag  = N_diag(LOC_2a);
  
      LOC_2b  = LOC_2b(LOC_2b~=0);
      M_site  = M_list(LOC_2b); 

      M_bars  = [M_trans M_long  M_diag];
  
      handle  = handle+1;
  
      figure(handle)

      h_M = bar(M_bars, 'stacked');
  
      if isfield(DENSITY, 'TransBarColor')
        h_M(1).FaceColor = DENSITY.TransBarColor;
      end
      if isfield(DENSITY, 'LongBarColor')
        h_M(2).FaceColor = DENSITY.LongBarColor;
      end
      if isfield(DENSITY, 'DiagBarColor')
        h_M(3).FaceColor = DENSITY.DiagBarColor;
      end  
  
      set(gca, 'XTick', 1:numel(M_list), 'XTickLabel', M_list);
      LegList_M = {'Transverse', 'Longitudinal', 'Diagonal'};
      legend(h_M(end:-1:1), LegList_M{end:-1:1}, 'Location', 'northeast')
      legend boxoff
      title({sprintf('MODEL : %s', ModelName), RuleStr, ParamStr{1:end}, BlankStr}, 'interpreter', 'none')
 
      print(handle,fullfile(dir_PLOT, 'Fig_M_bars.pdf'), '-dpdf')
  
      All_trans = [All_trans M_trans']';
      All_long  = [All_long  M_long']';
      All_diag  = [All_diag  M_diag']';
      All_names = [All_names M_site];
    end
    
    if ~isempty(All_names)
      for s=1:numel(All_names)
        CRACK(s).N_trans = All_trans(s);
        CRACK(s).N_long  = All_long(s);
        CRACK(s).N_diag  = All_diag(s);
        CRACK(s).N_sum   = All_trans(s)+All_long(s)+All_diag(s);
        CRACK(s).Site    = All_names{s};
      end
  
      Plot_trans = [FC_trans' NaN M_trans']';
      Plot_long  = [FC_long' NaN M_long']';
      Plot_diag  = [FC_diag' NaN M_diag']';
      Plot_bars  = [Plot_trans Plot_long Plot_diag];
      Plot_names = {FC_site{1:end} '' M_site{1:end}};
  
      handle = handle+1;
      figure(handle)
  
      h_bar = bar(Plot_bars, 'stacked');
  
      if isfield(DENSITY, 'TransBarColor')
        h_bar(1).FaceColor = DENSITY.TransBarColor;
      end
      if isfield(DENSITY, 'LongBarColor')
        h_bar(2).FaceColor = DENSITY.LongBarColor;
      end
      if isfield(DENSITY, 'DiagBarColor')
        h_bar(3).FaceColor = DENSITY.DiagBarColor;
      end  
  
      set(gca, 'XTick', 1:numel(Plot_names), 'XTickLabel', Plot_names);
      x_line = numel(FC_list)+1;
      hold on
      Y_lim = ylim;
      plot([x_line x_line], Y_lim, 'k')
      LegList = {'Transverse', 'Longitudinal', 'Diagonal'};
      legend(h_bar(end:-1:1), LegList{end:-1:1}, 'Location', 'northeast');
      legend boxoff
      title({sprintf('MODEL : %s', ModelName), RuleStr, ParamStr{1:end}, BlankStr}, 'interpreter', 'none')
  
      print(handle, fullfile(dir_PLOT, 'Fig_ALL_bars.pdf'), '-dpdf')

      save(fullfile(dir_OUT, 'CrackCount.mat'), 'MODEL', 'CRACK', 'CRACK_MODEL', 'file_source', 'file_date') 
    end      
  otherwise
    error('main_calculate_tabulate_and_plot_crack_density(): Unprogrammed case')
end    
  
fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);  
  
