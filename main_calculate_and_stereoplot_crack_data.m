function main_calculate_and_stereoplot_crack_data(ModelName)

% Read enhanced track file, perform and even analysis and generate event stereoplots

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
else
  clear global
  clearvars -except ModelName
  ModelDirName = fullfile(pwd, sprintf('Archive-%s', ModelName));
end

global RHO RHO_w g VERBOSITY

POLARITY  = 'BIPOLAR';

VERBOSITY   = 1;
PLOT_TITLE  = 1;

if VERBOSITY>0
  switch POLARITY
    case 'UNIPOLAR'
      fprintf(1,'\nPlotting crack orientation data for Model %s - Rosette plots are UNIPOLAR\n', ModelName);
    case 'BIPOLAR'
      fprintf(1,'\nPlotting crack orientation data for Model %s - Rosette plots are BIPOLAR\n', ModelName);
    otherwise
      error('main_calculate_and_stereoplot_crack_data(): Unprogrammed "POLARITY"=%s', POLARITY)
  end    
end  

CRACK_MODEL = GetCrackModelParameters(fullfile(ModelDirName, 'CrackModel.dat'));

RULE        = CRACK_MODEL.RULE;
RuleStr     = CRACK_MODEL.RuleStr;
ParamStr    = CRACK_MODEL.ParamStr;

BLANKING_RULE  = CRACK_MODEL.BLANKING_RULE;
BLANKING_MODEL = GetBlankingModelParameters(BLANKING_RULE, fullfile(ModelDirName, 'BlankingModel.dat'));
BlankStr   = BLANKING_MODEL.ParamStr;

YRSEC      = 3600*24*365.25;
dir_ARC    = ModelDirName;
file_DAT   = fullfile(dir_ARC, 'DOWN_TRACKS.mat');

DIR        = dir(file_DAT);

fprintf(1,'\n1. Loading down-tracking data (%.2f GB) from "%s"\n', DIR.bytes/10^9, file_DAT);

load(file_DAT, 'MODEL', 'P_dn')

RHO   = MODEL.par.RHO;
RHO_w = MODEL.par.RHO_w;
g     = MODEL.par.g;

% Select surface sites and ignore the remainder

PLOT    = InitializePlotStructure(fullfile(dir_ARC, 'PlotList.dat'), MODEL);
STEREO  = PLOT.Stereo.Crack;

TargetGroup = STEREO.Filter;  
L_match     = true(1, numel(P_dn));

for s=1:numel(P_dn)
  GroupText = P_dn(s).Group.Text;
  L_match(s) = strcmp(TargetGroup, GroupText);
end

P_dn  = P_dn(L_match);

if isfield(STEREO, 'NoPlotList') && ~isempty(STEREO.NoPlotList)
  P_dn = SitePurge(P_dn, STEREO.NoPlotList);
end

ns       = numel(P_dn);

handle   = 0;

dir_PLOT = fullfile(ModelDirName, 'PLOTS');   % If specified this will store hardcopy. If no hardcopy desires then set dir_PLOTS ={};
dir_PLOT = fullfile(dir_PLOT, 'TriptychPlots', sprintf('CrackRule_#%02d', RULE), sprintf('%s', POLARITY)); 


OPTIONS.dir_PLOT        = dir_PLOT;
OPTIONS.Polarity        = POLARITY;
OPTIONS.CrackRule       = RULE;
OPTIONS.ShowNet         = 0;
OPTIONS.PlotTitle       = 0;

OPTIONS.HardCopy        = STEREO.HARDCOPY;
OPTIONS.AzGrid          = 0;
OPTIONS.AzTicks         = 1;
OPTIONS.AzTickInterval  = 90;
OPTIONS.AzTickRfract    = 0.07;
OPTIONS.AzLabelRfract   = 0.25;
OPTIONS.AzLabels        = 'alphabet';
OPTIONS.AzLabelFontScale = 1.0;
OPTIONS.VelArrow        = 1;
OPTIONS.Monochrome      = 1;
OPTIONS.TitleFontScale  = 1.20;
OPTIONS.TextFontScale   = 1.0;
OPTIONS.MarkerScale     = 1.0;
OPTIONS.SurgeMarkers    = 1;
OPTIONS.CullPoints      = 1;
OPTIONS.CullTgtPoints   = 50;


% The files file_tbl and file_vec grow by appending results for each site
% while executing this script. If these files already exist (from previous
% runs) they should be deleted here

file_tbl = fullfile(dir_PLOT, 'EigenvalueTable.dat');
file_vec = fullfile(dir_PLOT, 'EigenpropertiesTable.dat');

if exist(file_tbl, 'file')
  eval(sprintf('delete %s', file_tbl))
end
if exist(file_vec, 'file')
  eval(sprintf('delete %s', file_vec))
end

if ~exist(dir_PLOT, 'dir')
  mkdir(dir_PLOT)
end  

% Eigenvalue table file grows by "appending". It should be destroyed before starting a new run

file_tbl = fullfile(dir_PLOT, sprintf('EigenvalueTable_%s.dat', ModelName));

if exist(file_tbl, 'file')
  delete(file_tbl)
end  
 
for s=1:ns
  
  t_pt     = P_dn(s).t;
  X_pt     = P_dn(s).X;
  Y_pt     = P_dn(s).Y;
  
  dX_pt    = X_pt(2:end)-X_pt(1:end-1);
  dY_pt    = Y_pt(2:end)-Y_pt(1:end-1);
  ds_xy_pt = [0 sqrt(dX_pt.^2 + dY_pt.^2)];
  s_xy_pt  = cumsum(ds_xy_pt);
  
  [L_crack, theta_c, c_type] = ApplyCrackRule(CRACK_MODEL, P_dn(s));
     
  L_use   = L_crack;
  t_use   = t_pt(L_use);
  
  % Eigenvectors D,  s, sigmma and R are identical (but not their eigenvalues)
  
  X_vecs = P_dn(s).vecs_D;
      
  Ext_pts = squeeze(X_vecs(:,1,L_use));   % Select first eigenvector (i.e., the one associated with principal axis of extension)

  % R_post is the rotation that occurs during transport from the site of a
  % crack event to the downstream sampling site. It allows crevasse traces
  % to be rotated from their place (and orientation) of origin to their
  % place of observation
  
  R_Ext   = P_dn(s).R_post(:,:,L_use);
  
  for k=1:sum(L_use)
    Ext_rot(:,k) = squeeze(R_Ext(:,:, k))*squeeze(Ext_pts(:,k));   % Rotate the eigenvector for maximum extension (set at the crack formation point) to that obtaining at the sampling site
    Ext_rot(:,k) = Ext_rot(:,k)/norm(Ext_rot(:,k));                % Ensure that the resulting rotated vector is a still unit vector 
  end
  
  % Assume trajectory downflow end-point is at the sampling site
  % and calculate 2D flow direction vector at that point
    
  n_flow_xy_site(1) = P_dn(s).u(end)/sqrt(P_dn(s).u(end)^2 + P_dn(s).v(end)^2);
  n_flow_xy_site(2) = P_dn(s).v(end)/sqrt(P_dn(s).u(end)^2 + P_dn(s).v(end)^2);
    
  if sum(L_use)>0  
        
    Ext_rot = ProjectToLowerHemisphere(Ext_rot);

    t_p    = t_use;
    x_p    = Ext_rot(1,:);           % Coordinates (x_p, y_p, z_p) are for a unit vector orthogonal to a rotated fault plane
    y_p    = Ext_rot(2,:);
    z_p    = Ext_rot(3,:);
    
    s_xy_p = s_xy_pt(L_use);
       
    P_p = [x_p; y_p; z_p];
    P_p = ProjectToLowerHemisphere(P_p);
    
    x_p = P_p(1,:);
    y_p = P_p(2,:);
    z_p = P_p(3,:);
    
    if BLANKING_RULE>0
      [x_p, y_p, z_p, t_p] = ApplyBlanking(BLANKING_MODEL, x_p, y_p, z_p, t_p, s_xy_p);
    end  
    
    N_pts    = numel(t_p);  % Number of points that survive the blanking rule
    
    % Use random selection of candidate points to obtain target number
    % Also use same random number seed to allow consistent reruns
     
    if OPTIONS.CullPoints    
      N_tgt    = OPTIONS.CullTgtPoints;
       
      if N_pts>N_tgt
        rng('default')
        R      = rand(1, N_pts);
        R_sort = sort(R, 'ascend');
        R_cut  = R_sort(N_tgt);
        L_keep = R<=R_cut;
        
        t_p    = t_p(L_keep);
        x_p    = x_p(L_keep);
        y_p    = y_p(L_keep);
        z_p    = z_p(L_keep);
        s_xy_p = s_xy_p(L_keep);
      end   
    end
       
    if PLOT_TITLE
       PlotTitle = {sprintf('%s', ModelName), RuleStr, ParamStr{1:end}, BlankStr};
    else
       PlotTitle = {};
    end
    handle = Triptych(handle, MODEL, P_dn(s), N_pts, x_p, y_p, z_p, n_flow_xy_site, PlotTitle, OPTIONS, VERBOSITY);
    SaveEigenvalueTable(dir_PLOT, P_dn(s), ModelName, x_p, y_p, z_p, n_flow_xy_site, OPTIONS.Polarity, VERBOSITY)

  else
    fprintf(1,'\nmain_stereoplot_crack_data(): No plot for site %s. No cracks found on this trajectory\n', char(P_dn(s).Name_site));
  end  
  
  clear Ext_rot n_flow_xy_site

end  

fprintf(1,'\nALL DONE\n');  
  
