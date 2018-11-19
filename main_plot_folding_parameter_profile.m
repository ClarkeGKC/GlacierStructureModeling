function main_plot_folding_parameter_profile(ModelName)

% This script plots the folding parameter

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
  error('main_calculate_and_plot_folding_parameter_profile(): Unprogrammed number of input arguments')
end  

t_start = tic;

dir_ARC     = ModelDirName;
dir_SNAP    = fullfile(dir_ARC, 'Snapshots');
dir_PLOT    = fullfile(dir_ARC, 'PLOTS', 'FoldingParameter');

file_DAT    = fullfile(dir_ARC, 'FoldingParameter', 'FOLD_PAR.mat');

DIR         = dir(file_DAT);

fprintf(1,'\n1. Loading down-tracking data (%.2f GB) from "%s"\n', DIR.bytes/10^9, file_DAT);

load(file_DAT, 'MODEL', 'F_par')

ns      = numel(F_par);
nx      = MODEL.nx;
ny      = MODEL.ny;
n_JJ    = 11;        % Number of quadratic invariants

PLOT    = InitializePlotStructure(fullfile(dir_ARC, 'PlotList.dat'), MODEL);
PROFILE = PLOT.Profile.FoldingParameter;

n_PROFILES = numel(PROFILE);
yr_site = [F_par.t_site];

if any(yr_site~=yr_site(1))
  error('main_calculate_and_plot_folding_parameter_profile(): All strain sites must have same observation year')
else
  yr_site       = yr_site(1);
  yr_last_surge = F_par(1).yr_last_surge;
end

DIR = dir(fullfile(dir_SNAP, '*.mat'));

ModelFileList = {DIR.name};
ModelYearList = ExtractYearListFromFileList(ModelFileList);

[D, ~, ~, ~, ~, YearSwap] = GetAliasYearData(yr_site, yr_last_surge, dir_SNAP, ModelFileList, ModelYearList, MODEL, 0);

H            = zeros(ny,nx);
H(D.ic_jc_I) = D.H;
handle       = 0;

for p=1:n_PROFILES 
  L_profile  = false(1, ns);
  
  NoPlotList = PROFILE(p).NoPlotList;
  TextSearch = PROFILE(p).Filter;

  for s=1:ns
    t_min(s) = F_par(s).t(1);
    OK = 1;
    
    if numel(NoPlotList)>0
      for q=1:numel(NoPlotList)
        if ~isempty(strfind(char(F_par(s).Name_site), NoPlotList{q}))
          OK = 0;
        end
      end
    end
    
    if OK
      ng       = numel(F_par(s).Group);
      for gg=1:ng
        if strcmp(F_par(s).Group(gg).Text, TextSearch)
          L_profile(s) = 1;
        end
      end
    end
  end
  
  if sum(L_profile)==0
    fprintf(1,'main_calculate_and_plot_folding_parameter_profile(): No group points found for Profile p=%d with search text "%s"', p, TextSearch)
    continue
  end
  
  FF_par = F_par(L_profile);
  t_min  = t_min(L_profile);

  X_all  = [FF_par.X_site];
  Y_all  = [FF_par.Y_site];
  XI_all = [FF_par.XI_site];
  Name_all = [FF_par.Name_site];
  
  J_use  = false;
  JJ     = false(1,11);
  
  for par=1:numel(PROFILE(p).Parameter)
    if strcmp(PROFILE(p).Parameter{par}, 'J')
      J_use = 1;
    else  
      eval(sprintf('%s = 1;', PROFILE(p).Parameter{par}) )
    end
  end
  
  JJ_cnt  = 1:n_JJ;
  JJ_list = JJ_cnt(JJ); 

  n_PT   = numel(X_all);

  PT     = [X_all; Y_all; XI_all];

  x_plot_axis = PROFILE(p).x_plot_axis';
  y_plot_axis = PROFILE(p).y_plot_axis';
  n           = cross(x_plot_axis, y_plot_axis);  % Outward normal to plane of x-y plot
  
  PT_0    = PROFILE(p).Origin';
  PT_proj = PT - dot(PT - PT_0, repmat(n, 1, n_PT), 1).*repmat(n, 1, n_PT);
  
  X_c_plot  = dot(PT_proj, repmat(x_plot_axis, 1, n_PT), 1);
  Y_c_plot  = dot(PT_proj, repmat(y_plot_axis, 1, n_PT), 1);
  
  B_fold = NaN(n_PT, 1);
  H_fold = NaN(n_PT, 1);
       
  for ii=1:n_PT
     B_fold(ii)  = FF_par(ii).B(end);
     H_fold(ii)  = FF_par(ii).H(end);
  end
  
  [x_grid, y_grid, index] = GridPointRowColumnIndices(X_c_plot, Y_c_plot);
  
  n_y = numel(y_grid);
  n_x = numel(x_grid);
  
  B_fold_grid = reshape(B_fold(index), n_y, n_x);
  H_fold_grid = reshape(H_fold(index), n_y, n_x);
    
  B_profile = B_fold_grid(1,:);
  H_profile = H_fold_grid(1,:);
  S_profile = B_profile+H_profile;
  
  XX = repmat(x_grid, n_y, 1);
  YY = B_fold_grid + repmat(y_grid', 1, n_x).*H_fold_grid;
  
  if J_use
    J_fold  = NaN(n_PT, 1);     % linear invariant

    handle = handle+1;
    figure(handle)
    
    for ii=1:n_PT
      J_fold(ii) = FF_par(ii).J(end);
    end
    J_fold_grid = reshape(J_fold(index), n_y, n_x);
    
    if PROFILE(p).FillContours
      contourf(XX, YY, J_fold_grid, PROFILE(p).ContourLevels, 'LineStyle', PROFILE(p).ContourLineStyle), axis tight, colorbar
    else
      contour(XX, YY, J_fold_grid, PROFILE(p).ContourLevels, 'LineStyle', PROFILE(p).ContourLineStyle), axis tight, colorbar
    end  
    hold on  
       
    yL = ylim;
    
    if isfield(PROFILE, 'Gray') && ~isempty(PROFILE(p).Gray)
      GrayPatch = [x_grid     x_grid(end) x_grid(1) x_grid(1); ...
                   B_profile    yL(1)      yL(1)     B_profile(1)];
      patch(GrayPatch(1,:), GrayPatch(2,:), PROFILE(p).Gray)           
    end
    
    plot(x_grid, B_profile, 'k')
    plot(x_grid, S_profile, 'k')
      
    xlabel(PROFILE(p).xlabel)
    ylabel(PROFILE(p).ylabel)
    title_str = {sprintf('MODEL %s', MODEL.name), PROFILE(p).title, 'Folding parameter J'};
    title(title_str, 'interpreter', 'none')
    hold off
  
    if PROFILE(p).HARDCOPY
      if ~exist(dir_PLOT, 'dir')
        mkdir(dir_PLOT)
      end  
      print(handle, fullfile(dir_PLOT, sprintf('Profile(%d)_Contour_J.pdf', p)), '-dpdf')
    end  
  end
  
  if numel(JJ_list)>0
    for j=1:numel(JJ_list)
      J_fold = NaN(n_PT, 1); 

      handle = handle+1;
      figure(handle)
      
      for ii=1:n_PT
        J_fold(ii) = FF_par(ii).JJ(JJ_list(j), end);
      end
      
      if any(J_fold<0)
        fprintf(1,'WARNING one or more JJ(%d) values is negative\n', JJ_list(j))
      end
      
      J_fold      = sqrt(abs(J_fold));
      J_fold_grid = reshape(J_fold(index), n_y, n_x);
    
      if PROFILE(p).FillContours
        contourf(XX, YY, log10(J_fold_grid), PROFILE(p).ContourLevels, 'LineStyle', PROFILE(p).ContourLineStyle), axis tight, colorbar
      else
        contour(XX, YY, log10(J_fold_grid), PROFILE(p).ContourLevels, 'LineStyle', PROFILE(p).ContourLineStyle), axis tight, colorbar
      end  
      hold on  
      
      if isfield(PROFILE, 'cmin') && isfield(PROFILE, 'cmax') && ~isempty(PROFILE(p).cmin) && ~isempty(PROFILE(p).cmax)
        [cmin, cmax] = caxis;
        if cmin<PROFILE(p).cmin || cmax>PROFILE(p).cmax
          fprintf(1,'main_plot_folding_parameter_profile(): Reqested cmin cmax range does not bracket the range of the data ****** APPLY ANYWAY\n')'
          fprintf(1,'  Data cmin = %.3f\n', cmin);
          fprintf(1,'  Data cmax = %.3f\n\n', cmax);
        end
        caxis([PROFILE(p).cmin PROFILE(p).cmax])
      end
      
      yL = ylim;
    
      if isfield(PROFILE, 'Gray') && ~isempty(PROFILE(p).Gray)
        GrayPatch = [x_grid     x_grid(end) x_grid(1) x_grid(1); ...
                     B_profile    yL(1)      yL(1)     B_profile(1)];
        patch(GrayPatch(1,:), GrayPatch(2,:), PROFILE(p).Gray)           
      end
          
      plot(x_grid, B_profile, 'k')
      plot(x_grid, S_profile, 'k')

      xlabel(PROFILE(p).xlabel)
      ylabel(PROFILE(p).ylabel)
      title_str = {sprintf('MODEL %s', MODEL.name), PROFILE(p).title, sprintf('Folding parameter JJ(%d)', JJ_list(j))};
      title(title_str, 'interpreter', 'none')
      hold off
  
      if PROFILE(p).HARDCOPY
        if ~exist(dir_PLOT, 'dir')
          mkdir(dir_PLOT)
        end  
        print(handle, fullfile(dir_PLOT, sprintf('Profile(%d)_Contour_JJ(%d).pdf', p, JJ_list(j))), '-dpdf')
      end  
    end
  end
end

fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);







