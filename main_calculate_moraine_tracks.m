function main_calculate_moraine_tracks(ModelName)

% This script uses the moraine debris emergence archive as its input and performs moraine tracking

close all

dir_ROOT = pwd;

if nargin==0
  clear global
  clearvars -except dir_ROOT
  [ModelName, ModelNum, ModelDirName] = PickModel;
  ModelDirName   = fullfile(dir_ROOT, ModelDirName);
elseif nargin==1
  clear global
  clearvars -except dir_ROOT ModelName 
  ModelDirName = fullfile(dir_ROOT, sprintf('Archive-%s', ModelName));
end

fprintf(1,'\nLAUNCHING :: main_calculate_moraine_tracks(''%s'') - Supraglacial particle tracking of moraine debris\n', ModelName);

setpref('Internet', 'E_mail', 'clarke@eos.ubc.ca');
setpref('Internet', 'SMTP_Server', 'smtp.eos.ubc.ca');

MoraineMapList = fullfile(ModelDirName, 'MoraineMapList.dat');

t_start     = tic;

file_date   = datestr(now);
file_source = fullfile(pwd, mfilename);

file_inp     = fullfile(ModelDirName, 'MoraineTracks', 'DEBRIS_EXIT.mat');
file_dat     = MoraineMapList;
dir_arch_out = fullfile(ModelDirName, 'Snapshots');

OK = CheckExistence({file_inp, dir_arch_out});
if ~OK
  return
end  

[yr_start_list, yr_display_list] = GetMoraineDataList(file_dat);

yr_start = min(yr_start_list);
yr_stop  = max(yr_display_list);

dir_moraine_out = fullfile(ModelDirName, 'MoraineTracks');

if ~exist(dir_moraine_out, 'dir')
  mkdir(dir_moraine_out)
end  

load(file_inp, 'MODEL', 'SOURCE')

ns      = numel(SOURCE);    % number of moraine debris sources
dt      = MODEL.dt;
t_CYCLE = MODEL.t_CYCLE;
[X,Y]   = meshgrid(MODEL.x_grid, MODEL.y_grid);

n_XI    = MODEL.n_XI;

DIR = dir(fullfile(dir_arch_out, 'Out(*).mat'));

FileList     = {DIR.name};
FileYearList = ExtractYearListFromFileList(FileList);

N_maps       = numel(yr_start_list);

YearList     = yr_start:dt:yr_stop;
nt           = numel(YearList);

t_last_surge = SOURCE(1).t_last_surge;

X_e   = [SOURCE.X_e];
Y_e   = [SOURCE.Y_e];
t_e   = [SOURCE.t_e];

N_max   = 10*numel(X_e);   % Allows 10 surge cycles of run-time

X_trk   = NaN(N_max, 1);
Y_trk   = NaN(N_max, 1);
H_trk   = NaN(N_max, 1);
s_trk   = NaN(N_max, 1);
t_trk   = NaN(N_max, 1);

ne    = zeros(ns, 1);

for s=1:ns
  ne(s) = numel(SOURCE(s).X_e);
end

cumsum_ne = cumsum(ne);
ns_lo     = [1 cumsum_ne(1:end-1)'+1]';
ns_hi     = [cumsum_ne];

s_e   = zeros(1, numel(X_e));
for s=1:ns
  s_e(1,ns_lo(s):ns_hi(s)) = s;       
end

n_map      = 1;
yr_display = yr_display_list(n_map);
PLOT       = InitializeMorainePlotStructure(yr_start, SOURCE);
handle     = 0;

for s=1:ns
  PLOT(s).yr_display = yr_display;
end

len_old = 0;

for t=1:nt
  Year    = YearList(t);
  if mod(Year, 100)==0
    fprintf(1,'%d', Year);
  elseif mod(Year,50)==0
    fprintf(1,'50');
  elseif mod(Year, 10)==0
    fprintf(1,'.');
  end
  if t==numel(YearList)
    fprintf(1,'\n');
  end   
  
  [DD, V, ~, ~] = GetAliasYearData(YearList(t), t_last_surge, dir_arch_out, FileList, FileYearList, MODEL, 0); 
  V_S_x_I = V.x(n_XI,:);
  V_S_y_I = V.y(n_XI,:);
  ic_jc_I = DD.ic_jc_I;
  
  H          = zeros(DD.ny, DD.nx);
  H(ic_jc_I) = DD.H;
  
  U_s     = zeros(DD.ny, DD.nx);
  V_s     = zeros(DD.ny, DD.nx);
  
  U_s(ic_jc_I) = V_S_x_I;
  V_s(ic_jc_I) = V_S_y_I;
    
  if ~exist('H_now', 'var')
    U_now = U_s;
    V_now = V_s;
    H_nxt = H;
    U_nxt = U_s;
    V_nxt = V_s;
  else  
    H_now = H_nxt;
    U_now = U_nxt;
    V_now = V_nxt;
    H_nxt = H;
    U_nxt = U_s;
    V_nxt = V_s;
  end
    
  H_mid   = H_nxt;
  U_mid   = 0.5*(U_now+U_nxt);
  V_mid   = 0.5*(V_now+V_nxt);
  
  % I am not fully confident in the logic of the following line
    
  L_yr    = mod(t_e, t_CYCLE)==mod(Year, t_CYCLE);
  
  % Introduce newly emerged makers (if any)
    
  if any(L_yr)
    len_new = len_old+sum(L_yr);
    L_old   = true(len_old,1);
    L_new   = true(len_new,1);
 
    X_e_yr  = X_e(L_yr);
    Y_e_yr  = Y_e(L_yr);
    s_e_yr  = s_e(L_yr);
    
    t_e_yr  = Year*ones(1, sum(L_yr));
  
    X_trk(L_new)   = [X_trk(L_old)' X_e_yr]';
    Y_trk(L_new)   = [Y_trk(L_old)' Y_e_yr]';
    s_trk(L_new)   = [s_trk(L_old)' s_e_yr]';
    t_trk(L_new)   = [t_trk(L_old)' t_e_yr]';
    H_trk(L_new)   = interp2(X, Y, H_mid, X_trk(L_new), Y_trk(L_new));
  else  
    len_new        = len_old;
    L_new          = true(len_old, 1);  
  end
  
  % Note that L_new is necessary because X_trk etc have been initialized to large arrays
  % of NaN values and we only wish to calculate the "live" markers
  
  % The following is the crudest approach to calculating particle
  % trajectories but probably adequate for the purpose

  dX_trk         = dt*interp2(X,Y,U_mid, X_trk(L_new), Y_trk(L_new));
  dY_trk         = dt*interp2(X,Y,V_mid, X_trk(L_new), Y_trk(L_new));
  
  X_trk(L_new)   = X_trk(L_new)+dX_trk;
  Y_trk(L_new)   = Y_trk(L_new)+dY_trk;
   
  L_del  = H_trk(L_new)<=0;
  L_tail = isnan(X_trk);
   
  if any(L_del)
    X_trk   = [X_trk(~L_del)' X_trk(L_tail)']';
    Y_trk   = [Y_trk(~L_del)' Y_trk(L_tail)']';
    s_trk   = [s_trk(~L_del)' s_trk(L_tail)']';
    t_trk   = [t_trk(~L_del)' t_trk(L_tail)']';
    H_trk   = [H_trk(~L_del)' H_trk(L_tail)']';
    len_new = len_new-sum(L_del);
  end
  
  if abs(YearList(t)-yr_display)<0.5*dt
    for s=1:ns
      PLOT(s).X_trk   = X_trk(s_trk==s);
      PLOT(s).Y_trk   = Y_trk(s_trk==s);
      PLOT(s).H_trk   = H_trk(s_trk==s);
      PLOT(s).t_trk   = t_trk(s_trk==s);
    end
    if ~exist(dir_moraine_out, 'dir')
      mkdir(dir_moraine_out);
    end
 
    handle = handle+1;

    figure(handle)
    imagesc(MODEL.x_grid, MODEL.y_grid, double(H>0)), colorbar, axis equal, axis tight   % Ice mask plot. Hence no NaN values
    set(gca, 'YDir', 'normal')

    colormap(gray)
    hold on

    PenList   = {'k.', 'b.', 'r.', 'g.', 'c.', 'm.', 'y.'};
    ColorList = {'k',  'b',  'r',  'g',  'c',  'm',  'y'};

    np      = numel(PenList);

    for s=1:ns
      plot([PLOT(s).X_trk], [PLOT(s).Y_trk], PenList{1+mod(s-1, np)})
      hold on
      plot(SOURCE(s).X_source, SOURCE(s).Y_source, 'k+')
      plot(SOURCE(s).X_source, SOURCE(s).Y_source, 'Marker', 'o', 'MarkerEdgeColor', ColorList{1+mod(s-1,np)})  
      
      if floor(yr_display)==yr_display
        title({sprintf('MODEL %s at %04d', ModelName, yr_display), 'All sources'}, 'Interpreter', 'None')
      else
        title({sprintf('MODEL %s at %.2f', ModelName, yr_display), 'All sources'}, 'Interpreter', 'None')
      end
    end

    hold off
    
    if floor(yr_display)==yr_display
      print(handle, fullfile(dir_moraine_out, sprintf('Fig_Moraine_tracks(%04d).pdf',yr_display)), '-dpdf')
      file_out = fullfile(dir_moraine_out, sprintf('MoraineTracks(%d).mat', yr_display));
    else
      print(handle, fullfile(dir_moraine_out, sprintf('Fig_Moraine_tracks(%.2f).pdf',yr_display)), '-dpdf')
      file_out = fullfile(dir_moraine_out, sprintf('MoraineTracks(%.2f).mat', yr_display));
    end
    
    save(file_out, 'MODEL', 'SOURCE', 'DD', 'PLOT', 'file_date', 'file_source');
    
    fprintf(1,'01. Saving moraine tracks for yr_display=%.2f\n', yr_display);
    fprintf(1,'    Used space in track arrays is %d (N_max=%d)\n', len_new, N_max); 
      
    n_map = n_map+1;
    if n_map>numel(yr_display_list)
      break
    else
     yr_display = yr_display_list(n_map);
     for s=1:ns
       PLOT(s).yr_display = yr_display;
     end  
    end  
  end
  
  len_old        = len_new; 
end

fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);

sendmail('clarke@eos.ubc.ca', sprintf('main_calculate_moraine_tracks(): Moraine tracking run "%s" completed', ModelName), ...
  sprintf('Model %s completed at %s. Elapsed time %.2f hr.', ModelName, datestr(now), toc(t_start)/3600));
