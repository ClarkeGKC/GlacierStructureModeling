function main_track_closure_check(ModelName)

close all

if nargin==0
  clear all
  [ModelName, ModelNum, ModelDirName] = PickModel;
elseif nargin==1
  clear global
  clearvars -except ModelName
end  

fprintf(1,'\nLAUNCHING :: main_track_closure_check(''%s'') - Compare results of bottom-up with top-down particle tracking (as an accuracy check)\n', ModelName);

dir_ROOT = fullfile(pwd, sprintf('Archive-%s', ModelName));

file_UP  = fullfile(dir_ROOT, 'UP_TRACKS.mat');
file_DN  = fullfile(dir_ROOT, 'DOWN_TRACKS.mat');
 
load(file_UP, 'MODEL', 'SITE', 'P_up')

X_site_up    = [SITE.X_site];
Y_site_up    = [SITE.Y_site];
Z_site_up    = [SITE.Z_site];
t_site_up    = [SITE.t_site];

Name_site_up = [SITE.Name_site];

clear MODEL SITE

load(file_DN, 'MODEL', 'P_dn');

X_site    = [P_dn.X_site];
Y_site    = [P_dn.Y_site];
Z_site    = [P_dn.Z_site];
t_site    = [P_dn.t_site];

Name_site = [P_dn.Name_site];

% Compare the sites to ensure that up-tracking and down-tracking sites are identical

A   = setdiff(X_site, X_site_up);
B   = setdiff(Y_site, Y_site_up);
C   = setdiff(Z_site, Z_site_up);
D   = setdiff(t_site, t_site_up);
E   = setdiff(Name_site, Name_site_up);

L_test = isempty(A) && isempty(B) && isempty(C) && isempty(D) && isempty(D);

if ~L_test
  error('main_track_closure_check(): The up-tracking and down-tracking site lists differ')
end

for s=1:numel(X_site)
  X_start(s) = P_up(s).X(1);
  Y_start(s) = P_up(s).Y(1);
  Z_start(s) = P_up(s).Z(1);
  t_start(s) = P_up(s).t(1);
  
  X_end(s)   = P_dn(s).X(end);
  Y_end(s)   = P_dn(s).Y(end);
  Z_end(s)   = P_dn(s).Z(end);
  t_end(s)   = P_dn(s).t(end);
end

X_err = X_end-X_start;
Y_err = Y_end-Y_start;
Z_err = Z_end-Z_start;
t_err = t_end-t_start;
  

str = {};

str{end+1} = sprintf('Table 1. Closure checks for bottom-up vs. top-down particle tracking for "%s"', ModelName);
str{end+1} = ' ';
str{end+1} = '       ------- Start of bottom-up track -----     ------ End of top-down track ----------      ---------- Closure errors ---------';
str{end+1} = 'Site   X_start    Y_start   Z_start   t_site           X_end     Y_end    Z_end    t_chk        err(X)    err(Y)   err(Z)   err(t)';
str{end+1} = '         (m)       (m)        (m)      (yr)             (m)       (m)      (m)     (yr)          (m)       (m)      (m)      (yr)';

for s=1:numel(X_site)
  str{end+1} = sprintf('%-4s  %8.3f  %8.3f  %8.3f  %8.3f      %8.3f  %8.3f  %8.3f  %8.3f     %8.3f  %8.3f  %8.3f  %8.3f', ...
    char(Name_site{s}), X_start(s), Y_start(s), Z_start(s), t_start(s), ...
      X_end(s), Y_end(s), Z_end(s), t_end(s), X_err(s), Y_err(s), Z_err(s), t_err(s));
end

file_out = fullfile(dir_ROOT, 'TrackClosureCheck.txt');

f_out = fopen(file_out, 'w');

for is=1:numel(str)
  fprintf(1,'%s\n', str{is});
  fprintf(f_out, '%s\n', str{is});
end

fclose(f_out);




 

