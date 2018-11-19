function main_calculate_and_save_grad_L(ModelName)

% This calculates grad_L which is required to calculate grad_F and other quantities
% relevant to deformation structures and folding

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
  error('main_calculate_and_save_grad_L(): Unprogrammed number of input arguments')
end  

fprintf(1,'\nLAUNCHING :: main_calculate_and_save_grad_L(''%s'') - Calculate grad_L (rank-3 gradient velocity gradient tensor) from archived L tensor field\n', ModelName);

global X2 Y2 X3 Y3 Z3 XI IC_JC

t_start = tic;

file_date   = datestr(now);
file_source = fullfile(pwd, mfilename);

setpref('Internet', 'E_mail', 'clarke@eos.ubc.ca');
setpref('Internet', 'SMTP_Server', 'smtp.eos.ubc.ca');

dir_ARC     = ModelDirName;
dir_SNAP    = fullfile(dir_ARC, 'Snapshots');
dir_GRAD    = fullfile(dir_ARC, 'grad_L');

load(fullfile(dir_ARC, 'Model.mat'), 'MODEL')

DIR = dir(fullfile(dir_SNAP, '*.mat'));

ModelFileList = {DIR.name};
ModelYearList = ExtractYearListFromFileList(ModelFileList);

if ~exist(dir_GRAD, 'dir')
  mkdir(dir_GRAD)
end

for f=1:numel(ModelFileList)
  if mod(f, 1000)==0
    fprintf(1,'%d', f);
  elseif mod(f,500)==0
    fprintf(1,'500');
  elseif mod(f, 100)==0
    fprintf(1,'.');
  end
  if f==numel(ModelYearList)
    fprintf(1,'\n');
  end 
  load(fullfile(dir_SNAP, ModelFileList{f}), 'DD', 'V')
  grad_L = Calculate_grad_L(MODEL, DD, V);
  file_L = fullfile(dir_GRAD, sprintf('grad_L(%08.3f).mat', ModelYearList(f)));
  save(file_L, 'grad_L', 'file_source', 'file_date')
end     

fprintf(1,'\nALL DONE. Elapsed time %.2f min\n', toc(t_start)/60);

sendmail('clarke@eos.ubc.ca', sprintf('main_calculate_and_save_grad_L(): Gradient of L tensor evolution run "%s" completed', ModelName), ...
  sprintf('Model %s completed at %s. Elapsed time %.2f hr.', ModelName, datestr(now), toc(t_start)/3600));





