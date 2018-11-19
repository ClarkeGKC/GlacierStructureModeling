function [ModelName, ModelNum, ModelDirName] = PickModel(ModelFilterString)

D = dir(fullfile(pwd,'Archive-*'));
n_len  = length('Archive-');

D_list    = {D.name};

if nargin==1
  D_temp = D_list;
  D_list = {};
  for l=1:numel(D_temp)
    if ~isempty(strfind(D_temp{l},ModelFilterString))
      D_list{end+1} = D_temp{l};
    end
  end
end

for d=1:numel(D_list)
  ModelList{d} = D_list{d}(n_len+1:end);
end

if numel(D_list)==0
  error('PickModel(): No matching models found')
elseif numel(ModelList)==1
  ModelName = ModelList{1};
  ModelNum  = 1;
  ModelDirName = D_list{1};
  return
end  

while 1
  fprintf(1,'\nSELECT MODEL FROM LIST\n\n');
  for d=1:numel(ModelList)
    fprintf(1,'%2d - %s\n', d, ModelList{d});
  end
  key = input(sprintf('\nEnter selection (1,2,3 ... %d) ... ', numel(D_list)), 's');
  if str2num(key)>0 && str2num(key)<=numel(ModelList)
    NumKey = str2num(key);
    break
  end
end

ModelName    = ModelList{NumKey};
ModelNum     = NumKey;
ModelDirName = D_list{NumKey};



  