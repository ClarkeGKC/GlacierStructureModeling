function [FileSwap, YearSwap] = GetYearSwapFileName(MODEL, t_last_surge, ModelFileList, ModelYearList, YR_tgt)

% Calculate alias year that matches target year and identify that file

t_CYCLE = MODEL.t_CYCLE;
IPREC   = GetTimeStepPrecision(MODEL.dt);

if numel(t_CYCLE)>1
  t_CYCLE  = t_CYCLE(1);
end  

L_start       = mod(ModelYearList, t_CYCLE)==0;
ModelYearStartList = ModelYearList(L_start);

% The following is the stored Model year that matches to the start of the last modelled surge

YearsSinceLastStart = YR_tgt-t_last_surge.start;  % Time since the start of last surge (can be negative)
YearsSinceLastStart = mod(YearsSinceLastStart, t_CYCLE);

DeltaYears = mod(round(YearsSinceLastStart-ModelYearList, IPREC), t_CYCLE);
cnt        = 1:numel(ModelYearList);
i_match    = cnt(DeltaYears==0);

if numel(i_match)>0
  YearSwap = ModelYearList(max(i_match));
  FileSwap = ModelFileList{max(i_match)};
else
  [~, i_min] = min(DeltaYears);
  
  YearAlt    = ModelYearList(i_min);
  FileAlt    = ModelFileList{i_min};
  fprintf(1,'*** WARNING ***\n');
  fprintf(1,'Desired Swap file %s (for YearGet=%.3f and YearSwap=%.3f) does not exist. Best substitute file (for %.3f) is applied\n', ...
    FileAlt, YR_tgt, YearAlt);
  YearSwap  = YearAlt;
  FileSwap  = FileAlt;
end
