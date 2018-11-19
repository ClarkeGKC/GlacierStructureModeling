function [Grid, X_NAD27_source, Y_NAD27_source, X_WGS84_source, Y_WGS84_source, Name_source, Start_year, Stop_year, t_last_surge] = GetSourceList(SourceFileName, MODEL)

f_dat = fopen(SourceFileName, 'r');

while 1
  str = fgetl(f_dat);
  if ~ischar(str) || numel(str)<2
    break
  end
  if ~strcmp(str(end), ';')
    str = strcat(str, ';');
  end
  eval(str)
end

if ~exist('X_source', 'var') || ~exist('Y_source', 'var')
  X_source = MODEL.x_grid(i_source);
  Y_source = MODEL.y_grid(j_source);
end

if strcmp(Grid, 'WGS84')
  X_WGS84_source = X_source;
  Y_WGS84_source = Y_source;
           
  X_NAD27_source = X_WGS84_source + 110.557;
  Y_NAD27_source = Y_WGS84_source -165.460;  
else
  X_NAD27_source = X_source;
  Y_NAD27_source = Y_source;

  X_WGS84_source = X_source - 110.557;
  Y_WGS84_source = Y_source + 165.460;
end

if ~exist('Name_source', 'var')
  Name_source = repmat({}, 1, numel(X_source));
end  

if ~exist('Start_year', 'var')
  error('GetSourceList(): Start_year not specified in file "%s"', SourceFileName)
end  

if ~exist('Stop_year', 'var')
  error('GetSourceList(): Stop_year not specified in file "%s"', SourceFileName)
end  

if ~exist('t_last_surge', 'var')
  error('GetSiteList(): No values assigned to structure "t_last_surge"')
elseif isempty(t_last_surge.start) || isempty(t_last_surge.stop)
  error('GetSiteList(): No value(s) assigned to one or all of "t_last_surge.start", "t_last_surge.stop"')
end  

% Verify that the last surge times are consistent with the MODEL-specified
% surge cycle

if t_last_surge.stop-t_last_surge.start ~=MODEL.dt_SURGE
  error('GetSiteList(): "t_last_surge.stop-t_last_surge.start=%.2f " is not consistent with MODEL-specified "dt_SURGE=%.2f"', ...
    t_last_surge.stop-t_last_surge.start, MODEL.dt_SURGE);
end

if ~strcmp(Grid, 'NAD27')
  error('GetSourceList(): Grid for source point(s) is WGS84. Must be converted to NAD27')
end


