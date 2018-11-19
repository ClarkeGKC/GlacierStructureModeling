function MODEL = GetModelData(file_dat)

% Set default values of MODEL structure

f_dat = fopen(file_dat, 'r');

if f_dat==-1
  error('GetModelData(): Unable to open run data file "%s"', file_dat)
end 

while 1
  str = fgetl(f_dat);
  if ~ischar(str) || length(str)<2
    break
  end

  if ~strcmp(str(end), ';')
    str = strcat(str, ';');
  end
  
  if contains(str, '../')
    i_match = findstr(str,'../');
    str = fullfile(str(1:i_match-1), pwd, str(i_match+3:end));
  end
  
  eval(str)
end  
 
if exist('name', 'var')
  MODEL.name = name;
end

if exist('descriptor', 'var')
  MODEL.descriptor = descriptor;
else
  MODEL.descriptor = {};
end

if exist(file_dem, 'file')
  MODEL.file_dem = file_dem;
end

% Note that multiple sliding masks can be used to generate complex surge patterns

for f=1:numel(file_slide)
  if exist(file_slide{f}, 'file')
    MODEL.file_slide{f} = file_slide{f};
  end
end

if exist(file_balance, 'file')
  MODEL.file_balance = file_balance;
end

if exist('file_exclude', 'var')
  MODEL.file_exclude = file_exclude;
else
  MODEL.file_exclude = {};
end  

if exist('dir_bootstrap', 'var')
  MODEL.dir_bootstrap = dir_bootstrap;
  MODEL.BOOTSTRAP     = true;
else
  MODEL.dir_bootstrap = {};
  MODEL.BOOTSTRAP     = false;
end

if exist('dx', 'var')
  MODEL.dx = dx;
else
  error('GetModelData(): dx not defined in data file')
end
  
if exist('dy', 'var')
  MODEL.dy = dy;
else
  error('GetModelData(): dy not defined in data file')
end

if exist('nx', 'var')
  MODEL.nx = nx;
else
  error('GetModelData(): nx not defined in data file')
end
if exist('ny', 'var')
  MODEL.ny = ny;
else
  error('GetModelData(): ny not defined in data file')
end

MODEL.N_xy = nx*ny;

if exist('dt', 'var')
  MODEL.dt = dt;
end
if exist('n_sub_steps', 'var')
  MODEL.n_sub_steps = n_sub_steps;
end
if exist('t_START', 'var')
  MODEL.t_START = t_START;
end
if exist('t_STOP', 'var')
  MODEL.t_STOP  = t_STOP;
end
if exist('t_archive_START', 'var')
  MODEL.t_archive_START = t_archive_START;
end 
if exist('t_tracer_START', 'var')
  MODEL.t_tracer_START = t_tracer_START;
else
  MODEL.t_tracer_START = tg_START;
end 
if exist('m_SLIDE', 'var')
  MODEL.m_SLIDE = m_SLIDE;
end 
if exist('C_SLIDE', 'var')
  MODEL.C_SLIDE = C_SLIDE;
end
if exist('C_SURGE', 'var')
  MODEL.C_SURGE = C_SURGE;
end           
if exist('t_CYCLE', 'var')
  MODEL.t_CYCLE = t_CYCLE;
end
if exist('dt_SURGE', 'var')
  MODEL.dt_SURGE = dt_SURGE;
end
if exist('dt_PEAK', 'var')
  MODEL.dt_PEAK = dt_PEAK;
end
if exist('dL_front', 'var')
  MODEL.dL_front = dL_front;
end
if exist('X_mask_min', 'var')
  MODEL.X_mask_min = X_mask_min;
end
if exist('X_mask_max', 'var')
  MODEL.X_mask_max = X_mask_max;
end

if exist('Y_mask_min', 'var')
  MODEL.Y_mask_min = Y_mask_min;
end
if exist('X_mask_max', 'var')
  MODEL.Y_mask_max = Y_mask_max;
end

if exist('THERMAL', 'var')
  MODEL.THERMAL = THERMAL;
else
  MODEL.THERMAL = 0;
end

if exist('TRACER', 'var')
  MODEL.TRACER = TRACER;
else
  MODEL.TRACER = 0;
end

if exist('TRACER_METHOD', 'var')
  MODEL.TRACER_METHOD = TRACER_METHOD;
else
  MODEL.TRACER_METHOD = '1S2T';
end  
  
if exist('MB_CHECK', 'var')
  MODEL.MB_CHECK = MB_CHECK;
else
  MODEL.MB_CHECK = 0;
end

if exist('VERBOSITY', 'var')
  MODEL.VERBOSITY = VERBOSITY;
else
  MODEL.VERBOSITY = 1;
end  

if exist('GRID_shift', 'var')
  MODEL.GRID_shift = GRID_shift;
else
  MODEL.GRID_shift = 0;
end  

if exist('PLOT', 'var')
  MODEL.PLOT = PLOT;
else
  MODEL.PLOT = 1;
end  

if exist('n_XI', 'var')
  MODEL.n_XI = n_XI;
else
  error('GetModelData(): n_XI not specified in data file');
end

MODEL.XI   = linspace(0, 1, n_XI);
MODEL.d_XI = 1/(n_XI-1);

if exist('n_ZETA', 'var')
  MODEL.n_ZETA = n_ZETA;
else
  error('GetModelData(): n_ZETA not specified in data file')
end

MODEL.ZETA   = linspace(0, 1, n_ZETA);
MODEL.d_ZETA = 1/(n_ZETA-1);


