function   [B, S, T, T_S, q_G, xd, yd, td] = LoadBootstrap(MODEL, BootYear)

global nx ny n_XI X2 Y2 X3 Y3 Z3

dir_ARCH = fullfile(pwd, sprintf('Archive-%s', MODEL.name));

load(MODEL.file_dem, 'B')

[Y3, Z3, X3] = meshgrid(MODEL.y_grid, MODEL.XI, MODEL.x_grid);
[X2, Y2]     = meshgrid(MODEL.x_grid, MODEL.y_grid); 

F_name = fullfile(dir_ARCH, 'Snapshots', sprintf('Out(%08.3f).mat', BootYear));

if ~exist(F_name, 'file')
  DIR = dir(fullfile(dir_ARCH, 'Snapshots', 'Out(*).mat'));
  D_name = {DIR.name};
  error('LoadBootstrap(): No output found for run year %.3f: First file is "%s"; last file is "%s"', BootYear, char(D_name(1,1)), char(D_name(1,end)))
end  

load(F_name, 'DD', 'T')

H = zeros(DD.ny, DD.nx);
H(DD.ic_jc_I) = DD.H;
S = B+H;

if MODEL.THERMAL
  T_S = MODEL.T_S;
  q_G = MODEL.q_G;
else  
  T.A     = [];
  T.grid = {};
  T_S    = [];
  q_G    = [];
end

if isfield(DD, 'xd') && ~isempty(DD.xd)
  xd = X3;
  yd = Y3;
  td = ones(size(X3))*BootYear;    
  
  xd(:,DD.ic_jc_I) = DD.xd;
  yd(:,DD.ic_jc_I) = DD.yd;
  td(:,DD.ic_jc_I) = DD.td;
else
  xd = [];
  yd = [];
  td = [];
end
