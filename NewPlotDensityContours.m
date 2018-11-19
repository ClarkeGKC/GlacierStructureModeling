function handle_out = NewPlotDensityContours(handle, MODEL, SITE, x_p, y_p, z_p, n_vel_site, OPTIONS)

% (x_p, y_p, z_p) are coordinates of the unit vector n_p that is
% associated with the stereographic projection of crack "p"

% This version does not employ histogram binning but assigns to each point
% a distribution function (e.g., 2D Gaussian) and the summed result is
% contoured

if isfield(OPTIONS, {'dir_PLOT'})
  dir_PLOT = OPTIONS.dir_PLOT;
else
  dir_PLOT = [];
end  

if isfield(OPTIONS,{'ShowNet'})
  ShowNet = OPTIONS.ShowNet;
else
  ShowNet = true;
end  

if isfield(OPTIONS, {'PlotTitle'})
  PlotTitle = OPTIONS.PlotTitle;
else
  PlotTitle = true;
end  

if isfield(OPTIONS, {'HardCopy'})
  HardCopy = OPTIONS.HardCopy;
else
  HardCopy = false;
end

if isfield(OPTIONS, {'VelArrow'})
  VelArrow = OPTIONS.VelArrow;
else
  VelArrow = false;
end  

if isfield(OPTIONS, {'Monochrome'})
  Monochrome = OPTIONS.Monochrome;
else
  Monochrome = false;
end

phi_grid      = linspace(0, 90, 46);
lambda_grid   = linspace(0, 360, 73);

phi_grid_m    = 0.5*(phi_grid(1:end-1) + phi_grid(2:end));
lambda_grid_m = 0.5*(lambda_grid(1:end-1) + lambda_grid(2:end));

phi_grid_m    = pi*phi_grid_m/180;
lambda_grid_m = pi*lambda_grid_m/180;

[LAMBDA_grid_m, PHI_grid_m] = meshgrid(lambda_grid_m, phi_grid_m);

[X_grid_m, Y_grid_m, Z_grid_m] = DipAndAzimuth_to_xyz(PHI_grid_m, LAMBDA_grid_m, 'lower');

C_plot = zeros(size(X_grid_m));

k = 100;
k = 50;
k = 10;

for p=1:numel(x_p)
  cos_arc_p = cos_arcangle(x_p(p),y_p(p),z_p(p), X_grid_m, Y_grid_m, Z_grid_m);
  C_plot    = C_plot + exp(k*(cos_arc_p-1));
end 

if ~isempty(handle)
  handle = handle+1;
  figure(handle)
end

n_circ = 361;

X_circ = cos(linspace(0, 2*pi, n_circ));
Y_circ = sin(linspace(0, 2*pi, n_circ));

geographical_polar([], [], OPTIONS);

patch(X_circ, Y_circ, 'w')
hold on

% Extend arrays to connect 0 and 360 degrees

X_grid_m = [X_grid_m X_grid_m(:,1)];
Y_grid_m = [Y_grid_m Y_grid_m(:,1)];

C_plot = [C_plot C_plot(:,1)];

if Monochrome  
  contour(X_grid_m, Y_grid_m, C_plot, 11, 'LineColor', 'black')
else
  contour(X_grid_m, Y_grid_m, C_plot, 11)
end

PlotTitleStr = {sprintf('Density plot for site %s', char(SITE.Name_site)), sprintf('Model %s', MODEL.name)}; 

if VelArrow && ~isempty(n_vel_site)
  if abs(norm(n_vel_site)-1)>1.0e-08
    error('NewPlotDensityContours(): Site velocity unit vector must be normalized')
  end  
  R_arrow = 1.25;
  x_1     = R_arrow*n_vel_site(1);
  y_1     = R_arrow*n_vel_site(2);
  x_2     = 1.10*R_arrow*n_vel_site(1);
  y_2     = 1.10*R_arrow*n_vel_site(2);
    
  arrow([x_1 y_1], [x_2 y_2], 'length', 7, 'width', 0.25, 'tipangle', 15, 'baseangle', 60, 'facecolor', [0.6 0.6 0.6])

end

plot(X_circ, Y_circ, 'k')

axis equal
if PlotTitle
  title(PlotTitleStr, 'interpreter', 'none');
end

hold off

handle_out = handle;

 
