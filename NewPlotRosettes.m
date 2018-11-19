function  handle_out = NewPlotRosettes(handle, MODEL, SITE, x_p, y_p, z_p, n_vel_site, OPTIONS)

[phi, lambda] = xyz_to_DipAndAzimuth(x_p, y_p, z_p, 'upper');

% Plot rosettes showing dip directions

switch OPTIONS.Polarity
  case 'BIPOLAR'
    lambda_minus  = lambda-pi;
    lambda_plus   = lambda+pi;
    L_plus_ok     = lambda_plus<2*pi;
    lambda_comp   = zeros(size(lambda));
    lambda_comp(L_plus_ok)  = lambda_plus(L_plus_ok);
    lambda_comp(~L_plus_ok) = lambda_minus(~L_plus_ok);
    lambda = [lambda lambda_comp];
  case 'UNIPOLAR'
  otherwise
    error('NewPlotRosettes(): Unprogrammed "POLARITY"=%s', POLARITY)
end

n_lambda_bins = 36;   

if ~isempty(handle)
  handle = handle+1;
  figure(handle);
end

PlotTitleStr = {sprintf('Rosette plot for site %s', char(SITE.Name_site)), sprintf('Model %s', MODEL.name)}; 

% rho and theta are polar angles with "theta=0" corresponding to the "x"
% axis and theta is a "counter-clockwise angle"

% Our Rosette plots to "0" to align with true north and measure angles in
% the clockwise direction

if isfield(OPTIONS, 'PlotTitle')
  PlotTitle = OPTIONS.PlotTitle;
else
  PlotTitle = true;
end  

if isfield(OPTIONS, 'dir_PLOT')
  dir_PLOT = OPTIONS.dir_PLOT;
else
  dir_PLOT = {};
end  

if isfield(OPTIONS, 'HardCopy')
  HardCopy = OPTIONS.HardCopy;
else
  HardCopy = false;
end

if isfield(OPTIONS, 'VelArrow')
  VelArrow = OPTIONS.VelArrow;
else
  VelArrow = false;
end  

if isfield(OPTIONS, 'Monochrome')
  Monochrome = OPTIONS.Monochrome;
else
  Monochrome = false;
end  
  
[theta, rho] = rose(lambda, n_lambda_bins);

v    = geographical_polar(theta, rho, OPTIONS);

rmin = 0;
rmax = v(4);

if VelArrow && ~isempty(n_vel_site)
  if abs(norm(n_vel_site)-1)>1.0e-08
    error('NewPlotRosettes(): Site velocity unit vector must be normalized')
  end  
  R_arrow = 1.25*rmax;
  x_1     = R_arrow*n_vel_site(1);
  y_1     = R_arrow*n_vel_site(2);
  x_2     = 1.10*R_arrow*n_vel_site(1);
  y_2     = 1.10*R_arrow*n_vel_site(2);
    
  arrow([x_1 y_1], [x_2 y_2], 'length', 7, 'width', 0.25, 'tipangle', 15, 'baseangle', 60, 'facecolor', [0.6 0.6 0.6])

end

if PlotTitle
  title(PlotTitleStr, 'interpreter', 'none');
end

if HardCopy && ~isempty(handle)
  if isempty(dir_PLOT)
    dir_PLOT = pwd;
  end  
  if ~exist(fullfile(dir_PLOT, 'RosettePlots'), 'dir')
    mkdir(fullfile(dir_PLOT, 'RosettePlots'))
  end  
  print(handle, fullfile(dir_PLOT, 'RosettePlots', sprintf('Fig_%s.pdf', char(SITE.Name_site))), '-dpdf')
end  

handle_out = handle;
