function handle_out = PlotSchmidt_S1(handle, MODEL, SITE, t_p, x_p, y_p, z_p, Color_p, n_vel_site, OPTIONS)

% (x_p, y_p, z_p) are coordinates of the unit vector n_p that is
% associated with the stereographic projection of crack "p"

if isfield(OPTIONS, 'dir_PLOT')
  dir_PLOT = OPTIONS.dir_PLOT;
else
  dir_PLOT = [];
end  

if isfield(OPTIONS, 'ShowNet')
  ShowNet = OPTIONS.ShowNet;
else
  ShowNet = true;
end  

if isfield(OPTIONS, 'PlotTitle')
  PlotTitle = OPTIONS.PlotTitle;
else
  PlotTitle = true;
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

if isfield(OPTIONS, 'MarkerScale')
  MarkerScale = OPTIONS.MarkerScale;
else
  MarkerScale = 1.0;
end

if isfield(OPTIONS, 'SurgeMarkers')
  SurgeMarkers = OPTIONS.SurgeMarkers;
else
  SurgeMarkers = 1;
end  
  
phi_grid      = linspace(0, 90, 10);
lambda_grid   = linspace(0, 360, 25);
lambda_fine_grid = linspace(0, 360, 361);

if isempty(t_p)
  HAMBREY_DATA = true;
  N_pts        = numel(x_p);
else
  HAMBREY_DATA = false;
end

if ~HAMBREY_DATA
  t_p_min      = min(t_p);
  t_p_max      = max(t_p);    % SITE.yr_site

  N_pts        = numel(t_p);

  t_last_stop  = SITE.yr_last_surge.stop;
  t_last_start = SITE.yr_last_surge.start;

  if numel(MODEL.t_CYCLE)>1
    t_cycle    = MODEL.t_CYCLE(1);
    dt_SURGE   = MODEL.dt_SURGE(1);
  else  
    t_cycle      = MODEL.t_CYCLE;
    dt_SURGE     = MODEL.dt_SURGE;
  end

  if t_last_stop>=t_p_min
    t_stop     = t_last_stop:-t_cycle:t_p_min;
  else
    t_stop     = [];
  end

  if t_last_start>=t_p_min
    t_start    = t_last_start:-t_cycle:t_p_min;
  else
    t_start    = [];
  end

  if t_stop(end)-dt_SURGE<=t_p_min
    t_start(end+1) = t_p_min;
  end

  n_surge     = numel(t_start);

  if numel(t_start)~=numel(t_stop)
    error('NewPlotSchmidt(): Tilt at Site %s', char(SITE.Name_site));
  end  

  t_stop  = t_stop(end:-1:1);
  t_start = t_start(end:-1:1);

  pen_nonsurge   = {'k+'};  
  marker_nonsurge = {'+'};
  color_nonsurge   = {'k'};
  
  if SurgeMarkers
    pen_surge = {'r+', 'g+', 'b+', 'm+', 'c+', 'y+'};  % Surges (most recent to oldest)
    color_surge  = {'r', 'g', 'b', 'm', 'c', 'y'}; 
    marker_surge = {'+', '+', '+', '+', '+', '+'};
  else
    pen_surge    = {'k.', 'k.', 'k.', 'k.', 'k.', 'k.'};
    color_surge  = {'k', 'k', 'k', 'k', 'k', 'k'};
    marker_surge = {'.',  '.', '.', '.', '.', '.'};
  end
  
  pen_p    = repmat(pen_nonsurge, 1, numel(x_p));         % Non-surge marker
  color_p  = repmat(color_nonsurge, 1, numel(x_p));
  marker_p = repmat(marker_nonsurge, 1, numel(x_p));
  
  np       = numel(color_surge);

  for n=n_surge:-1:1
    p_use = 1+mod(n_surge-n, np);
    L_p = t_p>=t_start(n) & t_p<=t_stop(n);
    pen_p(L_p) = pen_surge(p_use);
    marker_p(L_p) = marker_surge(p_use);
    color_p(L_p)  = color_surge(p_use);
  end
end

pen_default    = {'k.'};
marker_default = 'o';
marker_face_color_default = 'k';
marker_size_default       = MarkerScale*5.0;

PlotTitleStr = {sprintf('Schmidt plot for site %s', char(SITE.Name_site)), sprintf('Model %s', MODEL.name)}; 

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

% Placing arrow at this point of the program suppresses complaint about
% axis resizing

if VelArrow && ~isempty(n_vel_site)
  if abs(norm(n_vel_site)-1)>1.0e-08
    error('NewPlotSchmidt(): Site velocity unit vector must be normalized')
  end  
  R_arrow = 1.25;
  x_1     = R_arrow*n_vel_site(1);
  y_1     = R_arrow*n_vel_site(2);
  x_2     = 1.10*R_arrow*n_vel_site(1);
  y_2     = 1.10*R_arrow*n_vel_site(2);
    
  arrow([x_1 y_1], [x_2 y_2], 'length', 7, 'width', 0.25, 'tipangle', 15, 'baseangle', 60, 'facecolor', [0.6 0.6 0.6])

end

[x_grid,y_grid,z_grid] = DipAndAzimuth_to_xyz(pi/2*ones(1,numel(lambda_fine_grid)), pi*lambda_fine_grid/180, 'upper');
[X_s,Y_s] = ProjectDown(x_grid, y_grid, z_grid);
plot(X_s, Y_s, 'k')
axis equal
hold on

[X_s,Y_s] = ProjectDown(x_p, y_p, z_p);

if ~HAMBREY_DATA
  LegendStr = {};

  pen_log   = {};
  
  s_plot  = [];
  s_num   = 0;
  non_s_init = 1;

  for i=1:numel(X_s)
    if non_s_init && strcmp(color_p{i}, color_nonsurge) && strcmp(marker_p{i}, marker_nonsurge)
      non_s_init = 0;
      s_plot(end+1) = plot(X_s(i), Y_s(i), 'LineStyle', 'none', 'MarkerFaceColor', char(color_p(i)), 'Marker', char(marker_p(i)));
      LegendStr{end+1} = sprintf('Non-surge');
    elseif ~strcmp(color_p{i}, color_nonsurge) && ~strcmp(marker_p{i}, marker_nonsurge) && ~ismember(pen_p{i}, pen_log); 
      s_num  = s_num+1;
      s_plot(end+1) = plot(X_s(i), Y_s(i), 'LineStyle', 'none', 'MarkerFaceColor', char(color_p(i)), 'Marker', char(marker_p(i)));
      pen_log{end+1} = strcat(color_p{i}, marker_p{i});
      LegendStr{end+1} = sprintf('%04d-%04d', round(t_start(s_num)), round(t_stop(s_num)));    
    else  
      plot(X_s(i), Y_s(i), 'LineStyle', 'none', 'MarkerFaceColor', char(color_p(i)), 'Marker', char(marker_p(i)), 'MarkerSize', marker_size_default);
    end  
    hold on
  end
else
  plot(X_s(Color_p==1), Y_s(Color_p==1), 'LineStyle', 'none', 'MarkerFaceColor', OPTIONS.Color_1, 'MarkerEdgeColor', marker_face_color_default, ...
    'Marker', marker_default, 'MarkerSize', marker_size_default)
  plot(X_s(Color_p==2), Y_s(Color_p==2), 'LineStyle', 'none', 'MarkerFaceColor', OPTIONS.Color_2, 'MarkerEdgeColor', marker_face_color_default, ...
    'Marker', marker_default, 'MarkerSize', marker_size_default)
  plot(X_s(Color_p==3), Y_s(Color_p==3), 'LineStyle', 'none', 'MarkerFaceColor', OPTIONS.Color_3, 'MarkerEdgeColor', marker_face_color_default, ...
    'Marker', marker_default, 'MarkerSize', marker_size_default)
  plot(X_s(Color_p==4), Y_s(Color_p==4), 'LineStyle', 'none', 'MarkerFaceColor', OPTIONS.Color_4, 'MarkerEdgeColor', marker_face_color_default, ...
    'Marker', marker_default, 'MarkerSize', marker_size_default)
end

if ~HAMBREY_DATA && SurgeMarkers
  legend(s_plot, LegendStr, 'location', 'southeast');
  legend('boxoff')
end

% axis equal

if PlotTitle
  title(PlotTitleStr, 'interpreter', 'none');
end

hold off

if ~isempty(handle) && ~isempty(dir_PLOT)
  if ~exist(fullfile(dir_PLOT, 'SchmidtPlots'), 'dir')
    mkdir(fullfile(dir_PLOT, 'SchmidtPlots'))
  end  
  print(handle, fullfile(dir_PLOT, 'SchmidtPlots', sprintf('Fig_%s.pdf', char(SITE.Name_site))), '-dpdf')
end

handle_out = handle;
