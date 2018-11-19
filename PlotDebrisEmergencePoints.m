function  handle_out = PlotDebrisEmergencePoints(handle, MODEL, SOURCE, DD, dir_Moraine, HARDCOPY)

dx            = MODEL.dx;
dy            = MODEL.dy;
nx            = MODEL.nx;
ny            = MODEL.ny;
x_grid        = MODEL.x_grid;
y_grid        = MODEL.y_grid;

year_cal      = DD.cal_year;
year_model    = DD.t;

H             = zeros(ny, nx);
H(DD.ic_jc_I) = DD.H;

c_map       = colormap;
[nr, ~ ]    = size(c_map);
c_map(1,:)  = [0 0 0];
c_map(nr,:) = [1 1 1];
colormap(c_map)

X_min = x_grid(1)-0.5*dx;
X_max = x_grid(end)+0.5*dx;
Y_min = y_grid(end)-0.5*dy;
Y_max = y_grid(1)+0.5*dy;

X_e   = [SOURCE.X_e];
Y_e   = [SOURCE.Y_e];
t_e   = [SOURCE.t_e];

handle = handle+1;
  
figure(handle)
plot([X_min X_max X_max X_min X_min], [Y_min Y_min Y_max Y_max Y_min], 'k')
hold on

imagesc(x_grid, y_grid, H>0) 
set(gca, 'YDir', 'normal')

plot(X_e, Y_e, 'k.')
plot([SOURCE.X_source], [SOURCE.Y_source], 'r+')
hold off
title(sprintf('Model "%s"', MODEL.name), 'interpreter', 'none')

if HARDCOPY
  print(handle, fullfile(dir_Moraine, 'DebrisEmergencePointMap.pdf'), '-dpdf')
end

handle_out = handle;


