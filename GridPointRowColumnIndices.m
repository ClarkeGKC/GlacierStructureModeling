function  [x_grid, y_grid, index] = GridPointRowColIndices(X_c_pts, Y_c_pts)

x_grid = sort(unique(X_c_pts), 'ascend');
y_grid = sort(unique(Y_c_pts), 'descend');

nx     = numel(x_grid);
ny     = numel(y_grid);

N_xy   = nx*ny;

if N_xy~=numel(X_c_pts)
  error('GridPoints(): Product of nx*ny should equal number of points. Probable rounding trouble.')
end

x_min = x_grid(1);
y_max = y_grid(1);

dx    = (x_grid(nx)-x_grid(1))/(nx-1);
dy    = (y_grid(1)-y_grid(ny))/(ny-1);

i_pts = 1+round((X_c_pts-x_min)/dx);
j_pts = 1+round((y_max-Y_c_pts)/dy);

index = zeros(ny, nx);

for k=1:N_xy
  index(j_pts(k), i_pts(k)) = k;
end  

index = reshape(index, 1, N_xy);


