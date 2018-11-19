function [x_ELA, y_ELA] = GetELA(x_grid, y_grid, b_dot)

BW = zeros(size(b_dot));
BW(b_dot>=0) = 1;

B = bwboundaries(BW, 8);

x_ELA = [];
y_ELA = [];

[nr_x, nc_x] = size(x_grid);
[nr_y, nc_y] = size(y_grid);

if nc_x==1
  x_grid = x_grid';
end

if nc_y==1
  y_grid = y_grid';
end  

for k=1:length(B)
  rc = B{k};
  r  = rc(:,1);
  c  = rc(:,2);
  x_add = x_grid(c);
  y_add = y_grid(r);
  x_ELA = [x_ELA x_add NaN];
  y_ELA = [y_ELA y_add NaN];
end



