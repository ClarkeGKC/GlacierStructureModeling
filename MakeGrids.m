function [XX, YY, ZZ] = MakeGrids(MODEL)

x_grid = MODEL.x_grid;
y_grid = MODEL.y_grid;
XI     = MODEL.XI;
n_XI   = MODEL.n_XI;

nx     = MODEL.nx;
ny     = MODEL.ny;

x_grid_IC = 0.5*(x_grid(1:nx-1)+x_grid(2:nx));
y_grid_JC = 0.5*(y_grid(1:ny-1)+y_grid(2:ny));

[XX.IC_jc, YY.IC_jc] = meshgrid(x_grid_IC, y_grid);
[XX.ic_JC, YY.ic_JC] = meshgrid(x_grid, y_grid_JC);

[XX.ic_jc, YY.ic_jc] = meshgrid(x_grid, y_grid);
[YY.ic_jc_kc, ZZ.ic_jc_kc, XX.ic_jc_kc] = meshgrid(y_grid, XI, x_grid);

XX.ic_jc_kc = reshape(XX.ic_jc_kc, n_XI, ny, nx);
YY.ic_jc_kc = reshape(YY.ic_jc_kc, n_XI, ny, nx);
ZZ.ic_jc_kc = reshape(ZZ.ic_jc_kc, n_XI, ny, nx);


