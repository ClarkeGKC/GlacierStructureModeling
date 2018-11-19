function [X_s,Y_s] = ProjectDown(x_p, y_p, z_p)

% Project points (x_p,y_p,z_p) onto lower hemisphere and thence onto lower plane

% Check that all points have been normalized

n_p = sqrt(x_p.^2 + y_p.^2  + z_p.^2);

if abs(n_p-1)>1.e-09
  error('plot_xyz): Input data are not normalized')
end

L_up = z_p>0;

if any(L_up)
  % fprintf(1,'plot_xyz(): Upper hemisphere points encountered and converted to lower hemisphere\n');
  z_p(L_up) = -z_p(L_up);
  x_p(L_up) = -x_p(L_up);
  y_p(L_up) = -y_p(L_up);
end

X_p = x_p;
Y_p = y_p;
Z_p = z_p+1;
R_p = sqrt(X_p.^2+Y_p.^2 + Z_p.^2);

r_xy  = sqrt(x_p.^2 + y_p.^2);
L_xy  = r_xy==0;

X_s  = R_p.*x_p./r_xy;
X_s(L_xy) = R_p(L_xy);

Y_s  = R_p.*y_p./r_xy;
Y_s(L_xy) = R_p(L_xy);

% Note that for equal-area projection the maximum R value is sqrt(2)
% Normalize the projected data to sqrt(2) so that the radius of the plotted data
% circle is 1 so that it matched the contour plot and the radial histogram
% plot in "tryptich"

X_s = X_s/sqrt(2);
Y_s = Y_s/sqrt(2);

