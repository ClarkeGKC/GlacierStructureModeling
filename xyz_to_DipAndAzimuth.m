function [phi, lambda, phi_deg, lambda_deg] = xyz_to_DipAndAzimuth(x_p, y_p, z_p, hemisphere)

% Script to convert from unit vector (x_p,y_p,z_p) representation to (phi, lambda)
% Noting that phi is angle from vertical and lambda is clockwise rotation from N

% Check that all points have been normalized

n_p = sqrt(x_p.^2 + y_p.^2  + z_p.^2);

if abs(n_p-1)>1.e-09
  error('plot_xyz): Input data are not normalized')
end

switch lower(hemisphere)
  case 'lower'     % convert all data to lower hemisphere
    L_up = z_p>0;
    
    z_p(L_up) = -z_p(L_up);
    x_p(L_up) = -x_p(L_up);
    y_p(L_up) = -y_p(L_up);
    
  case 'upper'    % Convert all data to upper hemisphere 
    L_dn = z_p<0;

    z_p(L_dn) = -z_p(L_dn);
    x_p(L_dn) = -x_p(L_dn);
    y_p(L_dn) = -y_p(L_dn);
  otherwise
    error('xyz_to_DipAndAzimuth(): Unprogrammed value of "hemisphere"')
end

theta  = atan2(y_p, x_p);
lambda = pi/2-theta;

lambda = mod(lambda, 2*pi);

phi    = acos(z_p);

lambda_deg = 180*lambda/pi;
phi_deg    = 180*phi/pi;


