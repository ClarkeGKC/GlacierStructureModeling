function [n_x,n_y,n_z] = DipAndAzimuth_to_xyz(phi, lambda, hemisphere)

[nr, nc] = size(phi);
phi      = reshape(phi, numel(phi), 1);
lambda   = reshape(lambda, numel(lambda), 1);

if any(phi>pi/2) || any(lambda<0 | lambda>2*pi)
  error('DipAndAzimuth_to_xyz(): phi and/or lambda lie outside expected range')
end  
  
n_z = cos(phi);

theta = pi/2-lambda;

n_x = cos(theta).*sin(phi);
n_y = sin(theta).*sin(phi);

switch hemisphere
  case 'lower'    % Finally convert to lower hemisphere projection if n_z is positive
    L_up = n_z>0;

    n_z(L_up) = -n_z(L_up);
    n_x(L_up) = -n_x(L_up);
    n_y(L_up) = -n_y(L_up);
  case 'upper'   % Convert to upper hemisphere if n_z is negative
    L_dn = n_z<0;
    n_z(L_dn) = -n_z(L_dn);
    n_x(L_dn) = -n_x(L_dn);
    n_y(L_dn) = -n_y(L_dn);
  case 'full'
    
  otherwise
    error('DipAndAzimuth_to_xyz(): Unprogrammed value of "hemisphere"')
end    

n_x = reshape(n_x, nr, nc);
n_y = reshape(n_y, nr, nc);
n_z = reshape(n_z, nr, nc);