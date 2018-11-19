function [x_p, y_p, z_p, phi, lambda, Ratio_23] = StrainFoliationStereoplotPoints(X_c, Y_c, R, E_U, W_U, METHOD)

% Check to ensure that vecs_U forms a right-handed system of basis vectors
% [THis should actually be done in the top-down tracking runs

% An alternative approach working directly with E_V and W_V is more direct
% except that E_V and W_B not computed in P_dn structure

ns         = numel(X_c);
phi        = zeros(1, ns);
lambda     = zeros(1, ns);
x_DIRECT_p = zeros(1, ns);
y_DIRECT_p = zeros(1, ns);
z_DIRECT_p = zeros(1, ns);
dip_DIRECT = zeros(1, ns);

for s=1:ns
  RR    = R(:,:,s);
  J_U   = det(W_U(:,:,s));

  if abs(J_U-1)<0.00001
    W_U_rh      = W_U(:,:,s);
    W_U_lh      = W_U(:,:,s);
    W_U_lh(:,3) = -W_U(:,3,s);
  else  
    W_U_lh      = W_U(:,:,s);
    W_U_rh      = W_U(:,:,s);
    W_U_rh(:,3) = W_U(:,3,s);
  end

  % REMINDER: All operations performed on U are not in the final (R-rotated) system
  %           Much of the following is thus nonsense
  
  % Foliation plane is defined by the two principle eigenvector directions
  % Outward normal to F plane is the 3rd principle eigenvector direction

  W_rh   = RR*W_U_rh;
  W_lh   = RR*W_U_lh;

  n_rh_a = W_rh(:,1);    % Unit vector aligned with major semi-axis direction
  n_rh_b = W_rh(:,2);    % Unit vector aligned with intermediate semi-axis direction
  n_rh_c = W_rh(:,3);    % Unit vector aligned with minor semi-axis direction
  
  x_DIRECT_p(s) = n_rh_c(1);
  y_DIRECT_p(s) = n_rh_c(2);
  z_DIRECT_p(s) = n_rh_c(3);
  
  if z_DIRECT_p(s)>0
    x_DIRECT_p(s) = -x_DIRECT_p(s);
    y_DIRECT_p(s) = -y_DIRECT_p(s);
    z_DIRECT_p(s) = -z_DIRECT_p(s);
  end  
  
  dip_DIRECT(s) = pi/2 - acos(sqrt(x_DIRECT_p(s)^2 + y_DIRECT_p(s)^2));

  n_lh_a = W_lh(:,1);   % As above but LHS coordinate system
  n_lh_b = W_lh(:,2);
  n_lh_c = W_lh(:,3);

  cos_f_dip = dot(n_rh_c, [0 0 1]');   % f_dip is foliation dip angle

  n_x_f_dipdir = n_rh_c(1)/sqrt(n_rh_c(1)^2 + n_rh_c(2)^2);
  n_y_f_dipdir = n_rh_c(2)/sqrt(n_rh_c(1)^2 + n_rh_c(2)^2);

  if cos_f_dip<0
    cos_f_dip    = -cos_f_dip;
    n_x_f_dipdir = -n_x_f_dipdir;
    n_y_f_dipdir = -n_y_f_dipdir;
  end

  f_dip_deg = 180*acos(cos_f_dip)/pi;
  phi(s)    = acos(cos_f_dip);
  theta     = atan(n_y_f_dipdir/n_x_f_dipdir);
  lambda(s) = pi/2 - theta; 
  Ratio_23(s) = E_U(2,s)/E_U(3,s);
end

switch METHOD
  case 'STRIKE-DIP'
    [x_p, y_p, z_p] = DipAndAzimuth_to_xyz(phi, lambda, 'lower');
  case 'DIRECT'
    lambda = [];
    x_p = x_DIRECT_p;
    y_p = y_DIRECT_p;
    z_p = z_DIRECT_p;
    phi = dip_DIRECT;
  otherwise
    error('StrainFoliationStereoplotPoints(): Unprogrammed case')
end    

