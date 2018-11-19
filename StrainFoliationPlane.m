function StrainFoliationPlane(P_dn, L_use, MAP)

R_factor = MAP.R_factor;
Ratio    = MAP.Ratio;
TextSize = MAP.TextSize;

X_c    = P_dn.X(L_use);
Y_c    = P_dn.Y(L_use);

R      = P_dn.R_pre(:,:,L_use);
E_U    = P_dn.E_U(:,L_use);
W_U    = P_dn.vecs_U(:,:,L_use);

% Check to ensure that vecs_U forms a right-handed system of basis vectors
% [THis should actually be done in the top-down tracking runs

J_U   = det(W_U);

if det(W_U)==1
  W_U_rh = W_U;
  W_U_lh = W_U;
  W_U_lh(:,3) = -W_U(:,3);
else  
  W_U_lh = W_U;
  W_U_rh = W_U;
  W_U_rh(:,3) = W_U(:,3);
end

% REMINDER: All operations performed on U are not in the final (R-rotated) system
%           Much of the following is thus nonsense
  
% Foliation plane is defined by the two principle eigenvector directions
% Outward normal to F plane is the 3rd principle eigenvector direction

W_rh   = R*W_U_rh;
W_lh   = R*W_U_lh;

n_rh_a = W_rh(:,1);    % Unit vector aligned with major semi-axis direction
n_rh_b = W_rh(:,2);    % Unit vector aligned with intermediate semi-axis direction
n_rh_c = W_rh(:,3);    % Unit vector aligned with minor semi-axis direction

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

% Equation of foliation plane is   n_F(1) x + n_F(2) y  + n_F(3)*z = 0;
% Intersection of foliation plane with z=0 plane is the strike direction (or we could refer to the plane of local ice surface)

% n_F(1)*x + n_F(2)*y = 0

theta_f_strikedir = atan(-n_rh_c(1)/n_rh_c(2));

n_x_f_strikedir = cos(theta_f_strikedir);
n_y_f_strikedir = sin(theta_f_strikedir);

% Note that "hold on" is already in effect from calling program

R_l    = R_factor;
R_p    = Ratio*R_factor;

plot(X_c, Y_c, 'k')
plot([X_c-R_l*n_x_f_strikedir  X_c+R_l*n_x_f_strikedir],  [Y_c-R_l*n_y_f_strikedir   Y_c+R_l*n_y_f_strikedir], 'k-') 

patch([X_c-0.2*R_p*n_x_f_strikedir  X_c+R_p*n_x_f_dipdir  X_c+0.2*R_p*n_x_f_strikedir], [Y_c-0.2*R_p*n_y_f_strikedir  Y_c+R_p*n_y_f_dipdir  Y_c+0.2*R_p*n_y_f_strikedir], 'k')

n_1 = round(E_U(1)/E_U(3));
n_2 = round(E_U(2)/E_U(3));
n_3 = round(E_U(3)/E_U(3));

e_text = sprintf('%d:%d:%d', n_1, n_2, n_3);

e_TextSize = MAP.EllipsoidTextSize;

R_d   = 2.00*R_p;
R_t   = 1.50*R_p;

text(X_c+n_x_f_dipdir*R_d,  Y_c+n_y_f_dipdir*R_d, sprintf('%.0f', f_dip_deg), 'FontSize', TextSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 180*atan(n_y_f_strikedir/n_x_f_strikedir)/pi)
text(X_c-n_x_f_dipdir*R_t, Y_c-n_y_f_dipdir*R_t, e_text, 'FontSize', e_TextSize, 'HorizontalAlignment', 'center', 'Rotation', 180*atan(n_y_f_strikedir/n_x_f_strikedir)/pi);
