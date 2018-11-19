function  [L_crack, theta_p, c_type_p] = ApplyCrackRule(CRACK_MODEL, P_dn)

% Note that n_x, n_y are components of the unit vector aligned with 2D flow direction
% and that nc_x, nc_y are the components of unit vector perpendicular to
% the crack plane (or failing a crack they are aligned with the maximum
% tensile stress direction)

% The instantaneous flow direction (or the direction of maximum downslope -- which is
% identical is required for the binning of crack directions into
% longitudinal, transverse and diagonal as in Hambrey's analysis

% Note that theta_p and c_type_p refer to the on-track orientations not to the final orientations at the observation site

global RHO RHO_w g VERBOSITY
 
DEBUG = 0;

% GKCC has changed his point of view on the following. With the new (x-y)
% based stress analysis, e.g., using the hydrostatic approximation and
% reducing stress tensor to 2d, the correct angle comparison is between
% the 2D flow direction and the 2D crack plane (which is always vertical when it forms)

vel_2d_mag = sqrt(P_dn.u.^2 + P_dn.v.^2);

n_x = P_dn.u./vel_2d_mag;   % 2D flow directions (not used here but could be)
n_y = P_dn.v./vel_2d_mag;

% n_z = zeros(size(n_x));

nt  = numel(P_dn.t);

switch CRACK_MODEL.RULE
  case 1    
    
    [V_1, V_2, V_3, E_1, E_2] = DecomposeCrackGeneratingTensor(P_dn.D);
    
  case 2
    
    [V_1, V_2, V_3, E_1, E_2] = DecomposeCrackGeneratingTensor(P_dn.sigma);
   
  case {3,4,5,6,7,8}
    
    R      = P_dn.sigma;
    R_plus = kron(reshape(eye(3), numel(eye(3)), 1), P_dn.sigma_zz);
    R      = R - reshape(R_plus, 3, 3, nt);
    
   [V_1, V_2, V_3, E_1, E_2] = DecomposeCrackGeneratingTensor(R);
    
  case 9

    [V_1, V_2, V_3, E_1, E_2] = DecomposeCrackGeneratingTensor(P_dn.sigma);
         
  otherwise
    error('ApplyCrackRule(): Unprogrammed RULE')
end

% outward normal to crack plane 
   
nc_x         = squeeze(V_1(1,:,:));
nc_y         = squeeze(V_1(2,:,:));        
nc_z         = squeeze(V_1(3,:,:));

% 2D outward normal (in map plane) of the above crack plane.

nc_x         = nc_x./sqrt(nc_x.^2 + nc_y.^2);
nc_y         = nc_y./sqrt(nc_x.^2 + nc_y.^2);

% Angle between ice flow direction at point along the flow trajectory and the potential orientation of the crack plane direction at the same site

% NOTE: This differs from Hambrey's field measurements which compare the ice flow direction (inferred) at the measuirement site with the
%       the evolved orientations of crack planes at that sight.

% Here  I would expect the orientation of the cracks and flow to be either co-aligned or orthogonal 

% It does not seem to matter whether nc is 2D or 3D because nc_z=0 with
% the crack models that I use.

% cos_theta_p = n_x.*nc_x + n_y.*nc_y + n_z.*nc_z;
cos_theta_p = n_x.*nc_x + n_y.*nc_y;

theta_p = 180*acos(cos_theta_p)/pi;   % Maps into region 0-180 deg (0-pi rad)

% convert all angles in quadrant 2 to angles in quadrant 4

theta_p(theta_p>90) = theta_p(theta_p>90)-180;

switch CRACK_MODEL.RULE
  case 1
    D_eig_max = E_1;
    L_crack   = D_eig_max>CRACK_MODEL.D_Ext_threshold;
  case 2   
    sigma_eig_max = E_1;
    L_crack       = sigma_eig_max>CRACK_MODEL.SIGMA_Ext_threshold;
  case 3 
    R_xx    = E_1;
    Depth   = P_dn.S - P_dn.Z;
    L_neg   = Depth<0;
    if any(L_neg)
      fprintf(1,'ApplyCrackRule(): WARNING - negative particle depth encountered for site %s and CRACK_MODEL.RULE=%d   ******\n', char(P_dn.Name_site), CRACK_MODEL.RULE);
    end  
    K_I     = 1.12*R_xx.*sqrt(pi*Depth) - 0.683*RHO*g*Depth.^1.5;
    L_crack = K_I>CRACK_MODEL.K_threshold;
  case 4
    R_xx    = E_1; 
    H       = P_dn.H;
    Z       = P_dn.Z;
    Depth   = P_dn.S - P_dn.Z;   
    L_neg   = Depth<0;
    if any(L_neg)
      fprintf(1,'ApplyCrackRule(): WARNING - negative particle depth encountered for site %s and CRACK_MODEL.RULE=%d   ******\n', char(P_dn.Name_site), CRACK_MODEL.RULE);
    end      
    h_w     = CRACK_MODEL.f_0*RHO*H/RHO_w;
    z_w     = P_dn.B + h_w;   
    b_w     = max(z_w-Z, 0);
      
    K_I     = 1.12*R_xx.*sqrt(pi*Depth) - 0.683*RHO*g*Depth.^1.5 + 0.683*RHO_w*g*b_w.^1.5;
    L_crack = K_I>CRACK_MODEL.K_threshold;
  case 5
    R_xx    = E_1;
    M       = P_dn.M;
    f       = CRACK_MODEL.f_min + (CRACK_MODEL.f_max-CRACK_MODEL.f_min)*M;
    H       = P_dn.H;
    Z       = P_dn.Z;
    Depth   = P_dn.S - P_dn.Z;   
    L_neg   = Depth<0;
    if any(L_neg)
      fprintf(1,'ApplyCrackRule(): WARNING - negative particle depth encountered for site %s and CRACK_MODEL.RULE=%d   ******\n', char(P_dn.Name_site), CRACK_MODEL.RULE);
    end     
    h_w     = (RHO/RHO_w)*f.*H;
    z_w     = P_dn.B + h_w;
    b_w     = max(z_w-Z, 0);
    
    K_I     = 1.12*R_xx.*sqrt(pi*Depth) - 0.683*RHO*g*Depth.^1.5 + 0.683*RHO_w*g*b_w.^1.5;
    L_crack = K_I>CRACK_MODEL.K_threshold;
  case 6
    R_xx    = E_1;
    Depth   = P_dn.S -P_dn.Z; 
    b_term  = 0.5*CRACK_MODEL.spacing;
    S_term  = b_term + Depth;
    b_S     = b_term./S_term;
    
    F_1     = 1 + 1/2*b_S + 3/8*b_S.^2 + 5/16*b_S.^3 + 35/128*b_S.^4 + 63/256*b_S.^5 + 231/1024*b_S.^6;
    F_2     = 22.501*b_S.^7 -63.502*b_S.^8 + 58.045*b_S.^9 -17.577*b_S.^10;
    F       = F_1/sqrt(pi) + F_2;
    
    K_I     = F.*sqrt(b_S).*R_xx.*sqrt(pi*Depth);
    L_crack = K_I>CRACK_MODEL.K_threshold; 
  case 7
    error('ApplyCrackRule(7): RULE %d unprogrammed', RULE)    
  case 8
    error('ApplyCrackRule(8): RULE %d unprogrammed', RULE)    
  case 9
    sigma_xx = E_1;
    Depth    = P_dn.S - P_dn.Z;
    L_crack  = sigma_xx>0; 
    
    if VERBOSITY>0
      fprintf(1,'\n\nCRACK RESULTS FOR SITE "%s" TRAJECTORY\n', char(P_dn.Name_site));
      fprintf(1,'  Total points along trajectory                        = %d\n', nt);
      fprintf(1,'  Total points passing fracture test (before blanking) = %d\n', sum(L_crack));
          
      figure(900)

      subplot(2,1,1)
      plot(P_dn.t, P_dn.XI, 'k:', P_dn.t(L_crack), P_dn.XI(L_crack), 'r.')
      txL = text(P_dn.t(1), P_dn.XI(1)-0.01, char(P_dn.Name_site));
      txL.Rotation = 60;
      txL.HorizontalAlignment = 'right';
      hold on
      
      subplot(2,1,2)
      plot(P_dn.X, P_dn.XI, 'k:', P_dn.X(L_crack), P_dn.XI(L_crack), 'r.')
      txL = text(P_dn.X(1), P_dn.XI(1)-0.01, char(P_dn.Name_site));
      txL.Rotation = 60;
      txL.HorizontalAlignment = 'right';
      txR = text(P_dn.X(nt), P_dn.XI(nt)-0.01, char(P_dn.Name_site));
      txR.Rotation = -60;

      hold on
      
      figure(901)
      plot(P_dn.X, P_dn.Y, 'k'), axis equal
      txL = text(P_dn.X(1)-50, P_dn.Y(1), char(P_dn.Name_site));
      txl.HorizontalAlignment = 'right';
      txR = text(P_dn.X(nt)+10, P_dn.Y(nt), char(P_dn.Name_site));
      txR.HorizontalAlignment = 'left';
      hold on
    end
  otherwise
    error('ApplyCrackRule(): Unprogrammed RULE')
end

c_type_p      = cell(1, nt);
[c_type_p{:}] = deal(' ');

switch CRACK_MODEL.BIN_DEF
  case 'Hambrey'
    theta = [20 90-20];  % Consistent with Hambrey's transverse and longitudinal bins 
  case 'Equal'
    theta = [22.5 90-22.5];  % All bins have equal angle span
  otherwise
    error('main_crack_density_plots_and_table(): Unprogrammed CRACK_BIN_DEF')
end

% Longitudinal crack plane is parallel to flow and the normal to crack
% plane is orthogonal to flow vector
% Transverse crack is perpendicular to flow

% These orientations apply AT the crack formation point and will differ
% from the advected and rotated orientations at the measurements sites

L_trans = L_crack & abs(theta_p)<theta(1);
L_long  = L_crack & abs(theta_p)>theta(2);
L_diag  = L_crack & abs(theta_p)>=theta(1) & abs(theta_p)<=theta(2);

c_type_p(L_trans) = {'T'};
c_type_p(L_long)  = {'L'};
c_type_p(L_diag)  = {'D'};
