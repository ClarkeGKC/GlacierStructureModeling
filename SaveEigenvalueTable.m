function SaveEigenvalueTable(dir_OUT, SITE, TableString, x_p, y_p, z_p, n_vel_site, POLARITY, VERBOSITY)

if ~exist(dir_OUT, 'dir')
  mkdir(dir_OUT)
end  

if isempty(POLARITY)
  file_tbl = fullfile(dir_OUT, sprintf('EigenvalueTable(%s).dat', char(SITE.Name_site)));
  file_vec = fullfile(dir_OUT, sprintf('EigenpropertiesTable(%s).dat', char(SITE.Name_site)));
else  
  switch POLARITY
    case 'UNIPOLAR'
      file_tbl = fullfile(dir_OUT, sprintf('EigenvalueTable-UNIPOLAR(%s).dat', char(SITE.Name_site)));
      file_vec = fullfile(dir_OUT, sprintf('EigenpropertiesTable-UNIPOLAR(%s).dat', char(SITE.Name_site)));
    case 'BIPOLAR'
      file_tbl = fullfile(dir_OUT, sprintf('EigenvalueTable-BIPOLAR(%s).dat', char(SITE.Name_site)));
      file_vec = fullfile(dir_OUT, sprintf('EigenpropertiesTable-BIPOLAR(%s).dat', char(SITE.Name_site)));
    otherwise
      error('SaveEigenvalueTable(): Unprogrammed value of "POLARITY"=%s', POLARITY)
  end
end

[V, E] = OrientationTensor(x_p,y_p,z_p);
N_pts  = numel(x_p);
   
E_1 = E(1,1);   % smallest
E_2 = E(2,2);
E_3 = E(3,3);   % largest
  
V_1 = V(:,1);
V_2 = V(:,2);
V_3 = V(:,3);

% If the eigenvector is aligned with the plunge direction then it will have
% a NEGATIVE z component. If not then we must change the pointing direction
% to make this true

if V_3(3)>0
  V_3 = -V_3;
  if VERBOSITY>0  
    fprintf(1,'SaveEigenvalueTable(): "z" axis of V_3 is pointing upward. Applying correction to plunge azimuth\n');
  end  
end  

% There is a potential problem here. It seems that Hambrey and I both start with the assumption that
% The first two eigenvalues of the orientation distribution are x- and y-like. Probably this is true but not necessarily bomb-proof
% However the main objective is to determine the spatial orientation of the ellipsoid, This might be better characterized by 
% simply working with the three eigenvectors (only two are needed since they are mutually orthogonal) -- I think we are actually OK (GKCC 2017-05-09)

Theta_xy  = 180*atan2(V_3(2), V_3(1))/pi;
Lambda_xy = mod(90-Theta_xy, 360);

Delta_xy  = 180*atan2(V_3(3), sqrt(V_3(1).^2 + V_3(2).^2))/pi;

sin_Alpha_V2 = norm(cross([0 1 0], V_2));
Alpha_V2     = 180*asin(sin_Alpha_V2)/pi;

if ~isempty(n_vel_site)
  Azim_v_S  = 90-180*atan2(n_vel_site(2), n_vel_site(1))/pi;
  Azim_v_S  = mod(Azim_v_S, 360);
else
  Azim_v_S  = [];
end  

f_tbl = fopen(file_tbl, 'w');    % Append to the eigenvalue and vector file within any given run
f_vec = fopen(file_vec, 'w');


fprintf(f_tbl,'Table 1. Schmidt plot eigenvalues for "%s" sites\n', TableString);
fprintf(f_tbl,'         Notes: Azimuth for V_1 is plunge azimuth\n');
fprintf(f_tbl,'                "Dip" angle for V_1 is plunge angle relative to horizontal\n');
fprintf(f_tbl,'                "Rot(V_2)" angle is anti-clockwise rotation of V2 about axis of V1 relative to "+y" axis\n\n');
 
fprintf(f_tbl, 'Site     S_1     S_2     S_3      Az(V_1)   Dip(V_1)  Rot(V_2)  Az(v_S)\n');
fprintf(f_tbl, '                                   (deg)     (deg)     (deg)     (deg)\n\n');
  
fprintf(f_vec,'Table 2. Eigenproperties for model\n\n');
fprintf(f_vec, 'Site     Vector      E         V_x       V_y       V_z\n\n');  

if ~isempty(Azim_v_S)
  fprintf(f_tbl,'%-4s   %6.4f  %6.4f  %6.4f     %6.1f   %6.1f    %6.1f     %6.1f\n', char(SITE.Name_site), E_3, E_2, E_1, Lambda_xy, Delta_xy, Alpha_V2, Azim_v_S);
else
  fprintf(f_tbl,'%-4s   %6.4f  %6.4f  %6.4f     %6.1f   %6.1f    %6.1f       --\n', char(SITE.Name_site), E_3, E_2, E_1, Lambda_xy, Delta_xy, Alpha_V2);
end
fprintf(f_vec,'%-8s  V_1  %10.6f %10.6f %10.6f %10.6f\n', char(SITE.Name_site), E_3, V_3(1), V_3(2), V_3(3));
fprintf(f_vec,'          V_2  %10.6f %10.6f %10.6f %10.6f\n', E_2, V_2(1), V_2(2), V_2(3));
fprintf(f_vec,'          V_3  %10.6f %10.6f %10.6f %10.6f\n\n', E_1, V_1(1), V_1(2), V_1(3));

fclose(f_tbl);
fclose(f_vec);

