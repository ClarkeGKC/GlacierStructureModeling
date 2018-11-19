function P_dn = Compute_D_tensor_etc(P_dn)

global RHO g YRSEC

ns = numel(P_dn);

for ss=1:ns
  yr     = P_dn(ss).t;
  L      = P_dn(ss).L;
  ETA    = P_dn(ss).ETA;
  s      = P_dn(ss).s;
  
  nt     = numel(yr);
  
  D      = 0.5*(L + permute(L, [2 1 3]));
  W      = 0.5*(L - permute(L, [2 1 3]));
  
  E_D    = NaN(3,nt);       % eigenproperties of the deformation rate tensor
  E_R    = NaN(3,nt);
  vecs_D = NaN(3,3,nt);
    
  sigma_zz = -RHO*g*(P_dn(ss).S - P_dn(ss).Z);   % glaciostatic pressure  
  sigma    = P_dn(ss).s;                         % Start by setting full-stress tensor to deviatoric tensor then insert diagonal values
  
  sigma(1,1,:) = 2*squeeze(s(1,1,:))' + squeeze(s(2,2,:))' + sigma_zz;
  sigma(2,2,:) = 2*squeeze(s(2,2,:))' + squeeze(s(1,1,:))' + sigma_zz;
  sigma(3,3,:) = sigma_zz;
  
  P_dn(ss).sigma_zz  = sigma_zz;
  
  p                  = -sigma_zz - squeeze(s(1,1,:))' - squeeze(s(2,2,:))';   % e.g., Greve & Blatter, p. 73
  P_dn(ss).p         = p;
  P_dn(ss).sigma     = sigma;
  
  R_zz               = 0;  
  R                  = s;  % Start by setting resistive stress tensor to deviatoric stress then insert diagonal values
  
  R(1,1,:)           = 2*s(1,1,:) + s(2,2,:) + R_zz;
  R(2,2,:)           = 2*s(2,2,:) + s(1,1,:) + R_zz;
  R(3,3,:)           = R_zz;
       
  for t=1:nt
            
    % Eigenproperties of D (strain rate)
              
    [VVecs_D, EE_D] = eig(D(:,:,t), 'vector');
    
    % Sort the eigenvalues and eigenvectors so that E1>E2>E3 
    
    [~, i_sort]   = sort(EE_D, 'descend');    
    E_D(:,t)      = EE_D(i_sort);
    vecs_D(:,:,t) = VVecs_D(:,i_sort);

    [~, EE_R]     = eig(R(:,:,t), 'vector');
    [~, i_sort]   = sort(EE_R, 'descend');
    E_R(:,t)      = EE_R(i_sort);
  end
  
  % Units of D are 1/yr
  % Must convert to 1/sec for stress calculations from viscosity
  
  P_dn(ss).D       = D;
  P_dn(ss).W       = W;
  P_dn(ss).E_D     = E_D;
  P_dn(ss).vecs_D  = vecs_D;  
  E_s              = 2*repmat(ETA, 3, 1).*E_D/YRSEC;
  E_sigma          = E_s-p;
  
  P_dn(ss).E_s     = E_s;
  P_dn(ss).E_sigma = E_sigma;
  P_dn(ss).E_R     = E_R;           % No algebraic shortcut here
end  