function P_dn = Compute_F_tensor_etc(P_dn)

global dt

ns = numel(P_dn);

for s=1:ns
  yr     = P_dn(s).t;
  L      = P_dn(s).L;

  nt     = numel(yr);
  
  F      = zeros(3,3,nt);
  G      = zeros(3,3,nt);
  R_pre  = NaN(3,3,nt);
  R_post = NaN(3,3,nt);
  U      = NaN(3,3,nt);
  V      = NaN(3,3,nt);
    
  E_U    = NaN(3,nt);
  
  vecs_U   = NaN(3,3,nt);  
  R_vecs_U = NaN(3,3,nt);
  
  F_det   = NaN(nt, 1);
      
  for t=1:nt
    
    A      = eye(3) - 0.5*dt*squeeze(L(:,:,t));
    B      = eye(3) + 0.5*dt*squeeze(L(:,:,t));
    if t==1
      F(:,:,1) = eye(3);
    else
      F(:,:,t) = A\B*squeeze(F(:,:,t-1));  % Same as  F = A_inv*B*squeeze(F(:,:,t-1));
    end  
    
    F_det(t) = det(F(:,:,t));    
 
    G(:,:,t) = 0.5*(F(:,:,t)'*F(:,:,t) - eye(3));  % Green strain tensor
    
    % R_pre is the rotation from point of deposition to point along track

    [R_pre(:,:,t), U(:,:,t), V(:,:,t)] = poldecomp(F(:,:,t));    
        
    % Eigenproperties of U
              
    [VVecs_U, EE_U] = eig(U(:,:,t), 'vector');
    
    % Sort the eigenvalues and eigenvectors so that E1>E2>E3 
    
    [~, i_sort]   = sort(EE_U, 'descend');
    
    E_U(:,t)      = EE_U(i_sort);
    vecs_U(:,:,t) = VVecs_U(:,i_sort);
    R_vecs_U(:,:,t) = R_pre(:,:,t)*vecs_U(:,:,t);       
  end
  
  % R_post is the rotation from point along track to point of arrival at sampling site
    
  for t=1:nt
    R_post(:,:,t) = R_pre(:,:,nt)*squeeze(R_pre(:,:,t))';
  end
  
  P_dn(s).F      = F;
  P_dn(s).det_F  = F_det;
  P_dn(s).R_pre  = R_pre;
  P_dn(s).R_post = R_post;
  P_dn(s).U      = U;
  P_dn(s).V      = V;
  
  P_dn(s).E_U      = E_U;  
  P_dn(s).vecs_U   = vecs_U;
  P_dn(s).R_vecs_U = R_vecs_U;  % It is not clear whether this is a useful quantity
  
end  