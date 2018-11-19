function [V_1, V_2, V_3, E_1, E_2] = DecomposeCrackGeneratingTensor(D)

[~,~,nt] = size(D);

D_2x2    = D(1:2,1:2,:);
V_2x2    = zeros(2,2, nt);
E_2x2    = zeros(2,nt);

for t=1:nt
  [V2, E2] = eig(D_2x2(:,:,t), 'vector');
  [E_2x2(:,t), i_sort] = sort(E2, 'descend');
  V_2x2(:,:,t) = V2(:,i_sort);
end
           
E_1        = E_2x2(1,:);  % Principal 2D eigenvalue
E_2        = E_2x2(2,:);  %

V_1        = zeros(3,nt);  % Principal direction
V_1(1:2,:) = squeeze(V_2x2(:,1,:));
    
V_2        = zeros(3,nt);
V_2(1:2,:) = squeeze(V_2x2(:,2,:)); 
V_3        = cross(V_1, V_2, 1);