function   [V, E] = OrientationTensor(n_x,n_y,n_z)

N    = numel(n_x);
n_11 = sum(n_x.^2)/N;
n_12 = sum(n_x.*n_y)/N;
n_13 = sum(n_x.*n_z)/N;
n_22 = sum(n_y.^2)./N;
n_23 = sum(n_y.*n_z)/N;
n_33 = sum(n_z.^2)/N;

A = [n_11 n_12 n_13;
     n_12 n_22 n_23;
     n_13 n_23 n_33];
   
[V, E] = eig(A);   
   
   
   
     