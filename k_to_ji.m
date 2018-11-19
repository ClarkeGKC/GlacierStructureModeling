function [j i]=k_to_ji(k, ny)

i = 1+floor((k-1)/ny);
j = k-(i-1)*ny;

