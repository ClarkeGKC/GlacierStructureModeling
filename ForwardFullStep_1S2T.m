function [X_p, Y_p, Z_p, XI_p, L_ok_p] = ForwardFullStep_1S2T(X, Y, Z, B, H, b_dot, u, v, w, t, L_ok)

% Note "Forward" steps are downflow

global dt X2 Y2 X3 Y3 Z3 

d_max       = 0.0001;  % maximum distance error for convergence

X_p         = NaN(size(L_ok));
Y_p         = NaN(size(L_ok));
Z_p         = NaN(size(L_ok));
XI_p        = NaN(size(L_ok));

X_0         = X(L_ok);
Y_0         = Y(L_ok);
Z_0         = Z(L_ok);

B_0         = interp2(X2, Y2, B, X_0, Y_0);
H_0         = interp2(X2, Y2, H, X_0, Y_0);
XI_0        = min((Z_0-B_0)./H_0, 1);

counter     = 1:numel(L_ok);

i_alive     = counter(L_ok);

% This was added to prevent NaN velocity values from being encountered
% while exploring potential step endpoints

u(isnan(u)) = 0;
v(isnan(v)) = 0;
w(isnan(w)) = 0;

cnt = 0;

a   = 0.5*dt*interp3(Y3, Z3, X3, u, Y_0, XI_0, X_0);
b   = 0.5*dt*interp3(Y3, Z3, X3, v, Y_0, XI_0, X_0);
c   = 0.5*dt*interp3(Y3, Z3, X3, w, Y_0, XI_0, X_0);

while 1
  a_last = a;
  b_last = b;
  c_last = c;
  
  cnt    = cnt+1;
  
  X_p05  = X_0+a;
  Y_p05  = Y_0+b;
  Z_p05  = Z_0+c;
  
  B_p05  = interp2(X2, Y2, B, X_p05, Y_p05);
  H_p05  = interp2(X2, Y2, H, X_p05, Y_p05);
  XI_p05 = (Z_p05-B_p05)./H_p05;
  
  b_dot_p05 = interp2(X2, Y2, b_dot, X_p05, Y_p05);
  
  L_B    = XI_p05<0;
  
  if any(L_B)
    fprintf(1,'ForwardFullStep_1S2T(1): Negative elevation (relative to bed). Possibly zero ice thickness here');
    XI_p05(L_B) = 0;
    Z_p05(L_B)  = B_p05(L_B);
  end

  L_S   = XI_p05>1;
  L_abl = b_dot_p05<0;
  
  if any(L_S)
    L_done               = L_S & L_abl;    % flow trajectory emerges in ablation zone
    XI_p05(L_S & ~L_abl) = 1;              % flow trajectory emerges in accumulation zone (not legal--hence correct accuracy problem)
        
    if any(L_done)
      L_alive = ~L_done;        
      i_done  = i_alive(L_done);
      i_alive = i_alive(~L_done);
    
      L_ok(i_done) = 0;
    
      X_p(i_done)  = X_p05(L_done);
      Y_p(i_done)  = Y_p05(L_done);
      Z_p(i_done)  = B_p05(L_done) + H_p05(L_done);
      XI_p(i_done) = 1;
    
      if ~any(L_ok)
        L_ok_p = L_ok;
        return
      end
    
      a_last = a_last(L_alive);
      b_last = b_last(L_alive);
      c_last = c_last(L_alive);
    
      X_p05  = X_p05(L_alive);
      Y_p05  = Y_p05(L_alive);
      Z_p05  = Z_p05(L_alive);
      XI_p05 = XI_p05(L_alive);
    
      B_p05  = B_p05(L_alive);
      H_p05  = H_p05(L_alive);
    
      X_0    = X_0(L_alive);
      Y_0    = Y_0(L_alive);
      Z_0    = Z_0(L_alive);
      XI_0   = XI_0(L_alive);
    end                          
  end  

  a  = 0.50*dt*interp3(Y3, Z3, X3, u, Y_p05, XI_p05, X_p05);
  b  = 0.50*dt*interp3(Y3, Z3, X3, v, Y_p05, XI_p05, X_p05);
  c  = 0.50*dt*interp3(Y3, Z3, X3, w, Y_p05, XI_p05, X_p05);
  
  if any(isnan(a))
    error('ForwardFullStep_1S2T(): %d NaN values encountered in a at t=%.2f', sum(isnan(a)), t)
  end

  L_conv = abs(a-a_last)<d_max & abs(b-b_last)<d_max & abs(c-c_last)<d_max;
  if all(L_conv)
    break
  end
end

X_p1     = X_0 + 2*a;
Y_p1     = Y_0 + 2*b;
Z_p1     = Z_0 + 2*c;

B_p1     = interp2(X2, Y2, B, X_p1, Y_p1);
H_p1     = interp2(X2, Y2, H, X_p1, Y_p1);
b_dot_p1 = interp2(X2, Y2, b_dot, X_p1, Y_p1);

XI_p1     = (Z_p1-B_p1)./H_p1;

L_B        = XI_p1<0;
XI_p1(L_B) = 0;
Z_p1(L_B)  = B_p1(L_B);

L_S   = XI_p1>1;
L_abl = b_dot_p1<0;

if any(L_S)
  L_done              = L_S & L_abl;       % Flow trajectory emerges in ablation zone
  XI_p1(L_S & ~L_abl) = 1;                 % Flow trajector emerges in accumulation zone -- correct accuracy error
  
  if any(L_done)
    L_alive = ~L_done;
    i_done  = i_alive(L_done);
    i_alive = i_alive(L_alive);
  
    L_ok(i_done) = 0;
  
    X_p(i_done)  = X_p05(L_done);
    Y_p(i_done)  = Y_p05(L_done);
    Z_p(i_done)  = B_p05(L_done) + H_p05(L_done);
    XI_p(i_done) = 1;
  
    X_p1  = X_p1(L_alive);
    Y_p1  = Y_p1(L_alive);
    Z_p1  = Z_p1(L_alive);
    XI_p1 = XI_p1(L_alive);
  end
end

X_p(i_alive)  = X_p1;
Y_p(i_alive)  = Y_p1;
Z_p(i_alive)  = Z_p1;
XI_p(i_alive) = XI_p1;

L_ok_p = L_ok;
