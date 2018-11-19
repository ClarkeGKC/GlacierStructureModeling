function [X_m, Y_m, Z_m, XI_m, S_m, L_ok_m] = BackwardFullStep_1S2T(X, Y, Z, B, H, b_dot, u, v, w, t, L_ok)

% Note "Backward" half steps are steps upflow whereas foreward halfsteps are downflow
% This is relevant to treatment of points that emerge prematurely  and are
% diagnosed basis of sign of b_dot_p

% Note that inputs X, Y, Z, L_ok and outputs X_m, Y_m, Z_m, XI_m, S_m, L_ok_m never change size
% The number of true elements in L_ok decreases with time until all work is done

% Inside this function X_m, Y_m, XI_m etc are redefined to include only the
% currently active particles that are being traced

global dt X2 Y2 X3 Y3 Z3

d_max     = 0.0001;   % maximum distance error for convergence
ITMAX     = 25;       % Maximum iteration count before point treated as STALLED

X_m       = NaN(size(L_ok));
Y_m       = NaN(size(L_ok));
Z_m       = NaN(size(L_ok));
XI_m      = NaN(size(L_ok));
S_m       = NaN(size(L_ok)); 

X_0       = X(L_ok);
Y_0       = Y(L_ok);
Z_0       = Z(L_ok);

B_0       = interp2(X2, Y2, B, X_0, Y_0);
H_0       = interp2(X2, Y2, H, X_0, Y_0);
XI_0      = min((Z_0-B_0)./H_0, 1);

S_0       = B_0 + H_0;
Z_0       = min(Z_0, S_0);

counter   = 1:numel(L_ok);

i_alive   = counter(L_ok);

% This was added to prevent NaN velocity values from being encountered
% while exploring potential step endpoints

u(isnan(u)) = 0;
v(isnan(v)) = 0;
w(isnan(w)) = 0;

cnt     = 0;

a       = 0.50*dt*interp3(Y3, Z3, X3, u, Y_0, XI_0, X_0);
b       = 0.50*dt*interp3(Y3, Z3, X3, v, Y_0, XI_0, X_0);
c       = 0.50*dt*interp3(Y3, Z3, X3, w, Y_0, XI_0, X_0);

while 1
  a_last = a;
  b_last = b;
  c_last = c;
  
  cnt    = cnt+1;
 
  X_m05  = X_0-a;
  Y_m05  = Y_0-b;
  Z_m05  = Z_0-c;
  
  B_m05     = interp2(X2, Y2, B, X_m05, Y_m05);
  H_m05     = interp2(X2, Y2, H, X_m05, Y_m05);  
  S_m05     = B_m05 + H_m05;
  XI_m05    = (Z_m05-B_m05)./H_m05;
  
  if any(isnan(XI_m05))
    fprintf(1,'BackwardFullStep_1S2T(1): NaN value(s) of XI_m05 (probably because H_m05=0) set to XI_m05=0.5 at t=%.3f\n', t);
    L_nan = isnan(XI_m05);
    XI_m05(L_nan) = 0.5;
  end
      
  L_B         = XI_m05<0;
  
  XI_m05(L_B) = 0;
  Z_m05(L_B)  = B_m05(L_B);
   
  L_S         = XI_m05>1;
  L_H         = H_m05<=0;

  b_dot_m05   = interp2(X2, Y2, b_dot, X_m05, Y_m05);
  L_acc       = b_dot_m05>0;
  
  L_test      = (L_S & L_acc) | L_H;    
  
  if any(L_test)
    L_done       = L_test;
    L_alive      = ~L_done;
    i_done       = i_alive(L_done);
    i_alive      = i_alive(~L_done);
       
    L_ok(i_done) = 0;
    
    X_m(i_done)  = X_m05(L_done);
    Y_m(i_done)  = Y_m05(L_done);
    XI_m(i_done) = 1;
    S_m(i_done)  = B_m05(L_done) + H_m05(L_done);
    Z_m(i_done)  = S_m(i_done);
    
    if ~any(L_ok)
      L_ok_m = L_ok;
      return
    end
    
    a_last       = a_last(L_alive);
    b_last       = b_last(L_alive);
    c_last       = c_last(L_alive);
   
    X_m05        = X_m05(L_alive);
    Y_m05        = Y_m05(L_alive);
%   Z_m05        = Z_m05(L_alive);
    S_m05        = S_m05(L_alive);
    XI_m05       = XI_m05(L_alive);
    
    X_0          = X_0(L_alive);
    Y_0          = Y_0(L_alive);
    Z_0          = Z_0(L_alive);
    XI_0         = XI_0(L_alive);
  end
  
  L_tweak = L_S & ~L_acc;
  XI_m05(L_tweak) = 1;
  
  a  = 0.50*dt*interp3(Y3, Z3, X3, u, Y_m05, XI_m05, X_m05);
  b  = 0.50*dt*interp3(Y3, Z3, X3, v, Y_m05, XI_m05, X_m05);
  c  = 0.50*dt*interp3(Y3, Z3, X3, w, Y_m05, XI_m05, X_m05);
    
  if any(isnan(a))
    fprintf(1,'\nBackwardFullStep_1S2T(): %d NaN values encountered in a at t=%.2f\n', sum(isnan(a)), t);
    fprintf(1,'\nOffending points are as follows:\n');
    fprintf(1,'  X_m05    Y_m05   XI_m05   H_m05\n\n');
    cnt   = 1:numel(X_m05);
    i_nan      = cnt(isnan(a));
    X_m05_nan  = X_m05(i_nan);
    Y_m05_nan  = Y_m05(i_nan);
    XI_m05_nan = XI_m05(i_nan);
    H_m05_nan  = interp2(X2, Y2, H, X_m05_nan, Y_m05_nan);
    for ii = 1:numel(i_nan)
      fprintf(1,' %9.3f  %9.3f %6.3f %8.3f\n', X_m05_nan(ii), Y_m05_nan(ii), XI_m05_nan(ii), H_m05_nan(ii));
    end
    error('BackwardFullStep_1S2T(): %d NaN values encountered in a at t=%.2f', sum(isnan(a)), t)
  end

  L_conv = abs(a-a_last)<d_max & abs(b-b_last)<d_max & abs(c-c_last)<d_max;
  if all(L_conv)
    break
  end
  
  if cnt>=ITMAX
    L_alive   = L_conv;
    L_stalled = ~L_conv;

    fprintf(1,'\nBackwardFullStep_1S2T(): Maximum iteration %d reached at time=%.2f for %d point(s). Integration of these trajectories is terminated\n\n', cnt, t, sum(L_stalled));
    
    i_stalled    = i_alive(L_stalled);
    i_alive      = i_alive(L_alive);
    
    L_ok(i_stalled) = 0;
    
    X_0_stalled   = X_0(L_stalled);
    Y_0_stalled   = Y_0(L_stalled);
    XI_0_stalled  = XI_0(L_stalled);
    b_dot_stalled = b_dot_m05(L_stalled);
    H_stalled     = H_m05(L_stalled);
    
    for p=1:numel(i_stalled)
      fprintf('\n  Point #%d:  At X_p=%.2f; Y_p=%.2f; XI_p=%.2f; b_dot_p=%.3f; H_p=%.3f\n', i_stalled(p), X_0_stalled(p), Y_0_stalled(p), XI_0_stalled(p), b_dot_stalled(p), H_stalled(p));
    end
    fprintf(1,'\n');   
   
    X_0          = X_0(L_alive);
    Y_0          = Y_0(L_alive);
    Z_0          = Z_0(L_alive);
    XI_0         = XI_0(L_alive);    
    
    a  = 0.50*dt*interp3(Y3, Z3, X3, u, Y_0, XI_0, X_0);
    b  = 0.50*dt*interp3(Y3, Z3, X3, v, Y_0, XI_0, X_0);
    c  = 0.50*dt*interp3(Y3, Z3, X3, w, Y_0, XI_0, X_0);
    
    if any(isnan(a))
      error('BackwardFullStep_1S2T(): %d NaN values encountered in a at t=%.2f', sum(isnan(a)), t)
    end
    
    break
  end  
end

X_m1 = X_0 - 2*a;
Y_m1 = Y_0 - 2*b;
Z_m1 = Z_0 - 2*c;

B_m1 = interp2(X2, Y2, B, X_m1, Y_m1);
H_m1 = interp2(X2, Y2, H, X_m1, Y_m1);
S_m1 = B_m1 + H_m1;

XI_m1 = (Z_m1-B_m1)./H_m1;

if any(isnan(XI_m1))
  fprintf(1,'BackwardFullStep_1S2T(2): NaN value(s) of XI_m1 (probably because H_m15=0) set to XI_m1=0.5 at t=%.3f\n', t);
  L_nan        = isnan(XI_m1);
  XI_m1(L_nan) = 0.5;
end

L_B  = XI_m1<0;
XI_m1(L_B) = 0;
Z_m1(L_B)  = B_m1(L_B);

L_S = XI_m1>1; 
L_H = H_m1<=0;

b_dot_m1    = interp2(X2, Y2, b_dot, X_m1, Y_m1);
L_acc       = b_dot_m1>0;

L_test      = (L_S & L_acc) | L_H;

if any(L_test)
  L_done       = L_test;
  L_alive      = ~L_done;
  i_done       = i_alive(L_done);
  i_alive      = i_alive(L_alive);
  
  L_ok(i_done) = 0;
    
  X_m(i_done)  = X_m05(L_done);
  Y_m(i_done)  = Y_m05(L_done);
  S_m(i_done)  = S_m05(L_done);
  Z_m(i_done)  = S_m05(L_done);
  XI_m(i_done) = 1;
  
  X_m1  = X_m1(L_alive);
  Y_m1  = Y_m1(L_alive);
  Z_m1  = Z_m1(L_alive);
  XI_m1 = XI_m1(L_alive);
  S_m1  = S_m1(L_alive);
end

% Re-establish the original sizes to match input (X,Y,Z) but do not change names

X_m(i_alive)  = X_m1;
Y_m(i_alive)  = Y_m1;
Z_m(i_alive)  = Z_m1;
XI_m(i_alive) = XI_m1;
S_m(i_alive)  = S_m1;

L_ok_m        = L_ok;

