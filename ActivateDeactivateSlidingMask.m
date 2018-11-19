function [M_t, MODEL] = ActivateDeactivateSlidingMask(p, M, X, Y, MODEL, t)

nx    = MODEL.nx;
ny    = MODEL.ny;

[nr,nc] = size(M);

if nr~=1 && nc~=1
  M     = reshape(M(p,:), ny, nx);    % Select the mask corresponding to "p" indexs
else
  M = reshape(M, ny, nx);
end  

if ~isfield(MODEL, 'THETA') || isempty(MODEL.THETA)
  
  % New origin is at upstream limit of sliding zone (M=1)
  
  X_mask_min = MODEL.X_mask_min-MODEL.x_0;
  X_mask_max = MODEL.X_mask_max-MODEL.x_0;
  Y_mask_min = MODEL.Y_mask_min-MODEL.y_0;
  Y_mask_max = MODEL.Y_mask_max-MODEL.y_0;

  MODEL.THETA = atan2(Y_mask_max-Y_mask_min, X_mask_max-X_mask_min);

  MODEL.X_rot = +(X-X_mask_min)*cos(MODEL.THETA) + (Y-Y_mask_min)*sin(MODEL.THETA);
  MODEL.Y_rot = -(X-X_mask_min)*sin(MODEL.THETA) + (Y-Y_mask_min)*cos(MODEL.THETA);

  MODEL.DL_act  = sqrt((X_mask_max-X_mask_min)^2+(Y_mask_max-Y_mask_min)^2);  % length of active zone'
  
  if isnan(MODEL.dt_PEAK(p))
    MODEL.V_act   = NaN;
    MODEL.V_deact = NaN;
  else  
    MODEL.V_act   = MODEL.DL_act/MODEL.dt_PEAK(p);                       % velocity of activation wave
    MODEL.V_deact = MODEL.DL_act/(MODEL.dt_SURGE(p)-MODEL.dt_PEAK(p));   % velocity of deactivation wave
  end
end  

if isnan(MODEL.dt_SURGE(p))
  state = '0';
else
  [t_PHASE, t_MODEL, state] = PhaseClock(t, MODEL.t_CYCLE(p), MODEL.dt_SURGE(p), MODEL.dt_PEAK(p));
end

switch state
  case '0'
    M_t        = zeros(size(M));
  case '+'
    L_front    = MODEL.V_act*t_MODEL;                              % location of surge front
    arg_X_rot  = 2.5*(MODEL.X_rot-L_front)/MODEL.dL_front(p);
    arg_X0_rot = 2.5*MODEL.X_rot/MODEL.dL_front(p);
    X_rot_act  = 0.5*(tanh(arg_X0_rot)-tanh(arg_X_rot));  
    M_t        = M.*X_rot_act;      
  case '-'
    L_front    = MODEL.DL_act - MODEL.V_deact*(t_MODEL-MODEL.dt_PEAK(p));
    arg_X_rot  = 2.5*(MODEL.X_rot-L_front)/MODEL.dL_front(p);
    arg_X0_rot = 2.5*MODEL.X_rot/MODEL.dL_front(p);    
    X_rot_act  = 0.5*(tanh(arg_X0_rot)-tanh(arg_X_rot));
    M_t        = M.*X_rot_act;
  otherwise
    error('main_activation_waves(): Unprogrammed glacier state')
end
 



