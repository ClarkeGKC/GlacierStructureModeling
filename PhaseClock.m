function [t_PHASE, t_SURGE, state] = PhaseClock(yr, t_CYCLE, dt_SURGE, dt_PEAK)

% Determine timing into surge phase

t_PHASE = mod(yr, t_CYCLE);

if t_PHASE>=0 && t_PHASE<=dt_PEAK
  state   = '+';
  t_SURGE = t_PHASE;
elseif t_PHASE>=0 && t_PHASE<=dt_SURGE
  state   = '-';
  t_SURGE = t_PHASE;
else
  state   = '0';
  t_SURGE = 0;
end  