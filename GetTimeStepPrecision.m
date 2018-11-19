function IPREC = GetTimeStepPrecision(dt)

dt_str = num2str(dt);

i_dec  = strfind(dt_str, '.');
dt_str = dt_str(i_dec+1:end);

IPREC = numel(dt_str);
