function F_zeta = XI_to_ZETA(F_xi)

% Convert data from the XI grid to the ZETA grid using linear interpolation

% Note: It is very important to use interp1q (quick interp1)
%       because it runs VERY much faster than interp1 (at least in Apr 2014 MATLAB version)

global XI ZETA

if strcmp(F_xi.grid, 'XI')==0
  error('XI_to_ZETA(): Incorrect input grid')
end

XI_of_ZETA = ZETA;

F_zeta.A = interp1q(XI, F_xi.A, XI_of_ZETA);

F_zeta.grid = 'ZETA';

