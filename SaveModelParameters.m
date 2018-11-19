function SaveModelParameters(file_name, MODEL, par, file_source, file_date)

MODEL.N_xy     = MODEL.nx*MODEL.ny;
MODEL.N_all    = MODEL.N_xy*MODEL.n_XI;

MODEL.par      = par;

MODEL.OMEGA    = par.OMEGA;
MODEL.A_GLEN   = par.A_GLEN;
MODEL.n_GLEN   = par.n_GLEN;
MODEL.A_0      = par.A_0;
MODEL.Q_creep  = par.Q_creep;
MODEL.R_gas    = par.R_gas;
MODEL.T_KELVIN = par.T_KELVIN;
MODEL.BETA     = par.BETA;
MODEL.ETA_0    = par.ETA_0;
MODEL.YRSEC    = par.YRSEC;

FullList = fieldnames(MODEL);

OrderedFieldList = {'name', 'descriptor', 'file_dem', 'file_slide', 'file_balance', 'file_exclude', 'nx', 'ny', 'n_ZETA', 'n_XI', 'N_xy', 'N_all', 'dx', ...
  'dy', 'd_ZETA', 'd_XI', 'dt', 'n_sub_steps', 'OMEGA', 't_START', 't_STOP', 't_tracer_START', 't_archive_START', 'RunTime_hr', 't_CYCLE', 'dt_SURGE', 'dt_PEAK', ...
  'dL_front', 'X_mask_min', 'X_mask_max', 'Y_mask_min', 'Y_mask_max', 'n_GLEN', 'A_GLEN', 'A_0', 'Q_creep', 'R_gas', 'T_KELVIN', 'BETA', 'm_SLIDE', 'C_SLIDE', ...
  'C_SURGE', 'ETA_0', 'YRSEC', 'x_grid', 'y_grid', 'ZETA', 'XI', 'THERMAL', 'GRID_shift', 'PLOT', 'TRACER', 'TRACER_METHOD', 'MB_CHECK', 'VERBOSITY', ...
  'x_0', 'y_0', 'z_0', 'T_S', 'q_G', 'par'};

% Remove fields that are not on ordered field list (but may have been
% introduced en passant for internal calculations)

% These actions are local to this script and do not influence the main code

for f=1:numel(FullList)
  if ~any(strcmp(FullList{f}, OrderedFieldList))
    MODEL = rmfield(MODEL, FullList{f});
  end  
end

% Set non-defined but required fields to empty

for f=1:numel(OrderedFieldList)
  if ~isfield(MODEL, OrderedFieldList{f})
    eval(sprintf('MODEL.%s = [];', OrderedFieldList{f}))
  end
end

MODEL = orderfields(MODEL, OrderedFieldList);

save(file_name, 'MODEL', 'file_source', 'file_date')
