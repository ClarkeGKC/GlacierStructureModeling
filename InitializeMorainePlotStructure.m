function PLOT = InitializeMorainePlotStructure(yr_start, SOURCE)

ns = numel(SOURCE);
PLOT(ns) = struct('yr_start', [], 'yr_display', [], 'X_source', [], 'Y_source', [], 'X_e', [], 'Y_e', [], 't_e', [], 't_src', [], ...
  'X_trk', [], 'Y_trk', [], 't_trk', [], 't_e_trk', [], 'H_trk', []);

for s=1:ns
  PLOT(s).X_source   = SOURCE(s).X_source;
  PLOT(s).Y_source   = SOURCE(s).Y_source;
  PLOT(s).X_e        = SOURCE(s).X_e;
  PLOT(s).Y_e        = SOURCE(s).Y_e;
  PLOT(s).t_e        = SOURCE(s).t_e;
end

  
