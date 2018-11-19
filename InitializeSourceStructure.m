function SOURCE = InitializeSourceStructure(X_source, Y_source, Name_source, t_last_surge, MODEL)

ns      = numel(X_source);

t_CYCLE = MODEL.t_CYCLE;
x_grid  = MODEL.x_grid;
y_grid  = MODEL.y_grid;

load(MODEL.file_balance, 'b_dot');

b_dot_source = interp2(x_grid, y_grid, b_dot, X_source, Y_source);
L_abl        = b_dot_source<=0;

if any(L_abl)
  ns    = numel(X_source);
  cnt   = 1:ns;
  i_abl = cnt(L_abl);
  fprintf(1,'\nInitializeSourceStructure(): %d source sites are located in the ablation zone (i.e., no englacial pathway:\n\n', numel(i_abl));
  for ii=1:numel(i_abl)
    fprintf(1,'%2d. Site %s\n', ii, Name_source{i_abl(ii)});
  end
  error('Run terminated')
end   

SOURCE(ns) = struct('X_source', [], 'Y_source', [], 'Name_source', [], 't_last_surge', [], 'X_e', [], 'Y_e', [], 'Z_e', [], 't_e', []);

for s=1:ns
  SOURCE(s).X_source     = X_source(s);
  SOURCE(s).Y_source     = Y_source(s);
  SOURCE(s).b_dot_source = b_dot_source(s);
  SOURCE(s).Name_source  = Name_source{s};
  SOURCE(s).t_last_surge = t_last_surge;
  
  SOURCE(s).B_p    = [];
  SOURCE(s).H_p    = [];
  SOURCE(s).S_p    = [];
  SOURCE(s).X_p    = [];
  SOURCE(s).Y_p    = [];
  SOURCE(s).Z_p    = [];
  SOURCE(s).XI_p   = [];
  SOURCE(s).u_p    = [];
  SOURCE(s).v_p    = [];
  SOURCE(s).w_p    = [];
end
