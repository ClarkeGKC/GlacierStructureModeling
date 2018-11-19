function  [S_out, H_out, yr_end, LAMBDA_max, T, V, L, DD, E] = ...
  step_multi_year(ModelName, S_inp, B, b_dot, M, C_t, T, T_S, q_G, dt_multi_yr, n_sub_steps, yr_start)

global VERBOSITY
global PatchStruct

PatchStruct.t_explicit = 0;

t_step = tic;

t       = yr_start;
dt      = dt_multi_yr/n_sub_steps;

PatchStruct.calls      = 0;
PatchStruct.e_cnt      = 0;
PatchStruct.t_implicit = 0;

for i=1:n_sub_steps
  [S_out, t_n, LAMBDA_max_step, T, V, L, DD, E] = step(S_inp, B, b_dot, C_t, T, T_S, q_G, dt, t);

  t  = t_n;
  S_inp = S_out;
  if i==1 || LAMBDA_max_step>LAMBDA_max
    LAMBDA_max = LAMBDA_max_step;
  end
end

H_out  = S_out-B;    % This output is not essential but it helps reduce confusion when isostatic adjustment is applied

v_max  = max(max(sqrt(V.u_ic_jc.^2 + V.v_ic_jc.^2 + V.w_ic_jc.^2)));
v_surf = sqrt(V.u_ic_jc(1, :).^2 + V.v_ic_jc(1, :).^2 + V.w_ic_jc(1,:).^2);

[dim_1,~] = size(M);

if dim_1==numel(B)
  v_max_M = max(v_surf(M==1));
else
  v_max_M = 0;
  for d=1:dim_1
    v_max_M = max(max(v_surf(M(d,:)==1)), v_max_M);
  end
end

if VERBOSITY>=1
  fprintf(1,'MODEL %s(%.0f s) => step_multi_year(%.2f): yr_start=%.2f; max(H)=%.2f m; max(b_dot)=%.4f m/yr (i.e.); max(v)=%.2f m/yr; max[v(M==1)]=%.2f m/yr\n', ...
    ModelName, toc(t_step), t_n, yr_start, max(H_out), max(b_dot), v_max, v_max_M);
end  

PatchStruct.e_rate = PatchStruct.e_cnt/n_sub_steps;

yr_end = yr_start+dt_multi_yr;
