function [S_out dt_min] = SolvePatch(S, B, b_dot, F, yr_inp, yr_out)

global CFL nx_sub_default ny_sub_default cumsum_patch_cnt

% persistent sub_default
global sub_default

cumsum_patch_cnt = cumsum_patch_cnt+1;

if F.nx~=nx_sub_default || F.ny~=ny_sub_default
  [~, use] = SetupSubgridIndexArrays(F.nx, F.ny);
  sub    = use;
elseif isempty(sub_default)
  [~, sub_default] = SetupSubgridIndexArrays(nx_sub_default, ny_sub_default);
  sub = sub_default;
else
  sub = sub_default;
end

dt_yr  = yr_out-yr_inp;
t_done = yr_inp;

cnt    = 0;

dt_min = Inf;
dt_note = 1.0e-4;

while 1
  cnt = cnt+1;
  [D_IC_jc, D_IP_jc, D_ic_JC, D_ic_JP] = diffusion_subgrid(S, B, sub);
  
  try

  div_q  = -(D_IP_jc.*(S(sub.ip_jc)-S(sub.ic_jc)) - D_IC_jc.*(S(sub.ic_jc)-S(sub.im_jc)) ...
            +D_ic_JP.*(S(sub.ic_jp)-S(sub.ic_jc)) - D_ic_JC.*(S(sub.ic_jc)-S(sub.ic_jm)));
              
  catch
    fprintf(1,'DAMMIT at %.2f\n', yr_inp);  
  end
  
  fprintf(1,'SolvePatch(): Success at %.2f\n', yr_inp);
          
  dt_CFL = CFL/max([abs(D_ic_JP)' abs(D_ic_JC)' abs(D_IP_jc)' abs(D_IC_jc)']);
  
  if dt_CFL<dt_min
    dt_min = dt_CFL;
    if cnt==1 && dt_min<1.0e-4
      fprintf(1,' *** CFL time step dt_CFL=%.3e yr is very small\n', dt_CFL);
      fprintf(1,'     Frame details: j_mid = %d; i_mid=%d\n', round(0.5*(F.j_LO+F.j_HI)), ...
        round(0.5*(F.i_LO+F.i_HI)));
    end
  end
  
  dt_use = min([dt_CFL  yr_out-t_done]);
  t_done = t_done + dt_use;
  S_out  = S + (b_dot - div_q)*dt_use;
  
  S_out(S_out<=B) = B(S_out<=B);
     
  if t_done==yr_out 
    break
  end
  
end



