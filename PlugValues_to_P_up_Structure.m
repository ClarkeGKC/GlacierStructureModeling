function P_up = PlugValues_to_P_up_Structure(P_up, t_pt, X_pt, Y_pt, Z_pt, XI_pt, S_pt)

s_seq  = 1:numel(X_pt);
s_list = s_seq(~isnan(X_pt));

for s=s_list
  P_up(s).t  = [P_up(s).t  t_pt];
  P_up(s).X  = [P_up(s).X  X_pt(s)];
  P_up(s).Y  = [P_up(s).Y  Y_pt(s)];
  P_up(s).Z  = [P_up(s).Z  Z_pt(s)];
  P_up(s).XI = [P_up(s).XI  XI_pt(s)];
  P_up(s).S  = [P_up(s).S  S_pt(s)];
end  






