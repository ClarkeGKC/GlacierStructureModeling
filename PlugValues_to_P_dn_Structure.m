function P_dn = PlugValues_to_P_dn_Structure(P_dn, t_pt, X_pt, Y_pt, Z_pt, XI_pt,  B_pt, H_pt, S_pt, b_dot_pt, M_pt, ...
  u_pt, v_pt, w_pt, L_pt, s_pt, A_pt, ETA_pt, TAU_pt, T_pt)

s_seq  = 1:numel(X_pt);
s_list = s_seq(~isnan(X_pt));

for s=s_list
  P_dn(s).t  = [P_dn(s).t t_pt];
  P_dn(s).X  = [P_dn(s).X X_pt(s)];
  P_dn(s).Y  = [P_dn(s).Y Y_pt(s)];
  P_dn(s).Z  = [P_dn(s).Z Z_pt(s)];
  P_dn(s).XI = [P_dn(s).XI XI_pt(s)];
  P_dn(s).B  = [P_dn(s).B B_pt(s)];
  P_dn(s).H  = [P_dn(s).H H_pt(s)];
  P_dn(s).S  = [P_dn(s).S S_pt(s)];
  P_dn(s).b_dot = [P_dn(s).b_dot b_dot_pt(s)];
  P_dn(s).M  = [P_dn(s).M M_pt(s)];
  P_dn(s).u  = [P_dn(s).u u_pt(s)];
  P_dn(s).v  = [P_dn(s).v v_pt(s)];
  P_dn(s).w  = [P_dn(s).w w_pt(s)];
  P_dn(s).L  = [P_dn(s).L L_pt(:,:,s)];
  P_dn(s).s  = [P_dn(s).s s_pt(:,:,s)];
  P_dn(s).A  = [P_dn(s).A A_pt(s)];
  P_dn(s).ETA = [P_dn(s).ETA ETA_pt(s)];
  P_dn(s).TAU = [P_dn(s).TAU TAU_pt(s)];
  P_dn(s).T   = [P_dn(s).T T_pt(s)];
end  






