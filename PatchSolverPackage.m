function [S_out, F] = PatchSolverPackage(ic_jc_err, S, S_out_implicit, B, b_dot, t, dt)

global      IC_JC
global      nx_sub_default ny_sub_default nx_pad ny_pad

ic_jc_2d = IC_JC.ic_jc_2d;

nx    = IC_JC.nx;
ny    = IC_JC.ny;
N     = IC_JC.N;

S_out = S_out_implicit;    % Use the superimplicit solution except within the patches
 
F = DefineDefaultSubframes(ic_jc_err, nx_sub_default, ny_sub_default);
 
if numel(F)>1
  F  = EditSubframes(F, ic_jc_err, nx, ny, N);
  OK = CheckFrames(F, nx, ny);  % Check for overlapping frames -- this could undermine the subsequent patch grafting
  if ~OK
    error('step(): Subframe overlap of future patch grafts has been detected')
  end
   
  cnt_sum = sum([F.cnt]);

  if cnt_sum~=numel(ic_jc_err)     % Belt & suspenders test
    error('step(): sum([F.cnt])=%d but should equal numel(ic_jc_err)=%d', cnt_sum, numel(ic_jc_err))
  end  
end
 
n_frames = numel(F);

for n=1:n_frames
  if F(n).cnt>0
    [ny_2d, nx_2d] = size(ic_jc_2d);
    if F(n).j_HI>ny_2d || F(n).i_HI>nx_2d
      fprintf(1,'PatchSolverPackage(): Array size error for ic_jc_2d for frame F(%d)\n', n);
      fprintf(1,'  ny_2d = %d\n', ny_2d);
      fprintf(1,'  nx_2d = %d\n', nx_2d);
      fprintf(1,'  F(n).j_LO = %d\n', F(n).j_LO);
      fprintf(1,'  F(n).j_HI = %d\n', F(n).j_HI);
      fprintf(1,'  F(n).i_LO = %d\n', F(n).i_LO);
      fprintf(1,'  F(n).i_HI = %d\n', F(n).i_HI);
    end  
    ic_jc_sub  = reshape(ic_jc_2d(F(n).j_LO:F(n).j_HI, F(n).i_LO:F(n).i_HI), F(n).ny*F(n).nx, 1);
    S_inp_sub  = S(ic_jc_sub);       % Save the input to the implicit solver (as platform for adding)
    B_sub      = B(ic_jc_sub);
    b_dot_sub  = b_dot(ic_jc_sub);
      
    [S_out_sub, dt_min] = SolvePatch(S_inp_sub, B_sub, b_dot_sub, F(n), t, t+dt);
      
    nx_sub     = F(n).nx;
    ny_sub     = F(n).ny;
    
    % CheckFrames has established that no grafts will overlap hence they
    % can be patched in without special attention to this matter
      
    if F(n).i_LO==1
      i_LO_g = 1;
      i_LO_p = 1;
    else
      i_LO_g = F(n).i_LO+nx_pad;
      i_LO_p = 1+nx_pad;
    end
    if F(n).i_HI==nx
      i_HI_g = nx;
      i_HI_p = nx_sub;
    else
      i_HI_g = F(n).i_HI-nx_pad;
      i_HI_p = nx_sub-nx_pad;
    end
    if F(n).j_LO==1
      j_LO_g = 1;
      j_LO_p = 1;
    else
      j_LO_g = F(n).j_LO+ny_pad;
      j_LO_p = 1+ny_pad;
    end
    if F(n).j_HI==ny
      j_HI_g = ny;
      j_HI_p = ny_sub;
    else
      j_HI_g = F(n).j_HI-ny_pad;
      j_HI_p = ny_sub-ny_pad;
    end
      
    nx_g = i_HI_g-i_LO_g+1;
    ny_g = j_HI_g-j_LO_g+1;
    N_g  = nx_g*ny_g;
           
    ic_jc_graft  = reshape(ic_jc_2d(j_LO_g:j_HI_g, i_LO_g:i_HI_g), N_g, 1);
   
    S_out_sub           = reshape(S_out_sub, ny_sub, nx_sub);     
    S_out_graft         = S_out_sub(j_LO_p:j_HI_p, i_LO_p:i_HI_p);   
    S_out(ic_jc_graft)  = reshape(S_out_graft, N_g, 1);
  end   
end

