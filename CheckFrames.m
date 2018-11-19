function OK = CheckFrames(F, nx, ny)

% function checks for overlap of patch grafts

global nx_pad ny_pad

L = false(ny, nx);

OK = true;

for n=1:numel(F)
  
  if F(n).cnt>0                   % only investigate occupied frames
    i_lo = F(n).i_LO+nx_pad;
    i_hi = F(n).i_HI-nx_pad;
    j_lo = F(n).j_LO+ny_pad;
    j_hi = F(n).j_HI-ny_pad;
  
    if any(reshape(L(j_lo:j_hi, i_lo:i_hi), j_hi-j_lo+1, i_hi-i_lo+1))
      OK = false;
      break
    end
  
    L(j_lo:j_hi,i_lo:i_hi) = true;
  end  
end  