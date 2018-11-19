function frame = DefineDefaultSubframes(ic_jc_err, nx_sub_default, ny_sub_default)

% function scans through the problem points and defines locations and size
% of subframes as well as the location of the frame centre point

% when two or more points are located within the same mini subframe the
% frame size is redefined and recentred. Thus the number of subframes <=
% number of probloem points

global IC_JC nx_skirt ny_skirt

nx = IC_JC.nx;
ny = IC_JC.ny;

L = false(ny, nx);

L(ic_jc_err) = true;

n_err = numel(ic_jc_err);

frame(n_err) = struct('k', [], 'ic', [],  'jc', [], 'i_LO', [], 'i_HI', [], ...
  'j_LO', [], 'j_HI', [], 'nx', [], 'ny', [], 'cnt', []);

for n=1:numel(ic_jc_err);
  L(ic_jc_err(n)) = false;
  
  [jc, ic]  = k_to_ji(ic_jc_err(n), ny);
  
  frame(n).k   = ic_jc_err(n);
  frame(n).ic  = ic;
  frame(n).jc  = jc;
  frame(n).cnt = 1;

  % Set the limits of the rectangular patch to be extracted from the maxi grid
  % If error source point is near boundaries of the maxi grid these limits
  % must take this into account. The ideal graming places the error source
  % at the centre of the mini grid but his is not possible if the patch is
  % near boundaries of the maxi grid in which case we use a non-centred 
  % frame that touches one or two edges of the maxi-grid

  if ic-nx_skirt-1<0 
    frame(n).i_LO  = 1;
    frame(n).i_HI  = frame(n).i_LO+nx_sub_default-1;
  elseif ic+nx_skirt-nx>0
    frame(n).i_HI  = nx;
    frame(n).i_LO  = frame(n).i_HI-nx_sub_default+1;
  else  
    frame(n).i_LO  = max(1,  ic-nx_skirt);
    frame(n).i_HI  = min(nx, ic+nx_skirt);
  end

  if jc-ny_skirt-1<0
    frame(n).j_LO   = 1;
    frame(n).j_HI   = frame(n).j_LO+ny_sub_default-1;
  elseif jc+ny_skirt-ny>0
    frame(n).j_HI  = ny;
    frame(n).j_LO  = frame(n).j_HI-ny_sub_default+1; 
  else  
    frame(n).j_LO  = max(1,  jc-ny_skirt);
    frame(n).j_HI  = min(ny, jc+ny_skirt);
   end
  
   frame(n).nx = frame(n).i_HI-frame(n).i_LO+1;
   frame(n).ny = frame(n).j_HI-frame(n).j_LO+1;
end

