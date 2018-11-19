function F_edt = EditSubframes(F, ic_jc_err, nx, ny, N)

% Detect edge-frames and trim their size. Detect multiply-occupied subframes, 
% modify their dimensions and edit the frame list to remove duplications

global nx_skirt ny_skirt nx_pad ny_pad

% Step 1: Trim dimensions of edge frames

for n=1:numel(F)
  if F(n).ic<=nx_skirt 
    F(n).i_LO = 1;
    F(n).i_HI = F(n).ic+nx_skirt;
    F(n).nx   = F(n).i_HI-F(n).i_LO+1;
  end
  if F(n).ic+nx_skirt>nx
    F(n).i_HI = nx;
    F(n).i_LO = F(n).ic-nx_skirt;
    F(n).nx   = F(n).i_HI-F(n).i_LO+1;
  end  
  if F(n).jc<=ny_skirt
    F(n).j_LO = 1;
    F(n).j_HI = F(n).jc+ny_skirt;
    F(n).ny   = F(n).j_HI-F(n).j_LO+1;
  end
  if F(n).jc+ny_skirt>ny
    F(n).j_HI = ny;
    F(n).j_LO = F(n).jc-ny_skirt;
    F(n).ny   = F(n).j_HI-F(n).j_LO+1;
  end
end

% Step 2: At this point every point should be uniquely associated with a frame 
%         but frame overlap is possible so a single point could lie within the
%         boundaries of multiple frames. Search for multiply occupied frames and
%         merge them into a superframe containing more than one point

ic_jc_loc = zeros(N,1);
ic_jc_loc(ic_jc_err) = ic_jc_err;
ic_jc_loc = reshape(ic_jc_loc, ny, nx);

EXPAND = false;

n = 1;

while 1
  ic_jc_sub = ic_jc_loc(F(n).j_LO:F(n).j_HI, F(n).i_LO:F(n).i_HI); 
  ic_jc_sub = reshape(ic_jc_sub, F(n).nx*F(n).ny, 1);
  ic_jc_pts = ic_jc_sub(ic_jc_sub~=0);
  
  if numel(ic_jc_pts)>1
    
    EXPAND = true;
    
    OK = false;
    
    while ~OK && F(n).cnt>0
      
      [jc_pts, ic_pts] = k_to_ji(ic_jc_pts, ny);
         
      i_LO = max(1,  min(ic_pts)-nx_skirt);
      i_HI = min(nx, max(ic_pts)+nx_skirt);
      j_LO = max(1,  min(jc_pts)-ny_skirt);
      j_HI = min(ny, max(jc_pts)+ny_skirt);
    
      ic_jc_new = ic_jc_loc(j_LO:j_HI, i_LO:i_HI);
      ic_jc_new = reshape(ic_jc_new, numel(ic_jc_new), 1);
      ic_jc_new_pts = ic_jc_new(ic_jc_new~=0);
    
      if numel(ic_jc_new_pts)==numel(ic_jc_pts)
        OK = true;
      else
        ic_jc_pts = ic_jc_new_pts;
      end
    end
    
    F(n).cnt  = numel(ic_jc_pts);
    
    F(n).k    = ic_jc_pts;
    F(n).ic   = ic_pts;
    F(n).jc   = jc_pts;
    
    F(n).i_LO = i_LO;
    F(n).i_HI = i_HI;
    F(n).j_LO = j_LO;
    F(n).j_HI = j_HI;
      
    F(n).nx   = F(n).i_HI-F(n).i_LO+1;
    F(n).ny   = F(n).j_HI-F(n).j_LO+1;
      
    % Set point count to zero in frames that have been absorbed
    
    n_done = 0;
       
    for r=[1:n-1 n+1:numel(F)]
      for l=1:numel(F(r).k)
        if any(F(n).k==F(r).k(l))
          
          F(r).cnt  = F(r).cnt-1;
          n_done    = n_done+1;
        end
        if F(r).cnt==0
          break
        end
      end
       
      if F(r).cnt==0        
        F(r) = struct('k', [], 'ic', [],  'jc', [], 'i_LO', [], 'i_HI', [], ...
                      'j_LO', [], 'j_HI', [], 'nx', [], 'ny', [], 'cnt', 0);
                    
        if n_done==numel(F(n).k)-1     % Original frame-defining point remains
          break
        end
      end              
    end
    
  end
  if n==numel(F)
    break;
  end
  n = n+1;
end

% Step 3: Search for any cells which overlap in their unpadded "core" and
%         combine to form a single super cell. Overlaps should be unusual
%         unless the unpadded core of skirt is large. (If no cells were
%         expanded then this step is skipped.)

if EXPAND
  DONE = false;

  while ~DONE
    DONE = true;
  
    for n=1:numel(F)
      if F(n).cnt~=0
        for l=n+1:numel(F)
          if F(l).cnt~=0
            
            i_LO_A = F(n).i_LO+nx_pad;    % Not tidy but at least clear
            i_HI_A = F(n).i_HI-nx_pad;
            j_LO_A = F(n).j_LO+ny_pad;
            j_HI_A = F(n).j_HI-ny_pad;
          
            i_LO_B = F(l).i_LO+nx_pad;
            i_HI_B = F(l).i_HI-nx_pad;
            j_LO_B = F(l).j_LO+ny_pad;
            j_HI_B = F(l).j_HI-ny_pad;
          
            % Establish conditions for vertex of rectangle A in B or vertex of rectangle B in A

            L_nw   = (i_LO_B>=i_LO_A  &  i_LO_B< i_HI_A  &  j_LO_B>=j_LO_A  &  j_LO_B< j_HI_A) | ...
                     (i_LO_A>=i_LO_B  &  i_LO_A< i_HI_B  &  j_LO_A>=j_LO_B  &  j_LO_A< j_HI_B);
                 
            L_ne   = (i_HI_B> i_LO_A  &  i_HI_B<=i_HI_A  &  j_LO_B< j_HI_A  &  j_LO_B>=j_LO_A) | ...
                     (i_HI_A> i_LO_B  &  i_HI_A<=i_HI_B  &  j_LO_A< j_HI_B  &  j_LO_A>=j_LO_B);
                  
            L_se   = (i_HI_B> i_LO_A  &  i_HI_B<=i_HI_A  &  j_HI_B<=j_HI_A  &  j_HI_B> j_LO_A) | ...
                     (i_HI_A> i_LO_B  &  i_HI_A<=i_HI_B  &  j_HI_A<=j_HI_B  &  j_HI_A> j_LO_B);
          
            L_sw   = (i_LO_B< i_HI_A  &  i_LO_B>=i_LO_A  &  j_HI_B> j_LO_A  &  j_HI_B<=j_HI_A) | ...
                     (i_LO_A< i_HI_B  &  i_LO_A>=i_LO_B  &  j_HI_A> j_LO_B  &  j_HI_A<=j_HI_B);
          
            if L_nw || L_ne || L_se || L_sw
              F(n).cnt  = F(n).cnt + F(l).cnt;
          
              F(n).i_LO = min(F(n).i_LO, F(l).i_LO);     % Reset frame limits to new padded extents
              F(n).i_HI = max(F(n).i_HI, F(l).i_HI);
              F(n).j_LO = min(F(n).j_LO, F(l).j_LO);
              F(n).j_HI = max(F(n).j_HI, F(l).j_HI);
            
              F(n).k    = [];
              F(n).ic   = [];
              F(n).jc   = [];
              F(n).nx   = F(n).i_HI-F(n).i_LO+1;
              F(n).ny   = F(n).j_HI-F(n).j_LO+1;
            
              F(l) = struct('k', [], 'ic', [],  'jc', [], 'i_LO', [], 'i_HI', [], ...
                            'j_LO', [], 'j_HI', [], 'nx', [], 'ny', [], 'cnt', 0);

              DONE = false;
              break
            end  
          end
        end
      end
      if ~DONE
        break
      end
    end
  end
end

% Step 4: Perform consistency check

cnt_sum = 0;
for k=1:numel(F)
  cnt_sum = F(k).cnt+cnt_sum;
end

if cnt_sum~=numel(ic_jc_err)
  fprintf(1,'EditSubframes(*): cnt_sum = %d; numel(ic_jc_err)=%d; numel(F)=%d\n', cnt_sum, numel(ic_jc_err), numel(F));
end  
F_edt = F;