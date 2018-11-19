function   [x_p, y_p, z_p, t_p] = ApplyBlanking(BLANKING_MODEL, x_p, y_p, z_p, t_p, s_xy_p)

switch BLANKING_MODEL.Rule
  case 1  % Point blanking
    BlankArc = pi*BLANKING_MODEL.BlankArcDeg/180;
    cos_max  = cos(BlankArc);

    i    = 0;

    x_n  = -x_p;
    y_n  = -y_p;
    z_n  = -z_p;

    N_i  = numel(x_p);

    L_blank = false(N_i, 1);

    k         = 1;   % keep

    for i=2:numel(x_p)
  
      % Note also that we make no allowance for proximity of points that have
      % been separated by hemispheric split
   
      cos_arc_p   = x_p(k)*x_p(i) + y_p(k)*y_p(i) + z_p(k)*z_p(i);
      cos_arc_n   = x_p(k)*x_n(i) + y_p(k)*y_n(i) + z_p(k)*z_n(i);
  
      if cos_arc_p>cos_max || cos_arc_n>cos_max
        L_blank(i) = 1;
      else
        k = i;
      end
    end
  case 2  % Space blanking
    s_xy        = s_xy_p;
    Delta_s_min = BLANKING_MODEL.Distance;
    N_i         = numel(x_p);
    L_keep      = [true false(1, N_i-1)];
    cnt         = 1:N_i;
    i_keep      = cnt(L_keep);
    while 1
      s_keep_max = s_xy(max(i_keep));
      L_drop     = s_xy>s_keep_max & s_xy<s_keep_max+Delta_s_min;
      L_res      = ~L_drop & cnt>i_keep;
      if ~any(L_res)
        break
      end  
      i_keep     = min(cnt(L_res));
      L_keep(i_keep) = 1;
    end 
    L_blank = ~L_keep;
  case 3  % TIme blanking
    Delta_t_min = BLANKING_MODEL.Time; 
    N_i         = numel(x_p);
    L_keep      = [true false(1, N_i-1)];
    cnt         = 1:N_i;
    i_keep      = cnt(L_keep);
    while 1
      t_keep_max = t_p(max(i_keep));
      L_drop     = t_p>t_keep_max & t_p<t_keep_max + Delta_t_min;
      L_res      = ~L_drop & cnt>i_keep;
      if ~any(L_res)
        break
      end
      i_keep     = min(cnt(L_res));
      L_keep(i_keep) = 1;
    end
    L_blank = ~L_keep;
  case 4  % Space OR time blanking
    s_xy        = s_xy_p;
    Delta_s_min = BLANKING_MODEL.Distance;
    Delta_t_min = BLANKING_MODEL.Time;
    N_i         = numel(x_p);
    L_keep      = [true  false(1, N_i-1)];
    cnt         = 1:N_i;
    i_keep      = cnt(L_keep);
    while 1
      s_keep_max = s_xy(max(i_keep));
      t_keep_max = t_p(max(i_keep));
      L_drop_s   = s_xy>s_keep_max & s_xy<s_keep_max+Delta_s_min;
      L_drop_t   = t_p>t_keep_max & t_p<t_keep_max + Delta_t_min;
      L_drop     = L_drop_s & L_drop_t;   % Most realistic condition: EITHER the transit distance or healing time condition must be satisified
      
      try
        L_res      = ~L_drop & cnt>i_keep;
      catch
        disp('OOPS')
      end  
      if ~any(L_res)
        break
      end
      i_keep = min(cnt(L_res));
      L_keep(i_keep) = 1;
    end
    L_blank = ~L_keep;
  case 5  % Space AND time blanking
    s_xy        = s_xy_p;
    Delta_s_min = BLANKING_MODEL.Distance;
    Delta_t_min = BLANKING_MODEL.Time; 
    N_i         = numel(x_p);
    L_keep      = [true  false(1, N_i-1)];
    cnt         = 1:N_i;
    i_keep      = cnt(L_keep);
    while 1
      s_keep_max = s_xy(max(i_keep));
      t_keep_max = t_p(max(i_keep));
      L_drop_s   = s_xy>s_keep_max & s_xy<s_keep_max+Delta_s_min;
      L_drop_t   = t_p>t_keep_max & t_p<t_keep_max + Delta_t_min;
      L_drop     = L_drop_s | L_drop_t;    % Most stringent condition: BOTH both transit distance and healing time condition must be satisifed
      L_res      = ~L_drop & cnt>i_keep;
      if ~any(L_res)
        break
      end
      i_keep = min(cnt(L_res));
      L_keep(i_keep) = 1;
    end
    L_blank = ~L_keep;
    
  otherwise
    error('ApplyBlanking(): Unprogrammed case')
end

x_p = x_p(~L_blank);
y_p = y_p(~L_blank);
z_p = z_p(~L_blank);

t_p = t_p(~L_blank);

