function v = geographical_polar(theta, rho, OPTIONS)
    %POLAR  Polar coordinate plot.
    %   GKCC modification to convert MATLAB polar.m script to geographical
    %   orientation (North=0 with clockwise rotation positive)
    %   POLAR(THETA, RHO) makes a plot using polar coordinates of
    %   the angle THETA, in radians, versus the radius RHO.
    %   POLAR(THETA, RHO, S) uses the linestyle specified in string S.
    %   See PLOT for a description of legal linestyles.
    %
    %   POLAR(AX, ...) plots into AX instead of GCA.
    %
    %   H = POLAR(...) returns a handle to the plotted object in H.
    %
    %   Example:
    %      t = 0 : .01 : 2 * pi;
    %      polar(t, sin(2 * t) .* cos(2 * t), '--r');
    %
    %   See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.
    
    %   Copyright 1984-2015 MathWorks, Inc.
    
    if isfield(OPTIONS, 'AzTicks')
      AzTicks = OPTIONS.AzTicks;
    else
      AzTicks = false;
    end  
    
    if isfield(OPTIONS, 'AzTickInterval')
      AzTickInterval = OPTIONS.AzTickInterval;
    else
      AzTickInterval = 30;
    end  
    
    if isfield(OPTIONS, 'AzGrid')
      AzGrid = OPTIONS.AzGrid;
    else
      AzGrid = true;
    end  

    n_TickHalf = (360/AzTickInterval)/2;
    
    if isfield(OPTIONS, 'AzLabels')
      AzLabels = OPTIONS.AzLabels;
    else
      AzLabels = 'numeric';
    end
    
    if isfield(OPTIONS, 'AzLabelFontScale')
      AzLabelFontScale = OPTIONS.AzLabelFontScale;
    else
      AzLabelFontScale = 1;
    end 
    
    if isfield(OPTIONS, 'AzLabelRfract')
      AzLabelRfract = OPTIONS.AzLabelRfract;
    else
      AzLabelRfract = 0.10;
    end
    
    if isfield(OPTIONS, 'AzTickRfract')
      AzTickRfract = OPTIONS.AzTickRfract;
    else
      AzTickRfract = 0;
    end
    
    if isfield(OPTIONS, 'AzLabels')
      AzLabels = OPTIONS.AzLabels;
    else
      AzLabels = 'numeric';
    end  
    
    RLabel = 1 + AzLabelRfract;
    RTick  = 1 + AzTickRfract;

    DATA = true;

    if isempty(rho)
       DATA = false;
    end  
    
    cax = newplot;
    axis off
    
    next = lower(get(cax, 'NextPlot'));
    hold_state = ishold(cax);
    
    % get x-axis text color so grid is in same color
    % get the axis gridColor
    
    axColor     = get(cax, 'Color');
    gridAlpha   = get(cax, 'GridAlpha');
    axGridColor = get(cax,'GridColor').*gridAlpha + axColor.*(1-gridAlpha);
    tc          = axGridColor;
    if any(tc~=0)
      tc = [0 0 0];    % GKCC override to obtain black circle outline for Rosette plots
    end
    
    ls          = get(cax, 'GridLineStyle');
    
    % Hold on to current Text defaults, reset them to the
    % Axes' font attributes so tick marks use them.
    
    fAngle = get(cax, 'DefaultTextFontAngle');
    fName = get(cax, 'DefaultTextFontName');
    fSize = get(cax, 'DefaultTextFontSize');
    fWeight = get(cax, 'DefaultTextFontWeight');
    fUnits = get(cax, 'DefaultTextUnits');
    
    set(cax, ...
        'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
        'DefaultTextFontName', get(cax, 'FontName'), ...
        'DefaultTextFontSize', AzLabelFontScale*get(cax, 'FontSize'), ...
        'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
        'DefaultTextUnits', 'data');
    
    % only do grids if hold is off
    if ~hold_state
        
      if DATA
        % make a radial grid
        hold(cax, 'on');
        % ensure that Inf values don't enter into the limit calculation.
        arho = abs(rho(:));
        maxrho = max(arho(arho ~= Inf));
        if AzGrid
          linestyle = '-';
        else  
          linestyle = 'none';
        end 
        hhh = line([-maxrho, -maxrho, maxrho, maxrho], [-maxrho, maxrho, maxrho, -maxrho], 'linestyle', linestyle, 'Parent', cax);
        
        XLim = get(cax, 'XLim');
        YLim = get(cax, 'YLim');
        
        if -XLim(1) ~=XLim(2) || -YLim(1)~=YLim(2) || XLim(2)~=YLim(2)
          error('geographical_polar(): Expecting plot limits to be square')
        end  
        
        if XLim(2)~=maxrho
          % fprintf(1,'... geographical_polar(m): resizing XLim and/or YLim\n');
          XLim(1) = -maxrho;
          XLim(2) = maxrho;
          YLim(1) = -maxrho;
          YLim(2) = maxrho;
          set(cax, 'XLim', XLim, 'YLim', YLim)
        end
        
        set(cax, 'DataAspectRatio', [1, 1, 1], 'PlotBoxAspectRatioMode', 'auto');
        v = [get(cax, 'XLim') get(cax, 'YLim')];
        ticks = sum(get(cax, 'YTick') >= 0);
        delete(hhh);
        % check radial limits and ticks
        rmin = 0;
        rmax = v(4);
        rticks = max(ticks - 1, 2);
        if rticks > 5   % see if we can reduce the number
          if rem(rticks, 2) == 0
            rticks = rticks / 2;
          elseif rem(rticks, 3) == 0
            rticks = rticks / 3;
          end
        end
      else
        v = [get(cax, 'XLim') get(cax, 'YLim')];
        rmin = 0;
        rmax = v(4);
      end
        
      % define a circle
      th = 0 : pi / 50 : 2 * pi;
      xunit = cos(th);
      yunit = sin(th);
      % now really force points on x/y axes to lie on them exactly
      inds = 1 : (length(th) - 1) / 4 : length(th);
      xunit(inds(2 : 2 : 4)) = zeros(2, 1);
      yunit(inds(1 : 2 : 5)) = zeros(3, 1);
      % plot background if necessary
      if ~ischar(get(cax, 'Color'))
        patch('XData', xunit * rmax, 'YData', yunit * rmax, ...
              'EdgeColor', tc, 'FaceColor', get(cax, 'Color'), ...
              'HandleVisibility', 'off', 'Parent', cax);
      end
        
      % draw radial circles
        
      if DATA && AzGrid
        c82 = cos(82 * pi / 180);
        s82 = sin(82 * pi / 180);
        rinc = (rmax - rmin) / rticks;
        for i = (rmin + rinc) : rinc : rmax
          hhh = line(xunit * i, yunit * i, 'LineStyle', ls, 'Color', tc, 'LineWidth', 1, ...
              'HandleVisibility', 'off', 'Parent', cax);
        end
        set(hhh, 'LineStyle', '-'); % Make outer circle solid
        
        % plot spokes
        th = (1 : n_TickHalf) * 2 * pi / (2*n_TickHalf);
        
        th  = -th + pi/2;     % Here is where the theta plot annotation is converted to geographical azimuth
        cst = cos(th);
        snt = sin(th);
        cs = [-cst; cst];
        sn = [-snt; snt];
        R_spoke = rmax*RTick;
        line(rmax * cs, rmax * sn, 'LineStyle', ls, 'Color', tc, 'LineWidth', 1, ...
            'HandleVisibility', 'off', 'Parent', cax);
      else
        th  = (1 : n_TickHalf) * 2 * pi /(2*n_TickHalf);
        th  = -th + pi/2;
        cst = cos(th);
        snt = sin(th);    
      end
                
      % annotate spokes in degrees
      rt = RLabel * rmax;
      
      if strcmp(AzLabels, 'numeric') || n_TickHalf~=2        
        for i = 1 : length(th)
          text(rt * cst(i), rt * snt(i), int2str(i * AzTickInterval),...
            'HorizontalAlignment', 'center', ...
            'HandleVisibility', 'off', 'Parent', cax);
          if i == length(th)
            loc = int2str(0);
          else
            loc = int2str(180 + i * AzTickInterval);
          end
          text(-rt * cst(i), -rt * snt(i), loc, 'HorizontalAlignment', 'center', ...
              'HandleVisibility', 'off', 'Parent', cax);
        end
      else
        cardinal = {'E', 'S', 'W', 'N'};
        for i = 1 : length(th)
          text(rt * cst(i), rt * snt(i), cardinal{i},...
            'HorizontalAlignment', 'center', ...
            'HandleVisibility', 'off', 'FontWeight', 'bold', 'Parent', cax);
          text(-rt * cst(i), -rt * snt(i), cardinal{2+i}, 'HorizontalAlignment', 'center', ...
              'HandleVisibility', 'off', 'FontWeight', 'bold', 'Parent', cax);
        end
      end   
      
      % Plot external aximuthal ticks
        
      if AzTicks  
        th = (1 : 2*n_TickHalf) * 2 * pi / (2*n_TickHalf);        
        th  = -th + pi/2;     % Here is where the theta plot annotation is converted to geographical azimuth
        cst = cos(th);
        snt = sin(th);
        cs = [-cst; cst];
        sn = [-snt; snt];
        R_spoke = rmax*RTick;
          
        R_cs = [rmax*cst; R_spoke*cst];
        R_sn = [rmax*snt; R_spoke*snt];
          
        line(R_cs, R_sn, 'LineStyle', ls, 'Color', 'k', 'LineWidth', 1, ...
            'HandleVisibility', 'off', 'Parent', cax);
      end
        
      % set view to 2-D
      view(cax, 2);
      % set axis limits
      %axis(cax, rmax * [-1, 1, -1.15, 1.15]);
      if exist('R_spoke', 'var')
        R_max = R_spoke;
      else  
        R_max = rmax;
      end  
      axis(cax, R_max*[-1.15, 1.15, -1.15, 1.15]);
    end
    
    % Reset defaults.
    set(cax, ...
      'DefaultTextFontAngle', fAngle , ...
      'DefaultTextFontName', fName , ...
      'DefaultTextFontSize', fSize, ...
      'DefaultTextFontWeight', fWeight, ...
      'DefaultTextUnits', fUnits );
    
    % transform data to Cartesian coordinates.
   
    % Here convert from polar coordinates to geographical polar coords (N=0
    % deg with theta increasing in clock sense
    
    if DATA
      theta = -theta + pi/2;
  
      xx = rho .* cos(theta);
      yy = rho .* sin(theta);
    
      p_lo = 1:4:numel(xx);
      p_hi = 4:4:numel(xx);

      for p=1:numel(p_lo)
        q = patch(xx(p_lo(p):p_hi(p)), yy(p_lo(p):p_hi(p)), 'k', 'Parent', cax);
      end
    end
    
    if ~hold_state
      set(cax, 'DataAspectRatio', [1, 1, 1]), axis(cax, 'off');
      set(cax, 'NextPlot', next);
    end
    set(get(cax, 'XLabel'), 'Visible', 'on');
    set(get(cax, 'YLabel'), 'Visible', 'on');
end
