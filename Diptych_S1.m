function handle_out = Diptych_S1(handle, MODEL, SITE, N_pts, x_p, y_p, z_p, L_use_p, Color_p, n_vel_site, PlotTitle, OPTIONS, VERBOSITY)

if isfield(OPTIONS, 'HardCopy')
  HardCopy = OPTIONS.HardCopy;
else
  HardCopy = false;
end

if isfield(OPTIONS, 'dir_PLOT')
  dir_PLOT = OPTIONS.dir_PLOT;
else
  dir_PLOT = {};
end  

if isfield(OPTIONS, 'AzLabelFontScale')
  OPTIONS.AzLabelFontScale = 0.65*OPTIONS.AzLabelFontScale;
end

if isfield(OPTIONS, 'TitleFontScale')
  OPTIONS.TitleFontScale = 1.00*OPTIONS.TitleFontScale;
end

if isfield(OPTIONS, 'TextFontScale')
  OPTIONS.TextFontScale = 0.75*OPTIONS.TextFontScale;
end
if isfield(OPTIONS, 'MarkerScale')
  OPTIONS.MarkerScale = 0.5*OPTIONS.MarkerScale;
end

if isfield(OPTIONS, 'SurgeMarkers')
  OPTIONS.SurgeMarkers = 0;
end  

np      = numel(x_p);

x_p = x_p(L_use_p);
y_p = y_p(L_use_p);
z_p = z_p(L_use_p);

Color_p  = Color_p(L_use_p);

handle = handle+1;
figure(handle);

set(handle, 'units', 'inches')
fig_pos_start = get(handle, 'position');

hp =  uipanel('BorderType', 'line', 'BorderWidth', 1, 'HighLightColor', 'black', 'BackgroundColor','white','Position', [0.25 .15 .55 .30]);

set(hp, 'units', 'inches')

if ~isempty(PlotTitle)
  title(PlotTitle, 'interpreter', 'none')
  set(gca, 'XTick', [], 'YTick', [], 'box', 'off')
end  

pos_start = get(hp, 'position');

% Order of panel plotting has been re-ordered to eliminate blanking by overprinting

g9 = subplot(3,3,9,'Parent', hp);

axis off

TAB = 0.115;

p9 = get(g9, 'position');
p9(1) = 0.86;
p9(2) = 0.035;
p9(3) = 0.75;
P9(4) = 0.05; % 10;
set(g9, 'position', p9)

ax9 = gca;
set(ax9, 'Title', [])

text(0, 0.25, sprintf('N=%d', numel(x_p)), 'FontSize', OPTIONS.TitleFontScale*9)

g3 = subplot(3,3,3, 'Parent', hp);

axis off

p3 = get(g3, 'position');
p3(1) = 0.86;
p3(2) = 0.90;
p3(3) = 0.75;
p3(4) = 0.05;

set(g3, 'position', p3)
ax3   = gca;
set(ax3, 'Title', [])

text(0, 0.25, sprintf('(n=%d)', N_pts), 'FontSize', OPTIONS.TitleFontScale*6)

g1 = subplot(3,3,1, 'Parent', hp);
axis off

p1 = get(g1, 'position');

p1(1) = p1(1)-TAB;  % shift left
p1(2) = 0.89;
p1(3) = p1(3)+TAB;  % donate shift space to panel width
p1(4) = 0.10;

set(g1, 'position', p1)

text(0, 0.25, sprintf('%s', char(SITE.Name_site)), 'FontWeight', 'bold', 'FontSize', OPTIONS.TitleFontScale*10)

g4 = subplot(3,3,4, 'Parent', hp);
axis off

p4 = get(g4, 'position');

p4(1)   = 0.01;
p4(2)   = 0.15;
p4(3)   = 0.35;
p4(4)   = 0.70;

set(g4, 'position', p4)

PlotSchmidt_S1([], MODEL, SITE, [], x_p, y_p, z_p, Color_p, n_vel_site, OPTIONS);
ax4 = gca;
set(ax4,'Title',[])

g5 = subplot(3,3,5, 'Parent', hp);
axis off

p5 = get(g5, 'position');
p5(1) = 0.333;
p5(2) = p4(2);
p5(3) = p4(3);
p5(4) = p4(4);
set(g5, 'position', p5);

OPTIONS.Monochrome = 1;

NewPlotDensityContours([], MODEL, SITE, x_p, y_p, z_p, n_vel_site, OPTIONS);
ax5 = gca;
set(ax5,'Title',[])

if any(fig_pos_start-get(handle, 'position'))
  fprintf(1,'Diptych_2(1): Figure changes size for Site %s\n', char(SITE.Name_site))
end

if HardCopy && ~isempty(handle)
  if isempty(dir_PLOT)
    dir_PLOT = pwd;
  end  
  
  if ~exist(dir_PLOT, 'dir')
    mkdir(dir_PLOT)
  end  

  set(handle, 'PaperPositionMode', 'auto')
  print(handle, fullfile(dir_PLOT, sprintf('Fig_%s.pdf', char(SITE.Name_site))), '-dpdf')    
end  

handle_out = handle;

