function [Grid, t_site, t_last_surge, X_NAD27_site, Y_NAD27_site, X_WGS84_site, Y_WGS84_site, XI_site, dZ_tweak_site, Name_site, Type_site, Group] = GetSiteList(SiteFileName, MODEL)

f_dat = fopen(SiteFileName, 'r');

while 1
  if feof(f_dat)
    break
  end  
  str = fgetl(f_dat);
  if ischar(str) && numel(str)>0
    if ~strcmp(str(end), ';')
      str = strcat(str, ';');
    end
    try
      eval(str)
    catch
      fprintf(1, '%s\n',str);
      error('FAIL')
    end
  end
end

if ~exist('t_site', 'var')
  error('GetSiteList(): No value assigned to "t_site"')
end

if ~exist('t_last_surge', 'var')
  t_last_surge = [];
elseif isempty(t_last_surge.start) || isempty(t_last_surge.stop)
  error('GetSiteList(): No value(s) assigned to one or all of "t_last_surge.start", "t_last_surge.peak", "t_last_surge.stop"')
end  

% Verify that the last surge times are consistent with the MODEL-specified surge cycle

if ~isempty(t_last_surge)
  if t_last_surge.stop-t_last_surge.start ~=MODEL.dt_SURGE
    error('GetSiteList(): "t_last_surge.stop-t_last_surge.start=%.2f " is not consistent with MODEL-specified "dt_SURGE=%.2f"', ...
      t_last_surge.stop-t_last_surge.start, MODEL.dt_SURGE);
  end
end

nl = numel(X_site); 

XX_site = cell2mat(X_site{1});
YY_site = cell2mat(Y_site{1});

for n=2:nl
  XX_site = [XX_site cell2mat(X_site{n})];
end

for n=2:nl
  YY_site = [YY_site cell2mat(Y_site{n})];
end

X_site = XX_site;
Y_site = YY_site;

ns      = numel(X_site);

if exist('XI_site', 'var')
  xi_site = cell2mat(XI_site{1});
  for n=2:nl
    xi_site = [xi_site cell2mat(XI_site{n})];
  end
  XI_site = xi_site;
else
  XI_site = ones(1, ns);
end

if exist('dZ_tweak_site', 'var')
  dz_tweak_site = cell2mat(dZ_tweak_site{1});
  for n=2:nl
    dz_tweak_site = [dz_tweak_site  cell2mat(dZ_tweak_site{n})];
  end
  dZ_tweak_site = dz_tweak_site;
else
  dZ_tweak_site = zeros(1, ns);
end

if strcmp(Grid, 'WGS84')
  X_WGS84_site = X_site;
  Y_WGS84_site = Y_site;
           
  X_NAD27_site = X_WGS84_site + 110.557;
  Y_NAD27_site = Y_WGS84_site -165.460;  
else
  X_NAD27_site = X_site;
  Y_NAD27_site = Y_site;

  X_WGS84_site = X_site - 110.557;
  Y_WGS84_site = Y_site + 165.460;
end

if exist('Name_site', 'var')
  NN_site = Name_site{1};
  for n=2:nl
   NN_site = [NN_site Name_site{n}];
  end
  Name_site = NN_site;
else
  Name_site = cell(1, ns);
  for s=1:ns
    Name_site{s} = num2str(s);
  end
end

if exist('Type_site', 'var')
  TT_site = Type_site{1};
  for n=2:numel(Type_site)
    TT_site = [TT_site Type_site{n}];
  end
  Type_site = TT_site;
else
  Type_site = cell(1, ns);
  for s=1:ns
    Type_site{s} = 'A';
   end
end  
