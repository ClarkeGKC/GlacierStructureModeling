function SITE = InitializeSiteStructure(t_site, t_last_surge, X_site, Y_site, XI_site, dZ_tweak_site, Name_site, Type_site, Group, b_dot, FileSwap, MODEL)

global X2 Y2

% Note that B is in the original NAD27 coordinate system when extracted from DEM file

load(MODEL.file_dem, 'B');

if MODEL.GRID_shift
  B = B - MODEL.z_0;
end

load(FileSwap, 'DD')

H   = zeros(MODEL.ny, MODEL.nx);
H(DD.ic_jc_I) = DD.H;
S   = B+H;

nx     = MODEL.nx;
ny     = MODEL.ny;

N_xy   = MODEL.N_xy;

XI     = MODEL.XI;
n_XI   = MODEL.n_XI;

dx     = MODEL.dx;
dy     = MODEL.dy;

load(MODEL.file_balance, 'b_dot')

NaNCnt  = [];
ArgList = {};

B_site    = interp2(X2, Y2, B, X_site, Y_site);  
H_site    = interp2(X2, Y2, H, X_site, Y_site);
S_site    = B_site + H_site;

Z_site    = B_site + XI_site.*H_site;
L_S       = Z_site>S_site;
L_B       = Z_site<B_site;
Z_site(L_S) = S_site(L_S);
Z_site(L_B) = B_site(L_B);
D_site      = S_site - Z_site;

b_dot_site = interp2(X2, Y2, b_dot, X_site, Y_site);

L_bad_site = b_dot_site>0 & D_site==0;

if any(L_bad_site)
  error('InitializeSiteStructure(): Illogical site data for Site %s. Surface site is in accumulation zone', Name_site{L_bad_site}) 
end

if any(Z_site<B_site)
  error('InitializeSiteStructure(): Illogical site data for Site %s. Subsurface site is below bed elevation', Name_site{Z_site<B_site}) 
end

for s=1:numel(X_site)
  SITE(s).Name_site = Name_site(s);  
  SITE(s).t_site    = t_site;  
  SITE(s).X_site    = X_site(s);
  SITE(s).Y_site    = Y_site(s);
  SITE(s).Z_site    = Z_site(s);
  SITE(s).XI_site   = XI_site(s);
  
  if ~exist('dZ_tweak_site', 'var') || isempty(dZ_tweak_site(s))
    dZ_tweak_site(s) = 0;
    Z_init_site(s)  = Z_site(s);
    XI_init_site(s) = XI_site(s);
  else
    Z_init_site(s)  = Z_site(s) + dZ_tweak_site(s);
    XI_init_site(s) = (Z_init_site(s)-B_site(s))/H_site(s);
  end
 
  SITE(s).Z_init_site  = Z_init_site(s);
  SITE(s).XI_init_site = XI_init_site(s);
  
  SITE(s).B_site    = B_site(s);
  SITE(s).H_site    = H_site(s);
  SITE(s).S_site    = S_site(s);
  SITE(s).b_dot_site = b_dot_site(s);  
  SITE(s).Type_site = Type_site(s);
  SITE(s).yr_last_surge = t_last_surge;
  Type_list        = Type_site{s}';
  Group_list       = [Group.ID]';
  [~, LOC]         = ismember(Type_list, Group_list);
  SITE(s).Group    = Group(LOC);
end

