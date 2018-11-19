function CRACK_MODEL = GetCrackModelParameters(CrackFileName)

f_dat = fopen(CrackFileName, 'r');

while 1
  str = fgetl(f_dat);
  if ~ischar(str) || numel(str)<2
    break
  end
  if ~strcmp(str(end), ';')
    str = strcat(str, ';');
  end
  eval(str)
end

CRACK_MODEL.Text       = CrackModelText;
CRACK_MODEL.BIN_DEF    = CRACK_BIN_DEF;

if exist('min_cos_n_S', 'var')
  CRACK_MODEL.min_cos_n_S = min_cos_n_S;
else
  CRACK_MODEL.min_cos_n_S = [];
end  

if exist('BLANKING_RULE', 'var')
  CRACK_MODEL.BLANKING_RULE = BLANKING_RULE;
else  
  CRACK_MODEL.BLANKING_RULE = 0;
end

switch CRACK_BIN_DEF
  case 'Hambrey'
    CRACK_MODEL.BinStr = 'Hambrey crack bin definitions';
  case 'Equal'
    CRACK_MODEL.BinStr = 'Crack bins have equal angular span';
  otherwise
    error('GetCrackModelParameters(): Unexpected value CRACK_BIN_DEF')
end

CrackString            = CrackString{CRACK_RULE};
CRACK_MODEL.RULE       = CRACK_RULE;
CRACK_MODEL.RuleStr    = CrackString;

CRACK_MODEL.D_Ext_threshold     = [];  % strain-rate threshold for crack formation ("Ext" refers to extension/tensile crack)
CRACK_MODEL.SIGMA_Ext_threshold = [];
CRACK_MODEL.K_threshold         = [];
CRACK_MODEL.f_0                 = [];     
CRACK_MODEL.f_min               = [];
CRACK_MODEL.f_max               = [];
CRACK_MODEL.spacing             = [];

switch CRACK_RULE
  case 1
    CRACK_MODEL.ParamStr            = {sprintf('D_Ext_threshold=%.3f 1/yr', D_Ext_threshold)};
    CRACK_MODEL.D_Ext_threshold     = D_Ext_threshold;
  case 2    
    CRACK_MODEL.ParamStr            = {sprintf('SIGMA_Ext_theshold=%.0f Pa', SIGMA_Ext_threshold)};
    CRACK_MODEL.SIGMA_Ext_threshold = SIGMA_Ext_threshold;
  case 3
    CRACK_MODEL.ParamStr            = {sprintf('K_threshold=%.0f Pa m^(1/2)', K_threshold)};
    CRACK_MODEL.K_threshold         = K_threshold;
  case 4
    CRACK_MODEL.ParamStr            = {sprintf('K_threshold(f_0)=%.0f Pa m^(1/2)', K_threshold), sprintf('f_0=%.3f',f_0)};
    CRACK_MODEL.K_threshold         = K_threshold;
    CRACK_MODEL.f_0                 = f_0;
  case 5
    CRACK_MODEL.ParamStr            = {sprintf('K_threshold(f)=%.0f Pa m^(1/2)', K_threshold), sprintf('f_min=%.3f',f_min), sprintf('f_max=%.3f',f_max)};
    CRACK_MODEL.K_threshold         = K_threshold;
    CRACK_MODEL.f_min               = f_min;
    CRACK_MODEL.f_max               = f_max;
  case 6
    CRACK_MODEL.ParamStr            = {sprintf('K_threshold=%.0f Pa m^(1/2)', K_threshold)};
    CRACK_MODEL.K_threshold         = K_threshold;
    CRACK_MODEL.spacing             = spacing;
  case 7
    CRACK_MODEL.ParamStr            = {sprintf('K_threshold(f_0)=%.0f Pa m^(1/2)', K_threshold), sprintf('f_0=%.3f',f_0)};
    CRACK_MODEL.K_threshold         = K_threshold;
    CRACK_MODEL.f_0                 = f_0;
    CRACK_MODEL.spacing             = spacing;    
  case 8
    CRACK_MODEL.ParamStr            = {sprintf('K_threshold(f)=%.0f Pa m^(1/2)', K_threshold), sprintf('f_min=%.3f',f_min), sprintf('f_max=%.3f',f_max)};
    CRACK_MODEL.K_threshold         = K_threshold;
    CRACK_MODEL.f_min               = f_min;
    CRACK_MODEL.f_max               = f_max;
    CRACK_MODEL.spacing             = spacing;    
  case 9 
    CRACK_MODEL.ParamStr            = {'None'};       % This is the Nye crevasse model. It has no free parameters
  otherwise
    error('GetCrackModelParameters(): Unprogrammed CRACK_RULE')
end

