function BLANKING_MODEL = GetBlankingModelParameters(BlankingRule, BlankingFileName)

if BlankingRule==0
  BLANKING_MODEL.Rule     = 0;
  BLANKING_MODEL.str      = 'No Blanking';
  BLANKING_MODEL.BlankDeg = [];
  BLANKING_MDOEL.Distance = [];
  BLANKING_MODEL.Time     = [];
  BLANKING_MODEL.ParamStr = {};
else
  BLANKING_MODEL.Rule = BlankingRule;
  f_dat = fopen(BlankingFileName, 'r');
  
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
  
  switch BlankingRule
    case 1  % Point blanking
      BLANKING_MODEL.str         = BlankingString{1};
      BLANKING_MODEL.BlankArcDeg = BlankingArcDeg(1);
      BLANKING_MODEL.Distance    = [];
      BLANKING_MODEL.Time        = [];
      BLANKING_MODEL.ParamStr    = sprintf('Blanking angular radius=%.2f deg', BlankArcDeg(1));
    case 2  % Distance blanking
      BLANKING_MODEL.str         = BlankingString{2};
      BLANKING_MODEL.BlankArcDeg = [];
      BLANKING_MODEL.Distance    = BlankingDistance(2);
      BLANKING_MODEL.Time        = [];
      BLANKING_MODEL.ParamStr    = sprintf('Blanking distance=%.2f m', BlankingDistance(2));
    case 3  % Time blanking
      BLANKING_MODEL.str         = BlankingString{3};
      BLANKING_MODEL.BlankArcDeg = [];
      BLANKING_MODEL.Distance    = [];
      BLANKING_MODEL.Time        = BlankingTime(3);
      BLANKING_MODEL.ParamStr    = sprintf('Blanking time=%.2f yr', BlankingTime(3));
    case 4  % Space OR time blanking
      BLANKING_MODEL.str         = BlankingString{4};
      BLANKING_MODEL.BlankArcDeg = [];
      BLANKING_MODEL.Distance    = BlankingDistance(4);
      BLANKING_MODEL.Time        = BlankingTime(4);
      BLANKING_MODEL.ParamStr    = sprintf('Blanking distance=%.2f m; blanking time=%.2f yr', BlankingDistance(4), BlankingTime(4));     
    case 5  % Space AND time blanking
      BLANKING_MODEL.str         = BlankingString{5};
      BLANKING_MODEL.BlankArcDeg = [];
      BLANKING_MODEL.Distance    = BlankingDistance(5);
      BLANKING_MODEL.Time        = BlankingTime(5);
      BLANKING_MODEL.ParamStr    = sprintf('Blanking distance=%.2f m; blanking time=%.2f yr', BlankingDistance(5), BlankingTime(5));      
    otherwise
      error('GetBlankingModelParameters(): Unprogrammed blanking rule')
  end
end      

