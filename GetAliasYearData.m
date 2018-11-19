function   [DD, V, LL, EE, T, YearSwap] = GetAliasYearData(YearGet, t_last_surge, dir_arch, ModelFileList, ModelYearList, MODEL, VERBOSITY)

% Use alias years to match stop year with an output file

if nargin==7
  VERBOSE = VERBOSITY;
else
  VERBOSE = 1;
end

[FileSwap, YearSwap] = GetYearSwapFileName(MODEL, t_last_surge, ModelFileList, ModelYearList, YearGet);

YearShift   = YearGet-YearSwap;

FullFileSwap = fullfile(dir_arch, FileSwap);

load(FullFileSwap, 'DD', 'V', 'LL', 'EE', 'T')

DD.cal_year = YearGet;

if VERBOSE
  fprintf(1,'GetAliasYearData(): YearGet (t) = %.2f yr; YearSwap (cal_year) = %.2f yr; YearShift = %.2f yr\n', YearGet, YearSwap, YearShift);
end  
