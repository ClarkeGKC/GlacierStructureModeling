function YearList = ExtractYearListFromFileList(FileList)

N = numel(FileList);

YearList = zeros(N, 1);

for n=1:N
  
  FileName = FileList{n};
  L_par = strfind(FileName, '(');
  R_par = strfind(FileName, ')');
  
  YearList(n) = str2num(FileName(L_par+1:R_par-1)); 
end  