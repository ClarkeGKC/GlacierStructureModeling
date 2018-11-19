function OK = CheckExistence(FileList)

OK = true;

for f=1:numel(FileList)
  if exist(FileList{f}, 'file')==0 && exist(FileList{f}, 'dir')==0
    fprintf(1,'CheckExistence(): Required data file or directory "%s" does not exist\n', FileList{f})
    OK = false;
  end
end


