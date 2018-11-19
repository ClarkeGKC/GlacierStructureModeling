function PLOT = InitializePlotStructure(PlotFileName, MODEL)

f_dat = fopen(PlotFileName, 'r');

if f_dat==-1
  error('InitializePlotStructure(): Failed to open "%s"\n', PlotFileName)
end  

while 1
  if feof(f_dat)
    break
  end  
  str = fgetl(f_dat);
  if ischar(str) && numel(str)>0
    if ~strcmp(str(end), ';')
      str = strcat(str, ';');
    end
    eval(str)
  end
end

PLOT = struct('Map', [], 'Section', [], 'Profile', [], 'Density', []);

if exist('Map', 'var')
  PLOT.Map = Map;
end

if exist('Section', 'var')
  PLOT.Section = Section;
end

if exist('Profile', 'var')
  PLOT.Profile = Profile;
end

if exist('Density', 'var')
  PLOT.Density = Density;
end

if exist('Stereo', 'var')
  PLOT.Stereo  = Stereo;
end  

if exist('Ellipsoid', 'var')
  PLOT.Ellipsoid = Ellipsoid;  
end

if exist('TimeSeries', 'var')
  PLOT.TimeSeries = TimeSeries;
end   

if exist('Graph', 'var')
  PLOT.Graph = Graph;
end 

if exist('Check', 'var')
  PLOT.Check = Check;
end  

if exist('Trajectory', 'var')
  PLOT.Trajectory = Trajectory;
end  
  
  
