function [yr_start, yr_display] = GetMoraineDataList(file_dat)

f_dat = fopen(file_dat, 'r');

C     = textscan(f_dat, '%f%f', 'headerlines', 1);
yr_start   = C{1};
yr_display = C{2};

fclose(f_dat);