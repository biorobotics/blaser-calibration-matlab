function [n] = count_lines(fname)
%COUNT_LINES Count number of lines in a file
%   Detailed explanation goes here

n = 0;
fid = fopen(fname,'r');
tline = fgetl(fid);
while ischar(tline)
  tline = fgetl(fid);
  n = n+1;
end
fclose(fid);
end

