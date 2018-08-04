clear 
close all

% Create the correct path
base_path = pwd;
if exist([base_path '\Data'], 'dir')
      addpath([base_path '\Data']) %for Microsoft
else
      addpath([base_path '/Data']) %for Mac
end

filename = 'GZ_Spread_Monthly';
sheet    = 'Sheet1';
range    = 'A1:C525';
[data, ~] = xlsread(filename, sheet, range);
j = 1;
for i = 1:length(data)
      if floor(((i - 1)/3)) == (i - 1)/3 && length(data) - i > 3
      GZ_Spread(j,1) = mean(data(i:i+3,2));
      j = j + 1;
      end
end
plot(GZ_Spread)