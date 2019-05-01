
clear
close all

% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:DX300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);

% Assess name to each variable
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end


omeg = linspace(0,pi,50)';
Z    = exp(-1i*omeg);
lags = 25;
X    = UnempRate;
X    = truncate_data(X);
for iZ = 1:length(omeg)
      [GxZ(iZ,1), Cx] = autocov_GenFunc(X,lags,Z(iZ));
end
periodicity = (2*pi./omeg);
figure(1)
plot(omeg,GxZ)
figure(2)
plot(periodicity,GxZ)

