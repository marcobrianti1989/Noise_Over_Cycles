clear
close all

% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:CT300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);

nNaN                        = 20; % adding some NaN at the end to have space for leads
dataset                     = [dataset; NaN(nNaN,size(dataset,2))];
% Assess name to each variable
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

system = [TFP PC1 PC2 PC3 PC4 PC5 PC6];
nlags  = 8;
[system, truncation_point, truncation_point2] = truncate_data(system);
[A,B,res,sigma] = sr_var(system,nlags);