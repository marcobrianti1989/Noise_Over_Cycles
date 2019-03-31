clear
close all
Starting_date = 1982;
trend = 'quad';
nlags  = 8; %8
% Read main dataset
filename                    = 'main_file';
sheet                       = 'Sheet1';
range                       = 'B1:CT300';
do_truncation               = 0; %Do not truncate data. You will have many NaN
[dataset, var_names]        = read_data2(filename, sheet, range, do_truncation);

%nNaN                        = 20; % adding some NaN at the end to have space for leads
%dataset                     = [dataset; NaN(nNaN,size(dataset,2))];
% Assess name to each variable
for i = 1:size(dataset,2)
      eval([var_names{i} ' = dataset(:,i);']);
end

start_time = find(Time==Starting_date);

%system = [TFP RealGDP RealCons RealInvestment Hours]; 
%VAR_VEC = [TFP RealGDP RealCons RealInvestment];  
VAR_VEC = [TFP PC1];  %[PC4]
[system, first, last] = truncate_data(VAR_VEC);
VAR_VEC = VAR_VEC(start_time:last,:); 
dTFP = detrend_func(VAR_VEC(:,1),trend);
VAR_VEC(:,1) = dTFP;
[A,B,res,sigma] = sr_var(VAR_VEC,nlags);
TechSurprise = [NaN(start_time-1+nlags,1); res(:,1); NaN(length(Time)-last,1)]; 
save('TechShock_identification.mat','TechSurprise');