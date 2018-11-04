
function [Data, codeData,tsc] = ReadDataNews_EJ_data

[x,text] = xlsread('dataset_EJ_data.xlsx','Foglio1','B6:DH213');
tsc = xlsread('dataset_EJ_data.xlsx','Foglio1','B4:DH4');% Row 4 describes if series are in levels or in growth rates

for i=1:size(x,2)
    if     tsc(i) == 0;% levels
        data(:,i) = x(1:end,i);
        codeData(i)= 1;
    elseif tsc(i) == 1;% logs
        data(:,i) = 100*log(x(1:end,i));
        codeData(i) = 0;
    end
end
Data = data;