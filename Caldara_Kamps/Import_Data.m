clear all; 
close all;
																	
DATASET.TSERIES=xlsread('CK_RESTUD_DATASET.xlsx','FINAL');

DATASET.LABEL     = {'DATES','GDP','TAX','G','TB3MS','CPI_PIQ4','MUNI1Y','PDVMILY','DTFP_UTIL','HAMILTON3YP','RESID08','TAXNARRATIVE'};
DATASET.VALUE	= [  1	  ,    2,    3,	 4,	5     ,	   6     ,	 7    ,	       8,	       9,	  10	  , 11      ,   12	  ];	  
DATASET.UNIT	= [  0	  ,    2,    2,	 2,	1     ,	   1     ,	 1    ,	       1,	       1,	   1	  , 1       ,	1	  ];	    

% Detrend Non-stationary Series

TSERIESDETREND   = detrend(DATASET.TSERIES(:,DATASET.UNIT==2));
DATASET.TSERIES  = [DATASET.TSERIES TSERIESDETREND];

DATASET.LABEL     = [DATASET.LABEL,{'GDP_S','TAX_S','G_S'}];
DATASET.VALUE	= [DATASET.VALUE,  13,        14 , 15	];
DATASET.UNIT	= [DATASET.UNIT,	1,         1 ,  1	];

DATASET.MAP = containers.Map(DATASET.LABEL,DATASET.VALUE);

temp = DATASET.TSERIES(:,DATASET.UNIT==2);
DATASET.GDPRATIO = zeros(size(temp,2),1);
for ii = 1:size(temp,2);
    DATASET.GDPRATIO(ii,1) = mean(exp(temp(:,ii))./exp(temp(:,1)));
end

DATASET.FIGLABELS{1,1}  = 'Year';
DATASET.FIGLABELS{2,1}  = 'Real GDP';                                          
DATASET.FIGLABELS{3,1}  = 'Tax Revenue';                                       
DATASET.FIGLABELS{4,1}  = 'Government Spending';                               
DATASET.FIGLABELS{5,1}  = '3-Month Yield';                                     
DATASET.FIGLABELS{6,1}  = 'CPI Inflation (Q4/Q4)';                             
DATASET.FIGLABELS{7,1}  = 'Tax News Shocks';					    
DATASET.FIGLABELS{8,1}  = 'Defense Spending News';					   
DATASET.FIGLABELS{9,1}  = 'Utilization-adjusted TFP';					    
DATASET.FIGLABELS{10,1} = 'Oil Price Shocks';					  
DATASET.FIGLABELS{11,1} = 'Narrative Tax Shocks';	  
DATASET.FIGLABELS{12,1} = 'Romer and Romer Monetary Policy Shocks';   
DATASET.FIGLABELS{13,1} = 'Real GDP (Detrended)';					     
DATASET.FIGLABELS{14,1} = 'Tax Revenue (Detrended)';			     
DATASET.FIGLABELS{15,1} = 'Government Spending (Detrended)';		     

clear TSERIESDETREND

save('DATASET','DATASET');
