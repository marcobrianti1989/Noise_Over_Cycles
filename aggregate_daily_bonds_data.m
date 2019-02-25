

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Marco Brianti and Vito Cormun - February 22, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add path to get Data
base_path = pwd;
if exist([base_path '\Data'], 'dir')
      addpath([base_path '\Data']) %for Microsoft
else
      addpath([base_path '/Data']) %for Mac
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Rating Datasets - There are 3 of them since they are too long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file                     = 'rating';
sheet                    = 'WRDS';
range                    = 'A2:C1039062';
[data_ratings1, rating]   = xlsread(file,sheet,range);
code_rating              = data_ratings1(:,1);
date_rating              = data_ratings1(:,3);

file                     = 'rating2';
sheet                    = 'WRDS';
range                    = 'A2:C1048266';
[data_ratings2, rating2] = xlsread(file,sheet,range);
code_rating              = [code_rating; data_ratings2(:,1)];
date_rating              = [date_rating; data_ratings2(:,3)];
rating                   = [rating; rating2];

file                     = 'rating3';
sheet                    = 'WRDS';
range                    = 'A2:C949563';
[data_ratings3, rating3] = xlsread(file,sheet,range);
code_rating              = [code_rating; data_ratings3(:,1)];
date_rating              = [date_rating; data_ratings3(:,3)];
rating                   = [rating; rating3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Principal Datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file                     = 'principal';
sheet                    = 'WRDS';
range                    = 'A2:C414775';
[data_principal, ~]      = xlsread(file,sheet,range);
code_principal           = data_principal(:,1);
date_principal           = data_principal(:,2);
principal                = data_principal(:,3);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distinguish between Junk and Safe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JUNK_CODE = {'BB+','BB','BB-','B+','B','B-','CCC+','CCC','CCC-','CC','C',...
      'RD','SD','D','DDD','DD',...
      'Ba1','Ba2','Ba3','B1','B2','B3','Caa1','Caa2','Caa3','Ca'};
NOT_RATED = {'NR'};
JUNK = zeros(size(code_principal,1),1);
SAFE = zeros(size(code_principal,1),1);
for i = 1:length(code_principal)
      codei = code_principal(i);
      locs  = find(codei == code_rating);
      if isempty(locs) == 1
            JUNK(i,1) = 0;
            SAFE(i,1) = 0;
      else
            for il = 1:length(locs)
                  locsi                 = locs(il);
                  ratinglocsi           = rating{locsi};
                  JUNK_flagil           = strcmp(JUNK_CODE,ratinglocsi);
                  NR_flagil             = strcmp(NOT_RATED,ratinglocsi);
                  JUNK_flag(il,1)       = sum(JUNK_flagil);
                  NR_flag(il,1)         = sum(NR_flagil);
            end
            JUNK_flag_SUM      = sum(JUNK_flag);
            NR_flag_SUM        = sum(NR_flag);
            if JUNK_flag_SUM > 0 
                  JUNK(i,1) = principal(i,1);
                  SAFE(i,1) = 0;
                  if i/100 == ceil(i/100)
                        disp(['Observation ',num2str(i),' is pure junk.'])
                  end
            elseif JUNK_flag_SUM == 0 && NR_flag_SUM == length(locs)
                  JUNK(i,1) = 0;
                  SAFE(i,1) = 0;
                  if i/100 == ceil(i/100)
                        disp(['Observation ',num2str(i),' is not rated.'])
                  end
            else
                  JUNK(i,1) = 0;
                  SAFE(i,1) = principal(i,1);
                  if i/100 == ceil(i/100)
                        disp(['Observation ',num2str(i),' is pure safe.'])
                  end
            end
      end
      clear JUNK_flag NR_flag JUNK_flag_SUM NR_flag_SUM
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aggregate over Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bundle_daily = [date_principal JUNK SAFE];
idcount = 1;
imcount = 1;
iycount = 1;
for iy = 1970:2018
      JUNK_month_sum = 0;
      SAFE_month_sum = 0;
      for im = 1:12
            JUNK_day_sum = 0;
            SAFE_day_sum = 0;
            for id = 1:31
                  if id < 10 && im < 10
                        day_name = [num2str(iy),'0',num2str(im),'0',num2str(id)];
                        eval(['day = ',day_name,';']);
                  elseif id < 10 && im >=10
                        day_name = [num2str(iy),'0',num2str(im),num2str(id)];
                        eval(['day = ',day_name,';']);
                  elseif id >= 10 && im < 10
                        day_name = [num2str(iy),num2str(im),'0',num2str(id)];
                        eval(['day = ',day_name,';']);
                  else
                        day_name = [num2str(iy),num2str(im),num2str(id)];
                        eval(['day = ',day_name,';']);
                  end
                  locsday = find(date_principal == day);
                  if isempty(locsday) == 1
                        JUNK_day = 0;
                        SAFE_day = 0;
                  else
                        JUNK_day = sum(JUNK(locsday));
                        SAFE_day = sum(SAFE(locsday));
                  end
                  JUNK_day_sum                = JUNK_day_sum + JUNK_day;
                  SAFE_day_sum                = SAFE_day_sum + SAFE_day;
                  JUNK_day_STORE(idcount,1)   = JUNK_day;
                  SAFE_day_STORE(idcount,1)   = SAFE_day;
                  idcount                     = idcount + 1;
            end
            JUNK_month                    = JUNK_day_sum;
            SAFE_month                    = SAFE_day_sum;
            JUNK_month_STORE(imcount,1)   = JUNK_month;
            SAFE_month_STORE(imcount,1)   = SAFE_month;
            JUNK_month_sum                = JUNK_month_sum + JUNK_month;
            SAFE_month_sum                = SAFE_month_sum + SAFE_month;
            imcount                       = imcount + 1;
      end
      JUNK_year                    = JUNK_month_sum;
      JUNK_year_STORE(iycount,1)   = JUNK_year;
      SAFE_year                    = SAFE_month_sum;
      SAFE_year_STORE(iycount,1)   = SAFE_year;
      iycount                      = iycount + 1
end
% Evaluate monthly Junk Share
JUNK_SHARE_month  = log(JUNK_month_STORE)./(log(JUNK_month_STORE) + log(SAFE_month_STORE));
monthly           = [1970:(1/12):(2019-(1/12))]';
% Evaluate annual Junk Share
JUNK_SHARE_year   = log(JUNK_year_STORE)./(log(JUNK_year_STORE) + log(SAFE_year_STORE));
yearly            = [1970:1:2018]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aggregate monthly to Quarterly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iquarter = 1;
im       = 1;
sumim_JUNK = 0;
sumim_SAFE = 0;
for im = 1:length(JUNK_month_STORE)
      sumim_JUNK = sumim_JUNK + JUNK_month_STORE(im);
      sumim_SAFE = sumim_SAFE + SAFE_month_STORE(im);
      if im/3 == ceil(im/3)
            JUNK_quarter_STORE(iquarter,1) = sumim_JUNK;
            SAFE_quarter_STORE(iquarter,1) = sumim_SAFE;
            iquarter                       = iquarter + 1;
            sumim_JUNK                     = 0;
            sumim_SAFE                     = 0;
      end
end
% Evaluate quarterly Junk Share
quarterly            = [1970:0.25:2018.75]';
JUNK_SHARE_quarter   = log(JUNK_quarter_STORE)./(log(JUNK_quarter_STORE) + log(SAFE_quarter_STORE));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plut Junk Shares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(3,1,1)
plot(monthly,JUNK_SHARE_month)
subplot(3,1,2)
plot(quarterly,JUNK_SHARE_quarter)
subplot(3,1,3)
plot(yearly,JUNK_SHARE_year)








