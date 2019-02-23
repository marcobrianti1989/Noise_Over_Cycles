% clear
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          Marco Brianti and Vito Cormun - February 22, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Add path to get Data
% base_path = pwd;
% if exist([base_path '\Data'], 'dir')
%       addpath([base_path '\Data']) %for Microsoft
% else
%       addpath([base_path '/Data']) %for Mac
% end
% 
% % Read Rating Datasets - There are 3 of them since they are too long for
% % Excel
% file                     = 'rating';
% sheet                    = 'WRDS';
% range                    = 'A2:C1039062';
% [data_ratings1, rating]   = xlsread(file,sheet,range);
% code_rating              = data_ratings1(:,1);
% date_rating              = data_ratings1(:,3);
% 
% file                     = 'rating2';
% sheet                    = 'WRDS';
% range                    = 'A2:C1048266';
% [data_ratings2, rating2] = xlsread(file,sheet,range);
% code_rating              = [code_rating; data_ratings2(:,1)];
% date_rating              = [date_rating; data_ratings2(:,3)];
% rating                   = [rating; rating2];
% 
% file                     = 'rating3';
% sheet                    = 'WRDS';
% range                    = 'A2:C949563';
% [data_ratings3, rating3] = xlsread(file,sheet,range);
% code_rating              = [code_rating; data_ratings3(:,1)];
% date_rating              = [date_rating; data_ratings3(:,3)];
% rating                   = [rating; rating3];
% 
% % Read Principal Datasets
% file                     = 'principal';
% sheet                    = 'WRDS';
% range                    = 'A2:C414775';
% [data_principal, ~]      = xlsread(file,sheet,range);
% code_principal           = data_principal(:,1);
% date_principal           = data_principal(:,2);
% principal                = data_principal(:,3);
tic
JUNK_CODE = {'BB+','BB','BB-','B+','B','B-','CCC+','CCC','CCC-','CC','C',...
      'RD','SD','D','DDD','DD',...
      'Ba1','Ba2','Ba3','B1','B2','B3','Caa1','Caa2','Caa3','Ca'};

JUNK = zeros(size(code_principal,1),1);
SAFE = zeros(size(code_principal,1),1);
for i = 1:length(code_principal)
      codei = code_principal(i);
      locs  = find(codei == code_rating);
      for il = 1:length(locs)
         locsi           = locs(il);
         ratinglocsi     = rating{locsi};
         JUNK_flag(il)   = sum(strcmp(JUNK_CODE,ratinglocsi));
      end
      BOND_flag          = sum(JUNK_flag);     
      if BOND_flag >= 1
            JUNK(i,1) = principal(i,1);
            SAFE(i,1) = 0;
      else
            JUNK(i,1) = 0;
            SAFE(i,1) = principal(i,1);
      end  
      if i/100 == ceil(i/100)
            disp(['Iteration number ', num2str(i)])
      end
end
toc

boundle = [JUNK SAFE]
















