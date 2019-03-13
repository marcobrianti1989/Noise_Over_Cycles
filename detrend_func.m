function system = detrend_func(system,which_trend)
% System is a matrix of (T,N) where T is number of observations and N
% number of variables.
% which_trend cvan be BP, HP, lin, quad, none, demean
% For both HP and BP check the parameters are correct
for i = 1:size(system,2)
      switch which_trend
            case 'BP'
                  system(:,i) = bpass(system(:,i),2,36);
%                   if i == 1
%                         disp('Variables are BPfiltered. Check Frequencies.')
%                         fprintf('\n')
%                   end
            case 'HP'
                  [~,system(:,i)] = hpfilter(system(:,i),1600);
%                   if i == 1
%                         disp('Variables are HPfiltered. Check Frequencies.')
%                         fprintf('\n')
%                   end
            case 'lin'
                  quadratic = 0;
                  system(:,i) = lin_timetrend(system(:,i),quadratic);
%                   if i == 1
%                         disp('Variables are linearly detrended.')
%                         fprintf('\n')
%                   end
            case 'quad'
                  quadratic = 1;
                  system(:,i) = lin_timetrend(system(:,i),quadratic);
%                   if i == 1
%                         disp('Variables are quadratically detrended.')
%                         fprintf('\n')
%                   end
            case 'diff'
                  system(:,i) = [NaN; diff(system(:,i))];
%                   if i == 1
%                         disp('Variables are differentiated.')
%                         fprintf('\n')
%                   end
            case 'none'
                  system(:,i) = system(:,i);
%                   if i == 1
%                         disp('Variables are not detrended')
%                         fprintf('\n')
%                   end
            case 'demean'
                  means       = mean(system(:,i));
                  system(:,i) = system(:,i) - means;
%                   if i == 1
%                         disp('Variables are demeaned')
%                         fprintf('\n')
%                   end
      end
end



end