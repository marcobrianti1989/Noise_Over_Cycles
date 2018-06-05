%provs
function [y, yss] = generate_AR1(rho,sigm,T,drift,tt)
%This function generates an AR(1) process with persistence rho, 
%a possible drift different from zero, a possible linear time trand (tt) and
%disturbances with variance sigm. Length of the sample is T and initial
%values is y0.

y0 = drift/(1-rho); %Initial value is y in steady state
yss = y0; %Set as an output also y in steady state
y = zeros(1,T); 
errors = sigm*randn(1,T);

y(1) = rho*y0 + errors(1);
for i = 2:T
    y(i) = drift + tt*i + rho*y(i-1) + errors(i);
end

end

