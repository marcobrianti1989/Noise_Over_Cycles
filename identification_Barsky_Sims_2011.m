function [impact, gam] = ...
      identification_Barsky_Sims_2011(A,B,horizon,TFPposition)

% A chol impact matrix (nvar,nvar)
% Reduce-form regressors coefficient (1+nvar*nlags,nvar)
% horizon is which horizon to max news shock
% TFPposition for key variable to max

% Technical values
B                     = B(2:end,:); % remove the constant
[nvarlags, nvar]      = size(B);
nlags                 = nvarlags/nvar;

% Defining initial values
D              = eye(nvar);
gam1_zero      = D(:,1); %financial shock impact vector (initial value)
gam2_zero      = D(:,2); %uncertainty shock impact vector (initial value)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First Step - Identifying gam1 - News TFP Shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the objective function of Step1
obj1 = @(gam1) objective_barskysims(TFPposition,horizon,B,A,gam1);

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

%Constraint that News shocks have no contemporaneous effect on TFP
% Me3       = 1; %  Me = no. of equality constraints
% Beq3      = zeros(Me3,1); % Beq is (Me x 1) where
% Aeq3      = zeros(Me3,1*nvar); % Aeq is (Me x (nshock*nvar)) - nshock is 1 at this step
% Aeq3(1,1) = 1; %zero-impact of news on TFP

% Optimization
[gam1_opt] = fmincon(obj1, gam1_zero,[],[],[],[],[],[],...
      @(gam1) constraint_orthogonality(gam1),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

if (gam1_opt'*gam1_opt - 1)^2 > 10^(-10)
      warning('The problem is not consistent with the constraints.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Step - Identifying gam2 - Surprise TFP Shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the objective function of Step1
obj2 = @(gam2) objective_max_impact(TFPposition,A,gam2);

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

%Optimization - Notice the contraint for gamNews and gamTFP
gam2_opt = fmincon(obj2, gam2_zero,[],[],[],[],[],[],...
      @(gam2) constraint_orthogonality([gam1_opt gam2]),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

if (gam2_opt'*gam2_opt - 1)^2 > 10^(-10) || sum(sum(([gam1_opt gam2_opt]'*[gam1_opt gam2_opt] - eye(2)).^2)) > 10^(-10)
      warning('The problem is not consistent with the constraints.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Step - Obtain impact matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gam          = [gam1_opt gam2_opt];
impact1      = A*gam;
impact       = [impact1 zeros(size(impact1,1),size(impact1,1)-size(impact1,2))];

end
