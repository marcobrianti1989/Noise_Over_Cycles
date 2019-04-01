function [residual, g1, g2, g3] = Beaudry_Portier_2004_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
T18 = exp(y(13))*params(3)*params(8)*y(7)^(params(8)-1);
T35 = params(8)*params(3)*exp(y(18))*y(15)^(params(8)-1);
T46 = exp(y(17))*params(4)*y(16)^params(7);
T54 = params(5)*T46^params(6)+(1-params(5))*y(6)^params(6);
T57 = (1-params(5))*exp(y(14))*params(2)*T54^(-1);
T59 = y(6)^(params(6)-1);
T61 = params(9)*exp(y(14))*params(2)*(1-params(1))*T35^(-1)+T57*T59;
T68 = params(4)*exp(y(12))*y(8)^params(7);
T74 = params(5)*T68^params(6)+(1-params(5))*y(1)^params(6);
T76 = params(5)*T74^(-1);
T77 = T68^(params(6)-1);
T83 = y(8)^(params(7)-1);
lhs =params(9)*T18^(-1);
rhs =T61;
residual(1)= lhs-rhs;
lhs =params(9);
rhs =params(7)*params(4)*exp(y(12))*T76*T77*T83;
residual(2)= lhs-rhs;
lhs =y(5);
rhs =T74^(1/params(6));
residual(3)= lhs-rhs;
lhs =y(10);
rhs =y(7)+y(8);
residual(4)= lhs-rhs;
lhs =y(6);
rhs =(1-params(1))*y(1)+y(11);
residual(5)= lhs-rhs;
lhs =y(11);
rhs =params(3)*y(7)^params(8);
residual(6)= lhs-rhs;
lhs =y(9);
rhs =T68;
residual(7)= lhs-rhs;
lhs =y(12);
rhs =params(10)*y(2)+x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(13);
rhs =params(11)*y(3)+x(it_, 2);
residual(9)= lhs-rhs;
lhs =y(14);
rhs =params(12)*y(4)+x(it_, 3);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 21);

  %
  % Jacobian matrix
  %

T121 = (1-params(5))*getPowerDeriv(y(1),params(6),1);
T122 = getPowerDeriv(T74,(-1),1);
T131 = getPowerDeriv(T74,1/params(6),1);
T137 = getPowerDeriv(T54,(-1),1);
T148 = getPowerDeriv(T18,(-1),1);
T156 = getPowerDeriv(T35,(-1),1);
T161 = params(4)*exp(y(12))*getPowerDeriv(y(8),params(7),1);
T162 = getPowerDeriv(T68,params(6),1);
T167 = getPowerDeriv(T68,params(6)-1,1);
T185 = getPowerDeriv(T46,params(6),1);
  g1(1,6)=(-(T59*(1-params(5))*exp(y(14))*params(2)*(1-params(5))*getPowerDeriv(y(6),params(6),1)*T137+T57*getPowerDeriv(y(6),params(6)-1,1)));
  g1(1,7)=params(9)*exp(y(13))*params(3)*params(8)*getPowerDeriv(y(7),params(8)-1,1)*T148;
  g1(1,15)=(-(params(9)*exp(y(14))*params(2)*(1-params(1))*params(8)*params(3)*exp(y(18))*getPowerDeriv(y(15),params(8)-1,1)*T156));
  g1(1,16)=(-(T59*(1-params(5))*exp(y(14))*params(2)*T137*params(5)*exp(y(17))*params(4)*getPowerDeriv(y(16),params(7),1)*T185));
  g1(1,17)=(-(T59*(1-params(5))*exp(y(14))*params(2)*T137*params(5)*T46*T185));
  g1(1,13)=params(9)*T18*T148;
  g1(1,18)=(-(params(9)*exp(y(14))*params(2)*(1-params(1))*T35*T156));
  g1(1,14)=(-T61);
  g1(2,1)=(-(T83*params(7)*params(4)*exp(y(12))*T77*params(5)*T121*T122));
  g1(2,8)=(-(T83*params(7)*params(4)*exp(y(12))*(T77*params(5)*T122*params(5)*T161*T162+T76*T161*T167)+params(7)*params(4)*exp(y(12))*T76*T77*getPowerDeriv(y(8),params(7)-1,1)));
  g1(2,12)=(-(T83*params(7)*params(4)*(exp(y(12))*T76*T77+exp(y(12))*(T77*params(5)*T122*params(5)*T68*T162+T76*T68*T167))));
  g1(3,5)=1;
  g1(3,1)=(-(T121*T131));
  g1(3,8)=(-(T131*params(5)*T161*T162));
  g1(3,12)=(-(T131*params(5)*T68*T162));
  g1(4,7)=(-1);
  g1(4,8)=(-1);
  g1(4,10)=1;
  g1(5,1)=(-(1-params(1)));
  g1(5,6)=1;
  g1(5,11)=(-1);
  g1(6,7)=(-(params(3)*getPowerDeriv(y(7),params(8),1)));
  g1(6,11)=1;
  g1(7,8)=(-T161);
  g1(7,9)=1;
  g1(7,12)=(-T68);
  g1(8,2)=(-params(10));
  g1(8,12)=1;
  g1(8,19)=(-1);
  g1(9,3)=(-params(11));
  g1(9,13)=1;
  g1(9,20)=(-1);
  g1(10,4)=(-params(12));
  g1(10,14)=1;
  g1(10,21)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,441);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,9261);
end
end
end
end
