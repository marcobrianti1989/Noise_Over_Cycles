function [residual, g1, g2, g3] = JQ_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(13, 1);
T42 = params(3)*exp(y(17))*y(7)^(params(3)-1);
T45 = y(2)^(1-params(3));
T63 = y(7)^params(3);
T66 = y(2)^(-params(3));
lhs =y(5)/y(6);
rhs =params(1)/(1-y(7));
residual(1)= lhs-rhs;
lhs =y(18)/y(6);
rhs =params(5)*(y(8)-params(2))/(1-params(2));
residual(2)= lhs-rhs;
residual(3) = y(5)*y(7)+y(1)-y(9)/y(8)+y(10)-y(6);
lhs =T42*T45;
rhs =y(5)*1/(1-y(12)*y(13));
residual(4)= lhs-rhs;
lhs =y(14)*(1-params(4)+exp(y(17))*(1-params(3))*(1-y(19)*y(20))*T63*T66)+y(13)*y(12)*exp(y(16));
rhs =1;
residual(5)= lhs-rhs;
lhs =y(12)*exp(y(16))+y(8)*y(14)+y(13)*y(8)*(1-params(2))/(y(8)-params(2));
rhs =1;
residual(6)= lhs-rhs;
residual(7) = y(9)/y(8)+y(2)*(1-params(4))+T45*exp(y(17))*T63-y(5)*y(7)-y(1)-y(11)-y(15);
lhs =exp(y(16))*(y(11)-(1-params(2))*y(9)/(y(8)-params(2)));
rhs =T45*exp(y(17))*T63;
residual(8)= lhs-rhs;
lhs =y(14);
rhs =params(5)*y(6)/y(18)*y(13)/y(20);
residual(9)= lhs-rhs;
lhs =y(15);
rhs =y(10)+params(8)*(y(10)-params(6))^2;
residual(10)= lhs-rhs;
lhs =y(13);
rhs =1+(y(10)-params(6))*2*params(8);
residual(11)= lhs-rhs;
lhs =exp(y(17));
rhs =exp(params(9)*y(4)+x(it_, 1));
residual(12)= lhs-rhs;
lhs =exp(y(16));
rhs =exp(params(10)*y(3)+params(7)*(1-params(10))+x(it_, 2));
residual(13)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(13, 22);

  %
  % Jacobian matrix
  %

T159 = getPowerDeriv(y(7),params(3),1);
T193 = getPowerDeriv(y(2),1-params(3),1);
  g1(1,5)=1/y(6);
  g1(1,6)=(-y(5))/(y(6)*y(6));
  g1(1,7)=(-(params(1)/((1-y(7))*(1-y(7)))));
  g1(2,6)=(-y(18))/(y(6)*y(6));
  g1(2,18)=1/y(6);
  g1(2,8)=(-(params(5)/(1-params(2))));
  g1(3,5)=y(7);
  g1(3,6)=(-1);
  g1(3,7)=y(5);
  g1(3,8)=(-((-y(9))/(y(8)*y(8))));
  g1(3,1)=1;
  g1(3,9)=(-(1/y(8)));
  g1(3,10)=1;
  g1(4,5)=(-(1/(1-y(12)*y(13))));
  g1(4,7)=T45*params(3)*exp(y(17))*getPowerDeriv(y(7),params(3)-1,1);
  g1(4,2)=T42*T193;
  g1(4,12)=(-(y(5)*y(13)/((1-y(12)*y(13))*(1-y(12)*y(13)))));
  g1(4,13)=(-(y(5)*y(12)/((1-y(12)*y(13))*(1-y(12)*y(13)))));
  g1(4,17)=T42*T45;
  g1(5,7)=y(14)*T66*exp(y(17))*(1-params(3))*(1-y(19)*y(20))*T159;
  g1(5,2)=y(14)*exp(y(17))*(1-params(3))*(1-y(19)*y(20))*T63*getPowerDeriv(y(2),(-params(3)),1);
  g1(5,12)=y(13)*exp(y(16));
  g1(5,19)=y(14)*T66*T63*exp(y(17))*(1-params(3))*(-y(20));
  g1(5,13)=y(12)*exp(y(16));
  g1(5,20)=y(14)*T66*T63*exp(y(17))*(1-params(3))*(-y(19));
  g1(5,14)=1-params(4)+exp(y(17))*(1-params(3))*(1-y(19)*y(20))*T63*T66;
  g1(5,16)=y(13)*y(12)*exp(y(16));
  g1(5,17)=y(14)*exp(y(17))*(1-params(3))*(1-y(19)*y(20))*T63*T66;
  g1(6,8)=y(14)+y(13)*((y(8)-params(2))*(1-params(2))-y(8)*(1-params(2)))/((y(8)-params(2))*(y(8)-params(2)));
  g1(6,12)=exp(y(16));
  g1(6,13)=y(8)*(1-params(2))/(y(8)-params(2));
  g1(6,14)=y(8);
  g1(6,16)=y(12)*exp(y(16));
  g1(7,5)=(-y(7));
  g1(7,7)=T45*exp(y(17))*T159-y(5);
  g1(7,8)=(-y(9))/(y(8)*y(8));
  g1(7,1)=(-1);
  g1(7,9)=1/y(8);
  g1(7,2)=1-params(4)+exp(y(17))*T63*T193;
  g1(7,11)=(-1);
  g1(7,15)=(-1);
  g1(7,17)=T45*exp(y(17))*T63;
  g1(8,7)=(-(T45*exp(y(17))*T159));
  g1(8,8)=exp(y(16))*(-((-((1-params(2))*y(9)))/((y(8)-params(2))*(y(8)-params(2)))));
  g1(8,9)=exp(y(16))*(-((1-params(2))/(y(8)-params(2))));
  g1(8,2)=(-(exp(y(17))*T63*T193));
  g1(8,11)=exp(y(16));
  g1(8,16)=exp(y(16))*(y(11)-(1-params(2))*y(9)/(y(8)-params(2)));
  g1(8,17)=(-(T45*exp(y(17))*T63));
  g1(9,6)=(-(y(13)/y(20)*params(5)*1/y(18)));
  g1(9,18)=(-(y(13)/y(20)*params(5)*(-y(6))/(y(18)*y(18))));
  g1(9,13)=(-(params(5)*y(6)/y(18)*1/y(20)));
  g1(9,20)=(-(params(5)*y(6)/y(18)*(-y(13))/(y(20)*y(20))));
  g1(9,14)=1;
  g1(10,10)=(-(1+params(8)*2*(y(10)-params(6))));
  g1(10,15)=1;
  g1(11,10)=(-(2*params(8)));
  g1(11,13)=1;
  g1(12,4)=(-(params(9)*exp(params(9)*y(4)+x(it_, 1))));
  g1(12,17)=exp(y(17));
  g1(12,21)=(-exp(params(9)*y(4)+x(it_, 1)));
  g1(13,3)=(-(params(10)*exp(params(10)*y(3)+params(7)*(1-params(10))+x(it_, 2))));
  g1(13,16)=exp(y(16));
  g1(13,22)=(-exp(params(10)*y(3)+params(7)*(1-params(10))+x(it_, 2)));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],13,484);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],13,10648);
end
end
end
end
