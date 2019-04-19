function [residual, g1, g2, g3] = BP2004_extension_dynamic(y, x, params, steady_state, it_)
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
T16 = y(19)^(1-params(6));
T27 = y(7)^(params(6)-1);
T34 = params(10)/2;
T36 = y(20)^2;
T37 = y(7)^2;
T42 = 1+T34*(T36/T37-params(1)^2);
T45 = exp(y(18))*params(2)*(T16/(y(19)-params(9)*y(6))*(1-params(5))*T27+(1-params(1))*y(21)*T42);
T48 = y(6)^(1-params(6));
T55 = y(10)^(params(6)-1);
T64 = y(9)^params(15);
T69 = y(9)^params(7);
T70 = exp(y(16))*T64*params(4)*T69;
T76 = params(5)*T70^params(6)+(1-params(5))*y(2)^params(6);
T90 = (y(12)/y(2)-params(1))^2;
T97 = y(8)^params(16);
T102 = y(8)^params(8);
T110 = y(8)^(params(8)-1);
T124 = y(9)^(params(7)-1);
lhs =y(13);
rhs =T45;
residual(1)= lhs-rhs;
lhs =params(11);
rhs =params(5)*T48/(y(6)-params(9)*y(1))*T55*y(15);
residual(2)= lhs-rhs;
lhs =y(6);
rhs =T76^(1/params(6));
residual(3)= lhs-rhs;
lhs =y(11);
rhs =y(9)+y(8);
residual(4)= lhs-rhs;
lhs =y(7);
rhs =(1-params(1))*y(2)+y(12)-T34*y(2)*T90;
residual(5)= lhs-rhs;
lhs =y(12);
rhs =exp(y(17))*T97*params(3)*T102;
residual(6)= lhs-rhs;
lhs =y(14);
rhs =params(3)*T97*exp(y(17))*params(8)*T110;
residual(7)= lhs-rhs;
lhs =y(13);
rhs =params(11)*(y(14)*(1-params(10)*(y(12)/y(2)-params(1))))^(-1);
residual(8)= lhs-rhs;
lhs =y(10);
rhs =T70;
residual(9)= lhs-rhs;
lhs =y(15);
rhs =params(4)*T64*exp(y(16))*params(7)*T124;
residual(10)= lhs-rhs;
lhs =y(16);
rhs =params(12)*y(3)+x(it_, 1);
residual(11)= lhs-rhs;
lhs =y(17);
rhs =params(13)*y(4)+x(it_, 2);
residual(12)= lhs-rhs;
lhs =y(18);
rhs =params(14)*y(5)+x(it_, 3);
residual(13)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(13, 24);

  %
  % Jacobian matrix
  %

T180 = getPowerDeriv(T76,1/params(6),1);
T196 = getPowerDeriv(y(14)*(1-params(10)*(y(12)/y(2)-params(1))),(-1),1);
T212 = getPowerDeriv(y(8),params(16),1);
T227 = getPowerDeriv(y(9),params(15),1);
T233 = T69*params(4)*exp(y(16))*T227+exp(y(16))*T64*params(4)*getPowerDeriv(y(9),params(7),1);
T234 = getPowerDeriv(T70,params(6),1);
  g1(1,6)=(-(exp(y(18))*params(2)*T27*(1-params(5))*(-(T16*(-params(9))))/((y(19)-params(9)*y(6))*(y(19)-params(9)*y(6)))));
  g1(1,19)=(-(exp(y(18))*params(2)*T27*(1-params(5))*((y(19)-params(9)*y(6))*getPowerDeriv(y(19),1-params(6),1)-T16)/((y(19)-params(9)*y(6))*(y(19)-params(9)*y(6)))));
  g1(1,7)=(-(exp(y(18))*params(2)*(T16/(y(19)-params(9)*y(6))*(1-params(5))*getPowerDeriv(y(7),params(6)-1,1)+(1-params(1))*y(21)*T34*(-(T36*2*y(7)))/(T37*T37))));
  g1(1,20)=(-(exp(y(18))*params(2)*(1-params(1))*y(21)*T34*2*y(20)/T37));
  g1(1,13)=1;
  g1(1,21)=(-(exp(y(18))*params(2)*(1-params(1))*T42));
  g1(1,18)=(-T45);
  g1(2,1)=(-(y(15)*T55*params(5)*(-(T48*(-params(9))))/((y(6)-params(9)*y(1))*(y(6)-params(9)*y(1)))));
  g1(2,6)=(-(y(15)*T55*params(5)*((y(6)-params(9)*y(1))*getPowerDeriv(y(6),1-params(6),1)-T48)/((y(6)-params(9)*y(1))*(y(6)-params(9)*y(1)))));
  g1(2,10)=(-(y(15)*params(5)*T48/(y(6)-params(9)*y(1))*getPowerDeriv(y(10),params(6)-1,1)));
  g1(2,15)=(-(params(5)*T48/(y(6)-params(9)*y(1))*T55));
  g1(3,6)=1;
  g1(3,2)=(-((1-params(5))*getPowerDeriv(y(2),params(6),1)*T180));
  g1(3,9)=(-(T180*params(5)*T233*T234));
  g1(3,16)=(-(T180*params(5)*T70*T234));
  g1(4,8)=(-1);
  g1(4,9)=(-1);
  g1(4,11)=1;
  g1(5,2)=(-(1-params(1)-(T34*T90+T34*y(2)*(-y(12))/(y(2)*y(2))*2*(y(12)/y(2)-params(1)))));
  g1(5,7)=1;
  g1(5,12)=(-(1-T34*y(2)*2*(y(12)/y(2)-params(1))*1/y(2)));
  g1(6,8)=(-(T102*params(3)*exp(y(17))*T212+exp(y(17))*T97*params(3)*getPowerDeriv(y(8),params(8),1)));
  g1(6,12)=1;
  g1(6,17)=(-(exp(y(17))*T97*params(3)*T102));
  g1(7,8)=(-(T110*params(3)*exp(y(17))*params(8)*T212+params(3)*T97*exp(y(17))*params(8)*getPowerDeriv(y(8),params(8)-1,1)));
  g1(7,14)=1;
  g1(7,17)=(-(params(3)*T97*exp(y(17))*params(8)*T110));
  g1(8,2)=(-(params(11)*y(14)*(-(params(10)*(-y(12))/(y(2)*y(2))))*T196));
  g1(8,12)=(-(params(11)*T196*y(14)*(-(params(10)*1/y(2)))));
  g1(8,13)=1;
  g1(8,14)=(-(params(11)*(1-params(10)*(y(12)/y(2)-params(1)))*T196));
  g1(9,9)=(-T233);
  g1(9,10)=1;
  g1(9,16)=(-T70);
  g1(10,9)=(-(T124*params(4)*exp(y(16))*params(7)*T227+params(4)*T64*exp(y(16))*params(7)*getPowerDeriv(y(9),params(7)-1,1)));
  g1(10,15)=1;
  g1(10,16)=(-(params(4)*T64*exp(y(16))*params(7)*T124));
  g1(11,3)=(-params(12));
  g1(11,16)=1;
  g1(11,22)=(-1);
  g1(12,4)=(-params(13));
  g1(12,17)=1;
  g1(12,23)=(-1);
  g1(13,5)=(-params(14));
  g1(13,18)=1;
  g1(13,24)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],13,576);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],13,13824);
end
end
end
end
