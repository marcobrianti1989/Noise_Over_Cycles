function [residual, g1, g2, g3] = Beaudry_Portier_2004_habit_consumption_dynamic(y, x, params, steady_state, it_)
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
T19 = (1-y(13)+params(11)*y(3))^params(9);
T20 = exp(y(17))*params(12)/T19;
T28 = exp(y(16))*y(10)^params(17)*params(3);
T32 = y(10)^(params(8)-1);
T34 = (T28*params(8)*T32)^(-1);
T47 = (1-y(21)+y(13)*params(11))^params(9);
T48 = params(12)*(1-params(1))*exp(y(23))/T47;
T55 = params(8)*params(3)*exp(y(22))*y(20)^params(17);
T56 = y(20)^(params(8)-1);
T58 = (T55*T56)^(-1);
T63 = y(18)^(1-params(6));
T74 = y(19)^(params(6)-1);
T83 = y(8)^(1-params(6));
T85 = params(5)*1/(y(8)-params(10)*y(1))*T83;
T90 = y(11)^params(16);
T95 = y(11)^params(7);
T96 = exp(y(15))*T90*params(4)*T95;
T97 = T96^(params(6)-1);
T99 = exp(y(15))*T85*T97;
T102 = params(7)*params(4)*T90*T99;
T104 = y(11)^(params(7)-1);
T113 = params(5)*y(12)^params(6)+(1-params(5))*y(9)^params(6);
T125 = y(10)^params(8);
lhs =T20*T34;
rhs =params(2)*(T48*T58+T63/(y(18)-params(10)*y(8))*(1-params(5))*T74);
residual(1)= lhs-rhs;
lhs =T20;
rhs =T102*T104;
residual(2)= lhs-rhs;
lhs =y(8);
rhs =T113^(1/params(6));
residual(3)= lhs-rhs;
lhs =y(13);
rhs =y(10)+y(11);
residual(4)= lhs-rhs;
lhs =y(9);
rhs =(1-params(1))*y(2)+y(4);
residual(5)= lhs-rhs;
lhs =y(14);
rhs =T28*T125;
residual(6)= lhs-rhs;
lhs =y(12);
rhs =T96;
residual(7)= lhs-rhs;
lhs =y(15);
rhs =params(13)*y(5)+x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(16);
rhs =params(14)*y(6)+x(it_, 2);
residual(9)= lhs-rhs;
lhs =y(17);
rhs =params(15)*y(7)-x(it_, 3);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 26);

  %
  % Jacobian matrix
  %

T191 = getPowerDeriv(T113,1/params(6),1);
T200 = params(3)*exp(y(16))*getPowerDeriv(y(10),params(17),1);
T206 = getPowerDeriv(T28*params(8)*T32,(-1),1);
T222 = getPowerDeriv(T55*T56,(-1),1);
T227 = getPowerDeriv(y(11),params(16),1);
T233 = T95*params(4)*exp(y(15))*T227+exp(y(15))*T90*params(4)*getPowerDeriv(y(11),params(7),1);
T234 = getPowerDeriv(T96,params(6)-1,1);
T253 = getPowerDeriv(1-y(13)+params(11)*y(3),params(9),1);
T265 = getPowerDeriv(1-y(21)+y(13)*params(11),params(9),1);
  g1(1,8)=(-(params(2)*T74*(1-params(5))*(-(T63*(-params(10))))/((y(18)-params(10)*y(8))*(y(18)-params(10)*y(8)))));
  g1(1,18)=(-(params(2)*T74*(1-params(5))*((y(18)-params(10)*y(8))*getPowerDeriv(y(18),1-params(6),1)-T63)/((y(18)-params(10)*y(8))*(y(18)-params(10)*y(8)))));
  g1(1,19)=(-(params(2)*T63/(y(18)-params(10)*y(8))*(1-params(5))*getPowerDeriv(y(19),params(6)-1,1)));
  g1(1,10)=T20*(T32*params(8)*T200+T28*params(8)*getPowerDeriv(y(10),params(8)-1,1))*T206;
  g1(1,20)=(-(params(2)*T48*(T56*params(8)*params(3)*exp(y(22))*getPowerDeriv(y(20),params(17),1)+T55*getPowerDeriv(y(20),params(8)-1,1))*T222));
  g1(1,3)=T34*(-(exp(y(17))*params(12)*params(11)*T253))/(T19*T19);
  g1(1,13)=T34*(-(exp(y(17))*params(12)*(-T253)))/(T19*T19)-params(2)*T58*(-(params(12)*(1-params(1))*exp(y(23))*params(11)*T265))/(T47*T47);
  g1(1,21)=(-(params(2)*T58*(-(params(12)*(1-params(1))*exp(y(23))*(-T265)))/(T47*T47)));
  g1(1,16)=T20*T28*params(8)*T32*T206;
  g1(1,22)=(-(params(2)*T48*T55*T56*T222));
  g1(1,17)=T20*T34;
  g1(1,23)=(-(params(2)*T48*T58));
  g1(2,1)=(-(T104*params(7)*params(4)*T90*exp(y(15))*T97*params(5)*T83*params(10)/((y(8)-params(10)*y(1))*(y(8)-params(10)*y(1)))));
  g1(2,8)=(-(T104*params(7)*params(4)*T90*exp(y(15))*T97*params(5)*(T83*(-1)/((y(8)-params(10)*y(1))*(y(8)-params(10)*y(1)))+1/(y(8)-params(10)*y(1))*getPowerDeriv(y(8),1-params(6),1))));
  g1(2,11)=(-(T104*params(7)*params(4)*(T99*T227+T90*exp(y(15))*T85*T233*T234)+T102*getPowerDeriv(y(11),params(7)-1,1)));
  g1(2,3)=(-(exp(y(17))*params(12)*params(11)*T253))/(T19*T19);
  g1(2,13)=(-(exp(y(17))*params(12)*(-T253)))/(T19*T19);
  g1(2,15)=(-(T104*params(7)*params(4)*T90*(T99+exp(y(15))*T85*T96*T234)));
  g1(2,17)=T20;
  g1(3,8)=1;
  g1(3,9)=(-((1-params(5))*getPowerDeriv(y(9),params(6),1)*T191));
  g1(3,12)=(-(T191*params(5)*getPowerDeriv(y(12),params(6),1)));
  g1(4,10)=(-1);
  g1(4,11)=(-1);
  g1(4,13)=1;
  g1(5,2)=(-(1-params(1)));
  g1(5,9)=1;
  g1(5,4)=(-1);
  g1(6,10)=(-(T125*T200+T28*getPowerDeriv(y(10),params(8),1)));
  g1(6,14)=1;
  g1(6,16)=(-(T28*T125));
  g1(7,11)=(-T233);
  g1(7,12)=1;
  g1(7,15)=(-T96);
  g1(8,5)=(-params(13));
  g1(8,15)=1;
  g1(8,24)=(-1);
  g1(9,6)=(-params(14));
  g1(9,16)=1;
  g1(9,25)=(-1);
  g1(10,7)=(-params(15));
  g1(10,17)=1;
  g1(10,26)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,676);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,17576);
end
end
end
end
