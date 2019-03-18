function [residual, g1, g2, g3] = BC_19March2019_model3_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(11, 1);
T16 = params(1)*exp(y(15))*y(9)^params(5);
T19 = y(10)^(1-params(5));
T20 = T16*T19;
T68 = params(6)*((y(8)-params(7)*y(2))/(y(17)-y(8)*params(7)))^params(9);
lhs =y(7);
rhs =T20*params(2)-params(4);
residual(1)= lhs-rhs;
lhs =y(6);
rhs =params(14)*y(1)-params(3)*log(y(3));
residual(2)= lhs-rhs;
lhs =y(7);
rhs =y(10)*y(11)+y(9)*y(14);
residual(3)= lhs-rhs;
lhs =params(10)*y(10)^params(11);
rhs =y(11);
residual(4)= lhs-rhs;
lhs =y(8)+y(13);
rhs =y(12);
residual(5)= lhs-rhs;
lhs =y(11)*params(5)*y(10);
rhs =y(9)*(1-params(5))*y(14);
residual(6)= lhs-rhs;
lhs =1;
rhs =T68*(1+y(19)-params(8))*exp(y(16));
residual(7)= lhs-rhs;
lhs =y(12);
rhs =T20;
residual(8)= lhs-rhs;
lhs =y(18);
rhs =y(13)+y(9)*(1-params(8));
residual(9)= lhs-rhs;
lhs =y(15);
rhs =params(12)*y(4)+x(it_, 1);
residual(10)= lhs-rhs;
lhs =y(16);
rhs =params(13)*y(5)-x(it_, 2);
residual(11)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(11, 21);

  %
  % Jacobian matrix
  %

T99 = getPowerDeriv((y(8)-params(7)*y(2))/(y(17)-y(8)*params(7)),params(9),1);
T126 = T19*params(1)*exp(y(15))*getPowerDeriv(y(9),params(5),1);
T134 = T16*getPowerDeriv(y(10),1-params(5),1);
  g1(1,7)=1;
  g1(1,9)=(-(params(2)*T126));
  g1(1,10)=(-(params(2)*T134));
  g1(1,15)=(-(T20*params(2)));
  g1(2,1)=(-params(14));
  g1(2,6)=1;
  g1(2,3)=params(3)*1/y(3);
  g1(3,7)=1;
  g1(3,9)=(-y(14));
  g1(3,10)=(-y(11));
  g1(3,11)=(-y(10));
  g1(3,14)=(-y(9));
  g1(4,10)=params(10)*getPowerDeriv(y(10),params(11),1);
  g1(4,11)=(-1);
  g1(5,8)=1;
  g1(5,12)=(-1);
  g1(5,13)=1;
  g1(6,9)=(-((1-params(5))*y(14)));
  g1(6,10)=params(5)*y(11);
  g1(6,11)=params(5)*y(10);
  g1(6,14)=(-(y(9)*(1-params(5))));
  g1(7,2)=(-(exp(y(16))*(1+y(19)-params(8))*params(6)*(-params(7))/(y(17)-y(8)*params(7))*T99));
  g1(7,8)=(-(exp(y(16))*(1+y(19)-params(8))*params(6)*T99*(y(17)-y(8)*params(7)-(y(8)-params(7)*y(2))*(-params(7)))/((y(17)-y(8)*params(7))*(y(17)-y(8)*params(7)))));
  g1(7,17)=(-(exp(y(16))*(1+y(19)-params(8))*params(6)*T99*(-(y(8)-params(7)*y(2)))/((y(17)-y(8)*params(7))*(y(17)-y(8)*params(7)))));
  g1(7,19)=(-(T68*exp(y(16))));
  g1(7,16)=(-(T68*(1+y(19)-params(8))*exp(y(16))));
  g1(8,9)=(-T126);
  g1(8,10)=(-T134);
  g1(8,12)=1;
  g1(8,15)=(-T20);
  g1(9,9)=(-(1-params(8)));
  g1(9,18)=1;
  g1(9,13)=(-1);
  g1(10,4)=(-params(12));
  g1(10,15)=1;
  g1(10,20)=(-1);
  g1(11,5)=(-params(13));
  g1(11,16)=1;
  g1(11,21)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],11,441);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],11,9261);
end
end
end
end
