function [residual, g1, g2, g3] = BC_19March2019_model2_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(8, 1);
T20 = y(2)^params(5);
T50 = params(6)*((y(6)-params(7)*y(1))/(y(13)-y(6)*params(7)))^params(9);
T65 = exp(y(11))*params(1)*0.5*(1-y(5)^2);
T67 = y(7)^params(5);
lhs =y(10)*y(2);
rhs =params(3)*params(1)*exp(y(11))*y(5)*T20-params(4);
residual(1)= lhs-rhs;
lhs =T20*y(5)*params(1)*exp(y(11))-y(10)*y(2);
rhs =params(2);
residual(2)= lhs-rhs;
lhs =y(6)+y(9);
rhs =y(8)-params(2)*(1-y(5));
residual(3)= lhs-rhs;
lhs =1;
rhs =T50*(1+y(14)-params(8))*exp(y(12));
residual(4)= lhs-rhs;
lhs =y(8);
rhs =T65*T67;
residual(5)= lhs-rhs;
lhs =y(7);
rhs =y(9)+y(2)*(1-params(8));
residual(6)= lhs-rhs;
lhs =y(11);
rhs =params(10)*y(3)+x(it_, 1);
residual(7)= lhs-rhs;
lhs =y(12);
rhs =params(11)*y(4)-x(it_, 2);
residual(8)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(8, 16);

  %
  % Jacobian matrix
  %

T97 = getPowerDeriv((y(6)-params(7)*y(1))/(y(13)-y(6)*params(7)),params(9),1);
T119 = getPowerDeriv(y(2),params(5),1);
  g1(1,5)=(-(params(3)*params(1)*exp(y(11))*T20));
  g1(1,2)=y(10)-params(3)*params(1)*exp(y(11))*y(5)*T119;
  g1(1,10)=y(2);
  g1(1,11)=(-(params(3)*params(1)*exp(y(11))*y(5)*T20));
  g1(2,5)=T20*params(1)*exp(y(11));
  g1(2,2)=y(5)*params(1)*exp(y(11))*T119-y(10);
  g1(2,10)=(-y(2));
  g1(2,11)=T20*y(5)*params(1)*exp(y(11));
  g1(3,5)=(-params(2));
  g1(3,6)=1;
  g1(3,8)=(-1);
  g1(3,9)=1;
  g1(4,1)=(-(exp(y(12))*(1+y(14)-params(8))*params(6)*(-params(7))/(y(13)-y(6)*params(7))*T97));
  g1(4,6)=(-(exp(y(12))*(1+y(14)-params(8))*params(6)*T97*(y(13)-y(6)*params(7)-(y(6)-params(7)*y(1))*(-params(7)))/((y(13)-y(6)*params(7))*(y(13)-y(6)*params(7)))));
  g1(4,13)=(-(exp(y(12))*(1+y(14)-params(8))*params(6)*T97*(-(y(6)-params(7)*y(1)))/((y(13)-y(6)*params(7))*(y(13)-y(6)*params(7)))));
  g1(4,14)=(-(T50*exp(y(12))));
  g1(4,12)=(-(T50*(1+y(14)-params(8))*exp(y(12))));
  g1(5,5)=(-(T67*exp(y(11))*params(1)*0.5*(-(2*y(5)))));
  g1(5,7)=(-(T65*getPowerDeriv(y(7),params(5),1)));
  g1(5,8)=1;
  g1(5,11)=(-(T65*T67));
  g1(6,2)=(-(1-params(8)));
  g1(6,7)=1;
  g1(6,9)=(-1);
  g1(7,3)=(-params(10));
  g1(7,11)=1;
  g1(7,15)=(-1);
  g1(8,4)=(-params(11));
  g1(8,12)=1;
  g1(8,16)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,256);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,4096);
end
end
end
end
