function [residual, g1, g2, g3] = BGP_Jan2019_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(5, 1);
T28 = exp(y(8))*y(6)^params(8);
T57 = exp(y(3))*y(1)^params(8);
T60 = params(7)*(y(7)+T28)-params(7)*params(2)*(y(2)+T57);
T73 = (y(6)+(1-y(6))*params(5))*params(11)*params(10)*exp(y(9))/exp(y(4));
T79 = exp(y(12))*y(10)^params(8);
T83 = params(7)*(y(11)+T79)-(y(7)+T28)*params(7)*params(2);
T85 = T83^(1-params(1))-1;
lhs =y(5);
rhs =(1+(1-y(6))*params(5)*params(6))/(y(6)+(1-y(6))*params(5));
residual(1)= lhs-rhs;
lhs =y(11);
rhs =(1-params(9))*y(7)+T28;
residual(2)= lhs-rhs;
lhs =y(8);
rhs =params(12)*y(3)+params(13)*x(it_, 1);
residual(3)= lhs-rhs;
lhs =y(9);
rhs =params(14)*y(4)+params(15)*x(it_, 2);
residual(4)= lhs-rhs;
lhs =(T60^(1-params(1))-1)/(1-params(1));
rhs =y(5)*T73*T85/(1-params(1));
residual(5)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(5, 14);

  %
  % Jacobian matrix
  %

T96 = getPowerDeriv(T60,1-params(1),1);
T109 = exp(y(8))*getPowerDeriv(y(6),params(8),1);
T118 = getPowerDeriv(T83,1-params(1),1);
  g1(1,5)=1;
  g1(1,6)=(-(((y(6)+(1-y(6))*params(5))*params(6)*(-params(5))-(1+(1-y(6))*params(5)*params(6))*(1-params(5)))/((y(6)+(1-y(6))*params(5))*(y(6)+(1-y(6))*params(5)))));
  g1(2,6)=(-T109);
  g1(2,7)=(-(1-params(9)));
  g1(2,11)=1;
  g1(2,8)=(-T28);
  g1(3,3)=(-params(12));
  g1(3,8)=1;
  g1(3,13)=(-params(13));
  g1(4,4)=(-params(14));
  g1(4,9)=1;
  g1(4,14)=(-params(15));
  g1(5,5)=(-(T73*T85/(1-params(1))));
  g1(5,1)=(-(params(7)*params(2)*exp(y(3))*getPowerDeriv(y(1),params(8),1)))*T96/(1-params(1));
  g1(5,6)=T96*params(7)*T109/(1-params(1))-(T85*y(5)*params(11)*params(10)*exp(y(9))/exp(y(4))*(1-params(5))+y(5)*T73*(-(params(7)*params(2)*T109))*T118)/(1-params(1));
  g1(5,10)=(-(y(5)*T73*T118*params(7)*exp(y(12))*getPowerDeriv(y(10),params(8),1)/(1-params(1))));
  g1(5,2)=T96*(-(params(7)*params(2)))/(1-params(1));
  g1(5,7)=params(7)*T96/(1-params(1))-y(5)*T73*T118*(-(params(7)*params(2)))/(1-params(1));
  g1(5,11)=(-(y(5)*T73*params(7)*T118/(1-params(1))));
  g1(5,3)=T96*(-(params(7)*params(2)*T57))/(1-params(1));
  g1(5,8)=T96*T28*params(7)/(1-params(1))-y(5)*T73*T118*(-(T28*params(7)*params(2)))/(1-params(1));
  g1(5,12)=(-(y(5)*T73*T118*params(7)*T79/(1-params(1))));
  g1(5,4)=(-(T85*y(5)*(y(6)+(1-y(6))*params(5))*(-(params(11)*params(10)*exp(y(9))*exp(y(4))))/(exp(y(4))*exp(y(4)))/(1-params(1))));
  g1(5,9)=(-(y(5)*T73*T85/(1-params(1))));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],5,196);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],5,2744);
end
end
end
end
