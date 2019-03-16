function [residual, g1, g2, g3] = BC_19March2019_RBC_with_K_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(7, 1);
T13 = params(2)*y(6)^(params(2)-1);
T31 = params(3)*((y(5)-params(4)*y(1))/(y(12)-y(5)*params(4)))^params(6);
T48 = 0.5*params(1)*exp(y(10))*y(2)^params(2);
lhs =y(9);
rhs =T13;
residual(1)= lhs-rhs;
lhs =y(5)+y(8);
rhs =y(7);
residual(2)= lhs-rhs;
lhs =1;
rhs =T31*(1+T13-params(5))*exp(y(11));
residual(3)= lhs-rhs;
lhs =y(7);
rhs =T48;
residual(4)= lhs-rhs;
lhs =y(6);
rhs =y(8)+y(2)*(1-params(5));
residual(5)= lhs-rhs;
lhs =y(10);
rhs =params(7)*y(3)+x(it_, 1);
residual(6)= lhs-rhs;
lhs =y(11);
rhs =params(8)*y(4)-x(it_, 2);
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 14);

  %
  % Jacobian matrix
  %

T68 = getPowerDeriv((y(5)-params(4)*y(1))/(y(12)-y(5)*params(4)),params(6),1);
T95 = params(2)*getPowerDeriv(y(6),params(2)-1,1);
  g1(1,6)=(-T95);
  g1(1,9)=1;
  g1(2,5)=1;
  g1(2,7)=(-1);
  g1(2,8)=1;
  g1(3,1)=(-(exp(y(11))*(1+T13-params(5))*params(3)*(-params(4))/(y(12)-y(5)*params(4))*T68));
  g1(3,5)=(-(exp(y(11))*(1+T13-params(5))*params(3)*T68*(y(12)-y(5)*params(4)-(y(5)-params(4)*y(1))*(-params(4)))/((y(12)-y(5)*params(4))*(y(12)-y(5)*params(4)))));
  g1(3,12)=(-(exp(y(11))*(1+T13-params(5))*params(3)*T68*(-(y(5)-params(4)*y(1)))/((y(12)-y(5)*params(4))*(y(12)-y(5)*params(4)))));
  g1(3,6)=(-(exp(y(11))*T31*T95));
  g1(3,11)=(-(T31*(1+T13-params(5))*exp(y(11))));
  g1(4,2)=(-(0.5*params(1)*exp(y(10))*getPowerDeriv(y(2),params(2),1)));
  g1(4,7)=1;
  g1(4,10)=(-T48);
  g1(5,2)=(-(1-params(5)));
  g1(5,6)=1;
  g1(5,8)=(-1);
  g1(6,3)=(-params(7));
  g1(6,10)=1;
  g1(6,13)=(-1);
  g1(7,4)=(-params(8));
  g1(7,11)=1;
  g1(7,14)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,196);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,2744);
end
end
end
end
