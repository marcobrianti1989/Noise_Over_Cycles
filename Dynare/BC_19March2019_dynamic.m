function [residual, g1, g2, g3] = BC_19March2019_dynamic(y, x, params, steady_state, it_)
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
T12 = params(7)*(1-y(5)^2);
T16 = y(9)^(1+params(8));
T17 = T12*T16;
T26 = y(9)^params(5);
T42 = params(1)*0.5*(1-y(1)^2);
T44 = y(4)^params(5);
T45 = T42*T44;
lhs =T17;
rhs =y(5)*params(3)*params(1)*exp(y(7))*T26-params(4)*exp(y(8));
residual(1)= lhs-rhs;
residual(2) = y(5)*T26-T17-params(2)/T45;
lhs =y(6);
rhs =T26*(1-y(5));
residual(3)= lhs-rhs;
lhs =y(7);
rhs =params(9)*y(2)+x(it_, 1);
residual(4)= lhs-rhs;
lhs =y(8);
rhs =params(10)*y(3)-x(it_, 2);
residual(5)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(5, 11);

  %
  % Jacobian matrix
  %

T93 = T12*getPowerDeriv(y(9),1+params(8),1);
T94 = getPowerDeriv(y(9),params(5),1);
  g1(1,5)=T16*params(7)*(-(2*y(5)))-params(3)*params(1)*exp(y(7))*T26;
  g1(1,7)=(-(y(5)*params(3)*params(1)*exp(y(7))*T26));
  g1(1,8)=params(4)*exp(y(8));
  g1(1,9)=T93-y(5)*params(3)*params(1)*exp(y(7))*T94;
  g1(2,1)=(-((-(params(2)*T44*params(1)*0.5*(-(2*y(1)))))/(T45*T45)));
  g1(2,5)=T26-T16*params(7)*(-(2*y(5)));
  g1(2,4)=(-((-(params(2)*T42*getPowerDeriv(y(4),params(5),1)))/(T45*T45)));
  g1(2,9)=y(5)*T94-T93;
  g1(3,5)=T26;
  g1(3,6)=1;
  g1(3,9)=(-((1-y(5))*T94));
  g1(4,2)=(-params(9));
  g1(4,7)=1;
  g1(4,10)=(-1);
  g1(5,3)=(-params(10));
  g1(5,8)=1;
  g1(5,11)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],5,121);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],5,1331);
end
end
end
end
