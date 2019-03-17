function [residual, g1, g2, g3] = RBC_CME_basic_dynamic(y, x, params, steady_state, it_)
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
Rbar__ = 1/params(1);
T13 = params(1)/y(10);
T19 = y(5)^(params(3)-1);
T23 = params(3)*y(11)*T19+1-params(2);
T33 = y(1)^params(3);
T34 = y(6)*T33;
lhs =1/y(4);
rhs =T13*T23;
residual(1)= lhs-rhs;
lhs =1/y(4);
rhs =T13*y(8)/y(12);
residual(2)= lhs-rhs;
lhs =T34;
rhs =y(4)+y(5)-(1-params(2))*y(1);
residual(3)= lhs-rhs;
lhs =y(3);
rhs =T34;
residual(4)= lhs-rhs;
lhs =y(8)/Rbar__;
rhs =(y(9)/params(6))^params(5);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =exp(y(7));
residual(6)= lhs-rhs;
lhs =y(7);
rhs =params(4)*y(2)+x(it_, 1);
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 13);

  %
  % Jacobian matrix
  %

T68 = y(6)*getPowerDeriv(y(1),params(3),1);
  g1(1,4)=(-1)/(y(4)*y(4));
  g1(1,10)=(-(T23*(-params(1))/(y(10)*y(10))));
  g1(1,5)=(-(T13*params(3)*y(11)*getPowerDeriv(y(5),params(3)-1,1)));
  g1(1,11)=(-(T13*params(3)*T19));
  g1(2,4)=(-1)/(y(4)*y(4));
  g1(2,10)=(-(y(8)/y(12)*(-params(1))/(y(10)*y(10))));
  g1(2,8)=(-(T13*1/y(12)));
  g1(2,12)=(-(T13*(-y(8))/(y(12)*y(12))));
  g1(3,4)=(-1);
  g1(3,1)=T68-(-(1-params(2)));
  g1(3,5)=(-1);
  g1(3,6)=T33;
  g1(4,3)=1;
  g1(4,1)=(-T68);
  g1(4,6)=(-T33);
  g1(5,8)=1/Rbar__;
  g1(5,9)=(-(1/params(6)*getPowerDeriv(y(9)/params(6),params(5),1)));
  g1(6,6)=1;
  g1(6,7)=(-exp(y(7)));
  g1(7,2)=(-params(4));
  g1(7,7)=1;
  g1(7,13)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,169);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,2197);
end
end
end
end
