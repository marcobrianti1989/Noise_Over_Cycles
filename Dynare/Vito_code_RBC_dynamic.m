function [residual, g1, g2, g3] = Vito_code_RBC_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(9, 1);
lhs =y(8);
rhs =(1-params(3))*y(2)+params(3)*y(7);
residual(1)= lhs-rhs;
lhs =y(11);
rhs =params(8)*y(4)+params(7)*y(5)+y(7)*(1-params(8)-params(7));
residual(2)= lhs-rhs;
lhs =y(11);
rhs =y(2)*params(1)+(1-params(1))*y(12)+(1-params(1))*y(6);
residual(3)= lhs-rhs;
lhs =y(10);
rhs =y(2)*params(1)+(1-params(1))*y(12)-params(1)*y(6);
residual(4)= lhs-rhs;
lhs =y(9);
rhs =(1-params(1))*y(6)+(1-params(1))*y(12)+y(2)*(-(1-params(1)));
residual(5)= lhs-rhs;
lhs =y(4)+y(6)*params(4);
rhs =y(10);
residual(6)= lhs-rhs;
lhs =y(13)-y(4);
rhs =(1-params(2)+params(3)*params(2))*y(14);
residual(7)= lhs-rhs;
lhs =y(12);
rhs =params(6)*y(3)+x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(5);
rhs =params(5)*y(1)+x(it_, 2);
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 16);

  %
  % Jacobian matrix
  %

  g1(1,7)=(-params(3));
  g1(1,2)=(-(1-params(3)));
  g1(1,8)=1;
  g1(2,4)=(-params(8));
  g1(2,5)=(-params(7));
  g1(2,7)=(-(1-params(8)-params(7)));
  g1(2,11)=1;
  g1(3,6)=(-(1-params(1)));
  g1(3,2)=(-params(1));
  g1(3,11)=1;
  g1(3,12)=(-(1-params(1)));
  g1(4,6)=params(1);
  g1(4,2)=(-params(1));
  g1(4,10)=1;
  g1(4,12)=(-(1-params(1)));
  g1(5,6)=(-(1-params(1)));
  g1(5,2)=1-params(1);
  g1(5,9)=1;
  g1(5,12)=(-(1-params(1)));
  g1(6,4)=1;
  g1(6,6)=params(4);
  g1(6,10)=(-1);
  g1(7,4)=(-1);
  g1(7,13)=1;
  g1(7,14)=(-(1-params(2)+params(3)*params(2)));
  g1(8,3)=(-params(6));
  g1(8,12)=1;
  g1(8,15)=(-1);
  g1(9,1)=(-params(5));
  g1(9,5)=1;
  g1(9,16)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,256);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,4096);
end
end
end
end
