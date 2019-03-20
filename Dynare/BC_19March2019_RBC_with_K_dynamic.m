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
T16 = exp(y(11))*params(2)*y(7)^(params(2)-1);
T34 = params(3)*((y(6)-params(4)*y(1))/(y(13)-y(6)*params(4)))^params(6);
T49 = exp(y(11))*0.5*params(1)*y(7)^params(2);
lhs =y(10);
rhs =T16;
residual(1)= lhs-rhs;
lhs =y(6)+y(9);
rhs =y(8);
residual(2)= lhs-rhs;
lhs =1;
rhs =T34*(1+y(14)-params(5))*exp(y(12));
residual(3)= lhs-rhs;
lhs =y(8);
rhs =T49;
residual(4)= lhs-rhs;
lhs =y(7);
rhs =(1-params(5))*y(2)+y(3);
residual(5)= lhs-rhs;
lhs =y(11);
rhs =params(8)*y(4)+x(it_, 1);
residual(6)= lhs-rhs;
lhs =y(12);
rhs =params(9)*y(5)-x(it_, 2);
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 16);

  %
  % Jacobian matrix
  %

T71 = getPowerDeriv((y(6)-params(4)*y(1))/(y(13)-y(6)*params(4)),params(6),1);
  g1(1,7)=(-(exp(y(11))*params(2)*getPowerDeriv(y(7),params(2)-1,1)));
  g1(1,10)=1;
  g1(1,11)=(-T16);
  g1(2,6)=1;
  g1(2,8)=(-1);
  g1(2,9)=1;
  g1(3,1)=(-(exp(y(12))*(1+y(14)-params(5))*params(3)*(-params(4))/(y(13)-y(6)*params(4))*T71));
  g1(3,6)=(-(exp(y(12))*(1+y(14)-params(5))*params(3)*T71*(y(13)-y(6)*params(4)-(y(6)-params(4)*y(1))*(-params(4)))/((y(13)-y(6)*params(4))*(y(13)-y(6)*params(4)))));
  g1(3,13)=(-(exp(y(12))*(1+y(14)-params(5))*params(3)*T71*(-(y(6)-params(4)*y(1)))/((y(13)-y(6)*params(4))*(y(13)-y(6)*params(4)))));
  g1(3,14)=(-(T34*exp(y(12))));
  g1(3,12)=(-(T34*(1+y(14)-params(5))*exp(y(12))));
  g1(4,7)=(-(exp(y(11))*0.5*params(1)*getPowerDeriv(y(7),params(2),1)));
  g1(4,8)=1;
  g1(4,11)=(-T49);
  g1(5,2)=(-(1-params(5)));
  g1(5,7)=1;
  g1(5,3)=(-1);
  g1(6,4)=(-params(8));
  g1(6,11)=1;
  g1(6,15)=(-1);
  g1(7,5)=(-params(9));
  g1(7,12)=1;
  g1(7,16)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,256);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,4096);
end
end
end
end
