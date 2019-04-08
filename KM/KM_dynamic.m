function [residual, g1, g2, g3] = KM_dynamic(y, x, params, steady_state, it_)
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
T17 = params(1)/(y(7)+params(3)-y(17)/params(4));
T63 = params(7)*exp(y(18))*((y(10)-y(6))/params(6))^(params(7)-1)/params(4);
T74 = exp(y(9))*params(6)^(1-params(7));
T78 = T74*(y(4)-y(1))^params(7);
lhs =y(6);
rhs =T17*((y(7)+params(5)*exp(y(9))+params(3)*params(2))*y(1)-params(4)*y(8))+y(1)*params(2)*(1-params(1));
residual(1)= lhs-rhs;
lhs =y(8)*exp(y(11));
rhs =params(4)*y(2)+y(7)*(y(6)-y(1))+params(3)*(y(6)-params(2)*y(1))-params(5)*exp(y(9))*y(1);
residual(2)= lhs-rhs;
lhs =T63;
rhs =y(7)-y(17)/params(4);
residual(3)= lhs-rhs;
lhs =y(10);
rhs =2*exp(y(9)*10);
residual(4)= lhs-rhs;
lhs =y(13);
rhs =T78;
residual(5)= lhs-rhs;
lhs =y(14);
rhs =y(13)+y(1)*exp(y(9))*(params(5)+params(12));
residual(6)= lhs-rhs;
lhs =y(15);
rhs =params(5)*exp(y(9))*y(1)+y(13)-params(3)*(y(6)-params(2)*y(1))+y(1)*exp(y(9))*params(12);
residual(7)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(params(9)*y(3)+x(it_, 1));
residual(8)= lhs-rhs;
lhs =exp(y(11));
rhs =exp(params(10)*y(5)+x(it_, 2));
residual(9)= lhs-rhs;
lhs =y(12);
rhs =y(7)-y(17)/params(4);
residual(10)= lhs-rhs;
lhs =y(16);
rhs =params(3)*(y(6)-params(2)*y(1));
residual(11)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(11, 20);

  %
  % Jacobian matrix
  %

T121 = getPowerDeriv(y(4)-y(1),params(7),1);
T132 = getPowerDeriv((y(10)-y(6))/params(6),params(7)-1,1);
T138 = (y(7)+params(3)-y(17)/params(4))*(y(7)+params(3)-y(17)/params(4));
  g1(1,1)=(-(params(2)*(1-params(1))+T17*(y(7)+params(5)*exp(y(9))+params(3)*params(2))));
  g1(1,6)=1;
  g1(1,7)=(-(((y(7)+params(5)*exp(y(9))+params(3)*params(2))*y(1)-params(4)*y(8))*(-params(1))/T138+T17*y(1)));
  g1(1,17)=(-(((y(7)+params(5)*exp(y(9))+params(3)*params(2))*y(1)-params(4)*y(8))*(-(params(1)*(-(1/params(4)))))/T138));
  g1(1,8)=(-(T17*(-params(4))));
  g1(1,9)=(-(T17*params(5)*exp(y(9))*y(1)));
  g1(2,1)=(-((-y(7))+params(3)*(-params(2))-params(5)*exp(y(9))));
  g1(2,6)=(-(y(7)+params(3)));
  g1(2,7)=(-(y(6)-y(1)));
  g1(2,2)=(-params(4));
  g1(2,8)=exp(y(11));
  g1(2,9)=params(5)*exp(y(9))*y(1);
  g1(2,11)=y(8)*exp(y(11));
  g1(3,6)=params(7)*exp(y(18))*(-1)/params(6)*T132/params(4);
  g1(3,7)=(-1);
  g1(3,17)=1/params(4);
  g1(3,18)=T63;
  g1(3,10)=params(7)*exp(y(18))*T132*1/params(6)/params(4);
  g1(4,9)=(-(2*10*exp(y(9)*10)));
  g1(4,10)=1;
  g1(5,1)=(-(T74*(-T121)));
  g1(5,9)=(-T78);
  g1(5,4)=(-(T74*T121));
  g1(5,13)=1;
  g1(6,1)=(-(exp(y(9))*(params(5)+params(12))));
  g1(6,9)=(-(y(1)*exp(y(9))*(params(5)+params(12))));
  g1(6,13)=(-1);
  g1(6,14)=1;
  g1(7,1)=(-(exp(y(9))*params(12)+params(5)*exp(y(9))-params(3)*(-params(2))));
  g1(7,6)=params(3);
  g1(7,9)=(-(params(5)*exp(y(9))*y(1)+y(1)*exp(y(9))*params(12)));
  g1(7,13)=(-1);
  g1(7,15)=1;
  g1(8,3)=(-(params(9)*exp(params(9)*y(3)+x(it_, 1))));
  g1(8,9)=exp(y(9));
  g1(8,19)=(-exp(params(9)*y(3)+x(it_, 1)));
  g1(9,5)=(-(params(10)*exp(params(10)*y(5)+x(it_, 2))));
  g1(9,11)=exp(y(11));
  g1(9,20)=(-exp(params(10)*y(5)+x(it_, 2)));
  g1(10,7)=(-1);
  g1(10,17)=1/params(4);
  g1(10,12)=1;
  g1(11,1)=(-(params(3)*(-params(2))));
  g1(11,6)=(-params(3));
  g1(11,16)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],11,400);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],11,8000);
end
end
end
end
