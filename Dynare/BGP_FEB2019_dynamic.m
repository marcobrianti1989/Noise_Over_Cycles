function [residual, g1, g2, g3] = BGP_FEB2019_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(6, 1);
T30 = y(1)^params(8);
T54 = exp(y(9))*y(7)^params(8);
T55 = y(8)+T54;
T72 = y(6)*(y(7)+(1-y(7))*params(5))*params(11)*params(10)*exp(y(10));
T78 = exp(y(13))*y(11)^params(8);
T82 = params(7)*(y(12)+T78)-T55*params(7)*params(2);
T83 = T82^(-params(1));
T86 = y(11)^params(4);
lhs =y(6);
rhs =(1+(1-y(7))*params(5)*params(6))/(y(7)+(1-y(7))*params(5));
residual(1)= lhs-rhs;
lhs =y(8);
rhs =(1-params(9))*y(2)+params(3)*exp(y(3))*T30;
residual(2)= lhs-rhs;
lhs =y(9);
rhs =y(3)*params(12)+params(13)*x(it_, 1);
residual(3)= lhs-rhs;
lhs =y(10);
rhs =params(14)*y(4)-params(15)*x(it_, 2);
residual(4)= lhs-rhs;
lhs =(params(7)*T55-params(7)*params(2)*(y(2)+exp(y(3))*T30))^(-params(1));
rhs =T72*T83*T86;
residual(5)= lhs-rhs;
lhs =y(5);
rhs =T55;
residual(6)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(6, 15);

  %
  % Jacobian matrix
  %

T94 = getPowerDeriv(y(1),params(8),1);
T100 = getPowerDeriv(params(7)*T55-params(7)*params(2)*(y(2)+exp(y(3))*T30),(-params(1)),1);
T112 = exp(y(9))*getPowerDeriv(y(7),params(8),1);
T119 = getPowerDeriv(T82,(-params(1)),1);
  g1(1,6)=1;
  g1(1,7)=(-(((y(7)+(1-y(7))*params(5))*params(6)*(-params(5))-(1+(1-y(7))*params(5)*params(6))*(1-params(5)))/((y(7)+(1-y(7))*params(5))*(y(7)+(1-y(7))*params(5)))));
  g1(2,1)=(-(params(3)*exp(y(3))*T94));
  g1(2,2)=(-(1-params(9)));
  g1(2,8)=1;
  g1(2,3)=(-(params(3)*exp(y(3))*T30));
  g1(3,3)=(-params(12));
  g1(3,9)=1;
  g1(3,14)=(-params(13));
  g1(4,4)=(-params(14));
  g1(4,10)=1;
  g1(4,15)=params(15);
  g1(5,6)=(-(T86*(y(7)+(1-y(7))*params(5))*params(11)*params(10)*exp(y(10))*T83));
  g1(5,1)=(-(params(7)*params(2)*exp(y(3))*T94))*T100;
  g1(5,7)=T100*params(7)*T112-T86*(T83*y(6)*params(11)*params(10)*exp(y(10))*(1-params(5))+T72*(-(params(7)*params(2)*T112))*T119);
  g1(5,11)=(-(T86*T72*T119*params(7)*exp(y(13))*getPowerDeriv(y(11),params(8),1)+T72*T83*getPowerDeriv(y(11),params(4),1)));
  g1(5,2)=T100*(-(params(7)*params(2)));
  g1(5,8)=params(7)*T100-T86*T72*T119*(-(params(7)*params(2)));
  g1(5,12)=(-(T86*T72*params(7)*T119));
  g1(5,3)=T100*(-(params(7)*params(2)*exp(y(3))*T30));
  g1(5,9)=T100*params(7)*T54-T86*T72*T119*(-(T54*params(7)*params(2)));
  g1(5,13)=(-(T86*T72*T119*params(7)*T78));
  g1(5,10)=(-(T72*T83*T86));
  g1(6,5)=1;
  g1(6,7)=(-T112);
  g1(6,8)=(-1);
  g1(6,9)=(-T54);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,225);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,3375);
end
end
end
end
