function [residual, g1, g2, g3] = BGP_Jan2019_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 5, 1);

%
% Model equations
%

T27 = exp(y(4))*y(2)^params(8);
T55 = (params(7)*(y(3)+T27)-(y(3)+T27)*params(7)*params(2))^(1-params(1))-1;
T63 = (y(2)+(1-y(2))*params(5))*params(11)*params(10)*exp(y(5))/exp(y(5));
lhs =y(1);
rhs =(1+(1-y(2))*params(5)*params(6))/(y(2)+(1-y(2))*params(5));
residual(1)= lhs-rhs;
lhs =y(3);
rhs =y(3)*(1-params(9))+T27;
residual(2)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(12)+params(13)*x(1);
residual(3)= lhs-rhs;
lhs =y(5);
rhs =y(5)*params(14)+params(15)*x(2);
residual(4)= lhs-rhs;
lhs =T55/(1-params(1));
rhs =T55*y(1)*T63/(1-params(1));
residual(5)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(5, 5);

  %
  % Jacobian matrix
  %

T81 = exp(y(4))*getPowerDeriv(y(2),params(8),1);
T86 = getPowerDeriv(params(7)*(y(3)+T27)-(y(3)+T27)*params(7)*params(2),1-params(1),1);
T87 = (params(7)*T81-params(7)*params(2)*T81)*T86;
  g1(1,1)=1;
  g1(1,2)=(-(((y(2)+(1-y(2))*params(5))*params(6)*(-params(5))-(1+(1-y(2))*params(5)*params(6))*(1-params(5)))/((y(2)+(1-y(2))*params(5))*(y(2)+(1-y(2))*params(5)))));
  g1(2,2)=(-T81);
  g1(2,3)=1-(1-params(9));
  g1(2,4)=(-T27);
  g1(3,4)=1-params(12);
  g1(4,5)=1-params(14);
  g1(5,1)=(-(T55*T63/(1-params(1))));
  g1(5,2)=T87/(1-params(1))-(y(1)*T63*T87+T55*y(1)*params(11)*params(10)*exp(y(5))/exp(y(5))*(1-params(5)))/(1-params(1));
  g1(5,3)=T86*(params(7)-params(7)*params(2))/(1-params(1))-y(1)*T63*T86*(params(7)-params(7)*params(2))/(1-params(1));
  g1(5,4)=T86*(T27*params(7)-T27*params(7)*params(2))/(1-params(1))-y(1)*T63*T86*(T27*params(7)-T27*params(7)*params(2))/(1-params(1));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],5,25);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],5,125);
end
end
end
end
