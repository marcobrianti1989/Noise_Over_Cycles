function [residual, g1, g2, g3] = BGP_FEB2019_static(y, x, params)
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

residual = zeros( 6, 1);

%
% Model equations
%

T28 = y(3)^params(8);
T57 = (params(7)*(y(4)+exp(y(5))*T28)-(y(4)+exp(y(5))*T28)*params(7)*params(2))^(-params(1));
T67 = y(3)^params(4);
lhs =y(2);
rhs =(1+(1-y(3))*params(5)*params(6))/(y(3)+(1-y(3))*params(5));
residual(1)= lhs-rhs;
lhs =y(4);
rhs =y(4)*(1-params(9))+params(3)*exp(y(5))*T28;
residual(2)= lhs-rhs;
lhs =y(5);
rhs =y(5)*params(12)+params(13)*x(1);
residual(3)= lhs-rhs;
lhs =y(6);
rhs =y(6)*params(14)-params(15)*x(2);
residual(4)= lhs-rhs;
lhs =T57;
rhs =T57*y(2)*(y(3)+(1-y(3))*params(5))*params(11)*params(10)*exp(y(6))*T67;
residual(5)= lhs-rhs;
lhs =y(1);
rhs =y(4)+exp(y(5))*T28;
residual(6)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(6, 6);

  %
  % Jacobian matrix
  %

T84 = getPowerDeriv(y(3),params(8),1);
T91 = getPowerDeriv(params(7)*(y(4)+exp(y(5))*T28)-(y(4)+exp(y(5))*T28)*params(7)*params(2),(-params(1)),1);
T92 = (params(7)*exp(y(5))*T84-params(7)*params(2)*exp(y(5))*T84)*T91;
  g1(1,2)=1;
  g1(1,3)=(-(((y(3)+(1-y(3))*params(5))*params(6)*(-params(5))-(1+(1-y(3))*params(5)*params(6))*(1-params(5)))/((y(3)+(1-y(3))*params(5))*(y(3)+(1-y(3))*params(5)))));
  g1(2,3)=(-(params(3)*exp(y(5))*T84));
  g1(2,4)=1-(1-params(9));
  g1(2,5)=(-(params(3)*exp(y(5))*T28));
  g1(3,5)=1-params(12);
  g1(4,6)=1-params(14);
  g1(5,2)=(-(T67*T57*(y(3)+(1-y(3))*params(5))*params(11)*params(10)*exp(y(6))));
  g1(5,3)=T92-(T67*(y(2)*(y(3)+(1-y(3))*params(5))*params(11)*params(10)*exp(y(6))*T92+T57*y(2)*params(11)*params(10)*exp(y(6))*(1-params(5)))+T57*y(2)*(y(3)+(1-y(3))*params(5))*params(11)*params(10)*exp(y(6))*getPowerDeriv(y(3),params(4),1));
  g1(5,4)=T91*(params(7)-params(7)*params(2))-T67*y(2)*(y(3)+(1-y(3))*params(5))*params(11)*params(10)*exp(y(6))*T91*(params(7)-params(7)*params(2));
  g1(5,5)=T91*(params(7)*exp(y(5))*T28-exp(y(5))*T28*params(7)*params(2))-T67*y(2)*(y(3)+(1-y(3))*params(5))*params(11)*params(10)*exp(y(6))*T91*(params(7)*exp(y(5))*T28-exp(y(5))*T28*params(7)*params(2));
  g1(5,6)=(-(T57*y(2)*(y(3)+(1-y(3))*params(5))*params(11)*params(10)*exp(y(6))*T67));
  g1(6,1)=1;
  g1(6,3)=(-(exp(y(5))*T84));
  g1(6,4)=(-1);
  g1(6,5)=(-(exp(y(5))*T28));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,36);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,216);
end
end
end
end
