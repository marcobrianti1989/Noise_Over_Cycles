function [residual, g1, g2, g3] = JQ_static(y, x, params)
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

residual = zeros( 13, 1);

%
% Model equations
%

T39 = params(3)*exp(y(13))*y(3)^(params(3)-1);
T42 = y(7)^(1-params(3));
T56 = y(3)^params(3);
T59 = y(7)^(-params(3));
lhs =y(1)/y(2);
rhs =params(1)/(1-y(3));
residual(1)= lhs-rhs;
lhs =1;
rhs =params(5)*(y(4)-params(2))/(1-params(2));
residual(2)= lhs-rhs;
residual(3) = y(1)*y(3)+y(5)-y(5)/y(4)+y(6)-y(2);
lhs =T39*T42;
rhs =y(1)*1/(1-y(8)*y(9));
residual(4)= lhs-rhs;
lhs =y(10)*(1-params(4)+exp(y(13))*(1-params(3))*(1-y(8)*y(9))*T56*T59)+y(9)*y(8)*exp(y(12));
rhs =1;
residual(5)= lhs-rhs;
lhs =y(8)*exp(y(12))+y(4)*y(10)+y(9)*y(4)*(1-params(2))/(y(4)-params(2));
rhs =1;
residual(6)= lhs-rhs;
residual(7) = y(5)/y(4)+y(7)*(1-params(4))+T42*exp(y(13))*T56-y(1)*y(3)-y(5)-y(7)-y(11);
lhs =exp(y(12))*(y(7)-(1-params(2))*y(5)/(y(4)-params(2)));
rhs =T42*exp(y(13))*T56;
residual(8)= lhs-rhs;
lhs =y(10);
rhs =params(5);
residual(9)= lhs-rhs;
lhs =y(11);
rhs =y(6)+params(8)*(y(6)-params(6))^2;
residual(10)= lhs-rhs;
lhs =y(9);
rhs =1+(y(6)-params(6))*2*params(8);
residual(11)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(13)*params(9)+x(1));
residual(12)= lhs-rhs;
lhs =exp(y(12));
rhs =exp(y(12)*params(10)+params(7)*(1-params(10))+x(2));
residual(13)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(13, 13);

  %
  % Jacobian matrix
  %

T133 = getPowerDeriv(y(3),params(3),1);
T168 = getPowerDeriv(y(7),1-params(3),1);
  g1(1,1)=1/y(2);
  g1(1,2)=(-y(1))/(y(2)*y(2));
  g1(1,3)=(-(params(1)/((1-y(3))*(1-y(3)))));
  g1(2,4)=(-(params(5)/(1-params(2))));
  g1(3,1)=y(3);
  g1(3,2)=(-1);
  g1(3,3)=y(1);
  g1(3,4)=(-((-y(5))/(y(4)*y(4))));
  g1(3,5)=1-1/y(4);
  g1(3,6)=1;
  g1(4,1)=(-(1/(1-y(8)*y(9))));
  g1(4,3)=T42*params(3)*exp(y(13))*getPowerDeriv(y(3),params(3)-1,1);
  g1(4,7)=T39*T168;
  g1(4,8)=(-(y(1)*y(9)/((1-y(8)*y(9))*(1-y(8)*y(9)))));
  g1(4,9)=(-(y(1)*y(8)/((1-y(8)*y(9))*(1-y(8)*y(9)))));
  g1(4,13)=T39*T42;
  g1(5,3)=y(10)*T59*exp(y(13))*(1-params(3))*(1-y(8)*y(9))*T133;
  g1(5,7)=y(10)*exp(y(13))*(1-params(3))*(1-y(8)*y(9))*T56*getPowerDeriv(y(7),(-params(3)),1);
  g1(5,8)=y(10)*T59*T56*exp(y(13))*(1-params(3))*(-y(9))+y(9)*exp(y(12));
  g1(5,9)=y(8)*exp(y(12))+y(10)*T59*T56*exp(y(13))*(1-params(3))*(-y(8));
  g1(5,10)=1-params(4)+exp(y(13))*(1-params(3))*(1-y(8)*y(9))*T56*T59;
  g1(5,12)=y(9)*y(8)*exp(y(12));
  g1(5,13)=y(10)*exp(y(13))*(1-params(3))*(1-y(8)*y(9))*T56*T59;
  g1(6,4)=y(10)+y(9)*((y(4)-params(2))*(1-params(2))-y(4)*(1-params(2)))/((y(4)-params(2))*(y(4)-params(2)));
  g1(6,8)=exp(y(12));
  g1(6,9)=y(4)*(1-params(2))/(y(4)-params(2));
  g1(6,10)=y(4);
  g1(6,12)=y(8)*exp(y(12));
  g1(7,1)=(-y(3));
  g1(7,3)=T42*exp(y(13))*T133-y(1);
  g1(7,4)=(-y(5))/(y(4)*y(4));
  g1(7,5)=1/y(4)-1;
  g1(7,7)=1-params(4)+exp(y(13))*T56*T168-1;
  g1(7,11)=(-1);
  g1(7,13)=T42*exp(y(13))*T56;
  g1(8,3)=(-(T42*exp(y(13))*T133));
  g1(8,4)=exp(y(12))*(-((-((1-params(2))*y(5)))/((y(4)-params(2))*(y(4)-params(2)))));
  g1(8,5)=exp(y(12))*(-((1-params(2))/(y(4)-params(2))));
  g1(8,7)=exp(y(12))-exp(y(13))*T56*T168;
  g1(8,12)=exp(y(12))*(y(7)-(1-params(2))*y(5)/(y(4)-params(2)));
  g1(8,13)=(-(T42*exp(y(13))*T56));
  g1(9,10)=1;
  g1(10,6)=(-(1+params(8)*2*(y(6)-params(6))));
  g1(10,11)=1;
  g1(11,6)=(-(2*params(8)));
  g1(11,9)=1;
  g1(12,13)=exp(y(13))-params(9)*exp(y(13)*params(9)+x(1));
  g1(13,12)=exp(y(12))-params(10)*exp(y(12)*params(10)+params(7)*(1-params(10))+x(2));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],13,169);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],13,2197);
end
end
end
end
