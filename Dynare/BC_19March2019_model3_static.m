function [residual, g1, g2, g3] = BC_19March2019_model3_static(y, x, params)
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

residual = zeros( 11, 1);

%
% Model equations
%

T16 = params(1)*exp(y(10))*y(4)^params(5);
T19 = y(5)^(1-params(5));
T20 = T16*T19;
lhs =y(2);
rhs =T20*params(2)-params(4);
residual(1)= lhs-rhs;
lhs =y(1);
rhs =y(1)*params(14)-params(3)*log(y(4));
residual(2)= lhs-rhs;
lhs =y(2);
rhs =y(5)*y(6)+y(4)*y(9);
residual(3)= lhs-rhs;
lhs =params(10)*y(5)^params(11);
rhs =y(6);
residual(4)= lhs-rhs;
lhs =y(3)+y(8);
rhs =y(7);
residual(5)= lhs-rhs;
lhs =y(6)*params(5)*y(5);
rhs =y(4)*(1-params(5))*y(9);
residual(6)= lhs-rhs;
lhs =1;
rhs =params(6)*(1+y(9)-params(8))*exp(y(11));
residual(7)= lhs-rhs;
lhs =y(7);
rhs =T20;
residual(8)= lhs-rhs;
lhs =y(4);
rhs =y(8)+y(4)*(1-params(8));
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(12)+x(1);
residual(10)= lhs-rhs;
lhs =y(11);
rhs =y(11)*params(13)-x(2);
residual(11)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(11, 11);

  %
  % Jacobian matrix
  %

T86 = T19*params(1)*exp(y(10))*getPowerDeriv(y(4),params(5),1);
T97 = T16*getPowerDeriv(y(5),1-params(5),1);
  g1(1,2)=1;
  g1(1,4)=(-(params(2)*T86));
  g1(1,5)=(-(params(2)*T97));
  g1(1,10)=(-(T20*params(2)));
  g1(2,1)=1-params(14);
  g1(2,4)=params(3)*1/y(4);
  g1(3,2)=1;
  g1(3,4)=(-y(9));
  g1(3,5)=(-y(6));
  g1(3,6)=(-y(5));
  g1(3,9)=(-y(4));
  g1(4,5)=params(10)*getPowerDeriv(y(5),params(11),1);
  g1(4,6)=(-1);
  g1(5,3)=1;
  g1(5,7)=(-1);
  g1(5,8)=1;
  g1(6,4)=(-((1-params(5))*y(9)));
  g1(6,5)=params(5)*y(6);
  g1(6,6)=params(5)*y(5);
  g1(6,9)=(-(y(4)*(1-params(5))));
  g1(7,9)=(-(params(6)*exp(y(11))));
  g1(7,11)=(-(params(6)*(1+y(9)-params(8))*exp(y(11))));
  g1(8,4)=(-T86);
  g1(8,5)=(-T97);
  g1(8,7)=1;
  g1(8,10)=(-T20);
  g1(9,4)=1-(1-params(8));
  g1(9,8)=(-1);
  g1(10,10)=1-params(12);
  g1(11,11)=1-params(13);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],11,121);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],11,1331);
end
end
end
end
