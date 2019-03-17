function [residual, g1, g2, g3] = BC_19March2019_model2_static(y, x, params)
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

residual = zeros( 8, 1);

%
% Model equations
%

T20 = y(3)^params(5);
T57 = exp(y(7))*params(1)*0.5*(1-y(1)^2);
lhs =y(6)*y(3);
rhs =params(3)*params(1)*exp(y(7))*y(1)*T20-params(4);
residual(1)= lhs-rhs;
lhs =T20*y(1)*params(1)*exp(y(7))-y(6)*y(3);
rhs =params(2);
residual(2)= lhs-rhs;
lhs =y(2)+y(5);
rhs =y(4)-params(2)*(1-y(1));
residual(3)= lhs-rhs;
lhs =1;
rhs =params(6)*(1+y(6)-params(8))*exp(y(8));
residual(4)= lhs-rhs;
lhs =y(4);
rhs =T20*T57;
residual(5)= lhs-rhs;
lhs =y(3)*(1-y(1));
rhs =y(5)+(1-y(1))*y(3)*(1-params(8));
residual(6)= lhs-rhs;
lhs =y(7);
rhs =y(7)*params(10)+x(1);
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(8)*params(11)-x(2);
residual(8)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(8, 8);

  %
  % Jacobian matrix
  %

T88 = getPowerDeriv(y(3),params(5),1);
  g1(1,1)=(-(params(3)*params(1)*exp(y(7))*T20));
  g1(1,3)=y(6)-params(3)*params(1)*exp(y(7))*y(1)*T88;
  g1(1,6)=y(3);
  g1(1,7)=(-(params(3)*params(1)*exp(y(7))*y(1)*T20));
  g1(2,1)=T20*params(1)*exp(y(7));
  g1(2,3)=y(1)*params(1)*exp(y(7))*T88-y(6);
  g1(2,6)=(-y(3));
  g1(2,7)=T20*y(1)*params(1)*exp(y(7));
  g1(3,1)=(-params(2));
  g1(3,2)=1;
  g1(3,4)=(-1);
  g1(3,5)=1;
  g1(4,6)=(-(params(6)*exp(y(8))));
  g1(4,8)=(-(params(6)*(1+y(6)-params(8))*exp(y(8))));
  g1(5,1)=(-(T20*exp(y(7))*params(1)*0.5*(-(2*y(1)))));
  g1(5,3)=(-(T57*T88));
  g1(5,4)=1;
  g1(5,7)=(-(T20*T57));
  g1(6,1)=(-y(3))-(-(y(3)*(1-params(8))));
  g1(6,3)=1-y(1)-(1-y(1))*(1-params(8));
  g1(6,5)=(-1);
  g1(7,7)=1-params(10);
  g1(8,8)=1-params(11);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,64);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,512);
end
end
end
end
