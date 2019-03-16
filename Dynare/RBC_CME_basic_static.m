function [residual, g1, g2, g3] = RBC_CME_basic_static(y, x, params)
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

residual = zeros( 7, 1);

%
% Model equations
%

Rbar__ = 1/params(1);
T12 = params(1)/y(2);
T18 = y(3)^(params(3)-1);
T22 = params(3)*y(4)*T18+1-params(2);
T30 = y(3)^params(3);
T31 = y(4)*T30;
lhs =1/y(2);
rhs =T12*T22;
residual(1)= lhs-rhs;
lhs =1/y(2);
rhs =T12*y(6)/y(7);
residual(2)= lhs-rhs;
lhs =T31;
rhs =y(2)+y(3)-y(3)*(1-params(2));
residual(3)= lhs-rhs;
lhs =y(1);
rhs =T31;
residual(4)= lhs-rhs;
lhs =y(6)/Rbar__;
rhs =(y(7)/params(6))^params(5);
residual(5)= lhs-rhs;
lhs =y(4);
rhs =exp(y(5));
residual(6)= lhs-rhs;
lhs =y(5);
rhs =y(5)*params(4)+x(1);
residual(7)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(7, 7);

  %
  % Jacobian matrix
  %

T66 = y(4)*getPowerDeriv(y(3),params(3),1);
  g1(1,2)=(-1)/(y(2)*y(2))-T22*(-params(1))/(y(2)*y(2));
  g1(1,3)=(-(T12*params(3)*y(4)*getPowerDeriv(y(3),params(3)-1,1)));
  g1(1,4)=(-(T12*params(3)*T18));
  g1(2,2)=(-1)/(y(2)*y(2))-y(6)/y(7)*(-params(1))/(y(2)*y(2));
  g1(2,6)=(-(T12*1/y(7)));
  g1(2,7)=(-(T12*(-y(6))/(y(7)*y(7))));
  g1(3,2)=(-1);
  g1(3,3)=T66-(1-(1-params(2)));
  g1(3,4)=T30;
  g1(4,1)=1;
  g1(4,3)=(-T66);
  g1(4,4)=(-T30);
  g1(5,6)=1/Rbar__;
  g1(5,7)=(-(1/params(6)*getPowerDeriv(y(7)/params(6),params(5),1)));
  g1(6,4)=1;
  g1(6,5)=(-exp(y(5)));
  g1(7,5)=1-params(4);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,49);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,343);
end
end
end
end
