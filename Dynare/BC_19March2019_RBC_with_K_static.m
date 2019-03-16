function [residual, g1, g2, g3] = BC_19March2019_RBC_with_K_static(y, x, params)
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

T13 = params(2)*y(2)^(params(2)-1);
T40 = 0.5*params(1)*exp(y(6))*y(2)^params(2);
lhs =y(5);
rhs =T13;
residual(1)= lhs-rhs;
lhs =y(1)+y(4);
rhs =y(3);
residual(2)= lhs-rhs;
lhs =1;
rhs =params(3)*(1+T13-params(5))*exp(y(7));
residual(3)= lhs-rhs;
lhs =y(3);
rhs =T40;
residual(4)= lhs-rhs;
lhs =y(2);
rhs =y(4)+y(2)*(1-params(5));
residual(5)= lhs-rhs;
lhs =y(6);
rhs =y(6)*params(7)+x(1);
residual(6)= lhs-rhs;
lhs =y(7);
rhs =y(7)*params(8)-x(2);
residual(7)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(7, 7);

  %
  % Jacobian matrix
  %

T57 = params(2)*getPowerDeriv(y(2),params(2)-1,1);
  g1(1,2)=(-T57);
  g1(1,5)=1;
  g1(2,1)=1;
  g1(2,3)=(-1);
  g1(2,4)=1;
  g1(3,2)=(-(exp(y(7))*params(3)*T57));
  g1(3,7)=(-(params(3)*(1+T13-params(5))*exp(y(7))));
  g1(4,2)=(-(0.5*params(1)*exp(y(6))*getPowerDeriv(y(2),params(2),1)));
  g1(4,3)=1;
  g1(4,6)=(-T40);
  g1(5,2)=1-(1-params(5));
  g1(5,4)=(-1);
  g1(6,6)=1-params(7);
  g1(7,7)=1-params(8);
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
