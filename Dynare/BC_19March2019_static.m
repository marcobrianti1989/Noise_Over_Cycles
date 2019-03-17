function [residual, g1, g2, g3] = BC_19March2019_static(y, x, params)
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

T11 = 1-y(1)^2;
T16 = y(5)^(1+params(8));
T17 = params(7)*T11*T16;
T26 = y(5)^params(5);
lhs =T17;
rhs =y(1)*params(3)*params(1)*exp(y(3))*T26-params(4)*exp(y(4));
residual(1)= lhs-rhs;
residual(2) = y(1)*T26-T17-params(2)/(T26*T11*params(1)*0.5);
lhs =y(2);
rhs =T26*(1-y(1));
residual(3)= lhs-rhs;
lhs =y(3);
rhs =y(3)*params(9)+x(1);
residual(4)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(10)-x(2);
residual(5)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(5, 5);

  %
  % Jacobian matrix
  %

T78 = params(7)*T11*getPowerDeriv(y(5),1+params(8),1);
T79 = getPowerDeriv(y(5),params(5),1);
  g1(1,1)=T16*params(7)*(-(2*y(1)))-params(3)*params(1)*exp(y(3))*T26;
  g1(1,3)=(-(y(1)*params(3)*params(1)*exp(y(3))*T26));
  g1(1,4)=params(4)*exp(y(4));
  g1(1,5)=T78-y(1)*params(3)*params(1)*exp(y(3))*T79;
  g1(2,1)=T26-T16*params(7)*(-(2*y(1)))-(-(params(2)*T26*params(1)*0.5*(-(2*y(1)))))/(T26*T11*params(1)*0.5*T26*T11*params(1)*0.5);
  g1(2,5)=y(1)*T79-T78-(-(params(2)*T11*params(1)*0.5*T79))/(T26*T11*params(1)*0.5*T26*T11*params(1)*0.5);
  g1(3,1)=T26;
  g1(3,2)=1;
  g1(3,5)=(-((1-y(1))*T79));
  g1(4,3)=1-params(9);
  g1(5,4)=1-params(10);
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
