function [residual, g1, g2, g3] = Beaudry_Portier_2004_static(y, x, params)
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

residual = zeros( 10, 1);

%
% Model equations
%

T21 = exp(y(9))*params(3)*params(8)*y(3)^(params(8)-1);
T22 = T21^(-1);
T40 = exp(y(8))*params(4)*y(4)^params(7);
T47 = params(5)*T40^params(6)+(1-params(5))*y(2)^params(6);
T48 = T47^(-1);
T52 = y(2)^(params(6)-1);
T57 = T40^(params(6)-1);
T63 = y(4)^(params(7)-1);
lhs =exp(y(10))*params(9)*T22;
rhs =T22*params(9)*params(2)*(1-params(1))+(1-params(5))*exp(y(10))*params(2)*T48*T52;
residual(1)= lhs-rhs;
lhs =exp(y(10))*params(9);
rhs =params(7)*params(4)*exp(y(8))*params(5)*T48*T57*T63;
residual(2)= lhs-rhs;
lhs =exp(y(10))*y(1);
rhs =T47^(1/params(6));
residual(3)= lhs-rhs;
lhs =y(6);
rhs =y(3)+y(4);
residual(4)= lhs-rhs;
lhs =y(2);
rhs =(1-params(1))*y(2)+y(7);
residual(5)= lhs-rhs;
lhs =y(7);
rhs =params(3)*y(3)^params(8);
residual(6)= lhs-rhs;
lhs =y(5);
rhs =T40;
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(8)*params(10)+x(1);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(11)+x(2);
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(12)-x(3);
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

  %
  % Jacobian matrix
  %

T99 = (1-params(5))*getPowerDeriv(y(2),params(6),1);
T100 = getPowerDeriv(T47,(-1),1);
T101 = T99*T100;
T116 = getPowerDeriv(T47,1/params(6),1);
T122 = getPowerDeriv(T21,(-1),1);
T123 = exp(y(9))*params(3)*params(8)*getPowerDeriv(y(3),params(8)-1,1)*T122;
T131 = exp(y(8))*params(4)*getPowerDeriv(y(4),params(7),1);
T132 = getPowerDeriv(T40,params(6),1);
T135 = T100*params(5)*T131*T132;
T141 = getPowerDeriv(T40,params(6)-1,1);
  g1(1,2)=(-(T52*(1-params(5))*exp(y(10))*params(2)*T101+(1-params(5))*exp(y(10))*params(2)*T48*getPowerDeriv(y(2),params(6)-1,1)));
  g1(1,3)=exp(y(10))*params(9)*T123-params(9)*params(2)*(1-params(1))*T123;
  g1(1,4)=(-(T52*(1-params(5))*exp(y(10))*params(2)*T135));
  g1(1,8)=(-(T52*(1-params(5))*exp(y(10))*params(2)*T100*params(5)*T40*T132));
  g1(1,9)=exp(y(10))*params(9)*T21*T122-params(9)*params(2)*(1-params(1))*T21*T122;
  g1(1,10)=exp(y(10))*params(9)*T22-(1-params(5))*exp(y(10))*params(2)*T48*T52;
  g1(2,2)=(-(T63*params(7)*params(4)*exp(y(8))*T57*params(5)*T101));
  g1(2,4)=(-(T63*params(7)*params(4)*exp(y(8))*(T57*params(5)*T135+params(5)*T48*T131*T141)+params(7)*params(4)*exp(y(8))*params(5)*T48*T57*getPowerDeriv(y(4),params(7)-1,1)));
  g1(2,8)=(-(T63*params(7)*params(4)*(exp(y(8))*params(5)*T48*T57+exp(y(8))*(T57*params(5)*T100*params(5)*T40*T132+params(5)*T48*T40*T141))));
  g1(2,10)=exp(y(10))*params(9);
  g1(3,1)=exp(y(10));
  g1(3,2)=(-(T99*T116));
  g1(3,4)=(-(T116*params(5)*T131*T132));
  g1(3,8)=(-(T116*params(5)*T40*T132));
  g1(3,10)=exp(y(10))*y(1);
  g1(4,3)=(-1);
  g1(4,4)=(-1);
  g1(4,6)=1;
  g1(5,2)=1-(1-params(1));
  g1(5,7)=(-1);
  g1(6,3)=(-(params(3)*getPowerDeriv(y(3),params(8),1)));
  g1(6,7)=1;
  g1(7,4)=(-T131);
  g1(7,5)=1;
  g1(7,8)=(-T40);
  g1(8,8)=1-params(10);
  g1(9,9)=1-params(11);
  g1(10,10)=1-params(12);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,100);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,1000);
end
end
end
end
