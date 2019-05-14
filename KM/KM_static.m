function [residual, g1, g2, g3] = KM_static(y, x, params)
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

residual = zeros( 14, 1);

%
% Model equations
%

T16 = params(1)/(y(2)+params(3)-y(2)/params(4));
T54 = exp(y(4))*params(7)*((y(5)-y(1))/params(6))^(params(7)-1)/params(4);
T61 = exp(y(4))*params(6)^(1-params(7));
T63 = T61*(y(5)-y(1))^params(7);
lhs =y(1);
rhs =T16*(y(1)*(y(2)+params(5)*exp(y(4))+params(3)*params(2))-params(4)*y(3))+y(1)*params(2)*(1-params(1));
residual(1)= lhs-rhs;
lhs =y(3)*exp(y(6));
rhs =params(4)*y(3)+params(3)*(y(1)-y(1)*params(2))-y(1)*params(5)*exp(y(4));
residual(2)= lhs-rhs;
lhs =T54;
rhs =y(2)-y(2)/params(4);
residual(3)= lhs-rhs;
lhs =y(5);
rhs =2;
residual(4)= lhs-rhs;
lhs =y(8);
rhs =T63;
residual(5)= lhs-rhs;
lhs =y(9);
rhs =y(8)+y(1)*exp(y(4))*(params(5)+params(12));
residual(6)= lhs-rhs;
lhs =y(10);
rhs =y(1)*params(5)*exp(y(4))+y(8)-params(3)*(y(1)-y(1)*params(2))+y(1)*exp(y(4))*params(12);
residual(7)= lhs-rhs;
lhs =y(12);
rhs =y(10)-y(1)*exp(y(4))*params(12);
residual(8)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(4)*params(9)+x(1));
residual(9)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(6)*params(10)+x(2));
residual(10)= lhs-rhs;
lhs =y(7);
rhs =y(2)-y(2)/params(4);
residual(11)= lhs-rhs;
lhs =y(14);
rhs =y(3)/y(1);
residual(12)= lhs-rhs;
lhs =y(11);
rhs =params(3)*(y(1)-y(1)*params(2));
residual(13)= lhs-rhs;
lhs =y(13);
rhs =y(1)*y(2)-params(4)*y(3);
residual(14)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(14, 14);

  %
  % Jacobian matrix
  %

T113 = getPowerDeriv((y(5)-y(1))/params(6),params(7)-1,1);
T117 = getPowerDeriv(y(5)-y(1),params(7),1);
  g1(1,1)=1-(params(2)*(1-params(1))+T16*(y(2)+params(5)*exp(y(4))+params(3)*params(2)));
  g1(1,2)=(-((y(1)*(y(2)+params(5)*exp(y(4))+params(3)*params(2))-params(4)*y(3))*(-(params(1)*(1-1/params(4))))/((y(2)+params(3)-y(2)/params(4))*(y(2)+params(3)-y(2)/params(4)))+y(1)*T16));
  g1(1,3)=(-(T16*(-params(4))));
  g1(1,4)=(-(T16*y(1)*params(5)*exp(y(4))));
  g1(2,1)=(-(params(3)*(1-params(2))-params(5)*exp(y(4))));
  g1(2,3)=exp(y(6))-params(4);
  g1(2,4)=y(1)*params(5)*exp(y(4));
  g1(2,6)=y(3)*exp(y(6));
  g1(3,1)=exp(y(4))*params(7)*(-1)/params(6)*T113/params(4);
  g1(3,2)=(-(1-1/params(4)));
  g1(3,4)=T54;
  g1(3,5)=exp(y(4))*params(7)*T113*1/params(6)/params(4);
  g1(4,5)=1;
  g1(5,1)=(-(T61*(-T117)));
  g1(5,4)=(-T63);
  g1(5,5)=(-(T61*T117));
  g1(5,8)=1;
  g1(6,1)=(-(exp(y(4))*(params(5)+params(12))));
  g1(6,4)=(-(y(1)*exp(y(4))*(params(5)+params(12))));
  g1(6,8)=(-1);
  g1(6,9)=1;
  g1(7,1)=(-(exp(y(4))*params(12)+params(5)*exp(y(4))-params(3)*(1-params(2))));
  g1(7,4)=(-(y(1)*params(5)*exp(y(4))+y(1)*exp(y(4))*params(12)));
  g1(7,8)=(-1);
  g1(7,10)=1;
  g1(8,1)=exp(y(4))*params(12);
  g1(8,4)=y(1)*exp(y(4))*params(12);
  g1(8,10)=(-1);
  g1(8,12)=1;
  g1(9,4)=exp(y(4))-params(9)*exp(y(4)*params(9)+x(1));
  g1(10,6)=exp(y(6))-params(10)*exp(y(6)*params(10)+x(2));
  g1(11,2)=(-(1-1/params(4)));
  g1(11,7)=1;
  g1(12,1)=(-((-y(3))/(y(1)*y(1))));
  g1(12,3)=(-(1/y(1)));
  g1(12,14)=1;
  g1(13,1)=(-(params(3)*(1-params(2))));
  g1(13,11)=1;
  g1(14,1)=(-y(2));
  g1(14,2)=(-y(1));
  g1(14,3)=params(4);
  g1(14,13)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],14,196);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],14,2744);
end
end
end
end
