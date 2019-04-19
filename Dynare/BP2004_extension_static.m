function [residual, g1, g2, g3] = BP2004_extension_static(y, x, params)
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

T16 = y(1)^(1-params(6));
T20 = T16/(y(1)-y(1)*params(9));
T26 = y(2)^(params(6)-1);
T32 = params(10)/2;
T34 = y(7)^2;
T35 = y(2)^2;
T40 = 1+T32*(T34/T35-params(1)^2);
T48 = y(5)^(params(6)-1);
T57 = y(4)^params(15);
T62 = y(4)^params(7);
T63 = exp(y(11))*T57*params(4)*T62;
T68 = params(5)*T63^params(6)+(1-params(5))*y(2)^params(6);
T81 = (y(7)/y(2)-params(1))^2;
T88 = y(3)^params(16);
T93 = y(3)^params(8);
T101 = y(3)^(params(8)-1);
T115 = y(4)^(params(7)-1);
lhs =y(8);
rhs =exp(y(13))*params(2)*(T20*(1-params(5))*T26+y(8)*(1-params(1))*T40);
residual(1)= lhs-rhs;
lhs =params(11);
rhs =T20*params(5)*T48*y(10);
residual(2)= lhs-rhs;
lhs =y(1);
rhs =T68^(1/params(6));
residual(3)= lhs-rhs;
lhs =y(6);
rhs =y(4)+y(3);
residual(4)= lhs-rhs;
lhs =y(2);
rhs =y(7)+y(2)*(1-params(1))-y(2)*T32*T81;
residual(5)= lhs-rhs;
lhs =y(7);
rhs =exp(y(12))*T88*params(3)*T93;
residual(6)= lhs-rhs;
lhs =y(9);
rhs =params(3)*T88*exp(y(12))*params(8)*T101;
residual(7)= lhs-rhs;
lhs =y(8);
rhs =params(11)*(y(9)*(1-params(10)*(y(7)/y(2)-params(1))))^(-1);
residual(8)= lhs-rhs;
lhs =y(5);
rhs =T63;
residual(9)= lhs-rhs;
lhs =y(10);
rhs =params(4)*T57*exp(y(11))*params(7)*T115;
residual(10)= lhs-rhs;
lhs =y(11);
rhs =y(11)*params(12)+x(1);
residual(11)= lhs-rhs;
lhs =y(12);
rhs =y(12)*params(13)+x(2);
residual(12)= lhs-rhs;
lhs =y(13);
rhs =y(13)*params(14)+x(3);
residual(13)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(13, 13);

  %
  % Jacobian matrix
  %

T139 = ((y(1)-y(1)*params(9))*getPowerDeriv(y(1),1-params(6),1)-T16*(1-params(9)))/((y(1)-y(1)*params(9))*(y(1)-y(1)*params(9)));
T162 = getPowerDeriv(T68,1/params(6),1);
T178 = getPowerDeriv(y(9)*(1-params(10)*(y(7)/y(2)-params(1))),(-1),1);
T182 = getPowerDeriv(y(3),params(16),1);
T197 = getPowerDeriv(y(4),params(15),1);
T203 = T62*params(4)*exp(y(11))*T197+exp(y(11))*T57*params(4)*getPowerDeriv(y(4),params(7),1);
T204 = getPowerDeriv(T63,params(6),1);
  g1(1,1)=(-(exp(y(13))*params(2)*T26*(1-params(5))*T139));
  g1(1,2)=(-(exp(y(13))*params(2)*(T20*(1-params(5))*getPowerDeriv(y(2),params(6)-1,1)+y(8)*(1-params(1))*T32*(-(T34*2*y(2)))/(T35*T35))));
  g1(1,7)=(-(exp(y(13))*params(2)*y(8)*(1-params(1))*T32*2*y(7)/T35));
  g1(1,8)=1-exp(y(13))*params(2)*(1-params(1))*T40;
  g1(1,13)=(-(exp(y(13))*params(2)*(T20*(1-params(5))*T26+y(8)*(1-params(1))*T40)));
  g1(2,1)=(-(y(10)*T48*params(5)*T139));
  g1(2,5)=(-(y(10)*T20*params(5)*getPowerDeriv(y(5),params(6)-1,1)));
  g1(2,10)=(-(T20*params(5)*T48));
  g1(3,1)=1;
  g1(3,2)=(-((1-params(5))*getPowerDeriv(y(2),params(6),1)*T162));
  g1(3,4)=(-(T162*params(5)*T203*T204));
  g1(3,11)=(-(T162*params(5)*T63*T204));
  g1(4,3)=(-1);
  g1(4,4)=(-1);
  g1(4,6)=1;
  g1(5,2)=1-(1-params(1)-(T32*T81+y(2)*T32*(-y(7))/(y(2)*y(2))*2*(y(7)/y(2)-params(1))));
  g1(5,7)=(-(1-y(2)*T32*2*(y(7)/y(2)-params(1))*1/y(2)));
  g1(6,3)=(-(T93*params(3)*exp(y(12))*T182+exp(y(12))*T88*params(3)*getPowerDeriv(y(3),params(8),1)));
  g1(6,7)=1;
  g1(6,12)=(-(exp(y(12))*T88*params(3)*T93));
  g1(7,3)=(-(T101*params(3)*exp(y(12))*params(8)*T182+params(3)*T88*exp(y(12))*params(8)*getPowerDeriv(y(3),params(8)-1,1)));
  g1(7,9)=1;
  g1(7,12)=(-(params(3)*T88*exp(y(12))*params(8)*T101));
  g1(8,2)=(-(params(11)*y(9)*(-(params(10)*(-y(7))/(y(2)*y(2))))*T178));
  g1(8,7)=(-(params(11)*T178*y(9)*(-(params(10)*1/y(2)))));
  g1(8,8)=1;
  g1(8,9)=(-(params(11)*(1-params(10)*(y(7)/y(2)-params(1)))*T178));
  g1(9,4)=(-T203);
  g1(9,5)=1;
  g1(9,11)=(-T63);
  g1(10,4)=(-(T115*params(4)*exp(y(11))*params(7)*T197+params(4)*T57*exp(y(11))*params(7)*getPowerDeriv(y(4),params(7)-1,1)));
  g1(10,10)=1;
  g1(10,11)=(-(params(4)*T57*exp(y(11))*params(7)*T115));
  g1(11,11)=1-params(12);
  g1(12,12)=1-params(13);
  g1(13,13)=1-params(14);
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
