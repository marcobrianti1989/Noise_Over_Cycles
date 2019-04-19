function [residual, g1, g2, g3] = Beaudry_Portier_2004_habit_consumption_static(y, x, params)
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

T18 = (1-y(6)+y(6)*params(11))^params(9);
T19 = exp(y(10))*params(12)/T18;
T27 = exp(y(9))*y(3)^params(17)*params(3);
T31 = y(3)^(params(8)-1);
T33 = (T27*params(8)*T31)^(-1);
T40 = params(12)*exp(y(10))*(1-params(1))/T18;
T45 = y(1)^(1-params(6));
T55 = y(2)^(params(6)-1);
T64 = y(4)^params(16);
T69 = params(5)*T45*1/(y(1)-y(1)*params(10));
T72 = y(4)^params(7);
T73 = params(4)*T64*exp(y(8))*T72;
T74 = T73^(params(6)-1);
T76 = exp(y(8))*T69*T74;
T79 = params(7)*params(4)*T64*T76;
T81 = y(4)^(params(7)-1);
T89 = params(5)*y(5)^params(6)+(1-params(5))*y(2)^params(6);
T99 = y(3)^params(8);
lhs =T19*T33;
rhs =params(2)*(T33*T40+T45/(y(1)-y(1)*params(10))*(1-params(5))*T55);
residual(1)= lhs-rhs;
lhs =T19;
rhs =T79*T81;
residual(2)= lhs-rhs;
lhs =y(1);
rhs =T89^(1/params(6));
residual(3)= lhs-rhs;
lhs =y(6);
rhs =y(3)+y(4);
residual(4)= lhs-rhs;
lhs =y(2);
rhs =(1-params(1))*y(2)+y(7);
residual(5)= lhs-rhs;
lhs =y(7);
rhs =T27*T99;
residual(6)= lhs-rhs;
lhs =y(5);
rhs =T73;
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(8)*params(13)+x(1);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(14)+x(2);
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(15)-x(3);
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

  %
  % Jacobian matrix
  %

T118 = getPowerDeriv(y(1),1-params(6),1);
T148 = getPowerDeriv(T89,1/params(6),1);
T154 = params(3)*exp(y(9))*getPowerDeriv(y(3),params(17),1);
T160 = getPowerDeriv(T27*params(8)*T31,(-1),1);
T161 = (T31*params(8)*T154+T27*params(8)*getPowerDeriv(y(3),params(8)-1,1))*T160;
T171 = getPowerDeriv(y(4),params(16),1);
T177 = T72*params(4)*exp(y(8))*T171+params(4)*T64*exp(y(8))*getPowerDeriv(y(4),params(7),1);
T178 = getPowerDeriv(T73,params(6)-1,1);
T199 = ((-1)+params(11))*getPowerDeriv(1-y(6)+y(6)*params(11),params(9),1);
  g1(1,1)=(-(params(2)*T55*(1-params(5))*((y(1)-y(1)*params(10))*T118-T45*(1-params(10)))/((y(1)-y(1)*params(10))*(y(1)-y(1)*params(10)))));
  g1(1,2)=(-(params(2)*T45/(y(1)-y(1)*params(10))*(1-params(5))*getPowerDeriv(y(2),params(6)-1,1)));
  g1(1,3)=T19*T161-params(2)*T40*T161;
  g1(1,6)=T33*(-(exp(y(10))*params(12)*T199))/(T18*T18)-params(2)*T33*(-(params(12)*exp(y(10))*(1-params(1))*T199))/(T18*T18);
  g1(1,9)=T19*T27*params(8)*T31*T160-params(2)*T40*T27*params(8)*T31*T160;
  g1(1,10)=T19*T33-params(2)*T33*T40;
  g1(2,1)=(-(T81*params(7)*params(4)*T64*exp(y(8))*T74*params(5)*(1/(y(1)-y(1)*params(10))*T118+T45*(-(1-params(10)))/((y(1)-y(1)*params(10))*(y(1)-y(1)*params(10))))));
  g1(2,4)=(-(T81*params(7)*params(4)*(T76*T171+T64*exp(y(8))*T69*T177*T178)+T79*getPowerDeriv(y(4),params(7)-1,1)));
  g1(2,6)=(-(exp(y(10))*params(12)*T199))/(T18*T18);
  g1(2,8)=(-(T81*params(7)*params(4)*T64*(T76+exp(y(8))*T69*T73*T178)));
  g1(2,10)=T19;
  g1(3,1)=1;
  g1(3,2)=(-((1-params(5))*getPowerDeriv(y(2),params(6),1)*T148));
  g1(3,5)=(-(T148*params(5)*getPowerDeriv(y(5),params(6),1)));
  g1(4,3)=(-1);
  g1(4,4)=(-1);
  g1(4,6)=1;
  g1(5,2)=1-(1-params(1));
  g1(5,7)=(-1);
  g1(6,3)=(-(T99*T154+T27*getPowerDeriv(y(3),params(8),1)));
  g1(6,7)=1;
  g1(6,9)=(-(T27*T99));
  g1(7,4)=(-T177);
  g1(7,5)=1;
  g1(7,8)=(-T73);
  g1(8,8)=1-params(13);
  g1(9,9)=1-params(14);
  g1(10,10)=1-params(15);
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
