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

T18 = exp(y(9))*params(3)*params(8)*y(3)^(params(8)-1);
T19 = T18^(-1);
T39 = exp(y(8))*params(4)*y(4)^params(7);
T46 = params(5)*T39^params(6)+(1-params(5))*y(2)^params(6);
T47 = T46^(-1);
T51 = y(2)^(params(6)-1);
T56 = T39^(params(6)-1);
T62 = y(4)^(params(7)-1);
lhs =params(9)*T19;
rhs =T19*params(9)*exp(y(10))*params(2)*(1-params(1))+(1-params(5))*exp(y(10))*params(2)*T47*T51;
residual(1)= lhs-rhs;
lhs =params(9);
rhs =params(7)*params(4)*exp(y(8))*params(5)*T47*T56*T62;
residual(2)= lhs-rhs;
lhs =y(1);
rhs =T46^(1/params(6));
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
rhs =T39;
residual(7)= lhs-rhs;
lhs =y(8);
rhs =y(8)*params(10)+x(1);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(11)+x(2);
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(12)+x(3);
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

  %
  % Jacobian matrix
  %

T97 = (1-params(5))*getPowerDeriv(y(2),params(6),1);
T98 = getPowerDeriv(T46,(-1),1);
T99 = T97*T98;
T114 = getPowerDeriv(T46,1/params(6),1);
T120 = getPowerDeriv(T18,(-1),1);
T121 = exp(y(9))*params(3)*params(8)*getPowerDeriv(y(3),params(8)-1,1)*T120;
T129 = exp(y(8))*params(4)*getPowerDeriv(y(4),params(7),1);
T130 = getPowerDeriv(T39,params(6),1);
T133 = T98*params(5)*T129*T130;
T139 = getPowerDeriv(T39,params(6)-1,1);
  g1(1,2)=(-(T51*(1-params(5))*exp(y(10))*params(2)*T99+(1-params(5))*exp(y(10))*params(2)*T47*getPowerDeriv(y(2),params(6)-1,1)));
  g1(1,3)=params(9)*T121-params(9)*exp(y(10))*params(2)*(1-params(1))*T121;
  g1(1,4)=(-(T51*(1-params(5))*exp(y(10))*params(2)*T133));
  g1(1,8)=(-(T51*(1-params(5))*exp(y(10))*params(2)*T98*params(5)*T39*T130));
  g1(1,9)=params(9)*T18*T120-params(9)*exp(y(10))*params(2)*(1-params(1))*T18*T120;
  g1(1,10)=(-(T19*params(9)*exp(y(10))*params(2)*(1-params(1))+(1-params(5))*exp(y(10))*params(2)*T47*T51));
  g1(2,2)=(-(T62*params(7)*params(4)*exp(y(8))*T56*params(5)*T99));
  g1(2,4)=(-(T62*params(7)*params(4)*exp(y(8))*(T56*params(5)*T133+params(5)*T47*T129*T139)+params(7)*params(4)*exp(y(8))*params(5)*T47*T56*getPowerDeriv(y(4),params(7)-1,1)));
  g1(2,8)=(-(T62*params(7)*params(4)*(exp(y(8))*params(5)*T47*T56+exp(y(8))*(T56*params(5)*T98*params(5)*T39*T130+params(5)*T47*T39*T139))));
  g1(3,1)=1;
  g1(3,2)=(-(T97*T114));
  g1(3,4)=(-(T114*params(5)*T129*T130));
  g1(3,8)=(-(T114*params(5)*T39*T130));
  g1(4,3)=(-1);
  g1(4,4)=(-1);
  g1(4,6)=1;
  g1(5,2)=1-(1-params(1));
  g1(5,7)=(-1);
  g1(6,3)=(-(params(3)*getPowerDeriv(y(3),params(8),1)));
  g1(6,7)=1;
  g1(7,4)=(-T129);
  g1(7,5)=1;
  g1(7,8)=(-T39);
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
