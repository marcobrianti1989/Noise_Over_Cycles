function [residual, g1, g2, g3] = Beaudry_Portier_2004_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
T21 = exp(y(13))*params(3)*params(8)*y(7)^(params(8)-1);
T23 = exp(y(14))*params(9)*T21^(-1);
T35 = params(8)*params(3)*exp(y(18))*y(15)^(params(8)-1);
T47 = exp(y(17))*params(4)*y(16)^params(7);
T55 = params(5)*T47^params(6)+(1-params(5))*y(6)^params(6);
T58 = (1-params(5))*exp(y(14))*params(2)*T55^(-1);
T60 = y(6)^(params(6)-1);
T69 = params(4)*exp(y(12))*y(8)^params(7);
T75 = params(5)*T69^params(6)+(1-params(5))*y(1)^params(6);
T77 = params(5)*T75^(-1);
T78 = T69^(params(6)-1);
T84 = y(8)^(params(7)-1);
lhs =T23;
rhs =params(9)*params(2)*(1-params(1))*T35^(-1)+T58*T60;
residual(1)= lhs-rhs;
lhs =exp(y(14))*params(9);
rhs =params(7)*params(4)*exp(y(12))*T77*T78*T84;
residual(2)= lhs-rhs;
lhs =exp(y(14))*y(5);
rhs =T75^(1/params(6));
residual(3)= lhs-rhs;
lhs =y(10);
rhs =y(7)+y(8);
residual(4)= lhs-rhs;
lhs =y(6);
rhs =(1-params(1))*y(1)+y(11);
residual(5)= lhs-rhs;
lhs =y(11);
rhs =params(3)*y(7)^params(8);
residual(6)= lhs-rhs;
lhs =y(9);
rhs =T69;
residual(7)= lhs-rhs;
lhs =y(12);
rhs =params(10)*y(2)+x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(13);
rhs =params(11)*y(3)+x(it_, 2);
residual(9)= lhs-rhs;
lhs =y(14);
rhs =params(12)*y(4)-x(it_, 3);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 21);

  %
  % Jacobian matrix
  %

T123 = (1-params(5))*getPowerDeriv(y(1),params(6),1);
T124 = getPowerDeriv(T75,(-1),1);
T133 = getPowerDeriv(T75,1/params(6),1);
T139 = getPowerDeriv(T55,(-1),1);
T150 = getPowerDeriv(T21,(-1),1);
T158 = getPowerDeriv(T35,(-1),1);
T163 = params(4)*exp(y(12))*getPowerDeriv(y(8),params(7),1);
T164 = getPowerDeriv(T69,params(6),1);
T169 = getPowerDeriv(T69,params(6)-1,1);
T187 = getPowerDeriv(T47,params(6),1);
  g1(1,6)=(-(T60*(1-params(5))*exp(y(14))*params(2)*(1-params(5))*getPowerDeriv(y(6),params(6),1)*T139+T58*getPowerDeriv(y(6),params(6)-1,1)));
  g1(1,7)=exp(y(14))*params(9)*exp(y(13))*params(3)*params(8)*getPowerDeriv(y(7),params(8)-1,1)*T150;
  g1(1,15)=(-(params(9)*params(2)*(1-params(1))*params(8)*params(3)*exp(y(18))*getPowerDeriv(y(15),params(8)-1,1)*T158));
  g1(1,16)=(-(T60*(1-params(5))*exp(y(14))*params(2)*T139*params(5)*exp(y(17))*params(4)*getPowerDeriv(y(16),params(7),1)*T187));
  g1(1,17)=(-(T60*(1-params(5))*exp(y(14))*params(2)*T139*params(5)*T47*T187));
  g1(1,13)=exp(y(14))*params(9)*T21*T150;
  g1(1,18)=(-(params(9)*params(2)*(1-params(1))*T35*T158));
  g1(1,14)=T23-T58*T60;
  g1(2,1)=(-(T84*params(7)*params(4)*exp(y(12))*T78*params(5)*T123*T124));
  g1(2,8)=(-(T84*params(7)*params(4)*exp(y(12))*(T78*params(5)*T124*params(5)*T163*T164+T77*T163*T169)+params(7)*params(4)*exp(y(12))*T77*T78*getPowerDeriv(y(8),params(7)-1,1)));
  g1(2,12)=(-(T84*params(7)*params(4)*(exp(y(12))*T77*T78+exp(y(12))*(T78*params(5)*T124*params(5)*T69*T164+T77*T69*T169))));
  g1(2,14)=exp(y(14))*params(9);
  g1(3,5)=exp(y(14));
  g1(3,1)=(-(T123*T133));
  g1(3,8)=(-(T133*params(5)*T163*T164));
  g1(3,12)=(-(T133*params(5)*T69*T164));
  g1(3,14)=exp(y(14))*y(5);
  g1(4,7)=(-1);
  g1(4,8)=(-1);
  g1(4,10)=1;
  g1(5,1)=(-(1-params(1)));
  g1(5,6)=1;
  g1(5,11)=(-1);
  g1(6,7)=(-(params(3)*getPowerDeriv(y(7),params(8),1)));
  g1(6,11)=1;
  g1(7,8)=(-T163);
  g1(7,9)=1;
  g1(7,12)=(-T69);
  g1(8,2)=(-params(10));
  g1(8,12)=1;
  g1(8,19)=(-1);
  g1(9,3)=(-params(11));
  g1(9,13)=1;
  g1(9,20)=(-1);
  g1(10,4)=(-params(12));
  g1(10,14)=1;
  g1(10,21)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,441);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,9261);
end
end
end
end
