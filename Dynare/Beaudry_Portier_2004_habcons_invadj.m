%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'Beaudry_Portier_2004_habcons_invadj';
M_.dynare_version = '4.5.4';
oo_.dynare_version = '4.5.4';
options_.dynare_version = '4.5.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('Beaudry_Portier_2004_habcons_invadj.log');
M_.exo_names = 'eps_X';
M_.exo_names_tex = 'eps\_X';
M_.exo_names_long = 'eps_X';
M_.exo_names = char(M_.exo_names, 'eps_K');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_K');
M_.exo_names_long = char(M_.exo_names_long, 'eps_K');
M_.exo_names = char(M_.exo_names, 'eps_S');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_S');
M_.exo_names_long = char(M_.exo_names_long, 'eps_S');
M_.endo_names = 'C';
M_.endo_names_tex = 'C';
M_.endo_names_long = 'C';
M_.endo_names = char(M_.endo_names, 'K');
M_.endo_names_tex = char(M_.endo_names_tex, 'K');
M_.endo_names_long = char(M_.endo_names_long, 'K');
M_.endo_names = char(M_.endo_names, 'LK');
M_.endo_names_tex = char(M_.endo_names_tex, 'LK');
M_.endo_names_long = char(M_.endo_names_long, 'LK');
M_.endo_names = char(M_.endo_names, 'LX');
M_.endo_names_tex = char(M_.endo_names_tex, 'LX');
M_.endo_names_long = char(M_.endo_names_long, 'LX');
M_.endo_names = char(M_.endo_names, 'X');
M_.endo_names_tex = char(M_.endo_names_tex, 'X');
M_.endo_names_long = char(M_.endo_names_long, 'X');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, 'L');
M_.endo_names_long = char(M_.endo_names_long, 'L');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, 'I');
M_.endo_names_long = char(M_.endo_names_long, 'I');
M_.endo_names = char(M_.endo_names, 'MU');
M_.endo_names_tex = char(M_.endo_names_tex, 'MU');
M_.endo_names_long = char(M_.endo_names_long, 'MU');
M_.endo_names = char(M_.endo_names, 'dI');
M_.endo_names_tex = char(M_.endo_names_tex, 'dI');
M_.endo_names_long = char(M_.endo_names_long, 'dI');
M_.endo_names = char(M_.endo_names, 'dX');
M_.endo_names_tex = char(M_.endo_names_tex, 'dX');
M_.endo_names_long = char(M_.endo_names_long, 'dX');
M_.endo_names = char(M_.endo_names, 'LOGTHETX');
M_.endo_names_tex = char(M_.endo_names_tex, 'LOGTHETX');
M_.endo_names_long = char(M_.endo_names_long, 'LOGTHETX');
M_.endo_names = char(M_.endo_names, 'LOGTHETK');
M_.endo_names_tex = char(M_.endo_names_tex, 'LOGTHETK');
M_.endo_names_long = char(M_.endo_names_long, 'LOGTHETK');
M_.endo_names = char(M_.endo_names, 'LOGSENT');
M_.endo_names_tex = char(M_.endo_names_tex, 'LOGSENT');
M_.endo_names_long = char(M_.endo_names_long, 'LOGSENT');
M_.endo_partitions = struct();
M_.param_names = 'del';
M_.param_names_tex = 'del';
M_.param_names_long = 'del';
M_.param_names = char(M_.param_names, 'bet');
M_.param_names_tex = char(M_.param_names_tex, 'bet');
M_.param_names_long = char(M_.param_names_long, 'bet');
M_.param_names = char(M_.param_names, 'thetk');
M_.param_names_tex = char(M_.param_names_tex, 'thetk');
M_.param_names_long = char(M_.param_names_long, 'thetk');
M_.param_names = char(M_.param_names, 'thetx');
M_.param_names_tex = char(M_.param_names_tex, 'thetx');
M_.param_names_long = char(M_.param_names_long, 'thetx');
M_.param_names = char(M_.param_names, 'a');
M_.param_names_tex = char(M_.param_names_tex, 'a');
M_.param_names_long = char(M_.param_names_long, 'a');
M_.param_names = char(M_.param_names, 'v');
M_.param_names_tex = char(M_.param_names_tex, 'v');
M_.param_names_long = char(M_.param_names_long, 'v');
M_.param_names = char(M_.param_names, 'alp');
M_.param_names_tex = char(M_.param_names_tex, 'alp');
M_.param_names_long = char(M_.param_names_long, 'alp');
M_.param_names = char(M_.param_names, 'gam');
M_.param_names_tex = char(M_.param_names_tex, 'gam');
M_.param_names_long = char(M_.param_names_long, 'gam');
M_.param_names = char(M_.param_names, 'sig');
M_.param_names_tex = char(M_.param_names_tex, 'sig');
M_.param_names_long = char(M_.param_names_long, 'sig');
M_.param_names = char(M_.param_names, 'tau');
M_.param_names_tex = char(M_.param_names_tex, 'tau');
M_.param_names_long = char(M_.param_names_long, 'tau');
M_.param_names = char(M_.param_names, 'v0');
M_.param_names_tex = char(M_.param_names_tex, 'v0');
M_.param_names_long = char(M_.param_names_long, 'v0');
M_.param_names = char(M_.param_names, 'rhox');
M_.param_names_tex = char(M_.param_names_tex, 'rhox');
M_.param_names_long = char(M_.param_names_long, 'rhox');
M_.param_names = char(M_.param_names, 'rhok');
M_.param_names_tex = char(M_.param_names_tex, 'rhok');
M_.param_names_long = char(M_.param_names_long, 'rhok');
M_.param_names = char(M_.param_names, 'rhos');
M_.param_names_tex = char(M_.param_names_tex, 'rhos');
M_.param_names_long = char(M_.param_names_long, 'rhos');
M_.param_names = char(M_.param_names, 'extX');
M_.param_names_tex = char(M_.param_names_tex, 'extX');
M_.param_names_long = char(M_.param_names_long, 'extX');
M_.param_names = char(M_.param_names, 'extK');
M_.param_names_tex = char(M_.param_names_tex, 'extK');
M_.param_names_long = char(M_.param_names_long, 'extK');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 13;
M_.param_nbr = 16;
M_.orig_endo_nbr = 13;
M_.aux_vars = [];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
