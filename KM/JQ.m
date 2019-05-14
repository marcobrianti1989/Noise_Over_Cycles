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
M_.fname = 'JQ';
M_.dynare_version = '4.5.6';
oo_.dynare_version = '4.5.6';
options_.dynare_version = '4.5.6';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('JQ.log');
M_.exo_names = 'ez';
M_.exo_names_tex = 'ez';
M_.exo_names_long = 'ez';
M_.exo_names = char(M_.exo_names, 'exi');
M_.exo_names_tex = char(M_.exo_names_tex, 'exi');
M_.exo_names_long = char(M_.exo_names_long, 'exi');
M_.endo_names = 'w';
M_.endo_names_tex = 'w';
M_.endo_names_long = 'w';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'b');
M_.endo_names_tex = char(M_.endo_names_tex, 'b');
M_.endo_names_long = char(M_.endo_names_long, 'b');
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'mu');
M_.endo_names_tex = char(M_.endo_names_tex, 'mu');
M_.endo_names_long = char(M_.endo_names_long, 'mu');
M_.endo_names = char(M_.endo_names, 'phip');
M_.endo_names_tex = char(M_.endo_names_tex, 'phip');
M_.endo_names_long = char(M_.endo_names_long, 'phip');
M_.endo_names = char(M_.endo_names, 'm');
M_.endo_names_tex = char(M_.endo_names_tex, 'm');
M_.endo_names_long = char(M_.endo_names_long, 'm');
M_.endo_names = char(M_.endo_names, 'phi');
M_.endo_names_tex = char(M_.endo_names_tex, 'phi');
M_.endo_names_long = char(M_.endo_names_long, 'phi');
M_.endo_names = char(M_.endo_names, 'xi');
M_.endo_names_tex = char(M_.endo_names_tex, 'xi');
M_.endo_names_long = char(M_.endo_names_long, 'xi');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_partitions = struct();
M_.param_names = 'ALPHA';
M_.param_names_tex = 'ALPHA';
M_.param_names_long = 'ALPHA';
M_.param_names = char(M_.param_names, 'TAU');
M_.param_names_tex = char(M_.param_names_tex, 'TAU');
M_.param_names_long = char(M_.param_names_long, 'TAU');
M_.param_names = char(M_.param_names, 'THETA');
M_.param_names_tex = char(M_.param_names_tex, 'THETA');
M_.param_names_long = char(M_.param_names_long, 'THETA');
M_.param_names = char(M_.param_names, 'DELTA');
M_.param_names_tex = char(M_.param_names_tex, 'DELTA');
M_.param_names_long = char(M_.param_names_long, 'DELTA');
M_.param_names = char(M_.param_names, 'BETA');
M_.param_names_tex = char(M_.param_names_tex, 'BETA');
M_.param_names_long = char(M_.param_names_long, 'BETA');
M_.param_names = char(M_.param_names, 'DBAR');
M_.param_names_tex = char(M_.param_names_tex, 'DBAR');
M_.param_names_long = char(M_.param_names_long, 'DBAR');
M_.param_names = char(M_.param_names, 'XIBAR');
M_.param_names_tex = char(M_.param_names_tex, 'XIBAR');
M_.param_names_long = char(M_.param_names_long, 'XIBAR');
M_.param_names = char(M_.param_names, 'KAPPA');
M_.param_names_tex = char(M_.param_names_tex, 'KAPPA');
M_.param_names_long = char(M_.param_names_long, 'KAPPA');
M_.param_names = char(M_.param_names, 'RHOZ');
M_.param_names_tex = char(M_.param_names_tex, 'RHOZ');
M_.param_names_long = char(M_.param_names_long, 'RHOZ');
M_.param_names = char(M_.param_names, 'RHOXI');
M_.param_names_tex = char(M_.param_names_tex, 'RHOXI');
M_.param_names_long = char(M_.param_names_long, 'RHOXI');
M_.param_names = char(M_.param_names, 'VOLZ');
M_.param_names_tex = char(M_.param_names_tex, 'VOLZ');
M_.param_names_long = char(M_.param_names_long, 'VOLZ');
M_.param_names = char(M_.param_names, 'VOLXI');
M_.param_names_tex = char(M_.param_names_tex, 'VOLXI');
M_.param_names_long = char(M_.param_names_long, 'VOLXI');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 13;
M_.param_nbr = 12;
M_.orig_endo_nbr = 13;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('JQ_static');
erase_compiled_function('JQ_dynamic');
M_.orig_eq_nbr = 13;
M_.eq_nbr = 13;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 5 0;
 0 6 18;
 0 7 0;
 0 8 0;
 1 9 0;
 0 10 0;
 2 11 0;
 0 12 19;
 0 13 20;
 0 14 0;
 0 15 0;
 3 16 0;
 4 17 0;]';
M_.nstatic = 6;
M_.nfwrd   = 3;
M_.npred   = 4;
M_.nboth   = 0;
M_.nsfwrd   = 3;
M_.nspred   = 4;
M_.ndynamic   = 7;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(13, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(12, 1);
M_.NNZDerivatives = [64; -1; -1];
M_.params( 1 ) = 1.8834;
ALPHA = M_.params( 1 );
M_.params( 2 ) = 0.35;
TAU = M_.params( 2 );
M_.params( 3 ) = 0.64;
THETA = M_.params( 3 );
M_.params( 4 ) = 0.025;
DELTA = M_.params( 4 );
M_.params( 5 ) = 0.9825;
BETA = M_.params( 5 );
M_.params( 7 ) = (-1.811554096556235);
XIBAR = M_.params( 7 );
M_.params( 6 ) = 0;
DBAR = M_.params( 6 );
M_.params( 8 ) = 0.146;
KAPPA = M_.params( 8 );
M_.params( 9 ) = .9;
RHOZ = M_.params( 9 );
M_.params( 10 ) = .9;
RHOXI = M_.params( 10 );
M_.params( 11 ) = 1;
VOLZ = M_.params( 11 );
M_.params( 12 ) = 1;
VOLXI = M_.params( 12 );
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(11))^2;
M_.Sigma_e(2, 2) = (M_.params(12))^2;
M_.Sigma_e(1, 2) = 0;
M_.Sigma_e(2, 1) = M_.Sigma_e(1, 2);
M_.sigma_e_is_diagonal = 0;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 25;
options_.order = 1;
options_.periods = 200;
var_list_ = char('k');
info = stoch_simul(var_list_);
save('JQ_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('JQ_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('JQ_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('JQ_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('JQ_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('JQ_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('JQ_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
