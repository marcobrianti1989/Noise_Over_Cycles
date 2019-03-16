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
M_.fname = 'Vito_code_RBC';
M_.dynare_version = '4.5.4';
oo_.dynare_version = '4.5.4';
options_.dynare_version = '4.5.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('Vito_code_RBC.log');
M_.exo_names = 'ez';
M_.exo_names_tex = 'ez';
M_.exo_names_long = 'ez';
M_.exo_names = char(M_.exo_names, 'eg');
M_.exo_names_tex = char(M_.exo_names_tex, 'eg');
M_.exo_names_long = char(M_.exo_names_long, 'eg');
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names_long = 'c';
M_.endo_names = char(M_.endo_names, 'g');
M_.endo_names_tex = char(M_.endo_names_tex, 'g');
M_.endo_names_long = char(M_.endo_names_long, 'g');
M_.endo_names = char(M_.endo_names, 'h');
M_.endo_names_tex = char(M_.endo_names_tex, 'h');
M_.endo_names_long = char(M_.endo_names_long, 'h');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_partitions = struct();
M_.param_names = 'ALPHA';
M_.param_names_tex = 'ALPHA';
M_.param_names_long = 'ALPHA';
M_.param_names = char(M_.param_names, 'BETA');
M_.param_names_tex = char(M_.param_names_tex, 'BETA');
M_.param_names_long = char(M_.param_names_long, 'BETA');
M_.param_names = char(M_.param_names, 'DELTA');
M_.param_names_tex = char(M_.param_names_tex, 'DELTA');
M_.param_names_long = char(M_.param_names_long, 'DELTA');
M_.param_names = char(M_.param_names, 'PHI');
M_.param_names_tex = char(M_.param_names_tex, 'PHI');
M_.param_names_long = char(M_.param_names_long, 'PHI');
M_.param_names = char(M_.param_names, 'RHOG');
M_.param_names_tex = char(M_.param_names_tex, 'RHOG');
M_.param_names_long = char(M_.param_names_long, 'RHOG');
M_.param_names = char(M_.param_names, 'RHOZ');
M_.param_names_tex = char(M_.param_names_tex, 'RHOZ');
M_.param_names_long = char(M_.param_names_long, 'RHOZ');
M_.param_names = char(M_.param_names, 'SG');
M_.param_names_tex = char(M_.param_names_tex, 'SG');
M_.param_names_long = char(M_.param_names_long, 'SG');
M_.param_names = char(M_.param_names, 'SC');
M_.param_names_tex = char(M_.param_names_tex, 'SC');
M_.param_names_long = char(M_.param_names_long, 'SC');
M_.param_names = char(M_.param_names, 'VOLG');
M_.param_names_tex = char(M_.param_names_tex, 'VOLG');
M_.param_names_long = char(M_.param_names_long, 'VOLG');
M_.param_names = char(M_.param_names, 'VOLZ');
M_.param_names_tex = char(M_.param_names_tex, 'VOLZ');
M_.param_names_long = char(M_.param_names_long, 'VOLZ');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 9;
M_.param_nbr = 10;
M_.orig_endo_nbr = 9;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('Vito_code_RBC_static');
erase_compiled_function('Vito_code_RBC_dynamic');
M_.orig_eq_nbr = 9;
M_.eq_nbr = 9;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 4 13;
 1 5 0;
 0 6 0;
 0 7 0;
 2 8 0;
 0 9 14;
 0 10 0;
 0 11 0;
 3 12 0;]';
M_.nstatic = 4;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(9, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(10, 1);
M_.NNZDerivatives = [31; -1; -1];
M_.params( 1 ) = 0.3333333333333333;
ALPHA = M_.params( 1 );
M_.params( 2 ) = 0.99;
BETA = M_.params( 2 );
M_.params( 3 ) = 0.025;
DELTA = M_.params( 3 );
M_.params( 4 ) = 0.25;
PHI = M_.params( 4 );
M_.params( 5 ) = 0.90;
RHOG = M_.params( 5 );
M_.params( 6 ) = 0.99;
RHOZ = M_.params( 6 );
M_.params( 7 ) = 0.20;
SG = M_.params( 7 );
M_.params( 8 ) = 0.65;
SC = M_.params( 8 );
M_.params( 9 ) = M_.params(7)^(-1)/100;
VOLG = M_.params( 9 );
M_.params( 10 ) = (1-M_.params(1))^(-1)/100;
VOLZ = M_.params( 10 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = 0;
oo_.steady_state( 2 ) = 0;
oo_.steady_state( 3 ) = 0;
oo_.steady_state( 4 ) = 0;
oo_.steady_state( 5 ) = 0;
oo_.steady_state( 6 ) = 0;
oo_.steady_state( 7 ) = 0;
oo_.steady_state( 8 ) = 0;
oo_.steady_state( 9 ) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(10))^2;
M_.Sigma_e(2, 2) = (M_.params(9))^2;
M_.Sigma_e(1, 2) = 0;
M_.Sigma_e(2, 1) = M_.Sigma_e(1, 2);
M_.sigma_e_is_diagonal = 0;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 40;
options_.nocorr = 1;
options_.order = 1;
options_.periods = 200;
var_list_ = char('y','c','i','h','w','k','r','z','g');
info = stoch_simul(var_list_);
save('Vito_code_RBC_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('Vito_code_RBC_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('Vito_code_RBC_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('Vito_code_RBC_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('Vito_code_RBC_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('Vito_code_RBC_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('Vito_code_RBC_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
