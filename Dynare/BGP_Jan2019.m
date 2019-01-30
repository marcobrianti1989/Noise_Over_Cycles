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
M_.fname = 'BGP_Jan2019';
M_.dynare_version = '4.5.4';
oo_.dynare_version = '4.5.4';
options_.dynare_version = '4.5.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('BGP_Jan2019.log');
M_.exo_names = 'eps_thet';
M_.exo_names_tex = 'eps\_thet';
M_.exo_names_long = 'eps_thet';
M_.exo_names = char(M_.exo_names, 'eps_zeta');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_zeta');
M_.exo_names_long = char(M_.exo_names_long, 'eps_zeta');
M_.endo_names = 'rp';
M_.endo_names_tex = 'rp';
M_.endo_names_long = 'rp';
M_.endo_names = char(M_.endo_names, 'e');
M_.endo_names_tex = char(M_.endo_names_tex, 'e');
M_.endo_names_long = char(M_.endo_names_long, 'e');
M_.endo_names = char(M_.endo_names, 'x');
M_.endo_names_tex = char(M_.endo_names_tex, 'x');
M_.endo_names_long = char(M_.endo_names_long, 'x');
M_.endo_names = char(M_.endo_names, 'logthet');
M_.endo_names_tex = char(M_.endo_names_tex, 'logthet');
M_.endo_names_long = char(M_.endo_names_long, 'logthet');
M_.endo_names = char(M_.endo_names, 'logzeta');
M_.endo_names_tex = char(M_.endo_names_tex, 'logzeta');
M_.endo_names_long = char(M_.endo_names_long, 'logzeta');
M_.endo_partitions = struct();
M_.param_names = 'OMEG';
M_.param_names_tex = 'OMEG';
M_.param_names_long = 'OMEG';
M_.param_names = char(M_.param_names, 'GAM');
M_.param_names_tex = char(M_.param_names_tex, 'GAM');
M_.param_names_long = char(M_.param_names_long, 'GAM');
M_.param_names = char(M_.param_names, 'PSI');
M_.param_names_tex = char(M_.param_names_tex, 'PSI');
M_.param_names_long = char(M_.param_names_long, 'PSI');
M_.param_names = char(M_.param_names, 'PHIE');
M_.param_names_tex = char(M_.param_names_tex, 'PHIE');
M_.param_names_long = char(M_.param_names_long, 'PHIE');
M_.param_names = char(M_.param_names, 'PHI');
M_.param_names_tex = char(M_.param_names_tex, 'PHI');
M_.param_names_long = char(M_.param_names_long, 'PHI');
M_.param_names = char(M_.param_names, 'PHIBIG');
M_.param_names_tex = char(M_.param_names_tex, 'PHIBIG');
M_.param_names_long = char(M_.param_names_long, 'PHIBIG');
M_.param_names = char(M_.param_names, 'S');
M_.param_names_tex = char(M_.param_names_tex, 'S');
M_.param_names_long = char(M_.param_names_long, 'S');
M_.param_names = char(M_.param_names, 'ALP');
M_.param_names_tex = char(M_.param_names_tex, 'ALP');
M_.param_names_long = char(M_.param_names_long, 'ALP');
M_.param_names = char(M_.param_names, 'DEL');
M_.param_names_tex = char(M_.param_names_tex, 'DEL');
M_.param_names_long = char(M_.param_names_long, 'DEL');
M_.param_names = char(M_.param_names, 'THET');
M_.param_names_tex = char(M_.param_names_tex, 'THET');
M_.param_names_long = char(M_.param_names_long, 'THET');
M_.param_names = char(M_.param_names, 'BET');
M_.param_names_tex = char(M_.param_names_tex, 'BET');
M_.param_names_long = char(M_.param_names_long, 'BET');
M_.param_names = char(M_.param_names, 'RHO_THET');
M_.param_names_tex = char(M_.param_names_tex, 'RHO\_THET');
M_.param_names_long = char(M_.param_names_long, 'RHO_THET');
M_.param_names = char(M_.param_names, 'SIGMA_THET');
M_.param_names_tex = char(M_.param_names_tex, 'SIGMA\_THET');
M_.param_names_long = char(M_.param_names_long, 'SIGMA_THET');
M_.param_names = char(M_.param_names, 'RHO_ZETA');
M_.param_names_tex = char(M_.param_names_tex, 'RHO\_ZETA');
M_.param_names_long = char(M_.param_names_long, 'RHO_ZETA');
M_.param_names = char(M_.param_names, 'SIGMA_ZETA');
M_.param_names_tex = char(M_.param_names_tex, 'SIGMA\_ZETA');
M_.param_names_long = char(M_.param_names_long, 'SIGMA_ZETA');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 5;
M_.param_nbr = 15;
M_.orig_endo_nbr = 5;
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
erase_compiled_function('BGP_Jan2019_static');
erase_compiled_function('BGP_Jan2019_dynamic');
M_.orig_eq_nbr = 5;
M_.eq_nbr = 5;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 5 0;
 1 6 10;
 2 7 11;
 3 8 12;
 4 9 0;]';
M_.nstatic = 1;
M_.nfwrd   = 0;
M_.npred   = 1;
M_.nboth   = 3;
M_.nsfwrd   = 3;
M_.nspred   = 4;
M_.ndynamic   = 4;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(5, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(15, 1);
M_.NNZDerivatives = [24; -1; -1];
M_.params( 1 ) = 0.2408;
OMEG = M_.params( 1 );
M_.params( 2 ) = 0.5876;
GAM = M_.params( 2 );
M_.params( 3 ) = 0.2994;
PSI = M_.params( 3 );
M_.params( 4 ) = 0.0467;
PHIE = M_.params( 4 );
M_.params( 5 ) = 0.8827;
PHI = M_.params( 5 );
M_.params( 6 ) = 0.0458;
PHIBIG = M_.params( 6 );
M_.params( 7 ) = 1;
S = M_.params( 7 );
M_.params( 8 ) = 0.6666666666666666;
ALP = M_.params( 8 );
M_.params( 9 ) = 0.05;
DEL = M_.params( 9 );
M_.params( 10 ) = 1;
THET = M_.params( 10 );
M_.params( 11 ) = 0.99;
BET = M_.params( 11 );
M_.params( 12 ) = 0.9;
RHO_THET = M_.params( 12 );
M_.params( 13 ) = 1;
SIGMA_THET = M_.params( 13 );
M_.params( 14 ) = 0.9;
RHO_ZETA = M_.params( 14 );
M_.params( 15 ) = 1;
SIGMA_ZETA = M_.params( 15 );
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
options_.solve_algo = 2;
options_.steady.maxit = 1000;
steady;
options_.irf = 50;
options_.order = 1;
var_list_ = char('logthet','e','x','rp');
info = stoch_simul(var_list_);
save('BGP_Jan2019_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('BGP_Jan2019_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('BGP_Jan2019_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('BGP_Jan2019_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('BGP_Jan2019_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('BGP_Jan2019_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('BGP_Jan2019_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
