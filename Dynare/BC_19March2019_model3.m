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
M_.fname = 'BC_19March2019_model3';
M_.dynare_version = '4.5.4';
oo_.dynare_version = '4.5.4';
options_.dynare_version = '4.5.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('BC_19March2019_model3.log');
M_.exo_names = 'epsz';
M_.exo_names_tex = 'epsz';
M_.exo_names_long = 'epsz';
M_.exo_names = char(M_.exo_names, 'epsF');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsF');
M_.exo_names_long = char(M_.exo_names_long, 'epsF');
M_.endo_names = 'f';
M_.endo_names_tex = 'f';
M_.endo_names_long = 'f';
M_.endo_names = char(M_.endo_names, 'b');
M_.endo_names_tex = char(M_.endo_names_tex, 'b');
M_.endo_names_long = char(M_.endo_names_long, 'b');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names_long = char(M_.endo_names_long, 'n');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'logz');
M_.endo_names_tex = char(M_.endo_names_tex, 'logz');
M_.endo_names_long = char(M_.endo_names_long, 'logz');
M_.endo_names = char(M_.endo_names, 'logF');
M_.endo_names_tex = char(M_.endo_names_tex, 'logF');
M_.endo_names_long = char(M_.endo_names_long, 'logF');
M_.endo_partitions = struct();
M_.param_names = 'Z';
M_.param_names_tex = 'Z';
M_.param_names_long = 'Z';
M_.param_names = char(M_.param_names, 'EPS');
M_.param_names_tex = char(M_.param_names_tex, 'EPS');
M_.param_names_long = char(M_.param_names_long, 'EPS');
M_.param_names = char(M_.param_names, 'ETA');
M_.param_names_tex = char(M_.param_names_tex, 'ETA');
M_.param_names_long = char(M_.param_names_long, 'ETA');
M_.param_names = char(M_.param_names, 'FF');
M_.param_names_tex = char(M_.param_names_tex, 'FF');
M_.param_names_long = char(M_.param_names_long, 'FF');
M_.param_names = char(M_.param_names, 'ALP');
M_.param_names_tex = char(M_.param_names_tex, 'ALP');
M_.param_names_long = char(M_.param_names_long, 'ALP');
M_.param_names = char(M_.param_names, 'BET');
M_.param_names_tex = char(M_.param_names_tex, 'BET');
M_.param_names_long = char(M_.param_names_long, 'BET');
M_.param_names = char(M_.param_names, 'GAM');
M_.param_names_tex = char(M_.param_names_tex, 'GAM');
M_.param_names_long = char(M_.param_names_long, 'GAM');
M_.param_names = char(M_.param_names, 'DEL');
M_.param_names_tex = char(M_.param_names_tex, 'DEL');
M_.param_names_long = char(M_.param_names_long, 'DEL');
M_.param_names = char(M_.param_names, 'SIG');
M_.param_names_tex = char(M_.param_names_tex, 'SIG');
M_.param_names_long = char(M_.param_names_long, 'SIG');
M_.param_names = char(M_.param_names, 'PSI');
M_.param_names_tex = char(M_.param_names_tex, 'PSI');
M_.param_names_long = char(M_.param_names_long, 'PSI');
M_.param_names = char(M_.param_names, 'CHI');
M_.param_names_tex = char(M_.param_names_tex, 'CHI');
M_.param_names_long = char(M_.param_names_long, 'CHI');
M_.param_names = char(M_.param_names, 'RHOZ');
M_.param_names_tex = char(M_.param_names_tex, 'RHOZ');
M_.param_names_long = char(M_.param_names_long, 'RHOZ');
M_.param_names = char(M_.param_names, 'RHOF');
M_.param_names_tex = char(M_.param_names_tex, 'RHOF');
M_.param_names_long = char(M_.param_names_long, 'RHOF');
M_.param_names = char(M_.param_names, 'RHOf');
M_.param_names_tex = char(M_.param_names_tex, 'RHOf');
M_.param_names_long = char(M_.param_names_long, 'RHOf');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 11;
M_.param_nbr = 14;
M_.orig_endo_nbr = 11;
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
erase_compiled_function('BC_19March2019_model3_static');
erase_compiled_function('BC_19March2019_model3_dynamic');
M_.orig_eq_nbr = 11;
M_.eq_nbr = 11;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 6 0;
 0 7 0;
 2 8 17;
 3 9 18;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 0;
 0 14 19;
 4 15 0;
 5 16 0;]';
M_.nstatic = 5;
M_.nfwrd   = 1;
M_.npred   = 3;
M_.nboth   = 2;
M_.nsfwrd   = 3;
M_.nspred   = 5;
M_.ndynamic   = 6;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(11, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(14, 1);
M_.NNZDerivatives = [39; -1; -1];
M_.params( 1 ) = 1;
Z = M_.params( 1 );
M_.params( 3 ) = 0.1;
ETA = M_.params( 3 );
M_.params( 2 ) = 0.8;
EPS = M_.params( 2 );
M_.params( 4 ) = 1;
FF = M_.params( 4 );
M_.params( 5 ) = 0.6666666666666666;
ALP = M_.params( 5 );
M_.params( 6 ) = 0.99;
BET = M_.params( 6 );
M_.params( 7 ) = 0.8;
GAM = M_.params( 7 );
M_.params( 8 ) = 0.05;
DEL = M_.params( 8 );
M_.params( 9 ) = 2;
SIG = M_.params( 9 );
M_.params( 10 ) = 1;
PSI = M_.params( 10 );
M_.params( 11 ) = 4;
CHI = M_.params( 11 );
M_.params( 12 ) = 0.95;
RHOZ = M_.params( 12 );
M_.params( 13 ) = 0.95;
RHOF = M_.params( 13 );
M_.params( 14 ) = 0.5;
RHOf = M_.params( 14 );
rss    = 1/BET - 1 + DEL;
css    = 0.0172;
nss    = 0.3704;
kss    = 1/rss*ALP/(1-ALP)*PSI*nss^(1+CHI);
wss    = PSI*nss^CHI;
yss    = Z*kss^ALP*nss^(1-ALP);
iss    = kss*DEL;
bss    = wss*nss + rss*kss;
logzss = 0;
logFss = 0;
fss    = 0;
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 4 ) = kss;
oo_.steady_state( 3 ) = css;
oo_.steady_state( 9 ) = rss;
oo_.steady_state( 7 ) = yss;
oo_.steady_state( 8 ) = iss;
oo_.steady_state( 2 ) = bss;
oo_.steady_state( 5 ) = nss;
oo_.steady_state( 6 ) = wss;
oo_.steady_state( 10 ) = 0;
oo_.steady_state( 11 ) = 0;
oo_.steady_state( 1 ) = 0;
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
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 60;
options_.order = 1;
var_list_ = char('y','c','i','k','r','w','n','b');
info = stoch_simul(var_list_);
save('BC_19March2019_model3_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('BC_19March2019_model3_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('BC_19March2019_model3_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('BC_19March2019_model3_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('BC_19March2019_model3_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('BC_19March2019_model3_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('BC_19March2019_model3_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
