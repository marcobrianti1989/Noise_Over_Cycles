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
M_.fname = 'KM';
M_.dynare_version = '4.5.6';
oo_.dynare_version = '4.5.6';
options_.dynare_version = '4.5.6';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('KM.log');
M_.exo_names = 'ez';
M_.exo_names_tex = 'ez';
M_.exo_names_long = 'ez';
M_.exo_names = char(M_.exo_names, 'es');
M_.exo_names_tex = char(M_.exo_names_tex, 'es');
M_.exo_names_long = char(M_.exo_names_long, 'es');
M_.endo_names = 'k';
M_.endo_names_tex = 'k';
M_.endo_names_long = 'k';
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'b');
M_.endo_names_tex = char(M_.endo_names_tex, 'b');
M_.endo_names_long = char(M_.endo_names_long, 'b');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_names = char(M_.endo_names, 'kbar');
M_.endo_names_tex = char(M_.endo_names_tex, 'kbar');
M_.endo_names_long = char(M_.endo_names_long, 'kbar');
M_.endo_names = char(M_.endo_names, 's');
M_.endo_names_tex = char(M_.endo_names_tex, 's');
M_.endo_names_long = char(M_.endo_names_long, 's');
M_.endo_names = char(M_.endo_names, 'dp');
M_.endo_names_tex = char(M_.endo_names_tex, 'dp');
M_.endo_names_long = char(M_.endo_names_long, 'dp');
M_.endo_names = char(M_.endo_names, 'yg');
M_.endo_names_tex = char(M_.endo_names_tex, 'yg');
M_.endo_names_long = char(M_.endo_names_long, 'yg');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'cg');
M_.endo_names_tex = char(M_.endo_names_tex, 'cg');
M_.endo_names_long = char(M_.endo_names_long, 'cg');
M_.endo_names = char(M_.endo_names, 'cons');
M_.endo_names_tex = char(M_.endo_names_tex, 'cons');
M_.endo_names_long = char(M_.endo_names_long, 'cons');
M_.endo_partitions = struct();
M_.param_names = 'PI';
M_.param_names_tex = 'PI';
M_.param_names_long = 'PI';
M_.param_names = char(M_.param_names, 'LAMBDA');
M_.param_names_tex = char(M_.param_names_tex, 'LAMBDA');
M_.param_names_long = char(M_.param_names_long, 'LAMBDA');
M_.param_names = char(M_.param_names, 'PHI');
M_.param_names_tex = char(M_.param_names_tex, 'PHI');
M_.param_names_long = char(M_.param_names_long, 'PHI');
M_.param_names = char(M_.param_names, 'R');
M_.param_names_tex = char(M_.param_names_tex, 'R');
M_.param_names_long = char(M_.param_names_long, 'R');
M_.param_names = char(M_.param_names, 'A');
M_.param_names_tex = char(M_.param_names_tex, 'A');
M_.param_names_long = char(M_.param_names_long, 'A');
M_.param_names = char(M_.param_names, 'M');
M_.param_names_tex = char(M_.param_names_tex, 'M');
M_.param_names_long = char(M_.param_names_long, 'M');
M_.param_names = char(M_.param_names, 'ALPHA');
M_.param_names_tex = char(M_.param_names_tex, 'ALPHA');
M_.param_names_long = char(M_.param_names_long, 'ALPHA');
M_.param_names = char(M_.param_names, 'VOLZ');
M_.param_names_tex = char(M_.param_names_tex, 'VOLZ');
M_.param_names_long = char(M_.param_names_long, 'VOLZ');
M_.param_names = char(M_.param_names, 'RHOZ');
M_.param_names_tex = char(M_.param_names_tex, 'RHOZ');
M_.param_names_long = char(M_.param_names_long, 'RHOZ');
M_.param_names = char(M_.param_names, 'RHOS');
M_.param_names_tex = char(M_.param_names_tex, 'RHOS');
M_.param_names_long = char(M_.param_names_long, 'RHOS');
M_.param_names = char(M_.param_names, 'VOLS');
M_.param_names_tex = char(M_.param_names_tex, 'VOLS');
M_.param_names_long = char(M_.param_names_long, 'VOLS');
M_.param_names = char(M_.param_names, 'C');
M_.param_names_tex = char(M_.param_names_tex, 'C');
M_.param_names_long = char(M_.param_names_long, 'C');
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
erase_compiled_function('KM_static');
erase_compiled_function('KM_dynamic');
M_.orig_eq_nbr = 13;
M_.eq_nbr = 13;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 6 0;
 0 7 19;
 2 8 0;
 3 9 20;
 4 10 0;
 5 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 0 15 0;
 0 16 0;
 0 17 0;
 0 18 0;]';
M_.nstatic = 7;
M_.nfwrd   = 1;
M_.npred   = 4;
M_.nboth   = 1;
M_.nsfwrd   = 2;
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
oo_.steady_state = zeros(13, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(12, 1);
M_.NNZDerivatives = [52; -1; -1];
M_.params( 5 ) = .5;
A = M_.params( 5 );
M_.params( 3 ) = 1.4;
PHI = M_.params( 3 );
M_.params( 4 ) = 1.025;
R = M_.params( 4 );
M_.params( 2 ) = 0.9875;
LAMBDA = M_.params( 2 );
M_.params( 1 ) = 0.60;
PI = M_.params( 1 );
M_.params( 7 ) = 0.1666666666666667;
ALPHA = M_.params( 7 );
M_.params( 6 ) = 1;
M = M_.params( 6 );
M_.params( 8 ) = .1;
VOLZ = M_.params( 8 );
M_.params( 11 ) = 0.01;
VOLS = M_.params( 11 );
M_.params( 9 ) = .92;
RHOZ = M_.params( 9 );
M_.params( 10 ) = .92;
RHOS = M_.params( 10 );
M_.params( 12 ) = .5;
C = M_.params( 12 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 2 ) = M_.params(4)/(M_.params(4)-1)*(M_.params(1)*M_.params(5)-(1-M_.params(2))*(1-M_.params(4)+M_.params(4)*M_.params(1))*M_.params(3))/((1-M_.params(2))*(1-M_.params(4)+M_.params(4)*M_.params(1))+M_.params(1)*M_.params(2));
oo_.steady_state( 3 ) = (M_.params(5)+M_.params(2)*M_.params(3)-M_.params(3))/(M_.params(4)-1)*oo_.steady_state(1);
oo_.steady_state( 4 ) = 0;
oo_.steady_state( 5 ) = 2;
oo_.steady_state( 1 ) = oo_.steady_state(5)-M_.params(6)*(M_.params(4)*oo_.steady_state(2)/M_.params(7)-oo_.steady_state(2)/M_.params(7))^(1/(M_.params(7)-1));
oo_.steady_state( 6 ) = 0;
oo_.steady_state( 8 ) = M_.params(6)^(1-M_.params(7))*(oo_.steady_state(5)-oo_.steady_state(1))^M_.params(7);
oo_.steady_state( 9 ) = oo_.steady_state(8)+oo_.steady_state(1)*(M_.params(5)+M_.params(12));
oo_.steady_state( 10 ) = oo_.steady_state(8)+M_.params(5)*oo_.steady_state(1)-M_.params(3)*(oo_.steady_state(1)-M_.params(2)*oo_.steady_state(1))+oo_.steady_state(1)*M_.params(12);
oo_.steady_state( 7 ) = oo_.steady_state(2)*(1-1/M_.params(4));
oo_.steady_state( 11 ) = M_.params(3)*(oo_.steady_state(1)-M_.params(2)*oo_.steady_state(1));
oo_.steady_state( 12 ) = oo_.steady_state(10)-oo_.steady_state(1)*M_.params(12);
oo_.steady_state( 13 ) = oo_.steady_state(1)*oo_.steady_state(2)-M_.params(4)*oo_.steady_state(3);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(8))^2;
M_.Sigma_e(2, 2) = (M_.params(11))^2;
M_.Sigma_e(1, 2) = 0;
M_.Sigma_e(2, 1) = M_.Sigma_e(1, 2);
M_.sigma_e_is_diagonal = 0;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 25;
options_.order = 1;
options_.periods = 200;
var_list_ = char('k','b','q','z','s','i','y','c','dp','cons');
info = stoch_simul(var_list_);
save('KM_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('KM_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('KM_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('KM_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('KM_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('KM_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('KM_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
