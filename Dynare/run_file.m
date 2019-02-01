clear
close all

addpath c:\dynare\4.5.4\matlab

%Run Dynare
dynare BGP_Jan2019
asd


zeta_shock = [rp_eps_zeta e_eps_zeta x_eps_zeta];

thet_shock = [rp_eps_thet e_eps_thet x_eps_thet];

IRFs_Dynare = [zeta_shock thet_shock];

save IRFs_Dynare IRFs_Dynare
