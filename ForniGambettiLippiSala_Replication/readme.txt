Readme file

This file explains how to reproduce the figures in ``No News in Business Cycles''
by Mario Forni, Luca Gambetti, Marco Lippi and Luca Sala
AEJMacro-2015-0359 - No News in Business Cycles


There are 2 folders (DSGE and kilian_bootstrap), put them in the Matlab path.

%%%%%%%%%%%%%%%%%%%
The main replication file is:  replication_Noisy_News_Business_Cycles.m
By running this, all the figures in the papers and all the tables are generated.
%%%%%%%%%%%%%%%%%%%

MainProgAEJRevision_1.m			: replicates the empirical results in the paper
MainProgAEJRevision_2.m			: replicates the empirical results in the paper

FAVARNewsNoiseChol.m			: estimates a VAR model, identifies a noise and a 						news shock (FAVARCholImp.m) and computes 							confidence bands (FAVARCholBoot.m).

DoFiguresBC_luca2_rev.m			: draw figures
DoFiguresBC_luca2_rev2.m			: draw figures

ReadDataNews_EJ_data.m			: loads data
dataset_EJ_data.xlsx			: data in excel format
population.xlsx				: population data in excel format
Robustness_mario3.m				: robustness analysis
Robustness_luca2.m				: robustness analysis
ortotest.m					: computes the orthogonality test in Forni- 							Gambetti, JME, 2014.

FAVARCholIdent.m				: identifies reduced form impulse reponses							 with a Cholesky order
FAVARCholBoot.m				: computes bootstrapped impulse responses from a 						Cholesky identification
FAVARNewsNoiseMax.m				: estimates a VAR model, identifies a noise and a 						news shock and computes confidence bands.
FAVARCholBootMax.m				: bootstrap analysis in the estimated VARs using 						the maximal impulse response  identification 							scheme

myols.m					: OLS estimation
FAVARNewsImp.m				: computes impulse responses using the maximal 						impulse response identification scheme
LongRunEffectZ.m				: function to be maximized to obtain the maximal 						impulse response identification scheme
ComputeIrfOneShock.m			: takes an identified shock and computes its 							impulse responses
GenerateNewSeries.m				: uses the estimated VAR generates new series
VarParameters.m				: takes VAR estimates and builds a companion form
aicbic.m					: computes information criteria
cfilter.m					: computes band passed time series
FAVARRaw.m					: computes impulse responses
FAVARCholImp.m				: computes impulse responses using a Cholesky 						identification scheme
VAR_str.m					: constructs the regressors for a VAR(k)
principalcomponents.m			: computes principal components 
polynomialmatricesproduct.m		: multiplies matrix polynomials
woldimpulse.m					: takes a VAR and computes the MA representation 
companion.m					: computes a companion form from a VAR
invertpolynomialmatrix.m			: inverts a matrix polynomial
myvar.m					: estimates a VAR(k)
center.m					: subtracts the mean for time series
standardize.m					: standardizes time series


DSGE folder

FGLS_experiments_news_and_noise_confidence_8shocks.m		: replicates the DSGE results in the 								paper
FGLS_4lags_generate_restricted_8shocks.m			: generates data from the DSGE 									model
fai_fill.m								: plots results


kilian_bootstrap folder 

These files implements Kilian's IRF bias correction 



