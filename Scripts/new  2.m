% Hausman Script (Based on Guggenberger 2010)
% June 26, 2010
% Nickolai M Riabov

% The purpose of this script is to simulate quantiles of the empirical
% distribution found in equation 20 of Guggenberger, 2010

cd 'C:\Users\nriabov\Dropbox\McCloskey RA\2014\Hausman\Code'
addpath 'C:\Users\nriabov\Dropbox\McCloskey RA\2014\Hausman\Code'
addpath 'C:\Users\nriabov\Dropbox\McCloskey RA\2014\Hausman\Code\Scripts'

% Enable parallel processings
distcomp.feature( 'LocalUseMpiexec', false );
matlabpool(feature('numCores')); %Use 8 nodes

%% Simulation (Multiple parameter values, single quantile of eta)

% Parameters
alpha = 0.05; 
beta = 0.05;
h_2 = 2;
R = 350000; % Sample size
g = 30;
alpha_g = alpha / 10;



