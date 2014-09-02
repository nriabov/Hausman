% Hausman Script (Based on Guggenberger 2010)
% June 26, 2010
% Nickolai M Riabov

% The purpose of this script is to simulate quantiles of the empirical
% distribution found in equation 20 of Guggenberger, 2010

cd '/users/nriabov/McCloskey2014-Sims/Hausman/'
addpath '/users/nriabov/McCloskey2014-Sims/Hausman/'
addpath '/users/nriabov/McCloskey2014-Sims/Hausman/Scripts'

% Enable parallel processings
distcomp.feature( 'LocalUseMpiexec', false );
matlabpool(feature('numCores')); %Use 8 nodes

%% Simulation (Level Adjustment)

% Parameters
alpha = 0.05; 
beta_1 = 0.1; % Parameter of eta_h_star
beta_2 = 0.05; % Parameter of the sup of the quantiles
h_2 = 2;
R = 500000; % Sample size
g = 40;
alpha_g = alpha / 10;

% Simulate the alpha_adj for the given parameter parameter values
alpha_adj = level_adjust( alpha, beta_1, beta_2, h_2, R, g, alpha_g);

% Save the data
savefile = 'adj_alpha.mat';
save(savefile, 'alpha_adj', 'alpha', 'beta_1', 'beta_2', 'h_2', ...
	'R', 'g', 'alpha_g');

		
%% Print the output
disp('**** Parameters ****')
disp('Given parameter values:')
disp(['alpha: ', num2str(alpha)])
disp(['beta_1: ', num2str(beta_1)])
disp(['beta_2: ', num2str(beta_2)])
disp(['h_2: ', num2str(h_2)])
disp(['R: ', num2str(R)])
disp(['g: ', num2str(g)])
disp(['alpha_g: ', num2str(alpha_g)])
disp('****Result****')
disp(['alpha_adj:', num2str(alpha_adj)])



%sprintf('**** Parameters ****')
%sprintf('Given parameter values:')
%sprintf('alpha: %1.2f', alpha);
%sprintf('beta_1: %1.2f', beta_1);
%sprintf('beta_2: %1.2f', beta_2);
%sprintf('h_2: %1.4f', h_2);
%sprintf('R: %0.1f', R);
%sprintf('g: %0.1f', g);
%sprintf('alpha_g: %1.4f', alpha_g);

%sprintf('**** Result ****');
%sprintf('alpha_adj: %1.4f', alpha_adj);


