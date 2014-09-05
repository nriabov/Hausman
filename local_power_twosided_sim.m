% Hausman Script (Based on Guggenberger 2010)
% June 26, 2010
% Nickolai M Riabov

% The purpose of this script is to simulate the local power of the Hausman
% test vs. the TSLS test under various parameter configurations

%cd 'Y:\McCloskey2014-Sims\Hausman'
%addpath 'Y:\McCloskey2014-Sims\Hausman'
%addpath 'Y:\McCloskey2014-Sims\Hausman\Scripts'
cd '/users/nriabov/McCloskey2014-Sims/Hausman/'
addpath '/users/nriabov/McCloskey2014-Sims/Hausman/'
addpath '/users/nriabov/McCloskey2014-Sims/Hausman/Scripts'

% Enable parallel processings
distcomp.feature( 'LocalUseMpiexec', false );
matlabpool(feature('numCores')); %Use all cores

%% Simulation (Level Adjustment)

% Parameters
alpha = 0.05;
alpha_half = alpha / 2;
num_h1 = 9; % Number of h_1 values in the range of h_1 to be simulated
num_lambda = 100; % Number of lambda values between 0 and lambda_max to be plotted for each h_2
R = 500000; % Sample size
%R = 2000000; % Sample size for small values of lambda
k_2 = 1;
p = 0;
alpha_g = alpha_half / 10;
g = 40;

h2_vec = [0.001, 0.017, 0.028, 0.1, 0.5, 1, 2, 10]; % h_2 values to be tested
h_1range = [-3000, 4000; % Range of h_1 values
            -180, 240;
            -90, 140;
            -30, 35;
            -5, 12;
            -5, 10;
            -4, 8;
            -6, 7; 
            -6, 6];
lambda_max = [4000, 200, 150, 45, 9, 5, 3, 0.45, 0.2];
beta1_vec = 0.1:0.1:0.9;
beta2_vec = 0.1:0.1:0.9; % CHANGE SPACING BACK TO 0.1!!!

% Save file output
tsls_vec_perm = zeros(num_h1, num_lambda, length(h2_vec), length(beta1_vec));
post_haus_vec_perm = zeros(num_h1, num_lambda, length(h2_vec), length(beta1_vec));

% Simulate R_loc draws of Psi for the sake of generating the quantiles
R_loc = R;
[psi_up_1, psi_vp_1, psi_uvp_1] = Psi_vec_gen(R_loc, p, k_2); % Make the Psi vectors in eq 15

% Save optimal betas and optimal alpha_adj values
beta1_optim = zeros(length(h2_vec), length(beta1_vec));
alpha_adj_optim_pos = beta1_optim;
alpha_adj_optim_neg = beta1_optim;

mkdir('/users/nriabov/McCloskey2014-Sims/Hausman/beta1search-twosided/');
cd '/users/nriabov/McCloskey2014-Sims/Hausman/beta1search-twosided/'
t_overall = tic;
simnum = 1;
% Simulate local power for different values of beta_1
for k = 1:length(beta1_vec)
%for k = length(beta1_vec):length(beta1_vec)
    beta_1 = beta1_vec(k);
    loc_dir = ['beta1-', num2str(beta_1)];
    % Simulate the local power for different values of h_1 and h_2
    %for j = 3:3
    %for j = length(h2_vec)-2:length(h2_vec)
    for j = 1:length(h2_vec)
        t_loc = tic;
        h_2 = h2_vec(j); % Strength of the instrument

        % Calculate the level-adjusted alpha for all the possible beta_2
        % associated with beta_1
        [alpha_pos_vec, alpha_neg_vec] = level_adjust_mod_twosided(alpha_half, beta_1, beta2_vec, h_2, R, g, alpha_g, psi_up_1, psi_vp_1, psi_uvp_1);
        % Find the alpha_vec closest to 0.045 and use it. If there's
        % multiple, just use the first one
        alpha_vec = alpha_pos_vec + alpha_neg_vec;
        [~, closest_ind] = min(abs(alpha_vec - 0.045));
        alpha_pos_adeq = alpha_pos_vec(closest_ind);
        alpha_neg_adeq = alpha_neg_vec(closest_ind);
        alpha_pos_adj = alpha_pos_adeq(1);
        alpha_neg_adj = alpha_neg_adeq(1);
        
        % Save the values of alpha_adj and beta_optimal for further analysis
        beta2_optim_loc = beta2_vec(closest_ind(1));
        beta2_optim(j, k) = beta2_optim_loc;
        alpha_adj_optim_pos(j, k) = alpha_pos_adj;
        alpha_adj_optim_neg(j, k) = alpha_neg_adj;


        % Generate the lambda
        lambda = linspace(-lambda_max(j), lambda_max(j), num_lambda);

        % Make matrices to store results
        tsls_vec = zeros(num_h1, num_lambda);
        post_haus_vec = zeros(num_h1, num_lambda);

        % Set boundaries of h1_left and h1_right
        h1_left = linspace(h_1range(j,1), 0, 4);
        h1_right = linspace(0, h_1range(j,2), 6);
        h_1 = [h1_left(1:end-1), h1_right];

        % Calculate the local power
        parfor i = 1:length(h_1)
            [ tsls_pow, post_haus_pow ] = local_power_mod_twosided(R, alpha_half, beta_1, beta2_optim_loc, h_1(i), h_2, lambda, alpha_pos_adj, alpha_neg_adj);
            tsls_vec(i,:) = tsls_pow;
            post_haus_vec(i,:) = post_haus_pow;
        end
        
        % Save the TSLS power values and the Post-Hausman power
        tsls_vec_perm(:,:,j,k) = tsls_vec;
        post_haus_vec_perm(:,:,j,k) = post_haus_vec;
        
        % Make the plots
        for a = 1:length(h_1)
            subplot(3, 3, a);
            plot(lambda, tsls_vec(a, :), '-r', lambda, post_haus_vec(a, :), '-b', 'LineWidth', 2);
            xlim([min(lambda), max(lambda)]);
            title_vec = [' $h_1 = ', num2str(h_1(a)), '$'];
            title(title_vec, 'Interpreter', 'latex', 'Fontsize', 13)
            if a == length(h_1)
                hL = legend('TSLS', 'Hausman', 'Location', 'Southeast', 'Orientation', 'horizontal');
                xlabel('$\lambda$', 'Interpreter', 'latex', 'Fontsize', 13);
            end
        end

        name = 'Power Curves, ';
        name = [name, '$\alpha = ', num2str(alpha), ...
            '$, $\beta_1 = ', num2str(beta_1), '$, $\beta_2 = ', num2str(beta2_optim_loc), ...
            '$, $ h_2 = ', num2str(h_2), '$'];
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off', ...
            'Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1, name ,'HorizontalAlignment', ...
            'center','VerticalAlignment', 'top', 'Interpreter', 'latex', 'Fontsize', 13)
		
        legend boxoff;
		newPosition = [0.5, 0.019, 0.06, 0.04];
		set(hL, 'Position', newPosition, 'Units', 'normalized');

        mkdir(loc_dir);
        cd(loc_dir);

        % Output the figure to pdf
        saveas(gcf, ['PowerCurves-h2-', num2str(j), '.m']); 
        %PrintFig(['PowerCurves-h2-', num2str(h2_vec(j))], gcf);
        cd ..
        t_loc = toc(t_loc);
        disp(['Simulation number: ',  num2str(simnum), ' / ', num2str(length(h2_vec) * length(beta1_vec))])
        simnum = simnum + 1;
        disp(['This simulation took: ', num2str(t_loc), ' secs'])
    end
end

t_overall = toc(t_overall);
disp(['Total Runtime: ', num2str(t_overall./60), ' mins']);

savefile = 'power_out.mat';
save(savefile, 'beta2_optim', 'alpha_adj_optim_pos', ...
    'alpha_adj_optim_neg', 'tsls_vec_perm', 'post_haus_vec_perm');