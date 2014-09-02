function [ alpha_adj, h1_min_val ] = level_adjust( alpha, beta_1, beta_2, h_2, R, g, alpha_g)
%level_adjust: This function takes in an alpha value between 0 and 1, a
% beta_1 parameter of eta_h_star between 0 and 1, a
% beta_2 parameter of the sup of the quantiles between 0 and 1, a positive value of h_2, 
% a number of simulations R, a number of draws from eta_h alpha_g, a number of
% iterations of alpha g, and a minimal value of alpha alpha_g. First, it
% generates a grid of candidate h_1 values. Then, based on this grid, it
% generates P(eta_h_star > sup_h1 \in I_b(xi) c(h_1, h_2, alpha, beta)) for
% the vector of h_1 values and alpha_j.

% Throw error if:
%   alpha not between 0 and 1
%   beta not between 0 and 1
%   h_2 less than 0
%   g is unreasonable
%   alpha_g is greater than alpha

if alpha < 0 || alpha > 1
    error('alpha must be between 0 and 1');
end
if beta_1 < 0 || beta_1 > 1 || beta_2 < 0 || beta_2 > 1
    error('beta must be between 0 and 1');
end
if h_2 < 0
    error('h_2 must be larger than 0');
end
if g < 0 || isnan(g) || isinf(g)
    error('Please use a valid value of g');
end
if alpha_g >= alpha
   error('alpha_g must be strictly less than alpha');
end


% 0.) Set all stationary parameter values (p, s_vec, k_2, etc...)
% Parameters
p = 0; % Correlation (Must be 0 for xi and eta to be computed according to equation 16 and 17)
k_2 = 1; % Dimension of instrument vector Z
s_vec = [1, zeros(1, k_2 - 1)]; % Vector of s_k2

% Divide (0, \alpha] into [0, alpha_g, alpha_g-1, ... , alpha], save and
% set aside
alpha_space = (alpha - alpha_g) / g;
alpha_vec = alpha : -alpha_space : alpha_g;

% 1.) Generate h_1

% If the parameters have nice values, then generate a grid of length 40
% which captures the point of the quantile of eta where it is not flat,
% along with some spots where it is flat around the sides

tic
if (beta_1 <= 0.1 && alpha_g > 0.001)
   if h_2 <= 0.75 % If h_2 is small, start with a wide range of potential h1 values
       h1_cand = -200:20:(10/h_2);
   elseif h_2 > 0.75 && h_2 < 5 % If h_2 is larger, start with a narrower range of h1 values
       h1_cand = -15:2:30/h_2;
   elseif h_2 >= 5 && h_2 <= 10 % Note: Don't make h_2 larger than 10!!!!
       h1_cand = -10:0.5:100/h_2;
   else 
       h1_cand = -10:0.5:400/h_2;
   end
   
   % Generate the quantiles of eta
   eta_quant_cand = Eta_quantile_gen_vec(R, h1_cand, h_2, [min(alpha_vec), max(alpha_vec)], beta_1);
   
   % Find the largest possible range for the unique values
   [~, uniq_ind_1, ~] = unique(eta_quant_cand(:,1));
   [~, uniq_ind_2, ~] = unique(eta_quant_cand(:,2));
   
   % Choose the indices of the h_vals to correspond to whichever one has
   % the longer range
   if length(uniq_ind_1) > length(uniq_ind_2)
       uniq_ind = sort(uniq_ind_1);
   else
       uniq_ind = sort(uniq_ind_2);
   end
   
   % Pad the end values of uniq_ind
   uniq_ind = uniq_ind(2:end);
   
   if(h_2 < 0.5) % Get a bigger chunk of the flat section for smaller values of h_2
       min_ind = min(uniq_ind) - 3;
       max_ind = max(uniq_ind) + 3;
   else
       min_ind = min(uniq_ind) - 1;
       max_ind = max(uniq_ind) + 1;
   end
   
   % Generate a finer grid of h1_vals over the discovered relevant indices,
   % including the padded values
   h1_vals = linspace(h1_cand(min_ind), h1_cand(max_ind), 30);
   
   % Test code
%    eta_quant_test = Eta_quantile_gen_vec(R, h1_vals, h_2, [min(alpha_vec), max(alpha_vec)], beta_1);
%    subplot(1, 2, 1);
%    plot(h1_vals, eta_quant_test(:,1))
%    subplot(1, 2, 2);
%    plot(h1_vals, eta_quant_test(:,2))
%    
else % If parameters are not nice, set h_1 by hand
    h1_vals = -100:20:300; % CHANGE ME LATER!!!
end

% 2.) For each h_1, 
% Simulate R draws of (eta_h_star, xi_2,h)
[psi_up, psi_vp, psi_uvp] = Psi_vec_gen(R, p, k_2); % Make the Psi vectors in eq 15

% Simulate R_loc draws of Psi for the sake of generating the quantiles
R_loc = 500000;
[psi_up_1, psi_vp_1, psi_uvp_1] = Psi_vec_gen(R_loc, p, k_2); % Make the Psi vectors in eq 15

% generate output vector for the entire loop
alpha_res_vec = zeros(length(h1_vals), 1);

tic
for h1_ind = 1:length(h1_vals)
    h1_loc = h1_vals(h1_ind);
    h_vec = [h1_loc, h_2]; % Make the local h_vec
    % Generate eta_{1:3},h as in equation 17
    [eta_1, eta_2, eta_3, xi_2] =  Eta_gen_comb(h_vec, psi_up, psi_vp, psi_uvp, s_vec, k_2);
    % Generate eta_h_star, as in equation 20
    eta_star = Eta_star_gen( eta_1, eta_2, eta_3, beta_1);
    
    %Find z_1-beta/2 in order to obtain the index function
    norm_quant = norminv(1 - beta_2/2, 0, 1);

    % Generate I_beta(xi)
    [I_lo, I_hi] = indicfn(xi_2, norm_quant, h_2);

    % Minimum and maximum possible values of xi_ind
    I_lo_min = min(I_lo);
    I_hi_max = max(I_hi);

    % Generate the grid of h1_local for the entire domain between
    % I_lo_min and I_hi_max
    h1_loc_val = linspace(I_lo_min, I_hi_max, 50);
    % Calculate and save quantiles c_(h1_loc, h_2, alpha_vec, beta)
    eta_quant_mat = Eta_quantile_gen_vec_psigiven(psi_up_1, psi_vp_1, ...
        psi_uvp_1, k_2, h1_loc_val, h_2, s_vec, alpha_vec, beta_1 ); 
    
    % Do parts 2.c through 2.e of the instruction

    % These are the indices of the closest value of h1_loc_val to the
    % boundaries of I_beta(xi) for each value of xi. They are used to find
    % the sup of the quantiles in h_1 for each value of xi_2 in the sample
    [~, xi_low_indices] = min(abs(bsxfun(@minus, I_lo', h1_loc_val)),[],2);
    [~, xi_hi_indices] = min(abs(bsxfun(@minus, I_hi', h1_loc_val)),[],2);
    xi_low_indices = cast(xi_low_indices, 'uint32');
    xi_hi_indices = cast(xi_hi_indices, 'uint32');
    
    % Initialize while loop conditions
    probvals = 1;
    j = 1;
    while probvals > alpha && j <= g % Stop if either probvals is <= alpha or you reach the end of alpha_vec
        % Save the sup_eta values
        sup_eta = zeros(length(xi_2),1);
        
        % Take only the quantiles of eta relevant for alpha_j
        eta_quant_loc = eta_quant_mat(:,j);
        
        % Using precalculated quantiles, for each value of xi, 
        % find the index of the closest h1_loc_val to the boundaries of 
        % I_beta(xi). Then, find the sup of c_(h1_loc, h_2) (1 - alpha_loc, 1 - beta) 
        % within those boundaries
        % New code: boundary indices xi_low_indices and xi_hi_indices are 
        % precomputed before the while loop.
        
        parfor xi_ind = 1:length(xi_2)
            % Take the values of eta_quant_loc that are between the above
            % boundaries
            eta_quant_rel = eta_quant_loc(xi_low_indices(xi_ind): xi_hi_indices(xi_ind));
            
            %Find the sup of the eta over the local boundaries, save the
            %value
            sup_eta(xi_ind) = max(eta_quant_rel);
        end
		
        % Calculate the probability values
        probvals = sum(eta_star > sup_eta')/length(eta_star);
        
        % Iterate j
        j = j + 1;
    end
    
    % If probvals <= alpha and j < g, set alpha_res_vec(h_1) = alpha_j
    if j < g
		alpha_loc = alpha_vec(j-1); % j-1 because last command in while loop is j = j + 1;
        alpha_res_vec(h1_ind) = alpha_loc; 
        % Display alpha value 
        disp(['alpha_h: ', num2str(alpha_loc)])
    end
    %parfor_progress;
    disp(['Simulation number: ',  num2str(h1_ind), ' / ', num2str(length(h1_vals))])
end
toc

% Output alpha_adj from part 3 of the instructions and possibly the 
% corresponding minimum value of h_1
[alpha_adj, h1_min] = min(alpha_res_vec);

if nargout > 1
	h1_min_val = h1_vals(h1_min);
end
	
end

