function [ tsls_pow, post_haus_pow ] = local_power_mod_twosided(R, alpha_half, beta_1, beta_2, h_1, h_2, lambda, alpha_pos_adj, alpha_neg_adj)
%local_power: This function takes in the arguments R, alpha, beta_1,
%beta_2, h_1, h_2 and lambda, where all are constants, R is the number of
%psi values to be simulated, the psi_up etc. are pregenerated vectors of psi
%used for calculating a_adj, 
%and the rest of the values are as given in the
%instructions. The output is the local power of the 2sls test and the
%post-Hausman test at the given parameter values.

%% Throw errors if parameter values are invalid
% DO this later

%% Arguments for the generation of psi vectors
k_2 = 1;
s_vec = [1, zeros(1, k_2 - 1)]; % Vector of s_k2
p = 0;
g = 40;
alpha_g = alpha_half / 10;

%% 1.) 2SLS power
% Calculate the standard normal quantile
z = norminv(1 - alpha_half, 0, 1);
% Calculate psi values and add lambda * h_2
[psi_u, psi_v, psi_uv] = Psi_vec_gen(R, p, k_2);
psi_sum = abs(bsxfun(@plus, lambda .* h_2, psi_u));
% Generate the two-stage least squares power
tsls_pow = sum(psi_sum > z, 1) ./ length(psi_sum);

%% 2.) Hausman power

% 2.1.) Terms to the left of the inequality
% Generate eta using pregenerated psi from previous section
[eta_1, eta_2, eta_3, xi_2] = Eta_gen_comb_vec(h_1, h_2, psi_u, psi_v, psi_uv, s_vec, k_2);

% Generate eta_star_mod (equivalent to the term to the left of the
% inequality in the second equation on page 2 of the instructions)
eta_lhs = Eta_star_mod(eta_1, eta_2, eta_3, beta_1, lambda, h_2);

% Calculate alpha_adj
if nargin < 8
    [alpha_pos_adj, alpha_neg_adj] = level_adjust_twosided(alpha_half, beta_1, beta_2, h_2, R, g, alpha_g);
end

% 2.2.) Terms to the right of the inequality
% Compute I_beta_2(xi_2,h) for all xi:

%Find z_1-beta/2 in order to obtain the index function
norm_quant = norminv(1 - beta_2/2, 0, 1);

% Generate I_beta(xi)
[I_lo, I_hi] = indicfn(xi_2, norm_quant, h_2);

% Minimum and maximum possible values of xi_ind
I_lo_min = min(I_lo);
I_hi_max = max(I_hi);

% Calculate the quantiles of h1_local for the entire domain between
% I_lo_min and I_hi_max
h1_loc_val = linspace(I_lo_min, I_hi_max, 50);

% Calculate c(h1_loc, h_2, alpha_adj, beta_1), save and set aside
alpha_adj = [alpha_pos_adj, alpha_neg_adj];
[eta_quant_pos] = Eta_quantile_gen_vec(R, h1_loc_val, h_2, alpha_adj, beta_1);
[eta_quant_neg] = Eta_quantile_gen_vec_neg(R, h1_loc_val, h_2, alpha_adj, beta_1);


% Generate and save the sup_eta values
sup_eta_pos = zeros(length(xi_2),1)';
sup_eta_neg = zeros(length(xi_2),1)';

% These are the indices of the closest value of h1_loc_val to the
% boundaries of I_beta(xi) for each value of xi. They are used to find
% the sup of the quantiles in h_1 for each value of xi_2 in the sample
[~, xi_low_indices] = min(abs(bsxfun(@minus, I_lo', h1_loc_val)),[],2);
[~, xi_hi_indices] = min(abs(bsxfun(@minus, I_hi', h1_loc_val)),[],2);
xi_low_indices = cast(xi_low_indices, 'uint32');
xi_hi_indices = cast(xi_hi_indices, 'uint32');

% Using precalculated quantiles, for each configuration of xi_ind, 
% find the sup of c_(h1_loc, h_2) (1 - alpha_loc, 1 - beta) 
% within the precalculated boundaries of \mathfrac{h}_1
parfor xi_ind = 1:length(xi_2)
    % Take the values of eta_quant_loc that are between the above
    % boundaries
    eta_quant_lpos = eta_quant_pos(xi_low_indices(xi_ind): xi_hi_indices(xi_ind));
    eta_quant_lneg = eta_quant_neg(xi_low_indices(xi_ind): xi_hi_indices(xi_ind));
    
    %Find the sup of the eta over the local boundaries, save the
    %value
    sup_eta_pos(xi_ind) = max(eta_quant_lpos);
    sup_eta_neg(xi_ind) = max(eta_quant_lneg);
    %parfor_progress;
end

% 2.3.) Find asymptotic power of the local Hausman test
cond1 = bsxfun(@ge, eta_lhs', sup_eta_pos'); % lhs > sup c
cond2 = bsxfun(@ge, -1*eta_lhs', sup_eta_neg'); % -1*lhs > sup c
post_haus_pow = sum(bsxfun(@or, cond1, cond2))./size(eta_lhs,2); % probability of cond1 or cond2 being satisfied

end

