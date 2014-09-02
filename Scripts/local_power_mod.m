function [ tsls_pow, post_haus_pow ] = local_power_mod(R, alpha, beta_1, beta_2, h_1, h_2, lambda, psi_up_1, psi_vp_1, psi_uvp_1 )
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
alpha_g = alpha / 10;

%% 1.) 2SLS power
% Calculate the standard normal quantile
z = norminv(1 - alpha, 0, 1);
% Calculate psi values and add lambda * h_2
[psi_u, psi_v, psi_uv] = Psi_vec_gen(R, p, k_2);
%psi_sum = lambda .* h_2 + psi_u;
psi_sum = bsxfun(@plus, lambda .* h_2, psi_u);
% Generate output
tsls_pow = sum(psi_sum > z, 1) ./ length(psi_sum);

%% 2.) Hausman power

% 2.1.) Terms to the left of the inequality
% Generate eta using pregenerated psi from previous section
[eta_1, eta_2, eta_3, xi_2] = Eta_gen_comb_vec(h_1, h_2, psi_u, psi_v, psi_uv, s_vec, k_2);

% Generate eta_star_mod (equivalent to the term to the left of the
% inequality in the second equation on page 2 of the instructions)
eta_lhs = Eta_star_mod(eta_1, eta_2, eta_3, beta_1, lambda, h_2);

% Calculate alpha_adj

alpha_adj = level_adjust(alpha, beta_1, beta_2, h_2, R, g, alpha_g);
%alpha_adj = level_adjust_mod(alpha, beta_1, beta_2, h_2, R, g, alpha_g, psi_up_1, psi_vp_1, psi_uvp_1);
%alpha_adj = 0.04;

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
%h1_loc_val = I_lo_min:1:I_hi_max; % Make the grid finer later for greater accuracy
h1_loc_val = linspace(I_lo_min, I_hi_max, 50);

% Calculate c(h1_loc, h_2, alpha_adj, beta_1), save and set aside
eta_quant_mat = Eta_quantile_gen_vec_psigiven(psi_up_1, psi_vp_1, ...
        psi_uvp_1, k_2, h1_loc_val, h_2, s_vec, alpha_adj, beta_1 ); 
%eta_quant_mat = Eta_quantile_gen_vec(R, h1_loc_val, h_2, alpha_adj, beta_1);

% Generate and save the sup_eta values
sup_eta = zeros(length(xi_2),1)';

% Using precalculated quantiles, for each configuration of xi_ind, 
% find the index of the closest h1_loc_val to the boundaries of 
% I_beta(xi) and find the sup of c_(h1_loc, h_2) (1 - alpha_loc, 1 - beta) 
% within those boundaries

parfor xi_ind = 1:length(xi_2)
    [~, xi_low_ind] = min(abs(I_lo(xi_ind) - h1_loc_val));
    [~, xi_hi_ind] = min(abs(I_hi(xi_ind) - h1_loc_val));

    % Take the values of eta_quant_loc that are between the above
    % boundaries
    eta_quant_rel = eta_quant_mat(xi_low_ind: xi_hi_ind);

    %Find the sup of the eta over the local boundaries, save the
    %value
    sup_eta(xi_ind) = max(eta_quant_rel);
    %parfor_progress;
end

% 2.3.) Find asymptotic power of the local Hausman test
post_haus_pow = sum(bsxfun(@ge, eta_lhs', sup_eta'))./size(eta_lhs,2);
%post_haus_pow = sum(eta_lhs > sup_eta) / length(eta_lhs);

end

