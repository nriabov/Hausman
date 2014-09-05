function [ eta_quant_mat ] = Eta_quantile_gen_vec_neg(N, h_1, h_2, alpha, beta )
%Eta_quantile_gen: This code takes in the sample size per simulation N
% the number of instruments k_2, the vector of
% h_1 values h_1, the scalar of h_2 values h_2, the vector of s_k2 s_vec, 
% the desired quantile of the distribution of eta_star alpha, and the 
% quantile of the chi-square distribution beta. It outputs an |h_1| *
% |alpha| matrix of quantiles of -1*eta_star computed for a given vector of 
% h_1 and other parameter values
% [ eta_quant_mat ] = Eta_quantile_gen(N, k_2, h_1, h_2, s_vec, alpha, beta )


% p value used in the generation of the Psi vector.
p = 0; 
k_2 = 1;
s_vec = [1, zeros(1, k_2 - 1)]; % Vector of s_k2

% Make the Psi vectors in eq 15
[psi_up, psi_vp, psi_uvp] = Psi_vec_gen(N, p, k_2);

% Generate eta_{1:3},h as in equation 17
[eta_1, eta_2, eta_3] =  Eta_gen_comb_vec(h_1, h_2, psi_up, psi_vp, psi_uvp, s_vec, k_2);

% Generate eta_h_star, as in equation 20
eta_star = -1*Eta_star_gen_vec(eta_1, eta_2, eta_3, beta);

% Generate the 1 - alpha quantile of the eta_h distribution for each h_1
%   The output is a |h_1| * |alpha| vector.
eta_quant_mat = quantile(eta_star, 1 - alpha, 2);

end

