function [ eta_quant_mat ] = Eta_quantile_gen(N, k_2, h_vec, s_vec, alpha, beta )
%Eta_quantile_gen: This code takes in the sample size per simulation N
% the number of instruments k_2, the vector of
% h values h_vec, the vector of s_k2 s_vec, the desired quantile of the
% distribution of eta_star alpha, and the quantile of the chi-square
% distribution beta.
% [ eta_quant_mat ] = Eta_quantile_gen(N, k_2, h_vec, s_vec, alpha, beta )

% p value used in the generation of the Psi vector.
p = 0; 

% Make the Psi vectors in eq 15
[psi_up, psi_vp, psi_uvp] = Psi_vec_gen(N, p, k_2);

% Generate eta_{1:3},h as in equation 17
[eta_1, eta_2, eta_3] =  Eta_gen_comb(h_vec, psi_up, psi_vp, psi_uvp, s_vec, k_2);

% Generate eta_h_star, as in equation 20
eta_star = Eta_star_gen( eta_1, eta_2, eta_3, beta);

% Generate the 1 - alpha quantile of the eta_h distribution
eta_quant_mat = quantile(eta_star, 1 - alpha);


end

