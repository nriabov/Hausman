function [ psi_up, psi_vp, psi_uvp ] = Psi_vec_gen( N, p, k_2)
%Psi_vec_gen: This function accepts as arguments a sample size N, a
%correlation coefficient p, and a dimension parameter k_2 and 
%generates N realizations of a multivariate normal random variable of the 
%structure found in eq. 15 of Guggenberger
%
% [ psi_up, psi_vp, psi_uvp ] = Psi_vec_gen( N, p, k_2)

% If p = 0, generate an identity matrix of variance
if p == 0
    cov_mat = eye(2*k_2 + 1);
else % Otherwise, calculate the covariances
    V_p = [1, p; p, 1];
    ul_sigma = kron(V_p, eye(k_2));
    num_zero = size(ul_sigma, 1);
    cov_mat = [ul_sigma, zeros(num_zero,1); zeros(1, num_zero), 1 + p^2];
    % Take the square root of the covariance matrix by Cholesky
    % factorization
    cov_mat = chol(cov_mat);
end

% Generate Standard Normal random variables, multiply them by the square
% root of the covariance matrix
psi_mat = randn(N, 2*k_2 + 1) * cov_mat;

psi_up = psi_mat(:,1:k_2);
psi_vp = psi_mat(:,k_2 + 1:2*k_2);
psi_uvp = psi_mat(:, end);

end

