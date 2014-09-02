function [ eta_1, eta_2, eta_3, xi_2 ] = Eta_gen_comb_vec( h_1, h_2, psi_up, psi_vp, psi_uvp, s_vec, k_2 )
%Eta_gen: This function takes in 2-dimensional vector of h values, the 
% three matrices of psi values made by Psi_vec_gen, and a k_2-dimensional 
% vector of s_k2 that has a norm equal to 1. It
%outputs the three nontrivial random variables eta_h given in equation 17
% [ eta_1, eta_2, eta_3 ] = Eta_gen(s_vec, psi_up, h_vec, xi_2 )

% Fail conditions:
%   h_vec isn't 2-dimensional
%   Psi_up or Psi_vp aren't N*k_2 dimensional
%   s_vec doesn't have a norm of 1

if length(h_2) ~= 1
    error('Please provide an h_2 vector of length 1');
end
if size(psi_up, 2) ~= k_2
    error('psi_up is supposed to have k_2 columns');
end
if size(psi_vp, 2) ~= k_2
    error('psi_vp is supposed to have k_2 columns');
end
if size(psi_uvp, 2) ~= 1
    error('psi_vp is supposed to have k_2 columns');
end
if norm(s_vec) ~= 1
    error('The norm of s_vec is supposed to be 1');
end
if size(s_vec,2) ~= k_2
    error('The length of s_vec is supposed to be k_2 with k_2 columns');
end

%h_1 = h_vec(1);
h_sum = 1 + h_2^2;

% If all conditions met
%   If nargout == 1, make eta_1 only
%   If nargout == 2, make eta_1 and eta_2
%   If nargout == 3, ...

% Generate the xi values from equation 16 that will be used in equation 17
xi_1 = h_2 * s_vec * psi_up';
xi_2 = bsxfun(@plus, xi_1 + psi_uvp', h_1');

eta_1 = s_vec * psi_up';
if nargout >= 2
    eta_2 = h_sum^(-1/2) .* xi_2;
    if nargout >= 3
        eta_3 = h_sum .* (bsxfun(@minus, eta_1,  h_2 * h_sum^(-1) * xi_2)).^2;
    end
end



end

