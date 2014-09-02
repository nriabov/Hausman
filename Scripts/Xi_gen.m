function [ xi_1, xi_2, xi_3, xi_4 ] = Xi_gen( h_vec, psi_up, psi_vp, psi_uvp, s_vec, k_2 )
%Xi_gen: This function takes in a 2-dimensional vector of h values, the 
% three matrices of psi values made by Psi_vec_gen, and a k_2-dimensional 
% vector of s_k2 that has a norm equal to 1. It outputs the xi_h vector
% found in equation 16 of Guggenberger(2010)

% Fail conditions:
%   h_vec isn't 2-dimensional
%   Psi_up or Psi_vp aren't N*k_2 dimensional
%   s_vec doesn't have a norm of 1

if length(h_vec) ~= 2
    error('Please provide an h vector of length 2');
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

% If all conditions met
%   If nargout == 1, make Xi_1 only
%   If nargout == 2, make Xi_2 and Xi_2
%   If nargout == 3, ...

h_1 = h_vec(1);
h_2 = h_vec(2);

xi_1 = h_2 * s_vec * psi_up';
if nargout >= 2
    xi_2 = xi_1 + psi_uvp' + h_1;
    if nargout >= 3
        xi_3 = h_2^2;
        if nargout == 4
            xi_4 = xi_3 + 1; 
        end
    end
end



end

