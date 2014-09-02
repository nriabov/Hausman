function [ eta_1, eta_2, eta_3 ] = Eta_gen(s_vec, psi_up, h_vec, xi_2 )
%Eta_gen: This function takes in the k_2 dimensional s_vec, the matrix
%psi_up generated by Psi_vec_gen, the vector h_vec with the limits of
%gamma_n_1 and gamma_n_2, and the vector xi_2 generated by Xi_gen. It
%outputs the three nontrivial random variables eta_h given in equation 17
% [ eta_1, eta_2, eta_3 ] = Eta_gen(s_vec, psi_up, h_vec, xi_2 )

% Generate h_2 and do precalculations
h_2 = h_vec(2);
h_sum = 1 + h_2^2;

% If all conditions met
%   If nargout == 1, make eta_1 only
%   If nargout == 2, make eta_1 and eta_2
%   If nargout == 3, ...


eta_1 = s_vec * psi_up';
if nargout >= 2
    eta_2 = h_sum^(-1/2) .* xi_2;
    if nargout == 3
        eta_3 = h_sum .* (eta_1 - h_2 * h_sum^(-1) * xi_2);
    end
end



end

