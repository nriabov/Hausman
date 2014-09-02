function [ I_lo, I_hi ] = indicfn( x, z, h_2 )
%indicfn: This function takes in x, z, and h_2 as given in I_b(x) in the
%instructions, and outputs the boundaries of the interval I_b(x)

in_sum = (1 + h_2^2)^(1/2);
I_lo = bsxfun(@minus, x', in_sum .* z)';
I_hi = bsxfun(@plus , x', in_sum .* z)';

end

