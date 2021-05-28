function [w] = haart_vec(x,lev)
% Implements the Haar wavelet transform for 1d signal
% It relies on the haart Matlab routine, but produces as an output a vector
% of the same size of the signal.
% It works only for signals whose size is a power of 2

% x: original signal
% lev: maximum level of resolution (default: log2(length(w)))
% w: coefficients of the wavelet transform, ordered as follows: first the
% detail coefficients of the largest levels, until the specified level
% lev, then the approximation coefficients 

N = length(x);
Mlev = round(log(N)/log(2));
if N ~= 2^Mlev
    error('If using wavelets, the signal size must be a power of 2')
end
if nargin == 1
    lev = Mlev;
end

[a,d] = haart(x,lev);
w = zeros(size(x));
index = 1;
for ll = 1:lev
    w(index:index+2^(Mlev-ll)-1)=d{ll};
    index = index+2^(Mlev-ll);
end
w(index:end) = a;

end