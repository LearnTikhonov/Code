function [x] = haart_vecT(w,lev)
% Implements the adjoint Haar wavelet transform for 1d signal
% It relies on the ihaart Matlab routine, but accepts as input a vector of
% the same size of the signal
% Works only for signals whose size is a power of 2

% w: coefficients of the wavelet transform, ordered as in haart_vec
% lev: maximum level of resolution (default: log2(length(w)))
% x: reconstructed signal

N = length(w);
Mlev = round(log(N)/log(2));
if N ~= 2^Mlev
    error('If using wavelets, the signal size must be a power of 2')
end

if nargin == 1
    lev = Mlev;
end

d = cell(lev,1);
index = 1;
for ll = 1:lev
    d{ll} = w(index:index+2^(Mlev-ll)-1);
    index = index+2^(Mlev-ll);
end
a = w(index:end);
x= ihaart(a,d);

end