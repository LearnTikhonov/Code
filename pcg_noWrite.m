function [x] = pcg_noWrite(A,b)
% Suppresses output when running pcg
    [x,~,~,~] = pcg(A,b);
end