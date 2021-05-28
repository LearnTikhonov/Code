function [tr] = trace_op(A,N)
% Computes the trace of a linear operator A, defined on R^N (in a matrix-free version)
% By definition, trace(A) = sum_{i=1}^N <A(e_i),e_i>, 
% where e_i are the vectors of the canonical basis.

tr = 0;
for i = 1:N
    e_i = zeros(N,1);
    e_i(i) = 1;
    A_e_i = A(e_i);
    tr = tr + A_e_i(i);  % note: A_e_i(i) = < A(e_i),e_i>
end
end