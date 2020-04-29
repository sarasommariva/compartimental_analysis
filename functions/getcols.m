function Acols = getcols(A,cols)
% Acols = GETCOLS(A,cols) get the columns of matrix A whose indices are given
% in the vector cols.
Acols=A(:,cols);