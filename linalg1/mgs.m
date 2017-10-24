function [A, R] = mgs(A)

% This function performs the reduced QR decomposition using the modified
% Gram-Schmidt algorithm.  The input A is overwritten to give the matrix Q.
% Input: A: complex mxn matrix with m>=n
% Outputs: A: mxn matrix containing the orthonormal vectors spaning range(A)
%          R: nxn upper-triangular matrix containing the coefficients to write
%             columns of A as a linear combination of the columns of Q.


[m, n] = size(A);
R = zeros(n,n);

for i = 1:n
    
    rii = norm(A(:,i)); %compute diagonal entry of R from ith column of A
    q = A(:,i)/rii; %compute the ith column of Q
    R(i,i) = rii; %store the entry ii in R
    A(:,i) = q; %store the ith column of Q as the ith column of A
    
    for j = i+1:n
        
        V = A(:,j); %set V as the jth column of A
        rij = q'*V; %take the inner product of q and V
        R(i,j) = rij; %store inner product in R
        V = V - rij*q; %subtract off component V in direction q
        A(:,j) = V; %store V as the jth column of A
        
    end
end