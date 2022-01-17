function [C_A, C_B, R_A, R_B, M, N] = cur_pair_deim(A, B, k, U, V, X)

%CUR_PAIR_DEIM  DEIM incurred CUR decomposition for matrix pairs
% [function [irow, iA, iB, M, N] = gcur_deim(A, B, k)
%
% CA = A(:,icol);  CB = B(:,icol);  RA = A(iA,:);  RB = B(iB,:);
%
% See also cur_deim
%
% (C) Perfect Gidisu, Michiel Hochstenbach 2021

if nargin < 3 || isempty(k), k = 2; end
[U, V, X, ~, ~] = gsvd(A,B,0);
U = fliplr(U(:,end-k+1:end));       % Select largest (C,S) pairs
V = fliplr(V(:,end-k+1:end));
X = fliplr(X(:,end-k+1:end));


for j = 1:k
  [~, icol(j)] = max(abs(X(:,j)));    % Iterative selection and projection
  [~, iA(j)] = max(abs(U(:,j)));
  [~, iB(j)] = max(abs(V(:,j)));
  X(:,j+1:k) = X(:,j+1:k) - X(:,1:j) * (X(icol,1:j) \ X(icol,j+1:k));
  U(:,j+1:k) = U(:,j+1:k) - U(:,1:j) * (U(iA,1:j) \ U(iA,j+1:k));
  V(:,j+1:k) = V(:,j+1:k) - V(:,1:j) * (V(iB,1:j) \ V(iB,j+1:k));
end
C_A = A(:,icol);  
C_B = B(:,icol);  
R_A = A(iA,:);  
R_B = B(iB,:);
M = C_A \ (A / R_A);
N = C_B \ (B / R_B);
