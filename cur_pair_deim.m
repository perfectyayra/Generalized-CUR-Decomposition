function [icol, iA, iB, M, N] = cur_pair_deim(A, B, k, U, V, X)

%CUR_PAIR_DEIM  DEIM incurred CUR decomposition for matrix pairs
% [function [irow, iA, iB, M, N] = cur_pair_deim(A, B, k, U, V, X)
%
% CA = A(:,icol);  CB = B(:,icol);  RA = A(iA,:);  RB = B(iB,:)
%
% See also CUR_DEIM
%
% Revision date: June 12, 2020
% (C) Perfect Gidisu, Michiel Hochstenbach 2020

if nargin < 3 || isempty(k), k = 2; end
if nargin < 4 || isempty(U)
  [U, V, X, ~, ~] = gsvd(A,B,0);
  U = fliplr(U(:,end-k+1:end));       % Select largest (C,S) pairs
  V = fliplr(V(:,end-k+1:end));
  X = fliplr(X(:,end-k+1:end));
end

for j = 1:k
  [~, icol(j)] = max(abs(X(:,j)));    % Iterative selection and projection
  [~, iA(j)] = max(abs(U(:,j)));
  [~, iB(j)] = max(abs(V(:,j)));
  X(:,j+1:end) = X(:,j+1:end) - X(:,1:j) * (X(icol,1:j) \ X(icol,j+1:end));
  U(:,j+1:end) = U(:,j+1:end) - U(:,1:j) * (U(iA,1:j) \ U(iA,j+1:end));
  V(:,j+1:end) = V(:,j+1:end) - V(:,1:j) * (V(iB,1:j) \ V(iB,j+1:end));
%   X(:,j+1:end) = X(:,j+1:end) - X(:,j) * (X(icol,j) \ X(icol,j+1:end));
%   U(:,j+1:end) = U(:,j+1:end) - U(:,j) * (U(iA,j) \ U(iA,j+1:end));
%   V(:,j+1:end) = V(:,j+1:end) - V(:,j) * (V(iB,j) \ V(iB,j+1:end));
end
M = A(:,icol) \ (A / A(iA,:));
N = B(:,icol) \ (B / B(iB,:));
