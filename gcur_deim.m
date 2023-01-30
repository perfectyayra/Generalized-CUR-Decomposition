function [icol, iA, iB, M, N] = gcur_deim(A, B, k)

%CUR_PAIR_DEIM  DEIM incurred CUR decomposition for matrix pairs
% function [irow, iA, iB, M, N] = gcur_deim(A, B, k)
% Matrix A and B should have same number of columns
% k is the desired rank of the approximation 
% icol contains the column indices selected
% iA contains the row indices of matrix A that has been selected
% iB contains the row indices of matrix B that has been selected
% CA = A(:,icol);  CB = B(:,icol);  RA = A(iA,:);  RB = B(iB,:);
% M and N are the middle matrix of the gcur decomposition of A and B, respectively
% See also cur_deim
%
% Reference: Gidisu and Hochstenbach 2022 (https://doi.org/10.1137/21M1432119)


[U, V, Y, ~, ~] = gsvd(A,B,0);  %matlab gsvd implementation
U = fliplr(U);       % Select largest (C,S) pairs
V = fliplr(V);
Y = fliplr(Y);

icol=zeros(1,k);
iA=zeros(1,k);
iB=zeros(1,k);
for j = 1:k
  [~, icol(j)] = max(abs(Y(:,j)));    % Iterative selection and projection
  [~, iA(j)] = max(abs(U(:,j)));
  [~, iB(j)] = max(abs(V(:,j)));
  if j<k
    Y(:,j+1) = Y(:,j+1) - Y(:,1:j) * (Y(icol(1:j),1:j) \ Y(icol(1:j),j+1));
    U(:,j+1) = U(:,j+1) - U(:,1:j) * (U(iA(1:j),1:j) \ U(iA(1:j),j+1));
    V(:,j+1) = V(:,j+1) - V(:,1:j) * (V(iB(1:j),1:j) \ V(iB(1:j),j+1));
  end
end
C_A = A(:,icol);  
C_B = B(:,icol);  
R_A = A(iA,:);  
R_B = B(iB,:);
M = pinv(C_A)*A*pinv(R_A);
N =pinv(C_B)*B*pinv(R_B);
