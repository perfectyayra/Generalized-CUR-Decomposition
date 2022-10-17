rng(0)
n=1000;
m=10000;
k=50;
C=toeplitz(0.99.^(0:n-1));
R=chol(C);
A_exact=zeros(m,n);

for i=1:10
    x=sprand(m, 1,0.025);
    y=sprand(n,1,0.025);
    A_exact=A_exact+(2/i)*x*y.';
end 
for i=11:50
    x=sprand(m, 1,0.025);
    y=sprand(n,1,0.025);
    A_exact=A_exact+(1/i)*x*y.';
end 
epsilon=0.2 %noise level

for j=1:5
    Correlated_noise= randn(m,n)*R ;
    E=epsilon*(norm(A_exact)/norm(Correlated_noise))*Correlated_noise;
    A=A_exact+E;
    [U,S,V]=svd(A,0);
    [U2,~,X2,C2,~]=gsvd(A,R,0);
    U2=fliplr(U2);
    X2=fliplr(X2);
    C2=rot90(C2,2);
    [icol, irow,~] = cur_deim(A,k);
    [icol2, iA,~,~,~] = gcur_deim(A,R,k);
    for i=1:k
         
        A_svd=U(:,1:i)*S(1:i,1:i)*V(:,1:i)';
        A_gsvd=U2(:,1:i)*C2(1:i,1:i)*X2(:,1:i)';
        M=A(:,icol(1:i))\A/A(irow(1:i),:);
        M2=A(:,icol2(1:i))\A/A(iA(1:i),:);
        CUR=A(:,icol(1:i))*M*A(irow(1:i),:);
        CUR2=A(:,icol2(1:i))*M2*A(iA(1:i),:);
        err(j,i)=norm(A_exact-CUR)/norm(A_exact);
        err1(j,i)=norm(A_exact-CUR2)/norm(A_exact);
        err2(j,i)=norm(A_exact-A_svd)/norm(A_exact);
        err3(j,i)=norm(A_exact-A_gsvd)/norm(A_exact);
     end

    
end

plot(1:50,mean(err),'r-');
hold on;
plot(1:50,mean(err1),'b-'); 
plot(1:50,mean(err2),'r-.');
plot(1:50,mean(err3),'b-.');

legend('DEIM-CUR','DEIM-GCUR','SVD','GSVD')
