rng(0)
N = 400; D = 30; gap=3;
% In B, all the data pts are from the same distribution, which has different variances in three subspaces.
B = zeros(N, D);
B(:,1:10) = normrnd(0,10,[N,10]);
B(:,11:20) = normrnd(0,3,[N,10]);
B(:,21:30) = normrnd(0,1,[N,10]);


% In A there are four clusters.
A = zeros(N, D);
A(:,1:10) = normrnd(0,10,[N,10]);
% group 1
A(1:100, 11:20) = normrnd(0,1,[100,10]);
A(1:100, 21:30) = normrnd(0,1,[100,10]);
% group 2
A(101:200, 11:20) = normrnd(0,1,[100,10]);
A(101:200, 21:30) = normrnd(gap,1,[100,10]);
% group 3
A(201:300, 11:20) = normrnd(2*gap,1,[100,10]);
A(201:300, 21:30) = normrnd(0,1,[100,10]);
% group 4
A(301:400, 11:20) = normrnd(2*gap,1,[100,10]);
A(301:400, 21:30) = normrnd(gap,1,[100,10]);

A_labels = zeros(N,1);
A_labels(101:200)=1;
A_labels(201:300)=2;
A_labels(301:400)=3;

A_center=A-mean(A);
B_center=B-mean(B);

[~,~,V]=svd(A_center);
[~,~,XX,~,~]=gsvd(A_center,B_center,0);

XX=fliplr(XX);
YY=inv(XX'); % right generalized singular vectors


sv=A_center*V(:,1:10); % projection onto the right singular vectors 
gsv=A_center*YY(:,1:10); % projection onto the right generalized singular vectors


icol=cur_deim(V(:,1:10),10);
icol1 = cur_deim(XX(:,1:10),10);


C_cur=A_center(:,icol);
C_gcur=A_center(:,icol1);


for i=0:3
    subplot(2,2,1)
    scatter(sv(A_labels==i,1),sv(A_labels==i,2),'.')
    xlabel('singular vector 1')
    ylabel('singular vector 2')
    title("Projected data onto the 2 leading right singular vectors",'FontSize',12)
    hold on;
    subplot(2,2,2)
    scatter(gsv(A_labels==i,1),gsv(A_labels==i,2),'.')
    xlabel('generalized singular vector 1')
    ylabel('generalized singular vector 2')
    title("Projected data onto the 2 leading right generalized singular vectors",'FontSize',12)
     hold on;
    subplot(2,2,3)
    scatter(C_cur(A_labels==i,1),C_cur(A_labels==i,2),'.')
    title("Visualizing data set A using CUR's first 2 important columns",'FontSize',12)
    xlabel('CUR 1st important column')
    ylabel('CUR 2nd important column')
    
     hold on;
    subplot(2,2,4)
    scatter(C_gcur(A_labels==i,1),C_gcur(A_labels==i,2),'.')
    xlabel('GCUR 1st important column')
    ylabel('GCUR 2nd important column')
    title("Visualizing data set A using GCUR's first 2 important columns",'FontSize',12)
     hold on;
   
    
end
