A=load('cell_patient_035.txt');
B=load('cell_background.txt');

A=full(spconvert(A(2:end,:)));
B=full(spconvert(B(2:end,:)));

A_labels=load('group_mem_035.txt');

A_center=normalize(A);
B_center=normalize(B);

[U,~,V]=svd(A_center);
[~,~,XX,~,~]=gsvd(A_center,B_center,0);

XX=fliplr(XX);
YY=inv(XX'); % generalized right singular vectors 

sv=A_center*V(:,1:3);
gsv=A_center*YY(:,1:3);

[icol,~,~]=cur_deim(A_center,3);
[icol1,~,~,~,~] = gcur_deim(A_center,B_cente,3);

C_cur=A_center(:,icol);
C_gcur=A_center(:,icol1);


for i=1:2
    subplot(2,2,1)
    scatter3(sv(A_labels==i,1),sv(A_labels==i,2),sv(A_labels==i,3),'.')
    title("Projected data onto the 3 leading right singular vectors",'FontSize',12)
    hold on;
    subplot(2,2,2)
    scatter3(gsv(A_labels==i,1),gsv(A_labels==i,2),gsv(A_labels==i,3),'.')
    title("Projected data onto the 3 leading right generalized singular vectors",'FontSize',12)
     hold on;
    subplot(2,2,3)
    scatter3(C_cur(A_labels==i,1),C_cur(A_labels==i,2),C_cur(A_labels==i,3),'.')
    title("Visualizing data set A using CUR's first 3 important columns",'FontSize',12)
     hold on;
    subplot(2,2,4)
    scatter3(C_gcur(A_labels==i,1),C_gcur(A_labels==i,2),C_gcur(A_labels==i,3),'.')
    title("Visualizing data set A using GCUR's first 3 important columns",'FontSize',12)
     hold on;
   
    
end
