%%Multi-criteria decision-making¡ª¡ªTOPSIS

clear

data=xlsread('D:\OneDrive\MyImportance\IO_optimization\results\eps_plot.xlsx','C1:C708');

dataplot=zeros(236,3);
for i=1:236
    dataplot(i,1)=data(3*(i-1)+1);
    dataplot(i,2)=data(3*(i-1)+2);
    dataplot(i,3)=data(3*(i-1)+3);
end

load epsdata_e.mat dataplot
A=cell2mat(dataplot);

W1=[0.55 0.15 0.15 0.15];
W2=[0.15 0.55 0.15 0.15];
W3=[0.15 0.15 0.55 0.15];
W4=[0.15 0.15 0.15 0.55];
W5=[0.25 0.25 0.25 0.25];

A=A';
A(:,2)=-A(:,2);
A(:,4)=-A(:,4);

[ score1 ] = TOPSIS(A,W1);%Calculate the TOPSIS values for all scenarios (30) for each d (d=11) taking values in the preference one scenario
[ score2 ] = TOPSIS(A,W2);%The smaller the TOPSIS value under the preference 2 scenario, the better the program.
[ score3 ] = TOPSIS(A,W3);
[ score4 ] = TOPSIS(A,W4);
[ score5 ] = TOPSIS(A,W5);
% a1=cell2mat(score1);
% a2=cell2mat(score2);
% a3=cell2mat(score3);
[max_a1,ind1]=max(score1,[],1);%Find the optimal solution for each value of d for the preference 1 scenario. ind is the number of rows where the optimal solution is located.
[max_a2,ind2]=max(score2,[],1);
[max_a3,ind3]=max(score3,[],1);
[max_a4,ind4]=max(score4,[],1);
[max_a5,ind5]=max(score5,[],1);
ind=[ind1' ind2' ind3' ind4' ind5'];


