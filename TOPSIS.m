function [ output_args ] = TOPSIS(A,W)

    %A is the decision matrix, W is the weight matrix, M is the column where positive indicators are located, N is the column where negative indicators are located.
    [ma,na]=size(A);          %ma is the number of rows of the A matrix and na is the number of columns of the A matrix
    for i=1:na
        A(:,i)=A(:,i)/norm(A(:,i),2);
    end
    for i=1:na
        B(:,i)=A(:,i)*W(i);     %Loop by column to get [weighted normalized matrix]
    end
    V1=zeros(1,na);            %Initialize positive and negative ideal solutions
    V2=zeros(1,na);
    BMAX=max(B);               %Take the maximum and minimum values of each column of the weighted normalized matrix
    BMIN=min(B);               %
    for i=1:na
        %if i<=size(M,2)        %Loop to get the ideal solution and negative ideal solution
        V1(i)=BMAX(i);
        V2(i)=BMIN(i);
        %end
        %if i=size(N,2)
        %V1(N(i))=BMIN(N(i));
        %V2(N(i))=BMAX(N(i));
        %end
    end
    
    for i=1:ma                 %Loop by row to find the closeness of each solution
        C1=B(i,:)-V1;
        S1(i)=norm(C1);        %S1,S2 are the distances from the positive and negative ideal points, respectively, using second-order paradigms
        
        C2=B(i,:)-V2;
        S2(i)=norm(C2);
        T(i)=S2(i)/(S1(i)+S2(i));     %T is the closeness
    end
    output_args=T';
end
