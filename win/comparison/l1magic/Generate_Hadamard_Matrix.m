function H = Generate_Hadamard_Matrix(m) 

% This builds the Hadamard matrix
H=zeros(m,m);
for k=1:1:m,
    disp(k/m);
    temp=zeros(m,1);
    temp(k)=1; 
    H(:,k)=fht(temp);
end;
