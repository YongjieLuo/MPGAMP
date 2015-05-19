function [Haar]=Generate_Haar_Matrix(n)

D1=sparse(n,n);
v=sparse([1 zeros(1,n-2), -1]/2);
for k=1:1:n
    D1(k,:)=v;
    v=[v(end),v(1:end-1)];
end;
D2=sparse(n,n);
v=[1 1 zeros(1,n-4), -1 -1]/4;
for k=1:1:n
    D2(k,:)=v;
    v=[v(end),v(1:end-1)];
end;
S1=abs(D1);
S2=abs(D2);
Haar=[kron(S2,S2),kron(S2,D2),kron(D2,S2),kron(D2,D2),...
                             kron(S1,D1),kron(D1,S1),kron(D1,D1)];
return;
