function xhat = omp(A, b)

% A: (n x m) matrix 
% b: m columns vector

% Apply OMP
thrOMP=1e-4;
[n, m] = size(A);
r=b;
SS=[];

disp('Creating another measurment matrix...');
A = randn(n,m);
A = orth(A')';
disp('omp');

while r'*r>thrOMP,
    Z=abs(A'*r);
    posZ=find(Z==max(Z));
    SS=sort([SS,posZ(1)]);
    r=b-A(:,SS)*pinv(A(:,SS))*b;    
end;
        
xOMP=zeros(m,1);
xOMP(SS)=pinv(A(:,SS))*b;

xhat = xOMP;


