load('E:\Luo\data\mnist_all')

% put key subdirectories in path if not already there
path(path, './Optimization');

close all;
clear A K N s x y yhat yp;
clc;

s = reshape(train1(2153,:), 28, 28)';
figure,
imshow(s);

x = reshape(s, 28*28, 1);
N = size(x, 1);

% number of observations to make
%K = round(N*4/10) + 40;
K = 20*20;

% measurement matrix
disp('Creating measurment matrix...');
A = randn(K,N);
A = orth(A')';
disp('Done.');
	
% observations
x = double(x);
y = A*x;

figure,
imshow(reshape(int8(y), 20, 20));

%% OMP recover
yhat = omp(A, y);
yp = int8(yhat);
figure,
imshow(reshape(yp, 28, 28));
figure; mesh(reshape(yhat, 28, 28));


%%
N = 20;
nums = zeros(28*28, N);    % every column record a number
csins = zeros(20*20, N); % every column record the compressed instance of that number
for i = 1:N
    nums(:, i) = train1(i,:)';
end

s = reshape(nums(:, 1), 28, 28)';
figure,
imshow(s);

% xbar = mean(nums');
% 
% x = reshape(int8(xbar), 28, 28)';
% figure,
% imshow(x);

% the same A
for i = 1:N
    csins(:, i ) = A * (double(nums(:, i))); 
end

y = reshape(int8(csins(:, 1)), 20, 20)';
figure,
imshow(y);

yhat = omp(A, csins(:, 1));
yp = int8(yhat);
figure,
imshow(reshape(yp, 28, 28)');
figure; mesh(reshape(yhat, 28, 28)');

%
ybar = mean(csins');
ys = reshape(int8(ybar), 20, 20)';
figure,
imshow(ys);
figure; mesh(reshape(ybar, 20, 20));

%% BP recover
% initial guess = min energy
% x0 = A'*y;
% 
% % solve the LP
% tic
% xp = l1eq_pd(x0, A, [], y, 1e-3);
% toc
% 
% xpp = int8(xp);
% figure,
% imshow(reshape(xpp, 28, 28));



%% caculate energy
% original signal
dfx = x - mean(x);
ndfx = dfx / sqrt(dfx' * dfx);
nx = reshape(ndfx, 28, 28);
[nx1, nx2 ,nx3] = svd(cov(nx));
ex = sum(sum(nx2(1:4, 1:4)));

% compressed signal
dfy = y - mean(y);
ndfy = dfy / sqrt(dfy' * dfy);
ny = reshape(ndfy, 20, 20);
[ny1, ny2 ,ny3] = svd(cov(ny));
ey = sum(sum(ny2(1:12, 1:12)));

disp(sprintf('the energy of Original signal difference = %.4f, the energy of Compressed signal difference = %.4f\n', ex, ey));

figure, pmusic(dfx, 4);
figure, pmusic(dfy, 12);

% large scale
% Afun = @(z) A*z;
% Atfun = @(z) A'*z;
% tic
% xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);
% toc





