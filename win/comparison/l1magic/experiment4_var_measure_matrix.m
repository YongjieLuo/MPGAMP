% put key subdirectories in path if not already there
path(path, './Optimization');

close all;
clear all;
clc;

load('E:\Luo\data\mnist_all');

s = reshape(train8(666,:), 28, 28)';
figure; imshow(s);
figure; mesh(double(s));

x = reshape(s, 28*28, 1);
N = size(x, 1);

% number of observations to make
%K = round(N*4/10) + 40;
N = 28*28;
M = 20;
K = M*M;

total = 20;
ys = zeros(K, total);
for i = 1:1:total
    % measurement matrix
    disp(sprintf('Creating measurment matrix A %d', i));
    A = randn(K,N);
    A = orth(A')';

    % observations
    x = double(x);
    y = A*x;
    ys(:, i) = y;
    %figure; imshow(reshape(int8(y), M, M));
    %figure; mesh(reshape(y, M, M));

    %% OMP recover
    %yhat = omp(A, y)';
    %yp = int8(yhat)';
    %figure; imshow(reshape(yp, 28, 28));
    %figure; mesh(reshape(yhat, 28, 28));
end
figure; 
for i = 1:total, 
    subplot(4,5,i), hist(ys(:, i), 40), 
end

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
% dfx = x - mean(x);
% ndfx = dfx / sqrt(dfx' * dfx);
% nx = reshape(ndfx, 28, 28);
% [nx1, nx2 ,nx3] = svd(cov(nx));
% ex = sum(sum(nx2(1:4, 1:4)));
% 
% % compressed signal
% dfy = y - mean(y);
% ndfy = dfy / sqrt(dfy' * dfy);
% ny = reshape(ndfy, 20, 20);
% [ny1, ny2 ,ny3] = svd(cov(ny));
% ey = sum(sum(ny2(1:12, 1:12)));
% 
% disp(sprintf('the energy of Original signal difference = %.4f, the energy of Compressed signal difference = %.4f\n', ex, ey));
% 
% figure, pmusic(dfx, 4);
% figure, pmusic(dfy, 12);

% large scale
% Afun = @(z) A*z;
% Atfun = @(z) A'*z;
% tic
% xp = l1eq_pd(x0, Afun, Atfun, y, 1e-3, 30, 1e-8, 200);
% toc





