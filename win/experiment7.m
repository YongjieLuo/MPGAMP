clear all; close all;
clc;

trialSetup;

trialNum = 20;

Alpha = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
alphaNum = length(Alpha);

del = 0.5; %ratio of m/n
rho = 0.1; %ratio of sparsity to number of measurements

% Setup and global options
%Specify problem size
N = 1000;
M = ceil(del*N);
K = floor(rho*M);

%Specify noise model.
%Model is (1-lambda) Nor(0,nu0) + lambda Nor(0,nu1)
nu0 = 1e-3;

%Set options for GAMP
GAMP_options = GampOpt;
MpGAMP_options = MpGAMPOpt;  %initialize the options object

methods = {'(a) OMP', '(b) BP', '(c) AMP', '(d) MPGAMP'};
[~,J] = size(methods);

err_absolute = cell(1,J);
err_relative = cell(1,J);
for j=1:J
    err_absolute{j} = zeros(alphaNum, trialNum);
    err_relative{j} = zeros(alphaNum, trialNum);
end

for n = 1:alphaNum
    % reset random number generator
    rand('state', n);
    randn('state', 2*n);
    
    for i = 1:trialNum

        alpha = Alpha(n);
        R = ceil(alpha*N);        
        
        P = randn(M,R);
        Q = randn(R,N);
        
        A = (1/N)*P*Q;

        % Generate the true signal
        L = 200;
        X = randn(N,L);
        for ll = 1:L
            yada = randperm(N);
            yada2 = zeros(N,1);
            yada2(yada(1:K)) = 1;
            X(:,ll) = X(:,ll) .* yada2;
        end

        x = X(:,50);
        %Generate the uncorrupted measurements
        z = A*x;

        % Generate noisy signal
        y = z + sqrt(nu0)*randn(size(z));        
        
        % OMP
        [xhatOMP,itersOMP,usedOMP] = SolveOMP(A,y,N,30,1e-6,0,0,1e-8);
        err_absolute{1}(n,i) = norm(x-xhatOMP);                
        err_relative{1}(n,i) = norm(x-xhatOMP)/norm(x);
        
        % BP
        % Determine noise variance
        SNR = 20;
        Z = A*X;
        sigma = norm(reshape(Z,[],1))^2/M/L*10^(-SNR/10);
        % initial guess = min energy
        x0 = A'*y;
        % take epsilon a little bigger than sigma*sqrt(K)
        epsilon =  sigma*sqrt(M)*sqrt(1 + 2*sqrt(2)/sqrt(M));
        %xhatBP = l1qc_logbarrier(x0,A,[],y,epsilon,1e-3);
        % large scale
        Afun = @(z) A*z;
        Atfun = @(z) A'*z;
        xhatBP = l1qc_logbarrier(x0,Afun,Atfun,y,epsilon,1e-3);
        err_absolute{2}(n,i) = norm(x-xhatBP);        
        err_relative{2}(n,i) = norm(x-xhatBP)/norm(x);
        
        % AMP
        T = 1000;  % Number of iterations
        tol = 0.001;  % Tolerance
        xhatAMP = reconstructAmp(A,y,T,tol,x,0);
        err_absolute{3}(n,i) = norm(x-xhatAMP);        
        err_relative{3}(n,i) = norm(x-xhatAMP)/norm(x);
        
        % GAMP with Gaussian noise model
        %Input channel
        inputEst = AwgnEstimIn(0, 1);
        inputEst = SparseScaEstim(inputEst,K/N);
        %Output channel
        outputEst = AwgnEstimOut(y, nu0);
        
        % MPGAMP
        [result, history] = ...
            MpGAMP(inputEst, outputEst, A, y, x, K, MpGAMP_options);
        C = (result.C)';
        xhatMPGAMP = result.xhat;
        err_absolute{4}(n,i) = norm(x-xhatMPGAMP);        
        err_relative{4}(n,i) = norm(x-xhatMPGAMP)/norm(x);
        
    end  % end for 
    
    fprintf('Alpha(%d) done\n', n)
end  % end for

disp('Done!');

OMP_mean = mean(err_relative{1},2);
OMP_std = std(err_relative{1},0,2);
BP_mean = mean(err_relative{2},2);
BP_std = std(err_relative{2},0,2);
AMP_mean = mean(err_relative{3},2);
AMP_std = std(err_relative{3},0,2);
MPGAMP_mean = mean(err_relative{4},2);
MPGAMP_std = std(err_relative{4},0,2);

figure
plot(Alpha,OMP_mean,'k-+','LineWidth',1);
hold on;
plot(Alpha,BP_mean,'g-s','LineWidth',1);
hold on;
plot(Alpha,AMP_mean,'b-x','LineWidth',1);
hold on;
plot(Alpha,MPGAMP_mean,'r-^','LineWidth',1);
xlabel('\alpha','FontSize',12); ylabel('Average Relative Error','FontSize',12);
box on; grid on;
axis([0.18,1.02,-0.1,2.0]);
legend('OMP','BP','AMP','MPGAMP',1);

%save DataExp7_N1000.mat err_absolute err_relative trialNum N M K Alpha nu0 methods
