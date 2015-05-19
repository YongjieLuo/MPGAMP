clear all; close all;
clc;

trialSetup;

trialNum = 20;

% Specify problem size
N = 1000;
K = 50;

baseM = 100;
addM = 150;

%Specify noise model.
%Model is (1-lambda) Nor(0,nu0) + lambda Nor(0,nu1)
nu0 = 1e-3;

%Set options for GAMP
GAMP_options = GampOpt;
MpGAMP_options = MpGAMPOpt; %initialize the options object

methods = {'(a) MP', '(b) OMP', '(c) BP', '(d) AMP', '(e) GAMP', '(f) MPGAMP'};
[~,J] = size(methods);

err_absolute = cell(1,J);
err_relative = cell(1,J);
for j=1:J
    err_absolute{j} = zeros(addM+1, trialNum);
    err_relative{j} = zeros(addM+1, trialNum);
end

for n = 1:addM+1
    % reset random number generator
    rand('state', n);
    randn('state', 3*n);

    M = baseM+n-1;    
    
    % Other trials
    for i = 1:trialNum
        
        A = randn(M,N);
        A = randn(M,N) / sqrt(N);
        %A = A*diag(1 ./ sqrt(diag(A'*A)));      
        
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
        
        % MP
        [xhatMP,itersMP,usedMP] = SolveMP(A,y,N,30,1e-6,0,0,1e-8);
        err_absolute{1}(n,i) = norm(x-xhatMP);        
        err_relative{1}(n,i) = norm(x-xhatMP)/norm(x);
                
        % OMP
        [xhatOMP,itersOMP,usedOMP] = SolveOMP(A,y,N,30,1e-6,0,0,1e-8);
        err_absolute{2}(n,i) = norm(x-xhatOMP);                
        err_relative{2}(n,i) = norm(x-xhatOMP)/norm(x);
        
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
        err_absolute{3}(n,i) = norm(x-xhatBP);        
        err_relative{3}(n,i) = norm(x-xhatBP)/norm(x);
        
        % AMP
        T = 1000;  % Number of iterations
        tol = 0.001;  % Tolerance
        xhatAMP = reconstructAmp(A,y,T,tol,x,0);
        err_absolute{4}(n,i) = norm(x-xhatAMP);        
        err_relative{4}(n,i) = norm(x-xhatAMP)/norm(x);
        
        % GAMP with Gaussian noise model
        %Input channel
        inputEst = AwgnEstimIn(0, 1);
        inputEst = SparseScaEstim(inputEst,K/N);
        %Output channel
        %nu0 = sigma;
        outputEst = AwgnEstimOut(y, nu0);
        %Run GAMP
        [resGAMP,~,~,~,~,~,~,~, estHistGAMP] = ...
            gampEst(inputEst, outputEst, A, GAMP_options);
        xhatGAMP = resGAMP;
        err_absolute{5}(n,i) = norm(x-xhatGAMP);        
        err_relative{5}(n,i) = norm(x-xhatGAMP)/norm(x);
        
        % MPGAMP
        [result, history] = ...
            MpGAMP(inputEst, outputEst, A, y, x, K, MpGAMP_options);
        C = (result.C)';
        xhatMPGAMP = result.xhat;
        err_absolute{6}(n,i) = norm(x-xhatMPGAMP);        
        err_relative{6}(n,i) = norm(x-xhatMPGAMP)/norm(x);
                      
    end
    
    fprintf('measurement number %d done\n', baseM+n-1)
end

disp('Done!');

MP_mean = mean(err_relative{1},2);
MP_std = std(err_relative{1},0,2);
OMP_mean = mean(err_relative{2},2);
OMP_std = std(err_relative{2},0,2);
BP_mean = mean(err_relative{3},2);
BP_std = std(err_relative{3},0,2);
AMP_mean = mean(err_relative{4},2);
AMP_std = std(err_relative{4},0,2);
GAMP_mean = mean(err_relative{5},2);
GAMP_std = std(err_relative{5},0,2);
MPGAMP_mean = mean(err_relative{6},2);
MPGAMP_std = std(err_relative{6},0,2);

figure
plot(baseM+(0:addM),MP_mean,'m-v','LineWidth',1);
hold on;
plot(baseM+(0:addM),OMP_mean,'k-+','LineWidth',1);
hold on;
plot(baseM+(0:addM),BP_mean,'c-s','LineWidth',1);
hold on;
plot(baseM+(0:addM),AMP_mean,'g-x','LineWidth',1);
hold on;
plot(baseM+(0:addM),GAMP_mean,'y-o','LineWidth',1);
hold on;
plot(baseM+(0:addM),MPGAMP_mean,'r-^','LineWidth',1);
xlabel('Number of Measurements','FontSize',12); ylabel('Average Relative Error','FontSize',14);
box on; grid on;
axis([95,250,-0.1,1.1]);
legend('MP','OMP','BP','AMP','GAMP','MPGAMP',1);

%save DataExp3_N1000.mat err_absolute err_relative trialNum baseM addM N K nu0 methods
