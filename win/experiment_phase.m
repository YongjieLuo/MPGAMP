clear all; close all;
clc;

trialSetup;

trialNum = 20;

% Specify problem size
N = 200;

delta = linspace(0.05,1,21);
rho = linspace(0.05,1,21);
delta = delta(1:end-1);
rho = rho(1:end-1);

%Specify noise model.
%Model is (1-lambda) Nor(0,nu0) + lambda Nor(0,nu1)
nu0 = 1e-3;

%Set options for GAMP
GAMP_options = GampOpt;
MpGAMP_options = MpGAMPOpt; %initialize the options object

methods = {'(a) AMP', '(b) GAMP', '(c) MPGAMP'};
[~,J] = size(methods);

errL0 = cell(1,J);
errL2_absolute = cell(1,J);
errL2_relative = cell(1,J);
for j=1:J
    errL0{j}          = zeros(length(delta), length(rho), trialNum);
    errL2_absolute{j} = zeros(length(delta), length(rho), trialNum);
    errL2_relative{j} = zeros(length(delta), length(rho), trialNum);
end


tic;
for p = 1:length(delta)
    % reset random number generator
    rand('state', p);
    randn('state', 3*p);

    M = floor(N*delta(p));
    
    for q = 1:length(rho)
        K = ceil(rho(q)*M);
        
        % Other trials
        for i = 1:trialNum
            %errL0{1}(p,q,i) = i;
            
            A = randn(M,N) / sqrt(N);
%             %A = randn(M,N);
%             %A = A*diag(1 ./ sqrt(diag(A'*A)));      
% 
            L = 20;
            X = randn(N,L);
            for ll = 1:L
                yada = randperm(N);
                yada2 = zeros(N,1);
                yada2(yada(1:K)) = 1;
                X(:,ll) = X(:,ll) .* yada2;
            end

            x = X(:,10);
            %Generate the uncorrupted measurements
            z = A*x;
            % Generate noisy signal
            y = z + sqrt(nu0)*randn(size(z));


            % AMP
            T = 1000;  % Number of iterations
            tol = 0.001;  % Tolerance
            xhatAMP = reconstructAmp(A,y,T,tol,x,0);
            errL0{1}(p,q,i)          = length(find(abs(x-xhatAMP) > 1e-2));
            errL2_absolute{1}(p,q,i) = norm(x-xhatAMP);        
            errL2_relative{1}(p,q,i) = norm(x-xhatAMP)/norm(x);

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
            errL0{2}(p,q,i)          = length(find(abs(x-xhatGAMP) > 1e-2));
            errL2_absolute{2}(p,q,i) = norm(x-xhatGAMP);        
            errL2_relative{2}(p,q,i) = norm(x-xhatGAMP)/norm(x);

            % MPGAMP
            [result, history] = ...
                MpGAMP(inputEst, outputEst, A, y, x, K, MpGAMP_options);
            C = (result.C)';
            xhatMPGAMP = result.xhat;
            errL0{3}(p,q,i)          = length(find(abs(x-xhatMPGAMP) > 1e-2));
            errL2_absolute{3}(p,q,i) = norm(x-xhatMPGAMP);        
            errL2_relative{3}(p,q,i) = norm(x-xhatMPGAMP)/norm(x); 
        end

        fprintf('K=%d, M=%d, N=%d done.\n', K, M, N)
    end
    
end
timeAll = toc;

disp('Done!');


%save DataExp9_phase_transition.mat errL0 errL2_absolute errL2_relative timeAll trialNum N delta rho nu0 methods
