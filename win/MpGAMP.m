% The MIT License (MIT)
% Copyright © 2015 <copyright holders>
%
% Permission is hereby granted, free of charge, to any person obtaining 
% a copy of this software and associated documentation files (the “Software”), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the Software 
% is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
% PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
% OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% Thanks for gamp software package
% http://sourceforge.net/projects/gampmatlab/
%
% This code is open to anyone interested in using it.
% Yongjie Luo, E.E., UESTC
% yongjie.luo@hotmail.com
% 2014.5

function [estFin, estHist] = MpGAMP(scaEstIn, scaEstOut, A, y, x, K, opt)

nit     = opt.nit;              % number of iterations
nitMin  = opt.nitMin;
verbose = opt.verbose;          % Print results in each iteration?
tol = opt.tol;                  % Convergence tolerance
pvarStep = opt.pvarStep;        % include stepsize in pvar?
pvarMin = opt.pvarMin;          % minimum value of pvar
zvarToPvarMax = opt.zvarToPvarMax;  % maximum zvar/pvar ratio
step    = opt.step;             % stepsize
stepMin = opt.stepMin;
stepMax = opt.stepMax;
stepIncr = opt.stepIncr;
stepDecr = opt.stepDecr;
stepWindow = opt.stepWindow;
adaptStep = opt.adaptStep;
adaptStepBethe = opt.adaptStepBethe;
maxBadSteps = opt.maxBadSteps;  % maximum number of allowed failed iterations
maxStepDecr = opt.maxStepDecr;  % amount to decrease maxStep after failures



% Get dimensions
[m,n] = size(A);

% Get default initialization values
[xhat,xvar,valIn] = scaEstIn.estimInit();
valIn = sum( valIn(:) );
valOpt = -Inf;
val = nan;
failCount = 0;

xhat = xhat*ones(n,1);
xvar = xvar*ones(n,1);

% Replace default initialization with user-provided values
if ~isempty(opt.xhat0)
    xhat = opt.xhat0;
end
if ~isempty(opt.xvar0)
    xvar = opt.xvar0;
end

% Declare variables, Continue with initialization
pvarOld = nan(m,1);     % will test for NaN later 
zvarOld = nan(m,1);     % will test for NaN later 

shatOld = zeros(m,1);
shatDamp = zeros(m,1);    % will test for NaN later
svarDamp = zeros(m,1);    % will test for NaN later
xhatDamp = nan(n,1);    % will test for NaN later

zhat0 = nan(m,1);
zvar0 = nan(m,1);

xhatFinal = nan(n,1);
xhatIFinal = nan;
uFinal = nan(m,1);

C = [];
maskIndics = 1:n;
xvarThresh = 1e-6;
xhatC = [];

pursuitDone = false;
k = 0;
e = norm(y);

A2 = (abs(A).^2);

while ~pursuitDone

    k = k + 1;
    
    if k > ceil(n/2)
        pursuitDone = true;
    end

    if k > K
        pursuitDone = true;
        %disp(['Indics of non zero xhat element are: ' num2str(C)]);
        break;
    end    
    
    % residual coherence
    if k == 1  % set C is empty 
        cohen = A' * y;
    else
        v = norm(u);

        AMaskC = A; AMaskC(:,C) = 0;
        cohen = AMaskC' * u;
    end

    [~, I] = max(abs(cohen));
    
    maskIndics = setdiff(maskIndics, I);

    AI = A(:,I);
    AI2 = A2(:,I);
    xhatI = xhat(I);
    xvarI = xvar(I);
    if k ~= 1  % set C is not empty
        AC = A(:,C);
        AC2 = A2(:,C);
        xhatC = xhat(C);
        xvarC = xvar(C);
    end
    
    % Control variables to terminate the iterations
    stop = false;
    it = 0;
        
    % Main iteration loop
    while ~stop

        % Iteration count
        it = it + 1;

        if k ~= 1  % set C is not empty
            % Output linear stage with no A uncertainty
            zvar = AC2*xvarC + AI2*xvarI;
            
            % Continued output linear stage
            zhat = AC*xhatC + AI*xhatI;
        else
            % Output linear stage with no A uncertainty
            zvar = AI2*xvarI;
            
            % Continued output linear stage
            zhat = AI*xhatI;
        end
            
        pvar = zvar;  % notation coincide with paper
        
        % Step in pvar
        if pvarStep
            if (it==1)
                if any(isnan(pvarOld)),    % if user didn't specify opt.pvarOld0
                    pvarOld = pvar;        % equivalent to making step=1
                end
                if any(isnan(zvarOld)),    % if user didn't specify opt.A2xvarOld0
                    zvarOld = zvar;    % equivalent to making step=1
                end
            end
            pvarDamp = (1-step)*pvarOld + step*pvar;
            zvarDamp = (1-step)*zvarOld + step*zvar;
        end
        pvarRobust = max(pvarDamp, pvarMin); % At very high SNR, use very small pvarMin!

        % Continued output linear stage
        phat = zhat - shatOld.*zvarDamp; % Note: uses A2xvar rather than pvar

        
        % Compute expected log-likelihood of the output and add to negative 
        % KL-divergence of the input, giving the current utility function 
        if ~adaptStepBethe
            valOut = sum(sum(scaEstOut.logLike(zhat,pvarRobust)));
        else
            valOut = sum(sum(scaEstOut.logScale(zhat,pvarRobust,phat)));
        end
        val = valOut + valIn;
        
        % An iteration "passes" if any of below is true: 
        % 1. Adaptive stepsizing is turned off
        % 2. Current stepsize is so small it can't be reduced 
        % 3. The current utility at least as large as the worst in the stepWindow  
        % Also, we force a pass on the first iteration else many quantities undefined
        stopInd = length(valOpt);
        startInd = max(1,stopInd - stepWindow);
        valMin = min(valOpt(startInd:stopInd));
        %pass = (~adaptStep) || (step <= stepMin) || isempty(valMin) || (val >= valMin);
        pass = (it==1) || (~adaptStep) || (step <= stepMin) || (val >= valMin);
        
        % If pass, set the optimal values and compute a new target shat and snew.
        if (pass)
            % Save states 
            zvarOld = zvarDamp;
            pvarOld = pvarRobust;
            shatOld = shatDamp;
            svarOld = svarDamp;
            xhatDampOld = xhatDamp; 
            xhatOld = xhat;
        
            xhatOldFinal = xhatFinal; % previous xhat 
            xhatIOldFinal = xhatIFinal;
            uOldFinal = uFinal;
    
            valOpt = [valOpt val];
        
            % residual coherence
            if k ~= 1  % C is not empty
                u = y - AI*xhatI - AC*xhatC;
            else
                u = y - AI*xhatI;
            end        
                
            uFinal = u;
        
            % Store variables for export
            xhatFinal = xhat;
            xhatIFinal = xhatI;
            xvarFinal = xvar;    
            zhatFinal = zhat;
            phatFinal = phat;
            pvarFinal = pvarRobust;
            zhat0Final = zhat0;
            zvar0Final = zvar0;
            shatFinal = shatDamp; 
            svarFinal = svarDamp; 

            if (it>nitMin) && (stop==false)
                if ( (norm(uOldFinal - uFinal) / norm(uFinal) < tol)...
                        && (abs(xhatIOldFinal - xhatIFinal) / abs(xhatIFinal) < tol) )
                    C = [C,I];              
                    stop = true;
                end
            end
        
            if (it >= nit) && (stop==false)
                C = [C,I];
                stop = true;
            end

            % Output nonlinear stage
            [zhat0, zvar0] = scaEstOut.estim(phat, pvarRobust); 
            shatNew = (1./pvarRobust).*(zhat0 - phat);
            svarNew = (1./pvarRobust).*(1 - min(zvar0./pvarRobust, zvarToPvarMax));
        
            % Increase stepsize, keeping within bounds
            step = min([stepIncr*max([step stepMin]) stepMax]);

        else
            % Automatically decrease stepMax (when opt.maxBadSteps<Inf)
            failCount = failCount + 1;
            if failCount > maxBadSteps
                failCount = 0;
                stepMax = max(stepMin,maxStepDecr*stepMax);
            end
        
            % Decrease stepsize, keeping within bounds
            step = min(max(stepMin, stepDecr*step),stepMax);
        
        end  % end of pass
        
        % Print results
        if (verbose)
            fprintf(1,'it=%3d |dx|/|x|=%12.4e\n', ...
                it, norm(xhatOldFinal(:) - xhatFinal(:)) / norm(xhatFinal(:)));
        end

        % Apply damping to shat, svar, and xhat
        if (it==1)
            if any(isnan(svarOld)) % if user didn't specify opt.svar0
                svarOld = svarNew;  % equivalent to making step=1
            end
            if any(isnan(xhatDampOld)) % if user didn't specify opt.xhatPrev0
                xhatDampOld = xhatOld; % equivalent to making step=1
            end
        end

        shatDamp = (1-step)*shatOld + step*shatNew;
        svarDamp = (1-step)*svarOld + step*svarNew;
        svarDamp = min(svarDamp, 1e5);
        xhatIDamp = (1-step)*xhatDampOld(I) + step*xhatOld(I);
        xhatDamp(I) = xhatIDamp;
        
        rvarI = 1./(AI2' * svarDamp);
        rvarIRobust = max(rvarI, xvarThresh);
        % Input linear stage
        rhatI = xhatIDamp + rvarIRobust*(AI' * shatDamp);

        % Input nonlinear stage
        % Compute mean and variance 
        [xhatI, xvarI, valIn] = scaEstIn.estim(rhatI, rvarI);
        valIn = sum( valIn(:) );
        
        xhat(I) = xhatI;
        xvar(I) = xvarI;
        
        if k ~= 1 % C is not empty
            if stop == false
                D = C;
            else
                D = C(1:end-1);
            end
            
            xhatCDamp = (1-step)*xhatDampOld(D) + step*xhatOld(D);
            xhatDamp(D) = xhatCDamp;
        
            % Step in rvar
            rvarC = 1./(AC2' * svarDamp);
            rvarCRobust = max(rvarC, xvarThresh);
            % Input linear stage
            rhatC = xhatCDamp + rvarCRobust.*(AC' * shatDamp);

            % Input nonlinear stage
            % Compute mean and variance 
            [xhatC, xvarC, valIn] = scaEstIn.estim(rhatC, rvarC);
            valIn = sum( valIn(:) );
            
            if stop == false
                xhat(C) = xhatC;
                xvar(C) = xvarC;
            else
                xhat(D) = xhatC;
                xvar(D) = xvarC;
            end
        end   
    end % main loop    
    
end  % end of pursuitDone

estFin.xhat = xhatFinal;
estFin.xvar = xvarFinal;
estFin.zhat = zhatFinal; 
estFin.phat = phatFinal;
estFin.pvar = pvarFinal;
estFin.zhat0 = zhat0Final;
estFin.zvar0 = zvar0Final;
estFin.shat = shatFinal;
estFin.svar = svarFinal;

estFin.nit = it;
estFin.C = C;

estHist.C = C;
