% Thanks for gamp software package
% http://sourceforge.net/projects/gampmatlab/
%
% This code is open to anyone interested in using it.
% Yongjie Luo, E.E., UESTC
% yongjie.luo@hotmail.com
% 2014.5

classdef MpGAMPOpt
    % Options for the matched pursuit GAMP optimizer.
    
    properties
        
        %***** General Options
        
        % Print progress
        verbose = false;
                
        %Number of iterations
        nit = 1200;
        
        %Minimum number of iterations- sometimes useful with warm starting
        %approaches
        nitMin = 50; %0 for no effect
        
        %Specify a convergence tolerance. Iterations are terminated when
        %the norm of the differnece in two iterates divided by the norm of
        %the current iterate falls below this threshold. Set to -1 to
        %disable
        tol = 1e-8;
        
        %Error function. This is a function handle that acts on Zhat to
        %compute an NMSE error.
        error_function = @(q) inf;
        
        %Error functions for X
        error_functionX = @(q) inf;
                
        % If desired, a custom stopping criterion can be used to terminate
        % GAMP and/or display useful quantities during runtime.  
        % There are two versions:
        stopFcn = [];  
        % called as stop = stopFcn(val, xhat, xhatPrev, Axhat)
        %   val: the current objective function value
        %   xhat, xhatPrev: the estimate of X at the current, previous iteration respectivelly.
        %   Axhat: the current transform vector (zhat(k) = A*xhat(k)). 
        % A "true" return value indicates GAMP should terminate.
        stopFcn2 = [];
        % called as:  stop = stopFcn2(gampstate);
        %   gampstate: a structure whose fields describe the current state.
        % A stop value of "true" indicates GAMP should terminate.
        % For full details, see gampEst.m.        
        
        %***** Initialization
        
        %Provide initial guesses for xhat0,xvar0,shat0. If these are set to
        %empty, then the appropriate method of the input estimator is
        %called to initialize these values. This functionality is useful
        %for warm starting the algorithm when not providing the full state.
        xhat0 = [];
        xvar0 = [];
        shat0 = [];
        
        
        % The following are used for warm-starting the algorithm:
        svar0 = [];
        pvarOld0 = [];
        rvarOld0 = [];
        zvarOld0 = [];
        xhatDampOld0 = [];

        
        %***** Step Size Control
        
        %Logical flag to include a step size in the pvar/zvar calculation.
        %This momentum term often improves numerical performance. On by
        %defualt.
        pvarStep = true;
        
        %Initial step size, or fixed size for non-adaptive steps
        %step = 0.05;
                
        %***** Variance Control
                
        %Minimum variances. See code for details of use.
        pvarMin = 1e-12;
        rvarMin = 1e-12;
        zvarToPvarMax = 0.99;   % prevents negative svar, should be near 1
        
        %Variance threshold for rvar and qvar, set large for no effect
        varThresh = 1e6;
       
        % History Interval: save memory by decimating the saved history
        % e.g. histIntvl=20 will save iteration 20,40,...
        histIntvl = 1; % defaults to save every iteration
               
        
        % Initial stepsize
        step = 1;
        
        % Minimum stepsize.  If stepsize is initialized below this value,
        % then the iteration will be automatically considered successful.
        stepMin = 0.01;
        
        % Maximum allowed stepsize
        stepMax = 1;
        
        % Multiplicative stepsize increase, when successful
        stepIncr = 1.1;
        
        % Multiplicative stepsize decrease, when unsuccessful
        stepDecr = 0.5;
        
        %Maximum number of allowed failed steps before we decrease stepMax,
        %inf for no effect; 10 is a decent value to play with
        maxBadSteps = inf;
        
        %Amount to decrease stepMax after maxBadSteps failed steps, 1 for
        %no effect
        maxStepDecr = 0.8;
        
        % Iterations are termined when the stepsize becomes smaller
        % than this value. Set to -1 to disable
        stepTol = 1e-10;
        
        
        % Enable adaptive stepsize
        adaptStep = true;
        
        % Disable computing the cost via Bethe free energy
        adaptStepBethe = false;
        
        %Create a window for the adaptive stepsize test. Setting this to
        %zero causes it to have no effect. For postive integer values,
        %creates a moving window of this length when checking the stepsize
        %acceptance criteria. The new value is only required to be better
        %than the worst in this window, i.e. the stepsize criteria is not
        %required to monotonicaly increase. As with other algorithms, this
        %modification tends to improve convergence speed
        stepWindow = 20;  %50
        
        %Set to true to use a Barzilai Borwein type rule to choose a new
        %stepsize after each succesful step. Otherwise, the stepsize is set to
        %the previous successful stepsize
        bbStep = false;
        
    end
    
    methods
        
        % Constructor with default options
        function opt = MpGAMPOpt(varargin)
            if nargin == 0
                % No custom parameters values, thus create default object
                return
            elseif mod(nargin, 2) == 0
                % User is providing property/value pairs
                names = fieldnames(opt);    % Get names of class properties
                
                % Iterate through every property/value pair, assigning
                % user-specified values.  Note that the matching is NOT
                % case-sensitive
                for i = 1 : 2 : nargin - 1
                    if any(strcmpi(varargin{i}, names))
                        % User has specified a valid property
                        propName = names(strcmpi(varargin{i}, names));
                        opt.(propName{1}) = varargin{i+1};
                    else
                        % Unrecognized property name
                        error('MpGAMPOpt: %s is an unrecognized option', ...
                            num2str(varargin{i}));
                    end
                end
                return
            else
                error(['The MpGAMPOpt constructor requires arguments ' ...
                    'be provided in pairs, e.g., MpGAMPOpt(''verbose'',' ...
                    ' false, ''nit'', 50)'])
            end
        end
    end
    
end
