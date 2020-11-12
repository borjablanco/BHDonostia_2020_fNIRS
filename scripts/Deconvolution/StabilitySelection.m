%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   StabilitySelection:
%   Makes different iterations to then calculate the stability path and AUC
%   of each timepoint for each of the voxels.
%
%   Inputs:
%   - input_signal: struct with data from the signal to be analysed.
%   - r2only: 0 for R2* only, 1 for model with S0.
%   - fista_params: struct with parameters needed for the FISTA.
%       - rho: rho value to weight L1 and L21 in the proximal method.
%       - lambda: lambda values of the proximal method.
%       - weights: weight of lambda values.
%       - itermax: maximum number of FISTA iterations.
%       - converged: convergence criteria to stop FISTA.
%       - proximal_mode: name of the proximal method to be used.
%   - stabilitySelection:
%       - nlambdas: number of lambdas to build the stability path on.
%       - maxLambdaPercentage: percentage of maximum lambda to set as
%       highest lambda value.
%       - minLambdaPercentage: percentage of maximum lambda to set as
%       lowest lambda value.
%       - isLambdaLog: 1 if lambda scale is log, 0 if it is linear.
%       - niter: number of iterations to run on each lambda.
%
%   Outputs:
%   - auc: matrix with AUC of each timepoint for each voxel.
%   - stability_cells: array of cells with stability path of each voxel on
%   each timepoint. Size: nscans x (nlambdas x nvoxels).
%
%   BCBL, August 2018
%   Eneko Urunuela
%   e.urunuela@bcbl.eu

function [auc, lambdarange, lambda_probs] = StabilitySelection(input_signal, r2only, fista_params, stabilitySelection, hrf)

% Argument checking
if (~isfield(stabilitySelection,'niter')) || isempty(stabilitySelection.niter); stabilitySelection.niter  = 100; end
if (~isfield(stabilitySelection,'nlambdas')) || isempty(stabilitySelection.nlambdas); stabilitySelection.nlambdas  = 50; end
if (~isfield(stabilitySelection,'maxpercentage')) || isempty(stabilitySelection.maxpercentage)
    stabilitySelection.maxpercentage  = 0.95;
end
if (~isfield(stabilitySelection,'minpercentage')) || isempty(stabilitySelection.minpercentage)
    stabilitySelection.minpercentage  = 0.05;
end
if (~isfield(stabilitySelection,'mode')) || isempty(stabilitySelection.mode); stabilitySelection.mode  = 1; end
if (~isfield(stabilitySelection,'lambdarange')) || isempty(stabilitySelection.lambdarange); stabilitySelection.lambdarange  = []; end

%% Parameterization for Regularization path

stabilityiterations = stabilitySelection.niter;
signal          = input_signal.data;
TE              = input_signal.te;
nscans          = size(signal,1)/length(TE);
nvoxels         = size(signal,2);

% HRF matrix
X_tilde = hrf.norm;
X_tilde_trans = X_tilde';

%% Lambda values parameterization
if isempty(stabilitySelection.lambdarange)
    maxlambda = round(max(max(abs(X_tilde_trans*signal))));
    
    nlambdas            = stabilitySelection.nlambdas;
    maxLambdaPercentage = stabilitySelection.maxpercentage;
    minLambdaPercentage = stabilitySelection.minpercentage;
    isLambdaLog         = stabilitySelection.mode;
    
    switch isLambdaLog
        case 0 % Linearly spaced values of lambda
            lambdarange = fliplr(linspace(minLambdaPercentage*maxlambda,maxLambdaPercentage*maxlambda,nlambdas));
        case 1 % Log spaced values of lambda
            lambdarange = fliplr(logspace(log10(minLambdaPercentage*maxlambda),log10(maxLambdaPercentage*maxlambda),nlambdas));
    end
else
    lambdarange = stabilitySelection.lambdarange;
    nlambdas = length(lambdarange);
end

%% Clearing variables
clear maxLambdaPercentage minLambdaPercentage isLambdaLog X_tilde X_tilde_trans TE_mean dt TE hrf_p

%% Iterations for Stability Path
if r2only
    lambda_probs = zeros(nscans,nvoxels,nlambdas);
else
    lambda_probs = zeros(2*nscans,nvoxels,nlambdas);
end

% parpool('ips_base',25);
parpool('local',2);
parallelpool = gcp;

fprintf('Stability selection progress:\n');
% fprintf(['\n' repmat('.',1,nlambdas) '\n\n']);
parfor_progress(nlambdas);

parfor lambdaidx = 1 : nlambdas
    temp_params = struct();
    temp_params = fista_params;
    temp_params.lambda = lambdarange(lambdaidx); % Proximal operator parameter
    
    if r2only
        temp = false(nscans,nvoxels,stabilityiterations);
    else
        temp = false(nscans*2,nvoxels,stabilityiterations);
    end
    
    for iter = 1 : stabilityiterations
        
        % Runs algorithm
        [beta, ~] = DeconvolutionAlgorithm(input_signal, r2only, temp_params, 1, hrf);
        
        temp(:,:,iter) = logical(beta);
    end
    
    prob = sum(int8(temp),3)/stabilityiterations;
    
    lambda_probs(:,:,lambdaidx) = prob;
    parfor_progress;
%     fprintf('\b|\n');
end
parfor_progress(0);

delete(parallelpool);

%% Calculates AUC

[auc] = StabilityPath(lambda_probs(1:nscans,:,:), lambdarange); % Only on R2*

end
