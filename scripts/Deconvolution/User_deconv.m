%% Insert data

% Add scripts path
path_scripts = ('/Users/irenearrietasagredo/Desktop/BCBL/Brainhack/Brainhack_2020/Fork_brainhack_2020/BHDonostia_2020_fNIRS/scripts');
addpath(genpath(path_scripts))

% Load data
path_data = ('/Users/irenearrietasagredo/Desktop/BCBL/Brainhack/Brainhack_2020/Fork_brainhack_2020/BHDonostia_2020_fNIRS/data_files');
addpath(genpath(path_data))

load('BH_data2.mat');
% data_2 = load('BH_data2.mat');

[~, ~, lstInc, ~] = glm_model(data);

%%

% seed = 19;

%%
% With GSR filtering, prepare for input signal format
% Oxy
yhbo = data.GSR_oxy;
% Deoxy
yhbr = data.GSR_deoxy;

% Clean artifacts
yhbo_clean = yhbo(lstInc);
yhbr_clean = yhbr(lstInc);

%% Inputs:
%   - input_signal: struct with data from the input_signal to be analysed.
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


input_signal = yhbo_clean;

%-- r2only : changes in R2* (fMRI stuff), not needed in our case.
r2only = 1;

fista_params = struct();

%-- Fista_params : name of the proximal method to be used.
fista_params.proximal_mode = {};

% Play a bit with the amount, so that if it does not converge, it stops at
% some time. Number of surrogates that u create.
fista_params.itermax = 100;

% When the changes are small enough (below this threshold) (play a bit)
fista_params.converged = 10^-6;

%-- stabilitySelection :

% How many regularization parameters you will explore.
stabilitySelection.nlambdas = 25;

% It will calculate the maximum lambda value (if u go higher, the result
% will be zero).
stabilitySelection.maxLambdaPercentage = 0.9;

% It will calculate the minimum lambda value (if u go higher, the result
% will be zero). If you go lower, it may be too noisy. (play after)
stabilitySelection.maxLambdaPercentage = 0.1;

% Lambda logarithmic better, in order to be more detailed in the interest
% lambda values. 'It's more probably to find  a solution when lambda = 0,
% it'd be least squares solution. Therefore, interesting to go step by step in more detail around those small lambda values'.

stabilitySelection.isLambdaLog = 1;

% Number of surrogates (subsamples). Probability to calculate/K.

stabilitySelection.niter = 50;


%%  Create HRF

tau = 0;
sigma = 6;

trange = [0 25];

% Takes the nT (inverse to sf) and creates a vector of ones which
% introduces in the exponential for creating the cannonical HRF.

t = data.time;
nT = length(t);
dt = t(2)-t(1);
nPre = round(trange(1)/dt);
nPost = round(trange(2)/dt);
nTpts = size(data.OD,1);
tHRF = (1*nPre*dt:dt:nPost*dt)';
ntHRF=length(tHRF);

tbasis = (exp(1)*(tHRF-tau).^2/sigma^2) .* exp( -(tHRF-tau).^2/sigma^2 );

% Make zero baseline values
lstNeg = find(tHRF<0);
tbasis(lstNeg,1) = 0;

% figure;
% plot(tHRF, tbasis)
% xlabel('Seconds')
% ylabel('Amplitude, normalized')
% title('Canonical HRF')

%% Create covariance matrix

tbasis_temp = tbasis;

l = length(yhbo_clean);
l_t = length(tbasis_temp);

H = zeros(l, l);

% quantity of zeros on top
q_0 = 0;
cut = l_t;

for k = 1:l
    term = q_0 + l_t;
    if (term>l) % If it already has to be cut
        cut = cut - 1;
        tbasis_temp = tbasis_temp([1:cut])
    end
    H([k:(k+cut-1)], k) = tbasis_temp;
    q_0 = q_0 + 1;
end


%%

hrf = H;

[auc, lambdarange, lambda_probs] = StabilitySelection(input_signal, r2only, fista_params, stabilitySelection, hrf)












