
clear; close all;clc

% Add scripts path
%path_scripts = ('/Users/irenearrietasagredo/Desktop/BCBL/Brainhack/Brainhack_2020/Fork_brainhack_2020/BHDonostia_2020_fNIRS/scripts');
path_scripts = ('/Users/borjablanco/Downloads/BrainHack_project/BHDonostia_2020_fNIRS/scripts/');
addpath(genpath(path_scripts))
% Load data
%path_data = ('/Users/irenearrietasagredo/Desktop/BCBL/Brainhack/Brainhack_2020/Fork_brainhack_2020/BHDonostia_2020_fNIRS/data_files');
path_data = ('/Users/borjablanco/Downloads/BrainHack_project/BHDonostia_2020_fNIRS/data_files');

cd(path_data)
load('BH_data2.mat')

%% Objective

% To find PPI at different areas of the brain (Check out the literature on : https://github.com/brainhackorg/global2020/issues/37#issue-727168940)

% Y_i_H = seed_H * beta_1 + FW_H * beta_2 + BW_H * beta_3 + [seed_n *
% FW]_H * beta_4 + [seed_n * BW]_H * beta_5 + e

% Objective day 1
% 1 - without deconvolution

% Objective day 2
% 2 - with deconvolution

% Objective 3. Interpretation of results, play with deconvolution
% parameters, preprocessing parameters, seed choice parameters...
% Final design matrix


% Solve GLM
%%

% How to start? What do we already have?

% FW = yes 
% BW = yes 
% HRF = yes

% seed_H = Which to choose? A priori one from borja's article?

% seed_n = NO, need to devonvolute. I can first just do it without
% deconvolution.



%% Check visually how the signals look like

% Understand data. If any doubts, we can go through tomorrow in 'person'.


%% Equation variables: FW and BW. One way of doing it, without taking into account the motion artifacts and global signal regression.
% % FW condition
% FW = data.s(:, 2);
% % BW condition
% BW = data.s(:, 6);
% 
% sf = data.sf;
% time = linspace(0, length(data.OD)/sf/60, length(data.OD)); % in minutes
% 
% % Plot FW and BW
% figure; 
% plot(time, FW, 'blue'); hold on;
% plot(time, BW, 'red');
% xlabel('Minutes')
% ylabel('Activation yes/no')
% title('FW and BW psychological vectors')

%%  Create HRF: A modified gamma function convolved with a square-wave of duration T
% Code retrieved from Homer2
% tau = 0;
% sigma = 6;
% 
% trange = [-5 25];
% 
% % Takes the nT (inverse to sf) and creates a vector of ones which
% % introduces in the exponential for creating the cannonical HRF.
% t = data.time;
% nT = length(t);
% dt = t(2)-t(1);
% nPre = round(trange(1)/dt);
% nPost = round(trange(2)/dt);
% nTpts = size(data.OD,1);
% tHRF = (1*nPre*dt:dt:nPost*dt)';
% ntHRF=length(tHRF);  
% 
% tbasis = (exp(1)*(tHRF-tau).^2/sigma^2) .* exp( -(tHRF-tau).^2/sigma^2 );
% 
% % Make zero baseline values
% lstNeg = find(tHRF<0);
% tbasis(lstNeg,1) = 0;
% 
% figure
% plot(tHRF, tbasis)
% xlabel('Seconds')
% ylabel('Amplitude, normalized')
% title('Canonical HRF')

%%  Create FW_H and BW_H (task dependent regressors)
% PARTE DE CESAR

% Take out artifacts and create FW_H and BW_H
% One way to do it
% FW_H = conv(FW, tbasis)
% BW_H = conv(BW, tbasis)

[FW_H, BW_H, lstInc, Am] = glm_model(data);

%% PPI_FW & PPI_BW [seed_n * FW]_H
% PPI_FW_hbo = yhbo_clean .* FW_H; % Oxy hemoglobin
% PPI_FW_hbr = yhbr_clean .* FW_H; % Deoxy hemoglobin
% 
% PPI_BW_hbo = yhbo_clean .* BW_H; % Oxy hemoglobin
% PPI_BW_hbr = yhbr_clean .* BW_H; % Deoxy hemoglobin
% 
% time_clean = linspace(0, length(PPI_FW_hbo)/data.sf, length(PPI_FW_hbo));
% 
% figure; 
% subplot(2, 1, 1);
% plot(time_clean, PPI_FW_hbo, 'red'); hold on;
% plot(time_clean, PPI_FW_hbr, 'blue');
% plot(time_clean, FW(lstInc), 'black')
% xlabel('Time (minutes)')
% ylabel('Amplitude')
% title('PPI FW')
% 
% subplot(2, 1, 2)
% plot(time_clean, PPI_BW_hbo, 'red'); hold on;
% plot(time_clean, PPI_BW_hbr, 'blue');
% plot(time_clean, BW(lstInc), 'black')
% xlabel('Time (minutes)')
% ylabel('Amplitude')
% title('PPI BW')

%% Solve GLM (whole brain version)
% It keeps the betas for each seed
% other_ch(:, 1) = yhbo_clean *beta_1 + FW_H * beta_2 + BW_H *beta_3 + PPI_FW * beta_4 + PPI_BW * beta_5 + e;
ch = data.nchannels/2;
betas_fw_hbo = zeros(ch, ch);
betas_fw_hbr = zeros(ch, ch);
betas_bw_hbo = zeros(ch, ch);
betas_bw_hbr = zeros(ch, ch);

% Data for analysis
yhbo = data.GSR_oxy(lstInc, :);
yhbr = data.GSR_deoxy(lstInc, :);

for seed = 1:ch
    
    % Define seed data
    % With GSR filtering
    seed_hbo = yhbo(:, seed);
    seed_hbr = yhbr(:, seed);
    
    % Define PPI terms
    PPI_FW_hbo = seed_hbo .* FW_H; % FW - HbO
    PPI_FW_hbr = seed_hbr .* FW_H; % FW - HbR
    
    PPI_BW_hbo = seed_hbo .* BW_H; % BW - HbO
    PPI_BW_hbr = seed_hbr .* BW_H; % BW - HbR

    % Create design matrix
    X_hbo = [seed_hbo, Am, PPI_FW_hbo, PPI_BW_hbo];
    X_hbo_inv = pinv(X_hbo);
    
    X_hbr = [seed_hbr, Am, PPI_FW_hbr, PPI_BW_hbr];
    X_hbr_inv = pinv(X_hbr);
    
    % Solve GLM
    betas_hbo = X_hbo_inv * yhbo;
    betas_hbr = X_hbr_inv * yhbr;
    
    % Store betas PPI 
    l_betas = size(betas_hbo,1);
    
    betas_fw_hbo(:,seed) = betas_hbo(l_betas-1,:); 
    betas_fw_hbr(:,seed) = betas_hbr(l_betas-1,:); 
    
    betas_bw_hbo(:,seed) = betas_hbo(l_betas,:); 
    betas_bw_hbr(:,seed) = betas_hbr(l_betas,:); 
    
end
    
figure; 
subplot (221)
imagesc(betas_fw_hbo, [-1 1]); title('betas PPI FW HbO')
subplot (222)
imagesc(betas_bw_hbo, [-1 1]);  title('betas PPI BW HbO')
subplot (223)
imagesc(betas_fw_hbr, [-1 1]);  title('betas PPI FW HbR')
subplot (224)
imagesc(betas_bw_hbr, [-1 1]);  title('betas PPI BW HbR')
colormap jet

%%  Deconvolution EU

% Inputs
%   - input_signal: struct with data from the signal to be analysed.
%   - r2only: 0 for R2* only, 1 for model with S0.
%   - fista_params: struct with parameters needed for the FISTA.
%       - rho: rho value to weight L1 and L21 in the proximal method.
%       - lambda: lambda values of the proximal method.
%       - weights: weight of lambda values.
%       - itermax: maximum number of FISTA iterations.
%       - converged: convergence criteria to stop FISTA.
%       - proximal_mode: name of the proximal method to be used.
%   - stability: 1 for stability selection purposes, 0 otherwise.
%   - hrf: struct with
%   Outputs:
%   - betastruct: struct with deconvolved activity signal.
%       - beta: non debiased deconvolved activity signal.
%       - debiased: debiased deconvolved activity signal.
%   - Ystruct: struct with bold signal of deconvolved activity signal.
%       - fit: non debiased bold signal.
%       - debiased: debiased bold signal.

% Tiene que empezar en 0
% osea hacer de 0 a 25 no de -5 a 25
% hrf.hrf hrf sin normalizar
% hrf.norm hrf normalizada a uno

% Define parameters
input_signal.data = seed_hbo;
input_signal.te = 1;
r2only = [];
fista_params.rho =  [];
fista_params.lambda = [];
fista_params.weights = [];
fista_params.itermax = 100;
fista_params.converged = [];
fista_params.proximal_mode = [];
fista_params.lambda = [];
stability = 1;
hrf = [];


% Perform deconvolution
[betastruct, Ystruct, s_zero] = DeconvolutionAlgorithm(input_signal, r2only, fista_params, stability, hrf);

% Generate PPI regressors + convolution


% Solve GLM













