
clear; close all;clc

% Add scripts path
path_scripts = ('/Users/irenearrietasagredo/Desktop/BCBL/Brainhack/BrainHack_2020/BHDonostia_2020_fNIRS/scripts');
addpath(genpath(path_scripts))
% Load data
path_data = ('/Users/irenearrietasagredo/Desktop/BCBL/Brainhack/BrainHack_2020/BHDonostia_2020_fNIRS/data_files');
cd(path_data)
load('BH_data2.mat')

%% Objective

% To find PPI at different areas of the brain.

% Y_i_H = seed_H * beta_1 + FW_H * beta_2 + BW_H * beta_3 + [seed_n *
% FW]_H * beta_4 + [seed_n * BW]_H * beta_5 + e


% 1 - without deconvolution


% 2 - with deconvolution


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

% Understand data 

% Optical density 
%plot(data.OD); hold on
% Movement parts
%plot(time, data.tInc_auto)

sf = data.sf;
time = linspace(0, length(data.OD)/sf/60, length(data.OD))


%% Equation variables: FW and BW
% FW
FW = data.s(:, 2)
FW_H = conv(FW, tbasis)

figure; 
plot(FW_H)

% BW
BW = data.s(:, 6)

sf = data.sf;
time = linspace(0, length(data.OD)/sf/60, length(data.OD))
xlabel('Minutes')

% Plot FW and BW
figure; 
plot(time, FW, 'blue'); hold on;
plot(time, BW, 'red'); hold on;
xlabel('Minutes')
ylabel('Activation yes/no')
title('FW and BW psychological vectors')

%%  Create HRF

tau = 0;
sigma = 6;

trange = [-5 25];

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

plot(tHRF, tbasis)
xlabel('Seconds')
ylabel('Amplitude, normalized')
title('Canonical HRF')
%% Create FW_H and BW_H

%% Seed signal
 
% Choose seed, by channel.

% Define seed signal
seed = 19; % choose seed channel

%% Filtering depth

% Without filtering by the global signal regression

% y_conc_hbo = data.conc(:, 1, :);
% y_conc_hbR = data.conc(:, 2, :);

% This is oxy
figure; 
plot(data.conc(:, 1, seed), 'blue');

% This is deoxy
% plot(data.conc(:, 2, seed), 'red'); hold on;


% Filtering by global signal regression
figure;
yhbo = data.GSR_oxy(:, seed);



plot(yhbo,'red')

yhbr = data.GSR_deoxy(:, seed);


%% 
% Create task dependent regressors
[tdr_fw, tdr_bw, lstInc] = glm_model(data);

% Keep only good time_points
yhbo_clean = yhbo(lstInc)

figure; 
plot(yhbo_clean); hold on; plot(yhbo, 'red')


%% [seed_n * FW]_H

PPI_FW = conv(yhbo_clean,  tdr_fw, 'same') % Oxy hemoglobin
PPI_BW = conv(data.conc(:, 1, seed), BW_H) 

figure; 
plot(PPI_FW)

%% 
% Create psychophysiological interaction terms


% 1 - without deconvolution


% 2 - with deconvolution


% Final design matrix


% Solve GLM












