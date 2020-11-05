
clear; close all;clc

% Add scripts path
path_scripts = ('/Users/borjablanco/Downloads/BrainHack_project/BHDonostia_2020_fNIRS/scripts');
addpath(genpath(path_scripts))
% Load data
path_data = ('/Users/borjablanco/Downloads/BrainHack_project/BHDonostia_2020_fNIRS/data_files');
cd(path_data)
load('BH_data2.mat')

plot(data.OD); hold on
plot(data.tInc_auto)

% Create signals of interest
% --> Remember analysis can be done in HbO, HbR data:
% with GSR (data.GSR_oxy)and without GSR (data.conc)

% Define seed signal
seed = 18; % choose seed channel
yhbo = data.GSR_oxy(:, seed);
yhbr = data.GSR_deoxy(:, seed);

% Create task dependent regressors
[tdr_fw, tdr_bw] = glm_model(data);

% Create psychophysiological interaction terms


% 1 - without deconvolution


% 2 - with deconvolution


% Final design matrix


% Solve GLM












