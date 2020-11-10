
%clear; close all;clc

% Add scripts path
path_scripts = ('/Users/irenearrietasagredo/Desktop/BCBL/Brainhack/Brainhack_2020/Fork_brainhack_2020/BHDonostia_2020_fNIRS/scripts');
addpath(genpath(path_scripts))
% Load data
path_data = ('/Users/irenearrietasagredo/Desktop/BCBL/Brainhack/Brainhack_2020/Fork_brainhack_2020/BHDonostia_2020_fNIRS/data_files');
cd(path_data)
load('BH_data1.mat')

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
% FW
FW = data.s(:, 2)


% BW
BW = data.s(:, 6)


sf = data.sf;
time = linspace(0, length(data.OD)/sf/60, length(data.OD))
xlabel('Minutes')

% Plot FW and BW
figure; 
plot(time, FW, 'blue'); hold on;
plot(time, BW, 'red');
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

%%  Create FW_H and BW_H (task dependent regressors)
% Take out artifacts and create FW_H and BW_H

% One way to do it
% FW_H = conv(FW, tbasis)
% BW_H = conv(BW, tbasis)

[FW_H, BW_H, lstInc, Am] = glm_model(data);

%% Seed signal
 
% Choose seed, by channel.

% Define seed signal
seed = 19; % choose seed channel

%% With and without GSR filtering

% Without filtering by the global signal regression
% This is oxy - better data.filt
conc_yhbo = data.conc(:, 1, seed);
% This is deoxy
conc_yhbr = data.conc(:, 2, seed);

% Clean artifacts
conc_yhbo_clean = conc_yhbo(lstInc);
conc_yhbr_clean = conc_yhbr(lstInc);


% With GSR filtering
% Oxy
yhbo = data.GSR_oxy(:, seed);
% Deoxy
yhbr = data.GSR_deoxy(:, seed);

% Clean artifacts
yhbo_clean = yhbo(lstInc)
yhbr_clean = yhbr(lstInc)

% Visual inspection

figure; 
subplot(4, 1, 1)
plot(conc_yhbo, 'blue'); hold on;
plot(conc_yhbr, 'red');hold on; 
title('no gsr, no clean')

subplot(4, 1, 2)
plot(conc_yhbo_clean, 'blue'); hold on; 
plot(conc_yhbr_clean, 'red');hold on; 
title('no gsr, clean')
 
subplot(4, 1, 3)
plot(yhbo, 'blue'); hold on;
plot(yhbr, 'red');hold on; 
title('gsr, no clean')

subplot(4, 1, 4)
plot(yhbo_clean, 'blue'); hold on; 
plot(yhbr_clean, 'red');
title('gsr, clean')



%% PPI_FW & PPI_BW [seed_n * FW]_H

PPI_FW = yhbo_clean .*  FW_H; % Oxy hemoglobin
PPI_FW_r = yhbr_clean .*  FW_H; % Deoxy hemoglobin

PPI_BW = yhbo_clean .*  BW_H; % Oxy hemoglobin
PPI_BW_r = yhbr_clean .*  BW_H; % Deoxy hemoglobin

time_clean = linspace(0, length(PPI_FW)/data.sf, length(PPI_FW))

figure; 
subplot(2, 1, 1);
plot(time_clean, PPI_FW, 'red'); hold on;
plot(time_clean, PPI_FW_r, 'blue'); hold on;
xlabel('Time (minutes)')
ylabel('Amplitude')
title('PPI FW')

subplot(2, 1, 2)
plot(time_clean, PPI_BW, 'red'); hold on;
plot(time_clean, PPI_BW_r, 'blue'); hold on;
xlabel('Time (minutes)')
ylabel('Amplitude')
title('PPI BW')


%% Solve the optimization problem through GLM, minimizing the error term.


other_ch = data.GSR_oxy(lstInc, :);
other_ch(:, seed) = zeros(1, length(other_ch(:, seed)));


% other_ch(:, 1) = yhbo_clean *beta_1 + FW_H * beta_2 + BW_H *beta_3 + PPI_FW * beta_4 + PPI_BW * beta_5 + e;

% Design matrix
figure;
X = [yhbo_clean, Am, PPI_FW, PPI_BW];
X_inv = pinv(X);

betas = X_inv * other_ch;

% [5 *15546] * [15646*23]

%% 
figure;
imagesc(betas([5, 6], :), [-1, 1])
colormap(jet)

yticks([1, 2])
yticklabels({'PPI FW','PPI BW'})
xticks([1:24])
title('Betas after Least Squares')


%% 
% Create psychophysiological interaction terms


% 1 - without deconvolution


% 2 - with deconvolution


% Final design matrix


% Solve GLM












