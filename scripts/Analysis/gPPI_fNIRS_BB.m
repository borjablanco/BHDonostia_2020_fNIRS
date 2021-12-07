
clear; close all;clc

% Add scripts path
path_scripts = ('/Users/borjablanco/Documents/GitHub/BHDonostia_2020_fNIRS/scripts/');
addpath(genpath(path_scripts))
% Load data
path_data = ('/Users/borjablanco/Documents/BCBL/Speech_data/Preprocessing/gamma_adjusted');

cd(path_data)
sub = dir('*preprocessed.mat');


for nsub = 1:length(sub)
    load(sub(nsub).name);
    % FW condition
    FW = data.s(:, 2);
    % BW condition
    BW = data.s(:, 6);
    
    sf = data.sf;
    time = linspace(0, length(data.OD)/sf/60, length(data.OD)); % in minutes
    
    %  Create FW_H and BW_H (task dependent regressors)
    % Take out artifacts and create FW_H and BW_H
    % One way to do it
    % FW_H = conv(FW, tbasis)
    % BW_H = conv(BW, tbasis)
    
    [FW_H, BW_H, lstInc, Am] = glm_model(data);
    
    % Solve GLM (whole brain version)
    % It keeps the betas for each seed
    % other_ch(:, 1) = yhbo_clean *beta_1 + FW_H * beta_2 + BW_H *beta_3 + PPI_FW * beta_4 + PPI_BW * beta_5 + e;
    ch = data.nchannels/2;
    betas_fw_hbo = zeros(ch, ch);
    betas_fw_hbr = zeros(ch, ch);
    betas_bw_hbo = zeros(ch, ch);
    betas_bw_hbr = zeros(ch, ch);
    for seed = 1:ch
        
        % Define seed data
        % With GSR filtering
        yhbo = data.GSR_oxy(:, seed);
        yhbr = data.GSR_deoxy(:, seed);
        
        % Censor motion artifacts
        yhbo_clean = yhbo(lstInc);
        yhbr_clean = yhbr(lstInc);
        
        % Define PPI terms
        PPI_FW_hbo = yhbo_clean .* FW_H; % Oxy hemoglobin
        PPI_FW_hbr = yhbr_clean .* FW_H; % Deoxy hemoglobin
        
        PPI_BW_hbo = yhbo_clean .* BW_H; % Oxy hemoglobin
        PPI_BW_hbr = yhbr_clean .* BW_H; % Deoxy hemoglobin
        
        % Create design matrix
        X_hbo = [yhbo_clean, Am, PPI_FW_hbo, PPI_BW_hbo];
        X_hbo_inv = pinv(X_hbo);
        
        X_hbr = [yhbr_clean, Am, PPI_FW_hbr, PPI_BW_hbo];
        X_hbr_inv = pinv(X_hbr);
        
        % Solve GLM and store betas for PPI terms
        betas_hbo = X_hbo_inv * data.GSR_oxy(lstInc, :);
        betas_hbr = X_hbr_inv * data.GSR_deoxy(lstInc, :);
        
        l_betas = size(betas_hbo,1);
        
        betas_fw_hbo(:,seed) = betas_hbo(l_betas-1,:);
        betas_fw_hbr(:,seed) = betas_hbr(l_betas-1,:);
        
        betas_bw_hbo(:,seed) = betas_hbo(l_betas,:);
        betas_bw_hbr(:,seed) = betas_hbr(l_betas,:);
        
    end
end


figure; 
subplot (221)
imagesc(betas_fw_hbo, [-1 1]); 
subplot (222)
imagesc(betas_bw_hbo, [-1 1]); 
subplot (223)
imagesc(betas_fw_hbr, [-1 1]); 
subplot (224)
imagesc(betas_bw_hbr, [-1 1]); 
colormap jet






