%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Debiasing:
%
%   Inputs:
%
%   Outputs:
%
%   BCBL, September 2018
%   Eneko Urunuela
%   e.urunuela@bcbl.eu

function [beta_debiased, s_zero, Y_fit] = Debiasing(input_signal, auc, r2only, mode, hrf)

% Default from "spfm.m"
signal          = input_signal.data;
TR              = input_signal.tr;
TE              = input_signal.te;
nscans          = size(signal,1)/length(TE);
nvoxels         = size(signal,2);
TE_mean         = mean(TE);
TE_norm         = TE/TE_mean;

% HRF matrix
X_hrf = hrf.hrf;
X_hrf_norm = hrf.norm;

%% Debiasing
if r2only
    beta_debiased_sum = zeros(nscans,nvoxels);
    beta_debiased_conv = zeros(nscans,nvoxels);
else
    beta_debiased_sum = zeros(2*nscans,nvoxels);
    beta_debiased_conv = zeros(2*nscans,nvoxels);
end

Y_fit_sum = zeros(size(signal));
Y_fit_conv = zeros(size(signal));

% Only columns (voxels) with non zero values need to be debiased
[~, col] = find(auc ~= 0);
voxelsToDebias = unique(col);
prog_old = 0;

for idxvox = 1 : length(voxelsToDebias)
    idx = find(auc(:,voxelsToDebias(idxvox)) ~= 0); % Indexes of non-zero coefficients
    
    if any(idx == nscans) 
        if length(idx) == 1 || idx(length(idx)-1) ~= (nscans - 1)
            idx(end) = [];
        end
    end
    
    if ~isempty(idx)
        X_sum = []; % X matrix of summed HRF signals
        X = X_hrf_norm(:,idx); % HRF columns of non-zero coefficients
        X_conv = []; % X matrix of convolved HRF signals
        flag_array = []; % Array with flags of consecutive non-zero coefficients
        pulse = zeros(nscans,1); % Array containing pulse to convolve with HRF
        
%         diff_old = 1; % Difference between non-zero coefficients indexes on previous iteration
        counter = 1; % Counts consecutive non-zero coefficients (used for flag_array)
        if length(idx) > 1 % If there is more than one coefficient
            if exist('diff_idx','var')
                clear diff_idx
            end
            for i = 1 : length(idx)-1
                diff_idx(i) = (idx(i+1) - idx(i));
            end
            
            sum_hrf_signal = X(:,1);
            pulse(idx(1)) = 1;
            
            for i = 1 : length(diff_idx)
                
                if diff_idx(i) > 1
                    X_sum = [X_sum sum_hrf_signal]; % Append summed hrf signal to X matrix
                    flag_array = [flag_array ones(1,counter)*size(X_sum,2)]; % Append flags of consecutive coefficients
                    temp = filter(spm_hrf(TR),1,pulse); % Convolve pulse signal with hrf
                    X_conv = [X_conv temp]; % Append convolved hrf signal to X matrix
                    
                    sum_hrf_signal = zeros(size(X(:,1))); % Resets sum of hrf signals array
                    pulse = zeros(nscans,1); % Resets pulse array
                    counter = 0; % Resets counter
                end
                
                sum_hrf_signal = sum_hrf_signal + X(:,i+1); % Sum both HRF signals
                pulse(idx(i+1)) = 1; % pulse
                counter = counter + 1; % Updates counter
            end
            
            X_sum = [X_sum sum_hrf_signal]; % Append summed hrf signal to X matrix
            flag_array = [flag_array ones(1,counter)*size(X_sum,2)]; % Append flags of consecutive coefficients
            temp = filter(spm_hrf(TR),1,pulse); % Convolve pulse signal with hrf
            X_conv = [X_conv temp]; % Append convolved hrf signal to X matrix

%             for i = 1 : length(idx)-1 % Iterate through coefficient indexes
%                 diff_idx = (idx(i+1) - idx(i)); % Difference between coefficient index and following one
%                 if i == 1 % First timepoint has a non-zero coefficient
%                     sum_hrf_signal = X(:,1);
%                     pulse(idx(1)) = 1; %auc(idx(1),idxvox);
%                 end
%                 if diff_idx == 1 % Coefficients are consecutive
%                     sum_hrf_signal = sum_hrf_signal + X(:,i+1); % Sum both HRF signals
%                     pulse(idx(i+1)) = 1; % pulse
%                 else % Coefficients are not consecutive
%                     if diff_old == 1 % Consecutive coefficients end here
%                         X_sum = [X_sum sum_hrf_signal]; % Append summed hrf signal to X matrix
%                         flag_array = [flag_array ones(1,counter)*size(X_sum,2)]; % Append flags of consecutive coefficients
%                     else % Previous coefficient was not consecutive either
%                         sum_hrf_signal = X(:,i+1); % No signal to sum, hrf stays the same
%                         X_sum = [X_sum sum_hrf_signal]; % Append summed hrf signal to X matrix
%                         flag_array = [flag_array ones(1,counter)*size(X_sum,2)]; % Append flags of consecutive coefficients
%                         pulse(idx(i+1)) = 1; %auc(idx(i+1),idxvox); % AUC value on index timepoint
%                     end
%                     
%                     temp = filter(spm_hrf(TR),1,pulse); % Convolve pulse signal with hrf
%                     X_conv = [X_conv temp]; % Append convolved hrf signal to X matrix
%                     
%                     sum_hrf_signal = zeros(size(X(:,1))); % Resets sum of hrf signals array
%                     pulse = zeros(nscans,1); % Resets pulse array
%                     counter = 0; % Resets counter
%                 end
%                 
%                 diff_old = diff_idx; % Saves current iteration value for next iteration
%                 counter = counter + 1; % Updates counter
%             end
%             if diff_old == 1 % In case there is a consecutive coefficient on last iteration
%                 X_sum = [X_sum sum_hrf_signal];
%                 temp = filter(spm_hrf(TR),1,pulse);
%                 X_conv = [X_conv temp];
%             end
%             if length(idx) > 1 % Updates flag_array with values from latest iterations
%                 flag_array = [flag_array ones(1,counter)*size(X_sum,2)];
%             elseif isempty(idx) % If there is no coefficient
%                 flag_array = 0;
%             else % If there is only one coefficient
%                 flag_array = 1;
%             end
            
            % Adds TE component to convolved hrf and normalizes the matrix
            X_conv_norm = [];
            for teidx = 1 : length(TE_norm)
                X_conv_norm = [X_conv_norm; -TE_norm(teidx)*X_conv];
            end
            
        else % If there is only one coefficient HRF stays as it is
            X_sum = X;
            pulse(idx(1)) = 1;
            X_conv = filter(spm_hrf(TR),1,pulse); % Convolve pulse signal with hrf
            X_conv_norm = [];
            for teidx = 1 : length(TE_norm)
                X_conv_norm = [X_conv_norm; -TE_norm(teidx)*X_conv];
            end
%             diff_old = 0;
            flag_array = 1;
            pulse = zeros(nscans,1); % Resets pulse array
        end
        
        X_conv_norm = bsxfun(@rdivide,-X_conv_norm,min(X_conv_norm,[],1));
        X_sum = bsxfun(@rdivide,-X_sum,min(X_sum,[],1));
        
        % Appends S0 term to both X matrices
        if ~r2only
            X_sum = [X_sum X_hrf_norm(:,nscans+1:end)];
            X_conv_norm = [X_conv_norm X_hrf_norm(:,nscans+1:end)];
        end

        X_sum_deb = zeros(size(X_hrf));
        X_conv_deb = zeros(size(X_hrf));
        for flag_idx = 1 : length(flag_array)
            if flag_array(flag_idx) ~= 0
                X_sum_deb(:,idx(flag_idx)) = X_sum(:,flag_array(flag_idx));
                X_conv_deb(:,idx(flag_idx)) = X_conv_norm(:,flag_array(flag_idx));
            else
                X_sum_deb(:,idx(flag_idx)) = zeros(size(X_sum_deb,1),1);
                X_conv_deb(:,idx(flag_idx)) = zeros(size(X_sum_deb,1),1);
            end
        end

        beta_debiased_sum(:,voxelsToDebias(idxvox)) = pinv(X_sum_deb'*X_sum_deb)*X_sum_deb'*signal(:,voxelsToDebias(idxvox));
        beta_debiased_conv(:,voxelsToDebias(idxvox)) = pinv(X_conv_deb'*X_conv_deb)*X_conv_deb'*signal(:,voxelsToDebias(idxvox));
        
        Y_fit_sum(:,voxelsToDebias(idxvox)) = X_sum_deb*beta_debiased_sum(:,voxelsToDebias(idxvox));
        Y_fit_conv(:,voxelsToDebias(idxvox)) = X_conv_deb*beta_debiased_conv(:,voxelsToDebias(idxvox));
        
    end
    prog = round(100*idxvox/length(voxelsToDebias));
    
    if prog > prog_old
        prog_old = prog;
        fprintf('Debiasing progress: %i %% \n',prog);
    end
end

if strcmp(mode, 'conv')
    beta_debiased = beta_debiased_conv;
    Y_fit = Y_fit_conv;
else
    beta_debiased = beta_debiased_sum;
    Y_fit = Y_fit_sum;
end

if ~r2only
    s_zero = beta_debiased(nscans+1:end,:);
    beta_debiased = beta_debiased(1:nscans,:);
else
    s_zero = [];
end


