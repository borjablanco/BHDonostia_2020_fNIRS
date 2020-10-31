%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   StabilityPath:
%   Calculates the stability path and AUC of each timepoint for each of the
%   voxels.
%
%   Inputs:
%   - probs: matrix with the probabilities for each lambda value.
%   - lambda_values: array with values of lambda used through the
%   iterations.
%
%   Outputs:
%   - auc: matrix with AUC of each timepoint for each voxel.
%   - stability_cells: array of cells with stability path of each voxel on
%   each timepoint. Size: nscans x (nlambdas x nvoxels).
%
%   BCBL, August 2018
%   Eneko Urunuela
%   e.urunuela@bcbl.eu

function [auc] = StabilityPath(probs, lambda_values)

tic

fprintf('Calculating Stability Path... \n');

temp = zeros(size(probs));

for lambdaidx = 1 : length(lambda_values)
    temp(:,:,lambdaidx) = probs(:,:,lambdaidx).*lambda_values(lambdaidx)./sum(lambda_values); % Calculates AUC of probabilities.
end

auc = sum(temp,3);

elapsed = toc;

% fprintf('Total Stability Path time was %.2f seconds \n',round(elapsed,2));
fprintf('Stability Path calculated... \n');
end