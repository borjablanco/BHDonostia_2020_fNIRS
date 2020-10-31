function x = proximal_operator_lasso(y,lambda,weights)
N = length(y);
if (~exist('weights','var')) || isempty(weights); weights = ones(size(y,1),size(y,2)); end
weights = reshape(weights,size(y,1),size(y,2));
x = y.*max(zeros(size(y,1),size(y,2)), (1-(weights*lambda./abs(y))));
x(abs(x) < eps)=0;
end
