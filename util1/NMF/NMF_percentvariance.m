function [pve] = NMF_percentvariance(X, FLs, W)
% function [pve] = NMF_percentvariance(X, FLs, W)
%
% Calculates percent variance explained using convex NMF factors
% where X ~ X*W*FLs
%
% Expects FLs and W as cell arrays with length equal to the number of
% cluster levels you try fitting
%
% 2015 Liberty Hamilton
%

% Initialize percent variance explained
pve = zeros(length(FLs), 1);

% Loop through all cluster levels
for sb=1:length(FLs)
    try
    % Get reconstruction of the data using NMF factors
    Xnew = X*W{sb}*FLs{sb};
    
    % This is using the residual sum of squares over the total sum of
    % squares.  Percent variance explained is 1 minus this.
    % See https://en.wikipedia.org/wiki/Coefficient_of_determination
    % The mean across all X is 0 because of z-scoring, so denominator is
    % simplified
    pve(sb)= 1-sum((X(:)-Xnew(:)).^2,1)./sum(X(:).^2,1);
    catch
        pve(sb) = 0;
    end
end