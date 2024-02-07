function [W,G] = NMF_convex(XX,K,W0,G0,Nit,X)
% function [W,G] = NMF_convex(XX,K,W0,G0,Nit,X)
%
%This does convex-NMF based on multiplicative update rules.
%It fits X = X*W*G' = F*G', i.e. cols of F are convex avgs of X.
%
%INPUTS:
% XX [P x P]: X'*X of data matrix X to factorize
% K [int]: number of clusters
% W0 [P x K]: initialize weight matrix (optional)
% G0 [K x P]: initialize indicator matrix (optional)
% Nit [int]: # of iterations (optional)
%
%OUTPUTS:
% W [P x K]: cluster weight matrix
% G [K x P]: cluster "indicator" matrix
%

P = size(XX,1);

if nargin<3 || isempty(W0)
    G = rand(P,K);
    % Multiply all rows in G by 1/sum(G) -- bsxfun does this faster than
    % using repmat
    W = bsxfun(@times,G,1./sum(G));
elseif nargin<4 || isempty(G0)
    W = W0;
    %This is based on Grindlay's code, following Ding et al. (2008)
    G0 = (W/(W'*W))';
    G = zeros(size(G0),'double');
    G(G0>0) = G0(G0>0) + ...
        0.2*ones(size(G0(G0>0)))*sum(G0(G0>0))./sum(G0(:)>0);
    G = max(G,10*eps)';
else
    W = W0;
    G = G0';
end

%Pos/neg parts of XX
XXp = .5*(abs(XX)+XX);
XXn = .5*(abs(XX)-XX);

if nargin<5 || isempty(Nit)
    Nit = 1000;
end
if nargin < 6 
    X = [];
end

meanerr = zeros(Nit,1);
for ii = 1:Nit
    XXpW = XXp*W;
    XXnW = XXn*W;
    GWt = G*W';
    G = G.*sqrt((XXpW + GWt*XXnW)./(XXnW + GWt*XXpW)); 
    G = max(G,eps); %This is from Li's code
    GtG = G'*G;
    W = W.*sqrt((XXp*G + XXnW*GtG)./(XXn*G + XXpW*GtG)); 
    W = max(W,eps);
    if ~isempty(X)
        meanerr(ii) = mean(sqrt(sum((X - X*W*G').^2))); 
    end
end
if ~isempty(X)
    figure(1);
    plot(meanerr,'.');
    xlabel('Iteration');
    ylabel('Mean error');
end
G = G';
