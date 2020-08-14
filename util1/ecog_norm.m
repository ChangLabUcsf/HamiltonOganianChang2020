function [v2] = ecog_norm(v, normdim)
if nargin<2
    [~,normdim] = max(size(v));
end
allmax =max(abs(v),[],normdim); 
allmax (allmax ==0)=1;
v2 = v./allmax;

