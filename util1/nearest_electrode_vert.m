function [vert_inds, coords] = nearest_electrode_vert(cortex, elecmatrix, matlab_order, yzplane)
%function [vert_inds, coords] = NEAREST_ELECTRODE_VERT(cortex, elecmatrix, matlab_order, yzplane)
%
% Find the vertex on a mesh that is closest to the given electrode
% coordinates
%
% Inputs:
%   cortex:       a struct containing cortex.tri and cortex.vert to construct the
%                 cortical surface mesh
%   elecmatrix:   nchans x 3 matrix of electrode coordinates in surface RAS
%                 space
%   matlab_order: 0 or 1, indicates whether to return 0-indexed or
%                 1-indexed vertices
%   yzplane:      0 or 1, whether to calculate distances from the yzplane only
%                 (disregarding left-right axis)
%
% Returns vertex indices and the new coordinates
%
% Written 2015 by Liberty Hamilton
%

if nargin<3
    matlab_order = 1;
end
if nargin<4
    yzplane = 0; % Use distance to the yz plane instead of using all coordinates
end

xyplane=0; % Doesn't work
% Find which hemisphere the electrodes are from using the x coordinate
if mean(elecmatrix(:,1)) < 0
    hem = 'lh';
else
    hem = 'rh';
end

cortex.vert
% Find distance between each electrode coordinate and the surface vertices
d=zeros(size(elecmatrix,1), size(cortex.vert,1));
for chan=1:size(elecmatrix,1)
    %fprintf(1,'Channel %d\n', chan);
    if yzplane
        d(chan,:)=sqrt((elecmatrix(chan,2)-cortex.vert(:,2)).^2+(elecmatrix(chan,3)-cortex.vert(:,3)).^2);
    elseif xyplane
        d(chan,:)=sqrt((elecmatrix(chan,1)-cortex.vert(:,1)).^2+(elecmatrix(chan,2)-cortex.vert(:,2)).^2);
    else
            d(chan,:)=sqrt((elecmatrix(chan,1)-cortex.vert(:,1)).^2+(elecmatrix(chan,2)-cortex.vert(:,2)).^2+(elecmatrix(chan,3)-cortex.vert(:,3)).^2);
    end
end
%keyboard;
[~,vert_inds]=min(d'); % Find the indices where the distance is minimized

coords = cortex.vert(vert_inds,:);

if matlab_order==0
    % Important!! Vertex indices are 0 indexed for freesurfer, so we have to
    % subtract 1 from them here!
    vert_inds = vert_inds - 1; % for freesurfer compatibility
end
