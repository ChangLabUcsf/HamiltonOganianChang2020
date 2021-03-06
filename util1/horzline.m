function [lh] = horzline(yloc,hax,color, linestyle, linewidth)
% function [lh] = vertline(xloc,hax,color, linestyle)
if nargin<1
    yloc =0;
end

if nargin<2 || isempty(hax)
    hax = gcf;
end

if nargin<3 || isempty(color)
    color =  'k';
end

if nargin<4 || isempty(linestyle)
    linestyle =  ':';
end


if nargin<5|| isempty(linewidth)
    linewidth =  1;
end
lh=line(get(gca, 'xlim'),[yloc yloc]', 'color', color, 'linestyle', linestyle,'linewidth',linewidth);