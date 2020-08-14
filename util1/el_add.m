function [h] = el_add(els, varargin)
% function [h] = EL_ADD(els,varargin)
%
% Optional arguments:
%   'color': string, or [electrodes x 3] matrix of colors
%   'msize': marker size (default = 6)
%   'mtype': marker type (default = 'o')
%   'edgecol': edge color, (default = 'none')
%   'numbers': whether to plot numbers next to the electrode markers or not
%   (default = [], no numbers)
%
% Edited by Liberty Hamilton 2015
% Original code (C) 2006  K.J. Miller

p = inputParser;
p.addOptional('color','r');
p.addOptional('msize',6);
p.addOptional('mtype','o');
p.addOptional('edgecol','none');
p.addOptional('numbers',[]);
p.addOptional('LineWidth',1);

p.parse(varargin{:});

elcol = p.Results.color;
msize = p.Results.msize;
mtype = p.Results.mtype;
edgecol = p.Results.edgecol;
numbers = p.Results.numbers;
linewidth = p.Results.LineWidth;

%hold on, plot3(els(:,1),els(:,2),els(:,3),'.','Color', elcol,'MarkerSize',msize)
hold on;

if ~isstr(elcol)
    if size(elcol,1)>1 && size(elcol,2)==3
        for i=1:size(els,1)
            % For plosive, vowel, nasal, fricative, tonic
            if strcmp(mtype, 'P') ||  strcmp(mtype, 'V') ||  strcmp(mtype, 'N') ||  strcmp(mtype, 'F') || strcmp(mtype,'T')
                h(i)=text(els(i,1),els(i,2),els(i,3), mtype);
                set(h(i), 'color', elcol(i,:),'fontsize',msize);
            else
                h(i)=plot3(els(i,1),els(i,2),els(i,3), 'Marker', mtype,'MarkerFaceColor', elcol(i,:),'MarkerSize',msize,'MarkerEdgeColor',edgecol,'LineWidth',linewidth,'Linestyle','none');
            end
        end
    else
        % For plosive, vowel, nasal, fricative, tonic
        if strcmp(mtype, 'P') ||  strcmp(mtype, 'V') ||  strcmp(mtype, 'N') ||  strcmp(mtype, 'F')|| strcmp(mtype,'T')
            h=text(els(:,1),els(:,2),els(:,3), mtype);
            set(h, 'color', elcol,'fontsize', msize, 'fontweight','bold');
        else
            h=plot3(els(:,1),els(:,2),els(:,3),'Marker',  mtype,'MarkerFaceColor', elcol,'MarkerSize',msize,'MarkerEdgeColor',edgecol,'LineWidth',linewidth,'Linestyle','none');
        end
    end
else
    h=plot3(els(:,1),els(:,2),els(:,3), 'Marker', mtype,'MarkerFaceColor', elcol,'MarkerSize',msize,'MarkerEdgeColor',edgecol,'LineWidth',linewidth,'Linestyle','none');
    
end

tt=[];
if ~isempty(numbers)
    for i=1:size(els,1)
        tt(i)=text(els(i,1)-1,els(i,2)-1,els(i,3)-1,num2str(numbers(i)));
        set(tt(i),'color','k','fontsize',12);
    end
end
for i=1:length(tt)
    uistack(tt(i), 'top');
end