function print_quality_fig(fighandle, fig_out, ftsize, wd, ht, unit_type, figtype)
% print_quality_fig(fighandle, fig_out, ftsize, wd, ht, unit_type)
%
% Inputs:
%   fighandle: the figure handle.  For example, if you have h=figure(1);, pass fighandle=h
%   fig_out: the file name of your output (default saves to a color .eps file, use .eps extension)
%   ftsize: the font size for labels (8 is a good choice for journal articles)
%   wd: the width of the figure
%   ht: the height of the figure
%   unit_type: units for width and height ('centimeters')
%
% Example usage:
%   fighandle = figure(1);
%   plot(randn(100,1)); % Plot some random data
%   fname = '/path/to/myfile.eps'; % Your output file name
%   fontsize = 8;
%   wd = 6;
%   ht = 3;
%   unit_type = 'centimeters';
%   figtype = 'epsc'; % Need the 'c' for color vector figures
%   print_quality_fig(fighandle, fname, fontsize, wd, ht, unit_type, figtype);
%
% Liberty Hamilton 2014
%

if nargin<3
    ftsize = 8;
end
if nargin<4
    wd = 8;
end
if nargin < 5
    ht = 5;
end
if nargin < 6
    unit_type = 'inches';
end
if nargin < 7
    figtype = 'epsc';
end
figure(fighandle);

% set font size and name
set(gca,'fontsize',ftsize,'fontname','Helvetica');
h=get(gca,'title');
set(h,'fontsize',ftsize,'fontname','Helvetica');
h=get(gca,'xlabel');
set(h,'fontsize',ftsize,'fontname','Helvetica');
h=get(gca,'ylabel');
set(h,'fontsize',ftsize,'fontname','Helvetica');

% turn off the box around the plot and put the tick marks on the outside
set(gca,'tickdir','out','box','off');

% make the figure the correct size
set(gcf, 'PaperPosition', [0 0 wd ht], 'PaperUnits', unit_type, 'PaperSize', [wd ht]);

% to make subplots work
figchild = get(gcf,'Children');
if ~isempty(figchild)
    for hndls = 1:length(figchild)
        try
            set(figchild(hndls), 'fontsize', ftsize, 'tickdir', 'out', 'box', 'off');
        catch
        end
    end
end

%whitebg(gcf,'k');
%set(gcf,'color','k');
if strcmp(figtype,'epsc') || strcmp(figtype, 'eps')
    set(gcf,'renderer','painters'); % For vector graphics
end
saveas(fighandle, fig_out, figtype);