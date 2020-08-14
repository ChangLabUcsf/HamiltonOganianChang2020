function plotstd(data_mean, data_std, plot_color, xdata)%, fill_color)
% function PLOTSTD(data_mean, data_std, plot_color, xdata)%, fill_color)
%
% Plots mean and standard error/deviation (whatever is supplied) with
% transparent fill and fill color of your choice.  
% 
% xdata contains the values for the x axis
%

if nargin<3
    plot_color = 'k';
end
if nargin<4
    xdata = [1:length(data_mean)];
end
if size(data_mean,2) == 1
    data_mean = data_mean';
end
if size(data_std,2) == 1
    data_std = data_std';
end

stdfill = [data_mean+data_std fliplr(data_mean-data_std)];
%ff=fill(xdata, stdfill, fill_color);
ff=fill([xdata fliplr(xdata)], stdfill, plot_color);
set(ff,'edgecolor','none');
set(ff,'facealpha',0.5,'edgecolor','none');
hold on;

plot(xdata, data_mean, 'color', plot_color);