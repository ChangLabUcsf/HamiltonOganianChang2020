function [c_h, vert_values] = ctmr_gauss_plot_2Dcmap(cortex,electrodes,weights,hemi,do_lighting, gsp)
% function [c_h]=ctmr_gauss_plot_2Dcmap(cortex,electrodes,weights,hemi,do_lighting)
%
% projects electrode locations onto their cortical spots in the
% left hemisphere and plots about them using a gaussian kernel
% for only cortex use:
% ctmr_gauss_plot(cortex,[0 0 0],0,'lh',1);
%
% You can also just call ctmr_gauss_plot(cortex); 
% and it will assume 
%
%     Copyright (C) 2009  K.J. Miller & D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   Version 1.1.0, released 26-11-2009
%
% Edited 2016 by Liberty Hamilton

if nargin<2
    electrodes = [0 0 0];
end
if nargin<3
    weights = 0;
end
if nargin<4
    hemi = 'lh';
end
if nargin<5
    do_lighting=1;
end
if nargin<6
    gsp = 20;
end

%load in colormap
% load('loc_colormap')
load('loc_colormap_thresh')
% load('BlWhRdYl_colormap')
% load('BlGyOrCp_colormap')
cm = flipud(cbrewer('div','RdBu',256));

brain=cortex.vert;
v='l';

if length(weights)~=length(electrodes(:,1))
    error('you sent a different number of weights than electrodes (perhaps a whole matrix instead of vector)')
end

%gaussian "cortical" spreading parameter - in mm, so if set at 10, its 1 cm
%- distance between adjacent electrodes
%gsp=20; %zg edited from 50
%gsp=10;
%gsp=20;
%gsp =4 ;
%gsp = 50;

c=zeros(length(cortex(:,1)),1);
for i=1:length(electrodes(:,1))
    b_z=abs(brain(:,3)-electrodes(i,3));
    b_y=abs(brain(:,2)-electrodes(i,2));
    b_x=abs(brain(:,1)-electrodes(i,1));
    %d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2).^.5)/gsp^.5); %exponential fall off
    d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian
    c=c+d';
end
vert_values = c;
%keyboard;
% c=(c/max(c));
c_h=tripatch(cortex, 'nofigure', c');

shading interp;
a=get(gca);

%%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
%%RELEVANCE OF CHANNEL FOR COMPARISONs ACROSS CORTICES
d=a.CLim;
set(gca,'CLim',[-max(abs(d)) max(abs(d))])
colormap(cm)
lighting phong; %play with lighting...
%material shiny;
material dull;
% material([.3 .8 .1 10 1]);
% material([.2 .9 .2 50 1]); %  BF: editing mesh viewing attributes
axis off

if do_lighting
    l=light;
    
    if strcmp(hemi,'lh')
        view(270, 0);
        set(l,'Position',[-1 0 0],'Color',[0.8 0.8 0.8]);
    elseif strcmp(hemi,'rh')
        view(90, 0);
        set(l,'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
    elseif strcmp(hemi, 'top')
        view(0, 90);
        set(l,'Position',[0 0 1],'Color',[0.8 0.8 0.8]);
    elseif strcmp(hemi, 'bottom')
        view(0, -90);
        set(l,'Position',[0 0 -1],'Color',[0.8 0.8 0.8]);
    end
end
