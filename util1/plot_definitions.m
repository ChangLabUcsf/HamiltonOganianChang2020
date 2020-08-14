onsetColor =  [85 184, 72]/256;
rampColor = [180, 71, 143]/256;
corpColors = [44 123 182;215 25 28]/255;
pitchColor = [71 180 163]/256;
pointsize = 60;
nColors = 100;
rbcm = [[repmat(linspace(0,1, nColors)',1,2),ones(nColors, 1)] ; [ones(nColors, 1), repmat(linspace(1,0, nColors)',1,2)]];
rcm = [ones(nColors, 1), repmat(linspace(1,0, nColors)',1,2)];
%%
figDir ='figures';
mkdir(figDir);

%% 
addpath(genpath('cbrewer/'));
strfcmap = flipud(cbrewer('div','PRGn',256));