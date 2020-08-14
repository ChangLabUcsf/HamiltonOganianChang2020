% AllFigures.m
%
% This code generates all of the figures for Hamilton, Oganian, and Chang 2020.
%
% You will first need to download the data from OSF. Save the location of this data
% directory and add it to the `config_paths.m` file so that the other scripts
% know which directory to reference.
%
% run figures

% Download the data

%% load data
addpath(genpath('../util1'));
heschl_load_data;
config_paths;
%% 
addpath(pwd);

%%
Figure1BC;

%%
Figure2A;
Figure2CDE;
Figure2BFGH_S2_S3;

%%
Figure3;

%%
Figure4_and_S4

%%
Figure5;

%%
Figure6_ABCDE;
Figure6FG;

%%
Figure7_and_S5;

%%
Figure2_S2_S3
