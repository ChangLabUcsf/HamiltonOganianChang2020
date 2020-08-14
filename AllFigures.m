% AllFigures.m
%
% Recreate all figures for Hamilton, Oganian, and Chang (2020). 
% Topography of speech acoustic and phonological feature encoding across the expanse of human auditory cortex
%
% You must download the data from OSF prior to running these scripts.
% See the README.md for more information
%
% Liberty Hamilton, Yulia Oganian, Edward Chang. 2020
% The University of San Francisco
% 

%% load data
clear all; close all;

addpath(genpath('util1'));
config_paths;
fprintf(1,'****Loading data from %s****\n', paper_data_dir);

if ~exist(paper_data_dir,'dir')
    error(sprintf('ERROR: %s does not exist. Have you downloaded the data and set your data directory by editing config_paths?\n', paper_data_dir));
end

heschl_load_data;

%% 
addpath(genpath(pwd));

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
Figure6ABCDE;
Figure6FG;

%%
Figure7_and_S5;

%%
FigureS1
