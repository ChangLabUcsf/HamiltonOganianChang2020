% AllFigures.m
%
% Recreate all figures for Hamilton, Oganian, Hall, and Chang (2021). 
% Topography of speech acoustic and phonological feature encoding across the expanse of human auditory cortex
%
% You must download the data from OSF prior to running these scripts.
% See the README.md for more information
%
% Liberty Hamilton, Yulia Oganian, Jeffery Hall, Edward Chang. 2021
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
Figure2ABC;
Figure2DE;

%%
Figure3A;
Figure3BFGH_S2_S3_7areas_MID;
Figure3CDE_MID;

%%
Figure4

%%
Figure5_and_S6;

%%
Figure6ABCDE;
Figure6FG;

%%
Figure7_and_S5;

%%
FigureS1
