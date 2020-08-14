% run figures

% Download the data

%% load data
config_paths;
fprintf(1,'****Loading data from %s****\n', paper_data_dir);

if ~exist(paper_data_dir,'dir')
    error(sprintf('ERROR: %s does not exist. Have you downloaded the data and set your data directory by editing config_paths?\n', paper_data_dir));
end

addpath(genpath('util1'));
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
Figure6_ABCDE;
Figure6FG;

%%
Figure7_and_S5;

%%
Figure2_S2_S3