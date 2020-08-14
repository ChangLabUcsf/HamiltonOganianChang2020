% Figure S1
% 
% Hamilton, Oganian, and Chang
%
% Electrodes from all participants on native brain

close all;
subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10','S11','S12','S13'};
SID = subjects;
allHemi={'lh', 'lh', 'lh',  'lh', 'lh', 'lh', 'lh','lh', 'lh','lh','lh','lh','lh'};

heschl_load_anatomy;
% Potential labels
planumpolare = {'planumpolare','ctx_lh_G_temp_sup-Plan_polar', 'ctx_rh_G_temp_sup-Plan_polar',...
    'ctx_lh_S_temp_sup-Plan_polar', 'ctx_rh_S_temp_sup-Plan_polar'};
planumtemporale = {'planumtemporale','ctx_lh_G_temp_sup-Plan_tempo', 'ctx_rh_G_temp_sup-Plan_tempo',...
    'ctx_lh_S_temp_sup-Plan_tempo', 'ctx_rh_S_temp_sup-Plan_tempo'};
transversetemporal = {'transversetemporal','ctx_lh_G_temporal_transverse','ctx_rh_G_temporal_transverse',...
    'ctx_lh_S_temporal_transverse','ctx_rh_S_temporal_transverse', 'ctx_rh_G_temp_sup-G_T_transv', 'ctx_lh_G_temp_sup-G_T_transv'};
superiortemporal = {'superiortemporal','ctx_rh_G_temp_sup-Lateral','ctx_lh_G_temp_sup-Lateral'};

% Colors for each area
pt_color = [0.17, 0.22, 0.58];
pp_color = [0.62, 0.22, 0.58];
mstg_color = [0.69, 0.70, 0.21];
pstg_color = [0.69, 0.12, 0.14];
hg_color = [0.06, 0.69, 0.29];

% Loop through the subjects, plot each temporal lobe and electrodes
for i=1:length(subjects)
    subj = subjects{i};
    figure(i)
    if ~strcmp(subj, 'S08')
        temporal = imgNative.(subj).temporal;
    else
        temporal = imgNative.(subj).temporal_lh;    
    end
    ctmr_gauss_plot(temporal, [0,0,0], 0, 'lh', 1);
    alpha(0.5);
    elecmatrix = imgNative.(subj).elecmatrix;
    anatomy = imgNative.(subj).newAnatomy;
    %anatomy = imgmni.(subj).anatomy;
    %load('/Applications/freesurfer_v5/subjects/EC75/elecs/TDT_elecs_all.mat');
    
    rh=find(elecmatrix(:,1) > 0);
    elecmatrix(rh,:) = [];
    anatomy(rh,:) = [];
    pt_elecs = find(ismember(anatomy(:,4), planumtemporale));
    hg_elecs = find(ismember(anatomy(:,4), transversetemporal));
    pp_elecs = find(ismember(anatomy(:,4), planumpolare));
    STG_elecs = find(ismember(anatomy(:,4), superiortemporal));
    
    pSTG_elecs = find(ismember(anatomy(:,4), 'pSTG'));
    mSTG_elecs = find(ismember(anatomy(:,4), 'mSTG'));
    
    el_add(elecmatrix(pt_elecs,:), 'color', pt_color, 'msize', 20);
    el_add(elecmatrix(hg_elecs,:), 'color', hg_color, 'msize', 20);
    el_add(elecmatrix(pp_elecs,:), 'color', pp_color, 'msize', 20);
    el_add(elecmatrix(pSTG_elecs,:), 'color', pstg_color, 'msize', 20);%, 'numbers', pSTG_elecs);
    el_add(elecmatrix(mSTG_elecs,:), 'color', mstg_color, 'msize', 20);%, 'numbers', mSTG_elecs);
    title(subj);
    set(gcf,'color','w');
    
    loc_view(-118, 40);
    %keyboard;
    
end