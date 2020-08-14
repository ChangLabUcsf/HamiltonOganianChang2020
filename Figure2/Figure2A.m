% Figure 2A
%
% Hamilton, Oganian, and Chang
% 

if 0 
    addpath(genpath('../util1'));
    heschl_load_data;
end
config_paths;

%% Show the whole brain (small inset)
subj = 'S03';
load(sprintf('%s/anatomy/S03/Meshes/S03_lh_temporal_pial.mat',paper_data_dir));
load(sprintf('%s/anatomy/S03/Meshes/S03_lh_pial.mat',paper_data_dir));

figure; 
c_h = ctmr_gauss_plot(temporal, [0 0 0], 0, 'lh'); alpha 1;
c_h2 = ctmr_gauss_plot(cortex, [0 0 0], 0, 'lh',0); alpha 0.5;

elecmatrix = imgNative.(subj).elecmatrix;
anatomy = imgNative.(subj).newAnatomy;
%load('/Applications/freesurfer_v5/subjects/EC85/elecs/TDT_elecs_all.mat');

anat_areas = {'pSTG','mSTG','planumtemporale','planumpolare','transversetemporal'};
pt_color = [0.17, 0.22, 0.58];
pp_color = [0.62, 0.22, 0.58];
stg_color = [0.69, 0.70, 0.21];
pstg_color = [0.69, 0.12, 0.14];
hg_color = [0.06, 0.69, 0.29];

h(1) = el_add(elecmatrix(strcmp(anatomy(:,4),'planumtemporale'),:), 'color', pt_color, 'msize', 10);
h(2) = el_add(elecmatrix(strcmp(anatomy(:,4),'planumpolare'),:), 'color', pp_color, 'msize', 10);
h(3) = el_add(elecmatrix(strcmp(anatomy(:,4),'transversetemporal'),:), 'color', hg_color, 'msize', 10);

stg_elecs = union(find(strcmp(anatomy(:,4),'pSTG')), find(strcmp(anatomy(:,4),'mSTG')));

pp=4;
for i=1:length(stg_elecs)
    if (elecmatrix(stg_elecs(i),2) < -9)
        h(pp)=el_add(elecmatrix(stg_elecs(i),:), 'color', pstg_color, 'msize', 10);
        pp=pp+1;
    else
        h(pp)=el_add(elecmatrix(stg_elecs(i),:), 'color', stg_color, 'msize', 10);
        pp=pp+1;
    end
end

% These are the receptive field electrodes
elecs = [265 266 267 268 53:56]; %257 258 259 260 
for e=1:length(elecs)
    if strcmp(anatomy(elecs(e),4), 'planumtemporale')
        clr = pt_color;
    end
    if strcmp(anatomy(elecs(e),4), 'planumpolare')
        clr = pp_color;
    end
    if strcmp(anatomy(elecs(e),4), 'transversetemporal')
        clr = hg_color;
    end
    if strcmp(anatomy(elecs(e),4), 'mSTG')
        clr = stg_color;
    end
    if strcmp(anatomy(elecs(e),4), 'pSTG')
            clr = pstg_color;
        
    end
    
    h(pp) = el_add(elecmatrix(elecs(e), :), 'color', clr, 'edgecol', 'k', 'msize', 10);
    pp=pp+1;
end

loc_view(-113,60);
set(gcf,'color','w');

%% Just show the temporal lobe

figure; 
c_h = ctmr_gauss_plot(temporal, [0 0 0], 0, 'lh'); alpha 1;

elecmatrix = imgNative.(subj).elecmatrix;
%anatomy = imgNative.(subj).anatomy;
%load('/Applications/freesurfer_v5/subjects/EC85/elecs/TDT_elecs_all.mat');

anat_areas = {'pSTG','mSTG','planumtemporale','planumpolare','transversetemporal'};
pt_color = [0.17, 0.22, 0.58];
pp_color = [0.62, 0.22, 0.58];
stg_color = [0.69, 0.70, 0.21];
pstg_color = [0.69, 0.12, 0.14];
hg_color = [0.06, 0.69, 0.29];

h(1) = el_add(elecmatrix(strcmp(anatomy(:,4),'planumtemporale'),:), 'color', pt_color, 'msize', 10);
h(2) = el_add(elecmatrix(strcmp(anatomy(:,4),'planumpolare'),:), 'color', pp_color, 'msize', 10);
h(3) = el_add(elecmatrix(strcmp(anatomy(:,4),'transversetemporal'),:), 'color', hg_color, 'msize', 10);

stg_elecs = union(find(strcmp(anatomy(:,4),'pSTG')), find(strcmp(anatomy(:,4),'mSTG')));

pp=4;
for i=1:length(stg_elecs)
    if (elecmatrix(stg_elecs(i),2) < -9)
        h(pp)=el_add(elecmatrix(stg_elecs(i),:), 'color', pstg_color, 'msize', 10);
        pp=pp+1;
    else
        h(pp)=el_add(elecmatrix(stg_elecs(i),:), 'color', stg_color, 'msize', 10);
        pp=pp+1;
    end
end

% These are the receptive field electrodes
elecs = [265 266 267 268 53:56]; %257 258 259 260 
for e=1:length(elecs)
    if strcmp(anatomy(elecs(e),4), 'planumtemporale')
        clr = pt_color;
    end
    if strcmp(anatomy(elecs(e),4), 'planumpolare')
        clr = pp_color;
    end
    if strcmp(anatomy(elecs(e),4), 'transversetemporal')
        clr = hg_color;
    end
    if strcmp(anatomy(elecs(e),4), 'mSTG')
        clr = stg_color;
    end
    if strcmp(anatomy(elecs(e),4), 'pSTG')
            clr = pstg_color;
        
    end
    
    h(pp) = el_add(elecmatrix(elecs(e), :), 'color', clr, 'edgecol', 'k', 'msize', 10);
    pp=pp+1;
end

loc_view(-113,60);
set(gcf,'color','w');



