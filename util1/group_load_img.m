function [img, mni_meshes] = group_load_img(SID, elect_space, datapath, hemifl)
if nargin<4
    switch elect_space
        case 'mni'
            hemifl=[];
        case 'native'
            error('need hemisphere flags to load meshes');
    end
end
%% load data
tic
for csid = 1:length(SID)
    cs = SID{csid};
    disp(cs);

    imgPath = fullfile(datapath, cs);
    mnipath = fullfile(datapath, 'cvs_avg35_inMNI152', 'Meshes');

    switch elect_space
        case 'native'
            imgfile = fullfile(imgPath, 'elecs/TDT_elecs_all_anatfix.mat');
            if exist(imgfile, 'file')
                img.(cs) = load(imgfile);
            end
            hemis = hemifl{csid};
            img.(cs).hemifl = hemifl{csid};
            imgfile = fullfile(imgPath, 'Meshes', [cs '_'  hemis '_pial.mat']);
            if exist(imgfile, 'file')
                img.(cs).cortex = getfield(load(imgfile, 'cortex'), 'cortex');
            end
            
        case 'mni'
            try
                img.(cs) = load(fullfile(imgPath, 'elecs/TDT_elecs_all_warped.mat'));
            catch
                warning('no warped electrode file for subject %s.' ,cs);
                img.(cs) = [];
            end
    end
end
if ~exist('img', 'var'), img = [];end
if strcmpi(elect_space, 'mni')
    mni_meshes.lh = getfield(load(fullfile(mnipath,'cvs_avg35_inMNI152_lh_pial.mat') , 'cortex'), 'cortex');
    mni_meshes.rh = getfield(load(fullfile(mnipath, 'cvs_avg35_inMNI152_rh_pial.mat'), 'cortex'), 'cortex');
end

toc
