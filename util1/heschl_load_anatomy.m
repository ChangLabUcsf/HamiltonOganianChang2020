%% load imaging data
config_paths;
imgpath = sprintf('%s/anatomy', paper_data_dir);
[imgmni, mniMesh] = group_load_img(SID, 'mni',imgpath,allHemi);
imgNative = group_load_img(SID, 'native',imgpath,allHemi);
%% load heschls mesh

for i = 1:length(SID)
    cs = SID{i};
    switch allHemi{i}
        case {'lh', 'rh'}
            try
                imgfile = fullfile(imgpath, cs,'Meshes', [cs '_'  allHemi{i} '_temporal_pial.mat']);
                imgNative.(cs).temporal = getfield(load(imgfile, 'temporal'), 'temporal');
            catch
                fprintf(2,'no heschl mesh for subject %s. \n', cs);
            end
    end
end

mniMesh.temporal.lh = getfield(load(fullfile(imgpath,'cvs_avg35_inMNI152','Meshes','cvs_avg35_inMNI152_lh_temporal_pial.mat')), 'temporal');
mniMesh.temporal.rh = getfield(load(fullfile(imgpath,'cvs_avg35_inMNI152','Meshes','cvs_avg35_inMNI152_rh_temporal_pial.mat')), 'temporal');

%% bilateral patient
%cs = 'EC180';
cs = 'S08';
imgfile = fullfile(imgpath, cs, 'Meshes', [cs '_lh_temporal_pial.mat']);
imgNative.(cs).temporal_lh = getfield(load(imgfile, 'temporal'), 'temporal');

imgfile = fullfile(imgpath, cs, 'Meshes', [cs '_rh_temporal_pial.mat']);
imgNative.(cs).temporal_rh = getfield(load(imgfile, 'temporal'), 'temporal');

imgfile = fullfile(imgpath, cs, 'Meshes', [cs '_lh_pial.mat']);
imgNative.(cs).cortex_lh = getfield(load(imgfile, 'cortex'), 'cortex');

imgfile = fullfile(imgpath, cs, 'Meshes', [cs '_rh_pial.mat']);
imgNative.(cs).cortex_rh = getfield(load(imgfile, 'cortex'), 'cortex');

%% manual corrections of anatomy
imgNative.S02.anatomy([14 22 29 30],4) = {'pSTG'};

if ismember('S07',SID)
    els = contains(imgNative.S07.anatomy(:,1), 'MG');
    imgNative.S07.anatomy(els,4) = {'planumpolare'};
end

if ismember('S09',SID)
    els = contains(imgNative.S09.anatomy(:,1), 'PMG ') & contains(imgNative.S09.anatomy(:,4), 'temp_sup-Lateral');
    imgNative.S09.anatomy(els,4) = {'planumpolare'};
end
%% anatomical areas for original anatomy
includeAreas = { {'planumtemporale','Plan_tempo'},...
    {'transverse', 'transv'},...
    {'Plan_polar','planumpolare'},...
    {'temp_sup-Lateral','superiortemporal','pSTG','mSTG'}};
includeAreasName = { 'PT', 'HG', 'PP', 'pSTG','mSTG'};
nAreas = length(includeAreasName);
areaCols = [0 0 179; 0 179 0;179 0 179 ;179 0 0;179 179 0]/256;
for csid = 1:length(SID)
    cs = SID{csid};
    imgNative.(cs).customAna = zeros(length(imgNative.(cs).anatomy),1);
    for carea = 1:length(includeAreas)
        for cname = 1:length(includeAreas{carea})
            cel = strfind(imgNative.(cs).anatomy(:,4), includeAreas{carea}{cname});
            cel = ~cellfun('isempty', cel);
            imgNative.(cs).customAna(cel) = carea;
        end
    end
    %% split stg into pstg and mid-ant stg
    if ~isempty(imgmni.(cs))
        imgNative.(cs).customAna(imgNative.(cs).customAna==4 & imgmni.(cs).elecmatrix(:,2)>=6) = 5;
        crosstab(imgNative.(SID{csid}).customAna)
    end
end

%% load electrode matrix with fixed STG division in native space
fixed_ana_path = sprintf('%s/anatomy/', paper_data_dir);

for csid = 1:length(SID)
    cs = SID{csid};
    try
        imgNative.(cs).newAnatomy = getfield(load(fullfile(fixed_ana_path, cs, 'elecs', 'TDT_elecs_all_anatfix.mat'), 'anatomy'), 'anatomy');
    catch
        imgNative.(cs).newAnatomy = imgNative.(cs).anatomy;
    end
end

imgNative.S02.newAnatomy([14 22 29 30],4) = {'pSTG'};

%% anatomical areas with native space division

includeAreas2 = { {'planumtemporale','Plan_tempo'},...
    {'transverse', 'transv'},...
    {'Plan_polar','planumpolare'},...
    {'pSTG'}, {'mSTG'}};
for csid = 1:length(SID)
    cs = SID{csid};
    imgNative.(cs).newCustomAna = imgNative.(cs).customAna;
    for carea = 1:length(includeAreas2)
        for cname = 1:length(includeAreas2{carea})
            cel = strfind(imgNative.(cs).newAnatomy(:,4), includeAreas2{carea}{cname});
            cel = ~cellfun('isempty', cel);
            imgNative.(cs).newCustomAna(cel) = carea;
        end
    end
end
