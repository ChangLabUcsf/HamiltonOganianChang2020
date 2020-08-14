% set path
%addpath(genpath('../../../ecog_scripts'))
%[datapath, dpath] = setDatapath;
config_paths;

%% patient list

SID = {'S01', 'S02', 'S03', 'S04','S05', 'S06', 'S07', 'S08', 'S09'};
allHemi={'lh', 'lh', 'lh',  'lh', 'lh', 'lh', 'lh','bil', 'lh'};
nSID = length(SID);

stimSID= { 'S12','S10', 'S11', 'S13'};%

%% plot definitions
plot_definitions;
%% load imaging data
heschl_load_anatomy;
%% load strfs
strfdatapath = sprintf('%s/STRFs',paper_data_dir);
clear strfall;
strfdirs = dir(fullfile(strfdatapath));
strfdirs = {strfdirs([strfdirs.isdir]==1).name};
strfdirs(strcmpi(strfdirs,'puretone'))=[];
strfdirs(contains(strfdirs,'perm'))=[];

for cdir = 3:length(strfdirs)
    strffiles = struct2cell(dir(fullfile(strfdatapath, strfdirs{cdir}, '*.hf5')));
    if ~isempty(strffiles)
        strffiles  = strffiles(1,:);
        %strffilenames = strrep(strffiles , 'allchans', 'allchans_');
        strffilenames = strffiles;
        for cf = 1:length(strffiles)
            cname = textscan(strffiles{cf}, '%s', 'Delimiter', '_');
            csid = cname{1}{1};
            if ismember(csid,SID)
                cpred(1) = strfind(strffilenames{cf}, 'STRF')+5;
                cpred(2) = strfind(strffilenames{cf}, '.hf5')-1;
                chf5 = h5info(fullfile(strfdatapath, strfdirs{cdir}, strffiles{cf}));
                for i = 1:length({chf5.Datasets.Name})
                    strfall((strcmpi(SID, csid)),cdir).(chf5.Datasets(i).Name) = h5read(fullfile(strfdatapath, strfdirs{cdir}, strffiles{cf}),['/' chf5.Datasets(i).Name]);
                end
                strfall((strcmpi(SID, csid)),cdir).sid = csid;
                strfall((strcmpi(SID, csid)),cdir).name = strffilenames{cf}(cpred(1):cpred(2));
                strfall((strcmpi(SID, csid)),cdir).fullname = strffilenames{cf};
                strfall((strcmpi(SID, csid)),cdir).dirname = strfdirs{cdir};
            end
        end
    else
        for i =1:length(SID)
            strfall(i,cdir).name = '';
            strfall(i,cdir).fullname = '';
            strfall(i,cdir).dirname = '';            
        end
    end    
end

strfall = strfall(:,~cellfun('isempty',{strfall(end,:).dirname}'));
strfnames = {strfall(end,:).dirname}';

% reshape weights
for i = 1:size(strfall,1)
    for j = 1:size(strfall,2)
        %         ntp = length(strfall(i,j).delays);
        nel = length(strfall(i,j).vcorrs);
        strfall(i,j).beta = reshape(strfall(i,j).wts(:,2:end), nel, [],60);
    end
end
%% load pure tone strfs

strfdatapath = sprintf('%s/STRFs/puretone', paper_data_dir);

strffiles = struct2cell(dir(fullfile(strfdatapath, '*.hf5')));
strffiles  = strffiles(1,:);
for cf = 1:length(strffiles)
    cname = textscan(strffiles{cf}, '%s', 'Delimiter', '_');
    csid = cname{1}{1};
    strffilenames = strffiles;
    if ismember(csid,SID)
        cpred(1) = strfind(strffilenames{cf}, 'STRF')+5;
        cpred(2) = strfind(strffilenames{cf}, '.hf5')-1;
        chf5 = h5info(fullfile(strfdatapath,  strffiles{cf}));
        for i = 1:length({chf5.Datasets.Name})
            strfTone((strcmpi(SID, csid)),1).(chf5.Datasets(i).Name) = h5read(fullfile(strfdatapath,  strffiles{cf}),['/' chf5.Datasets(i).Name]);
        end
        strfTone((strcmpi(SID, csid)),1).sid = csid;
        strfTone((strcmpi(SID, csid)),1).name = strffilenames{cf}(cpred(1):cpred(2));
        strfTone((strcmpi(SID, csid)),1).fullname = strffilenames{cf};
    end
end

% reshape weights
for i = 1:size(strfTone,1)
    for j = 1:size(strfTone,2)
        %         ntp = length(strfall(i,j).delays);
        nel = length(strfTone(i,j).vcorrs);
        strfTone(i,j).beta = reshape(strfTone(i,j).wts(:,2:end), nel, [],60);
    end
end
%% load strf permutations
strfdatapath = sprintf('%s/STRFs/', paper_data_dir);
strfdirs = dir(fullfile(strfdatapath,'*perm*'));
strfdirs = {strfdirs([strfdirs.isdir]==1).name};
for cdir = 1:length(strfdirs)
    strffiles = struct2cell(dir(fullfile(strfdatapath, strfdirs{cdir}, '*.hf5')));

    if ~isempty(strffiles)
        strffiles  = strffiles(1,:);
        strffilenames = strffiles;
        for cf = 1:length(strffiles)
            cname = textscan(strffiles{cf}, '%s', 'Delimiter', '_');
            csid = cname{1}{1};
            if ismember(csid,SID)
                cpred(1) = strfind(strffilenames{cf}, 'STRF')+5;
                cpred(2) = strfind(strffilenames{cf}, '.hf5')-1;
%             try
%                 chf5 = h5info(fullfile(strfdatapath, strfdirs{cdir}, strffiles{cf}), 'TextEncoding', 'UTF-8');
                strfperm((strcmpi(SID, csid)),cdir).perm_corrs = h5read(fullfile(strfdatapath, strfdirs{cdir}, strffiles{cf}),'/perm_corrs');
%                 for i = 1:length({chf5.Datasets.Name})
%                     strfperm((strcmpi(SID, csid)),cdir).(chf5.Datasets(i).Name) = h5read(fullfile(strfdatapath, strfdirs{cdir}, strffiles{cf}),['/' chf5.Datasets(i).Name]);
%                 end
                strfperm((strcmpi(SID, csid)),cdir).sid = csid;
                strfperm((strcmpi(SID, csid)),cdir).name = strffilenames{cf}(cpred(1):cpred(2));
                strfperm((strcmpi(SID, csid)),cdir).fullname = strffilenames{cf};
                strfperm((strcmpi(SID, csid)),cdir).dirname = strfdirs{cdir};
%             end
            end
        end   
    end    
end

%% load out
out_dir = sprintf('%s/TIMIT/', paper_data_dir);
for si = 1:length(SID)
    outfiles = dir(fullfile(out_dir,...
        [SID{si} '*mat']));
    outall.(SID{si}) = getfield(load(fullfile(outfiles(1).folder, outfiles(1).name)), 'out');
end

%% load feature names
features = load('Timit_features.mat'); 
%% load sentence info
sentdet = load(sprintf('%s/TIMIT/out_sentence_details_timit_all_loudness.mat', paper_data_dir));

%% variance of all models combined

allr.vcorrs=[];allr.vcorrs=[];allr.elecmatrix =[];allr.sid=[];
allr.elnum=[];allr.customAna=[];allr.strf=cell(1,size(strfall, 2)); allr.psi=[];
allr.tone_vcorrs=[];
allr.meansentresp=[];
allr.newCustomAna=[];
meanHGA = cell(length(SID),1);
for i =1:length(SID)
    nel = length(strfall(i,1).vcorrs);
    allr.vcorrs = [allr.vcorrs;[strfall(i,:).vcorrs]];%./max(max([strfall(i,:).vcorrs],[],2),[],1)];
    allr.customAna =[allr.customAna;imgNative.(SID{i}).customAna(1:nel)];
    allr.newCustomAna = [allr.newCustomAna; imgNative.(SID{i}).newCustomAna(1:nel)];
    
    allr.elecmatrix = [allr.elecmatrix; imgmni.(SID{i}).elecmatrix(1:nel,:)];
    
    allr.sid = [allr.sid; i*ones(nel,1)];
    allr.elnum = [allr.elnum; [1:nel]'];
    for cm = 1:size(strfall, 2)        
            allr.strf{cm} = cat(1, allr.strf{cm}, strfall(i,cm).beta);        
    end
%    allr.psi = [allr.psi ,Psiall(i).psimap(:,1:nel)];
    
    if ~isempty(strfTone(i).vcorrs)
        allr.tone_vcorrs = [allr.tone_vcorrs;strfTone(i).vcorrs(1:nel)];
    else
        allr.tone_vcorrs = [allr.tone_vcorrs;nan(size(strfall(i,1).vcorrs))];
    end
    
    %% out
    cout = outall.(SID{i});
    
    for csent =1:length(cout)
        cons = find(sum(cout(csent).phnmatonset,1)>0, 1, 'first');
        cmeanresp = mean(cout(csent).resp(1:nel  ,(cons-49):min(end,(cons + 200)),:),3);
        meanHGA{i}(1:nel,1:size(cmeanresp,2),csent) = cmeanresp;
    end
    allr.meansentresp = [allr.meansentresp; mean(meanHGA{i},3)];
end
allr.strfnames = strfnames;
[allr.maxr,allr.maxMod] = max(allr.vcorrs,[],2);


%% betas 
for cm = 1:length(allr.strf)
        [allr.betas.maxb{cm},allr.betas.maxbLat{cm}] = max(allr.strf{cm},[],3);
        allr.betas.maxbLat{cm}=allr.betas.maxbLat{cm}*10;        
        allr.betas.meanb(:,:,cm) = squeeze(mean(allr.strf{cm}, 2));
end

[allr.betas.maxmeanb,allr.betas.maxmeanbLat] = max(permute(allr.betas.meanb, [1 3 2]),[],3);
allr.betas.maxmeanbLat=allr.betas.maxmeanbLat*10;


%% unique variance

allr.uvNames = {'onset', 'peakRate', 'features', 'abs pitch', 'rel pitch','+spect'};

% onset
allr.uvar(:,1) = allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_peakRate')).^2 - allr.vcorrs(:,strcmpi(allr.strfnames, 'phnfeaturesonset_peakRate')).^2;
allr.uvar(:,2) = allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_peakRate')).^2 - allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset')).^2;
allr.uvar(:,3) = allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_peakRate')).^2 - allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_peakRate')).^2;
allr.uvar(:,4) = allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_abs_f0_peakRate')).^2 - ...
    allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_peakRate')).^2;
allr.uvar(:,5) = allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_abs_f0_peakRate')).^2 - ...
    allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_abs_f0_peakRate')).^2;
allr.uvar(:,6) = allr.vcorrs(:,strcmpi(allr.strfnames, 'spect_zscore')).^2 - ...
    allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_abs_f0_peakRate')).^2;


%% pitch significance
origr = allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_abs_f0_peakRate'));

allr.pitch.relRu =origr.^2 -  allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_abs_f0_peakRate')).^2;
allr.pitch.absRu = origr.^2 - allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_peakRate')).^2;

permrelr = [strfperm(:,strcmpi({strfperm(1,:).name}, 'onset_phnfeaturesonset_permrelative_log_f0_delta_relative_log_f0_abs_f0_peakRate')).perm_corrs]';
nperm = size(permrelr,2);
allr.pitch.relpval = 1-sum(origr>permrelr,2)/nperm;

permabsr = [strfperm(:,strcmpi({strfperm(1,:).name}, 'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_permabs_f0_peakRate')).perm_corrs]';
nperm = size(permrelr,2);
allr.pitch.abspval = 1-sum(origr>permabsr,2)/nperm;

%% feature significance

origr = allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_peakRate'));
allr.feat.peakRu =origr.^2 -  allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset')).^2;
allr.feat.featu = origr.^2 - allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_peakRate')).^2;

permr = [strfperm(:,strcmpi({strfperm(1,:).name}, 'onset_phnfeaturesonset_permpeakRate')).perm_corrs]';
nperm = size(permr,2);
allr.feat.peakRpval = 1-sum(origr>permr,2)/nperm;

permr2 = [strfperm(:,strcmpi({strfperm(1,:).name}, 'onset_permphnfeaturesonset_peakRate')).perm_corrs]';
nperm = size(permr2,2);
allr.feat.featpval = 1-sum(origr>permr2,2)/nperm;


permr3 = [strfperm(:,strcmpi({strfperm(1,:).name}, 'permonset_phnfeaturesonset_peakRate')).perm_corrs]';
nperm = size(permr2,2);
allr.feat.onsetpval = 1-sum(origr>permr2,2)/nperm;


%% ----  thresholds and electrode groups
%% definition of thresholds for electrodes
useAnaName = 'newCustomAna';
minr2  = .01; % used for figures in manuscript draft v7
minsign = .01;
AnalysisEl =(allr.(useAnaName)>0 & allr.maxr.^2>minr2);
AnalysisElnum =find(allr.(useAnaName)>0 & allr.maxr.^2>minr2);
crosstab(allr.(useAnaName)(AnalysisEl))
%% electrode groups
elGroups = cell(5,1);
% onset
elGroups{1} = find(allr.(useAnaName)>0 & allr.elecmatrix(:,1)<-30 & allr.maxr.^2>minr2 &...
    (allr.vcorrs(:,strcmpi(allr.strfnames, 'onset')).^2-allr.vcorrs(:,strcmpi(allr.strfnames, 'spect_zscore')).^2)>0 & ...
    allr.feat.onsetpval<minsign );
groupColor(1,:) = onsetColor;
% feature
elGroups{2} = setdiff(find(allr.(useAnaName)>0 & allr.elecmatrix(:,1)<-30 & allr.maxr.^2>minr2 &  ...
    allr.feat.featpval<minsign & allr.feat.peakRpval<minsign & ...
    (allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_peakRate')).^2-allr.vcorrs(:,strcmpi(allr.strfnames, 'spect_zscore')).^2)>0),elGroups{1});
groupColor(2,:) = rampColor;

% relative pitch
elGroups{3} = allr.(useAnaName)>0 & allr.pitch.relpval<minsign & allr.pitch.relRu>0 & allr.elecmatrix(:,1)<-30 & allr.maxr.^2>minr2;
groupColor(3,:) = corpColors(1,:);
% % absolute pitch
elGroups{4} = allr.(useAnaName)>0 & allr.pitch.abspval<minsign & allr.pitch.absRu>0 & allr.elecmatrix(:,1)<-30 & allr.maxr.^2>minr2;
groupColor(4,:) = corpColors(2,:);

% both pitch types
elGroups{5} = intersect(find(elGroups{4}),find(elGroups{3}));
groupColor(5,:) = [ 0.9290    0.6940    0.1250];

% features separately
elGroups{6} = setdiff(find(allr.(useAnaName)>0 & allr.elecmatrix(:,1)<-30 & allr.maxr.^2>minr2 &  ...
    allr.feat.featpval<minsign & ...
    (allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_peakRate')).^2-allr.vcorrs(:,strcmpi(allr.strfnames, 'spect_zscore')).^2)>0),elGroups{1});


% peakRate separately
elGroups{7} = setdiff(find(allr.(useAnaName)>0 & allr.elecmatrix(:,1)<-30 & allr.maxr.^2>minr2 &  ...
    allr.feat.peakRpval<minsign & ...
    (allr.vcorrs(:,strcmpi(allr.strfnames, 'onset_phnfeaturesonset_peakRate')).^2-allr.vcorrs(:,strcmpi(allr.strfnames, 'spect_zscore')).^2)>0),elGroups{1});




% % spectrotemporal is best model
% cels = find(allr.(useAnaName)>0 & allr.elecmatrix(:,1)<-30 & allr.maxr.^2>minr2 &...
%     (allr.maxr.^2-allr.vcorrs(:,strcmpi(allr.strfnames, 'spect_zscore')).^2)==0);
% cpl(1)=scatter3(allr.elecmatrix(cels,1),...
%     allr.elecmatrix(cels,2),allr.elecmatrix(cels,3),elsize,corpColors(1,:), 'filled');
groupNames = {'Onset', 'features+peakRate', 'rel. pitch', 'abs. pitch', 'rel & abs pitch', 'features only', 'peakRate only'};
%% list of electrodes by feature 
allr.elGroup = zeros(length(allr.elnum),7);
for i = 1:7
    allr.elGroup(elGroups{i},i) = 1;
end
allr.allEls = allr.(useAnaName)>0 & allr.maxr.^2>minr2;
allr.elGroupName = {'onset', 'feat+peakR', 'relP', 'absP', 'realP+absP', 'feat', 'peakR'};

elGroup = allr.elGroup;
elGroupName = allr.elGroupName;
