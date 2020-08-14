% Figure 7 and Figure S5
% 
% Hamilton, Oganian, and Chang
% 
% Stimulation effects and single subject reconstructions
%% set path
clc;
clear all;
addpath(genpath('../util1'));
heschl_load_data;
config_paths;

%% load anatomy data
SID = {'S06', 'S09', 'S08',  'S10', 'S11', 'S12','S13'};
allHemi={ 'lh', 'lh', 'bil', 'lh', 'lh', 'lh', 'lh'};%'lh',
nSID = length(SID);
heschl_load_anatomy;
%% load list of channels
xlsfile = sprintf('%s/stim/HG_stim_summary.xlsx', paper_data_dir);
chanlist = readtable(xlsfile, 'Sheet','quantitative');

chanlist.effect = sign(chanlist.passive_effect+chanlist.repetition_effect);
chanlist.rep_simple_eff=chanlist.repetition_effect==1;
chanlist.other_effect = chanlist.passive_effect ==2;
chanlist.distraction_effect = chanlist.repetition_effect ==3;
%% add anatomy to channel list

for cs = 1:length(SID)
    celidx = strcmpi(chanlist.patient, SID{cs}) & isfinite(chanlist.tdt_electrode1);
    cel = chanlist.tdt_electrode1(celidx);
    chanlist.customAna1(celidx)=imgNative.(SID{cs}).customAna(cel) ;
    celidx = strcmpi(chanlist.patient, SID{cs}) & isfinite(chanlist.tdt_electrode2);
    cel = chanlist.tdt_electrode2(celidx);
    chanlist.customAna2(celidx)=imgNative.(SID{cs}).customAna(cel) ;
    
end
chanlist.customAna = max(chanlist.customAna1, chanlist.customAna2);
%% load list of unique channels

chanlistU = readtable(xlsfile, 'Sheet','unique');
chanlistU.effect = sign(chanlistU.passive_effect+chanlistU.repetition_effect);
chanlistU.rep_simple_eff=chanlistU.repetition_effect==1;
chanlistU.other_effect = chanlistU.passive_effect ==2;
chanlistU.distraction_effect = chanlistU.repetition_effect ==3;
chanlistU.aura_effect = contains(chanlistU.passive_description, 'aura');
chanlistU.attenuation_effect = chanlistU.passive_effect ==3;

for cs = 1:length(SID)
    celidx = strcmpi(chanlistU.patient, SID{cs}) & isfinite(chanlistU.tdt_electrode1);
    cel = chanlistU.tdt_electrode1(celidx);
    chanlistU.customAna1(celidx)=imgNative.(SID{cs}).customAna(cel) ;
    celidx = strcmpi(chanlistU.patient, SID{cs}) & isfinite(chanlistU.tdt_electrode2);
    cel = chanlistU.tdt_electrode2(celidx);
    chanlistU.customAna2(celidx)=imgNative.(SID{cs}).customAna(cel) ;
    
end
chanlistU.customAna = max(chanlistU.customAna1, chanlistU.customAna2);
chanlistU.repetition_effect(isnan(chanlistU.repetition_effect))=0;
chanlistU.passive_effect(isnan(chanlistU.passive_effect))=0;
chanlistU.all_eff = sign(chanlistU.passive_effect) + 2*sign(chanlistU.repetition_effect);
chanlistU.all_eff(chanlistU.distraction_effect == 1) = 4;
chanlistU.all_eff(chanlistU.attenuation_effect == 1) = 5;
chanlistU.all_eff(chanlistU.aura_effect == 1) = 6;

% final electrode definition

chanlistU.simple_eff = chanlistU.all_eff;
chanlistU.simple_eff(chanlistU.simple_eff==6) = 0; % exclude aura channels
chanlistU.simple_eff(chanlistU.simple_eff==4) = 1; % add passive+distraction into the passive pool
% simple anatomy
chanlistU.simpleAna = chanlistU.customAna;
chanlistU.simpleAna(ismember(chanlistU.simpleAna, 1:3)) = 1;
chanlistU.simpleAna(ismember(chanlistU.simpleAna, 4:5)) = 2;
effNames = {'only passive percept','only repetition impaired', 'hearing affected, rep delay', ...
    'passive + distraction in rep', 'attenuated background sounds', 'aura-like'};

%% ---- plots
%% on MNI brain
% Figure 7A
hemifls = {'lh', 'rh'};
csiz = [30 60 60 60];
figure;
ca = 1; i = 1;
%subplot(2,2,ca*2 +i);
cla; hold on;
c_h = ctmr_gauss_plot(mniMesh.temporal.(hemifls{i}),[0 0 0],0,hemifls{i});
c_h.FaceAlpha = .5; loc_view((-1)^i*90,60);
inclels = isfinite(chanlistU.tdt_electrode1) & contains(chanlistU.hemisphere, hemifls{i}(1)) & ...
    sign(chanlistU.customAna)==ca;
condi=1; legends = cell(0);clegs=cell(0); effi=1;
effCols = 'wrkmgb';anaM ='odp<>';
for effu = [1 2 5]
    allels=[];cloc=[];allAna =[];
    for csid = 1:7
        cs = SID{csid};
        celmat = imgmni.(cs).elecmatrix;
        
        cel=[];ilocs=[];
        % both tested and no effect
        celidx = strcmpi(chanlistU.patient, cs) & inclels & chanlistU.simple_eff==effu;
        cel(:,1) = chanlistU.tdt_electrode1(celidx); cel(:,2) = chanlistU.tdt_electrode2(celidx);
        elu = unique(cel(:));
        [~, celmat(elu,:)] = nearest_electrode_vert(mniMesh.temporal.(hemifls{i}),celmat(elu,:));
        
        ti=chanlistU.customAna1(celidx,1)>0 & chanlistU.customAna2(celidx,1)>0;
        %         ilocs(ti,:) = (celmat(cel(ti,1),:) + celmat(cel(ti,2),:))/2;
        %         ti=chanlistU.customAna1(celidx,1)>0 & chanlistU.customAna2(celidx,1)==0;
        %         ilocs(ti,:) = celmat(cel(ti,1),:);
        %         ti=chanlistU.customAna1(celidx,1)==0 & chanlistU.customAna2(celidx,1)>0;
        %         ilocs(ti,:) = celmat(cel(ti,2),:);
        ilocs = (celmat(cel(:,1),:) + celmat(cel(:,2),:))/2;
        cloc = [cloc;ilocs];
        allels=[allels; cel];
        allAna =[allAna; chanlistU.customAna(celidx)];
    end
    if ~isempty(cloc)
        for j = 1:5
            % pl(effi) = scatter3(cloc(allAna==j,1),cloc(allAna==j,2), cloc(allAna==j,3),60,'filled', effCols(effu), 'Marker', anaM(j), 'MarkerEdgeColor', 'k');
            pl(effi) = scatter3(cloc(allAna==j,1),cloc(allAna==j,2), cloc(allAna==j,3),60,'filled', effCols(effu), 'Marker', 'o', 'MarkerEdgeColor', 'k');
        end
        %effpl(effi) = scatter3(0,0,0,40,'filled', effCols(effu), 'Marker', 'o','MarkerEdgeColor', 'k');
        %text(cloc(:,1),cloc(:,2), cloc(:,3), num2str(allels(:,1)));
        clegs{effi} = effNames{effu};
        effi =effi+1;
    end
end
legend(pl, clegs)
% for cm = 1:5
%     markerpl(cm) = scatter3(0,0,0,40,'filled', 'k', 'Marker', anaM(cm),'MarkerEdgeColor', 'k');
% end
% for effu = 1:6
%
% end
%
% cl=legend([markerpl, effpl], [includeAreasName, clegs]);
% cl.Location = 'north';
clear pl effpl clegs;


%% MNI/native single subjects
% Figure S5

hemifls = {'lh', 'rh'};
csiz = [30 60 60 60];
figure;
ca = 1; i = 1;
%subplot(2,2,ca*2 +i);


effCols = 'wrkmgb';
for csid = 1:7
    subplot(2,4,csid)
    cs =SID{csid};
    
    cla; hold on;
    cimg = imgNative.(cs);
    if ~strcmpi(cs, 'S08')
        c_h = ctmr_gauss_plot(cimg.temporal,[0 0 0],0,'lh');
        inclels = isfinite(chanlistU.tdt_electrode1) & contains(chanlistU.hemisphere, hemifls{1}(1)) & ...
            sign(chanlistU.customAna)==ca;
        
    else
        c_h = ctmr_gauss_plot(cimg.temporal_rh,[0 0 0],0,'rh');
        %c_h = ctmr_gauss_plot(cimg.cortex_rh,[0 0 0],0,'rh');
        
        inclels = isfinite(chanlistU.tdt_electrode1) & contains(chanlistU.hemisphere, hemifls{2}(1)) & ...
            sign(chanlistU.customAna)==ca;
        
    end
    %     cimg = imgmni.(cs);
    %     c_h = ctmr_gauss_plot(mniMesh.temporal.(hemifls{i}),[0 0 0],0,hemifls{i});
    c_h.FaceAlpha = .5; loc_view((-1)^i*90,60);
    for effu = [1 2 5]
        cs = SID{csid};
        cel=[];ilocs=[];
        % both tested and no effect
        celidx = strcmpi(chanlistU.patient, cs) & inclels & chanlistU.all_eff==effu;
        cel(:,1) = chanlistU.tdt_electrode1(celidx); cel(:,2) = chanlistU.tdt_electrode2(celidx);
        cloc = (cimg.elecmatrix(cel(:,1),:) + cimg.elecmatrix(cel(:,2),:))/2;
        
        if ~isempty(cloc)
            pl(effu) = scatter3(cloc(:,1),cloc(:,2), cloc(:,3),60,'filled', effCols(effu), 'Marker', 'o', 'MarkerEdgeColor', 'k');
            %text(cloc(:,1),cloc(:,2), cloc(:,3), num2str(cel(:,1)));
        end
    end
    title(SID{csid});
end
subplot(2,4,8)
cla; hold on;% legend
for effu = [1 2 5]
    effpl(effu) = scatter3(0,0,0,1,'filled',  effCols(effu), 'Marker', 'o','MarkerEdgeColor', 'k');
end
axis off;
cl=legend(effNames);

clear pl;

lgpaper = {'hallucination, no effect on repetition', 'no hallucination, repetition interrupted', 'attenuated background sounds'};
legend(lgpaper);

%% ---- counts
%% number of different electrodes stimulated
uSites = (unique(chanlist(:,[1 3 4,16]),'rows'));
crosstab(uSites.customAna, uSites.patient)
sum(uSites.customAna>0)
crosstab(uSites.customAna)
%% number of channels with effects
crosstab(chanlist.effect, chanlist.customAna)
%% channels with evoked sound percept

cel = find(chanlist.passive_effect==1 & chanlist.customAna>0);
height(unique(chanlist(cel,[1 3 4,16]),'rows'))
crosstab(chanlist.customAna(cel),sign(chanlist.repetition_effect(cel)))

crosstab(chanlist.patient(cel))

%% channels with repetition effect

cel = find(chanlist.repetition_effect==1 & chanlist.customAna>0);
unique(chanlist.customAna(cel))
crosstab(chanlist.customAna(cel),sign(chanlist.repetition_effect(cel)))

crosstab(chanlist.patient(cel))


%% barplot overall summary of counts - full anatomy
inclel = find(ismember(chanlistU.all_eff, [1 2 5]) & chanlistU.customAna>0);

figure
a=crosstab(chanlistU.all_eff(inclel),chanlistU.customAna(inclel));
areau = unique(chanlistU.customAna(inclel));
cb=bar(a);
xticklabels(effNames([1 2 5]));
set(gca, 'XTickLabelRotation', 60)
ylabel('# electrodes');
for i =1:3
    cb(i).FaceColor = areaCols(areau(i),:);
    cb(i).FaceAlpha = .5;
end
vertline((1.5:1:3.5)');
legend(includeAreasName(areau));

%% barplot overall summary of counts - simple anatomy: Figure 7A inset
pleff = [1 2 5];
inclel = find(ismember(chanlistU.all_eff, pleff) & chanlistU.customAna>0);
figure;
cla;hold on;
a=crosstab(chanlistU.all_eff(inclel),chanlistU.simpleAna(inclel));
cb=bar(a');
xticks([1, 2]);xticklabels({'temporal plane', 'STG'});

ylabel('# electrodes');
ylim([0 17]);
vertline((1.5)');
for i = 1:3
    cb(i).FaceColor = effCols(pleff(i));
    cb(i).FaceAlpha = .5;
    cb(i).LineWidth=2;
end

legend(effNames(pleff));

%% ------ plot sentence onset ERPs

%% load out structures
outdir = sprintf('%s/TIMIT',paper_data_dir);
for cs =1:length(SID)
    try
        coutfile = dir(fullfile(outdir, [SID{cs} '*.mat']));
        outall.(SID{cs}) = getfield(load(fullfile(coutfile(1).folder, coutfile(1).name)), 'out');
    end
end
outSID = fieldnames(outall);
%% meanresp for all electrodes

for cs = 1:length(outSID)
    ntp(cs) = min(arrayfun(@(x) size(x.resp,2), outall.(outSID{cs})));
    respall.(outSID{cs}).all = arrayfun(@(x) mean(x.resp(:,1:ntp(cs),:),3), outall.(outSID{cs}), 'UniformOutput', false);
end

for cs = 1:length(outSID)
    respall.(outSID{cs}).mean = nanmean(reshape(cell2mat(respall.(outSID{cs}).all), size(respall.(outSID{cs}).all{1},1), ntp(cs), []),3);
    respall.(outSID{cs}).sem = nansem(reshape(cell2mat(respall.(outSID{cs}).all), size(respall.(outSID{cs}).all{1},1), ntp(cs), []),3);
end

%% get meanresponses for all channels
chanlistResp = nan(height(chanlistU), min(ntp), 2);

for cchan = 1:height(chanlistU)
    cel=[];
    cs = chanlistU.patient{cchan};
    cel(1) = chanlistU.tdt_electrode1(cchan);
    cel(2) = chanlistU.tdt_electrode2(cchan);
    try
        chanlistResp(cchan,:,:) = permute(respall.(cs).mean(cel,1:min(ntp)), [3 , 2, 1]);
        chanlistSem(cchan,:,:) = permute(respall.(cs).sem(cel,1:min(ntp)), [3 , 2, 1]);
    catch
        fprintf(2,'error for %s,%d.. \n', cs, cchan);
    end
end
%% plot erps - Figure 7B

figure
condcols = 'kr';
for i = 1:2
    anaels  = isfinite(chanlistU.tdt_electrode1) & sign(chanlistU.customAna)==1;
    inclels = anaels & chanlistU.all_eff == i;
    cresp=chanlistResp(inclels,:,:);
    cresp(isnan(cresp(:,1,1)),:,:)=[];
    cresp = cresp - mean(cresp(:,1:50,:),2);
    size(cresp)

    hold on;
    plot((1:min(ntp))/100-.5,squeeze(mean(cresp,3))', condcols(i), 'linewidth',2);

end
vertline(0);
cax=gca;
legend([cax.Children(end), cax.Children(2)], effNames(1:2))
xlabel('time to sentence onset (s)'); xlim([-.25 1.25]);
xticks([0 .5 1]);
ylabel('HGA (z-score)'); yticks(0:1:3);
