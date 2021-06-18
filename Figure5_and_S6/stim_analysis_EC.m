%% set path
clc;
clear all;
addpath(genpath('../util1'));

%% load anatomy data
SID = {'EC143', 'EC193','EC180', 'EC118', 'EC151', 'EC219','EC220'};%
patienNum = {'S06', 'S07', 'S08',  'S10', 'S11', 'S12','S13'};
allHemi={ 'lh', 'lh', 'bil', 'lh', 'lh', 'lh', 'lh'};%'lh',
nSID = length(SID);
heschl_load_anatomy;
%% load list of channels
xlsfile = '../../data/stim/HG_stim_summary.xlsx';
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
%% plot electrodes with effect, without a customAna label
figure
i=0;
for sid = [6 2 4 7]
    i=i+1;
    subplot(2,2,i);
    cla;hold on;
    
    cs = SID{sid};
    celidx = strcmpi(chanlist.patient, cs) & isfinite(chanlist.tdt_electrode1) & chanlist.customAna==0 & ...
        (chanlist.repetition_effect==1 | chanlist.passive_effect==1);
    plel = chanlist.tdt_electrode1(celidx);
    if sid~=2
        c_h=ctmr_gauss_plot(imgNative.(cs).temporal ,[0 0 0],0,allHemi{sid});
    else
        c_h=ctmr_gauss_plot(imgNative.(cs).temporal_lh ,[0 0 0],0,'lh');
    end
    c_h.FaceAlpha = .5; loc_view(-90,90);
    pl(1) = scatter3(imgNative.(cs).elecmatrix(plel,1),...
        imgNative.(cs).elecmatrix(plel,2),imgNative.(cs).elecmatrix(plel,3), 'b', 'filled');
    title(cs)
end


a = chanlist(chanlist.customAna==0 & (chanlist.repetition_effect ==1 |chanlist.passive_effect==1),:);


%% ---- plots
%% on MNI brain
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


%% native space for single subject
cs = 'EC220';
chemi = 'lh';
figure; cla; hold on;
c_h = ctmr_gauss_plot(imgNative.(cs).temporal,[0 0 0],0,chemi);
c_h.FaceAlpha = .5; loc_view(-90,90);

inclels = isfinite(chanlist.tdt_electrode1) & contains(chanlist.hemisphere, chemi(1)) & strcmpi(chanlist.patient, cs);
i=1; legends = cell(0);
for ri = 0:2
    for pi = 0:2
        cel = chanlist.tdt_electrode1(inclels & ...
            chanlist.repetition_effect== ri & chanlist.passive_effect==pi);
        if sum(cel)>0
            pl(i) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
                imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),60, 'filled');
            text(imgNative.(cs).elecmatrix(cel,1),imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3), num2str(cel));
            legends{i} = sprintf('passive = %d, rep = %d', pi, ri);
            i=i+1;
            
        end
    end
end
title(cs);

%legend(pl, legends);
clear pl;
%% MNI/native single subjects

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
    if ~strcmpi(cs, 'EC180')
        c_h = ctmr_gauss_plot(cimg.temporal,[0 0 0],0,'lh');
        inclels = isfinite(chanlistU.tdt_electrode1) & contains(chanlistU.hemisphere, hemifls{1}(1)) & ...
            sign(chanlistU.customAna)==ca;
        
    else
        c_h = ctmr_gauss_plot(cimg.temporal_rh,[0 0 0],0,'rh');
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
    title(patienNum{csid});
end
subplot(2,4,8)
cla; hold on;% legend
for effu = [1 2 5]
    effpl(effu) = scatter3(0,0,0,1,'filled',  effCols(effu), 'Marker', 'o','MarkerEdgeColor', 'k');
end
axis off;
cl=legend(effNames);

if ca == 1 
    sgtitle('Stimulation in auditory cortex ROIs');
else
    sgtitle('Channel labels outside ROIs');
end
clear pl;

lgpaper = {'hallucination, no effect on repetition', 'no hallucination, repetition interrupted', 'attenuated background sounds'};
legend(lgpaper);
sgtitle('');

%% coverage on individual brains
hemifls = {'lh', 'rh'};
csiz = [30 60 60 60];
figure;
ca = 1; i = 1;
for csid = 1:7
    subplot(2,4,csid);
    cla; hold on;
    cs =SID{csid};
    cimg = imgNative.(cs);
    inclels = imgNative.(cs).customAna>0;
    if ~strcmpi(cs, 'EC180')
        c_h = ctmr_gauss_plot(cimg.temporal,[0 0 0],0,'lh');
        c_h(1).FaceAlpha = .5; 
    else
        c_h(1) = ctmr_gauss_plot(cimg.temporal_rh,[0 0 0],0,'rh');
        c_h(2) = ctmr_gauss_plot(cimg.temporal_lh,[0 0 0],0,'lh');
        c_h(1).FaceAlpha = .5; 
        c_h(2).FaceAlpha = .5; 
    end
    
    loc_view((-1)^i*90,60);
    for ai=1:5
        cel = find(imgNative.(cs).newCustomAna==ai);
        if ~isempty(cel)
        pl(ai) = scatter3(cimg.elecmatrix(cel,1),cimg.elecmatrix(cel,2),cimg.elecmatrix(cel,3),...
        40,'filled','MarkerFaceColor',areaCols(ai,:));
        end
    end
    title(SID{csid});
end

%% ---- counts
%% number of different electrodes stimulated
uSites = (unique(chanlist(:,[1 3 4,12]),'rows'));
crosstab(uSites.customAna, uSites.patient)
sum(uSites.customAna>0)
crosstab(uSites.customAna)
%% number of channels with effects
crosstab(chanlist.effect, chanlist.customAna)
%% channels with evoked sound percept

cel = find(chanlist.passive_effect==1 & chanlist.customAna>0);
height(unique(chanlist(cel,[1 3 4,12]),'rows'))
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

%% barplot overall summary of counts - simple anatomy
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
outdir = '/Users/yuliao/Dropbox (UCSF Department of Neurological Surgery)/proj_Heschls/Heschls_STRFs/out_structures';
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
%% plot erps

figure
condcols = 'kr';
for i = 1:2
    anaels  = isfinite(chanlistU.tdt_electrode1) & sign(chanlistU.customAna)==1;
    inclels = anaels & chanlistU.all_eff == i;
    cresp=chanlistResp(inclels,:,:);
    cresp(isnan(cresp(:,1,1)),:,:)=[];
    cresp = cresp - mean(cresp(:,1:50,:),2);
    size(cresp)
    
    %     for j = 1:2
    %         subplot(2,1,1); hold on;
    %         plot((1:min(ntp))/100-.5,squeeze(cresp(:,:,j))', condcols(i));
    %     end
    %subplot(2,1,1); 
    hold on;
    plot((1:min(ntp))/100-.5,squeeze(mean(cresp,3))', condcols(i));
    
    %subplot(2,1,2); hold on
    %cresp = mean(cresp,3);
    %cpl(i) = shadedErrorBar((1:min(ntp))/100-.5,nanmean(cresp,1),nansem(cresp,1), {'color', condcols(i), 'linewidth', 3});
end
vertline(0);
cax=gca;
%legend([cpl.mainLine], effNames(1:2))
legend([cax.Children(end), cax.Children(2)], effNames(1:2))
xlabel('time to sentence onset (s)'); xlim([-.25 1.25]);
xticks([0 .5 1]);
ylabel('HGA (z-score)'); yticks(0:1:3);

%% single electrodes

anaels  = isfinite(chanlistU.tdt_electrode1) & sign(chanlistU.customAna)==1 & chanlistU.all_eff<3;
for effi = 1:2
    figure
    inclels = find(anaels & chanlistU.all_eff == effi);
    cresp=chanlistResp(inclels,:,:);
    %cresp(isnan(cresp(:,1,1)),:,:)=[];
    cresp = cresp - mean(cresp(:,1:50,:),2);
    for eli = 1:size(cresp,1)
        subplot(2,ceil(size(cresp,1)/2),eli);
        cla;hold on;
        for tdti = 1:2
            plot((1:min(ntp))/100-.5,squeeze(cresp(eli,:,tdti))', condcols(effi));
        end
        title(chanlistU.patient(inclels(eli)));
        legend({num2str(chanlistU.tdt_electrode1(inclels(eli))),num2str(chanlistU.tdt_electrode2(inclels(eli)))});
    end
    sgtitle(effNames{effi});
end

%% ----------- Initial labeling etc
%% --- get tdt labels for all electrodes
%% plot lateral grid with tdt and clinical labels
cs='EC219';
figure; cla; hold on;
c_h = ctmr_gauss_plot(imgNative.(cs).cortex ,[0 0 0],0,'lh');
plel = reshape([1:32:225]+[0:2:14]',1,[]);
tdtel = 1:256;
pl(1) = scatter3(imgNative.(cs).elecmatrix(plel,1),...
    imgNative.(cs).elecmatrix(plel,2),imgNative.(cs).elecmatrix(plel,3), 'b', 'filled');
text(imgNative.(cs).elecmatrix(plel,1),...
    imgNative.(cs).elecmatrix(plel,2),imgNative.(cs).elecmatrix(plel,3), num2str((1:length(plel))'), 'fontsize', 14);

plel2 = setdiff(tdtel,plel);
pl(2) = scatter3(imgNative.(cs).elecmatrix(plel2,1),...
    imgNative.(cs).elecmatrix(plel2,2),imgNative.(cs).elecmatrix(plel2,3), 'r', 'filled');
text(imgNative.(cs).elecmatrix(plel2,1),...
    imgNative.(cs).elecmatrix(plel2,2),imgNative.(cs).elecmatrix(plel2,3), num2str(tdtel(plel2)'), 'fontsize', 10);
title(cs);
%% labels - EC219

a=find(contains(imgNative.EC219.anatomy(:,1), 'PII'));
imgNative.EC219.anatomy(a,:)

a=find(contains(imgNative.EC219.anatomy(:,1), 'PG29'))

a=find(strcmpi(imgNative.EC219.anatomy(:,1), 'AG4'))
imgNative.EC219.anatomy(a,:)

a=find(contains(imgNative.EC219.anatomy(:,1), 'AII'))
imgNative.EC219.anatomy(a,:)
%% labels - EC180
clc;
a=find(contains(imgNative.EC180.anatomy(:,1), 'LPISP5'));
imgNative.EC180.anatomy(a,:)
a'
%% labels - EC193

a=find(contains(imgNative.EC193.anatomy(:,1), 'PMGC'));
imgNative.EC193.anatomy(a,:)
a'

a=find(contains(imgNative.EC193.anatomy(:,1), 'TGB'));
imgNative.EC193.anatomy(a,:)
a'
%% labels - EC143

a=find(contains(imgNative.EC143.anatomy(:,1), 'SG13'));
imgNative.EC143.anatomy(a,:)
a'
%% labels - EC118
clc
a=find(contains(imgNative.EC118.anatomy(:,1), 'LPISP'));
imgNative.EC118.anatomy(a,:)
a'
%% labels EC151

clc
a=find(contains(imgNative.EC151.anatomy(:,1), 'HG'));
imgNative.EC151.anatomy(a,:)
a'
%% labels EC220
clc; cs ='EC220';
a=find(contains(imgNative.(cs).anatomy(:,1), 'LMTG1'));
imgNative.(cs).anatomy(a,:)
a'
%% --- plot single subject grids

%% plot in EC143 3  anterior channels with aura effects
celidx = find(strcmpi(chanlistU.patient, 'EC143') & ismember(chanlistU.tdt_electrode1, [263 264 269]));
cs = 'EC143';
figure
subplot(2,1,1) ; hold on
cimg = imgNative.(cs);
c_h = ctmr_gauss_plot(cimg.temporal,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
loc_view(-90,90);
for i = 1:3
    cel = [chanlistU.tdt_electrode1(celidx(i)), chanlistU.tdt_electrode2(celidx(i))];
    pl(i) = plot3(cimg.elecmatrix(cel,1),...
        cimg.elecmatrix(cel,2),cimg.elecmatrix(cel,3),'o-', 'MarkerFaceColor', 'w', 'linewidth', 2);
end
cel3 = find(~cellfun('isempty', strfind(cimg.anatomy(:,1),'SG')));
pl(3) = scatter3(cimg.elecmatrix(cel3,1),...
    cimg.elecmatrix(cel3,2),cimg.elecmatrix(cel3,3),30, 'filled', 'k', 'Marker', 'o');
subplot(2,1,2) ; hold on
cimg = imgmni.(cs);
c_h = ctmr_gauss_plot(mniMesh.temporal.lh,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
loc_view(-90,90);
for i = 1:3
    cel = [chanlistU.tdt_electrode1(celidx(i)), chanlistU.tdt_electrode2(celidx(i))];
    pl(i) = plot3(cimg.elecmatrix(cel,1),...
        cimg.elecmatrix(cel,2),cimg.elecmatrix(cel,3),'o-', 'MarkerFaceColor', 'w', 'linewidth', 2);
end
cel3 = find(~cellfun('isempty', strfind(cimg.anatomy(:,1),'SG')));
pl(3) = scatter3(cimg.elecmatrix(cel3,1),...
    cimg.elecmatrix(cel3,2),cimg.elecmatrix(cel3,3),30, 'filled', 'k', 'Marker', 'o');

%% EC143
cs = 'EC143';
figure
hold on
cimg = imgNative.(cs);
%c_h = ctmr_gauss_plot(cimg.temporal,[0 0 0],0,'lh');
cimg = imgmni.(cs);
c_h = ctmr_gauss_plot(mniMesh.temporal.lh,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
loc_view(-90,90);
cel3 = find(~cellfun('isempty', strfind(cimg.anatomy(:,1),'SG')));
pl(3) = scatter3(cimg.elecmatrix(cel3,1),...
    cimg.elecmatrix(cel3,2),cimg.elecmatrix(cel3,3),30, 'filled', 'k', 'Marker', 'o');
% text(imgNative.(cs).elecmatrix(cel3,1),...
%     imgNative.(cs).elecmatrix(cel3,2),imgNative.(cs).elecmatrix(cel3,3), num2str(cel3-cel3(1)+1));

cel1 = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.passive_effect==1);
pl(1) = scatter3(cimg.elecmatrix(cel1,1),...
    cimg.elecmatrix(cel1,2),cimg.elecmatrix(cel1,3),60, 'filled', 'r', 'Marker', 'o');
text(cimg.elecmatrix(cel1,1),...
    cimg.elecmatrix(cel1,2),cimg.elecmatrix(cel1,3), num2str(cel1), 'color', 'r');
cel2 = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.passive_effect==0);
pl(2) = scatter3(cimg.elecmatrix(cel2,1),...
    cimg.elecmatrix(cel2,2),cimg.elecmatrix(cel2,3),60, 'filled', 'b', 'Marker', 'o');
text(cimg.elecmatrix(cel2,1),...
    cimg.elecmatrix(cel2,2),cimg.elecmatrix(cel2,3), num2str(cel2), 'color', 'b');


legend(pl, {'sound', 'no sound', 'all channels'});
title('EC143');
%% EC193
cs = 'EC220';
figure
hold on
c_h = ctmr_gauss_plot(imgNative.(cs).temporal,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
loc_view(-90,90);
cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.passive_effect==1);
pl(1) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),60, 'filled', 'r', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.passive_effect==0);
pl(2) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),60, 'filled', 'b', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.repetition_effect==1);
pl(3) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),20,'filled', 'g', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.repetition_effect==0);
pl(4) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),20,'filled', 'y', 'Marker', 'o');
legend(pl, {'sound', 'no sound', 'no rep', 'rep intact'});
title(cs);
%% EC 180

cs = 'EC180';
figure
cla; hold on;
c_h = ctmr_gauss_plot(imgNative.(cs).temporal_lh,[0 0 0],0,'lh');
c_h.FaceAlpha = .5; loc_view(-90,90);
inclels = isfinite(chanlist.tdt_electrode1) & strcmpi(chanlist.hemisphere, 'left');
cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & ...
    chanlist.passive_effect==1 & inclels );
pl(1) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),60, 'filled', 'r', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.passive_effect==0 & inclels);
pl(2) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),60, 'filled', 'b', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.repetition_effect==1 & inclels);
pl(3) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),20,'filled', 'g', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.repetition_effect==0 &inclels);
pl(4) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),20,'filled', 'y', 'Marker', 'o');

% right
figure; cla; hold on;
c_h = ctmr_gauss_plot(imgNative.(cs).temporal_rh,[0 0 0],0,'rh');
c_h.FaceAlpha = .5; loc_view(-90,90);
inclels = isfinite(chanlist.tdt_electrode1) & strcmpi(chanlist.hemisphere, 'right');
cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & ...
    chanlist.passive_effect==1 & inclels );
pl(1) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),60, 'filled', 'r', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.passive_effect==0 & inclels);
pl(2) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),60, 'filled', 'b', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.repetition_effect==1 & inclels);
pl(3) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),20,'filled', 'g', 'Marker', 'o');

cel = chanlist.tdt_electrode1(strcmpi(chanlist.patient, cs) & chanlist.repetition_effect==0 &inclels);
pl(4) = scatter3(imgNative.(cs).elecmatrix(cel,1),...
    imgNative.(cs).elecmatrix(cel,2),imgNative.(cs).elecmatrix(cel,3),20,'filled', 'y', 'Marker', 'o');


legend(pl, {'sound', 'no sound', 'no rep', 'rep intact'});
sgtitle(cs);
