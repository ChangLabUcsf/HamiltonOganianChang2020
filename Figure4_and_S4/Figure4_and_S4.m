% Figure 4 and Figure S4
% 
% Hamilton, Oganian, and Chang
% 
% Speech selectivity and comparison of STRF models

%% load data
addpath(genpath('../util1'));
heschl_load_data;


%% -------- Figure 4  - Aug 2020

%% A - on brain: electrodes by best model

figure
cla;
hold on
c_h = ctmr_gauss_plot(mniMesh.temporal.lh,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
loc_view(-90,60);
hold on
elsize =60;% (allr.vcorrs(cel,plMod(2))+.0001).^2*300;%/max(strfall(i,plMod(1)).vcorrs(cel).^2);
clear cpl
for i = 1:5
    cpl(i)=scatter3(allr.elecmatrix(elGroups{i},1),...
        allr.elecmatrix(elGroups{i},2),allr.elecmatrix(elGroups{i},3),elsize,groupColor(i,:), 'filled', 'MarkerFaceAlpha', .8);
end

legend(cpl,allr.elGroupName(1:5));
title('Speech feature encoding')

% print(fullfile(figDir,'mniBrain_different_features_v3.jpg'), '-djpeg', '-painters', '-r600')
%% B - scatter: Unique variance for onset and features vs spect 
figure
hfig=gcf; hfig.Units = 'centimeters'; hfig.Position(3:4) = [42,10];

allmods = {{'onset','spect_zscore'} ...    
    {'onset_phnfeaturesonset_peakRate','spect_zscore'} ...    
%     {'ons+phnF+relF0+drelF0+absF0+pR','ons+phnF+pR'}  ...        
    %     {'ons+phnF+relF0+drelF0+absF0+pR','spectZ'}  ...
%     {'ons+phnF+relF0+drelF0+pR' 'relF0+drelF0'} ...
    %     {'ons+phnF+relF0+drelF0+absF0+pR','spectZ'} ...
    %     {'ons+phnF+pR','ons+phnF+relF0+drelF0+absF0+pR'} ...
    %     {'ons+phnF+relF0+drelF0+absF0+pR', 'ons+phnF+relF0+drelF0+pR'} ...
    %     {'ons+phnF+relF0+drelF0+absF0+pR', 'ons+phnF+absF0+pR'} ...
    };

modColors = [onsetColor; rampColor; pitchColor];

for i = 1:length(allmods)
    modn  = allmods{i};
    subplot(1,length(allmods),i)
    cla
    hold on
    plMod=[];
    for  j = 1:2
        plMod(j) = find(ismember(allr.strfnames,modn{j}));
    end   
    
    cxval = allr.vcorrs(:, plMod(1)).^2;
    cyval = allr.vcorrs(:, plMod(2)).^2; 
       
    inclel = AnalysisEl;%setdiff(find(AnalysisEl),find(allr.elGroup(:,1)));
    uarea = unique(allr.(useAnaName)(inclel));
    scatter(cxval(inclel), cyval(inclel),20, 'k', 'filled');
    
    inclel = AnalysisEl&logical(allr.elGroup(:,i));%    find(AnalysisEl & cxval > cyval);
    uarea = unique(allr.(useAnaName)(inclel));
    scatter(cxval(inclel), cyval(inclel),40, modColors(i,:), 'filled');
    
    axis equal; 
    xlim([0 .55]), xticks(0:.25:.5); xlabel(allr.strfnames(plMod(1)), 'Interpreter', 'none')
    ylim([0 .55]), yticks(0:.25:.5); ylabel(allr.strfnames(plMod(2)), 'Interpreter', 'none')    
    title('R^2')
    rfl=refline(1,0);     rfl.LineStyle = ':';     rfl.Color = 'k';     rfl.DisplayName='';     rfl.LineWidth = 2;    
end
% print(fullfile(figDir, [modn{1} '_VS_' modn{2}]), '-depsc');%, '-opengl', '-r600');%, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 5 10])
%% C histograms of pitch and feature electrodes
elloc=[];
cgr=1:5;
for i =1:length(cgr)
    uEl{i} = find(allr.elGroup(:,cgr(i)));
    [ccount,~,~,clabels] = crosstab(allr.(useAnaName)(uEl{i}));
    elloc(str2num(cell2mat(clabels)),i) = ccount;
end

figure
subplot(1,2,1)  % pitch
cpl=bar(elloc(:,3:5)./crosstab(allr.(useAnaName)(allr.allEls==1)));
for i =1:3, cpl(i).FaceColor = groupColor(i+2,:);end
legend(allr.elGroupName(3:5));

subplot(1,2,2)
cpl=bar(elloc(:,1:2)./crosstab(allr.(useAnaName)(allr.allEls==1)));
for i =1:2, cpl(i).FaceColor = groupColor(i,:);end
legend(allr.elGroupName(1:2));
for i = 1:2
    subplot(1,2,i)
    xticks(1:5); xticklabels(includeAreasName)
    ylabel('portion electrodes'); ylim([0 0.5]);
    set(gca, 'XTickLabelRotation', 30);
    set(gca, 'FontName', 'Helvetica');
    box off
    legend off;
end

%% D silhouette index

elloc = allr.elecmatrix(allr.allEls,:);
elClust = zeros(size(allr.allEls));
cnames = {'onset', 'features', 'relP', 'absP'};
for i= 1:4
    elClust(allr.elGroup(allr.allEls,i)==1)=i;
end

% calculate silhouette index
silind=cell(6,1);
idx=1;
% pairwise comparison of all 4 clusters
for i = 1:4
    for j = i+1:4
        silind{idx}(:,1) = silhouette(elloc(ismember(elClust, [i j]),:), elClust(ismember(elClust, [i j])));
        silind{idx}(:,2) =  elClust(ismember(elClust, [i j]));
        idx=idx+1;
    end
end
%% D shuffle electrodes to get significance
nperm = 1000;
tic
silindPerm = cell(6,1);
idx=1;
for i = 1:4
    for j = i+1:4
        celloc = elloc(ismember(elClust, [i j]),:);
        celClust  = elClust(ismember(elClust, [i j]));
        for cperm = 1:nperm            
            cord = randperm(length(celloc));
            silindPerm{idx}(:,cperm) = silhouette(celloc(cord,:) , celClust);
        end
        idx=idx+1;
    end
end
toc

% get significance

for idx = 1:6
    silindMean(idx,:) = grpstats(silind{idx}(:,1), silind{idx}(:,2), 'median');
    ccl = unique(silind{idx}(:,2));
    permmeans=[];
    for cc = 1:2
        permmeans(cc,:) = median(silindPerm{idx}(silind{idx}(:,2)==ccl(cc),:),1);
    end
    for cc = 1:2 % pvalues for first/second cluster in each comparison
        silindp(idx,cc) = 1-sum(permmeans(cc,:)<silindMean(idx,cc))/nperm;
    end
end
%% Figure S4 - plot silhouette index results
figure
for i = 1:6
    subplot(1,6,i)
    cla;hold on;
    a = grpstats(silind{i}(:,1),silind{i}(:,2));
    aClust = unique(silind{i}(:,2));
    for j = 1:2        
        bar(j, a(j), 'FaceColor', groupColor(aClust(j),:), 'FaceAlpha', .5);        
    end
    
    xticks(1:2), xticklabels(cnames(unique(silind{i}(:,2))));    set(gca, 'xticklabelrotation', 30);    
    for kl=1:2
        if silindp(i,kl)<.05
            text(kl, .2, '*', 'FontSize', 20, 'HorizontalAlignment', 'center');
        else
            text(kl, .2, sprintf('%.2f',silindp(i,kl)), 'FontSize', 10, 'HorizontalAlignment', 'center');
        end
    end
    set(gca, 'FontName', 'Helvetica');
    set(gca, 'FontSize', 10);
    yticks(0:.1:.3); ylim([-.01 .37]);
    if i >1, yticks([]);end
end
suptitle('Pairwise silhouette distances')


%% --------- further analyses - part of text but not plots. ----
%% --- peakRate vs phonetic features
%% feature and peakRate R^2 scatter plot
featureEl = find(allr.elGroup(:,2)==1);
figure; 
cla;hold on;
gscatter(allr.uvar(featureEl,2), allr.uvar(featureEl,3), allr.(useAnaName)(featureEl), ...
    areaCols(unique(allr.(useAnaName)(featureEl)),:))
cel = featureEl(allr.feat.featpval(featureEl)<minsign & allr.feat.featu(featureEl)>0);
scatter(allr.uvar(cel,2), allr.uvar(cel,3), 'k')

cel = featureEl(allr.feat.peakRpval(featureEl)<minsign & allr.feat.peakRu(featureEl)>0);
scatter(allr.uvar(cel,2), allr.uvar(cel,3), 'r')


refline(1,0);
horzline(0); vertline(0);
legend(includeAreasName(unique(allr.(useAnaName)(featureEl))))
xlabel([allr.uvNames{2} ' R^2'] )
ylabel([allr.uvNames{3} ' R^2'])
box off
title('Unique variance explained')
%% plot feature and peakRate electrodes on brain

figure
cla;
hold on
c_h = ctmr_gauss_plot(mniMesh.temporal.lh,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
loc_view(-90,60);
hold on
elsize =60;% (allr.vcorrs(cel,plMod(2))+.0001).^2*300;%/max(strfall(i,plMod(1)).vcorrs(cel).^2);
clear cpl


% % peakR electrodes
peakREl = intersect(featureEl, find(allr.(useAnaName)>0 & allr.elecmatrix(:,1)<-30 & allr.feat.peakRpval<minsign & allr.feat.peakRu>0));
featEl = intersect(featureEl, find(allr.(useAnaName)>0 & allr.elecmatrix(:,1)<-30 & allr.feat.featpval<minsign & allr.feat.featu>0));

% peakRate only
cel2 = setdiff(peakREl, featEl);
cpl(1)=scatter3(allr.elecmatrix(cel2,1),...
    allr.elecmatrix(cel2,2),allr.elecmatrix(cel2,3),elsize,[71 180 108]/256, 'filled', 'MarkerFaceAlpha', .8);

% feature only
cel2 = setdiff(featEl,peakREl);
cpl(2)=scatter3(allr.elecmatrix(cel2,1),...
    allr.elecmatrix(cel2,2),allr.elecmatrix(cel2,3),elsize,[71 143 180]/256, 'filled', 'MarkerFaceAlpha', .8);

% % shared
cel2 = intersect(featEl,peakREl);
cpl(3)=scatter3(allr.elecmatrix(cel2,1),...
    allr.elecmatrix(cel2,2),allr.elecmatrix(cel2,3),elsize,rampColor,'filled', 'MarkerFaceAlpha', .8);

legend(cpl,'peakRate only', 'feature only', 'both');

title('Phonetic feature and peakRate encoding')
%% --- unique variance correlations
figure, 
cla;hold on;
scatter(allr.uvar(AnalysisEl,3),allr.uvar(AnalysisEl,5));
scatter(allr.uvar(allr.feat.featpval<minsign,3),allr.uvar(allr.feat.featpval<minsign,5), 'filled');
scatter(allr.uvar(allr.pitch.relpval<minsign,3),allr.uvar(allr.pitch.relpval<minsign,5),'filled');
horzline(0); vertline(0);
lsline();
xlabel('features'); ylabel('rel pitch');
%% % joint and separate encoding

elcounts = crosstab(allr.feat.featpval(AnalysisEl)<minsign & allr.feat.featu(AnalysisEl)>0,...
    allr.pitch.relpval(AnalysisEl)<minsign& allr.pitch.relRu(AnalysisEl)>0);

elcounts./crosstab(allr.feat.featpval(AnalysisEl)<minsign & allr.feat.featu(AnalysisEl)>0)

elcounts./crosstab(allr.pitch.relpval(AnalysisEl)<minsign& allr.pitch.relRu(AnalysisEl)>0)'
%%
figure, 
scatter(allr.uvar(AnalysisEl,2),allr.uvar(AnalysisEl,5));
horzline(0); vertline(0);
lsline();
xlabel('peakRate'); ylabel('rel pitch');

figure, 
scatter(allr.uvar(AnalysisEl,2),allr.uvar(AnalysisEl,3));
horzline(0); vertline(0);
lsline();
xlabel('peakRate'); ylabel('features');

