% Figure4 (formerly figures 3 and 4)
%
% Hamilton, Oganian, and Chang
%
%% load data
addpath(genpath('../util1'));
heschl_load_data;
%% colors and anatomy definitions
useAnaName = 'CustomAna7area';
useAreaCols = area7Cols;
useAreaNames= new7AreaNames;
%% ----- Figure 4: feature matrix and example receptive fields
%% Figure 4A - feature matrix
exfeat = getfield(load('clockwork_features.mat'), 'out');
exfeat2 = load('pitch_matrices_clockwork.mat');
sentdet2 = sentdet;
sentdet2.sentdet = sentdet2.sentdet(13);
sentdur = length(exfeat.sound)/exfeat.soundf;
figure
subplot(6,1,1) % sound waveform
xax = (1:length(exfeat.sound))/exfeat.soundf-.5;
plot(xax, exfeat.sound,'color', .5*ones(1,3));
xlim([0, sentdur-1]);
box off,xticks([]),yticks([])
axis off
title('it had gone like clock work')
a=find(sum(sentdet2.sentdet.phnmatonset,1))/100-.5;
for cphn = 1:length(sentdet2.sentdet.phnnames)
    text(a(cphn), -1, sentdet2.sentdet.phnnames{cphn});
end
% phonemes

spl=subplot(6,1,2); % spectrogram
imagesc(exfeat.aud);
colormap(spl, cool);
set(gca, 'ydir', 'normal');
xlim([50, sentdur*100-50]);
box off,xticks([]),yticks([])
% axis off
xticks([])
yticks([1 80])
ylabel('F (kHz)')
yticklabels([.1 8])
spl=subplot(12,1,5); % onset
imagesc(exfeat.sentence_onset);
colormap(spl, [1 1 1; onsetColor])
box off,xticks([]),yticks([])
xlim([50, sentdur*100-50]);
% axis off
yticks(1)
yticklabels('onset')
spl=subplot(12,1,6:9); % features
imagesc([ecog_norm(exfeat.peakRate + [0 exfeat.peakRate(1:end-1)]);exfeat.phnfeaturesonset]);
colormap(spl, [1 1 1; rampColor])
box off,xticks([]),yticks([])
xlim([50, sentdur*100-50]);
% axis off
xticks([])
yticks(1:15)
yticklabels(['peakRate';features.featnames(1:14)]);

spl=subplot(12,1,10); % absolute pitch
imagesc(exfeat.abs_f0);
colormap(spl, [1 1 1; corpColors(2,:)])
box off,xticks([]),yticks([])
xlim([50, sentdur*100-50]);
% axis off
yticks([1 10]); yticklabels([90 250])
ylabel(sprintf('abs. \n pitch \n (Hz)'))
set(gca, 'YDir', 'normal')

spl=subplot(12,1,11); % relative pitch
imagesc([exfeat.relative_f0]);
colormap(spl, [1 1 1; corpColors(1,:)])
box off,xticks([]),yticks([])
xlim([50, sentdur*100-50]);
% axis off
yticks([1 5 10]); yticklabels([-1.9, 0, 1.9])
ylabel(sprintf('rel. \n pitch \n (z-score)'))
set(gca, 'YDir', 'normal')

spl=subplot(12,1,12); % d/dt relative pitch
imagesc(exfeat2.delta_relative_pitch');
colormap(spl, [1 1 1; corpColors(1,:)])
box off,xticks([]),yticks([])
xlim([0, sentdur*100-100]);
% axis off
yticks([1 5 10]); yticklabels([-.4, 0, .3])
ylabel(sprintf('d/dt \n rel. pitch \n (z-score)'))
set(gca, 'YDir', 'normal')
%% Figure 4B - example electrode receptive fields

inclelid = find(AnalysisEl);
[a,b] = max(allr.uvar,[],1);

% b(4) = 941;
% b(5) = 263;%897;%906;
% b(1)= 1024;%67;%317;
b(4) = 1100;
%b(1) = 67;
plotel = b(1:5);

%plotel = inclelid([202,28,336,272,247]);%49 ]);
%plotel = inclelid([100,28,335,276,325]);%49 ]); 55  #2 = 279 #5=325
plotm = [5 5 5 4 7];
plelsid = allr.sid(plotel);
plelnum = allr.elnum(plotel);
eltargetnames = {'onset', 'peakRate', 'features', 'absPitch', 'relPitch'};

figure
for i = 1:length(plotel)
    subplot(1,length(plotel),i)
    cstrf = squeeze(allr.strf{plotm(i)}(plotel(i),[1 end 2:end-1],:));
    imagesc(cstrf/max(abs(cstrf(:))),[-.5 .5]);
    %     horzline([1.5, 15.5,15.5, 25.5, 35.5,45.5]',[],'w','-');
    %title(sprintf('R^2=%0.2f \n %s', allr.vcorrs(plotel(i),plotm(i)).^2, eltargetnames{i}));
    title(sprintf('%s uR^2=%.2f \n model R^2=%0.2f \n %s', eltargetnames{i}, allr.uvar(plotel(i),i), allr.vcorrs(plotel(i),plotm(i)).^2, eltargetnames{i}));
    %     ylabel(allr.strfnames{plotm(i)})
    horzline([1.5,2.5, 16.5, 26.5, 36.5]', [],'k','-');
    if i == 5
        yticks([1, 2, 9, 21, 31, 41])
        yticklabels({'onset', 'peakRate', 'phn feat', 'absP', 'relP', '?relP'});
        set(gca,'YAxisLocation','Right');
    else
        yticks([]);
    end
    %     colorbar;
    ylim([.5 47.5])
    set(gca, 'XDir', 'reverse');
    set(gca, 'XTickLabel', '')
end
colormap(strfcmap)
%% Figure 4C - example electrodes on the brain
figure
cla;
hold on
c_h = ctmr_gauss_plot(mniMesh.temporal.lh,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
loc_view(-90,60);
hold on
elsize =60;% (allr.vcorrs(cel,plMod(2))+.0001).^2*300;%/max(strfall(i,plMod(1)).vcorrs(cel).^2);
clear cpl
elloc = allr.elecmatrix(plotel,:); % elloc(1,2) = elloc(1,2) - 2;  % move el1 by a tiny bit so that e1 and e5 do not overlap
scatter3(elloc(:,1),elloc(:,2),elloc(:,3), 60, useAreaCols(allr.(useAnaName)(plotel),:), 'filled')
text(elloc(:,1)+2,elloc(:,2),elloc(:,3), num2str([1:5]'),'FontSize',20,'FontWeight', 'normal', 'Color', 'k');
%print(fullfile(figDir,'mniBrain_exel_Fig3.jpg'), '-djpeg', '-painters', '-r600')


%% Figure 4C - on brain: electrodes by best model

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
%% D+E - scatter: Unique variance for onset and features vs spect 
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
%% F+G histograms of pitch and feature electrodes
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
    xticks(1:length(useAreaNames)); xticklabels(useAreaNames)
    ylabel('portion electrodes'); ylim([0 0.5]);
    set(gca, 'XTickLabelRotation', 30);
    set(gca, 'FontName', 'Helvetica');
    box off
    legend off;
end

%% F+G version 2
elloc=[];
cgr=1:5;
for i =1:length(cgr)
    uEl{i} = find(allr.elGroup(:,cgr(i)));
    [ccount,~,~,clabels] = crosstab(allr.(useAnaName)(uEl{i}));
    elloc(str2num(cell2mat(clabels)),i) = ccount;
end

figure
subplot(1,2,1)  % pitch
cpl=bar(elloc(:,3:5)'./crosstab(allr.(useAnaName)(allr.allEls==1))');
for i =1:7, cpl(i).FaceColor = useAreaCols(i,:);end
legend(useAreaNames);
set(gca, 'xticklabels', allr.elGroupName(3:5));

subplot(1,2,2)
cpl=bar(elloc(:,1:2)'./crosstab(allr.(useAnaName)(allr.allEls==1))');
for i =1:7, cpl(i).FaceColor = useAreaCols(i,:);end
legend(useAreaNames);
set(gca, 'xticklabels', allr.elGroupName(1:2));

for i = 1:2
    subplot(1,2,i)
    xticks(1:length(useAreaNames)); xticklabels(useAreaNames)
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
