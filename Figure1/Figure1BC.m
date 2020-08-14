%Figure 1BC
%
% Hamilton, Oganian, and Chang
% 
% Electrode counts from each area and electrodes on MNI brain

%% load data
addpath(genpath('../util1'));
heschl_load_data;
%% B - histogram of electrodes by area
figure
elcounts = crosstab(allr.(useAnaName)(AnalysisEl));
cla; hold on;
for i = 1:5
    bar(i, elcounts(i),  'FaceColor', areaCols(i,:));
end
xticks(1:5); xticklabels(includeAreasName); set(gca, 'xtickLabelRotation', 30);
yticks(0:50:200); ylabel('electrode count');ylim([0 max(elcounts)+25])
text(1:5, elcounts+20, num2str(elcounts))
box off

%% C - plot temporal electrodes on mni brains
figure

hold on
c_h = ctmr_gauss_plot(mniMesh.temporal.lh,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
cel = AnalysisEl & allr.elecmatrix(:,1)<-30;
celsize = min(max(5, (allr.maxr(cel).^2+10^(-10))*300), 75);

loc_view(-90,60)
scatter3(allr.elecmatrix(cel,1),...
    allr.elecmatrix(cel,2),allr.elecmatrix(cel,3),celsize, areaCols(allr.(useAnaName)(cel),:), 'filled')

hold on
gscatter(zeros(5,1),nan(5,1),includeAreasName', areaCols,'.',20);
% legend(['','',includeAreasName])
% cb=area_colorbar(areaCols,includeAreasName);
sgtitle('size ~ maximal receptive field model R^2 (min: .02, max: .5)')

%print(fullfile(figDir,'fig1E_brainExEl.jpg'), '-djpeg', '-painters', '-r800')
