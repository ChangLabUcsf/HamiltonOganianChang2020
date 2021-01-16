%Figure 1BC
%
% Hamilton, Oganian, and Chang
% 
% Electrode counts from each area and electrodes on MNI brain

%% load data
addpath(genpath('../util1'));
heschl_load_data;
%% B - histogram of electrodes by area
include7AreasName = {'PT','pmHG','alHG','PP','pSTG onset','pSTG non-onset','mSTG'};
figure
elcounts = crosstab(allr.(useAnaName)(AnalysisEl));
cla; hold on;
for i = 1:length(elcounts)
    bar(i, elcounts(i),  'FaceColor', area7Cols(i,:));
end
xticks(1:length(elcounts)); xticklabels(include7AreasName); set(gca, 'xtickLabelRotation', 30);
yticks(0:50:200); ylabel('electrode count');ylim([0 max(elcounts)+25])
text(1:length(elcounts), elcounts+20, num2str(elcounts))
box off

print_quality_fig(gcf,sprintf('%s/figures/elec_counts.eps',pwd), 8, 3, 3, 'inches','epsc');
%% C - plot temporal electrodes on mni brains
figure

hold on
c_h = ctmr_gauss_plot(mniMesh.temporal.lh,[0 0 0],0,'lh');
c_h.FaceAlpha = .5;
cel = AnalysisEl & allr.elecmatrix(:,1)<-30;
celsize = min(max(5, (allr.maxr(cel).^2+10^(-10))*300), 75);

loc_view(-90,60)
scatter3(allr.elecmatrix(cel,1),...
    allr.elecmatrix(cel,2),allr.elecmatrix(cel,3),celsize, area7Cols(allr.(useAnaName)(cel),:), 'filled')

hold on
gscatter(zeros(length(elcounts),1),nan(length(elcounts),1),include7AreasName', area7Cols,'.',20);

% Show the onset electrodes on top to make them more visible
onset_cel = allr.(useAnaName)==5 & cel;
onset_celsize = min(max(5, (allr.maxr(onset_cel).^2+10^(-10))*300), 75);

scatter3(allr.elecmatrix(onset_cel,1)-0.5, ...
    allr.elecmatrix(onset_cel,2),allr.elecmatrix(onset_cel,3),onset_celsize, area7Cols(allr.(useAnaName)(onset_cel),:), 'filled')

% legend(['','',includeAreasName])
% cb=area_colorbar(areaCols,includeAreasName);
suptitle('size ~ maximal receptive field model R^2 (min: .02, max: .5)')

print(fullfile(figDir,'fig1E_brainExEl_7areas.jpg'), '-djpeg', '-painters', '-r800')
