%% load data
addpath(genpath('../util1'));
heschl_load_data;


%% ----- Figure 3: feature matrix and example receptive fields
%% Figure 3A - feature matrix
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
%% Figure 3B - example electrode receptive fields

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
%% Figure 3B - example electrodes on the brain
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
scatter3(elloc(:,1),elloc(:,2),elloc(:,3), 60, areaCols(allr.(useAnaName)(plotel),:), 'filled')
text(elloc(:,1)+2,elloc(:,2),elloc(:,3), num2str([1:5]'),'FontSize',20,'FontWeight', 'normal', 'Color', 'k');
%print(fullfile(figDir,'mniBrain_exel_Fig3.jpg'), '-djpeg', '-painters', '-r600')
