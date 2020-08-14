% Figure2CDE
%
% Hamilton, Oganian, and Chang
%

if 0
    addpath(genpath('../util1'));
    heschl_load_data;
end

subj = 'S03';

strf_dir = '/Users/liberty/Documents/UCSF/Heschls/data/strfs';
load(sprintf('/Users/liberty/Dropbox/Heschls_STRFs/data/pure_tone/%s_RF.mat', subj),'all_RFs');

elecs = [268 267 266 265 53:56];

strf_file = '/Users/liberty/Dropbox/Heschls_STRFs/data/STRFs/spect_zscore/S03_STRF_0to600ms_100Hz_allchans_zscore_intercept_noedge_samealpha_consolidatedcode.hf5';

% Plot receptive fields (pure tone)
figure;
suptitle('Pure Tones');

for e=1:length(elecs)
    subplot(2,4,e);
    rfmat = all_RFs(:,:,elecs(e));
    imagesc(rfmat', [-max(abs(rfmat(:))) max(abs(rfmat(:)))]);
    text(1,70,num2str(e),'fontsize',12);

    set(gca,'yaxislocation','right');
    set(gca,'xtick',[1 2 3],'xticklabel',[0 10 20]);
    set(gca,'ytick',[12 44 80],'yticklabel',[0.5 2 8]);
    if e==length(elecs)
        
        ylabel('Freq (kHz)');
        xlabel('Int. (dB)');
    end
    axis xy;
end
colormap(flipud(cbrewer('div','PRGn',256)));

% Plot speech STRFs
figure;
suptitle('Speech STRF');
wts = h5read(strf_file,'/wts');
for e=1:length(elecs)
    subplot(2,4,e);
    strfmat = reshape(wts(elecs(e),2:end),[],60);
    imagesc(fliplr(strfmat), [-max(abs(strfmat(:))) max(abs(strfmat(:)))]);
    text(10,70,num2str(e),'fontsize',12);
    hold all;
    set(gca,'yaxislocation','right');
    set(gca,'xtick',[1 30 60],'xticklabel',[-0.6 -0.3 0]);
    set(gca,'ytick',[12 44 80],'yticklabel',[0.5 2 8]);

    if e==length(elecs)
        xlabel('Time');
        ylabel('Freq (kHz)');

    end
    axis xy;
end
colormap(flipud(cbrewer('div','PRGn',256)));

%%
pt_color = [0.17, 0.22, 0.58];
pp_color = [0.62, 0.22, 0.58];
hg_color = [0.06, 0.69, 0.29];
stg_color = [0.69, 0.70, 0.21];
pstg_color = [0.69, 0.12, 0.14];

%elecs = 53:56;
figure;
cs = [hg_color; hg_color; pt_color; pt_color; stg_color; stg_color; pstg_color; pstg_color];
for e=1:length(elecs)
    wts = h5read(strf_file,'/wts');
    subplot(2,4,e);
    rfmat = all_RFs(:,:,elecs(e));
    freq_tuning_pt = sum(rfmat);
    %freq_tuning_pt = medfilt1(freq_tuning_pt,11);
    freq_tuning_pt = sgolayfilt([zeros(10,1); freq_tuning_pt'; zeros(10,1) ], 3, 31);
	freq_tuning_pt = freq_tuning_pt(11:11+79); % Get rid of zero padding
    freq_tuning_pt = freq_tuning_pt/max(abs(freq_tuning_pt));
    
    plot(freq_tuning_pt,1:80,'color',cs(e,:),'linewidth',2);
    text(1,70,num2str(e));
    hold all;
    set(gca,'xdir','reverse');
    set(gca,'yaxislocation','right');
    
    strfmat = reshape(wts(elecs(e),2:end),[],60);
    freq_tuning = sum(strfmat,2);
    %freq_tuning = medfilt1(freq_tuning, 11);
    freq_tuning = sgolayfilt([zeros(10,1); freq_tuning; zeros(10,1) ], 3, 31);
	freq_tuning = freq_tuning(11:11+79); % Get rid of zero padding
    freq_tuning = freq_tuning/max(abs(freq_tuning));
    plot(freq_tuning,1:80,'color',cs(e,:),'linewidth',1);
    
    set(gca,'xtick',[]);
    if e==4 || e==8
        set(gca,'ytick',[12 44 80],'yticklabel',[0.5 2 8]);
    else
        set(gca,'yticklabel',[])
    end
    axis([-1 1.2 1 80]);
    axis square;
end

