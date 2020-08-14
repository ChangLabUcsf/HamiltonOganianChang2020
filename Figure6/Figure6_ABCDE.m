% Figure 6 Code panels A-E
% Hamilton, Oganian, and Chang (2020).

% Set the data directories
config_paths;
nmf_dir = sprintf('%s/NMF/',paper_data_dir);
strf_dir = sprintf('%s/STRFs/',paper_data_dir);

% Load the data for cluster assignments
load(sprintf('%s/allsubj_clustdata_avgtrial.mat', nmf_dir));
load(sprintf('%s/allsubj_NMF_avgtrial.mat', nmf_dir));
load(sprintf('%s/elecs_bestsubband_4bases.mat', nmf_dir));

addpath('/Users/liberty/Documents/UCSF/Heschls/code/cbrewer');
cmap  = cbrewer('qual','Dark2',9);

nF = FLs{4};

% Order subbands as onset/sust HG, onset only, sustained early, sustained
% late (the original NMF order is arbitrary)
sb_order= [1 3 4 2]; 

rng(1);
%% Figure 6A-B : latencies and average wts plots for each cluster

all_wts = {};
all_corrs = {};
for sb=1:4
    all_wts{sb} = [];
    all_corrs{sb} = [];
end
for s=1:length(subjects)    
    fname = sprintf('%s/spect_zscore/%s_STRF_spect_zscore.hf5',strf_dir, subjects{s}) %samealpha
    wts = h5read(fname,'/wts');
    vcorrs = h5read(fname,'/vcorrs');
    
    for sb=1:4
        chans_to_avg = elecs{sb}(elecs{sb}(:,1)==s,2);
        all_wts{sb} = [all_wts{sb}; wts(chans_to_avg,:)];
        all_corrs{sb} = [all_corrs{sb}; vcorrs(chans_to_avg)];
    end
end

cmap  = cbrewer('qual','Dark2',9);

figure;
all_lat = [];
all_sb = [];
for sb=1:length(sb_order)
    % Reshape weights to be nchans x nfeats x ndelays
    wtt=reshape(all_wts{sb_order(sb)}(:,2:end),size(all_wts{sb_order(sb)},1),[],60);
    
    % Plot latencies
    subplot(1,2,1);
    % Get the latency of each STRF
    lat = argmax(squeeze(median(wtt,2))');
    all_lat = [all_lat; lat'];
    all_sb = [all_sb; repmat(sb, length(lat),1)];
    jitt = 0.1*randn(length(lat),1); % Add some jitter
    plot(lat, sb+jitt, '.', 'color', cmap(sb_order(sb),:));
    hold all;
    
    set(gca,'xtick',[0:10:60], 'xticklabel', [0:0.1:0.6]);
    xlabel('Time Delay (s)','fontsize',10);
    set(gca,'ytick',[1 2 3 4]);
    ylabel('Cluster #','fontsize',10);
    
    subplot(1,2,2);
    onset_mean=squeeze(mean(mean(wtt,2),1));
    onset_err=squeeze(std(mean(wtt,2),[],1))/sqrt(size(wtt,1));
    transp=1;
    shadedErrorBar(1:60, onset_mean, onset_err,{'color',cmap(sb_order(sb),:)}, transp);
    hold all;
    set(gca,'xtick',[0:10:60], 'xticklabel', [0:0.1:0.6]);
    xlabel('Time Delay (s)','fontsize',10);
    axis tight;
    ylabel('Avg. \beta','fontsize',10);
end
subplot(1,2,1);
boxplot(all_lat,all_sb,'orientation','horizontal','plotstyle','compact','colors',cmap(sb_order,:)); hold all;
axis([0 60 0.5 4.5]);

%print_quality_fig(gcf,sprintf('/Users/liberty/Dropbox/Heschls_STRFs/data/figures/Figure6_AB.eps'),10,6,3,'inches','epsc');

% Statistics
[p,anovatab,stats] = kruskalwallis(all_lat,all_sb);
fprintf(1,'p = %2.2g, df=%.2f, chi-sq=%.2f\n', p, anovatab{2,3}, anovatab{2,5});
[c,m,h,nms] = multcompare(stats) %,'display','off'); % Compares average ranks of columns

%% Figure 6C - Spectrotemporal vs. Onset only R^2

fkeys = {'fricative', 'labial', 'nasal', 'back', 'sonorant', 'high', ...
    'voiced', 'coronal', 'front', 'dorsal', 'plosive', 'low',...
    'syllabic', 'obstruent'};
fkeys_order = [5 14 7 4 9 12 6 10 8 2 13 11 1 3];
fkeys = {fkeys{fkeys_order}};

cmap  = cbrewer('qual','Dark2',9);

feat_set1 = 'onset';
feat_set2 = 'spect_zscore';

all_corrs1 = {};
all_corrs2 = {};
for sb=1:4
    all_corrs1{sb} = [];
    all_corrs2{sb} = [];
end
for s=1:length(subjects)

    fname1 = sprintf('%s/onset/%s_STRF_%s.hf5',strf_dir, subjects{s}, feat_set1) %samealpha
    vcorrs1=h5read(fname1,'/vcorrs');
    
    fname2 = sprintf('%s/spect_zscore/%s_STRF_%s.hf5',strf_dir, subjects{s}, feat_set2) %samealpha
    vcorrs2=h5read(fname2,'/vcorrs');
    
    for sb=1:4
        chans_to_avg = elecs{sb}(elecs{sb}(:,1)==s,2);
        all_corrs1{sb} = [all_corrs1{sb}; vcorrs1(chans_to_avg)]
        all_corrs2{sb} = [all_corrs2{sb}; vcorrs2(chans_to_avg)];
    end
end

figure;
%subplot(1,2,1);
for sb=1:length(sb_order)
    plot(all_corrs1{sb_order(sb)}.^2, all_corrs2{sb_order(sb)}.^2, '.', 'markersize',20, 'color', cmap(sb_order(sb),:)); hold all;
    xlabel('R^2 onset only');
    ylabel('R^2 spectrotemporal');
end
plot([0 0.4], [0 0.4], 'k--');
axis square;
set(gcf,'color','w');
%print_quality_fig(gcf,'/Users/liberty/Dropbox/Heschls_STRFs/data/figures/Figure6C.eps',10,3,3,'inches','epsc');

%% Figure 6D: Get absolute vs. relative pitch vs both for clusters 1 - 4 as a bar chart
% bars show percentage of absolute pitch electrodes in a cluster

pthresh = 0.01; 
feat_set1 = 'absolute pitch';
feat_set2 = 'relative pitch';

% Unique absolute pitch
feat_set_pairs{1} = {'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_abs_f0_peakRate', 'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_peakRate','unique abs pitch'};

% Unique relative pitch
feat_set_pairs{2} = {'onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_abs_f0_peakRate', 'onset_phnfeaturesonset_abs_f0_peakRate','unique rel pitch'};

all_corrs1 = {};
all_corrs2 = {};
for sb=1:4
    all_corrs1{sb} = [];
    all_corrs2{sb} = [];
end
rel_abs_both = zeros(4,3);

for s=1:length(subjects)
    %try
    fname11 = sprintf('%s/%s/%s_STRF_%s.hf5',strf_dir, feat_set_pairs{1}{1},subjects{s}, feat_set_pairs{1}{1}); %samealpha
    vcorrs11=h5read(fname11,'/vcorrs');
    fname12 = sprintf('%s/%s/%s_STRF_%s.hf5',strf_dir, feat_set_pairs{1}{2}, subjects{s}, feat_set_pairs{1}{2}); %samealpha
    vcorrs12=h5read(fname12,'/vcorrs');
    vcorrs1 = vcorrs11.^2 - vcorrs12.^2;
    
    fname22 = sprintf('%s/%s/%s_STRF_%s.hf5',strf_dir, feat_set_pairs{2}{2}, subjects{s}, feat_set_pairs{2}{2}); %samealpha
    vcorrs22=h5read(fname22,'/vcorrs');
    vcorrs2 = vcorrs11.^2 - vcorrs22.^2;
    
    % Get permuted values
    permabs_dir = sprintf('%s/phnfeaturesonset_relative_log_f0_delta_relative_log_f0_permabs_f0_peakRate', strf_dir);
    permrel_dir = sprintf('%s/phnfeaturesonset_permrelative_log_f0_delta_relative_log_f0_abs_f0_peakRate', strf_dir);
    fnameperm_abs = sprintf('%s/%s_STRF_onset_phnfeaturesonset_relative_log_f0_delta_relative_log_f0_permabs_f0_peakRate.hf5', permabs_dir, subjects{s});
    fnameperm_rel = sprintf('%s/%s_STRF_onset_phnfeaturesonset_permrelative_log_f0_delta_relative_log_f0_abs_f0_peakRate.hf5', permrel_dir, subjects{s});
    
    perm_abs_corrs = h5read(fnameperm_abs,'/perm_corrs');
    perm_rel_corrs = h5read(fnameperm_rel,'/perm_corrs');
    
    pval_abs = 1-sum((repmat(vcorrs11', size(perm_abs_corrs,1),1)-perm_abs_corrs)>0)/1000;
    pval_rel = 1-sum((repmat(vcorrs11', size(perm_rel_corrs,1),1)-perm_rel_corrs)>0)/1000;

    for sb=1:4
        chans_to_avg = elecs{sb}(elecs{sb}(:,1)==s,2);
        for c=1:length(chans_to_avg)
            ch = chans_to_avg(c);
            if ((pval_abs(ch)<pthresh) && (pval_rel(ch)<pthresh))
                fprintf(1,'Subj %s, elec %d, pval_abs = %3.3d, pval_rel=%3.3d\n', subjects{s}, ch, pval_abs(ch), pval_rel(ch));
                rel_abs_both(sb,3) = rel_abs_both(sb,3)+1;
            elseif pval_abs(ch)<pthresh
                fprintf(1,'Subj %s, elec %d, pval_abs = %3.3d, pval_rel=%3.3d\n', subjects{s}, ch, pval_abs(ch), pval_rel(ch));
                rel_abs_both(sb,2) = rel_abs_both(sb,2)+1;
            elseif pval_rel(ch)<pthresh
                fprintf(1,'Subj %s, elec %d, pval_abs = %3.3d, pval_rel=%3.3d\n', subjects{s}, ch, pval_abs(ch), pval_rel(ch));
                rel_abs_both(sb,1) = rel_abs_both(sb,1)+1;
            else
                fprintf(1,'######Subj %s, elec %d, pval_abs = %3.3g, pval_rel=%3.3g\n', subjects{s}, ch, pval_abs(ch), pval_rel(ch));
            end
        end
    end
end

figure; 
bar(rel_abs_both./sum(rel_abs_both,1));
xlabel('Cluster');
ylabel('% electrodes');
legend('rel pitch','abs pitch','both');
%print_quality_fig(gcf,'/Users/liberty/Dropbox/Heschls_STRFs/data/figures/Figure6D.eps',10,3,3,'inches','epsc');

%% Figure 6E - Get onsets vs. peakRate tuning for clusters 1 - 4 as a bar chart
% bars show percentage of absolute pitch electrodes in a cluster

feat_set1 = 'onsets';
feat_set2 = 'peakRate+features';

% Unique onset
feat_set_pairs{1} = {'onset_phnfeaturesonset_peakRate', 'phnfeaturesonset_peakRate','unique onset'};

% Unique peakRate+features
feat_set_pairs{2} = {'onset_phnfeaturesonset_peakRate', 'onset', 'unique peakRate+features'};

all_corrs1 = {};
all_corrs2 = {};
for sb=1:4
    all_corrs1{sb} = [];
    all_corrs2{sb} = [];
end
phnpR_onset_both = zeros(4,2);

for s=1:length(subjects)
    fname11 = sprintf('%s/%s/%s_STRF_%s.hf5',strf_dir, feat_set_pairs{1}{1}, subjects{s}, feat_set_pairs{1}{1}); %samealpha
    vcorrs11=h5read(fname11,'/vcorrs');
    fname12 = sprintf('%s/%s/%s_STRF_%s.hf5',strf_dir, feat_set_pairs{1}{2}, subjects{s}, feat_set_pairs{1}{2}); %samealpha
    vcorrs12=h5read(fname12,'/vcorrs');
    vcorrs1 = vcorrs11.^2 - vcorrs12.^2;
    
    fname22 = sprintf('%s/%s/%s_STRF_%s.hf5',strf_dir, feat_set_pairs{2}{2}, subjects{s}, feat_set_pairs{2}{2}); %samealpha
    vcorrs22=h5read(fname22,'/vcorrs');
    vcorrs2 = vcorrs11.^2 - vcorrs22.^2;
    
    % Get permuted values
    permpR_dir = sprintf('%s/onset_phnfeaturesonset_permpeakRate', strf_dir);
    permphn_dir = sprintf('%s/onset_permphnfeaturesonset_peakRate', strf_dir);
    permonset_dir = sprintf('%s/permonset_phnfeaturesonset_peakRate', strf_dir);
    fnameperm_pR = sprintf('%s/%s_STRF_onset_phnfeaturesonset_permpeakRate.hf5', permpR_dir, subjects{s});
    fnameperm_phn = sprintf('%s/%s_STRF_onset_permphnfeaturesonset_peakRate.hf5', permphn_dir, subjects{s});
    fnameperm_onset = sprintf('%s/%s_STRF_permonset_phnfeaturesonset_peakRate.hf5', permonset_dir, subjects{s});
    
    perm_pR_corrs = h5read(fnameperm_pR,'/perm_corrs');
    perm_phn_corrs = h5read(fnameperm_phn,'/perm_corrs');
    perm_onset_corrs = h5read(fnameperm_onset,'/perm_corrs');
    
    pval_pR = 1-sum((repmat(vcorrs11', size(perm_pR_corrs,1),1)-perm_pR_corrs)>0)/1000;
    pval_phn = 1-sum((repmat(vcorrs11', size(perm_phn_corrs,1),1)-perm_phn_corrs)>0)/1000;
    pval_onset = 1-sum((repmat(vcorrs11', size(perm_onset_corrs,1),1)-perm_onset_corrs)>0)/1000;

    for sb=1:4
        chans_to_avg = elecs{sb}(elecs{sb}(:,1)==s,2);
        for c=1:length(chans_to_avg)
            ch = chans_to_avg(c);
            if pval_pR(ch)<pthresh && pval_phn(ch)<pthresh
                fprintf(1,'Subj %s, elec %d, pval_abs = %3.3g, pval_rel=%3.3g\n', subjects{s}, ch, pval_pR(ch), pval_onset(ch));
                phnpR_onset_both(sb,2) = phnpR_onset_both(sb,2)+1;
            elseif pval_onset(ch)<pthresh
                fprintf(1,'Subj %s, elec %d, pval_abs = %3.3g, pval_rel=%3.3g\n', subjects{s}, ch, pval_pR(ch), pval_onset(ch));
                phnpR_onset_both(sb,1) = phnpR_onset_both(sb,1)+1;
            else
                fprintf(1,'######Subj %s, elec %d, pval_abs = %3.3g, pval_rel=%3.3g\n', subjects{s}, ch, pval_pR(ch), pval_onset(ch));
            end
        end
    end
end

figure; 
barh=bar(phnpR_onset_both./sum(phnpR_onset_both,1));
xlabel('Cluster');
ylabel('proportion of electrodes');
cmap = [85 182 71;
       179 70 142]./255;
set(barh(1),'FaceColor',cmap(1,:));
set(barh(2),'FaceColor',cmap(2,:));

legend('onsets','peakRate+phnFeat');
%print_quality_fig(gcf,'/Users/liberty/Dropbox/Heschls_STRFs/data/figures/Figure6E.eps',10,3,3,'inches','epsc');
