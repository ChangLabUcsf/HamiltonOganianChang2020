% Figure 2BFGH, S2 and S3
%
% Hamilton, Oganian, and Chang
%
% Magnitude of pure tone response in each area

if 0
    addpath(genpath('../util1'));
    heschl_load_data;
end
%heschl_load_anatomy;
strf_dir = '/Users/liberty/Dropbox/Heschls_STRFs/data/STRFs/spect_zscore';

pure_tones_dir = '/Users/liberty/Dropbox/Heschls_STRFs/data/pure_tone';
pt_subjs = {'S03', 'S04', 'S07', 'S08', 'S09'};

auditory_anat = {'transversetemporal','planumpolare','planumtemporale','ctx_lh_G_temp_sup-Plan_polar',...
    'ctx_rh_G_temp_sup-Plan_polar','ctx_lh_S_temp_sup-Plan_polar','ctx_rh_S_temp_sup-Plan_polar',...
    'ctx_lh_G_temporal_transverse','ctx_rh_G_temporal_transverse',...
    'ctx_lh_S_temporal_transverse','ctx_rh_S_temporal_transverse',...
    'ctx_lh_G_temp_sup-Plan_tempo','ctx_rh_G_temp_sup-Plan_tempo',...
    'ctx_lh_S_temp_sup-Plan_tempo','ctx_rh_S_temp_sup-Plan_tempo','superiortemporal',...
    'pSTG','mSTG'};
planumpolare = {'planumpolare','ctx_lh_G_temp_sup-Plan_polar', 'ctx_rh_G_temp_sup-Plan_polar',...
    'ctx_lh_S_temp_sup-Plan_polar', 'ctx_rh_S_temp_sup-Plan_polar'};
planumtemporale = {'planumtemporale','ctx_lh_G_temp_sup-Plan_tempo', 'ctx_rh_G_temp_sup-Plan_tempo',...
    'ctx_lh_S_temp_sup-Plan_tempo', 'ctx_rh_S_temp_sup-Plan_tempo'};
transversetemporal = {'transversetemporal','ctx_lh_G_temporal_transverse','ctx_rh_G_temporal_transverse',...
    'ctx_lh_S_temporal_transverse','ctx_rh_S_temporal_transverse'};

pt_color = [0.17, 0.22, 0.58];
pp_color = [0.62, 0.22, 0.58];
mstg_color = [0.69, 0.70, 0.21];
pstg_color = [0.69, 0.12, 0.14];
hg_color = [0.06, 0.69, 0.29];

rf_dat = containers.Map;
all_pt_magnitudes = [];
all_sp_magnitudes = [];
all_pt_avg = [];
all_sp_avg = [];
all_anat = [];
all_elecs = [];
all_groups = [];
all_pvals_tones = [];
all_pvals_speech = [];
all_vcorrs = [];
all_rf = [];
all_rf_struct = struct();
all_tc_struct = struct();
for i=1:5
    all_rf_struct(i).rf = [];
    all_tc_struct(i).tc = [];
    all_tc_struct(i).tc_speech = [];
    all_tc_struct(i).vcorrs = [];
end

load('melfreq.mat')
freqs = melfreq(3:end);

for s=1:length(pt_subjs)
    subj = pt_subjs{s};
    fprintf(1,'Getting RF data for %s\n',subj)
    
    strf_file = sprintf('%s/%s_STRF_0to600ms_100Hz_allchans_zscore_intercept_noedge_samealpha_consolidatedcode.hf5', strf_dir, subj);
    vcorrs = h5read(strf_file, '/vcorrs');
    wts = h5read(strf_file, '/wts');
    
    sig_inds = find(vcorrs>0.1);
    
    rf_file = sprintf('%s/%s_RF.mat', pure_tones_dir, subj);
    this_rf = load(rf_file,'all_RFs');
    
    rfdat = load(rf_file,'rfdat');
    
    out = outall.(subj);
    for i=1:length(out)
        out(i).silence = out(i).resp(:, 1:50);
        out(i).speech = out(i).resp(:, 51:100); % 0.5 - 1 s
        out(i).resp = mean(out(i).resp(:,50:end-50,:),3);
    end
    resp = [out.resp];
    speech_pretrial_silence = [out.silence];
    speech_posttrial_speech = [out.speech];
    
    anatomy = imgNative.(subj).newAnatomy(:,4);
    elecmatrix = imgNative.(subj).elecmatrix;
    elecmatrix_warped = imgmni.(subj).elecmatrix;
    
    % Remove channels not in the anatomy file
    if strcmp(subj,'S07')
        anatomy = anatomy(1:378);
    elseif strcmp(subj, 'S08')
        anatomy = anatomy(1:160);
    elseif strcmp(subj, 'S09')
        anatomy = anatomy(1:274);
    end
    
    this_rf.all_RFs = this_rf.all_RFs(:,:,1:length(anatomy));
    
    %tdt_elecs = load(sprintf('/Applications/freesurfer_v5/subjects/%s/elecs/TDT_elecs_all_anatfix.mat',subj));
    %elecmatrix = tdt_elecs.elecmatrix;
    %anatomy = tdt_elecs.anatomy(:,4);
    
    for i = 1:length(anatomy)
        a = anatomy{i};
        g = 0;
        if ~all(all(this_rf.all_RFs(:,:,i)==0))
            if ismember(a, auditory_anat)
                %fprintf(1,'%s is part of auditory cortex\n', a);
                if ismember(a,'mSTG')
                    a = 'mSTG';
                    a_short = 'mSTG';
                    g=5;
                    clr = mstg_color;
                elseif ismember(a, 'pSTG')
                    a='pSTG';
                    a_short='pSTG';
                    clr = pstg_color;
                    g=4;
                elseif ismember(a,transversetemporal)
                    a = 'transversetemporal';
                    a_short = 'HG';
                    clr = hg_color;
                    g = 2;
                elseif ismember(a,planumpolare)
                    a = 'planumpolare';
                    clr = pp_color;
                    a_short = 'PP';
                    g = 3;
                elseif ismember(a, planumtemporale)
                    a = 'planumtemporale';
                    clr = pt_color;
                    a_short = 'PT';
                    g = 1;
                end
                
                mag = max(max(this_rf.all_RFs(:,:,i)));
                avg_mag_pt = mean(mean(this_rf.all_RFs(:,:,i)));
                mag_sp = max(resp(i,:));
                avg_mag_sp = mean(resp(i,:));
                
                %current_tuning_curve = squeeze(nanmean(this_rf.all_RFs(:,:,i),1));
                curr_rf = this_rf.all_RFs(:,:,i);
                %curr_rf = medfilt2(curr_rf, [3 3]); %LHLHLH
                current_tuning_curve = squeeze(nanmean(curr_rf,1));
                
                %current_tuning_curve = current_tuning_curve/max(current_tuning_curve);
                
                % Zero-pad and filter with Savitzky-Golay 3rd order filter
                tuning_curve = sgolayfilt([zeros(10,1); current_tuning_curve'; zeros(10,1) ], 3, 31);
                %tuning_curve = medfilt1([zeros(10,1); current_tuning_curve'; zeros(10,1) ], 11);
                %tuning_curve(tuning_curve<0) = 0;
                tuning_curve = tuning_curve(11:11+79); % Get rid of zero padding
                
                % Tuning curve for speech
                strfmat = reshape(wts(i,2:end),[],60);
                freq_tuning = sum(strfmat,2);
                %freq_tuning = medfilt1(freq_tuning, 11);
                freq_tuning = sgolayfilt([zeros(10,1); freq_tuning; zeros(10,1) ], 3, 31);
                freq_tuning = freq_tuning(11:11+79); % Get rid of zero padding
                freq_tuning(freq_tuning<0) = 0;
                freq_tuning = freq_tuning/max(abs(freq_tuning));
                
                debug = 0;
                [~, mm] = max(tuning_curve);
                
                if ~isnan(mag)
                    %if mag > 1
                    all_pt_magnitudes = [all_pt_magnitudes; mag];
                    all_pt_avg = [all_pt_avg; avg_mag_pt];
                    all_sp_magnitudes = [all_sp_magnitudes; mag_sp];
                    all_sp_avg = [all_sp_avg; avg_mag_sp];
                    
                    all_vcorrs = [all_vcorrs; vcorrs(i) i s];
                    
                    all_groups = [all_groups; g];
                    all_anat = {all_anat; a_short};
                    all_elecs = [all_elecs; elecmatrix_warped(i,:)];
                    
                    %fprintf(1,'Testing for significant pure tone response\n');
                    tt = rfdat.rfdat.tt;
                    % Take mean across repeats, intensities, and frequencies
                    curr_chan = nanmean(nanmean(nanmean(rfdat.rfdat.rf(i,tt>0.05 & tt<0.3,2:end,:,:),5),4),3);
                    silence = nanmean(nanmean(nanmean(rfdat.rfdat.rf(i,tt<=0.05,2:end,:,:),5),4),3);
                    [p,~] = ranksum(curr_chan(:), silence(:));
                    
                    all_pvals_tones = [all_pvals_tones; p];
                    
                    % Do the same thing for a comparable speech window
                    [p2,~] = ranksum(speech_posttrial_speech(i,:), speech_pretrial_silence(i,:));
                    all_pvals_speech = [all_pvals_speech; p2];
                    
                    
                    %fprintf(1,'group %d, anat %s, subj %s, channel %d\n', g, a, subj, i);
                    all_rf_struct(g).rf = cat(3,all_rf_struct(g).rf, this_rf.all_RFs(:,:,i));
                    all_tc_struct(g).tc = [all_tc_struct(g).tc tuning_curve];
                    all_tc_struct(g).tc_speech = [all_tc_struct(g).tc_speech freq_tuning];
                    all_tc_struct(g).vcorrs = [all_tc_struct(g).vcorrs vcorrs(i)];
                    all_rf = cat(3, all_rf, this_rf.all_RFs(:,:,i));
                    
                end
                
            end
        end
    end
end

%% Figure 2B
cmap = imread('RdBuPr_splinesqrt22.png');
cmap = flipdim(cmap,1);

all_pt_avg(all_pt_avg<0) = 0;
all_sp_avg(all_sp_avg<0) = 0;

all_pt_avg = all_pt_avg./max(all_pt_avg);
all_sp_avg = all_sp_avg./max(all_sp_avg);

inds2rm = find(all_elecs(:,1)>0);
all_sp_avg(inds2rm) = [];
all_pt_avg(inds2rm) = [];
all_elecs(inds2rm,:) = [];

x=round(all_sp_avg*255)+1; %red
y=round(all_pt_avg*255)+1; %blue

elec_colors = zeros(length(all_pt_avg),3);
for e = 1:length(all_pt_avg)
    elec_colors(e,:) = double(squeeze(cmap(x(e),y(e),:)))./256;
end

figure('units','normalized','outerposition',[0 0 1 1]);
temporal = mniMesh.temporal.lh;
ctmr_gauss_plot(temporal, [0 0 0], 1, 'lh');
el_add(all_elecs, 'color', elec_colors, 'msize', 15);
loc_view(-110,49);

figure;
imagesc(cmap)
axis square
axis xy
set(gca,'xtick',[1 255],'xticklabel',{'0','max'});
xlabel('Pure tone');
set(gca,'ytick',[1 255],'yticklabel',{'0','max'});
ylabel('Speech');
%print_quality_fig(gcf,'purplered_cmap.eps',10,4,4,'inches','epsc');


%%
% Figure 2F
anat_labels = {'PT','HG','PP','pSTG','mSTG'};

pt = ones(length(all_pt_magnitudes),1); % pure tone
sp = 1+ones(length(all_pt_magnitudes),1); % speech
measure_type = [pt; sp];

figure;
boxplot([all_pt_magnitudes; all_sp_magnitudes], [repmat(all_groups,2,1) measure_type], ...
    'colors', [pt_color; pt_color; hg_color;hg_color; pp_color; pp_color; ...
    pstg_color; pstg_color; mstg_color; mstg_color],'symbol','');

h = findobj(gca,'Tag','Box');
cs = [pt_color; pt_color; hg_color;hg_color; pp_color; pp_color; ...
    pstg_color; pstg_color; mstg_color; mstg_color];
for j=1:length(h)
    if mod(j,2)==0
        patch(get(h(j),'XData'),get(h(j),'YData'),cs(length(h)-(j-1),:),'FaceAlpha',1.00);
    else
        patch(get(h(j),'XData'),get(h(j),'YData'),cs(length(h)-(j-1),:),'FaceAlpha',.3);
    end
end

tags = {'Lower Whisker','Upper Whisker','Lower Adjacent Value',...
    'Upper Adjacent Value','Median','Box'};
for tt=1:length(tags)
    
    h = findobj(gca,'Tag',tags{tt});
    for j=1:length(h)
        set(h(j), 'linewidth',1,'color','k','linestyle','-');
    end
end
ylabel('Max Z-scored high gamma');

ylim([-0.1 12]);
set(gca,'ytick',[0:3:12]);
set(gca,'xtick',[1.5:2:9.5],'xticklabel', anat_labels);
%print_quality_fig(gcf,'/Users/liberty/Documents/UCSF/Heschls/data/Fig2F.eps',10, 3, 3, 'inches', 'epsc');

%%
% Stats for pure tone magnitudes
fprintf(1,'Pure tone magnitudes\n');
[P,ANOVATAB,STATS] = kruskalwallis(all_pt_magnitudes, all_groups)
[c,m,h,nms] = multcompare(STATS);

fprintf(1,'********************\n');
fprintf(1,'Pure tone magnitudes\n');
for i=1:size(c,1)
    fprintf(1,'%s vs %s: %.2f, [%.2f, %.2f], %.4f\n', anat_labels{c(i,1)}, anat_labels{c(i,2)}, c(i,4) ,c(i,3), c(i,5), c(i,6));
end

%%
% Stats for speech
[P,ANOVATAB,STATS] = kruskalwallis(all_sp_magnitudes, all_groups);
[c,m,h,nms] = multcompare(STATS);

fprintf(1,'********************\n');
fprintf(1,'Speech magnitudes\n');
for i=1:size(c,1)
    fprintf(1,'%s vs %s: %.2f, [%.2f, %.2f], %.4f\n', anat_labels{c(i,1)}, anat_labels{c(i,2)}, c(i,4) ,c(i,3), c(i,5), c(i,6));
end

%%
% percentages of each area with significant response to tones
anat_labels = {'PT','HG','PP','pSTG','mSTG'};
for i=1:5
    perc=sum(all_pvals_tones(all_groups==i)<0.05/length(all_pvals_tones))/length(find(all_groups==i)); % Bonferroni corrected
    fprintf(1,'%s : %2.2f percent significant response to tones \n', anat_labels{i}, perc*100);
end


%%
% percentages of each area with significant response to speech
anat_labels = {'PT','HG','PP','pSTG','mSTG'};
for i=1:5
    perc=sum(all_pvals_speech(all_groups==i)<0.05/length(all_pvals_speech))/length(find(all_groups==i)); % Bonferroni corrected
    fprintf(1,'%s : %2.2f percent significant response to speech\n', anat_labels{i}, perc*100);
end
%%

good_speech_inds=find(all_pvals_speech<0.05/length(all_pvals_speech));
good_tone_inds=find(all_pvals_tones<0.05/length(all_pvals_tones));
bad_speech_inds=find(all_pvals_speech>=0.05/length(all_pvals_speech));
bad_tone_inds=find(all_pvals_tones>=0.05/length(all_pvals_tones));

bad_speech_good_tone_inds = intersect(good_tone_inds, bad_speech_inds);
good_speech_good_tone_inds = intersect(good_tone_inds, good_speech_inds);
good_speech_bad_tone_inds = intersect(bad_tone_inds, good_speech_inds);

bad_sp_good_t=length(bad_speech_good_tone_inds)/length(all_pvals_speech)*100;
good_sp_good_t=length(good_speech_good_tone_inds)/length(all_pvals_speech)*100;
good_sp_bad_t=length(good_speech_bad_tone_inds)/length(all_pvals_speech)*100;

fprintf(1,'Not significant for speech, significant for tones: %3.1f percent\n', bad_sp_good_t);
fprintf(1,'Significant for speech, significant for tones: %3.1f percent\n', good_sp_good_t);
fprintf(1,'Significant for speech, not significant for tones: %3.1f percent\n', good_sp_bad_t);


%%
% Supplemental figures of RFs from each area
all_rf_mask = all_rf_struct;
sig_bw = [];
all_sig_grps = [];
group_names = {'PT','HG','PP','pSTG','mSTG'};
nrows = 11;
ncols = 13;
for grp=1:5
    plotnum = 1;
    figure;
    group_elecs = find(all_groups==grp);
    chans = size(all_rf_struct(grp).rf,3);
    
    ncols_to_use = ceil(chans/nrows);
    if nrows*ncols < chans
        nrows=nrows+1;
    end
    all_tc_struct(grp).pvals_rf_in_out = zeros(length(group_elecs),1);
    all_tc_struct(grp).zval_rf_in_out = zeros(length(group_elecs),1);
    all_tc_struct(grp).peaks = zeros(length(group_elecs),1);
    for c=1:length(group_elecs)
        tc_scaled = all_tc_struct(grp).tc(:,c)/max(abs(all_tc_struct(grp).tc(:,c)))*3;
        r_tc = round(tc_scaled);
        for f=1:80
            if r_tc(f)>0 & ~isnan(r_tc(f))
                all_rf_mask(grp).rf(1:r_tc(f),f,c) = 1;
                all_rf_mask(grp).rf(r_tc(f)+1:end,f,c) = NaN;
            else
                all_rf_mask(grp).rf(:,f,c) = NaN;
            end
        end
        
        subplot(nrows,ncols,plotnum);
        
        if mod(c, ncols_to_use) == 0 && c>1
            plotnum = plotnum + (ncols-ncols_to_use)+1;
        else
            plotnum = plotnum+1;
        end
        rf_data = flipud(all_rf_struct(grp).rf(:,:,c));
        imagesc(rf_data);
        caxis([-3 3]);
        
        hold all;
        plot(tc_scaled,'b','linewidth',3);
        [pks, locs, wd, prom] =findpeaks(tc_scaled,'MinPeakProminence',0.5,'MinPeakHeight',0.5);
        
        % Plot the peaks
        %         for pk=1:length(pks)
        %             plot([locs(pk) locs(pk)], [0 3.5], 'k');
        %         end
        
        all_tc_struct(grp).peaks(c) = length(pks);
        % Determine if inside RF is significantly > outside RF (spontaneous)
        in_rf = flipud(all_rf_struct(grp).rf(:,:,c)).*all_rf_mask(grp).rf(:,:,c);
        inside_rf = nanmean(nanmean(in_rf));
        
        outside_mask = isnan(all_rf_mask(grp).rf(:,:,c));
        out_rf = flipud(all_rf_struct(grp).rf(:,:,c)).*outside_mask;
        out_rf(out_rf==0) = NaN;
        outside_rf = nanmean(nanmean(out_rf));
        
        if all(all(isnan(out_rf))) || all(all(isnan(in_rf)))
            continue;
            p=1;
            stat=0;
        else
            [p,h,stat] = ranksum(out_rf(:), in_rf(:));
        end
        all_tc_struct(grp).pvals_rf_in_out(c) = p;
        all_tc_struct(grp).zval_rf_in_out(c) = stat.zval;
        
        if p<0.05/length(group_elecs)
            set(gca,'xcolor','r');
            set(gca,'ycolor','r');
            set(gca,'color','r');
            all_sig_grps = [all_sig_grps; grp];
        end
        
        set(gca,'xtick',[12 80], 'xticklabel', [0.5 8]);
        set(gca,'ytick',[]);
        colormap(flipud(cbrewer('div','PRGn',11)));
        
    end
    set(gca,'xtick',[12 80], 'xticklabel', [0.5 8]);
    set(gca,'ytick',[]);
    set(gcf,'color','w');
    suptitle(group_names{grp});
    %print_quality_fig(gcf,sprintf('RFs_%s.eps', group_names{grp}),10,8,8,'inches','epsc');
end

%%
% Figure 2G
% percentage of good RFs, and percentage of multi vs single peaked RFs
perc_good_rf = zeros(length(anat_labels),1);

for i=1:5
    perc_good_rf(i) = length(find(all_tc_struct(i).pvals_rf_in_out<0.05/length(all_groups)))/length(all_tc_struct(i).pvals_rf_in_out);
end

figure; bar(perc_good_rf);
set(gca,'xtick',[1:5],'xticklabels',anat_labels);

single_multi_peak = zeros(5,2);
figure;
for i=1:5
    single_multi_peak(i,1) = length(find(all_tc_struct(i).peaks(find(all_tc_struct(i).pvals_rf_in_out<0.05/length(all_groups)))==1))/length(all_tc_struct(i).pvals_rf_in_out);
    single_multi_peak(i,2) = length(find(all_tc_struct(i).peaks(find(all_tc_struct(i).pvals_rf_in_out<0.05/length(all_groups)))>=2))/length(all_tc_struct(i).pvals_rf_in_out);
    %single_multi_peak(i,3) = length(find(all_tc_struct(i).peaks(find(all_tc_struct(i).pvals_rf_in_out<0.05/length(all_groups)))>=3))/length(all_tc_struct(i).pvals_rf_in_out);
end

bh=bar(single_multi_peak,'stacked');
colorSet = [];
for i = 1:2
    myColors = [pt_color; hg_color; pp_color; pstg_color; mstg_color];
    bh(i).FaceColor = 'flat';
    bh(i).CData = myColors;
    bh(i).EdgeColor = 'none';
    if i==1
        bh(i).FaceAlpha = 1;
    else
        bh(i).FaceAlpha = 0.5;
    end
end
set(gca,'xtick',[1:5],'xticklabels',anat_labels);
set(gca,'ytick',[0:0.25:1],'yticklabels',[0:25:100]);
legend('single peak','multi-peak');
ylabel('% significant RF sites');
%print_quality_fig(gcf,'Figures/Fig2G.eps',10,4,4,'inches','epsc');

for i=1:5
    fprintf(1,'%s, total: %.2f, single peak: %.2f, multi: %.2f percent\n', ...
        anat_labels{i}, sum(single_multi_peak(i,:))*100, ...
        100*single_multi_peak(i,1)/sum(single_multi_peak(i,:)),...
        100*single_multi_peak(i,2)/sum(single_multi_peak(i,:)));
end

%%
peak_counts = [];
anat_mat = [];
%[all_tc_struct(4).peaks; all_tc_struct(4).peaks;
for i=1:5
    all_peaks = all_tc_struct(i).peaks;
    non_sig_peaks = find(all_tc_struct(i).pvals_rf_in_out>=0.05/length(all_groups));
    all_peaks(non_sig_peaks) = 0;
    anat_mat = [anat_mat; repmat(i, length(all_peaks), 1)];
    peak_counts = [peak_counts; all_peaks];
end
peak_counts(peak_counts>2) = 2;

[tbl, chi2,p] = crosstab(anat_mat, peak_counts);
df = (size(tbl,1)-1)*(size(tbl,2)-1);

fprintf(1,'Proportion of single vs multi-peaked RFs differed significantly:\n');
fprintf(1,'chi2 = %.1f, df=%d, n=%d electrodes, p=%2.2g\n', chi2, df, sum(tbl(:)), p);
% Do the same thing but count single- and multi-peaked together
%peak_counts(peak_counts>1) = 1;
%[tbl, chi2,p] = crosstab(anat_mat, peak_counts)

%% Figure 2H

all_rf_strf_corrs = [];
all_rf_strf_pvals = [];
grp = [];
for g=1:5
    [rho,p]=corr(all_tc_struct(g).tc, all_tc_struct(g).tc_speech, 'tail','right');
    rf_strf_corrs = diag(rho);
    rf_strf_pvals = diag(p);
    %rf_strf_corrs = rf_strf_corrs(all_tc_struct(g).vcorrs>0.1);
    all_rf_strf_corrs = [all_rf_strf_corrs; rf_strf_corrs];
    all_rf_strf_pvals = [all_rf_strf_pvals; rf_strf_pvals];
    grp = [grp; repmat(g, length(rf_strf_corrs),1)];
end

good_elecs = union(find(all_pvals_tones<(0.05/length(all_pvals_tones))),find(all_pvals_speech<0.05/length(all_pvals_speech)));
fh=figure;
boxplot(all_rf_strf_corrs(good_elecs), grp(good_elecs),'colors', [pt_color;hg_color; pp_color; ...
    pstg_color; mstg_color],'symbol','');

h = findobj(gca,'Tag','Box');
cs = [pt_color; hg_color; pp_color; ...
    pstg_color; mstg_color];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),cs(length(h)-(j-1),:),'FaceAlpha',0.5);
end

tags = {'Lower Whisker','Upper Whisker','Lower Adjacent Value',...
    'Upper Adjacent Value','Median','Box'};
for tt=1:length(tags)
    
    h = findobj(gca,'Tag',tags{tt});
    for j=1:length(h)
        set(h(j), 'linewidth',1,'color','k','linestyle','-');
    end
end
hold all;
plot([0.5 5.5],[0 0], 'k--');
set(gca,'xtick',[1:5],'xticklabel',{'PT','HG','PP','pSTG','mSTG'});
ylabel('r(speech, pure tone)');

[P,ANOVATAB,STATS] = kruskalwallis(all_rf_strf_corrs(good_elecs), grp(good_elecs))
[c,m,h,nms] = multcompare(STATS);

ht=1.2;
fprintf(1,'********************\n');
fprintf(1,'Correlation between STRF and pure tone RF\n');
for i=1:size(c,1)
    fprintf(1,'%s vs %s: %.2f, [%.2f, %.2f], %.4f\n', anat_labels{c(i,1)}, anat_labels{c(i,2)}, c(i,4) ,c(i,3), c(i,5), c(i,6));
    if c(i,6)<0.05
        ast = '*';
        if c(i,6) < 0.01
            ast='**';
            if c(i,6) < 0.001
                ast='***';
            end
        end
        figure(fh);
        plot([c(i,1) c(i,2)], [ht ht],'k-');
        plot([c(i,1) c(i,1)], [ht-0.02 ht], 'k-');
        plot([c(i,2) c(i,2)], [ht-0.02 ht], 'k-');
        text(mean([c(i,2),c(i,1)]), ht+0.02, ast,'fontsize',12);
        ht =ht+0.2;
    end
end
axis tight;
%print_quality_fig(gcf,'Figure2H.eps',8,3,3,'inches','epsc');
%%
pct_gt_zero = zeros(5,1);
for g=1:5
    pct_gt_zero(g) = length(find(all_rf_strf_pvals(good_elecs)<0.05 & grp(good_elecs)==g))/length(find(grp(good_elecs)==g));
    %length(find(all_rf_strf_corrs(good_elecs)>0 & grp(good_elecs)==g))/length(find(grp(good_elecs)==g));
    fprintf(1,'Percent greater than zero: %.2f\n', pct_gt_zero(g)*100)
end