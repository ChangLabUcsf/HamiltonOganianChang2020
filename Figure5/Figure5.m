% Figure 5
%
% Liberty Hamilton 2020
%
% Figure 5A Percent Variance Explained and basis functions

if 0
    addpath('../util1');
    heschl_load_data;
end

%%
% Load the NMF data
load('/Users/liberty/Dropbox/Heschls_STRFs/data/NMF/allsubj_clustdata_avgtrial.mat');
load('/Users/liberty/Dropbox/Heschls_STRFs/data/NMF/allsubj_NMF_avgtrial.mat');

%%
bandnum=4; % This corresponds to 4 clusters
sb_order = [1 4 2 3];

nF = FLs{bandnum};
nF = nF./max(sum(nF));

chgroup=[];
for s=1:length(subjects)
    chgroup = [chgroup; repmat(s, length(chs_by_subj{s}),1)];
end

pp=1;
for s=1:length(subjects)
    for ch=1:length(chs_by_subj{s})
        ch_sub(pp,2) = chs_by_subj{s}(ch);
        pp=pp+1;
    end
end

% Get the color map
cmap  = cbrewer('qual','Dark2',9);

% Figure 5A - Factor basis functions

figure(20);
legend_labels = {};
for sb=1:size(nF,1)
    plot(factor_TS{bandnum}(:,sb), 'color', cmap(sb,:),'linewidth',3); hold on;
end
xlabel('Time (s)');
set(gca,'xtick',[0 30:50:250], 'xticklabel', [-0.3 0:0.5:2.5]);
ylabel('Weighted Z-scored HG');
hold all;
plot([30 30],[-15 25],'color',[0.5 0.5 0.5]);
fontsize=8; wd=3; ht=3; unit_type='inches'; figtype='epsc';

%print_quality_fig(gcf, '~/Dropbox/Heschls_STRFs/data/figures/Figure5A.eps', fontsize, wd, ht, unit_type, figtype);

%%
% Anatomical labels to use
all_anat = {'transversetemporal','planumtemporale','planumpolare','superiortemporal',...
    'ctx_lh_S_temporal_sup','ctx_rh_S_temporal_sup','ctx_lh_G_temp_sup-Plan_tempo',...
    'ctx_rh_G_temp_sup-Plan_tempo','ctx_rh_S_temporal_transverse','ctx_lh_S_temporal_transverse',...
    'ctx_rh_G_temp_sup-Plan_polar','ctx_lh_G_temp_sup-Plan_polar','ctx_lh_G_temp_sup-Lateral',...
    'ctx_rh_G_temp_sup-Lateral','pSTG','mSTG'};

% Get the sentence data and do the NMF reconstruction
x=squeeze(nanmean(sent_matrix_all,2));
X = x(:,1:300)';

% Calculate percent variance explained using NMF reconstruction
pve=NMF_percentvariance(X, FLs, W);
Ms = 1:length(FLs);
Ms(find(pve==0)) = [];
pve(pve==0)=[];

figure;
plot(Ms, pve, 'linewidth',2);
hold all;
plot([4 4], [0.84 0.96], 'k--');
xlabel('Number of bases');
ylabel('Percent variance');
set(gca,'xtick',[1 4 8 12 16]);
axis([1 16 0.84 0.96]);
%print_quality_fig(gcf, '~/Dropbox/Heschls_STRFs/data/figures/Figure5A_subpanel.eps', fontsize, wd, ht, unit_type, figtype);

%%
temporal=mniMesh.temporal.lh;
%load('~/Dropbox/Heschls_STRFs/data/anatomy/cvs_avg35_inMNI152/Meshes/cvs_avg35_inMNI152_lh_temporal_pial.mat');
figure(2);
c_h = [];
g = struct(); % Always use the same camera angle so
% everything can be overlaid perfectly
for tF = 1:size(nF,1)
    %subplot(2,ceil(size(nF,1)/2),tF);
    subplot(4,1,tF);
    c_h(tF) = ctmr_gauss_plot(temporal, [0 0 0], 0, 'lh',1); alpha(1); hold all;
    loc_view(-106, 42);
    
    g(tF).CameraPosition     = get(gca, 'CameraPosition');
    g(tF).CameraTarget       = get(gca, 'CameraTarget');
    g(tF).CameraUpVector     = get(gca, 'CameraUpVector');
    g(tF).CameraUpVectorMode = get(gca, 'CameraUpVectorMode');
    g(tF).CameraViewAngle    = get(gca, 'CameraViewAngle');
end
e=1;
h=[];
elecs = {};
for sb=1:size(nF,1)
    elecs{sb} = [];
end
anat_labels = {};
for a=1:size(nF,1)
    anat_labels{a} = {};
end
for s=1:length(subjects)
    subj = subjects{s};
    subnum = find(strcmp(subj, subjects));
    elecmatrix = imgmni.(subj).elecmatrix;
    anatomy = imgNative.(subj).newAnatomy;
    anatomy(isnan(elecmatrix(:,1)),4) = {'NaN'};
    
    elecmatrix(elecmatrix(:,1)>0,:) = NaN;
    
    thisNF = nF(:,chgroup==subnum);
    typ_chans = zeros(size(thisNF,2),1);
    for c=1:length(chs_by_subj{subnum})
        ch = chs_by_subj{subnum}(c);
        
        best_sb = argmax(thisNF(:,c));
        
        x = thisNF(best_sb,c);
        
        xsort = sort(thisNF(:,c));
        
        alab=anatomy(ch,4)
        if xsort(end)>2*xsort(end-1) && ismember(alab{1}, all_anat) && xsort(end) > 0.1
            ecolor2 = [cmap(best_sb,:)];
            %subplot(2,ceil(size(nF,1)/2),best_sb);
            subplot(4,1,sb_order(best_sb));
            
            elecmatrix(ch,:)
            h(e) = el_add(elecmatrix(ch,:),'color',ecolor2, 'msize',6);
            elecs{best_sb} = [elecs{best_sb}; s ch x];
            
            anat_labels{best_sb} = [anat_labels{best_sb}; alab];
            e=e+1;
            
        end
    end
end

%save('/Users/liberty/Dropbox/Heschls_STRFs/data/NMF/elecs_bestsubband_4bases.mat','elecs','subjects','chs_by_subj','nF');

%%
% Figure 5B
% Plots bar plots with frequency
all_anat2 = {'planumtemporale','transversetemporal','planumpolare','pSTG','aSTG'};

heschls_labels = {'transversetemporal','ctx_rh_S_temporal_transverse',...
    'ctx_rh_G_temporal_transverse', 'ctx_lh_S_temporal_transverse',...
    'ctx_lh_G_temporal_transverse'};

stg_labels = {'superiortemporal','ctx_rh_G_temp_sup-Lateral','ctx_lh_G_temp_sup-Lateral','pSTG','mSTG'};
pstg_labels = {'pSTG'};
astg_labels = {'mSTG'};
pp_labels = {'planumpolare','ctx_rh_G_temp_sup-Plan_polar','ctx_lh_G_temp_sup-Plan_polar'};
pt_labels = {'planumtemporale','ctx_rh_G_temp_sup-Plan_tempo','ctx_lh_G_temp_sup-Plan_tempo'};

all_anat_bar = {pt_labels; heschls_labels; pp_labels; pstg_labels; astg_labels};

hists = zeros(length(all_anat_bar),size(nF,1));
for i=1:size(nF,1)
    for a=1:length(all_anat_bar)
        for b=1:length(all_anat_bar{a})
            hists(a,i) = hists(a,i) + length(find(strcmp(all_anat_bar{a}{b},anat_labels{i})));
        end
    end
end

pt_color = [0.17, 0.22, 0.58];
pp_color = [0.62, 0.22, 0.58];
stg_color = [0.69, 0.70, 0.21];
pstg_color = [0.69, 0.12, 0.14];
hg_color = [0.06, 0.69, 0.29];

figure;
sb_order = [1 3 4 2];
hists = hists(:,sb_order);

hist_vals = hists'./repmat(sum(hists'),size(nF,1),1);
for i=1:size(hist_vals,1)
    xcoords = linspace(0.7,1.3,5)+(i-1)
    alphas = [1, 0.8, 0.6, 0.4, 0.2];
    for j=1:5
        bb=bar(xcoords(j), hist_vals(i,j), 0.135);
        bb.FaceColor = cmap(sb_order(i),:);
        bb.FaceAlpha = alphas(j);
        hold all;
    end
    
    %pause();
end
set(gca,'xtick',[linspace(0.7,1.3,5) linspace(0.7,1.3,5)+1 linspace(0.7,1.3,5)+2 linspace(0.7,1.3,5)+3])
set(gca,'xticklabel',{'PT','HG','PP','pSTG','mSTG'})
set(gca,'XTickLabelRotation',90);
axis([0.5 4.5 0 1]);
xlabel('Basis function');
ylabel('Frequency');
fontsize=8; wd=3; ht=3; unit_type='inches'; figtype='epsc';
%print_quality_fig(gcf, '~/Dropbox/Heschls_STRFs/data/figures/Figure5B.eps', fontsize, wd, ht, unit_type, figtype);

%%
% Figure 5D
% Individual sentence responses S03

out = outall.S03;

[i,j]=get_timit5_inds(out);
out=out(i);

sent_nums = [2 5 7];
out=out(sent_nums);

mresp = [];
sresp = [];
sent_bounds = 1;
aud = [];
for i=1:length(out)
    aud = [aud out(i).aud];
    mresp = [mresp mean(out(i).resp,3)];
    sresp = [sresp std(out(i).resp,[],3)/sqrt(size(out(i).resp,3))];
    sent_onset = find(mean(out(i).phnmat,1)>0, 1);
    sent_offset = find(mean(out(i).phnmat,1)>0, 1, 'last');
    wav_offset = size(out(i).phnmatonset,2);
    
    sent_bounds = [sent_bounds sent_bounds(end)+sent_onset sent_bounds(end)+sent_offset sent_bounds(end)+wav_offset];
    
end

f=figure;
subplot(5,1,1);
imagesc(aud); axis xy; hold on;
set(gca,'xtick',[0:100:800],'xticklabel',[0:1:8]);
colormap(flipud(gray));

subplot(5,1,2);
plotstd(mresp(266,:)', sresp(266,:)', cmap(1,:));

subplot(5,1,5);
plotstd(mresp(38,:)', sresp(38,:)', cmap(2,:));

subplot(5,1,3);
plotstd(mresp(54,:)', sresp(54,:)', cmap(3,:));

subplot(5,1,4);
plotstd(mresp(267,:)', sresp(267,:)', cmap(4,:));

for sb=1:5
    subplot(5,1,sb);
    if sb==1
        yl = [1 80];
    else
        yl = [-2 3];
    end
    for b=1:length(sent_bounds)
        plot([sent_bounds(b) sent_bounds(b)], [yl(1) yl(2)], 'color', [0.5 0.5 0.5]);
        
    end
    set(gca,'xtick',[0:100:800],'xticklabel',[0:1:8]);
    if sb>1
        set(gca,'ytick',[-2 0 3]);
        ylabel('Z-scored HG');
    else
        set(gca,'ytick',[12 44 80], 'yticklabel', [0.5 2 8]);
        ylabel('Freq (kHz)');
    end
    axis([1 size(mresp,2) yl(1) yl(2)]);
    if sb==5
        xlabel('Time (s)');
    end
end

%print_quality_fig(gcf,sprintf('/Users/liberty/Dropbox/Heschls_STRFs/data/figures/Figure5D.eps'),8,3,6,'inches','epsc');