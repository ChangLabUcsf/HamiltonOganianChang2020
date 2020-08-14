function [FLs, W, factor_TS] = run_NMF_avgtrial_heschls(subjects, spatial_weighting, logflag, highpass_filt)
%function [FLs, factor_TS] = run_group_NMF_heschls(subj, hemshort, spatial_weighting)
%This makes oblique spatial factors for ECoG data as demo,
%using the distance-based weighting of spatial_dwFA.m.
%v.3 uses actim_contours.m for plotting
%

%close all;

% Set some parameters

use_sent_matrix = 1;
if nargin < 2
    spatial_weighting = 0;
end
if spatial_weighting
    wtname = '';
else
    wtname = '_nospatialweight';
end
if nargin < 3
    logflag = 1;
end
if logflag == 1
    logname = '_log';
    fprintf(1,'Using log high gamma\n');
else
    logname = '';
end
if nargin <4
    highpass_filt = 0;
end
if highpass_filt~=0
    highpass_name = sprintf('_sublowpass%1.1fHz', highpass_filt);
else
    highpass_name = '';
end
sort_by = 'none';
sort_by = 'sent_length';
do_plots=0;
do_norm = 1;
do_mean = 1;
if do_mean
    meannm = '_domean';
else
    meannm = '';
end

% Locate and load the data
data_dir = '~/Documents/UCSF/Heschls/data/';
timit_dir = '~/Documents/UCSF/data/timit/';

fig_dir = sprintf('%s/figures', data_dir);

disp('Load ECoG Dat...');

fprintf(1,'Multi-subject analysis\n');
out_file = sprintf('%s/allsubj_NMF_avgtrial_heschls_2019_consolidated.mat', data_dir);
out_file2 = sprintf('%s/allsubj_clustdata_avgtrial_heschls_2019_consolidated.mat', data_dir);

clustdata_file = out_file2

if 1%~exist(out_file, 'file')
    
    % Get sentences that were played to all subjects
    chs = [];
    get_out = 0;
    %name_file = '/Users/liberty/Documents/UCSF/Heschls/data/allsubj_out_overlap_heschls_allstg.mat';
    %name_file = '/Users/liberty/Documents/UCSF/Heschls/data/allsubj_out_overlap_heschlsonly.mat';
    name_file = '/Users/liberty/Documents/UCSF/Heschls/data/allsubj_out_overlap_heschlsonly_2019.mat';
    
    if length(subjects)>1
        if 1%~exist(name_file)
            [names] = get_overlapping_sentences(subjects, get_out, chs, logflag);
            save(name_file, 'names');
            names
        else
            fprintf(1,'Loading file with overlapping sentences\n');
            load(name_file);
        end
    else
        names = []
    end
    if 0%exist(clustdata_file)
        fprintf(1,'Loading cluster data\n');
        load(clustdata_file);
    else
        Xall = [];
        sent_matrix_all = [];
        chs_subj = zeros(length(subjects),1);
        chs_by_subj = {};
        for s=1:length(subjects)
            subj = subjects{s}
            
            if strcmp(subj, 'EC53') || strcmp(subj, 'EC89') || strcmp(subj, 'EC85') || strcmp(subj, 'EC143') 
                nchans = 288;
            elseif strcmp(subj, 'EC75') || strcmp(subj, 'EC124')
                nchans = 64;
            elseif strcmp(subj, 'EC157')
                nchans = 302;
            elseif strcmp(subj, 'EC193')
                nchans = 274;
            elseif strcmp(subj, 'EC180')
                nchans = 160; 
            else
                nchans = 256;
            end
            stg_elects = 1:nchans;
            %disp(stg_elects);
            
            [X, ~, ~, chs, sent_matrix, stim] = load_data_clustering_onerep(subj, names, stg_elects, highpass_filt, do_norm, logflag, sort_by, do_mean);
            
            if 1
                fprintf(1,'Getting electrodes with STRF correlations > 0.1\n');
                
                strf_dir = '/Users/liberty/Documents/UCSF/Heschls/data/strfs';
                strf_dir2 = sprintf('%s/%s/strfs/ecog%d_log', timit_dir, subj, nchans);
                
                %strf_file = sprintf('%s/%s_STRF_0to600ms_100Hz_allchans_zscore_intercept_noedge_samealpha_optimizeR2.hf5', strf_dir, subj)
                strf_file = sprintf('%s/%s_STRF_0to600ms_100Hz_allchans_zscore_intercept_noedge_samealpha_consolidatedcode.hf5', strf_dir, subj);
                vcorrs = h5read(strf_file, '/vcorrs');
                
                vcorrs = vcorrs(chs); % Get rid of the bad channels
                
                sig_inds = find(vcorrs>0.1);
            end
            
            
            chs = chs(sig_inds);
            
            chs_subj(s) = length(chs);
            chs_by_subj{s} = chs;
            
            length(sig_inds)
            
            X = X(:,sig_inds);
            sent_matrix = sent_matrix(sig_inds,:,:);
            size(Xall)
            size(X)
            Xall = [Xall X];
            sent_matrix_all = cat(1,sent_matrix_all, sent_matrix);
            
        end
        fprintf(1,'Saving the sentence matrix data... ');
        save(out_file2, 'sent_matrix_all', 'Xall', 'chs_subj', 'chs_by_subj', 'subjects', '-v7.3');
        fprintf(1,'Done.\n');
        
    end
    
    fprintf(1,'Calculating XX\n');
    tic;
    x = squeeze(nanmean(sent_matrix_all,2));
    x = x(:,1:300); % Cut off end because the trials are not this long and there are edge artifacts
    XX = x*x';
    toc;
    
    Q = length(size(XX,1));
    if ~spatial_weighting
        DM = ones(Q,Q);
    end
    
    disp('Calculating FLs...');
    tic;
    
    Ms = [1:16];%16];
    FLs = {};
    %factor_TS = {};
    W = {};
    fprintf(1,'Band: ');
    for sb = 1:length(Ms)
        fprintf(1,'%d ', Ms(sb));
        %subplot(1,1,sb);
        
        M = Ms(sb);
        %if length(chs) >= M
        try
            % Initialize to varimax solution
            W0 = spatial_dwFA(x',[],M,'w',DM);
            W0 = max(W0,0); %rectify
            
            G0 = []; % compute G0 in NMF_convex
            Nit = 5000;
            
            % Run convex NMF
            [W{sb}, FLs{sb}] = NMF_convex(XX, M, W0, G0, Nit);
            %FLs{sb} = FLs{sb}';
            
            factor_TS{sb} = x'*W{sb}; % If we wanted the whole time series
            % Get time series sorted by sentence
            %factor_TS{sb}=reshape(reshape(sent_matrix_all,size(sent_matrix_all,1),size(sent_matrix_all,2)*size(sent_matrix_all,3))'*W{sb},size(sent_matrix_all,2),size(sent_matrix_all,3),size(W{sb},2));
            
            %nrows = floor(sqrt(M));
            %ncols = ceil(M/nrows);
            %fprintf(1,'Spatially reordering factor loadings\n');
            %FLs{sb} = reorg_FLs_chs(FLs{sb}, chs, nrows, ncols);
        catch
            fprintf(1,'Did not converge\n');
            FLs{sb} = [];
            W{sb} = [];
            factor_TS{sb} = [];
        end
        %else
        %    fprintf(1,'\nThere are only %d channels, cannot cluster into %d levels\n', length(chs), M);
        %end
        
    end
    fprintf(1,'\n');
    toc;
    
    % Order factor loading clusters by x coordinate for plotting -DOESNT
    % WORK
    %fprintf(1,'Spatially reordering factor loadings\n');
    %[FLs, W] = order_FLs_NMF(FLs, W, xy(:,chs), hem);
    
    fprintf(1,'Saving factor loadings in %s\n', out_file);
    save(out_file, 'FLs', 'Ms', 'DM', 'chs_by_subj', 'W', 'subjects', 'factor_TS'); %, 'transient', 'sustained');
else
    fprintf(1,'Loading FL file %s\n', out_file);
    load(out_file);
    load(out_file2);
    x = squeeze(nanmean(sent_matrix_all,2));
    XX = x*x';
    
    for sb=1:length(Ms)
        factor_TS{sb}=x'*W{sb};
    end
end
%keyboard;

if do_plots
    figure;
    imagesc(XX);
    cs=cumsum(chs_subj);
    hold on;
    for i=1:length(chs_subj)
        plot([cs(i) cs(i)],[1 size(XX,1)],'k')
        plot([1 size(XX,1)],[cs(i) cs(i)],'k')
    end
    caxis([-0.4 0.4])
    set(gca,'xtick',cs,'ytick',cs)
    xlabel('Electrode')
    ylabel('Electrode');
    set(gca,'tickdir','out');
    colorbar
    axis square
    set(gcf,'PaperSize',[6 6],'PaperUnits','inches','PaperPosition',[0 0 6 6],'renderer','painters');
    %saveas(gcf,'groupXX.eps','epsc');
    %print_quality_fig(gcf, '~/Documents/UCSF/CorticalOrganization/data/groupXX_STGMTGHG_domean.eps', 20, 10, 10, 'inches');
    
    
    offsets = zeros(size(factor_TS{1},1),1);
    for i=1:size(factor_TS{1},1)
        if ~isempty(find(isnan(factor_TS{1}(i,:,2)),1))
            offsets(i)=find(isnan(factor_TS{1}(i,:,2)),1)-50;
        else
            offsets(i)=size(factor_TS{1},2)-50;
        end
    end
    
    figure;
    factor_TS{1}(isnan(factor_TS{1}))=0;
    for i=1:2
        subplot(1,2,i);
        imagesc(factor_TS{1}(:,:,i));
        hold on;
        plot([31 31],[0 size(factor_TS{1},1)],'k');
        %plot(offsets,1:length(offsets),'k');
        %set(gca,'xtick',[31:50:offsets(end)],'xticklabel',[0:0.5:0.01*(offsets(end)-30)]);
        xlabel('Time (s)');
        ylabel('Ordered sentence #');
        caxis([-50 50]);
    end
    
    if 1
        cs=[0; cs];
        for i=1:length(cs)-1
            XX_within(cs(i)+1:cs(i+1),cs(i)+1:cs(i+1))=1;
        end
        
        figure;
        plot([-0.04:0.02:1], histc(XX(find(XX_within==0)),[-0.04:0.02:1])/length(XX(find(XX_within==0))),'b')
        hold all; plot([-0.04:0.02:1], histc(XX(find(XX_within==1)),[-0.04:0.02:1])/length(XX(find(XX_within==1))),'r')
        xlabel('Pearson R');
        ylabel('Percentage');
        axis([-0.04 1 0 0.15]);
        set(gca,'xtick',[0:0.25:1],'ytick',[0:0.05:0.15],'yticklabel',[0:5:15]);
        legend('across subject','within subject');
        print_quality_fig(gcf, '~/Documents/UCSF/CorticalOrganization/data/groupXX_withinacross_all_domean_heschls_avgtrial.eps', 20, 4, 4, 'inches');
    end
    %keyboard;
end
