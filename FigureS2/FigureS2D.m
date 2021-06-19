%calc_info_MID_vs_STA()
% Need to replace the STA here with the ridge regression STRF we already
% fit because this is the comparison we want to do

addpath(genpath('/Users/jsh3653/Documents/UCSF/code/hamiltonoganianchang2020/'));
addpath(genpath('/Users/jsh3653/Dropbox/Heschls_STRFs/code_release/Liberty_MID_code'));
config_paths;
heschl_load_data;
%%
useAnaName = 'CustomAna7areas';

%%
subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09'};
mid_dir = '~/Dropbox/Heschls_STRFs/data/MID/';
figure(2);
all_info = [];
for s=1:length(subjects)
    subj = subjects{s};
    date = 20111004;
    site = 9;
    Nparts = 4;
    opt_level1 = 1;
    opt_level2 = 2;
    rw = 80;
    col = 60;
    
    % Get the true MID electrode numbers
    mid_feat_file = sprintf('%s/%s_TrainTest_MID.mat', mid_dir, subj);
    stimresp_dat = load(mid_feat_file);
    mid_elecs = stimresp_dat.elec_nums;
    features = stimresp_dat.features;
    test_stim = stimresp_dat.test_stim(features.melspec,:)'; % Get only spectrogram
%     test_resp = stimresp_dat.test_resp(cell_num,:)';
%     test_resp = test_resp./std(test_resp);
%     
%     maxmax = max(abs(test_resp));
%    xedges = linspace(-maxmax, maxmax, Nbins);
    
    % Create delayed stimulus matrix
    fprintf(1,'Creating delay matrix\n');
    [nt, ndim] = size(test_stim);
    dStims = [];
    delays = 0:59;
    ndelays = length(delays);
    for di = 1:ndelays
        d = delays(di);
        dstim = zeros(nt, ndim);
        if d<0
            dstim(1:d,:) = test_stim(end-d:end,:);
        elseif d>0
            dstim(d+1:end, :) = test_stim(1:end-d,:);
        else
            dstim = test_stim;
        end
        dStims = [dStims dstim];
    end
    
    anat = imgNative.(subj).(useAnaName);
    
    mid_files_dir = sprintf('%s/MID/All_MID1_MID2_output_files_for_%s', paper_data_dir, subj);
    cd(mid_files_dir);
    nsites = length(dir('rpx1pxpxt_sta*1.dat'))
    %nsites = 58;
    info_matrix = zeros(nsites,5); % ridge STRF, MID1, MID2, MID12, anatomy_number
    for c = 1:nsites
        [mtx_sta, mtx1, mtx2, stainf, mid1inf, mid2inf, mid12inf, strfinf] = get_MIDs_Heschls_LH(date, site, c, Nparts,rw,col,opt_level1,opt_level2,subj,paper_data_dir,dStims);
        %plot(mean(stainf), mean(mid12inf),'.'); hold all;
        size(mid1inf)
        info_matrix(c, 1) = mean(strfinf);
        info_matrix(c, 2) = (mid1inf(1)+sign(mid1inf(1)*mid1inf(2))*mid1inf(2)+sign(mid1inf(1)*mid1inf(3))*mid1inf(3)+sign(mid1inf(1)*mid1inf(4))*mid1inf(4))/Nparts;%mean(mid1inf);
        info_matrix(c, 3) = (mid2inf(1)+sign(mid2inf(1)*mid2inf(2))*mid2inf(2)+sign(mid2inf(1)*mid2inf(3))*mid2inf(3)+sign(mid2inf(1)*mid2inf(4))*mid2inf(4))/Nparts;%mean(mid2inf);
        info_matrix(c, 4) = (mid12inf(1)+sign(mid12inf(1)*mid12inf(2))*mid12inf(2)+sign(mid12inf(1)*mid12inf(3))*mid12inf(3)+sign(mid12inf(1)*mid12inf(4))*mid12inf(4))/Nparts;%mean(mid12inf);
        %info_matrix(c, 5) = strfinf;
        info_matrix(c, 5) = anat(mid_elecs(c));
    end
    all_info = [all_info; info_matrix];
    %
    
    labels = {'Ridge STRF','MID1','MID2','MID12'};

    pp=1;
    for i=1:4
        for j=1:4
            if i<j
                
                subplot(4,4,pp);
                for a = 1:7
                    chs = find(info_matrix(:,5) == a);
                    plot(info_matrix(chs,i), info_matrix(chs,j), '.', ...
                        'color', area7Cols(a,:), 'markersize',15);
                    hold all;
                end
                xlabel(sprintf('%s info',labels{i}));
                ylabel(sprintf('%s info',labels{j}));
                hold all;
                plot([0, 0.1], [0, 0.1], 'k');
            end
            pp=pp+1;
        end
    end
    %keyboard;
    %pause;
end

print_quality_fig(gcf,'/Users/jsh3653/Dropbox/Heschls_STRFs/data/figures/MID_vs_STRF.eps',10,8,8,'inches','epsc');

%% Average percent improvement
fprintf(1,'********\n');

for i=2:4
    p = mean((all_info(:,i)-all_info(:,1))./all_info(:,1));
    fprintf(1,'Percent improvement from STRF to %s: %.2f\n', labels{i}, p*100);
    p = mean(all_info(:,i)-all_info(:,1));
    fprintf(1,'Mean info improvement from STRF to %s: %.4f\n', labels{i}, p);
end

fprintf(1,'********\n');
for i=1:4
    p = mean((all_info(:,i)-all_info(:,2))./all_info(:,2));
    fprintf(1,'Percent improvement from MID1 to %s: %.2f\n', labels{i}, p*100);
    p = mean(all_info(:,i)-all_info(:,2));
    fprintf(1,'Mean info improvement from MID1 to %s: %.4f\n', labels{i}, p);
end
fprintf(1,'********\n');
