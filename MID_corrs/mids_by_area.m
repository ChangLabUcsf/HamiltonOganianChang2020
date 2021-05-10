% Plot MIDs by area

heschl_load_data
config_paths;


%%
mid2_subjects = {'S01','S02','S03','S05','S06','S07'};

mid_dir = '~/Dropbox/Heschls_STRFs/data/MID/';
strf_dir = sprintf('%s/STRFs/spect_zscore', paper_data_dir);

useAnaName = 'CustomAna7areas';
include7AreasName = {'PT','pmHG','alHG','PP','pSTG onset','pSTG other','mSTG'};

pp = ones(7,1);
ff = [1:7];
%%
mid1_by_area = {};
mid2_by_area = {};
strf_by_area = {};
for i=1:7
    mid1_by_area{i} = [];
    mid2_by_area{i} = [];
    strf_by_area{i} = [];
end
for s=1:length(mid2_subjects)
    subj = mid2_subjects{s};
    mid_file = sprintf('%s/%s_MID12.mat', mid_dir, subj);
    mid_data = load(mid_file);
    
    % Get the true MID electrode numbers
    mid_feat_file = sprintf('%s/%s_TrainTest_MID.mat', mid_dir, subj);
    mid_elecs = load(mid_feat_file, 'elec_nums');
    
    strf_file = sprintf('%s/%s_STRF_spect_zscore.hf5', strf_dir, subj);
    %vcorrs = h5read(strf_file, '/vcorrs');
    wts = h5read(strf_file, '/wts');
    
    anat = imgNative.(subj).(useAnaName);
    for ch=1:length(mid_data.data)
        a_num = anat(mid_elecs.elec_nums(ch));
        
        mid1_by_area{a_num} = [mid1_by_area{a_num} mid_data.data(ch).MID1{1}(:)];
        mid2_by_area{a_num} = [mid2_by_area{a_num} mid_data.data(ch).MID2{1}(:)];
        strf_by_area{a_num} = [strf_by_area{a_num} wts(mid_elecs.elec_nums(ch),2:end)'];
        if 0
            figure(ff(a_num));
            subplot(6,3,pp(a_num));
            imagesc(mid_data.data(ch).MID1{1}); axis xy;
            title(sprintf('MID1 ch%d %s', mid_elecs.elec_nums(ch), subj));
            pp(a_num) = pp(a_num) + 1;

            subplot(6,3,pp(a_num));
            imagesc(mid_data.data(ch).MID2{1}); axis xy;
            title(sprintf('MID2 ch%d %s', mid_elecs.elec_nums(ch), subj));
            pp(a_num) = pp(a_num) + 1;

            subplot(6,3,pp(a_num));
            strf = wts(mid_elecs.elec_nums(ch),2:end);
            strf = reshape(strf, [], 60);
            imagesc(strf); axis xy;
            title(sprintf('STRF ch%d %s %s', mid_elecs.elec_nums(ch), subj, include7AreasName{a_num}));
            pp(a_num) = pp(a_num) + 1;

            if pp(a_num) > 18
                ff(a_num) = ff(a_num) + 7;
                pp(a_num) = 1;
            end
        end
    end
end

%%
pp=1;
figure;
for i=1:7
    subplot(7,3,pp);
    strf = reshape(mean(strf_by_area{i},2),80,60);
    mmax = max(abs(strf(:)));
    imagesc(fliplr(strf), [-mmax mmax]); axis xy;
    set(gca,'xtick',[1 30 60], 'xticklabel', [-0.6 -0.3 0]);
    set(gca,'ytick',[1 44 80], 'yticklabel', [0.75 2 8]);
    set(gca,'yaxislocation','right');
    title(sprintf('%s STRF',include7AreasName{i}));
    pp=pp+1;
    
    subplot(7,3,pp);
    mid1 = reshape(mean(mid1_by_area{i},2),80,60);
    mmax = max(abs(mid1(:)));
    imagesc(fliplr(mid1), [-mmax mmax]); axis xy;
    set(gca,'xtick',[1 30 60], 'xticklabel', [-0.6 -0.3 0]);
    set(gca,'ytick',[1 44 80], 'yticklabel', [0.75 2 8]);
    set(gca,'yaxislocation','right');
    title('MID1');
    pp=pp+1;
    
    subplot(7,3,pp);
    mid2 = reshape(mean(mid2_by_area{i},2),80,60);
    mmax = max(abs(mid2(:)));
    imagesc(fliplr(mid2), [-mmax mmax]); axis xy;
    set(gca,'xtick',[1 30 60], 'xticklabel', [-0.6 -0.3 0]);
    set(gca,'ytick',[1 44 80], 'yticklabel', [0.75 2 8]);
    set(gca,'yaxislocation','right');
    title('MID2');
    pp=pp+1;
    
end
colormap(flipud(cbrewer('div','PRGn',256)));
for pp=19:21
subplot(7,3,pp);
xlabel('Time (s)');
end
ylabel('Freq (kHz)');

%% Plot MID1 only vs STRF for pmHG because they seem less correlated with the tones

mid2_subjects = {'S01','S02','S03','S04','S05','S06','S07','S08','S09'};
useAnaName = 'CustomAna7areas';
pp = ones(7,1);
ff = [1:7];
for s=1:length(mid2_subjects)
    subj = mid2_subjects{s};
    mid_file = sprintf('%s/%s_MID1.mat', mid_dir, subj);
    mid_data = load(mid_file);
    
    % Get the true MID electrode numbers
    mid_feat_file = sprintf('%s/%s_TrainTest_MID.mat', mid_dir, subj);
    mid_elecs = load(mid_feat_file, 'elec_nums');
    
    strf_file = sprintf('%s/%s_STRF_spect_zscore.hf5', strf_dir, subj);
    %vcorrs = h5read(strf_file, '/vcorrs');
    wts = h5read(strf_file, '/wts');
    
    anat = imgNative.(subj).(useAnaName);
    for ch=1:length(mid_data.data)
        a_num = anat(mid_elecs.elec_nums(ch));
        
        if a_num==2
            figure(ff(a_num));
            subplot(6,6,pp(a_num));
            mid=mid_data.data(ch).MID1{1};
            imagesc(mid); axis xy;
            title(sprintf('MID1 ch%d %s', mid_elecs.elec_nums(ch), subj));
            pp(a_num) = pp(a_num) + 1;
            
            subplot(6,6,pp(a_num));
            strf = wts(mid_elecs.elec_nums(ch),2:end);
            strf = reshape(strf, [], 60);
            imagesc(strf); axis xy;
            title(sprintf('STRF ch%d %s %s', mid_elecs.elec_nums(ch), subj, include7AreasName{a_num}));
            pp(a_num) = pp(a_num) + 1;
            
            if pp(a_num) > 36
                ff(a_num) = ff(a_num) + 7;
                pp(a_num) = 1;
            end
            
            figure(999);
            plot(sum(strf,2)/max(sum(strf,2)), 'r--'); hold all;
            plot(sum(mid,2)/max(sum(mid,2)), 'b');
            
        end
    end
end