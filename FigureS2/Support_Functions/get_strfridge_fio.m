function [fio] = get_strfridge_fio(subj, cell_num, paper_data_dir, Nbins, dStims)
%function [mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] =
%get_dat_sta_fio(folder, location, cellnum, nv, nh, coeff, varargin)
%
%[mean_firing, x_mean, ior_mean, ior_std, px_mean, px_std] = ...
%     get_dat_sta_fio(folder, location, cellnum, coeff, Nparts, Nbins, Nbins_short, minpx)
%
% summary of changes made to accomodate ripple stimuli(May 31 2005):
% downsampling in ior2_mean is done by summing pxt within a larger bin, summing px within the 
% larger bin, and dividing the two numbers, also introduced dx, and probability distribution is now shown as a probability density
%
%
% The function assumes that you have files of the form:
%
%     rpx1pxpxt_sta_707_1_1x25x20_1_1.dat            
%     rpx1pxpxt_sta_707_1_1x25x20_1_2.dat            
%     rpx1pxpxt_sta_707_1_1x25x20_1_3.dat            
%     rpx1pxpxt_sta_707_1_1x25x20_1_4.dat     
%
% Parts of the file name are variable, such as 707, 25 , 20. These are specified by
% input arguments.
%
% folder : location of data. Example: 'C:\MATLABR2007b\work\data\MID_Inferior_Colliculus\';
%
% location : Tanya's site identifier. Examples: 707, 708, 709, etc.
%
% cellnum : Tanya's global cell identifier. Examples: 1, 2, 3, ...
%
% coeffs : coefficient vector which tells the polarity of the filters
%
% Npart = 4 : num of test reps in mid code. Should be 4.
%
% Nbins = 21 : num bins used to obtain projection values in mid code
%
% Nbins_short = 14 : num bins used to get nonlinearity values
%
%
% If you call get_data_sta_fio.m as
%
% [fio] = get_dat_sta_fio(folder, location, cellnum, nv, nh, coeff)
%
% then Nparts = 4, Nbins = 21, and Nbins_short = 14 by default.
%
% caa 2/17/09

% Parameters that Tatyana uses in her programs
% Nbins = 21; % used for input/output function but not RF
% Nbins_short = 14; % used for input/output function but not RF
% Nparts = 4; % used to identify files, not sure what the significance
%             % is, though



%     files{1} = rpx1pxpxt_sta_707_1_1x25x20_1_1.dat            
%     files{2} = rpx1pxpxt_sta_707_1_1x25x20_1_2.dat            
%     files{3} = rpx1pxpxt_sta_707_1_1x25x20_1_3.dat            
%     files{4} = rpx1pxpxt_sta_707_1_1x25x20_1_4.dat     

Nparts=4;
minpx = 0;
mid_dir = '~/Dropbox/Heschls_STRFs/data/MID/';
strf_dir = sprintf('%s/STRFs/spect_zscore', paper_data_dir);

% Get the true MID electrode numbers
mid_feat_file = sprintf('%s/%s_TrainTest_MID.mat', mid_dir, subj);
stimresp_dat = load(mid_feat_file);
mid_elecs = stimresp_dat.elec_nums;
test_resp = stimresp_dat.test_resp(cell_num,:)';
test_resp = test_resp./std(test_resp);

maxmax = max(abs(test_resp));
xedges = linspace(-maxmax, maxmax, Nbins);

fprintf(1,'Loading ridge STRF data\n');
strf_file = sprintf('%s/%s_STRF_spect_zscore.hf5', strf_dir, subj);
wts = h5read(strf_file, '/wts');
wts = wts(mid_elecs,:);

fprintf(1,'Calculating projection\n');
proj = dStims*wts(cell_num,2:end)'+wts(cell_num,1);

% scale to Z-score
proj = proj - mean(proj);
proj = proj./std(proj);

px = histcounts(proj, xedges,'Normalization','probability');
%ind0 = find( px < minpx );
%xm = sum( x.*px );
%x = x - xm; % center the projection values to 0
%rbar = mean(test_resp);

% Scale to Z-scores
test_resp = test_resp - mean(test_resp);
test_resp = test_resp./std(test_resp);

% This is the probability of the responses (P(spk), or really P(lfp))
% These should sum to 1.
pspk = histcounts(test_resp, xedges,'Normalization','probability');

% Edges for "response" vs. not (>0 or not)... maybe this is not the right
% way to do it for LFPs....
yedges = [xedges(1) 0 xedges(end)];
[N,c] = histcounts2( proj,test_resp, xedges, yedges, 'Normalization', 'probability');

%pxs_mean = sum(N,1);%'; % for sum over 2
pxs_mean = N(:,2)';
pxs_mean = pxs_mean./sum(pxs_mean);
ind1 = find(pxs_mean>0 & px>0);
information = sum( pxs_mean(ind1) .* log2( pxs_mean(ind1) ./ px(ind1) ) )/Nparts;

fio.info = information;


