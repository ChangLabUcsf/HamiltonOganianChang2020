function [X, xy, img, chs, varargout] = load_data_clustering_onerep(subj, names, chs_to_use, highpass_filt, do_norm, logflag, sort_by, do_mean)
% function [X, xy, img, chs, varargout] = load_data_clustering_onerep(subj, names, chs_to_use, highpass_filt, do_norm, logflag, sort_by, do_mean)
%
% Load the data to be used in NMF clustering 
% 
% Inputs:
%   subj: the subject ID ('S01')
%   names: the names of the stimuli (from out.name)
%   chs_to_use: which channels to include in the final output matrix 
%   highpass_filt: whether to highpass filter the high gamma data or not
%   do_norm: Norm the response data
%   logflag: Use log-transformed high gamma
%   sort_by: 'none' or 'sent_length' to sort by sentence length
%   do_mean: whether to take the mean across repetitions of the same
%            stimulus
%
%

if nargin<2
    names = [];
end
if nargin<3
    chs_to_use=[];
end
if nargin < 4 || isempty(highpass_filt)
    highpass_filt = 0;
end
if nargin < 5
    do_norm =0;
end
if nargin <6
    logflag = 1;
end
if nargin < 7
    sort_by = 'sent_length';
end
if nargin<8
    do_mean = 0;
end

rm_first_phoneme = 0;

% output arguments
nout = max(nargout,1)-4;
if nout>=1
    get_sent_matrix = 1;
else
    get_sent_matrix = 0;
end


fprintf(1,'Loading data for subject %s\n', subj);
datadir = '/Users/liberty/Documents/UCSF/data/timit';
%datadir = '/Volumes/SeagateBackup/TIMITcontrols/TIMIT';
subdir = sprintf('%s/%s',datadir,subj);

if logflag
    outfiles = dir(sprintf('%s/%s_*out_resp_log.mat',subdir, subj));
else
    outfiles = dir(sprintf('%s/%s_*out_resp.mat',subdir, subj));
end
fprintf(1,'Loading outfile: %s\n', outfiles(1).name);
tic;
grid256=load(sprintf('%s/%s', subdir, outfiles(1).name));

alldat = grid256;
toc;

if ~isempty(names)
    fprintf(1,'Only using sentences in all subjects\n');
    [~,j]=intersect({alldat.out.name},names);
    alldat.out = alldat.out(j);
end

fprintf(1,'Getting bad channels\n');
badChans = getBadChannels(subdir);
badChans = [];
nchans = size(alldat.out(1).resp,1);
fprintf(1,'%d total channels, ', nchans);
chs = setdiff(1:nchans,badChans);
if ~isempty(chs_to_use)
    chs = intersect(chs, chs_to_use);
end
fprintf(1,'%d good channels\n', length(chs));

sentence_lengths = zeros(length(alldat.out),1);
for sentence = 1:length(alldat.out)
    sentence_lengths(sentence) = size(alldat.out(sentence).resp,2)-20;
end
max_sent_len = max(sentence_lengths);

switch sort_by
    case 'sent_length'
        fprintf(1,'Sorting by sentence length\n');
        [~, sentences_order] = sort(sentence_lengths);
    case 'none'
        sentences_order = 1:length(names); % Keep in the original order
end



if get_sent_matrix
    sent_matrix = nan(length(chs), length(alldat.out), max_sent_len); % 288 channels x reps x time
end
cnt=1;

tic;
X = [];
stim = [];
alldat.out = add_NF_to_out(alldat.out);

for i = 1:length(alldat.out)
    if do_mean
        r = mean(alldat.out(sentences_order(i)).resp(chs,21:end,:),3);
    else
        r = alldat.out(sentences_order(i)).resp(chs,21:end,1);
    end
    if highpass_filt>0
        fs = 100;
        nyq = fs/2; %Nyquist frequency
        cof = highpass_filt; %cut-off frequency of 0.5 Hz
        fprintf(1,'High pass filtering at %d Hz\n', cof);
        [b,a] = butter(3, cof/nyq, 'low'); %this designs a 3-pole low-pass filter
        r = (r' - filtfilt(b,a,r'))'; %filtfilt makes it non-causal (fwd/backward)
    end
    if get_sent_matrix
        sent_matrix(:, i, 1:size(r,2)) = r;
    end
    alldat.out(sentences_order(i)).resp = r;
    alldat.out(sentences_order(i)).phnmat = alldat.out(sentences_order(i)).phnmat(:,21:end);
    if rm_first_phoneme
        fprintf(1,'Marking the first phoneme\n');
        start_first = find(sum(alldat.out(sentences_order(i)).phnmat)==1,1);
        this_phn = find(alldat.out(sentences_order(i)).phnmat(:,start_first)==1)
        end_first = find(alldat.out(sentences_order(i)).phnmat(this_phn,start_first:end)==0,1)
        alldat.out(sentences_order(i)).phnmat(this_phn,start_first:end_first+start_first+1)=99;
    end
    alldat.out(sentences_order(i)).aud = alldat.out(sentences_order(i)).aud(:,21:end);
    alldat.out(sentences_order(i)).NFs = alldat.out(sentences_order(i)).NFs(:,21:end);
    
end
X = [alldat.out.resp]';

if do_norm
    fprintf(1,'Normalizing data\n');
    X = bsxfun(@minus, X, mean(X));
    normX = sqrt(sum(X.^2));
    X = bsxfun(@times, X, 1./normX);
end

stim = [alldat.out.phnmat]';
stimNF = [alldat.out.NFs]';

toc;

xy=[];
img=[];

if get_sent_matrix
    fprintf(1,'Getting sentence matrix\n');
    varargout(1) = {sent_matrix};
    varargout(2) = {stim};
    varargout(3) = {stimNF};
end