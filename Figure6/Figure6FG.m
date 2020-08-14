% Figure 6FG
%
% Hamilton, Oganian, and Chang 2020.
%
% Load data and get code dependencies

addpath(genpath('../util1'));
heschl_load_data;

%% ---------- Figure 6: word number analysis, cluster 4
elecs_by_cluster = load('elecs_by_cluster.mat');

%% word duration 
wrddur=[sentdet.sentdet.wordDur];
wrdnum = [sentdet.sentdet.wordOnsOff];
a2 = wrddur(wrddur>0);
a2num = wrdnum(1,wrddur>0);
figure
scatter(a2, a2num+rand(size(a2num))/5);
xlabel('word length');
ylabel('word number');
figure, histogram(a2(a2num==1));
%% get mean HG for each word
for csid = 1:length(SID)
    cs = SID{csid};
    for csent = 1:length(outall.(cs))        
        csentid = find(strcmpi(outall.(cs)(csent).name, {sentdet.sentdet.name}));
        outall.(cs)(csent).sentId = csentid;
        outall.(cs)(csent).meanWrdResp = (mean(outall.(cs)(csent).resp,3)*sentdet.sentdet(csentid).wordfeat')./sum(sentdet.sentdet(csentid).wordfeat,2)';
    end
end
%% across all electrodes by NMF cluster
wrdResp = cell(4,1);
for ccl = 1:4
    cels = elecs_by_cluster.elecs{ccl};
    allsubj = elecs_by_cluster.subjects(cels(:,1));
    wrdResp{ccl}=nan(length(cels),12,499);
    for cel = 1:length(cels)
        cs = allsubj{cel};
        for csent = 1:length(outall.(cs))
            cresps = outall.(cs)(csent).meanWrdResp(cels(cel,2),:);
            wrdResp{ccl}(cel,1:size(cresps,2),outall.(cs)(csent).sentId) = cresps;
        end
    end
end

%% ----------  F+G  grand mean by cluster
clustercols = [27 158 119 ; 117 112 179 ; 231 41 138;217 95 2]/255;
clear cpl;
figure,
subplot(1,2,1)
cla;hold on;
for ccl = 1:4   
    cpl(ccl) = shadedErrorBar(1:12, nanmean(nanmean(wrdResp{ccl}, 3), 1),  nansem(nanmean(wrdResp{ccl}, 3), 1),{'color', clustercols(ccl,:)});
    xlim([1 8])
    xlabel('Word #')
    ylabel('Average high gamma')    
end
legend([cpl.mainLine],'c1', 'c2', 'c3', 'c4');
% stats for each cluster
for ccl = 1:4   
    meanwrdResp{ccl} = nanmean(wrdResp{ccl},3);
    for cel = 1:size(meanwrdResp{ccl},1)
        modb{ccl}(cel) = regress((1:8)',zscore(meanwrdResp{ccl}(cel,1:8)'));
    end
end

% plot regression betas
% figure
% cla;hold on;
% for i = 1:4
%     scatter(i*ones(size(modb{i})),modb{i});
% end

allbetas = cell2mat(modb);
allbetacluster =[];for i = 1:4, allbetacluster =[allbetacluster ; i*ones(length(modb{i}),1)];end

% figure
subplot(1,2,2)
boxplot(allbetas',allbetacluster);
horzline(0);
xticklabels({'1','2','3','4'});
ylabel('mean word number beta');
xlabel('Cluster #');
title('electrode betas: word number ~ mean HGA')

%% ----------  stats for F+G
%% 1-way anova on betas
[p,tbl] = anovan(allbetas', allbetacluster);

group = allbetacluster==4;
anovan(allbetas', {group});

group = allbetacluster;%
anovan(allbetas(ismember(allbetacluster, [1 4]))', {group(ismember(allbetacluster, [1 4]))});

for i = 1:4
%     [h(i), p(i),~,stats(i)]=...
        ttest(allbetas(allbetacluster == i),0);
end


%% --- mixed model with word # and word durations 
%% word responses by word duration
wdurs=[];
for i = 1:length(sentdet.sentdet)
    wdurs(1:3,1:size(sentdet.sentdet(i).wordOnsTimes,2),i) = sentdet.sentdet(i).wordOnsTimes;
end
wdurs(wdurs==0)=nan;
%% histogram of word durations
figure,
histogram(diff(wdurs(1:2,:,:),1,1))

%% mixed model: 
% hga ~ word number x time since sentence onset
lmWnum = cell(4,1);
for cl = 1:4
    disp(cl);
    cnel = size(wrdResp{cl},1);
    tblsentnum = repmat(1:499,cnel,12,1);
    tblwnum = repmat(reshape(wdurs(3,:,:), 1,[]),cnel,1);
    tblwonstime = repmat(reshape(wdurs(1,:,:),1,[]),cnel,1);
    tblresp =squeeze(reshape(wrdResp{cl},cnel,1,[]));
    tblelnum = squeeze(reshape(repmat((1:cnel)',1,12,499), cnel, 1,[]));
        
    respdata=table(reshape(tblwonstime,[],1),reshape(tblwnum,[],1), reshape(tblresp,[],1),reshape(tblelnum,[],1),reshape(tblsentnum,[],1),...
        'VariableNames', {'onsTime', 'wordNum', 'hga', 'elnum', 'sentnum'});
    respdata(respdata.wordNum==0,:)=[];
%     grpstats(respdata.hga, respdata.wordNum)
    
    lmWnum{cl}=fitlme(respdata, 'hga~onsTime+wordNum +(1+onsTime+wordNum|elnum)+(1|sentnum)');
end

%% --- model on words 3 - 7 only
for ccl = 1:4   
    meanwrdResp{ccl} = nanmean(wrdResp{ccl},3);
    for cel = 1:size(meanwrdResp{ccl},1)
        modb2{ccl}(cel) = regress((3:7)',zscore(meanwrdResp{ccl}(cel,3:7)'));
    end
end
% plot regression betas

allbetas2 = cell2mat(modb2);
allbetacluster2 =[];for i = 1:4, allbetacluster2 =[allbetacluster2 ; i*ones(length(modb2{i}),1)];end

% figure
subplot(1,2,2)
boxplot(allbetas2',allbetacluster2);
horzline(0);
xticklabels({'1','2','3','4'});
xlabel('Cluster #');
ylabel('mean word number beta');
title('electrode betas: word number ~ mean HGA, words 3-7 only')
% 1-way anova on betas
[p,tbl] = anovan(allbetas2', allbetacluster2);

group = allbetacluster2==4;
anovan(allbetas2', {group});

group = allbetacluster2;%

anovan(allbetas2(ismember(allbetacluster2, [1 4]))', {group(ismember(allbetacluster2, [1 4]))});

for i = 1:4
    [h(i), p(i),~,stats(i)]=...
        ttest(allbetas2(allbetacluster2 == i),0);
end



%% --- more stuff
%% plot single electrodes in each cluster
figure,
for ccl = 1:4
    subplot(1,4,ccl)
    plot(nanmean(wrdResp{ccl}, 3)')
    xlim([1 10])
    xlabel('Word #')
    ylabel('Average high gamma')
    title(['cluster ' num2str(ccl)]);
end


