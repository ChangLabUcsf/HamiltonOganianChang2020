config_paths;
heschl_load_data;
heschl_load_anatomy;
useAnaName = 'CustomAna7area';
%%
rmin = outall.S01(1).dataf; % one second

include7AreasName = {'PT','pmHG','alHG','PP','pSTG onset','pSTG other','mSTG'};

maxlags = [];
resp_by_area = {};
resp_by_area_sent = {};
for i=1:7
    resp_by_area{i} = [];
    resp_by_area_sent{i} = [];
end

mni_elecs = {};
for i=1:7
    mni_elecs{i} = [];
end
all_subj = [];
for i = 1:length(subjects)
    subj = subjects{i}
    out = outall.(subj);
    
    elecs = allr.elnum(allr.sid==i & AnalysisEl);
    roi_areas = allr.(useAnaName)(allr.sid==i & AnalysisEl);
    
    all_resp = zeros(size(out(1).resp,1), rmin, length(out));
    
    [timit5inds, timit5names] = get_timit5_inds(out);
    out2 = out(timit5inds);
    
    all_resp_by_sent = [mean(out2(1).resp,3) mean(out2(7).resp,3)];
    snd = [out2(1).sound; out2(7).sound];
    for r=1:length(out)
        first_phn_time = get_first_phoneme_timit(out(r).name) + out(r).befaft(1);
        first_phoneme = round(first_phn_time * out(r).dataf);
        all_resp(:,:,r) = mean(out(r).resp(:,first_phoneme:first_phoneme+rmin-1,:),3); % start at first phoneme    
    end
        
    for area1=1:7
        elecs1 = elecs(roi_areas==area1);
        all_subj = [all_subj; repmat(i,length(elecs1),1)];
        mni_elecs{area1} = [mni_elecs{area1}; imgmni.(subj).elecmatrix(elecs1,:)];
        resp_by_area{area1} = [resp_by_area{area1}; mean(all_resp(elecs1,:,:),3)];
        resp_by_area_sent{area1} = [resp_by_area_sent{area1}; all_resp_by_sent(elecs1,:,:)];
    end
end
fclose('all');

%%  Figure 2A
plotnum2 = [3 1 4 5 2 6 7]+1;

subplot(8,1,1);
plot(snd,'color','k'); axis off; axis tight;
for i=1:7
    subplot(8,1,plotnum2(i));
    imagesc(resp_by_area_sent{i}(lat_sort_by_area{i}(:,1),:), [0 4]);
    hold all;
    
    first_phn = get_first_phoneme_timit(out2(1).name) + out2(1).befaft(1);
    plot(first_phn*out2(1).dataf + lat_sort_by_area{i}(:,2),1:size(resp_by_area_sent{i},1),'r--');
    
    first_phn = get_first_phoneme_timit(out2(7).name) + out2(7).befaft(1);
    plot(first_phn*out2(7).dataf + size(out2(1).resp,2) + lat_sort_by_area{i}(:,2),1:size(resp_by_area_sent{i},1),'r--');
    colormap(flipud(gray));
    if plotnum2(i)==8
        set(gca,'xtick',[0:100:500],'xticklabel',[0:1:5]);
        xlabel('Time (s)');
    else
        set(gca,'xtick',[0:100:500],'xticklabel',[]);
        xlabel('');
    end
    set(gca,'ytick',[1 size(resp_by_area_sent{i},1)]);
    ylabel({sprintf(include7AreasName{i}), 'Electrode #'});
end
colorbar;

%print_quality_fig(gcf,'figures/latency_twosentences_customAna.eps',8,4,8,'inches','epsc');
%% Figure 2B latencies on brain
rmin=0.50*100;

latency_colors = cbrewer('div','PiYG',rmin); 


figure;
c_h = ctmr_gauss_plot(mniMesh.temporal.lh, [0 0 0], 0, 'lh');
alpha 0.5

lats = {};
all_lat_sent = [];
all_area = [];
for area=1:7
    derp = diff(resp_by_area{area},1,2);
    lat = argmax(derp');
    all_lat_sent = [all_lat_sent; lat'];
    all_area = [all_area; repmat(area,length(lat),1)];
    lat(lat>50) = 50;
    rh = find(mni_elecs{area}(:,1)>0);
    lat(rh) = [];
    lats{area} = lat;
    mni_elecs{area}(rh,:) = [];
    if area == 1 || area == 2 || area == 5
        el_add(mni_elecs{area}, 'color', latency_colors(lat,:), ...
               'edgecol','k','msize', 10);
    else
        el_add(mni_elecs{area}, 'color', latency_colors(lat,:), 'msize', 10);
    end
    hold all;
end
onset_elecs = mni_elecs{5};
onset_elecs(:,1) = onset_elecs(:,1) - 0.5;
el_add(onset_elecs, 'color', latency_colors(lats{5},:), 'edgecol','k','msize', 10);

set(gcf,'color','w');
loc_view(-93,60);

%print(fullfile(figDir,'brain_latencies.jpg'), '-djpeg', '-painters', '-r800')

%% Figure 2C
hgpm_color = [0.06, 0.69, 0.29];
hgal_color = [0.06, 0.8, 0.8];
pt_color = [0.17, 0.22, 0.58];
pp_color = [0.62, 0.22, 0.58];
mstg_color = [0.69, 0.70, 0.21];
pstg_color = [0.69, 0.12, 0.14];
pstg_onset_color = [0, 0, 0];
hg_color = [0.06, 0.69, 0.29];

[sortedlag,sortedOrder] = sort(grpstats(all_lat_sent, all_area, @median));
pos = 1:7;
pos(sortedOrder) = 1:7;

cmap = [pt_color; hgpm_color; hgal_color; pp_color; pstg_onset_color; pstg_color; mstg_color; ];
figure;
%subplot(1,2,1);
boxplot(all_lat_sent, all_area, 'orientation','horizontal','plotstyle','compact',...
        'positions',pos,'colors',cmap(sortedOrder,:)); hold all;
set(gca,'xtick',[1:10:100], 'xticklabel', [0:0.1:1]);
set(gca,'xlim',[0 51]);
set(gca,'ytick',1:size(cmap,1),'yticklabel',{include7AreasName{sortedOrder}});
xlabel('Latency (s)');
axis ij;
%print_quality_fig(gcf,'figures/latency_boxplot_0-5s_rev.eps',8,4,3,'inches','epsc');
%%
[p,anovatab,stats ] =kruskalwallis(all_lat_sent,{include7AreasName{all_area}});

for i=1:7
    for j=1:7
        if i<j
            pval=ranksum(all_lat_sent(all_area==i), all_lat_sent(all_area==j));
            fprintf(1,'%s vs %s, p=%2.2g', include7AreasName{i}, ...
                include7AreasName{j}, pval);
            if pval>0.05
                fprintf(1,'****NS*****\n');
            else
                fprintf(1,'\n');
            end
        end
    end
end