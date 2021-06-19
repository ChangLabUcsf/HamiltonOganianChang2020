% Figure 2 subpanels D and E
% Plot the lead versus lag time of the latencies for different areas

heschl_load_data;
%heschl_load_anatomy;
useAnaName = 'CustomAna7area';
%%
rmin = outall.S01(1).dataf; % one second

include7AreasName = {'PT','pmHG','alHG','PP','pSTG onset','pSTG other','mSTG'};

maxlags = [];
minlags = [];
corramp = [];
for i = 1:length(subjects)
    subj = subjects{i}
    out = outall.(subj);
    
    elecs = allr.elnum(allr.sid==i & AnalysisEl);
    roi_areas = allr.(useAnaName)(allr.sid==i & AnalysisEl);
    
    for r=1:length(out)
        first_phn_time = get_first_phoneme_timit(out(r).name) + out(r).befaft(1);
        first_phoneme = round(first_phn_time * out(r).dataf);
        out(r).resp = mean(out(r).resp(:,first_phoneme:end,:),3); % start at first phoneme
    end
    
    all_resp = [out.resp];
    
    for area1 = 1:7
        elecs1 = elecs(roi_areas==area1);
        for area2 = 1:7
            elecs2 = elecs(roi_areas==area2);
            if ~isempty(elecs1) && ~isempty(elecs2)
                if area1 < area2
                    fprintf(1,'subj %s: %s to %s\n', subjects{i}, include7AreasName{area1}, include7AreasName{area2});
                    for e1=1:length(elecs1)
                        for e2=1:length(elecs2)
                            resp1 = all_resp(elecs1(e1),:);
                            resp2 = all_resp(elecs2(e2),:);
                            [acor, lag] = xcorr(resp1,resp2, 100,'coeff');
                            
                            if max(acor)>0.1
                                maxlags = [maxlags; lag(argmax(acor)) max(acor) area1 area2 i];
                            end
                            if min(acor)<-0.1
                                [mincor, mincorlag] = min(acor);
                                minlags = [minlags; lag(mincorlag) mincor area1 area2 i];
                            end
                        end
                    end
                end
            end
        end
    end
end

%%
% add the group type
grpnum = 1;
maxlags = [maxlags zeros(size(maxlags,1),1)];
grpname = {};
for grp1=1:7
    for grp2=1:7
        if grp1<grp2
            maxlags(maxlags(:,3)==grp1 & maxlags(:,4)==grp2,5)=grpnum;
            grpnum = grpnum+1;
            grpname = [grpname; sprintf('%s-%s', include7AreasName{grp1}, include7AreasName{grp2})];
        end
    end
end

lag_times = lag/100.0;

for g=1:length(unique(grpname))
    if mean(maxlags(maxlags(:,5)==g,1))>0
        disp(grpname{g})
    end
end

new_maxlags = maxlags;

% Plot by median lead/lag
pair_means = grpstats(maxlags(:,1), maxlags(:,5), @mean);
flip_pos = find(pair_means>0);
new_maxlags(ismember(maxlags(:,5), flip_pos),1) = -maxlags(ismember(maxlags(:,5), flip_pos),1);

[sortedlag,sortedOrder] = sort(grpstats(new_maxlags(:,1),new_maxlags(:,5),@mean));

pos = 1:length(unique(grpname));
pos(sortedOrder) = 1:length(unique(grpname));
%[sortedlagabs, sortedOrderabs] = sort(abs(sortedlag));
%pos(sortedOrder(sortedOrderabs)) = 1:length(unique(grpname));

figure;
boxplot(new_maxlags(:,1), new_maxlags(:,5), 'labels', grpname, 'plotstyle','compact',...
    'symbol','','orientation','horizontal',...
    'color',[0.5 0.5 0.5], 'positions',pos);
%
xlim([-20, 20]); hold all;
plot([0, 0], [1, length(unique(grpname))],'k');

split_labels = {};
left_label = {}; right_label = {};
for i=1:length(grpname)
    split_labels{i} = split(grpname{i},'-');
    left_label{i} = split_labels{i}{1};
    right_label{i} = split_labels{i}{2};
end

for i=1:length(flip_pos)
    tmp = left_label{flip_pos(i)};
    left_label{flip_pos(i)} = right_label{flip_pos(i)};
    right_label{flip_pos(i)} = tmp;
end


left_label = left_label(sortedOrder);
right_label = right_label(sortedOrder);

set(gca,'ytick',1:length(unique(grpname)),'yticklabel',left_label);
yyaxis right;
cla;
set(gca,'ytick',1:length(unique(grpname)),'yticklabel',right_label);
set(gca,'xtick',[-20:5:10], 'xticklabel', -0.2:0.05:0.2);
set(gca,'xlim',[-21 10]);
set(gca,'ylim',[0.5 length(unique(grpname))+0.5])
xlabel({'Time lag (s)','Leading - Lagging'})

print_quality_fig(gcf,'figures/Figure2D.eps',8,4,2.6,'inches','epsc');

% Find the non-significant ones (that overlap zero)
for i=1:length(grpname)
    [h,p]=ttest(new_maxlags(new_maxlags(:,5)==i,1));
    if p>0.05
        fprintf(1,'OVERLAPS ZERO %s p=%.4f\n', grpname{i}, p);
    else
        fprintf(1,'%s p=%.4f\n', grpname{i}, p);
    end
end

% OVERLAPS ZERO PT-pmHG p=0.6974
% OVERLAPS ZERO pmHG-pSTG onset p=0.1501
% OVERLAPS ZERO alHG-pSTG other p=0.8582

%% Connectivity diagram

lmax = round(max(abs(sortedlag)));

cmap = flipud(cbrewer('div','RdBu',lmax*2+1));
%cmap = flipud(cbrewer('seq','YlOrRd',lmax+1));
area_coords = [3 1;
    2 1.2;
    1 1;
    0 1;
    2.5 0;
    1.4 0;
    0.3 0;]
figure;
hold on;


lflag=0;
for grp1=1:7
    for grp2=1:7
        if grp1<grp2
            %fprintf(1,' %s to %s\n', include7AreasName{grp1}, include7AreasName{grp2});
            
            % Get line width based on the magnitude of this lagged corr
            acor_mean = mean(maxlags(maxlags(:,3)==grp1 & maxlags(:,4)==grp2,2));
            
            % Get the line color based on the max lag (leading vs lagging)
            leadlag = round(mean(maxlags(maxlags(:,3)==grp1 & maxlags(:,4)==grp2,1)));
            
            % Test whether the lags are significantly different from zero
            [h,p] = ttest(maxlags(maxlags(:,3)==grp1 & maxlags(:,4)==grp2,1));
            
            if leadlag > 0 % the second thing is leading instead of the first
                % flip the arrow around
                lflag = 1;
                grp2_tmp = grp2;
                grp1_tmp = grp1;
                grp2 = grp1_tmp;
                grp1 = grp2_tmp;
                leadlag = -leadlag;
            end
            
            if p<0.05
                fprintf(1,' SIGNIF. %s to %s\n', include7AreasName{grp1}, include7AreasName{grp2});
                %                 dp = area_coords(grp2,:) - area_coords(grp1,:);
                %                 quiver(area_coords(grp1,1), area_coords(grp1,2), ...
                %                     dp(1), dp(2), ...
                %                     'color', cmap(leadlag+lmax+1,:), ...
                %                     'linewidth', acor_mean*20);
                
                if (grp1==1 && grp2==3)
                    bez = bezier.eval([area_coords(1,:); area_coords(2,1) area_coords(2,2)+0.5; area_coords(3,:)])
                    plot(bez(:,1), bez(:,2), ...
                        'color', cmap(leadlag+lmax+1,:), ...
                        'linewidth', acor_mean*20);
                elseif (grp1==1 && grp2==4)
                    bez = bezier.eval([area_coords(1,:); area_coords(2,1)-0.5 area_coords(2,2)+0.5; area_coords(4,:)])
                    plot(bez(:,1), bez(:,2), ...
                        'color', cmap(leadlag+lmax+1,:), ...
                        'linewidth', acor_mean*20);
                elseif (grp1==5 && grp2==7)
                    bez = bezier.eval([area_coords(5,:); area_coords(6,1) area_coords(6,2)-0.5; area_coords(7,:)])
                    plot(bez(:,1), bez(:,2), ...
                        'color', cmap(leadlag+lmax+1,:), ...
                        'linewidth', acor_mean*20);
                else
                    plot([area_coords(grp1,1) area_coords(grp2,1)], ...
                        [area_coords(grp1,2) area_coords(grp2,2)], ...
                        'color', cmap(leadlag+lmax+1,:), ...
                        'linewidth', acor_mean*20);
                end
                if lflag
                    grp1 = grp1_tmp;
                    grp2 = grp2_tmp;
                    lflag=0;
                end
            else
                fprintf(1,'p=%.2f, NOT SIGNIF. %s to %s\n', p, include7AreasName{grp1}, include7AreasName{grp2});
            end
        end
    end
end
for i=1:7
    scatter(area_coords(i,1), area_coords(i,2), 100, 'k','filled'); hold all;
    text(area_coords(i,1)+0.1, area_coords(i,2)-0.1, include7AreasName{i});
end
axis([-0.5 3.5 -0.5 1.5]);
axis off;
colormap(cmap);
cb = colorbar;
set(cb, 'ytick', linspace(0, 0.5, 7), 'yticklabel',linspace(-0.06, 0, 7),'ylim',[0 0.5])
print_quality_fig(gcf,'figures/Figure2E.eps',8,4,3,'inches','epsc');

