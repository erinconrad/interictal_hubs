function dist_cosi_pc(whichPts,saved_out,out)

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
main_spike_results = [results_folder,'main_spikes/'];
if ~exist(main_spike_results,'dir')
    mkdir(main_spike_results);
end

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
spike_folder = [results_folder,'new_spikes/'];

if isempty(whichPts)
    whichPts = [20 103 106 107 35 109 110 111 94 97];
end

if saved_out == 1
    
    %out = load([main_spike_results,'out.mat']);
    %out = out.out;
    
else
    out = initialize_out_struct(length(whichPts));
    
    %% Get spike details
    fprintf('Getting spike details for pt...\n');
    for i = 1:length(whichPts)
        p = whichPts(i);
        fprintf('%d of %d\n',i,length(whichPts));
        out(i) = get_gdf_details(p);
    end
    save([main_spike_results,'out'],'out');
end

%% Colors
cols = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250];

%% Initialize figure
figure
set(gcf,'position',[29 100 1340 550])
tiledlayout(2,3,'TileSpacing','tight','padding','tight')

%% First row is example plots for a single patient
p = 1;
dist = out(p).dist;
cosi = out(p).cosi;
pc = out(p).metrics.added_pc;
labels = out(p).unchanged_labels;

% Correlate distance and co-spike index
nexttile
plot(dist,cosi,'o','linewidth',2,'color',cols(1,:))
if 0
plot(dist,cosi,'o','color',[1 1 1])
text(dist,cosi,labels,'horizontalalignment','center');
end
xlabel('Distance from added electrodes (mm)')
ylabel('Co-spike index')
set(gca,'fontsize',15)


% Correlate distance and pearson correlation with added electrodes
nexttile
plot(dist,pc,'o','linewidth',2,'color',cols(2,:))
xlabel('Distance from added electrodes (mm)')
ylabel('Functional connectivity')
set(gca,'fontsize',15)


% Correlate pearson correlation and co-spike index
nexttile
plot(pc,cosi,'o','linewidth',2,'color',cols(3,:))
xlabel('Functional connectivity')
ylabel('Co-spike index')
set(gca,'fontsize',15)

%% Second row is aggregate patient data
all_dist_cosi = [];
all_dist_pc = [];
all_cosi_pc = [];
for i = 1:length(whichPts)
    dist = out(i).dist;
    cosi = out(i).cosi;
    if isempty(out(i).metrics)
        pc = nan(length(dist),1);
    else
        pc = out(i).metrics.added_pc;
    end
    
    dist_cosi_r = corr(dist,cosi,'rows','pairwise');
    dist_pc_r = corr(dist,pc,'rows','pairwise');
    cosi_pc_r = corr(cosi,pc,'rows','pairwise');
    
    % Fisher r-to-z transformation, get z-score out
    [z,score,p] = fisher_transform(dist_cosi_r,sum(~isnan(dist) & ~isnan(cosi)));
    all_dist_cosi = [all_dist_cosi;dist_cosi_r z score p];
    
    [z,score,p] = fisher_transform(dist_pc_r,sum(~isnan(dist) & ~isnan(pc)));
    all_dist_pc = [all_dist_pc;dist_pc_r z score p];
    
    [z,score,p] = fisher_transform(cosi_pc_r,sum(~isnan(cosi) & ~isnan(pc)));
    all_cosi_pc = [all_cosi_pc;cosi_pc_r z score p];
    
end

%{
% Stouffer's method to combine the z scores
k = length(all_dist_cosi);
dist_cosi_z = nansum(all_dist_cosi(:,3))/sqrt(k);
dist_pc_z = nansum(all_dist_pc(:,3))/sqrt(k);
cosi_pc_z = nansum(all_cosi_pc(:,3))/sqrt(k);

% Do significance testing on the combined z scores
dist_cosi_p = 2*normcdf(-abs(dist_cosi_z));
dist_pc_p = 2*normcdf(-abs(dist_pc_z));
cosi_pc_p = 2*normcdf(-abs(cosi_pc_z));

% Back transform to get an average correlation coefficient
dist_cosi_r = tanh(nanmean(all_dist_cosi(:,2)));
dist_pc_r = tanh(nanmean(all_dist_pc(:,2)));
cosi_pc_r = tanh(nanmean(all_cosi_pc(:,2)));
%}

% Do an unpaired two-tailed t-test on the z's
[~,dist_cosi_p,~,stats] = ttest(all_dist_cosi(:,2));
dist_cosi_df = stats.df;
dist_cosi_tstat = stats.tstat;

[~,dist_pc_p,~,stats] = ttest(all_dist_pc(:,2));
dist_pc_df = stats.df;
dist_pc_tstat = stats.tstat;

[~,cosi_pc_p,~,stats] = ttest(all_cosi_pc(:,2));
cosi_pc_df = stats.df;
cosi_pc_tstat = stats.tstat;

% Plot dist_cosi for each patient
nexttile
plot(all_dist_cosi(:,1),'o','linewidth',2,'markersize',15,'color',cols(1,:))
hold on
%plot(xlim,[dist_cosi_r dist_cosi_r],'linewidth',2);
set(gca,'fontsize',15)
ylim([-1 1])
xlim([0.5 length(whichPts)+0.5])
plot(xlim,[0 0],'k--')
xticklabels([])
xlabel('Patient')
ylabel('Distance-CSI correlation')
xl = xlim;
yl = ylim;
text(xl(2),yl(2),sprintf('Mean r = %1.2f\n%s',mean(all_dist_cosi(:,1)),get_p_text(dist_cosi_p)),...
    'fontsize',15,'horizontalalignment','right','verticalalignment','top')


% Plot dist_pc for each patient
nexttile
plot(all_dist_pc(:,1),'o','linewidth',2,'markersize',15,'color',cols(2,:))
hold on
%plot(xlim,[dist_pc_r dist_pc_r],'linewidth',2);
set(gca,'fontsize',15)
ylim([-1 1])
xlim([0.5 length(whichPts)+0.5])

plot(xlim,[0 0],'k--')
xticklabels([])
xlabel('Patient')
ylabel('Distance-FC correlation')
xl = xlim;
yl = ylim;
text(xl(2),yl(2),sprintf('Mean r = %1.2f\n%s',mean(all_dist_pc(:,1)),get_p_text(dist_pc_p)),...
    'fontsize',15,'horizontalalignment','right','verticalalignment','top')

% Plot pc_cosi for each patient
nexttile
plot(all_cosi_pc(:,1),'o','linewidth',2,'markersize',15,'color',cols(3,:))
hold on
%plot(xlim,[cosi_pc_r cosi_pc_r],'linewidth',2);
set(gca,'fontsize',15)
xlim([0.5 length(whichPts)+0.5])

ylim([-1 1])
plot(xlim,[0 0],'k--')
xticklabels([])
xlabel('Patient')
ylabel('FC-CSI correlation')
xl = xlim;
yl = ylim;
text(xl(2),yl(2),sprintf('Mean r = %1.2f\n%s',mean(all_cosi_pc(:,1)),get_p_text(cosi_pc_p)),...
    'fontsize',15,'horizontalalignment','right','verticalalignment','top')

annotation('textbox',[0 0.91 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0.33 0.91 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0.67 0.91 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.42 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
annotation('textbox',[0.33 0.42 0.1 0.1],'String','E','fontsize',20,'linestyle','none')
annotation('textbox',[0.67 0.42 0.1 0.1],'String','F','fontsize',20,'linestyle','none')

print(gcf,[main_spike_results,'dist_cosi_pc_corr'],'-dpng');

%% Text
fprintf(['\nAggregated across patients, there was a significant negative '...
    'correlation between distance and co-spike index (mean r = %1.2f, '...
    'two-tailed unpaired t-test of individual patient Fisher''s r-to-z '...
    'transformed correlation coefficients: t(%d) = %1.2f, %s). There was also a '...
    'significant negative correlation between distance and functional '...
    'connectivity (mean r = %1.2f, t(%d) = %1.2f, %s), and a significant '...
    'positive correlation between co-spike index and functional connectivity '...
    '(mean r = %1.2f, t(%d) = %1.2f, %s) (Figure 1).\n'],...
    mean(all_dist_cosi(:,1)),dist_cosi_df,dist_cosi_tstat,get_p_text(dist_cosi_p),...
    mean(all_dist_pc(:,1)),dist_pc_df,dist_pc_tstat,get_p_text(dist_pc_p),...
    mean(all_cosi_pc(:,1)),cosi_pc_df,cosi_pc_tstat,get_p_text(cosi_pc_p));

end