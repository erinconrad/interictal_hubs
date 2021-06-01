function spikiest_elec(rate,labels,change,surround,results_folder,name,...
    block_dur,run_dur,locs)

%% Remove EKG and scalp electrodes
ekg = identify_ekg_scalp(labels);
rate(ekg,:) = [];
labels(ekg) = [];
locs(ekg,:) = [];

%% Redefine rate according to run_dur
rate = rate/run_dur;

%% Remove electrodes with too many nans
too_many_nans = (sum(~isnan(rate(:,change-surround:change-1)),2) < 0.5*length(change-surround:change-1) | ...
    sum(~isnan(rate(:,change+1:change+surround)),2) < 0.5*length(change+1:change+surround));
rate(too_many_nans,:) = [];
labels(too_many_nans) = [];
locs(too_many_nans,:) =  [];

if isempty(rate)
    return
end
 
%% Identify pre times and post times
pre = 1:change-1;
post = change + 1:size(rate,2);
 
%% Get mean pre- rate for each channel and sort channels by this
mean_rate = nanmean(rate(:,change-surround:change-1),2);
[sorted_mean_rate,sorted_chs] = sort(mean_rate,'descend');

%% Re-sort labels and rate by their pre-rate
labels = labels(sorted_chs);
rate = rate(sorted_chs,:);

%% Get the spikiest electrode pre-revision
[~,spikiest_pre] = max(nanmean(rate(:,change-surround:change-1),2));
[~,spikiest_post] = max(nanmean(rate(:,change+1:change+surround),2));

%% For each time period, get the spikiest electrode and whether it agrees with the pre-revision
[~,spikiest] = max(rate,[],1);

% Make a matrix of 0s for all electrodes except 1 for spikiest
spikiest_matrix = zeros(size(rate));
for i = 1:size(rate,2)
    spikiest_matrix(spikiest(i),i) = 1;
end

%% Significance testing for agreement
pval_spikiest = compare_agree(rate,change,surround);

%% For each time period, get the correlation between its rate order and this mean pre-rate order
corr_time = nan(size(rate,2),1);
for i = 1:size(rate,2)
    corr_time(i) = corr(rate(:,i),sorted_mean_rate,'Type','Spearman',...
        'Rows','pairwise');
end

%% Significance testing for rate order stability
pval_rho = compare_rhos(rate,change,surround);

%{
%% For each time period, get the identity of the spikiest electrode
[~,sp_ch] = max(rate,[],1);

%% For each time period, list whether it agrees with this mode or not
sp_mode_agg = nan(size(rate,2),1);
for i = 1:size(rate,2)
    sp_mode_agg(i) = sp_ch(i) == spikiest_pre;
end
pre_sp_mode_agg = nanmean(sp_mode_agg(pre));
post_sp_mode_agg = nanmean(sp_mode_agg(post));
%}

%% Prep figure
figure
set(gcf,'position',[100 100 1300 600]) 
t = tiledlayout(2,2,'TileSpacing','compact');

%% Show raster
nexttile(1,[1,2])
times = 1:size(rate,2);
times = times*block_dur;
h = turn_nans_white(rate);
set(h,'XData',[0:size(rate,2)*block_dur])
xlim([0 size(rate,2)*block_dur])
hold on
cp = plot([change*block_dur change*block_dur],ylim,'r--','linewidth',3);
%yticklabels([])
xlabel('Hour')
ylabel('Electrode')
h = colorbar; 

ylabel(h,'Spikes/min','fontsize',20)
set(gca,'fontsize',20)

%% Mark the pre- and post-max
%{
if spikiest_pre == spikiest_post
    yticks(spikiest_pre);
    yticklabels('pre & post\rightarrow');
else
    yticks([spikiest_pre,spikiest_post])
    yticklabels({'pre\rightarrow','post\rightarrow'});
end
%}

%% Show correlation with pre-revision rate-order
nexttile(3)
plot(times,corr_time,'linewidth',2);
hold on
cp = plot([change*block_dur change*block_dur],ylim,'r--','linewidth',3);
title(sprintf('Rate order stability (p = %1.3f by permutation test)',pval_rho))
ylabel({'Spike rate consistency'});
xlabel('Hour')
set(gca,'fontsize',20)

%
%% Compare pre- and post-revision spikiest channel
%{
nexttile(4)
scatter3(locs(:,1),locs(:,2),locs(:,3),'k')
hold on
if spikiest_pre == spikiest_post
    scatter3(locs(spikiest_pre,1),locs(spikiest_pre,2),locs(spikiest_pre,3),...
        'b','filled');
else
    scatter3(locs(spikiest_pre,1),locs(spikiest_pre,2),locs(spikiest_pre,3),...
        'g','filled');
    scatter3(locs(spikiest_post,1),locs(spikiest_post,2),locs(spikiest_post,3),...
        'r','filled');
end
xticklabels([])
yticklabels([])
zticklabels([])
dist = vecnorm(locs(spikiest_pre,:)-locs(spikiest_post,:));
title(sprintf('Distance = %1.2f mm',dist))
%}
%
ax4 = nexttile(4);
imagesc(spikiest_matrix)
colormap(ax4,[1 1 1;0, 0.4470, 0.7410])
hold on
cp = plot([change*block_dur change*block_dur],ylim,'r--','linewidth',3);
title(sprintf('Spikiest electrode agreement (p = %1.3f by permutation test)',...
    pval_spikiest));
xlabel('Hour')
ylabel('Spikiest electrode')
set(gca,'fontsize',20)
%}

title(t,name,'fontsize',20)

outfolder = [results_folder,'consistency/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end
print(gcf,[outfolder,name],'-dpng')
close(gcf)

end