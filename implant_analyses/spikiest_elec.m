function spikiest_elec(rate,labels,change,surround,results_folder,name,...
    block_dur,run_dur)

%% Remove EKG and scalp electrodes
ekg = identify_ekg_scalp(labels);
rate(ekg,:) = [];
labels(ekg) = [];

%% Redefine rate according to run_dur
rate = rate/run_dur;

%% Remove electrodes with all nans
all_nan = sum(~isnan(rate),2) == 0;
rate(all_nan,:) = [];
labels(all_nan) = [];
 
%% Identify pre times and post times
pre = 1:change-1;
post = change + 1:size(rate,2);
 
%% Get mean pre- rate for each channel and sort channels by this
mean_rate = nanmean(rate(:,pre),2);
[sorted_mean_rate,sorted_chs] = sort(mean_rate,'descend');


%% Re-sort labels and rate by their pre-rate
labels = labels(sorted_chs);
rate = rate(sorted_chs,:);

%% Get the spikiest electrode pre-revision and post revision
[~,spikiest_pre] = max(sorted_mean_rate);
[~,spikiest_post] = max(nanmean(rate(:,post),2));

%% For each time period, get the correlation between its rate order and this mean pre-rate order
corr_time = nan(size(rate,2),1);
for i = 1:size(rate,2)
    corr_time(i) = corr(rate(:,i),sorted_mean_rate,'Type','Spearman',...
        'Rows','pairwise');
end


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
t = tiledlayout(2,1,'TileSpacing','compact');

%% Show raster
nexttile(1)
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
if spikiest_pre == spikiest_post
    yticks(spikiest_pre);
    yticklabels('

%% Show correlation with pre-revision rate-order
nexttile(2)
plot(times,corr_time,'linewidth',2);
hold on
cp = plot([change*block_dur change*block_dur],ylim,'r--','linewidth',3);
ylabel({'Spike rate consistency'});
xlabel('Hour')
set(gca,'fontsize',20)

%{
%% Compare pre- and post-revision spikiest channel
nexttile(4)
bar([1 2],[pre_sp_mode_agg*100 post_sp_mode_agg*100])
xticks([1 2])
xlim([0.5 2.5])
xticklabels({'Pre-revision','Post-revision'})
ylabel('Spikiest electrode consistency (%%)');
set(gca,'fontsize',20);
%}

title(t,name,'fontsize',20)

outfolder = [results_folder,'consistency/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end
print(gcf,[outfolder,name],'-dpng')
close(gcf)

end