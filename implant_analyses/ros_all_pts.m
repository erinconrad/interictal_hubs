function ros_all_pts(whichPts,saved_out)

%{
Generates Figure 2 for implant effect paper -> raster plot for a single
patient and rate order stability for all patients and stability of spikiest
electrode for all patients

%}

%% Parameters
surround = 24*1;
do_save = 1;
nb = 1e4;
do_rel = 0;
ex_p = 5;
do_vec = 0;
type = 'Spearman';

%% Decide whether to do this!!
only_pre = 0; % for the MC analysis, compare to only the pre-revision times (in case the revision effect is delayed)


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
    
    out = load([main_spike_results,'out.mat']);
    out = out.out;
    
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

%% Initialize figure
figure
set(gcf,'position',[100 100 800 700])
tiledlayout(3,2,'TileSpacing','tight','padding','compact')

%% Example spike rate raster
nexttile
p = ex_p;
rate_raster = out(p).rate./out(p).run_dur;
% Remove EKG and scalp electrodes
ekg = identify_ekg_scalp(out(p).unchanged_labels);
rate_raster(ekg,:) = [];
curr_times = (1:size(rate_raster,2)) * out(p).block_dur;
curr_change = out(p).change_block*out(p).block_dur;
h = turn_nans_white(rate_raster);
set(h,'XData',[0:curr_times(end)]);
xlim([0 curr_times(end)])
hold on
cp = plot([curr_change curr_change],ylim,'r--','linewidth',3);
xlabel('Hour')
%legend(cp,'Revision','fontsize',20,'location','southeast')
yticklabels([])
ylabel('Electrode')
set(gca,'fontsize',20)
c = colorbar('location','northoutside');
ylabel(c,'Spikes/min','fontsize',20)

%% Example node strength raster
nexttile
p = ex_p;
ns_raster = out(p).metrics.ns_norm;

% Remove bad electrodes per ns
rate_raster = out(p).rate;
ns_raster(isnan(rate_raster)) = nan;

% Remove EKG and scalp electrodes
ekg = identify_ekg_scalp(out(p).unchanged_labels);
ns_raster(ekg,:) = [];



curr_times = (1:size(ns_raster,2)) * out(p).block_dur;
curr_change = out(p).change_block*out(p).block_dur;
h = turn_nans_white((ns_raster));
set(h,'XData',[0:curr_times(end)]);
xlim([0 curr_times(end)])
hold on
cp = plot([curr_change curr_change],ylim,'r--','linewidth',3);
xlabel('Hour')
%legend(cp,'Revision','fontsize',20,'location','southeast')
yticklabels([])
%ylabel('Electrode')
set(gca,'fontsize',20)
c = colorbar('location','northoutside');
ylabel(c,'Normalized NS','fontsize',20)

%% Histogram of rate change for all electrodes, all patients
nexttile
npts = length(whichPts);
all_rate_diff = [];
for i = 1:npts
    rate = out(i).rate./out(p).run_dur;
    cblock = out(i).change_block;
    % Get surround times, starting with first non nan
    [pre,post] = get_surround_times(rate,cblock,surround);
    
    %pre = cblock - surround:cblock-1;
    %post = cblock+1:cblock+surround;
    rate_pre = nanmean(rate(:,pre),2);
    rate_post = nanmean(rate(:,post),2);
    if do_rel
        rate_diff = (rate_post-rate_pre)./rate_pre;
    else
        rate_diff = rate_post-rate_pre;
    end
    unchanged_labels = out(i).unchanged_labels;
    ekg = identify_ekg_scalp(unchanged_labels);
    rate_diff = rate_diff(~ekg);
    all_rate_diff = [all_rate_diff;rate_diff];
end
histogram(all_rate_diff);
ylabel('# electrodes')
xlabel('Spike rate change (spikes/min)')
set(gca,'fontsize',20)

%% Histogram of node strength change for all electrodes, all patients
nexttile
all_ns_diff = [];
for i = 1:npts
    ns = out(i).metrics.ns_norm;
    cblock = out(i).change_block;
    rate = out(i).rate;
    %pre = cblock - surround:cblock-1;
    %post = cblock+1:cblock+surround;
    
    % Get surround times, starting with first non nan
    [pre,post] = get_surround_times(rate,cblock,surround);
    
    ns_pre = nanmean(ns(:,pre),2);
    ns_post = nanmean(ns(:,post),2);
    if do_rel
        ns_diff = (ns_post-ns_pre)./ns_pre;
    else
        ns_diff = ns_post-ns_pre;
    end
    unchanged_labels = out(i).unchanged_labels;
    ekg = identify_ekg_scalp(unchanged_labels);
    ns_diff = ns_diff(~ekg);
    all_ns_diff = [all_ns_diff;ns_diff];
end
histogram(all_ns_diff);
xlabel('NS change')
%ylabel('Number of electrodes')
set(gca,'fontsize',20)

%% Rate order stability
nexttile

% Get full surround time for plotting
max_time_before = 0;
max_time_after = 0;
for i = 1:length(whichPts)
    ntotal = size(out(i).rate,2);
    cblock = out(i).change_block;
    nbefore = cblock-1;
    nafter = ntotal-cblock;
    if nbefore > max_time_before
        max_time_before = nbefore;
    end
    if nafter > max_time_after
        max_time_after = nafter;
    end
end

%all_ros = nan(length(whichPts),pl_surround*2+1);
all_ros = nan(length(whichPts),max_time_before+max_time_after+1);
mid_pos = max_time_before+1;
all_ps = nan(length(whichPts),1);
for i = 1:length(whichPts)
    rate = out(i).rate;
    cblock = out(i).change_block;
    
    % Remove EKG and scalp electrodes
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    rate(ekg,:) = [];
    
    % Take mean over pre-revision period
    %first_period = max(1,cblock-pl_surround):cblock-1;
    first_period = 1:cblock-1;
    first_rate = nanmean(rate(:,first_period),2);
    %nblocks = pl_surround*2+1;
    
    ros = nan(size(rate,2),1);
    for h = 1:size(rate,2)
        ros(h) = corr(first_rate,rate(:,h),'Type',type,'rows','pairwise');
    end
    
    % Put it at the correct position in the larger array. If it's the one
    % with the latest cblock, then we can start filling it up at the first
    % position. If it has an earlier cblock, then we need to pad the
    % beginning
    pos_off = mid_pos - cblock;
    all_ros(i,1+pos_off:length(ros)+pos_off) = ros;
    
    % do a ros permutation test to see whether the pre-post change is
    % larger than expected for randomly chosen times
    if do_vec
        pval_curr = compare_vecs(rate,cblock,surround,nb);
    else
        pval_curr = compare_rhos(rate,cblock,surround,nb,only_pre,type);
    end
    all_ps(i) = pval_curr;
    
end

% Fisher test to combine pvalues
X_2 = -2 * sum(log(all_ps));
sum_p = 1-chi2cdf(X_2,2*length(all_ps));

% Find the times in which at least 2 are non nans
m = nanmean(all_ros,1);
st = nanstd(all_ros,[],1);

two_non_nans = sum(~isnan(all_ros),1)>=2;
times = -max_time_before:max_time_after;
times(~two_non_nans) = [];
m(~two_non_nans) = [];
st(~two_non_nans) = [];


[mp,stp] = shaded_error_bars(times,m,st,[]);
hold on
set(gca,'fontsize',20)
xlabel('Hours surrounding revision')
ylabel('Spike stability')
yl = get(gca,'ylim');
ylim([yl(1) yl(1) + 1.11*(yl(2)-yl(1))])
plot([0 0],[yl(1) yl(1)+0.98*(yl(2)-yl(1))],'r--','linewidth',3)
plot([-surround surround],[yl(1)+1.01*(yl(2)-yl(1)) yl(1)+1.01*(yl(2)-yl(1))],'k-','linewidth',2)
text(0,yl(1)+1.03*(yl(2)-yl(1)),get_asterisks(sum_p,1),...
    'horizontalalignment','center','fontsize',30)

xlim([times(1) times(end)])

%% NS order stability
nexttile
all_ros = nan(length(whichPts),max_time_before+max_time_after+1);
all_ps = nan(length(whichPts),1);
for i = 1:length(whichPts)
    if isempty(out(i).metrics)
        ns = nan(size(out(i).rate));
    else
        ns = out(i).metrics.ns_norm;
    end
    cblock = out(i).change_block;
    
    % Remove EKG and scalp electrodes
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    ns(ekg,:) = [];
    
    % Take mean over pre-revision period
    first_period = 1:cblock-1;
    first_ns = nanmean(ns(:,first_period),2);
    
    %{
    nblocks = pl_surround*2+1;
       
    % Get rate order stability
    ros = nan(nblocks,1);
    hcount = 0;
    for h = cblock-pl_surround:cblock+pl_surround
        hcount = hcount + 1;
        ros(hcount) = corr(first_ns,ns(:,h),'Type','Spearman','rows','pairwise');
    end
    
    % put it in the cell
    all_ros(i,:) = ros;
    %}
    
    ros = nan(size(ns,2),1);
    for h = 1:size(ns,2)
        ros(h) = corr(first_ns,ns(:,h),'Type',type,'rows','pairwise');
    end

    pos_off = mid_pos - cblock;
    all_ros(i,1+pos_off:length(ros)+pos_off) = ros;
    
    % do a ros permutation test to see whether the pre-post change is
    % larger than expected for randomly chosen times
    if do_vec
        pval_curr = compare_vecs(rate,cblock,surround,nb);
    else
        pval_curr = compare_rhos(rate,cblock,surround,nb,only_pre,type);
    end
    all_ps(i) = pval_curr;
    
end

% Fisher test to combine pvalues
X_2 = -2 * nansum(log(all_ps));
sum_p = 1-chi2cdf(X_2,2*sum(~isnan(all_ps)));

% Find the times in which at least 2 are non nans
m = nanmean(all_ros,1);
st = nanstd(all_ros,[],1);

two_non_nans = sum(~isnan(all_ros),1)>=2;
times = -max_time_before:max_time_after;
times(~two_non_nans) = [];
m(~two_non_nans) = [];
st(~two_non_nans) = [];

[mp,stp] = shaded_error_bars(times,m,st,[0.8500, 0.3250, 0.0980]);
hold on
set(gca,'fontsize',20)
xlabel('Hours surrounding revision')
ylabel('NS stability')
yl = get(gca,'ylim');
ylim([yl(1) yl(1) + 1.11*(yl(2)-yl(1))])
plot([0 0],[yl(1) yl(1)+0.98*(yl(2)-yl(1))],'r--','linewidth',3)
plot([-surround surround],[yl(1)+1.01*(yl(2)-yl(1)) yl(1)+1.01*(yl(2)-yl(1))],'k-','linewidth',2)
text(0,yl(1)+1.03*(yl(2)-yl(1)),get_asterisks(sum_p,1),...
    'horizontalalignment','center','fontsize',30)

xlim([times(1) times(end)])

%% Spikiest electrode stability
%{
nexttile
nblocks = pl_surround*2+1;
agree_all = nan(length(whichPts),nblocks);
pval_all = nan(length(whichPts),1);
for i = 1:length(whichPts)
    rate = out(i).rate;
    cblock = out(i).change_block;
    
    % Remove EKG and scalp electrodes
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    rate(ekg,:) = [];
    
    % For each time period, get the spikiest electrode
    [~,spikiest] = max(rate,[],1);
    
    % Get mean rate in pre-revision period
    
    mean_pre_rate = nanmean(rate(:,cblock-pl_surround:cblock-1),2);
    [~,spikiest_pre] = max(mean_pre_rate,[],1);
    
    % For each time period, get whether spikiest electrode agrees with pre
    agree = nan(nblocks,1);
    hcount = 0;
    for h = cblock-pl_surround:cblock+pl_surround
        hcount = hcount + 1;
        agree(hcount) = isequal(spikiest_pre,spikiest(h));
    end
    
    agree_all(i,:) = agree;
    
    % Significance testing
    pval = compare_agree(rate,cblock,surround,nb);
    pval_all(i) = pval;
    
end

% Fisher
X_2 = -2 * sum(log(pval_all));
sum_p = 1-chi2cdf(X_2,2*length(pval_all));

% Plot the agree
times = -pl_surround:pl_surround;
plot(times,nanmean(agree_all,1),'linewidth',3)
hold on
set(gca,'fontsize',20)
xlabel('Hours surrounding revision')
ylabel('Spike hub stability')
yl = get(gca,'ylim');
ylim([yl(1) yl(1) + 1.11*(yl(2)-yl(1))])
plot([0 0],[yl(1) 0.98*(yl(2)-yl(1))],'r--','linewidth',3)
plot([-surround surround],[1.01*(yl(2)-yl(1)) 1.01*(yl(2)-yl(1))],'k-','linewidth',2)
text(0,1.03*(yl(2)-yl(1)),get_asterisks(sum_p,1),...
    'horizontalalignment','center','fontsize',30)

xlim([times(1) times(end)])

%}

%% Add annotations
annotation('textbox',[0 0.9 0.1 0.1],'String','A','fontsize',30,'linestyle','none')
annotation('textbox',[0.53 0.9 0.1 0.1],'String','B','fontsize',30,'linestyle','none')
annotation('textbox',[0 0.55 0.1 0.1],'String','C','fontsize',30,'linestyle','none')
annotation('textbox',[0.53 0.55 0.1 0.1],'String','D','fontsize',30,'linestyle','none')
annotation('textbox',[0 0.25 0.1 0.1],'String','E','fontsize',30,'linestyle','none')
annotation('textbox',[0.53 0.25 0.1 0.1],'String','F','fontsize',30,'linestyle','none')



if do_save == 1
    print(gcf,[main_spike_results,'Fig2'],'-dpng')
    print(gcf,[main_spike_results,'Fig2'],'-depsc')
end

end