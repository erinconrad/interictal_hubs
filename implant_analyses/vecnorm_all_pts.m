function vecnorm_all_pts(whichPts,saved_out)

%% Parameters
surround = 24*2;
do_save = 1;
nb = 1e4;
do_rel = 0;
ex_p = 5;

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

%% Spike distribution stability



end