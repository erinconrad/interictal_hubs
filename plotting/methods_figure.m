function methods_figure(out)

%% Parameters
p = 1;
surround = 24;
csize = 200;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
main_spike_results = [results_folder,'main_spikes/'];
if ~exist(main_spike_results,'dir')
    mkdir(main_spike_results);
end

addpath(genpath(locations.script_folder));


%% Get basic info
unchanged_labels = out(p).unchanged_labels;
unchanged_locs = out(p).unchanged_locs;
ekg = identify_ekg_scalp(unchanged_labels);
unchanged_labels(ekg) = [];
unchanged_locs(ekg,:) = [];

added_labels = out(p).added_labels;
added_locs = out(p).added_locs;

rate = out(p).rate/out(p).run_dur;
rate(ekg,:) = [];
cblock = out(p).change_block;
ns = out(p).metrics.ns;
ns(ekg,:) = [];

% Get surround times, starting with first non nan
[pre,post] = get_surround_times(rate,cblock,surround);

% Get relative rate change
pre_rate = nanmean(rate(:,pre),2);
post_rate = nanmean(rate(:,post),2);
rel_rate_change = (post_rate-pre_rate)./abs(pre_rate);

% Get relative ns change
pre_ns = nanmean(ns(:,pre),2);
post_ns = nanmean(ns(:,post),2);
rel_ns_change = (post_ns-pre_ns)./abs(pre_ns);


%% Prep figure
figure
tiledlayout(3,2,'Padding','Compact','TileSpacing','Compact')

%% A: Brain pre-revision (top left)

%% B: Brain post-revision (top right)

%% C: Spikes (middle left)

%% D: spike rate change (middle right)
nexttile
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,rel_rate_change,'filled');
hold on
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,'k','linewidth',2);
%text(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),unchanged_labels);
scatter3(added_locs(:,1),added_locs(:,2),added_locs(:,3),csize,'rp','linewidth',2);
view([81,-6])
axis off
cs = colorbar;
title('Relative spike rate change','fontsize',20);


%% E: PC network (bottom left)


%% F: Node strength change (bottom right)
nexttile
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,rel_ns_change,'filled');
hold on
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,'k','linewidth',2);
%text(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),unchanged_labels);
scatter3(added_locs(:,1),added_locs(:,2),added_locs(:,3),csize,'rp','linewidth',2);
view([81,-6])
axis off
cns = colorbar;
title('Relative node strength change','fontsize',20);


end