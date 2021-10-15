function methods_figure(out)

%% Parameters
p = 6;
surround = 24;
csize = 200;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
out_folder = [results_folder,'job_talk/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder);
end

addpath(genpath(locations.script_folder));


%% Get basic info
unchanged_labels = out(p).unchanged_labels;
unchanged_locs = out(p).unchanged_locs;
ekg = identify_ekg_scalp(unchanged_labels);
unchanged_labels(ekg) = [];
unchanged_locs(ekg,:) = [];
dist = out(p).dist;
block_dur = out(p).block_dur;

added_labels = out(p).added_labels;
added_locs = out(p).added_locs;

rate = out(p).rate/out(p).run_dur;
rate(ekg,:) = [];
cblock = out(p).change_block;
ns = out(p).metrics.ns;
ns(ekg,:) = [];
dist(ekg,:) = [];

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

%% Get distance between a few elecs
[~,max_dist_elec] = max(dist);
% Get which added elec it's closest to
all_dists = vecnorm(unchanged_locs(max_dist_elec,:) - added_locs,2,2);
[~,closest_max_dist] = min(all_dists);

[~,min_dist_elec] = min(dist);
% Get which added elec it's closest to
all_dists = vecnorm(unchanged_locs(min_dist_elec,:) - added_locs,2,2);
[~,closest_min_dist] = min(all_dists);


%% test plot
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,'k','linewidth',2);
hold on
%text(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),unchanged_labels);
scatter3(added_locs(:,1),added_locs(:,2),added_locs(:,3),csize,'rp','linewidth',2);

% Plot some nearby ones
[~,I] = (sort(dist));
scatter3(unchanged_locs(I(1:3),1),unchanged_locs(I(1:3),2),unchanged_locs(I(1:3),3),csize,'b','linewidth',2);
all_dists = vecnorm(unchanged_locs(I(1),:) - added_locs,2,2);
[~,closest_min_dist] = min(all_dists);
scatter3(added_locs(closest_min_dist,1),added_locs(closest_min_dist,2),added_locs(closest_min_dist,3),csize,'bp','linewidth',2);

%% Spike rate over time
figure
set(gcf,'position',[317 452 1124 345])
plot(nanmean(rate,1),'k','linewidth',2)
hold on
nan_blocks = find(isnan(nanmean(rate,1)));
for b = 1:length(nan_blocks)
    bidx = [max(nan_blocks(b) - 0.5,1) min(nan_blocks(b)+0.5,size(rate,2))];
    bidx = bidx;
    yl = ylim;
    ap = fill([bidx(1),bidx(2),bidx(2),bidx(1)],...
        [yl(1),yl(1),yl(2),yl(2)],...
        [0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
end
xlim([0 500])
ylim([0 9])
xticklabels([])
yticklabels([])
xlabel('Time since implant')
ylabel('Spike rate')
set(gca,'fontsize',20)
print([out_folder,'example_rate'],'-dpng');

%{

%% UNchanged locs
figure
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,'k','linewidth',2);
hold on
view([81,-6])
axis off
print([out_folder,'unchanged_locs'],'-dpng');

%% Unchanged and added
figure
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,'k','linewidth',2);
hold on
scatter3(added_locs(:,1),added_locs(:,2),added_locs(:,3),csize,'rp','linewidth',2);
view([81,-6])
axis off
print([out_folder,'unchanged_added_locs'],'-dpng');

%% Fake correlation (FAKE DATA)
figure
x = dist;
y = 1./dist + 0.1*rand(length(dist),1);
plot(x,y,'o','markersize',10,'linewidth',2);
xticklabels([])
yticklabels([])
ylabel('Spike rate change')
xlabel('Distance')
set(gca,'fontsize',20)
print([out_folder,'fake_correlation'],'-dpng');
%}
close all

if 0
%% Prep figure
figure
tiledlayout(3,2,'Padding','Compact','TileSpacing','Compact')

%% A: Brain pre-revision (top left)
nexttile
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,'k','linewidth',2);
view([81,-6])
axis off
%cs = colorbar;
%title('Relative spike rate change','fontsize',20);

%% B: Brain post-revision (top right)
nexttile
scatter3(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),csize,'k','linewidth',2);
hold on
%text(unchanged_locs(:,1),unchanged_locs(:,2),unchanged_locs(:,3),unchanged_labels);
scatter3(added_locs(:,1),added_locs(:,2),added_locs(:,3),csize,'rp','linewidth',2);

% Plot some nearby ones
[~,I] = (sort(out(p).dist));
scatter3(unchanged_locs(I(1:3),1),unchanged_locs(I(1:3),2),unchanged_locs(I(1:3),3),csize,'b','linewidth',2);

view([81,-6])
axis off
%cs = colorbar;
%title('Relative spike rate change','fontsize',20);

%% C: spike rate change (bottom left)
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

%% Hypothetical correlation between distance and rate change

%{
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
%}
end

end