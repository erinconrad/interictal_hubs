function extreme_all_pts(whichPts,saved_out)


%% Parameters
surround = 48;
nb = 1e4;
min_rate = 0;
thresh_prc_dist = [20 80];
thresh_prc_ns = [20 80];
do_rel = 0;
do_med = 0;

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


%% Get percentiles of distance and connectivity
all_dist = [];
all_fc = [];
for i = 1:length(whichPts)
    dist = out(i).dist;
    fc = out(i).metrics.added_pc;
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    dist(ekg) = [];
    fc(ekg) = [];
    
    all_dist = [all_dist;dist];
    all_fc = [all_fc;fc];
end

%% Get the percentiles of most close and most connected
dist_prc = prctile(all_dist,thresh_prc_dist);
fc_prc = prctile(all_fc,thresh_prc_ns);
thresh_dist = dist_prc;
thresh_ns = fc_prc;
    
if 1
    figure
    tiledlayout(1,2)
    nexttile
    histogram(all_dist)
    hold on
    plot([dist_prc(1) dist_prc(1)],ylim)
    plot([dist_prc(2) dist_prc(2)],ylim)
    
    nexttile
    histogram(all_fc)
    hold on
    plot([fc_prc(1) fc_prc(1)],ylim)
    plot([fc_prc(2) fc_prc(2)],ylim)
end

npts = length(whichPts);
all_rclose = nan(npts,2);
all_rconn = nan(npts,2);
all_nsclose = nan(npts,2);
all_nsconn = nan(npts,2);

for i = 1:length(whichPts) 
    
    rate = out(i).rate./out(i).run_dur;
    if isempty(out(i).metrics)
        ns = nan(size(rate));
    else
        ns = out(i).metrics.ns_norm;
    end
    cblock = out(i).change_block;
    
    % Identify pre and post times
    if isempty(surround)
        pre = 1:cblock-1;
        post = cblock+1:size(rate,2);
    else
        pre = cblock-surround:cblock - 1;
        post = cblock + 1: cblock+surround;
    end
    
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    
    %% Get responses
    abs_rate_change = nanmean(rate(:,post),2) - nanmean(rate(:,pre),2);
    abs_ns_change = nanmean(ns(:,post),2) - nanmean(ns(:,pre),2);
    
    if do_rel
        abs_rate_change = abs_rate_change./nanmean(rate(:,pre),2);
        abs_ns_change = abs_ns_change./nanmean(ns(:,pre),2);
    end
    
    dist = out(i).dist;
    fc = out(i).metrics.added_pc;
 
    % Remove ekg
    abs_rate_change(ekg) = [];
    abs_ns_change(ekg) = [];
    dist(ekg) = [];
    fc(ekg) = [];
    
    % Find close and connected things
    close = dist < thresh_dist(1);
    conn = fc > thresh_ns(2);
    
    not_close = dist > thresh_dist(2);
    not_conn = fc < thresh_ns(1);
    
    % Get responses of these
    if do_med
        rclose = [nanmedian(abs_rate_change(close)) nanmedian(abs_rate_change(not_close))];
        rconn = [nanmedian(abs_rate_change(conn)) nanmedian(abs_rate_change(not_conn))];
        nsclose = [nanmedian(abs_ns_change(close)) nanmedian(abs_ns_change(not_close))];
        nsconn = [nanmedian(abs_ns_change(conn)) nanmedian(abs_ns_change(not_conn))];
    else
        rclose = [nanmean(abs_rate_change(close)) nanmean(abs_rate_change(not_close))];
        rconn = [nanmean(abs_rate_change(conn)) nanmean(abs_rate_change(not_conn))];
        nsclose = [nanmean(abs_ns_change(close)) nanmean(abs_ns_change(not_close))];
        nsconn = [nanmean(abs_ns_change(conn)) nanmean(abs_ns_change(not_conn))];
    end
    
    all_rclose(i,:) = rclose;
    all_rconn(i,:) = rconn;
    all_nsclose(i,:) = nsclose;
    all_nsconn(i,:) = nsconn;
    
    
    
end

figure
tiledlayout(2,2)

nexttile
plot(ones(npts,1)+0.05*rand(npts,1),all_rclose(:,1),'o')
hold on
plot(2*ones(npts,1)+0.05*rand(npts,1),all_rclose(:,2),'o')
[~,p] = ttest(all_rclose(:,1),all_rclose(:,2));
xl = xlim;
yl = ylim;
text(xl(2),yl(2),sprintf('p = %1.3f',p),...
    'horizontalalignment','right','verticalalignment','top')

nexttile
plot(ones(npts,1)+0.05*rand(npts,1),all_rconn(:,1),'o')
hold on
plot(2*ones(npts,1)+0.05*rand(npts,1),all_rconn(:,2),'o')
[~,p] = ttest(all_rconn(:,1),all_rconn(:,2));
xl = xlim;
yl = ylim;
text(xl(2),yl(2),sprintf('p = %1.3f',p),...
    'horizontalalignment','right','verticalalignment','top')

nexttile
plot(ones(npts,1)+0.05*rand(npts,1),all_nsclose(:,1),'o')
hold on
plot(2*ones(npts,1)+0.05*rand(npts,1),all_nsclose(:,2),'o')
[~,p] = ttest(all_nsclose(:,1),all_nsclose(:,2));
xl = xlim;
yl = ylim;
text(xl(2),yl(2),sprintf('p = %1.3f',p),...
    'horizontalalignment','right','verticalalignment','top')

nexttile
plot(ones(npts,1)+0.05*rand(npts,1),all_nsconn(:,1),'o')
hold on
plot(2*ones(npts,1)+0.05*rand(npts,1),all_nsconn(:,2),'o')
[~,p] = ttest(all_nsconn(:,1),all_nsconn(:,2));
xl = xlim;
yl = ylim;
text(xl(2),yl(2),sprintf('p = %1.3f',p),...
    'horizontalalignment','right','verticalalignment','top')

end