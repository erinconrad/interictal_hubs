function anatomy_analyses(whichPts,saved_out)

%% Parameters
surround = 48;

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

figure
tiledlayout(3,5)

%% Get anatomy loc names
[~,~,loc_names] = anatomy_grouper([]);

%% prep table
all_locs = array2table(nan(length(whichPts),length(loc_names)), 'VariableNames',loc_names);

for i = 1:length(whichPts)
    
    nexttile
    
    %% unchanegd labels, cblock
    unchanged_labels = out(i).unchanged_labels;
    cblock = out(i).change_block;
    rate = out(i).rate;
    
    %% Define pre and post
    [pre,post] = get_surround_times(rate,cblock,surround);
    
    %% Get ekg
    ekg = identify_ekg_scalp(unchanged_labels);
    
    %% Get anatomy
    unchanged_anatomy = out(i).unchanged_anatomy;
    
    %% Get rate change
    rate = out(i).rate;
    rate_change = nanmean(rate(:,post),2) - nanmean(rate(:,pre),2);
    
    %% Get ns change
    ns = out(i).metrics.ns_norm;
    ns_change = nanmean(ns(:,post),2) - nanmean(ns(:,pre),2);
    
    %% Remove ekg
    unchanged_anatomy(ekg) = [];
    rate_change(ekg) = [];
    ns_change(ekg) = [];

    %% Group by localization
    [~,ana_loc] = anatomy_grouper(unchanged_anatomy);
    % group groups
    [ana_loc_groups,~,ic] = unique(ana_loc);
    nLocs = length(ana_loc_groups);
    grouped_rates_loc = cell(nLocs,1);
    rates_table = [];
    ns_table = [];
    loc_table = {};
   
    
    % loop through anatomy groups
    for j = 1:nLocs

        % get the channels corresponding to that anatomy group
        chs = find(ic == j);

        % Get the average rate change for those chs
        grouped_rates_loc{j} = (rate_change(chs));
        rates_table = [rates_table;rate_change(chs)];
        loc_table = [loc_table;ana_loc(chs)];
        ns_table = [ns_table;ns_change(chs)];
    end
    
    % Prep table for anova
    [p_rate,~,~] = anova1(rates_table,loc_table,'off');
    [p_ns,~,~] = anova1(ns_table,loc_table,'off');
    
    for k = 1:length(grouped_rates_loc)
        plot(k+0.05*rand(length(grouped_rates_loc{k}),1),grouped_rates_loc{k},'o',...
            'markersize',15)
        hold on
    end
    xticks(1:nLocs)
    xticklabels(ana_loc_groups)
    xtickangle(45)
    xlim([0.5 length(grouped_rates_loc)+0.5])
    xl = xlim; yl = ylim;
    text(xl(2),yl(2),sprintf('p = %1.3f',p_rate),'horizontalalignment','right','verticalalignment','top');
    
    %% Get mean rate for each group
    for j = 1:length(loc_names)
        
        % Get the indices with that group
        curr_idx = strcmp(ana_loc,loc_names{j});
        
        % Mean rate change for that
        all_locs.(loc_names{j})(i) = nanmean(rate_change(curr_idx));
        
    end
    
end

nexttile([1,3])
all_locs_array = table2array(all_locs);
p = anova1(all_locs_array,[],'off');

for k = 1:size(all_locs_array,2)
    plot(k+0.05*rand(length(all_locs_array(:,k)),1),all_locs_array(:,k),'o',...
        'markersize',15)
    hold on
end

xticks(1:size(all_locs_array,2))
xticklabels(loc_names)
xtickangle(45)
xlim([0.5 length(loc_names)+0.5])

xl = xlim; yl = ylim;
text(xl(2),yl(2),sprintf('p = %1.3f',p),'horizontalalignment','right','verticalalignment','top');
    


end