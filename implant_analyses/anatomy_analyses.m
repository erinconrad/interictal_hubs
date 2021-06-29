function anatomy_analyses(whichPts,saved_out)

%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
do_norm = 0; % doesn't seem to make much difference so I will keep 0 for simplicity
main_surround = 3;

%% Other info
n_surrounds = length(all_surrounds);
all_metrics = {'spike rate','node strength'};
n_metrics = length(all_metrics);

% divide by 2 to get time in hours but then multiply by 2 to get the full
% pre+post surround
true_times = all_surrounds; 

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

%% Get anatomy loc names
[~,~,loc_names] = anatomy_grouper([]);

%% Prep final table
final_array = cell(length(all_metrics),n_surrounds);


for im = 1:n_metrics
    
    metric = all_metrics{im};
    
    for is = 1:n_surrounds
        
        surround = all_surrounds(is);
        
        %% prep table
        all_locs = array2table(nan(length(whichPts),length(loc_names)), 'VariableNames',loc_names);

        
        for i = 1:length(whichPts)


            %% unchanegd labels, cblock
            unchanged_labels = out(i).unchanged_labels;
            cblock = out(i).change_block;
            switch metric
                case 'spike rate'
                    rate = out(i).rate;
                case 'node strength'
                    if do_norm
                        rate = out(i).metrics.ns_norm;
                    else
                        rate = out(i).metrics.ns;
                    end
            end
            
            %% Get anatomy
            unchanged_anatomy = out(i).unchanged_anatomy;

            %% Remove EKG and scalp electrodes
            ekg = identify_ekg_scalp(out(i).unchanged_labels);
            rate(ekg,:) = [];
            unchanged_anatomy(ekg) = [];
            
            %% Define pre and post
            [pre,post] = get_surround_times(rate,cblock,surround);

            

            %% Get rate/ns change
            rate_change = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre));

            %% Group by localization
            [~,ana_loc] = anatomy_grouper(unchanged_anatomy);
            % group groups
            [ana_loc_groups,~,ic] = unique(ana_loc);
            nLocs = length(ana_loc_groups);
            grouped_rates_loc = cell(nLocs,1);
            rates_table = [];
            loc_table = {};


            % loop through anatomy groups
            for j = 1:nLocs

                % get the channels corresponding to that anatomy group
                chs = find(ic == j);

                % Get the average rate change for those chs
                grouped_rates_loc{j} = (rate_change(chs));
                rates_table = [rates_table;rate_change(chs)];
                loc_table = [loc_table;ana_loc(chs)];
            end

            %% Get mean rate for each group
            for j = 1:length(loc_names)

                % Get the indices with that group
                curr_idx = strcmp(ana_loc,loc_names{j});

                % Mean rate change for that
                all_locs.(loc_names{j})(i) = nanmean(rate_change(curr_idx));

            end

        end
        
        % Convert to array
        all_locs_array = table2array(all_locs);
        
        % for Friedman, remove unspecified category
        all_locs_minus_unspecified = all_locs_array(:,[2:end]);
        [p stats] = skillmack(all_locs_minus_unspecified,1);
        chi2 = stats.T;
        df = stats.df;
        
        % make a string out of this
        curr_string = sprintf('chi^2(%d) = %1.1f, p = %1.2f',df,chi2,p);
        final_array{im,is} = (curr_string);
        
        % Show full result if main surround
        if is == main_surround
            fprintf(['\nWhen the peri-implant surround period was defined to'...
                ' be %d hours, the mean (std) relative %s change was %1.1f (%1.1f)'...
                ' for %s, %1.1f (%1.1f) for %s, %1.1f (%1.1f) for %s, and %1.1f (%1.1f) for %s.'...
                ' The difference between groups was not significant (Skillings-Mack %s)\n'],...
                all_surrounds(is),all_metrics{im},...
                nanmean(all_locs_array(:,1)),nanstd(all_locs_array(:,1)),loc_names{2},...
                nanmean(all_locs_array(:,2)),nanstd(all_locs_array(:,2)),loc_names{3},...
                nanmean(all_locs_array(:,3)),nanstd(all_locs_array(:,3)),loc_names{4},...
                nanmean(all_locs_array(:,4)),nanstd(all_locs_array(:,4)),loc_names{5},...
                curr_string);
            
            % Plot it for my own understanding
            figure
            for k = 1:size(all_locs_array,2)
                plot(k+rand(size(all_locs_array,1),1)*0.05,all_locs_array(:,k),'o')
                hold on

            end
            title(all_metrics{im})
        end
        
    end
end

final_table = cell2table(final_array,...
    'VariableNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'RowNames',all_metrics);


writetable(final_table,[main_spike_results,'anatomy.csv'],'WriteRowNames',true)  


end