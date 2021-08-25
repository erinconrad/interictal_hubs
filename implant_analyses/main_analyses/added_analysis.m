function added_analysis(whichPts,saved_out,out)

%% User change parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
main_surround = 3; %24 hour peri-revision surround
main_metric = 1;
ex_p = 1;

%% Other info
n_surrounds = length(all_surrounds);
all_metrics = {'rate','ns'};
n_metrics = length(all_metrics);

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

%% Get names
names = cell(length(whichPts),1);
for i = 1:length(whichPts)
    names{i} = out(i).name;
end

% Initialize output
n_patients = length(whichPts);
all_rates = cell(n_patients,2); % unchanged, added

% Get individual patient data
for i = 1:length(whichPts)
    
    rate = out(i).rate;
    added_rate = out(i).rate_added;
    cblock = out(i).change_block;
    
    % Remove EKG and scalp electrodes
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    rate(ekg,:) = [];
    run_dur = out(i).run_dur;
    
    % overall rate
    overall_r_un = nansum(rate,1);
    overall_r_ad = nansum(added_rate,1);
    
    if 1
        plot(overall_r_un)
        hold on
        plot(overall_r_ad);
        plot([cblock cblock],ylim,'--')
        pause
        hold off
    end
    
    
end

end