function find_new_spikes(whichPts,saved_out)

%% Parameters
max_pre = 0.01;% in units of spikes/min
min_post = 5;

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

for i = 1:length(whichPts) 
    rate = out(i).rate./out(i).run_dur;
    cblock = out(i).change_block;
    chLabels = out(i).unchanged_labels;
    name = out(i).name;
    
    %% Find spikes meeting criteria for few pre and many post
    new_spikes = nanmean(rate(:,1:cblock-1),2) < max_pre & ...
        nanmean(rate(:,cblock+1:end),2) > min_post;
    
    name
    chLabels(new_spikes)
    
end

end