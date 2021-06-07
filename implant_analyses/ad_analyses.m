function ad_analyses(whichPts,saved_out)


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


%% Power spectrum of AD
% Frequency on x axis, power on y axis; all pts, look for peak at 1/24
% hours??? 

%% Single pt example spike rate vs AD

%% All pts correlation between AD and SR

%% AD over time all patients on graph, showing implant time in middle???

%% Single patient sleep ROS
% Either manually pick a sleep period, or automatically use all times with
% AD ratio below some threshold, and get rate order for those times and
% then do ROS over time for single patient

%% ALl patients, compare pre-revision ROS to post-revision ROS


end