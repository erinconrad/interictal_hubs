function ros_cat(whichPts,saved_out)

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


figure

for i = 1:length(whichPts)
    rate = out(i).rate;
    cblock = out(i).change_block;
    
    % Remove EKG and scalp electrodes
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    rate(ekg,:) = [];
    
   
    
    [true_cat,pseudo_cats] = mc_cat(rate,cblock,surround,nb);
    
    plot(true_cat)
    hold on
    shaded_error_bars(1:length(true_cat),mean(pseudo_cats,1),std(pseudo_cats,[],1));
    
end



end