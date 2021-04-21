function change_hour_to_block

% this is a one time use function to deal with the fact that I changed a
% field name from hour to block
%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
addpath(genpath(locations.ieeg_folder));
spike_folder = [results_folder,'old_april19/spikes_april19/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

whichPts = [];
if isempty(whichPts)
    listing = dir([spike_folder,'*.mat']);
    for i = 1:length(listing)
        C = listing(i).name;
        temp_name = strsplit(C,'_');
        temp_name = temp_name{1};
        for j = 1:length(pt)
            pt_name = pt(j).name;
            if strcmp(temp_name,pt_name)
                whichPts = [whichPts,j];
                break
            end
        end
    end
end

for p = whichPts
    pt_name = pt(p).name;


    %% Load spike file
    if exist([spike_folder,sprintf('%s_spikes.mat',pt_name)],'file') == 0
        continue;
    end
    
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
    name = spikes.name;
    nfiles = length(spikes.file);
    
    for f = 1:nfiles
        spikes.file(f).block = spikes.file(f).hour;
    end
    
    
    %% Save
    save([spike_folder,sprintf('%s_spikes.mat',pt_name)],'spikes');
end



end