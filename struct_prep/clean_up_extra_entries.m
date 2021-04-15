function clean_up_extra_entries

% This function just exists to fix a bug wherein for early runs I forgot to
% reinitialize and clear spikes between each patient, and so patients with
% more hours or more files earlier on would have their extra data appended
% to future patients.

%% Get file locs
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
data_folder = [scripts_folder,'data/'];
addpath(genpath(scripts_folder));
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
sp_folder = [results_folder,'spikes/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

% Loop over spike files
listing = dir([sp_folder,'*mat']);
for i = 1:length(listing)
    spikes = load([sp_folder,listing(i).name]);
    spikes = spikes.spikes;
    
    
    % Get corresponding patient
    name = spikes.name;
    for p = 1:length(pt)
        if strcmp(name,pt(p).name) == 1
            break
        end
    end
    
    
    nfiles = length(pt(p).ieeg.file);
    if nfiles > 1
        error('Geez do something');
    end
    
    nhours = length(pt(p).ieeg.file.block); % this is how many hours it should have
    n_spike_hours = length(spikes.file.hour); % this is how many hours it has
    
    % Trim extra hours from spike
    if n_spike_hours > nhours % if more hours than it should have
        
        % trim all hours after the number of hours it should have
        spikes.file.hour(nhours+1:n_spike_hours) = [];
        
    end
    
    % Save the spike structure
    save([sp_folder,listing(i).name],'spikes');
    
    
end




end