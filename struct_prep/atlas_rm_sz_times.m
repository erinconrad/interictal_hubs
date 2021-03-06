function atlas_rm_sz_times(overwrite)

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
sp_folder = [results_folder,'new_spikes_june6_2021/'];
%sp_folder = [results_folder,'new_spikes/'];
out_folder = [results_folder,'clean_atlas_spikes/'];

if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

%% LOad pt
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Loop through pt file
for i = 1:length(pt)
    name = pt(i).name;
    
    % Find corresponding spike file
    sp_file_name = [name,'_spikes.mat'];
    sp_path = [sp_folder,sp_file_name];
    
    % See if outfile already exists
    if overwrite == 0
        if exist([out_folder,sp_file_name],'file') ~= 0
            fprintf('\nAlready did %s, skipping...\n',name);
            continue
        end
    end
    
    % see if it exists
    if ~exist(sp_path,'file')
        fprintf('\nDid not find %s\n',name);
        continue
    end
    
    % Load it
    spikes = load(sp_path);
    spikes = spikes.spikes;
    
    old_spikes = spikes;
    
    % Get seizure times
    all_sz = [];
    for s = 1:length(pt(i).seizure_info.sz)
        eec = pt(i).seizure_info.sz(s).EEC;
        sz_end = pt(i).seizure_info.sz(s).End;
        if isempty(sz_end)
            sz_end = eec + 5*60; % assume five minutes;
        end
        all_sz = [all_sz;eec sz_end];
    end
    
    % Loop over files
    all_n_removed = 0;
    for f = 1:length(spikes.file)
        for h = 1:length(spikes.file(f).block)
            gdf = spikes.file(f).block(h).gdf;
            
            % Remove spikes in seizures
            [gdf,n_removed]= remove_spikes_in_sz(gdf,all_sz);
            all_n_removed = all_n_removed + n_removed;
            
            % replace gdf with the clean version
            spikes.file(f).block(h).gdf = gdf;
            
            % Empty this as I'm not using it
            spikes.file(f).block(h).details = [];
            
            
        end
    end
    spikes.n_removed_in_sz = all_n_removed;
    spikes.sz_times = all_sz;
    fprintf('\n\nFor %s, sz times are:\n',name);
    all_sz
    fprintf('\nFor %s, removed %d spikes\n\n',name,all_n_removed)
    
    
    
    % Save
    save([out_folder,sp_file_name],'spikes');
end

end