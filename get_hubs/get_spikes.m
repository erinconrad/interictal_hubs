function get_spikes(whichPts)

%% Parameters
overwrite = 0;

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];
scripts_folder = [locations.main_folder,'scripts/'];
addpath(genpath(scripts_folder));
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
sp_folder = [results_folder,'spikes/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Get which patients to run
if isempty(whichPts)
    whichPts = 1:length(pt);
end

% Loop over patients
for i = 1:length(whichPts)
    p = whichPts(i);
    name = pt(p).name;
    
    out_file = [name,'_spikes.mat'];
    
    %% Load the existing spike structure if it exists
    if overwrite == 0 && exist([sp_folder,out_file],'file') ~= 0
        spikes = load([sp_folder,out_file]);
        spikes = spikes.spikes;
        
        % find the last hour we have finished
        next_file = spikes.next_file;
        next_hour = spikes.next_hour;
        
        if isnan(next_file)
            fprintf('\nAlready finished %s, skipping...\n',name);
            continue
        end
    else
        % initialize it
        next_file = 1;
        next_hour = 1;
        spikes.name = name;
      
    end
    
    
    
    %% Pull spike detector parameters
    params = pull_detector_params(name);
    
    
    % Loop over ieeg files
    n_files = length(pt(p).ieeg.file);
    for f = next_file:n_files
        
        % sampling frequency
        fs = pt(p).ieeg.file(f).fs;
        
        % electrodes
        chLabels = pt(p).ieeg.file(f).chLabels(:,1);
        
        % filename
        fname = pt(p).ieeg.file(f).name;
        
        % Loop over hours
        n_hours = length(pt(p).ieeg.file(f).block);
        for h = next_hour:n_hours
            
            % get the run time (already randomly assigned minute
            run_times = pt(p).ieeg.file(f).block(h).run; 
            run_idx = run_times(1)*fs:run_times(2)*fs;
            
            %% Get the eeg data
            session = IEEGSession(fname, login_name, pwfile);
            values = session.data.getvalues(run_idx,':');
            
            %% Do pre-processing
            values = pre_process(values,chLabels);
            
            %% Designate electrodes over which to run spike detector
            which_chs = designate_chs(chLabels);
            
            %% Reject bad channels?????
            
            %% Spike detector
            gdf = detector(values,fs,which_chs,params);
            
            %% Example plot
            dur = diff(run_times);
            which_chs = 1:length(chLabels);
            show_eeg_and_spikes(values,which_chs,chLabels,gdf,dur,run_times(1),name,fs);
          
            
            %% Re-align gdf time
            gdf(:,2) = (gdf(:,2)+run_idx(1)-1)/fs;
            
            %% Add spikes to structure
            spikes.file(f).hour(h).run_times = run_times;
            spikes.file(f).hour(h).params = params;
            spikes.file(f).hour(h).gdf = gdf;
            
            %% Mark the next hour
            if h == n_hours
                next_hour = 1;
                next_file = f + 1;
                if f == n_files
                    next_file = nan;
                end
            else
                next_hour = h + 1;     
            end
            spikes.next_file = next_file;
            spikes.next_hour = next_hour;
            
            %% Save structure
            save([sp_folder,out_file],'spikes');
           
            
        end
    end
    
end



end