function get_spikes(whichPts)

%% Parameters
overwrite = 0;
test.do_test = 0;
do_plot = 0;

%% Test parameters
test.time = 20550; %193231.80;% - ok except RA4;  %1605.05 - super high variance LAF1; 5617.68 - flat; 28583.69 - RA4 bad
test.dur = 60;
test.pt = 1;
test.file = 1;
test.ch = [];

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
    
    if test.do_test == 1
        p = test.pt;
        name = pt(p).name;
        next_file = 1;
        next_hour = 1;
    else
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
    end
    
    
    
    %% Pull spike detector parameters
    params = pull_detector_params(name);
    
    
    % Loop over ieeg files
    n_files = length(pt(p).ieeg.file);
    for f = next_file:n_files
        
        if test.do_test == 1
            f = test.file;
        end
        
        % sampling frequency
        fs = pt(p).ieeg.file(f).fs;
        
        % electrodes
        chLabels = pt(p).ieeg.file(f).chLabels(:,1);
        
        % filename
        fname = pt(p).ieeg.file(f).name;
        
        % Loop over hours
        n_hours = length(pt(p).ieeg.file(f).block);
        for h = next_hour:n_hours
            tic;
            fprintf('\nDoing %s file %d hour %d\n',name,f,h);
            
            % get the run time (already randomly assigned minute)
            if test.do_test == 1
                run_times = [test.time,test.time+test.dur];
            else
                run_times = pt(p).ieeg.file(f).block(h).run; 
            
            end
            run_idx = run_times(1)*fs:run_times(2)*fs;
            
            
            %% Get the eeg data
            session = IEEGSession(fname, login_name, pwfile);
            values = session.data.getvalues(run_idx,':');
                       
            %% Do pre-processing
            [values,bipolar_labels] = pre_process(values,chLabels);
            
            %% Designate electrodes over which to run spike detector
            if test.do_test == 1 && ~isempty(test.ch)
                if iscell(test.ch)
                    which_chs = (cellfun(@(x) find(strcmp(chLabels,x)),test.ch));
                else
                    which_chs = test.ch;
                end
            else
                [which_chs,skip] = designate_chs(chLabels);
            end
            
            %% Reject bad channels
            bad = reject_bad_chs(values,which_chs,chLabels,fs);
            
            %% Spike detector
            run_chs = which_chs;
            run_chs(ismember(run_chs,bad)) = [];
            if ~isempty(run_chs)
                gdf = detector(values,fs,run_chs,params);
            else
                gdf = [];
            end
            
            %% Multi-channel requirements
            if ~isempty(gdf)
                gdf =  multi_channel_requirements(gdf,length(which_chs),fs);
            end
            
            %% Example plot
            dur = diff(run_times);
            t = toc;
            fprintf('\nTook %1.1f s and detected %d spikes\n',t,size(gdf,1));
            if do_plot
                show_eeg_and_spikes(values,bipolar_labels,gdf,dur,run_times(1),name,fs,bad,skip);
            end
          
            %% Re-align gdf time
            if ~isempty(gdf)
                gdf(:,2) = (gdf(:,2)+run_idx(1)-1)/fs;
            end
            
            %% Add spikes to structure
            spikes.file(f).hour(h).run_times = run_times;
            spikes.file(f).hour(h).params = params;
            spikes.file(f).hour(h).gdf = gdf;
            spikes.file(f).hour(h).chLabels = chLabels;
            spikes.file(f).hour(h).bipolar_labels = bipolar_labels;
            spikes.file(f).hour(h).bad = bad;
            spikes.file(f).hour(h).skip = skip;
            spikes.file(f).hour(h).run_chs = run_chs;
            
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
            if test.do_test == 0
                save([sp_folder,out_file],'spikes');
            end
           
            
        end
    end
    
end



end