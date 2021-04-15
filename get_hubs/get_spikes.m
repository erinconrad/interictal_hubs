function get_spikes(whichPts)

%% Parameters
overwrite = 0;
test.do_test = 1;
do_plot = 1;

%% Test parameters
test.pt = 21;
test.time = 1170334.89; %193231.80;% - ok except RA4;  %1605.05 - super high variance LAF1; 5617.68 - flat; 28583.69 - RA4 bad
test.dur = 15;
test.file = 1;
test.ch = [];

test.tmul = 17;
test.absthresh = 50;

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

if exist(sp_folder,'dir') == 0
    mkdir(sp_folder);
end

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load parameter file
param_table = readtable([data_folder,'detector_parameters.xlsx']);

%% Load localization structure
loc = load([data_folder,'patient_localization.mat']);
loc = loc.patient_localization;

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
            clear spikes % I must clear this 
            next_file = 1;
            next_hour = 1;
            spikes.name = name;

        end
    end
    
    %% identify the correct index of the localization structure
    for loc_p = 1:length(loc)
        if strcmp(loc(loc_p).patient,name) == 1
            break
        end
    end
    
    % Get clean loc labels
    clean_loc_labs = clean_labels_2(loc(loc_p).labels);
    
    %% Pull spike detector parameters
    params = pull_detector_params(name,param_table);
    
    if test.do_test == 1
        params.tmul = test.tmul;
        params.absthresh = test.absthresh;
    end
    
    % Loop over ieeg files
    if isempty(pt(p).ieeg)
        n_files = 0;
    else
        n_files = length(pt(p).ieeg.file);
    end
    for f = next_file:n_files
        
        if test.do_test == 1
            f = test.file;
        end
        
        % sampling frequency
        fs = pt(p).ieeg.file(f).fs;
        
        % electrodes
        chLabels = pt(p).ieeg.file(f).chLabels(:,1);
        
        %% Get cleaned labels
        clean_labs = clean_labels_2(chLabels);
        
        %% Reconcile cleaned labels with cleaned loc labels
        % For the purpose of knowing which electrodes to skip
        if ~isequal(clean_labs,clean_loc_labs)
            fprintf('\nWarning, localization labels do not match ieeg labels\n');
        end
        
        
        %% Designate which electrodes to skip
        [which_chs,skip] = designate_chs(chLabels,clean_labs,clean_loc_labs,loc(loc_p));
        non_skip = which_chs;
        
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
            dur = diff(run_times);
            
            
            %% Get the eeg data
            session = IEEGSession(fname, login_name, pwfile);
            values = session.data.getvalues(run_idx,':');
            session.delete;
                       
            %% Do pre-processing
            [values,bipolar_labels] = pre_process(values,clean_labs);
            
            %% Designate electrodes over which to run spike detector
            if test.do_test == 1 && ~isempty(test.ch)
                if iscell(test.ch)
                    which_chs = (cellfun(@(x) find(strcmp(chLabels,x)),test.ch));
                else
                    which_chs = test.ch;
                end
            else
                
                % I have already designated channels to run about at the
                % file level
            end
            
            %% Reject bad channels
            [bad,bad_details] = reject_bad_chs(values,which_chs,chLabels,fs);
            
            %% Skip the chunk entirely if enough channels are bad (suggests period of disconnection)
            if length(bad) >= 0.5*length(which_chs)
                fprintf('\nSkipping this run because %d of %d chs marked bad\n',length(bad),length(which_chs));
                gdf = [];
                run_chs = [];
                run_skip = 1;
            else
                run_skip = 0;
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
                    gdf =  multi_channel_requirements(gdf,length(run_chs),fs);
                end

                

                
            end
            
            %% Run details
            t = toc;
            fprintf('\nTook %1.1f s and detected %d spikes\n',t,size(gdf,1));
            fprintf('Of %d non-skipped chs, rejected %d for nans, %d for zeros,\n%d for variance, %d for noise, %d for std.\n',...
                length(non_skip),length(bad_details.nans),length(bad_details.zeros),length(bad_details.var),...
                length(bad_details.noisy),length(bad_details.higher_std));
            
            %% Example plot              
            if do_plot
                show_eeg_and_spikes(values,bipolar_labels,gdf,dur,run_times(1),name,fs,bad,skip,params);
            end
            
            %% Re-align gdf time
            if ~isempty(gdf)
                gdf(:,2) = (gdf(:,2)+run_idx(1)-1)/fs;
            end
            
            %% Add spikes to structure
            spikes.file(f).hour(h).run_times = run_times;
            spikes.file(f).hour(h).fs = fs;
            spikes.file(f).hour(h).params = params;
            spikes.file(f).hour(h).gdf = gdf;
            spikes.file(f).hour(h).chLabels = chLabels;
            spikes.file(f).hour(h).bipolar_labels = bipolar_labels;
            spikes.file(f).hour(h).bad = bad;
            spikes.file(f).hour(h).run_skip = run_skip;
            spikes.file(f).hour(h).bad_details = bad_details;
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