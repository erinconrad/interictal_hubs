function alt_spikes(whichPts)

%% Parameters
overwrite = 0;
test.do_test = 0;
do_plot = 0;
do_save = 1;

%% Test parameters
test.pt = 20;
test.time = 29800;
test.dur = 60;
test.file = 1;
test.ch =[];

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
sp_folder = [results_folder,'alt/spikes/'];

if exist(sp_folder,'dir') == 0
    mkdir(sp_folder);
end

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load localization structure
loc = load([data_folder,'patient_localization.mat']);
loc = loc.patient_localization;

%% Get which patients to run
if isempty(whichPts)
    % The implant revision patients
    whichPts = [20 103 104 105 106 107 108,...
        35 109 110 78 111 90 112 94 97];%1:length(pt);
end

% Loop over patients
for i = 1:length(whichPts)
    
    if test.do_test == 1
        p = test.pt;
        name = pt(p).name;
        next_file = 1;
        next_block = 1;
    else
        p = whichPts(i);
        name = pt(p).name;

        out_file = [name,'_spikes.mat'];
        
        %% Throw warning if there is missing seizure data
        all_missed = check_for_missing_szs(pt,p);
        if all_missed == 1, error('Missing seizure data, make sure you add this'); end

        %% Load the existing spike structure if it exists
        if overwrite == 0 && exist([sp_folder,out_file],'file') ~= 0
            spikes = load([sp_folder,out_file]);
            spikes = spikes.spikes;

            % find the last block we have finished
            next_file = spikes.next_file;
            next_block = spikes.next_block;

            if isnan(next_file)
                fprintf('\nAlready finished %s, skipping...\n',name);
                continue
            end
        else
            % initialize it
            clear spikes % I must clear this 
            next_file = 1;
            next_block = 1;
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
        orig_labels = clean_labs;
        
        %% Reconcile cleaned labels with cleaned loc labels
        % For the purpose of knowing which electrodes to skip
        if ~isequal(clean_labs,clean_loc_labs)
            fprintf('\nWarning, localization labels do not match ieeg labels\n');
        end
        
        
        %% Designate which electrodes to skip
        if ~strcmp(name,'HUP132')
            [which_chs,skip] = designate_chs(chLabels,clean_labs,clean_loc_labs,loc(loc_p));
            non_skip = which_chs;
        else
            which_chs = 1:length(chLabels); non_skip = which_chs; skip = [];
        end
        
        % filename
        fname = pt(p).ieeg.file(f).name;
        
        % Loop over blocks
        n_blocks = length(pt(p).ieeg.file(f).block);
        for h = next_block:n_blocks
            tic;
            fprintf('\nDoing %s file %d block %d of %d\n',name,f,h,n_blocks);
            
            % get the run time (already randomly assigned minute)
            if test.do_test == 1
                run_times = [test.time,test.time+test.dur];
            else
                run_times = pt(p).ieeg.file(f).block(h).run; 
            
            end
            if ~isempty(run_times)
             
                run_idx = run_times(1)*fs:run_times(2)*fs;
                dur = diff(run_times);

                %% Get the eeg data
                values = pull_ieeg_data(fname, login_name, pwfile, run_idx);
                
                %% Re-name chLabels for HUP132
                if strcmp(name,'HUP132')
                    clean_labs = fix_hup132(f,run_times,orig_labels,data_folder);
                    chLabels = clean_labs;
                end
                                
                %% Designate electrodes over which to run spike detector
                if test.do_test == 1 && ~isempty(test.ch)
                    if ischar(test.ch)
                        test.ch = {test.ch};
                    end
                    if iscell(test.ch)
                        which_chs = (cellfun(@(x) find(strcmp(clean_labs,x)),test.ch));
                    else
                        which_chs = test.ch;
                    end
                else

                    % I have already designated channels to run about at the
                    % file level
                end
                
                %% Find non-intracranial chs
                non_intracranial = find_non_intracranial(clean_labs);
                which_chs = find(~non_intracranial); % channels to do analysis on
                intracranial_chs = which_chs;

                %% Reject bad channels
                [bad,bad_details] = identify_bad_chs(values,which_chs,chLabels,fs);
                which_chs(ismember(which_chs,bad)) = [];
                
                %% Do pre-processing
                [car_values,car_labels] = car_montage(values,which_chs,clean_labs);
         
                is_run_car = ismember((1:length(clean_labs))',which_chs);

                %% Skip the chunk entirely if enough channels are bad (suggests period of disconnection)
                if length(bad) >= 0.5*length(intracranial_chs)
                    fprintf('\nSkipping this run because %d of %d chs marked bad\n',length(bad),length(intracranial_chs));
                    gdf = [];
                    which_chs = [];
                    run_skip = 1;
                    details = [];
                else
                    run_skip = 0;
                    %% Spike detector
                    
                    if ~isempty(which_chs)
                        % filters
                        car_values = notch_filter(car_values,fs);
                        car_values = bandpass_filter(car_values,fs);
                        
                        run_car_values  = car_values;
                        run_car_values(:,~is_run_car) = nan;
                        gdf = detector_alt(run_car_values,fs);
                    else
                        gdf = [];
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
                    ex_chs = [];
                    only_run = 0;

                    simple_plot(car_values,clean_labs,is_run_car,...
                        fs,ex_chs,gdf,...
                        only_run,bad,clean_labs)

                    pause
                    close(all)
                    %}
                end


                %% Re-align gdf time
                if ~isempty(gdf)
                    gdf(:,2) = (gdf(:,2)+run_idx(1)-1)/fs;
                end
                
            else
                gdf = [];
                bad = [];
                run_skip = 1;
                car_labels = [];
                bad_details =  [];
                skip = [];
                which_chs = [];
            end
            
            %% Add spikes to structure
            spikes.file(f).block(h).run_times = run_times;
            spikes.file(f).block(h).fs = fs;
            spikes.file(f).block(h).gdf = gdf;
            spikes.file(f).block(h).chLabels = chLabels;
            spikes.file(f).block(h).car_labels = car_labels;
            spikes.file(f).block(h).bad = bad;
            spikes.file(f).block(h).run_skip = run_skip;
            spikes.file(f).block(h).bad_details = bad_details;
            spikes.file(f).block(h).skip = skip;
            spikes.file(f).block(h).which_chs = which_chs;
            
            %% Mark the next block
            if h == n_blocks
                next_block = 1;
                next_file = f + 1;
                if f == n_files
                    next_file = nan;
                end
            else
                next_block = h + 1;     
            end
            spikes.next_file = next_file;
            spikes.next_block = next_block;
            
            %% Save structure
            if test.do_test == 0 && do_save == 1
                save([sp_folder,out_file],'spikes');
            end
           
            
        end
    end
    
end



end