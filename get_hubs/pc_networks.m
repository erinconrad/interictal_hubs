function pc_networks(whichPts)

%% Parameters
overwrite = 0;
test.do_test = 0;
do_plot = 0;
do_save = 1;

tw = 2; % 2 second calculations

%% Test parameters
test.pt = 20;%111;
test.time = 51384.29;%11725;
test.dur = 15;
test.file = 3;
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
out_folder = [results_folder,'pc/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load localization structure
loc = load([data_folder,'patient_localization.mat']);
loc = loc.patient_localization;

%% Get which patients to run
if isempty(whichPts)
    whichPts = [20 103 106 107 35 109 110 111 94 97];
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
        
        out_file = [name,'_pc.mat'];
        
        %% Throw warning if there is missing seizure data
        all_missed = check_for_missing_szs(pt,p);
        if all_missed == 1, fprintf('Missing seizure data, make sure you add this'); end

        %% Load the existing pc structure if it exists
        if overwrite == 0 && exist([out_folder,out_file],'file') ~= 0
            pc = load([out_folder,out_file]);
            pc = pc.pc;

            % find the last block we have finished
            next_file = pc.next_file;
            next_block = pc.next_block;

            if isnan(next_file)
                fprintf('\nAlready finished %s, skipping...\n',name);
                continue
            end
        else
            % initialize it
            clear pc % I must clear this 
            next_file = 1;
            next_block = 1;
            pc.name = name;

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

                %% Do pre-processing
                [bipolar_values,bipolar_labels,chs_in_bipolar] = pre_process(values,clean_labs);

                %% Designate electrodes over which to run
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
                
                %% Reject bad channels
                [bad,bad_details] = reject_bad_chs(bipolar_values,which_chs,chLabels,fs);
                
                if length(bad) >= 0.5*length(which_chs)
                    fprintf('\nSkipping this run because %d of %d chs marked bad\n',length(bad),length(which_chs));
                    run_chs = [];
                    run_skip = 1;
                    details = [];
                    pc_out = [];
                    run_labels = {};
                else
                    run_skip = 0;
                    
                    %% Do CAR on good chs
                    %{
                    I am adding something whereby I define the CAR to ONLY
                    include those electrodes that are pre-existing. This is
                    to avoid contaminating the pre-existing electrodes with
                    new electrode signal, which might falsely cause me to
                    see a change post-implantation (which is just the
                    signal from the new electrodes)
                    %}
                    run_chs = which_chs;
                    run_chs(ismember(run_chs,bad)) = [];
                    values = values(:,run_chs);
                    
                    % Get the labels of the pre-existing channels
                    [~,pre_existing_labels] = find_electrode_change_files(pt,p,0);
                    
                    % labels of channels I am currently planning to run it
                    % over
                    curr_labels = clean_labs(run_chs);
                    
                    % find those that match pre-existing labels 
                    pre_existing_idx = find(ismember(curr_labels,pre_existing_labels));
                    
                    old_values = values;
                    values = new_pre_process(values,pre_existing_idx); % define car according to preexisting only
                    %values = new_pre_process(values,1:size(values,2));
                    
                    %% Do filters
                    values = do_filters(values,fs);
                    
                    %% Do PC
                    pc_out = calc_pc(values,fs,tw);
                    run_labels = clean_labs(run_chs);
                    
                end
                
                %% Run details
                t = toc;
                fprintf('\nTook %1.1f s\n',t);
                fprintf('Of %d non-skipped chs, rejected %d for nans, %d for zeros,\n%d for variance, %d for noise, %d for std.\n',...
                    length(non_skip),length(bad_details.nans),length(bad_details.zeros),length(bad_details.var),...
                    length(bad_details.noisy),length(bad_details.higher_std));
                
                %% Example plot              
                if do_plot
                    pc_uw = wrap_or_unwrap_adjacency(pc_out);
                    figure
                    imagesc(pc_uw)
                    xticks(1:size(pc_uw,2));
                    yticks(1:size(pc_uw,2));
                    yticklabels(run_labels)
                    xticklabels(run_labels);
                    colorbar
                    pause
                    close(gcf)
                    
                end
                
                if 0
                    show_eeg_and_spikes(values,clean_labs,[],dur,run_times(1),name,fs,bad,skip,[]);
                end
            else
                
                bad = [];
                run_skip = 1;
                bad_details = [];
                skip = [];
                run_chs = [];
                run_labels = {};
                details = [];
                
            end
            
            %% Add pc matrix to structure
            pc.file(f).block(h).run_times = run_times;
            pc.file(f).block(h).fs = fs;
            pc.file(f).block(h).chLabels = chLabels;
            pc.file(f).block(h).bad = bad;
            pc.file(f).block(h).run_skip = run_skip;
            pc.file(f).block(h).bad_details = bad_details;
            pc.file(f).block(h).skip = skip;
            pc.file(f).block(h).run_chs = run_chs;
            pc.file(f).block(h).run_labels = run_labels;
            pc.file(f).block(h).pc = pc_out;
            
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
            pc.next_file = next_file;
            pc.next_block = next_block;
            
            %% Save structure
            if test.do_test == 0 && do_save == 1
                save([out_folder,out_file],'pc');
            end
            
            
        end
        
    end

end