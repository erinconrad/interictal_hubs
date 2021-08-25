function calculate_pause_plus_disconnect(out)

%% Temporarily add this codebase to the Matlab path
locations = interictal_hub_locations; % This is a script you will need to add to your own path (see ReadMe for details)
addpath(genpath(locations.script_folder)); % path to the interictal_hubs codebase
whichPts = [20 103 106 107 35 109 110 111 94 97];
surround = 24; % doesn't matter

%% Get individual patient data
all_pauses = nan(length(whichPts),1);
names = cell(length(whichPts),1);
for i = 1:length(whichPts)
    rate = out(i).rate;
    cblock = out(i).change_block;
            
    % Get file gap
    [~,buffer] = file_gaps(out(i).name);

     % Get surround times, starting with first non nan
    [~,~,pre_nans,post_nans] = get_surround_times(rate,cblock,surround);
    
    total_pause_blocks = buffer+(pre_nans+post_nans)/2; % convert nan blocks to hours
    pause_hours = total_pause_blocks;
    all_pauses(i) = pause_hours;
    names{i} = out(i).name;

end

table(names,all_pauses)

end