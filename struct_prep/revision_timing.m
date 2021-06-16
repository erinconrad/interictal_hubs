function revision_timing(whichPts)

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


%% Get which patients to run
if isempty(whichPts)
    whichPts = [20 103 106 107 35 109 110 111 94 97];
end

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

all_names = {};
all_gaps = [];

for i = [1:5,7:8]
    
    p = whichPts(i);

    name = pt(p).name;
    
    %% Get change info
    [change,no_change_ever] = find_electrode_change_files(pt,p,0);
    nchanges = length(change);
    pre_file = change(nchanges).files(1);

    %% This gives the metadata start times for the pre- and post-revision files
    switch name
        case 'HUP099'
            pre_start = '2000-01-04 16:40:51';
            post_start = '2000-01-07 11:00:56';
        case 'HUP100'   % might be wrong and not sure how to check; I guessed based on Natus EEG start times
            pre_start = '2000-01-01 16:10:51';
            post_start = '2000-01-08 14:57:26';
        case 'HUP111'
            pre_start = '2000-01-02 22:06:08';
            post_start = '2000-01-07 19:37:38';
        case 'HUP128'
            pre_start = '2000-01-01 13:23:04';
            post_start = '2000-01-11 13:40:01';
        case 'HUP132'
            pre_start = '2000-01-03 15:43:57';
            post_start = '2000-01-07 14:51:20';
        case 'HUP136'
            pre_start = '2000-01-01 12:15:46';
            post_start = '2000-01-09 06:45:31';
        case 'HUP152'
            pre_start = '2000-01-01 13:17:47';
            post_start = '2000-01-10 15:44:30';
        case 'HUP193'
            pre_start = '2000-01-05 11:59:16';
            post_start = '2000-01-09 10:44:24';
        case 'HUP201'
        case 'HUP209'

    end

    %% Calculates the difference in start times
    pre_t = datetime(pre_start);
    post_t = datetime(post_start);
    metadata_diff = seconds(post_t-pre_t); % get it in seconds

    %% Get the ieeg duration of the pre-revision file
    ieeg_dur = pt(p).ieeg.file(pre_file).duration/(1e6); %note that duration is in microseconds, convert to seconds
    
    %% The difference in metadata_diff and ieeg_dur is the gap in ieeg time
    gap = metadata_diff - ieeg_dur;
    if gap < 0, error('why is gap negative?'); end
    gap = gap/3600; % Convert to hours
    
    all_names = [all_names;name];
    all_gaps = [all_gaps;gap];

end

table(all_names,all_gaps)

end