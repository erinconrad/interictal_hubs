function measure_pause

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
timing_folder = [locations.main_folder,'data/'];
main_spike_results = [results_folder,'main_spikes/'];
if ~exist(main_spike_results,'dir')
    mkdir(main_spike_results);
end

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];

ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;


%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load file with timing
timing_file = [timing_folder,'reimplantation patients.xlsx'];
T = readtable(timing_file,'Sheet','pause details','ReadRowNames',true);

%% initialize arrays
prefile_duration = cell(size(T,1),1);
prefile_end = cell(size(T,1),1);
gap = cell(size(T,1),1);
all_num_sz = cell(size(T,1),1);

%% Get prefile duration and gap
% Loop over rows
for t = 1:size(T,1)
    
    % Get details
    name = T.Properties.RowNames{t};
    pre_start = datetime(T.PreFileStart(t));
    post_start = datetime(T.PostFileStart(t));
    f = T.Prefile(t);
    
    if isempty(T.PreFileStart{t})
        continue
    else
    
        % find corresponding pt index in pt struct
        for p = 1:length(pt)
            if strcmp(pt(p).name,name)
                % found it

                % get filename
                fname = pt(p).ieeg.file(f).name;

                % get duration
                dur = pt(p).ieeg.file(f).duration/1e6; % in microseconds, need to divide
                
                
                % Get pre-file end
                pre_end = pre_start + seconds(dur);
                
                % get gap
                pre_post_gap = post_start - pre_end;
                
                % fill up arrays
                prefile_duration{t} = dur;
                prefile_end{t} = pre_end;
                gap{t} = pre_post_gap;
                
                
                % get num szs
                num_sz = 0;
                for k = 1:length(pt(p).ieeg.file)
                    nsz = size(pt(p).ieeg.file(k).sz,1);
                    num_sz = num_sz + nsz;
                end
                all_num_sz{t} = num_sz;

                % break out of pt loop
                break

            end
        end
    end
    
end

%% Write this as columns in the table
new_col_T = cell2table([prefile_duration,prefile_end,gap],...
     'VariableNames',{'Prefile duration','Prefile end','Gap'});
writetable(new_col_T,timing_file,'Sheet','pause details','Range','E1:G17');

end