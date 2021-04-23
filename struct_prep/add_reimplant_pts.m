function add_reimplant_pts

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.script_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;



%% Load reimplantation struct
re = load([data_folder,'reimplant_pts.mat']);
re = re.pt;

next_p = length(pt) + 1;
% Loop through reimplantation patients
for i = 1:length(re)
    
    re_name = re(i).name;
    
    
    % see if it already exists
    done = 0;
    for j = 1:length(pt)
        pt_name = pt(j).name;
        if strcmp(pt_name,re_name)
            fprintf('\n%s already exists\n',pt_name);
            done = 1;
            break
        end
    end
    
    %% if doesn't exist, add stuff about it
    if done == 0
        pt(next_p).name = re_name;
        pt(next_p).clinical = [];
        
        pt(next_p).electrode_info.resected = 'missing';
        pt(next_p).electrode_info.ignore = [];
        
        %% Add ieeg stuff
        for f = 1:length(re(i).filename)
            pt(next_p).ieeg.file(f).fs = re(i).filename(f).fs;
            pt(next_p).ieeg.file(f).name = re(i).filename(f).ieeg_name;
            pt(next_p).ieeg.file(f).chLabels = re(i).filename(f).chLabels;
            pt(next_p).ieeg.file(f).duration = re(i).filename(f).duration*1e6; % convert to microseconds            
        end
        
        %% Add seizure time data
        % Need to do!!!
        
    end
    
    %% Add one to next index
    next_p = next_p + 1;
end


end