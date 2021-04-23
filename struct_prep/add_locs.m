%{
Note that this makes it so that the locs in patient's ieeg file correspond
exactly to the channel labels marked for that file

%}

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.script_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% load loc struct
loc = load([data_folder,'patient_localization.mat']);
loc = loc.patient_localization;

whichPts = 1:length(pt);

% Loop over patients
for i = 1:length(whichPts)
    p = whichPts(i);
    name = pt(p).name;
    
    if isempty(pt(p).ieeg), continue; end
    
    %% Get corresponding patient from localization thing
    pt_match = 0;
    for j = 1:length(loc)
        if strcmp(name,loc(j).patient)
            pt_match = 1;
            break
        end
    end
    if pt_match == 0
        fprintf('\nCannot find %s in localization file\n',name);
        continue;
        
    end
    
    match = 0;
    % Loop over files
    for f = 1:length(pt(p).ieeg.file)
        file = pt(p).ieeg.file(f);
        labels = clean_labels_2(file.chLabels(:,1));
        loc_labels = clean_labels_2(loc(j).labels);
        loc_coords = loc(j).coords;
        
        
        % see if labels are equal
        if isequal(labels,loc_labels)
            pt(p).ieeg.file(f).locs = loc_coords;
            match = 1;
        else
            % see if they are equal once I turn "RG" or "LG" into GRID
            temp_labels = convert_g_to_grid_label(labels,name);
            temp_loc_labels = convert_g_to_grid_label(loc_labels,name);
            if isequal(temp_labels,temp_loc_labels)
                pt(p).ieeg.file(f).locs = loc_coords;
                match = 1;
            end
        end
    end
    
    if match == 0
        error('electrode label mismatch');
    end
    
end

save([data_folder,'pt.mat'],'pt')
    