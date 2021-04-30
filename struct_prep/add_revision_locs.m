clear

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.script_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% load revision struct
re = load([data_folder,'reimplant_pts.mat']);
re = re.pt;

% Loop over revision patients
for r = 1:length(re)
    name = re(r).name;
    
    if ~isstruct(re(r).master_elecs)
        fprintf('\nWarning, no locs for %s any channel\n',name);
        continue
    end
    
    re_labels = clean_labels_2(re(r).master_elecs.master_labels);
    
    % Get re-locs in a non ridiculous format
    re_locs = nan(length(re_labels),3);
    
    if ~isfield(re(r).master_elecs,'locs')
        fprintf('\nWarning, no locs for %s any channel\n',name);
        continue
    end
    
    for e = 1:length(re(r).master_elecs.locs)
        
        if isempty(re(r).master_elecs.locs(e).system)
            fprintf('\nWarning, no locs for %s channel %d\n',name,e);
            continue
        end
        
        if ~isempty(re(r).master_elecs.locs(e).system(1).locs)
            re_locs(e,:) = re(r).master_elecs.locs(e).system(1).locs;
        elseif ~isempty(re(r).master_elecs.locs(e).system(2).locs)
            re_locs(e,:) = re(r).master_elecs.locs(e).system(2).locs;
        else
            fprintf('\nWarning, no locs for %s channel %d\n',name,e);
            continue
        end
    end
    
    pt_match = 0;
    % Loop over main pt struct
    for p = 1:length(pt)
        
        if strcmp(name,pt(p).name)
            pt_match = 1;
            break
        end
        
    end
    
    if pt_match == 0
        fprintf('\nNo corresponding pt\n');
        continue;
    end
    
    % Loop over files
    for f = 1:length(pt(p).ieeg.file)
        f_labels = clean_labels_2(pt(p).ieeg.file(f).chLabels(:,1));
        
        % Find the matching labels from the master list
        [~,locb] = ismember(f_labels,re_labels);
        
        matching_locs = re_locs(locb,:);
        pt(p).ieeg.file(f).locs = matching_locs;
        
    end
    
end


save([data_folder,'pt.mat'],'pt');