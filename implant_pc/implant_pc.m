function implant_pc(whichPts)

% probably should do
rm_sz = 1;
rm_dup = 1;
only_depth = 0;
clean_blocks = 1;


%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
pc_folder = [results_folder,'pc/'];
bct_folder = locations.bct;
addpath(genpath(bct_folder));

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

if ischar(whichPts)
    all_pt_names = cell(length(pt),1);
    for i = 1:length(pt)
        all_pt_names{i}=pt(i).name;
    end
    whichPts = find(strcmp(all_pt_names,whichPts));
end

if isempty(whichPts)
    whichPts = [20 103 106 107 35 109 110 111 94 97];
end

for p = whichPts
    name = pt(p).name;
    
    %% Load pc file
    if exist([pc_folder,sprintf('%s_pc.mat',name)],'file') == 0
        continue;
    end
      
    pc = load([pc_folder,sprintf('%s_pc.mat',name)]);
    pc = pc.pc;
    
    block_dur = diff(pt(p).ieeg.file(1).block_times(1,:))/3600;
    run_dur = diff(pc.file(1).block(1).run_times)/60;
    nfiles = length(pc.file);
    
    %% Identify files with a change in electrodes
    [change,no_change_ever] = find_electrode_change_files(pt,p,only_depth);
    nchanges = length(change);
    
    for c = nchanges
        if c < nchanges
            last_file = change(c+1).files(2)-1;
        else
            last_file = nfiles;
        end
        added = change(c).added;
        unchanged = no_change_ever;%change(c).unchanged;
        
        if isempty(added)
            continue
        end
    
        %% Ch labels
        chLabels = clean_labels_2(pc.file(change(c).files(2)).block(1).chLabels);
        [~,added_idx] = ismember(added,chLabels);
        [~,unchanged_idx] = ismember(unchanged,chLabels);
        unchanged_labels = chLabels(unchanged_idx);
        added_labels = chLabels(added_idx);
        
        %% Get locs
        if ~isfield(pt(p).ieeg.file(change(c).files(2)),'locs')
            dist_info = unchanged_labels;
            dist = nan(length(unchanged_labels),1);
        else
            added_locs = pt(p).ieeg.file(change(c).files(2)).locs(added_idx,:);
            unchanged_locs = pt(p).ieeg.file(change(c).files(2)).locs(unchanged_idx,:);
            

            %% For each unchanged electrode, get identity of and distance from nearest added electrode
            [dist,closest_added] = distance_from_closest_added(unchanged_locs,added_locs);
            closest_label = added_labels(closest_added);
            dist_info = cellfun(@(x,y,z) sprintf('%s %s %1.1f',x,y,z),unchanged_labels,closest_label,num2cell(dist),'UniformOutput',false);

        end
        
        %% Get anatomy
        if ~isfield(pt(p).ieeg.file(change(c).files(2)),'anatomy')
            unchanged_anatomy = cell(length(unchanged_labels),1);
            added_anatomy = cell(length(added_labels),1);
        else
            unchanged_anatomy = pt(p).ieeg.file(change(c).files(2)).anatomy(unchanged_idx);
            added_anatomy = pt(p).ieeg.file(change(c).files(2)).anatomy(added_idx);
        end
        
        findices = [];
        bindices = [];
        
        % Loop over files
        for f = 1:last_file
            
            flocs = pt(p).ieeg.file(f).locs;
            fdist = distance_from_closest_added(flocs,added_locs);
            
            chLabels = clean_labels_2(pc.file(f).block(1).chLabels);
           
             % Loop over blocks
            for h = 1:nblocks
                block = pc.file(f).block(h);
                fs = block.fs;
                findices = [findices,f];
                bindices = [bindices,h];
                run_labels = block.run_labels;
                
                if block.run_skip == 1
                    continue; % leave the whole block as nans
                end

                
                pc_uw = wrap_or_unwrap_adjacency(block.pc);
                
                
                if 1
                    figure
                    imagesc(pc_uw)
                    xticks(1:size(pc_uw,1))
                    yticks(1:size(pc_uw,1))
                    xticklabels(run_labels)
                    yticklabels(run_labels)
                    colorbar
                    pause
                    close(gcf)
                end
                
            end
            
        end
        
    
end


end