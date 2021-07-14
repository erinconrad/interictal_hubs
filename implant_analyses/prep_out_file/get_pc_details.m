function metrics = get_pc_details(p)

%{
Think more about how I define connectivity to the new elecs (I kind of like
my way)

Try this out as:
1) A thing that changes after revision (unchanged pc)
2) A predictor of change (unchanged-added pc)
3) A SAR weight matrix
%}

%% Parameters
only_depth = 0;

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
name = pt(p).name;

%% Load pc file
if exist([pc_folder,sprintf('%s_pc.mat',name)],'file') == 0
    metrics = [];
    return;
end

pc = load([pc_folder,sprintf('%s_pc.mat',name)]);
pc = pc.pc;

nfiles = length(pc.file);

%% Identify files with a change in electrodes
[change,no_change_ever] = find_electrode_change_files(pt,p,only_depth);
nchanges = length(change);
c = nchanges;

last_file = nfiles;

added = change(c).added;
unchanged = no_change_ever;%change(c).unchanged;

%% Get total number of blocks
nb = 0;
for f = 1:last_file
    nblocks_in_file = length(pc.file(f).block);
    if strcmp(name,'HUP136') && f == 1
        last_good_block = fix_hup136;
        nblocks_in_file = last_good_block;
    end
    
    nb = nb + nblocks_in_file;
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

%% initialize thingy
net = nan(length(unchanged_labels)*(length(unchanged_labels)-1)/2,nb);
net_post = [];
b_count = 0;
findices = [];
bindices = [];

% Loop over files
for f = 1:last_file

    nblocks = length(pc.file(f).block);
    % fix for hup136
    if strcmp(name,'HUP136') && f == 1
        last_good_block = fix_hup136;
        nblocks = last_good_block;
    end
    flocs = pt(p).ieeg.file(f).locs;

     % Loop over blocks
    for h = 1:nblocks
        block = pc.file(f).block(h);

        run_labels = block.run_labels;
        b_count = b_count + 1;

        if block.run_skip == 1
            continue; % leave the whole block as nans
        end


        pc_w = wrap_or_unwrap_adjacency(block.pc);
        pc_w_added = pc_w;


        %% Remove and pad as needed
        % Remove electrodes that are not unchanged
        [lia,locb] = ismember(run_labels,unchanged);
        new_labels = run_labels;
        new_labels(~lia) = [];
        pc_w(~lia,:) = [];
        pc_w(:,~lia) = [];

        % Pad with missing electrodes (primarily those that we excluded from run due to baddness)
        missing_idx = ~ismember(unchanged,run_labels);
        missing_labels = unchanged(missing_idx);
        pc_w = padarray(pc_w,...
            [length(missing_labels) length(missing_labels)],...
            nan,'post');
        new_labels = [new_labels;missing_labels];

        % Re-order as needed
        [lia,locb] = ismember(unchanged,new_labels);
        new_labels = new_labels(locb);
        pc_w = pc_w(locb,:);
        pc_w = pc_w(:,locb);

        if ~isequal(new_labels,unchanged)
            error('ruh roh');
        end
        
        %% Do stuff if it's a post-change file
        if f >= change(c).files(2)
            
            %% Also do matching for the matrix with unchanged and added
            
            % Remove
            lia = ismember(run_labels,unchanged);
            lib = ismember(run_labels,added);
            new_un_labels = run_labels;
            new_ad_labels = run_labels;
            new_un_labels(~lia) = [];
            new_ad_labels(~lib) = [];
            pc_w_added(~lia,:) = [];
            pc_w_added(:,~lib) = [];

            % Pad
            missing_un = ~ismember(unchanged,run_labels);
            missing_un_labels = unchanged(missing_un);
            missing_ad = ~ismember(added,run_labels);
            missing_ad_labels = added(missing_ad);
            pc_w_added = padarray(pc_w_added,...
                [sum(missing_un) sum(missing_ad)],...
                nan,'post');
            new_un_labels = [new_un_labels;missing_un_labels];
            new_ad_labels = [new_ad_labels;missing_ad_labels];

            % Reorder
            [~,loca] = ismember(unchanged,new_un_labels);
            new_un_labels = new_un_labels(loca);
            pc_w_added = pc_w_added(loca,:);

            [~,locb] = ismember(added,new_ad_labels);
            new_ad_labels = new_ad_labels(locb);
            pc_w_added = pc_w_added(:,locb);

            if ~isequal(new_un_labels,unchanged) || ~isequal(new_ad_labels,added)
                error('ruh roh');
            end
        
            %pc_uw_added = wrap_or_unwrap_adjacency(pc_w_added);
            
            if 0
                figure
                imagesc(pc_w_added)
                xticks(1:length(added))
                yticks(1:length(unchanged))
                xticklabels(added)
                yticklabels(unchanged)
            end
            
            net_post = cat(3,net_post,pc_w_added);
        
        end

        
        %% Fill up network array
        pc_uw =  wrap_or_unwrap_adjacency(pc_w);
        net(:,b_count) = pc_uw;

        if 0
            figure
            turn_nans_white(pc_w)
            xticks(1:size(pc_w,1))
            yticks(1:size(pc_w,1))
            xticklabels(new_labels)
            yticklabels(new_labels)
            colorbar
            title(sprintf('File %d block %d',f,h))
            pause
            close(gcf)
        end

    end

    if f + 1 == change(c).files(2)
        change_block = b_count;
    end

end


%% Get network metrics
metrics = get_network_metrics(net,unchanged_labels,1);
%metrics_old = get_network_metrics(net,unchanged_labels,0);

if 0
    fprintf('\n\n%s\n',name);
    r = corr(nanmean(metrics_rm_ekg.ns,2),nanmean(metrics_old.ns,2),'rows','pairwise')
end


if 0
    imagesc(metrics.avg_mat)
    xticks(1:size(metrics.avg_mat,1))
    yticks(1:size(metrics.avg_mat,1))
    xticklabels(unchanged_labels)
    yticklabels(unchanged_labels)
end

%% Get "ns" across added channels
mean_net_post = nanmean(nanmean(net_post,3),2);
metrics.added_pc = mean_net_post;

if 0
    figure
    imagesc(nanmean(net_post,3))
    xticks(1:length(added))
    yticks(1:length(unchanged))
    xticklabels(added)
    yticklabels(unchanged)
end

end