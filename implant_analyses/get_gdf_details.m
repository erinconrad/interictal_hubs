function out = get_gdf_details(p)

%% Parameters
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
spike_folder = [results_folder,'new_spikes/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
name = pt(p).name;

%% Load spike file
spikes = load([spike_folder,sprintf('%s_spikes.mat',name)]);
spikes = spikes.spikes;

%% Flip things that we think are bad to bad
if clean_blocks
    spikes = clean_missed_bad_blocks(spikes);
end

block_dur = diff(pt(p).ieeg.file(1).block_times(1,:))/3600;
run_dur = diff(spikes.file(1).block(1).run_times)/60;
nfiles = length(spikes.file);

%% Identify files with a change in electrodes
[change,no_change_ever] = find_electrode_change_files(pt,p,only_depth);
nchanges = length(change);

for c = nchanges % just do last one
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
    chLabels = clean_labels_2(spikes.file(change(c).files(2)).block(1).chLabels);
    [~,added_idx] = ismember(added,chLabels);
    [~,unchanged_idx] = ismember(unchanged,chLabels);
    unchanged_labels = chLabels(unchanged_idx);
    added_labels = chLabels(added_idx);


    %% Get locs
    if ~isfield(pt(p).ieeg.file(change(c).files(2)),'locs')
        dist = nan(length(unchanged_labels),1);
    else
        added_locs = pt(p).ieeg.file(change(c).files(2)).locs(added_idx,:);
        unchanged_locs = pt(p).ieeg.file(change(c).files(2)).locs(unchanged_idx,:);


        %% For each unchanged electrode, get identity of and distance from nearest added electrode
        dist = distance_from_closest_added(unchanged_locs,added_locs);

    end

    all_rate = [];
    rate_post = [];
    findices = [];
    bindices = [];
    all_cos = [];
    all_un_cos = [];

    % Loop over files
    for f = 1:last_file

        flocs = pt(p).ieeg.file(f).locs;
        fdist = distance_from_closest_added(flocs,added_locs);

        chLabels = clean_labels_2(spikes.file(f).block(1).chLabels);
        nchs = length(chLabels);
        nblocks = length(spikes.file(f).block);
        rate = nan(nchs,nblocks); % default should be nans
        cos = nan(nchs,nblocks);
        un_cos = nan(nchs,nchs,nblocks);
        
        
        sz_times = all_sz_times_in_file(pt,p,f);


        % Loop over blocks
        for h = 1:nblocks
            block = spikes.file(f).block(h);
            fs = block.fs;
            findices = [findices,f];
            bindices = [bindices,h];

            if block.run_skip == 1
                continue; % leave the whole block as nans
            end


            gdf = spikes.file(f).block(h).gdf;

            %% Spike rate
            if ~isempty(gdf)

                %% Remove duplicates
                if rm_dup 
                    [gdf,~] = remove_duplicates(gdf);
                end

                if isempty(gdf)
                    continue
                end


                %% Remove any spikes in sz
                if rm_sz
                    [gdf,n_removed] = remove_spikes_in_sz(gdf,sz_times);
                end

                for ich = 1:nchs
                    % If the channel was marked as bad or skip, keep it
                    % a nan
                    if ismember(ich,block.bad) ||...
                            ismember(ich,block.skip.all)
                        rate(ich,h) = nan;
                        continue
                    end
                    rate(ich,h) = sum(gdf(:,1)==ich);
                end

                if ~isempty(gdf)
                    %% Get cospike index
                    [cos(:,h),un_cos(:,:,h)] = co_spiking(gdf,fs,added_labels,chLabels);

                end
            else
                for ich = 1:nchs
                    % If the channel was marked as bad or skip, keep it
                    % a nan
                    if ismember(ich,block.bad) ||...
                            ismember(ich,block.skip.all)
                        rate(ich,h) = nan;
                        continue
                    end
                    rate(ich,h) = 0;
                end
            end



        end

        % This part is important to stitch together things properly

        %% Get indices of unchanged and remove changed
        [lia,locb] = ismember(chLabels,unchanged);
        new_labels = chLabels;
        new_labels(~lia) = [];
        rate(~lia,:) =[];
        cos(~lia,:) = [];
        un_cos(~lia,:,:) = [];
        un_cos(:,~lia,:) = [];

        %% Re-order as needed
        [lia,locb] = ismember(unchanged,new_labels);
        new_labels = new_labels(locb);
        rate = rate(locb,:);
        cos = cos(locb,:);
        un_cos = un_cos(lia,:,:);
        un_cos = un_cos(:,lia,:);
        if ~isequal(new_labels,unchanged)
            error('oh no');
        end
        %}

        all_rate = [all_rate,rate];
        all_un_cos = cat(3,all_un_cos,un_cos);


        %% If it's a post-change file, add the coa matrix
        if f >= change(c).files(2)
            all_cos = [all_cos,cos];
            rate_post = [rate_post,rate];

        end


        if f + 1 == change(c).files(2)
            change_block = size(all_rate,2);
        end

    end
    
    
end

out.rate = all_rate;
out.cos = all_cos;
out.un_cos = all_un_cos;
out.change = change;
out.unchanged_labels = unchanged_labels;
out.added_labels = added_labels;
out.dist = dist;
out.unchanged_locs = unchanged_locs;
out.added_locs = added_locs;
out.change_block = change_block;
out.findices = findices;
out.bindices = bindices;
out.fs = fs;
out.block_dur = block_dur;
out.run_dur = run_dur;

end