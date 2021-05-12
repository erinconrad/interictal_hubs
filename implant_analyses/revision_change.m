function revision_change(whichPts)

%% Parameters
filt = 2;
timing = 'peak';
only_depth = 1;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
ieeg_folder = locations.ieeg_folder;
addpath(genpath(ieeg_folder));
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;
addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
addpath(genpath(locations.ieeg_folder));
spike_folder = [results_folder,'new_spikes/'];
bct_folder = locations.bct;
addpath(genpath(bct_folder));

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

if isempty(whichPts)
    whichPts = [20 103 105 106 107 35 109 110 111 94 97];
end

%% Prep info for aggregate tables
all_names = {};
all_big_inc = {};
all_anat = {};
all_inc = [];

for p = whichPts
    pt_name = pt(p).name;

    %% Load spike file
    if exist([spike_folder,sprintf('%s_spikes.mat',pt_name)],'file') == 0
        continue;
    end
      
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
    name = spikes.name;
    block_dur = diff(pt(p).ieeg.file(1).block_times(1,:))/3600;
    run_dur = diff(spikes.file(1).block(1).run_times)/60;
    nfiles = length(spikes.file);
    
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
        chLabels = clean_labels_2(spikes.file(change(c).files(2)).block(1).chLabels);
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


            %% Show this
            if 0
                table(unchanged_labels,added_labels(closest_added),dist)

            end
        end
        
        %% Get anatomy
        if ~isfield(pt(p).ieeg.file(change(c).files(2)),'anatomy')
            unchanged_anatomy = cell(length(unchanged_labels),1);
        else
            unchanged_anatomy = pt(p).ieeg.file(change(c).files(2)).anatomy(unchanged_idx);
        end
        
        all_rate = [];
        coa_post = [];
        rate_post = [];
        all_nseq = [];
        all_seq = {};
        last_block = zeros(last_file-1,1);
        findices = [];
        bindices = [];
        
        % Loop over files
        for f = 1:last_file
            
            chLabels = clean_labels_2(spikes.file(f).block(1).chLabels);
            nchs = length(chLabels);
            nblocks = length(spikes.file(f).block);
            rate = nan(nchs,nblocks); % default should be nans
            coa = nan(nchs,nchs,nblocks);
            num_seq = nan(nchs,nblocks);
            
            
            
            % Loop over blocks
            for h = 1:nblocks
                block = spikes.file(f).block(h);
                
                findices = [findices,f];
                bindices = [bindices,h];

                if block.run_skip == 1
                    continue; % leave the whole block as nans
                end

                %% Get spike channels and times (whatever time you want)
                if isempty(block.details)
                    chs = [];
                    times = [];
                else
                    chs = block.details.filter(filt).gdf(:,1);
                    times = block.details.filter(filt).(timing);
                end
                fs = block.fs;

                %% Get sorted spike indices and chs
                [times,I] = sort(times);
                chs = chs(I);

                %% Construct gdf
                gdf = [chs,times];
                
                %% Spike rate
                if ~isempty(gdf)
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
                    
                     %% Get sequences, rl, coa
                    [seq,~,coa(:,:,h),num_seq(:,h)] = new_get_sequences(gdf,nchs,fs);
                    
                    if f >= change(c).files(2)
                        all_seq = [all_seq;seq];
                    end
                else
                    rate(:,h) = 0;
                end
                
                
                
            end
            
            % This part is important to stitch together things properly

            %% Get indices of unchanged and remove changed
            [lia,locb] = ismember(chLabels,unchanged);
            new_labels = chLabels;
            new_labels(~lia) = [];
            rate(~lia,:) =[];

            %% Re-order as needed
            %[lia,locb] = ismember(new_labels,unchanged);
            [lia,locb] = ismember(unchanged,new_labels);
            new_labels = new_labels(locb);
            rate = rate(locb,:);
            if ~isequal(new_labels,unchanged)
                error('oh no');
            end
            %}
            
            all_rate = [all_rate,rate];
            
            %% If it's a post-change file, add the coa matrix
            if f >= change(c).files(2)
                coa_post = cat(3,coa_post,coa);
                rate_post = [rate_post,rate];
                all_nseq = [all_nseq,num_seq];
                post_labels = chLabels;
            end
            
            if f<last_file
                last_block(f) = size(all_rate,2);
            end
            
            if f + 1 == change(c).files(2)
                change_block = size(all_rate,2);
            end
            
        end
        
        
    %% Text to designate new chs
    lia = ismember(post_labels,added);
    new_post_labels = post_labels;
    for j = 1:length(post_labels)
        if lia(j) == 1
            new_post_labels{j} = [new_post_labels{j},'***'];
        end
    end
        
    %% Perc before - perc after
    perc_diff = comp_sequences(all_seq,post_labels,added,post_labels);
    if 0
        figure
        imagesc(perc_diff)
        xticks(1:length(post_labels))
        yticks(1:length(post_labels))
        xticklabels(new_post_labels)
        xtickangle(90)
        yticklabels(new_post_labels)
        colorbar
    end
    
    
    
    if 0
        figure
        turn_nans_white(nanmean(all_coa(:,:,:),3))
        xticks(1:length(post_labels))
        yticks(1:length(post_labels))
        xticklabels(new_post_labels)
        xtickangle(90)
        yticklabels(new_post_labels)
    end
    
    if 0
        figure
        set(gcf,'position',[1 1 1400 800])
        turn_nans_white(all_rate)
        hold on
        for b = 1:length(last_block)
            plot([last_block(b) last_block(b)],ylim,'k','linewidth',3)
        end
        plot([change_block change_block],ylim,'r','linewidth',3)
        yticks(1:length(unchanged))
        yticklabels(dist_info)
        %title(added')
        
    end
    
    %% Overall change in spike rate
    if 0
        show_overall_rate(all_rate,block_dur,last_block,change_block,run_dur,name,results_folder)
    end
    
    %% clustered
    if 0
        if strcmp(name,'HUP128')
            clusts{1} = find(contains(unchanged,'L')); clusts{2} = find(contains(unchanged,'R'));
            clust_names = {'Left','Right'};
            show_clusters_rate(all_rate,block_dur,change_block,run_dur,name,results_folder,clusts,clust_names)
        end
    end
    
    %% Get rate increase of electrodes with minimum spike rate
    
    [rate_increase,spikey_labels,spikey_idx,mean_rate_post,...
        big_inc_labels,big_inc_rate_sorted,abs_increase] = find_rate_increase(all_rate,change_block,unchanged);
    % Get anatomy of the big increase electrodes
    [~,big_inc_chs] = ismember(big_inc_labels,unchanged_labels);
    big_inc_anatomy = unchanged_anatomy(big_inc_chs);
    
    all_big_inc = [all_big_inc;big_inc_labels];
    all_anat = [all_anat;big_inc_anatomy];
    all_names = [all_names;repmat(cellstr(pt_name),length(big_inc_anatomy),1)];
    all_inc = [all_inc; big_inc_rate_sorted];
    
    
    %% Are electrodes with bigger spike rate increase closer to the new electrodes than are other spikey electrodes?
    if 0
        
        % Get distance from closest new electrode of these spikey electrodes
        dist_spikey = dist(spikey_idx);
        
        
        plot_inc_as_fcn_of_dist(rate_increase,abs_increase,dist_spikey,...
            spikey_labels,name,results_folder,run_dur,big_inc_labels)

    end
    
    %% Are eletrodes with bigger spike rate increase more likely to co-spike with new electrodes?
    if 0
        blocks = 1:100;
        [cos,unchanged_spikey_labels] = ...
            co_spike_index(rate_post,spikey_idx,coa_post,...
            blocks,added,unchanged,post_labels,new_post_labels,mean_rate_post,...
            abs_increase,pt_name,run_dur,results_folder);
        %}

        
    
    end
    
    if 0
        raster_rate_chs(all_rate,last_block,block_dur,run_dur,change_block,...
    dist_info,unchanged,name,results_folder)
        
    end
    
    
    
    end
    
end

if 0
%% Show table with anatomical localizations
outfolder = [results_folder,'anatomy/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end
T = table(all_names,all_big_inc,all_inc,all_anat);
writetable(T,[outfolder,'anatomy.csv']);
end



end

