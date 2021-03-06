function revision_change(whichPts)

% for permutation test, restrict to post implant

%% Parameters
surround = 48; % Divide by 2 to get number of hours

% probably should do
rm_sz = 1;
rm_dup = 1;
only_depth = 0;
clean_blocks = 1;

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

%% Prep info for aggregate tables
all_names = {};
all_big_inc = {};
all_anat = {};
all_inc = [];
all_added_anat = {};
all_rate_change = [];
pt_names = {};
all_t = [];

for p = whichPts
    pt_name = pt(p).name;

    %% Load spike file
    if exist([spike_folder,sprintf('%s_spikes.mat',pt_name)],'file') == 0
        continue;
    end
      
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
    
    %% Flip things that we think are bad to bad
    if clean_blocks
        spikes = clean_missed_bad_blocks(spikes);
    end
    
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
            added_anatomy = cell(length(added_labels),1);
        else
            unchanged_anatomy = pt(p).ieeg.file(change(c).files(2)).anatomy(unchanged_idx);
            added_anatomy = pt(p).ieeg.file(change(c).files(2)).anatomy(added_idx);
        end
        all_added_anat{end+1} = added_anatomy;
        
        all_rate = [];
        all_rl = [];
        coa_post = [];
        rate_post = [];
        all_nseq = [];
        all_seq = {};
        last_block = zeros(last_file-1,1);
        findices = [];
        bindices = [];
        all_cos = [];
        all_dist = {};
        all_alt_cos = [];
        all_added_rate= [];
        all_un_cos = [];
        
        % Loop over files
        for f = 1:last_file
            
            flocs = pt(p).ieeg.file(f).locs;
            fdist = distance_from_closest_added(flocs,added_locs);
            
            chLabels = clean_labels_2(spikes.file(f).block(1).chLabels);
            nchs = length(chLabels);
            nblocks = length(spikes.file(f).block);
            rate = nan(nchs,nblocks); % default should be nans
            rl = nan(nchs,nblocks);
            coa = nan(nchs,nchs,nblocks);
            cos = nan(nchs,nblocks);
            un_cos = nan(nchs,nchs,nblocks);
            added_rate = nan(length(added_labels),nblocks);
            alt_cos = nan(length(added_labels),length(unchanged_labels),nblocks);
            num_seq = nan(nchs,nblocks);
            sdist = cell(nblocks,1);
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
                    
                    %% Get spike distances
                    
                    % find the channel indices corresponding to unchanged
                    % labels
                    unchanged_idx = find(ismember(chLabels,unchanged_labels));
                    
                    unchanged_spikes = ismember(gdf(:,1),unchanged_idx);
                    sdist{h} = fdist(gdf(unchanged_spikes,1));
                    
                     %% Get sequences, rl, coa
                     if ~isempty(gdf)
                        if f >= change(c).files(2)
                            is_post = 1;
                        else
                            is_post = 0;
                        end
                         
                        
                        %% Get cospike index
                        [cos(:,h),un_cos(:,:,h)] = co_spiking(gdf,fs,added_labels,chLabels);
                        [alt_cos(:,:,h),added_rate(:,h)] = alt_co_spiking(gdf,added_labels,unchanged_labels,chLabels);
                        if 0
                        [seq,rl(:,h),coa(:,:,h),num_seq(:,h),cos(:,h)] = new_get_sequences(gdf,nchs,fs,...
                            is_post,chLabels,added_labels,unchanged_labels);
                        else
                            seq = [];
                        end

                        if f >= change(c).files(2)
                            all_seq = [all_seq;seq];
                            
                        end
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
            rl(~lia,:) = [];
            cos(~lia,:) = [];
            un_cos(~lia,:,:) = [];
            un_cos(:,~lia,:) = [];

            %% Re-order as needed
            %[lia,locb] = ismember(new_labels,unchanged);
            [lia,locb] = ismember(unchanged,new_labels);
            new_labels = new_labels(locb);
            rate = rate(locb,:);
            rl = rl(locb,:);
            cos = cos(locb,:);
            un_cos = un_cos(lia,:,:);
            un_cos = un_cos(:,lia,:);
            if ~isequal(new_labels,unchanged)
                error('oh no');
            end
            %}
            
            all_rate = [all_rate,rate];
            all_rl = [all_rl,rl];
            all_dist = [all_dist;sdist];
            all_un_cos = cat(3,all_un_cos,un_cos);
            
            
            %% If it's a post-change file, add the coa matrix
            if f >= change(c).files(2)
                coa_post = cat(3,coa_post,coa);
                all_cos = [all_cos,cos];
                rate_post = [rate_post,rate];
                all_nseq = [all_nseq,num_seq];
                post_labels = chLabels;
                all_alt_cos = cat(3,all_alt_cos,alt_cos);
                all_added_rate = [all_added_rate,added_rate];
                
                if size(all_alt_cos,3) ~= size(all_added_rate,2)
                    error('what');
                end
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
        
    if 1
        %% Spike count
        fprintf('\n%s had %d spikes\n',pt_name,nansum(all_rate(:)));
    end
        
    end

    if 0
       change_dist_time(all_dist,block_dur,change_block,name,results_folder,surround,all_rate); 
        
    end
    
    %%
    if 0
    cosa = cosi_analysis(all_alt_cos,all_rate,change_block,surround,unchanged_labels,name,...
        all_added_rate);
    end
    
    
    %% Overall change in spike rate
    if 0
        show_overall_rate(all_rate,block_dur,change_block,run_dur,name,results_folder,surround)
    end
    
    %% clustered
    if 0
        %clust = {{'RB1','RB7','RB2'},{'LA4','LA2'}};
        clust = {{'LB8','LB9'};{'RB2','RB3'}};
        show_clusters_rate(all_rate,block_dur,change_block,run_dur,name,results_folder,clust,[],unchanged)
        %{
        if strcmp(name,'HUP128')
            clusts{1} = find(contains(unchanged,'L')); clusts{2} = find(contains(unchanged,'R'));
            clust_names = {'Left','Right'};
            show_clusters_rate(all_rate,block_dur,change_block,run_dur,name,results_folder,clusts,clust_names,unchanged)
        end
        %}
    end
    
    %% Overall rate increase
    mean_rate = nanmean(all_rate,1);
    mean_pre = nanmean(mean_rate(change_block-surround:change_block-1));
    mean_post = nanmean(mean_rate(change_block+1:change_block+surround));
    mean_change = (mean_post-mean_pre)./mean_pre;
    all_rate_change = [all_rate_change;mean_change];
    pt_names = [pt_names;pt_name];
    
    %% Get rate increase of electrodes with minimum spike rate
    [pre,post,cosi] = new_rate_increase(all_rate,change_block,all_cos,surround);
    abs_increase = (post-pre);
    rel_increase = (post-pre)./pre;
    
    %% Histogram of rate increase
    if 0
    histogram_rate_change(abs_increase,unchanged_labels);
    end
    
    %% AES plot
    if 0
        aes_plot(all_rate,block_dur,change_block,run_dur,unchanged_locs,added_locs,...
        name,results_folder,unchanged_labels)
    end
    
    %% Identity of spikiest electrode
    if 0
        spikiest_elec(all_rate,unchanged_labels,change_block,surround,...
            results_folder,name,block_dur,run_dur,unchanged_locs)
        
    end
   
    
    %% Restrict to minimum number of spikes
    min_num = 1;
    [spikey_idx,~] = find_spikey_elecs(all_rate,min_num,change_block,surround);

    spikey_labels = unchanged(spikey_idx);
    spikey_rate_inc = abs_increase(spikey_idx);
    coa_spikey_idx = ismember(post_labels,spikey_labels);
    coa_added_idx = ismember(post_labels,added_labels);
    
    %% SAR model
    if 0
        do_rel = 0;
        thing = cosi;
        ttext = 'cospike';
        fprintf('\n\n\n%s\n',name);
        outfolder = [results_folder,ttext,'/'];
        if ~exist(outfolder,'dir')
            mkdir(outfolder)
        end
        
        rprep(thing,all_rate,change_block,surround,unchanged_locs,spikey_idx,...
            outfolder,name,do_rel,ttext,unchanged_labels,all_un_cos)
        %{
        perm_sar(thing,all_rate,change_block,surround,unchanged_locs,spikey_idx,outfolder,...
    name,do_rel,ttext)
        %}
    end
    
    %{
    if 0
        
        thing = cosi;
        ttext = 'cospike';
        outfolder = [results_folder,ttext,'/'];
        if ~exist(outfolder,'dir')
            mkdir(outfolder)
        end
        do_abs = 0;
        do_simple = 1;
        tstat = corr_time_perm(thing,all_rate,change_block,surround,spikey_idx,do_abs,...
    ttext,outfolder,unchanged_labels,name,do_simple);
        all_t = [all_t;tstat];
        
    end
    %}
    
     %% Spatial clustering of rate increase
    if 0
    do_moran(abs_increase,unchanged_locs,spikey_idx,added_locs,name,results_folder,...
        run_dur,all_rate,change_block,surround)
    end
    
    %{
    [rate_increase,spikey_labels,spikey_idx,mean_rate_post,...
        big_inc_labels,big_inc_rate_sorted,abs_increase] = find_rate_increase(all_rate,change_block,unchanged);
    % Get anatomy of the big increase electrodes
    [~,big_inc_chs] = ismember(big_inc_labels,unchanged_labels);
    big_inc_anatomy = unchanged_anatomy(big_inc_chs);
    
    all_big_inc = [all_big_inc;big_inc_labels];
    all_anat = [all_anat;big_inc_anatomy];
    all_names = [all_names;repmat(cellstr(pt_name),length(big_inc_anatomy),1)];
    all_inc = [all_inc; big_inc_rate_sorted];
    
    if sum(spikey_idx) == 0
        continue;
    end
    %}
    
    if 0
        rate_order_stability(all_rate,all_rl,change_block,surround,block_dur,pt_name,results_folder)
    end
    
    %% Look at electrodes by spike rate change
    % A way to validate spikes
    if 0
    scatterplot_rate_change(spikey_rate_inc,spikey_labels,spikes,p,name)
    end
    
    %% Get anatomy of spikey electrodes
    spikey_anatomy = unchanged_anatomy(spikey_idx);
    
    % group anatomy
    [ana_lat,ana_loc] = anatomy_grouper(spikey_anatomy);
    
    %% Is there a difference in spike rate change based on anatomy?
    if 0
    rchange = rel_increase(spikey_idx);
    if sum(cell2mat(cellfun(@(x) ~isempty(x),spikey_anatomy,'uniformoutput',false))) ~= 0
        group_rate_change_by_anatomy(ana_lat,ana_loc,rchange,name,results_folder)
    end
    end
    

    

    %% Are electrodes with bigger spike rate increase closer to the new electrodes than are other spikey electrodes?
    %{
    dist_spikey = dist(spikey_idx);
    if 0
        
        % Get distance from closest new electrode of these spikey electrodes
        
        show_anatomy = 0;
        if show_anatomy
            plot_inc_as_fcn_of_dist([],spikey_rate_inc,dist_spikey,...
            spikey_anatomy,name,results_folder,run_dur,[])
        else
            plot_inc_as_fcn_of_dist([],spikey_rate_inc,dist_spikey,...
            spikey_labels,name,results_folder,run_dur,[])
        end

    end

   if 0
      new_distance(rel_increase,dist,spikey_idx,unchanged_labels,name,results_folder); 
       
   end
   
   if 0
   cosi_over_time(all_rate,cosi,spikey_idx,unchanged_labels,name,results_folder,...
    block_dur,change_block,run_dur,surround)
   end
    
    %% Are eletrodes with bigger spike rate increase more likely to co-spike with new electrodes?
    % NOTE WHETHER IM PASSING ABS INCREASE OR REL INCREASE
    if 0
        do_rel = 1;
        new_cos(post,cosi,abs_increase,spikey_idx,unchanged_labels,name,results_folder,dist,rel_increase,do_rel)
        %{
        blocks = 1:100;
        [cos,unchanged_spikey_labels] = ...
            co_spike_index(post_rate,spikey_idx,coa_post,...
            blocks,added,unchanged,post_labels,new_post_labels,...
            spikey_rate_inc,pt_name,run_dur,results_folder);
        %}

        
    
    end
    %}
    
    if 0
        raster_rate_chs(all_rate,block_dur,change_block,...
    unchanged,name,results_folder,findices,bindices,p,spikes,...
    spikey_idx,surround)
        
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

if 0
    multi_pt_added_info(all_rate_change,all_added_anat,results_folder,surround,pt_names)
end

end

