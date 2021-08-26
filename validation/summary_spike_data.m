function summary_spike_data(which_ver)

%% Parameters
rm_dup = 1;
rm_sz = 1;

%% Get file locs
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
data_folder = [scripts_folder,'data/'];
addpath(genpath(scripts_folder));

if which_ver == 1
    sp_folder = [results_folder,'spikes/'];
elseif which_ver == 2
    sp_folder = [results_folder,'new_spikes/'];
elseif which_ver == 3
    sp_folder = [results_folder,'nina_spikes/'];
elseif which_ver == 4
    sp_folder = [results_folder,'revision_spikes/'];
elseif which_ver == 5
    sp_folder = [results_folder,'alt/spikes/'];
end

% Loop over spike files
listing = dir([sp_folder,'*mat']);
for i = 1:length(listing)
    spikes = load([sp_folder,listing(i).name]);
    spikes = spikes.spikes;
    
    pt = load([data_folder,'pt.mat']);
    pt = pt.pt;
    
    % Get corresponding pt
    name = spikes.name;
    foundit = 0;
    for p = 1:length(pt)
        if strcmp(pt(p).name,name)
            foundit = 1;
            break
        end
    end
    if foundit == 0
        error('what');e
    end
    
    nbad = [];
    nskip = [];
    nspikes = 0;
    nblocks = 0;
    n_dup = 0;
    n_sz = 0;
    name = spikes.name;
    tmul1 = spikes.file(1).block(1).params.tmul;
    absthresh1 = spikes.file(1).block(1).params.absthresh;
    
    for f = 1:length(spikes.file)
        
        if rm_sz
            sz_times = all_sz_times_in_file(pt,p,f);
        end
            
        
        for h = 1:length(spikes.file(f).block)
            nbad = [nbad;length(spikes.file(f).block(h).bad)];
            if isempty(spikes.file(f).block(h).skip)
                nskip = [nskip;0];
            else
                nskip = [nskip;length(spikes.file(f).block(h).skip.all)];
            end
            
            gdf = spikes.file(f).block(h).gdf;
            if isempty(gdf)
                continue
            end
            
            if rm_dup
                [gdf,n_dup_temp] = remove_duplicates(gdf);
                n_dup = n_dup + n_dup_temp;
            end
            
            if rm_sz
                [gdf,n_sz_temp]= remove_spikes_in_sz(gdf,sz_times);
                n_sz = n_sz + n_sz_temp;
            end
            
            nspikes = nspikes + size(gdf,1);
            nblocks = nblocks + 1;
        end
    end
    
    nbad = mean(nbad);
    nskip = mean(nskip);
    
    fprintf(['\nFor %s, using tmul %d and absthresh %d\n'...
        'in %d blocks, %d spikes detected \nmean %1.1f contacts skipped'...
        ', %1.1f rejected as artifact, %d duplicates removed, %d in sz removed\n'],...
        name,tmul1,absthresh1,nblocks,nspikes,nskip,nbad,n_dup,n_sz);
    
end

end