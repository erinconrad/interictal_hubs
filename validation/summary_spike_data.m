function summary_spike_data(which_ver)

%% Get file locs
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

if which_ver == 1
    sp_folder = [results_folder,'spikes/'];
elseif which_ver == 2
    sp = [results_folder,'new_spikes/'];
end

% Loop over spike files
listing = dir([sp_folder,'*mat']);
for i = 1:length(listing)
    spikes = load([sp_folder,listing(i).name]);
    spikes = spikes.spikes;
    
    nbad = [];
    nskip = [];
    nspikes = 0;
    nblocks = 0;
    name = spikes.name;
    tmul1 = spikes.file(1).block(1).params.tmul;
    absthresh1 = spikes.file(1).block(1).params.absthresh;
    
    for f = 1:length(spikes.file)
        for h = 1:length(spikes.file(f).block)
            nbad = [nbad;length(spikes.file(f).block(h).bad)];
            if isempty(spikes.file(f).block(h).skip)
                nskip = [nskip;0];
            else
                nskip = [nskip;length(spikes.file(f).block(h).skip.all)];
            end
            nspikes = nspikes + size(spikes.file(f).block(h).gdf,1);
            nblocks = nblocks + 1;
        end
    end
    
    nbad = mean(nbad);
    nskip = mean(nskip);
    
    fprintf(['\nFor %s, using tmul %d and absthresh %d\n'...
        'in %d blocks, %d spikes detected, (mean %1.1f contacts skipped and %1.1f rejected as artifact)\n'],...
        name,tmul1,absthresh1,nblocks,nspikes,nskip,nbad);
    
end

end