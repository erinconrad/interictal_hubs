function summary_spike_data

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.main_folder,'data/'];
results_folder = [locations.main_folder,'results/'];
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));
sp_folder = [results_folder,'spikes/'];

% Loop over spike files
listing = dir([sp_folder,'*mat']);
for i = 1:length(listing)
    spikes = load([sp_folder,listing(i).name]);
    spikes = spikes.spikes;
    
    nbad = [];
    nskip = [];
    nspikes = 0;
    nhours = 0;
    name = spikes.name;
    tmul1 = spikes.file(1).hour(1).params.tmul;
    absthresh1 = spikes.file(1).hour(1).params.absthresh;
    
    for f = 1:length(spikes.file)
        for h = 1:length(spikes.file(f).hour)
            nbad = [nbad;length(spikes.file(f).hour(h).bad)];
            nskip = [nskip;length(spikes.file(f).hour(h).skip)];
            nspikes = nspikes + size(spikes.file(f).hour(h).gdf,1);
            nhours = nhours + 1;
        end
    end
    
    nbad = mean(nbad);
    nskip = mean(nskip);
    
    fprintf(['\nFor %s, using tmul %d and absthresh %d\n'...
        'in %d hours, %d spikes detected, (mean %1.1f contacts skipped and %1.1f rejected as artifact)\n'],...
        name,tmul1,absthresh1,nhours,nspikes,nskip,nbad);
    
end

end