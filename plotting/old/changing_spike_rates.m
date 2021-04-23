function changing_spike_rates(whichPts)

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
spike_folder = [results_folder,'old_april19/spikes_april19/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load file with info about which patients are good
T = readtable([data_folder,'detector_parameters.xlsx']);
goodness = T.num_real_outOf50__1; % CHANGE AS NEEDED
names = T.Patient;

if isempty(whichPts)
    listing = dir([spike_folder,'*.mat']);
    for i = 1:length(listing)
        C = listing(i).name;
        temp_name = strsplit(C,'_');
        temp_name = temp_name{1};
        for j = 1:length(pt)
            pt_name = pt(j).name;
            if strcmp(temp_name,pt_name)
                whichPts = [whichPts,j];
                break
            end
        end
    end
end

%% Only do good ones
skip = zeros(length(whichPts),1);
for i = 1:length(whichPts)
    name = pt(i).name;
    tb_idx = find(strcmp(names,name));
    if goodness(tb_idx) <= 30
        skip(i) = 1;
    end
end
whichPts(logical(skip)) = [];


for p = whichPts
    pt_name = pt(p).name;


    %% Load spike file
    if exist([spike_folder,sprintf('%s_spikes.mat',pt_name)],'file') == 0
        continue;
    end
      
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
    name = spikes.name;
    nfiles = length(spikes.file);
    
   
    for f = 1:nfiles
        chLabels = spikes.file(f).block(1).chLabels;
        nchs = length(chLabels);
        nblocks = length(spikes.file(f).block);
        blocks = 1:nblocks;
        raster = zeros(nchs,nblocks);
        
        for h = 1:nblocks
            
            if spikes.file(f).block(h).run_skip == 1
                raster(:,h) = nan;
            end
            
            
            gdf = spikes.file(f).block(h).gdf;
            if isempty(gdf), continue; end
            
            for ich = 1:nchs
                
                if ismember(ich,spikes.file(f).block(h).bad) ||...
                        ismember(ich,spikes.file(f).block(h).skip.all)
                    raster(ich,h) = nan;
                    continue
                end
                
                
                
                raster(ich,h) = sum(gdf(:,1)==ich);
            end
            

        end
        
        %% Plot as a bunch of lines
        % restrict to top 10 channels
        total_spikes = nansum(raster,2);
        [nspikes,I] = sort(total_spikes,'descend');
        for i = 1:10
        	plot(raster(I(i),:))
            hold on
            
        end
    end
    
end



end