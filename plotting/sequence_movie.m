function sequence_movie(whichPts)

%% Parameters
filt = 1;
timing = 'mid_rise';

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
spike_folder = [results_folder,'spikes/'];

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Load file with info about which patients are good
T = readtable([data_folder,'detector_parameters.xlsx']);
goodness = T.num_real_outOf50_; % CHANGE AS NEEDED
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
end

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
    
    
    % Loop over files
    for f = 1:nfiles
        
        chLabels = spikes.file(f).block(1).chLabels;
        nchs = length(chLabels);
        nblocks = length(spikes.file(f).block);
        
        rate = nan(nchs,nblocks);
        rl = nan(nchs,nblocks);
        coa = nan(nchs,nchs,nblocks);
        
        % Loop over blocks
        for h = 1:nblocks
            block = spikes.file(f).block(h);
            
            if block.run_skip == 1
                continue;
            end
            
            %% Get spike channels and times (whatever time you want)
            chs = block.details.filter(filt).gdf(:,1);
            times = block.details.filter(filt).(timing);
            fs = block.fs;
            
            %% Get sorted spike indices and chs
            [times,I] = sort(times);
            chs = chs(I);
            
            %% Construct gdf
            gdf = [chs,times];
            
            %% Spike rate
            for ich = 1:nchs
                if ismember(ich,block.bad) ||...
                        ismember(ich,block.skip.all)
                    rate(ich,h) = nan;
                    continue
                end
                rate(ich,h) = sum(gdf(:,1)==ich);
            end
            
            %% Get sequences, rl, coa
            [~,rl(:,h),coa(:,:,h)] = new_get_sequences(gdf,nchs,fs);
            
            
        end
        
        %% Raster plot of rate
        if 0
        r = corr(rate,nanmean(rate,2),'type','pearson','rows','pairwise');   
        figure
        turn_nans_white(rate)
        title(sprintf('Rate reliability: %1.1f',nanmean(r)))
        end
        
        %% Raster plot of rl
        if 1
        [~,I] = sort(nanmean(rl,2));
        r = corr(rl,nanmean(rl,2),'type','pearson','rows','pairwise');
        figure
        turn_nans_white(rl)
        title(sprintf('RL reliability: %1.1f',nanmean(r)))
        end
        
        %% gif of coa matrix
        if 0
        fig = figure;
        for h = 1:nblocks
            imagesc(coa(:,:,h))
            pause(0.2)
        end
        end
        
        %% Details of coa matrix
        ns = squeeze(sum(coa,1));
        
        
        %% Raster plot of ns
        if 0
        r = corr(ns,nanmean(ns,2),'type','pearson','rows','pairwise');   
        figure
        turn_nans_white(ns)
        title(sprintf('Node strength reliability: %1.1f',nanmean(r)))
        end
        
        
    end

end
end


function turn_nans_white(im)
    % white
    cmap = colormap;
    nanjet = [ 1,1,1; cmap  ];
    nanjetLen = length(nanjet); 
    pctDataSlotStart = 2/nanjetLen;
    pctDataSlotEnd   = 1;
    pctCmRange = pctDataSlotEnd - pctDataSlotStart;

    dmin = nanmin(im(:));
    dmax = nanmax(im(:));
    dRange = dmax - dmin;   % data range, excluding NaN

    cLimRange = dRange / pctCmRange;
    cmin = dmin - (pctDataSlotStart * cLimRange);
    cmax = dmax;
    imagesc(im);
    set(gcf,'colormap',nanjet);
    caxis([cmin cmax]);
end