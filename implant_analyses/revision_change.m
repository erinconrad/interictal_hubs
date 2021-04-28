function revision_change(whichPts)

%% Parameters
filt = 2;
timing = 'peak';

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
bct_folder = locations.bct;
addpath(genpath(bct_folder));

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

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
    
    %% Identify files with a change in electrodes
    change = find_electrode_change_files(pt,p);
    nchanges = length(change);
    
    for c = 1:nchanges
        if c < nchanges
            last_file = change(c+1).files(2)-1;
        else
            last_file = nfiles;
        end
        added = change(c).added;
        unchanged = change(c).unchanged;
        all_rate = [];
        last_block = zeros(last_file-1,1);
        
        % Loop over files
        for f = 1:last_file
            
            chLabels = clean_labels_2(spikes.file(f).block(1).chLabels);
            nchs = length(chLabels);
            nblocks = length(spikes.file(f).block);
            rate = nan(nchs,nblocks);
            
            % Loop over blocks
            for h = 1:nblocks
                block = spikes.file(f).block(h);

                if block.run_skip == 1
                    continue;
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
                        if ismember(ich,block.bad) ||...
                                ismember(ich,block.skip.all)
                            rate(ich,h) = nan;
                            continue
                        end
                        rate(ich,h) = sum(gdf(:,1)==ich);
                    end

                else
                    rate(:,h) = 0;
                end
                
            end
            
            %% Get indices of unchanged and remove changed
            [lia,locb] = ismember(chLabels,unchanged);
            new_labels = chLabels;
            new_labels(~lia) = [];
            rate(~lia,:) =[];
            
            %% Re-order as needed
            [lia,locb] = ismember(new_labels,unchanged);
            new_labels = new_labels(locb);
            rate = rate(locb,:);
            if ~isequal(new_labels,unchanged)
                error('oh no');
            end
            
            all_rate = [all_rate,rate];
            if f<last_file
                last_block(f) = size(all_rate,2);
            end
            
            if f + 1 == change(c).files(2)
                change_block = size(all_rate,2);
            end
            
        end
        
    
    
    if 1
        figure
        set(gcf,'position',[1 1 1400 800])
        turn_nans_white(all_rate)
        hold on
        for b = 1:length(last_block)
            plot([last_block(b) last_block(b)],ylim,'k','linewidth',3)
        end
        plot([change_block change_block],ylim,'r','linewidth',3)
        yticks(1:length(unchanged))
        yticklabels(unchanged)
        title(added)
        
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