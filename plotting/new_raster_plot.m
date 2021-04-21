function new_raster_plot(whichPts)

%% Parameters
cluster_time = 60*30;

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
    
   
    for f = 1:nfiles
        chLabels = spikes.file(f).block(1).chLabels;
        nchs = length(chLabels);
        nblocks = length(spikes.file(f).block);
        blocks = 1:nblocks;
        rate_raster = zeros(nchs,nblocks);
        rl = nan(nchs,nblocks);
        all_rl = [];
        seconds = [];
        nseq = nan(nblocks,1);
        all_chs = 1:nchs;
        skip_chs = spikes.file(f).block(1).skip.all;
        rate_rl_corr = nan(nblocks,1);
        for h = 1:nblocks
            
            if spikes.file(f).block(h).run_skip == 1
                rate_raster(:,h) = nan;
            end
            
            
            %% Get spike details as a table
            T = convert_details_to_table(spikes,f,h);

            %% Get sequences
            times = T.peak_idx./T.fs + T.run_start;
            [times,I] = sort(times);
            T = T(I,:);
            chs = T.ch;
            gdf = [chs,times];
            if isempty(gdf), continue; end
            [seq,rl_ind] = new_get_sequences(gdf,nchs);
            all_rl = [all_rl,rl_ind];
            seconds = [seconds,...
                repmat(spikes.file.block(h).run_times(1),...
                1,size(rl_ind,2))];
        
            for ich = 1:nchs
                
                if ismember(ich,spikes.file(f).block(h).bad) ||...
                        ismember(ich,spikes.file(f).block(h).skip.all)
                    rate_raster(ich,h) = nan;
                    continue
                end
                
                
                
                rate_raster(ich,h) = sum(gdf(:,1)==ich);
                
                if isempty(seq), continue; end
                
                if size(rl_ind,2) == 1
                    rl(ich,h) = rl_ind(ich);
                else
                    rl(ich,h) = nanmean(rl_ind(ich,2));
                end
            end
  
        end
        
        %% Cluster rl by time
        %clusters = cluster_by_time(all_rl,seconds,cluster_time);
        
        %% Take the median RL across all blocks
        mean_rl = nanmean(rl,2);
        
        %% Re-order channels basedon this (for display purposes)
        [~,I] = sort(mean_rl);
        
        %% Get sequence reliability
        sr = nan(size(rl,2),1);
        for s = 1:size(rl,2)
            r = corr(mean_rl,rl(:,s),'Type','Spearman','Rows','pairwise');
            sr(s) = r;
        end
        
        %}
        
        
        %% Plot recruitment latency
        if 1
            figure
            set(gcf,'position',[1 8 651 797])
            ha = tight_subplot(2,1,[0 0.01],[0.06 0.04],[0.07 0.01]);
            axes(ha(1))
            turn_nans_white(rl(I,:));
            yticklabels([])
            ylabel('Electrode')
            xlim([1 nblocks])
            title(sprintf('%s file %d',name,f))
            set(gca,'fontsize',15)
            
            axes(ha(2))
            plot(blocks(nseq>=median(nseq)),sr(nseq>=median(nseq)),'o','linewidth',2)
            hold on
            plot(blocks(nseq<median(nseq)),sr(nseq<median(nseq)),'x','linewidth',2)
            ylim([-1 1])
            xlim([1 nblocks])
            xlabel('block')
            set(gca,'fontsize',15)
 
        end
        
        %% Correlate recruitment latency with spike rate
        if 0
            figure
            plot(rate_rl_corr,'o')
            ylim([-1 1])
            
        end
        
        %% Show spike rate
        if 0
            figure
            set(gcf,'position',[721 4 720 797])
            turn_nans_white(rate_raster);
            yticks(1:nchs)
            yticklabels(chLabels)
            ylabel('Electrode')
            xlabel('block')

            title(sprintf('%s file %d',name,f))


            
        end
        
        pause
        close all
        
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