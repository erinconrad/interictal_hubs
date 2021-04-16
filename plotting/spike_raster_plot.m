function spike_raster_plot(whichPts)

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


for p = whichPts
    pt_name = pt(p).name;


    %% Load spike file
    spikes = load([spike_folder,sprintf('%s_spikes.mat',pt_name)]);
    spikes = spikes.spikes;
    name = spikes.name;
    nfiles = length(spikes.file);
    
   
    for f = 1:nfiles
        chLabels = spikes.file(f).hour(1).chLabels;
        nchs = length(chLabels);
        nhours = length(spikes.file(f).hour);
        hours = 1:nhours;
        raster = zeros(nchs,nhours);
        rl = nan(nchs,nhours);
        nseq = nan(nhours,1);
        all_chs = 1:nchs;
        skip_chs = spikes.file(f).hour(1).skip.all;
        rate_rl_corr = nan(nhours,1);
        for h = 1:nhours
            
            if spikes.file(f).hour(h).run_skip == 1
                raster(:,h) = nan;
            end
            
            
            gdf = spikes.file(f).hour(h).gdf;
            
            
            [rl(:,h),nseq(h)] = get_sequences(gdf,nchs);
            
            for ich = 1:nchs
                
                if ismember(ich,spikes.file(f).hour(h).bad) ||...
                        ismember(ich,spikes.file(f).hour(h).skip.all)
                    raster(ich,h) = nan;
                    continue
                end
                
                
                if isempty(gdf), continue; end
                raster(ich,h) = sum(gdf(:,1)==ich);
            end
            
            %% Correlate spike frequency with rl
            rate_rl_corr(h) = corr(raster(:,h),rl(:,h),...
                'Type','Spearman','rows','pairwise');
            
        end
        
        %% Get sequence reliability
        sr = rl_stability(rl,nseq);
        
        %% Take the median RL across all hours
        median_rl = nanmedian(rl,2);
        
        %% Re-order channels basedon this (for display purposes)
        [~,I] = sort(median_rl);
        
        
        
        
        %% Plot recruitment latency
        if 1
            figure
            set(gcf,'position',[440 1 835 797])
            ha = tight_subplot(2,1,[0 0.01],[0.06 0.04],[0.07 0.01]);
            axes(ha(1))
            turn_nans_white(rl(I,:));
            yticklabels([])
            ylabel('Electrode')
            xlim([1 nhours])
            title(sprintf('%s file %d',name,f))
            set(gca,'fontsize',15)
            
            axes(ha(2))
            plot(hours(nseq>=median(nseq)),sr(nseq>=median(nseq)),'o','linewidth',2)
            hold on
            plot(hours(nseq<median(nseq)),sr(nseq<median(nseq)),'x','linewidth',2)
            ylim([-1 1])
            xlim([1 nhours])
            xlabel('Hour')
            set(gca,'fontsize',15)
 
        end
        
        %% Correlate recruitment latency with spike rate
        if 1
            figure
            plot(rate_rl_corr,'o')
            ylim([-1 1])
            
        end
        
        %% Show spike rate
        if 1
            figure
            set(gcf,'position',[440 1 835 797])
            turn_nans_white(raster(I,:));
            yticks(1:nchs)
            yticklabels(chLabels(I))
            ylabel('Electrode')
            xlabel('Hour')

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