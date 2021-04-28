function sequence_movie(whichPts)

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
        locs = pt(p).ieeg.file(f).locs;
        resected = pt(p).electrode_info.resected;
        res = resected_ch_nums(resected,chLabels);
        clean_labs = clean_labels_2(chLabels);
        
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
                
                 %% Get sequences, rl, coa
                [~,rl(:,h),coa(:,:,h)] = new_get_sequences(gdf,nchs,fs);
            else
                rate(:,h) = 0;
            end
            
           
            
            
        end
        
        %% Raster plot of rate
        if 0
            if sum(sum(isnan(rate))) == size(rate,1)*size(rate,2)
                continue;
            end
        r = corr(rate,nanmean(rate,2),'type','pearson','rows','pairwise');   
        figure
        turn_nans_white(rate)
        yticks(1:nchs)
        yticklabels(chLabels)
        title(sprintf('%s file %d rate reliability: %1.1f',pt_name,f,nanmean(r)))
        pause
        close(gcf)
        end
        
        %% NMF rate
        if 1
        k = 2;
        A = rate;
        A(isnan(A)) = 0;
        [W,H] = rate_nmf(A,k);
        figure
        for j = 1:k
            subplot(1,k,j)
            scatter3(locs(:,1),locs(:,2),locs(:,3),100,W(:,j),'filled');
            hold on
            scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
            scatter3(locs(res,1),locs(res,2),locs(res,3),100,'rp');
            
            % Find main chs
            [~,I] = sort(W(:,j),'descend');
            main_chs = I(1:3);
            main_labs = chLabels(main_chs);
            title(main_labs)
        end
        
        figure
        plot(H');
        end
        
        
        %% Raster plot of rl
        if 0
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
            title(sprintf('Block %d of %d',h,size(coa,3)))
            yticks(1:nchs)
            xticks(1:nchs)
            yticklabels(chLabels)
            xticklabels(chLabels)
            pause
        end
        end
        
        %% GE/NS of coa matrix
        ns = squeeze(sum(coa,1));
        ge = nan(nblocks,1);
        for h = 1:nblocks
            ge(h) = efficiency_wei(coa(:,:,h));
        end
        if 0
        figure
        subplot(2,1,1)
        plot(mean(ns,1));
        subplot(2,1,2)
        plot(ge)
        end
        
        %% Raster plot of ns
        if 0
        r = corr(ns,nanmean(ns,2),'type','pearson','rows','pairwise');   
        figure
        turn_nans_white(ns)
        title(sprintf('Node strength reliability: %1.1f',nanmean(r)))
        end
        
        %% Gif of spike rate
        if 0 
        fig = figure;
        for h = 1:nblocks
            scatter3(locs(:,1),locs(:,2),locs(:,3),100,ns(:,h),'filled');
            hold on
            scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
            scatter3(locs(res,1),locs(res,2),locs(res,3),100,'rp');
            pause(0.2)
            hold off
        end
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