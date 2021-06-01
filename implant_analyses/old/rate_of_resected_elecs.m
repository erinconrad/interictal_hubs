function rate_of_resected_elecs(whichPts)

%% Parameters
filt = 2;
timing = 'peak';
do_smooth = 1;
sm_span = 1;
min_spike_rate = 1;

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
        %locs = pt(p).ieeg.file(f).locs;
        resected = pt(p).electrode_info.resected;
        res = resected_ch_nums(resected,chLabels);
        is_res = ismember(1:nchs,res);
        

        
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
        
        %% 
        avg_rate_by_ch = nanmean(rate,2);
        include_ch = avg_rate_by_ch > min_spike_rate;

        
        %% Plot rate over time according to resected vs non-resected status
        if 1
            figure
            set(gcf,'position',[100 280 1300 500])
            if do_smooth
                plot(smooth(nanmean(rate(is_res'&include_ch,:),1),sm_span))
                hold on
                plot(smooth(nanmean(rate(~is_res'&include_ch,:),1),sm_span))
            else
                
                plot(nanmean(rate(is_res,:),1))
                hold on
                plot(nanmean(rate(~is_res,:),1))
            end
            legend({'Resected','Not-resected'})
            xlabel('Block')
            ylabel('Average spike rate');
            title(sprintf('%s file %d',pt_name,f))
            set(gca,'fontsize',20)
            pause
            close(gcf)
        end
        
    end

end