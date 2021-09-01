function all_ad = get_ad_details(p)

only_depth = 0;

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
ad_folder = [results_folder,'ad/'];
bct_folder = locations.bct;
addpath(genpath(bct_folder));

%% Load pt file
pt = load([data_folder,'pt.mat']);
pt = pt.pt;
name = pt(p).name;

%% Load ad file
if exist([ad_folder,sprintf('%s_ad.mat',name)],'file') == 0
    out = [];
    return;
end

ad = load([ad_folder,sprintf('%s_ad.mat',name)]);
ad = ad.ad;

nfiles = length(ad.file);

%% Identify files with a change in electrodes
[change,no_change_ever] = alt_find_electrode_change_files(ad,only_depth);

%[change,no_change_ever] = find_electrode_change_files(pt,p,only_depth);
nchanges = length(change);
c = nchanges;

last_file = nfiles;

added = change(c).added;
unchanged = no_change_ever;%change(c).unchanged;

%% Get total number of blocks
nb = 0;
for f = 1:last_file
    nblocks_in_file = length(ad.file(f).block);
    if strcmp(name,'HUP136') && f == 1
        last_good_block = fix_hup136;
        nblocks_in_file = last_good_block;
    end
    
    nb = nb + nblocks_in_file;
end

%% Initialize things
b_count = 0;
all_ad = nan(length(unchanged),nb);
bdur = pt(p).ieeg.file(1).block(1).end - pt(p).ieeg.file(1).block(1).start;
bdur = bdur/3600; %convert to hours
% Loop over files
for f = 1:last_file

    nblocks = length(ad.file(f).block);
    % fix for hup136
    if strcmp(name,'HUP136') && f == 1
        last_good_block = fix_hup136;
        nblocks = last_good_block;
    end
    
    % Loop over blocks
    for h = 1:nblocks
        block = ad.file(f).block(h);

        run_labels = block.run_labels;
        b_count = b_count + 1;

        if block.run_skip == 1
            continue; % leave the whole block as nans
        end
        
        ad_rat = block.ad;
        
        %% Remove and pad as needed
        % Remove electrodes that are not unchanged
        [lia,locb] = ismember(run_labels,unchanged);
        new_labels = run_labels;
        new_labels(~lia) = [];
        ad_rat(~lia) = [];

        % Pad with missing electrodes (primarily those that we excluded from run due to baddness)
        missing_idx = ~ismember(unchanged,run_labels);
        missing_labels = unchanged(missing_idx);
        ad_rat = [ad_rat; nan(length(missing_labels),1)];
        new_labels = [new_labels;missing_labels];

        % Re-order as needed
        [lia,locb] = ismember(unchanged,new_labels);
        new_labels = new_labels(locb);
        ad_rat = ad_rat(locb);

        if ~isequal(new_labels,unchanged)
            error('ruh roh');
        end
        
        all_ad(:,b_count) = ad_rat;

    end
    
    if f + 1 == change(c).files(2)
        change_block = b_count;
    end


end

times = 1:size(all_ad,2);
times = times*bdur;
change_block = change_block*bdur;

if 0
    
    figure
    set(gcf,'position',[100 100 1000 600])
    tiledlayout(2,2)
    
    %% Alpha delta ratio
    nexttile([1 2])
    plot(times,nanmean(all_ad,1),'linewidth',2)
    hold on
    plot([change_block change_block],ylim,'r--','linewidth',3)
    ylabel('Alpha-delta ratio')
    xlabel('Hour');
    
    %% PSD
    nexttile
    mean_ad = nanmean(all_ad,1);
    mean_ad(isnan(mean_ad)) = nanmedian(mean_ad);
    mean_ad = mean_ad - mean(mean_ad);
    X = mean_ad;
    Y = fft(X);
    P = abs(Y).^2;
    fs = 1/bdur;
    freqs = linspace(0,fs,length(P)+1);
    freqs = freqs(1:end-1);
    P = P(1:ceil(length(P)/2)); % Take first half
    freqs = freqs(1:ceil(length(freqs)/2));
    P = P(1./freqs<100);
    freqs = freqs(1./freqs<100);
    plot(1./freqs,P,'linewidth',2);
    xlabel('Period (hours)')
    ylabel('Power');
    
    %% Wavelet power
    nexttile
    [cfs,periods] = cwt(X,hours(bdur));
    plot(periods,sqrt(nanmean((abs(cfs)).^2,2)),'linewidth',2)
    xlabel('Period (hours)')
    ylabel('Magnitude')
    pause
    close(gcf)
    
end

end