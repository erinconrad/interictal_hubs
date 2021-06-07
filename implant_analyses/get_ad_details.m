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
[change,no_change_ever] = find_electrode_change_files(pt,p,only_depth);
nchanges = length(change);
c = nchanges;

last_file = nfiles;

added = change(c).added;
unchanged = no_change_ever;%change(c).unchanged;

%% Get total number of blocks
nb = 0;
for f = 1:last_file
    nb = nb + length(ad.file(f).block);
end

%% Initialize things
b_count = 0;
all_ad = nan(length(unchanged),nb);
% Loop over files
for f = 1:last_file

    nblocks = length(ad.file(f).block);
    
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

if 1
    figure
    tiledlayout(2,1)
    nexttile
    turn_nans_white(all_ad)
    hold on
    plot([change_block change_block],ylim,'r--','linewidth',3)
    
    nexttile
    plot(nanmean(all_ad,1))
    hold on
    plot([change_block change_block],ylim,'r--','linewidth',3)
    pause
    close(gcf)
    
end

end