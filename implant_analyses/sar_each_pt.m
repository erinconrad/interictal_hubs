function sar_each_pt(whichPts,saved_out)

%% Parameters
surround = 48;
nb = 1e4;
which_resp = 'abs_rate';
which_pred = 'dist';
which_sar = 'pc';

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
main_spike_results = [results_folder,'main_spikes/'];
if ~exist(main_spike_results,'dir')
    mkdir(main_spike_results);
end

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
spike_folder = [results_folder,'new_spikes/'];

if isempty(whichPts)
    whichPts = [20 103 106 107 35 109 110 111 94 97];
end

if saved_out == 1
    
    out = load([main_spike_results,'out.mat']);
    out = out.out;
    
else
    out = initialize_out_struct(length(whichPts));
    
    %% Get spike details
    fprintf('Getting spike details for pt...\n');
    for i = 1:length(whichPts)
        p = whichPts(i);
        fprintf('%d of %d\n',i,length(whichPts));
        out(i) = get_gdf_details(p);
    end
    save([main_spike_results,'out'],'out');
end

%% Initialize figure
figure
set(gcf,'position',[1 100 1400 550])
ax = tiledlayout(2,5,'TileSpacing','compact','padding','compact');
all_zs = nan(length(whichPts),3);

for i = 1:length(whichPts) 
    
    nexttile
    
    rate = out(i).rate;
    if isempty(out(i).metrics)
        ns = nan(size(rate));
    else
        ns = out(i).metrics.ns_norm;
    end
    cblock = out(i).change_block;
    
    % Identify pre and post times
    if isempty(surround)
        pre = 1:cblock-1;
        post = cblock+1:size(rate,2);
    else
        [pre,post] = get_surround_times(rate,cblock,surround);
    end
    
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    
    % Calculate change in rate
    switch which_resp
        case 'abs_rate'
            resp = nanmean(rate(:,post),2) - nanmean(rate(:,pre),2);
            resp = resp./out(i).run_dur;
            rtext = 'Spike rate change (spikes/min)';     
            spikey_idx = ~ekg;
        case 'rel_rate'
            resp = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre),2);
            spikey_idx = nanmean(rate(:,[pre,post]),2) > 0.5 & ~ekg;
            rtext = 'Relative spike rate change';
        case 'ns'
            resp = nanmean(ns(:,post),2) - nanmean(ns(:,pre),2);
            rtext = 'Normalized node strength change';
            spikey_idx = ~ekg;
            
    end
    
    
    % get predictor
    switch which_pred
        case 'dist'
            predictor = out(i).dist;
            
            ptext = 'Distance from added electrodes (mm)';
        case 'cosi'
            %cos = out(i).cos;
            %cos_post = cos(:,1:end);
            %cosi = nanmean(cos_post./rate(:,cblock+1:end),2);
            predictor = out(i).cosi;
            ptext = 'Co-spike index with added electrodes';
        case 'ns'
            if isempty(out(i).metrics)
                predictor = nan(size(out(i).dist));
            else
                predictor = out(i).metrics.added_pc;
            end
            ptext = 'Functional connectivity with added electrodes';
    end
    
    
    % Get locs and pc
    locs = out(i).unchanged_locs;
    if isempty(out(i).metrics)
        pc = nan(size(locs,1),size(locs,1));
    else
        pc = out(i).metrics.avg_mat;
    end
    
    % Remove non-spikey
    resp(~spikey_idx) = [];
    predictor(~spikey_idx) = [];
    locs(~spikey_idx,:) = [];
    pc(~spikey_idx,:) = [];
    pc(:,~spikey_idx) = [];
    
    
    % Do SAR model
    [coeffs,pvals] = sar_model(predictor,resp,locs,pc,which_sar);
    [trad_rho,trad_p] = corr(predictor,resp,'type','pearson','rows','pairwise');
    
    % Plot
    plot(predictor,resp,'o')
    title(sprintf('Simple: r = %1.2f, p = %1.3f\nSAR: b = %1.2f, p = %1.3f',...
        trad_rho,trad_p,coeffs(1),pvals(1)))
    
end


end