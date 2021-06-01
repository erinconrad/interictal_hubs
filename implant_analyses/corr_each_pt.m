function corr_each_pt(whichPts,saved_out)

%{
I should incorporate some sort of SAR model as well

I should see if I can "whiten" the data for plotting

Try FC for correlation
%}

%% Parameters
surround = 48;
nb = 1e4;
which_pred = 'dist';
which_rate = 'abs';

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
all_zs = nan(length(whichPts),1);

for i = 1:length(whichPts) 
    
    nexttile
    
    rate = out(i).rate;
    cblock = out(i).change_block;
    
    % Identify pre and post times
    pre = cblock-surround:cblock - 1;
    post = cblock + 1: cblock+surround;
    
    % Calculate change in rate
    switch which_rate
        case 'abs'
            rate_change = nanmean(rate(:,post),2) - nanmean(rate(:,pre),2);
            rate_change = rate_change./out(i).run_dur;
            rtext = 'Spike rate change (spikes/min)';
            ekg = identify_ekg_scalp(out(i).unchanged_labels);
            spikey_idx = ~ekg;
        case 'rel'
            rate_change = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre),2);
            spikey_idx = nanmean(rate(:,[pre,post]),2) > 0.5;
            rtext = 'Relative spike rate change';
    end
    
    
    % get predictor
    switch which_pred
        case 'dist'
            predictor = out(i).dist;
            
            ptext = 'Distance from added electrodes (mm)';
        case 'cosi'
            cos = out(i).cos;
            cos_post = cos(:,1:end);
            cosi = nanmean(cos_post./rate(:,cblock+1:end),2);
            predictor = cosi;
            ptext = 'Co-spike index with added electrodes';
    end
    
    % Remove non-spikey if doing rel
    rate_change(~spikey_idx) = [];
    predictor(~spikey_idx) = [];
    
    % Monte carlo test to see how big correlation is compared to that at
    % different time points
    [rho,pval] = mc_corr(rate(spikey_idx,:),predictor,cblock,surround,nb,which_rate);
    
    % Fisher's R to z transformation
    z = atanh(rho);
    all_zs(i) = z;
    
    % Plot
    plot(predictor,rate_change,'o','linewidth',2);
    title(sprintf('rho = %1.2f, p = %1.3f',rho,pval))   
    set(gca,'fontsize',20)
    
end

ylabel(ax,rtext,'fontsize',20)
xlabel(ax,ptext,'fontsize',20)

%% Stouffer's Z-score
z = sum(all_zs)/sqrt(length(all_zs));

% get two sided p-value
pval = 2*normcdf(-abs(z));

% get r back
r = tanh(z);
fprintf('\nCombined r = %1.2f, p = %1.3f\n',r,pval);

end


function [true_rho,pval] = mc_corr(rate,predictor,cblock,surround,nb,which_rate)

nblocks = size(rate,2);

%% Get true corr
% Identify pre and post times
pre = cblock-surround:cblock - 1;
post = cblock + 1: cblock+surround;

% Calculate change in rate
switch which_rate
    case 'abs'
        rate_change = nanmean(rate(:,post),2) - nanmean(rate(:,pre),2);
    case 'rel'
        rate_change = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre),2);
end

true_rho = corr(rate_change,predictor,'Type','Pearson','rows','pairwise');

mc_rho = nan(nb,1);
for ib = 1:nb
    
    % Make a fake change time
    fchange = randi([surround+1,nblocks-surround]);
    
    % recalculate change around this time
    fpre = fchange-surround:fchange-1;
    fpost = fchange+1:fchange+surround;
    
    switch which_rate
        case 'abs'
            rate_change = nanmean(rate(:,fpost),2) - nanmean(rate(:,fpre),2);
        case 'rel'
            rate_change = (nanmean(rate(:,fpost),2) - nanmean(rate(:,fpre),2))./nanmean(rate(:,fpre),2);
    end
    
    mc_rho(ib) = corr(rate_change,predictor,'Type','Pearson','rows','pairwise');
    
end


%% Determine number with as or more extreme a correlation
num_as_sig = sum(abs(mc_rho) >= abs(true_rho)); % those more positive or more negative

pval = (num_as_sig+1)/(nb+1);


end