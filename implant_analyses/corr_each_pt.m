function corr_each_pt(whichPts,saved_out)


%% Parameters
surround = 24*0.5;
nb = 1e4;
which_resp = 'rel_rate';
which_pred = 'dist';
corr_type = 'Spearman';
min_rate = -1;
show_labels = 0;

%% Decide whether to do this!!
only_pre = 0; % for the MC analysis, compare to only the pre-revision times (in case the revision effect is delayed)

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
set(gcf,'position',[50 547 1000 800])
tiledlayout(4,3,'TileSpacing','compact','padding','compact')
tile_order = [1 2 4 5 7 8 10 11 3 6 9];

all_zs = nan(length(whichPts),7);
axesi = zeros(length(whichPts),1);

for i = 1:length(whichPts) 
    axesi(i) = nexttile(tile_order(i));

    rate = out(i).rate./out(i).run_dur;
    chLabels = out(i).unchanged_labels;
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
        %pre = cblock-surround:cblock - 1;
        %post = cblock + 1: cblock+surround;
        [pre,post] = get_surround_times(rate,cblock,surround);
    end
    
    ekg = identify_ekg_scalp(out(i).unchanged_labels);
    
    % Calculate change in rate
    switch which_resp
        case 'abs_rate'
            resp = nanmean(rate(:,post),2) - nanmean(rate(:,pre),2);
            rtext = 'Spike rate change (spikes/min)';     
            spikey_idx = ~ekg;
        case 'rel_rate'
            %error('Are you sure you want this, seems like a bad analysis');
            resp = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre),2);
            spikey_idx = nanmean(rate(:,[pre,post]),2) > min_rate & ~ekg;
            rtext = 'Spike rate change (relative)';
        case 'ns'
            resp = nanmean(ns(:,post),2) - nanmean(ns(:,pre),2);
            rtext = 'Node strength change';
            spikey_idx = ~ekg;
        case 'ns_rel'
            resp = (nanmean(ns(:,post),2) - nanmean(ns(:,pre),2))./nanmean(ns(:,pre),2);
            rtext = 'Relative node strength change';
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
    
    % Remove non-spikey
    resp(~spikey_idx) = [];
    predictor(~spikey_idx) = [];
    chLabels(~spikey_idx) = [];
    rate(~spikey_idx,:) = [];
    ns(~spikey_idx,:) = [];
    
    if 0
        figure
        turn_nans_white(rate)
        hold on
        plot([cblock cblock],ylim,'r--','linewidth',3)
        yticks(1:length(chLabels))
        yticklabels(chLabels);
        
        [~,I] = sort(predictor);
        table(chLabels(I),predictor(I),resp(I))
    end
    
    % Monte carlo test to see how big correlation is compared to that at
    % different time points
    if isempty(surround)
        rho = corr(resp,predictor,'Type',corr_type,'rows','pairwise');
        pval = nan;
    else
        [rho,pval] = mc_corr(rate,ns,predictor,...
            cblock,surround,nb,which_resp,corr_type,only_pre);
    end
    
    % Fisher's R to z transformation on the original rhos
    n = sum(~isnan(resp) & ~isnan(predictor));
    [z,z_score,pval2] = fisher_transform(rho,n);
    [~,alt_pval2] = corr(resp,predictor,'Type',corr_type,'rows','pairwise');
    
    % this is to make sure I am getting a reasonable z-score from my fisher
    % transformation.
    if abs(pval2-alt_pval2) > 0.05
        error('oh nos');
    end
    
    % Fill up all z's with info
    all_zs(i,:) = [z,z_score,rho,pval,n,pval2,alt_pval2];
    
    % Plot
    if show_labels
        plot(predictor,resp,'o','color',[1 1 1]);
        text(predictor,resp,chLabels,'horizontalalignment','center')
    else
        plot(predictor,resp,'o','linewidth',2);
    end
    % the rho is the simple rho, but the p-value is the p-value from the MC
    % test
    %title(sprintf('rho = %1.2f, MC p = %1.3f',rho,pval))   
    
    set(gca,'fontsize',15)
    %if i == 8
        xlabel(ptext)
    %end
    
    %if i == 1 || i == 6
        ylabel(rtext)
    %end

end

%% Stouffer's Z-score
% This tests whether the correlation goes in the same direction across
% patients and is significant in aggregate (regardless of whether it is
% larger than at other times)
z = nansum(all_zs(:,2))/sqrt(sum(~isnan(all_zs(:,2))));

% get two sided p-value from the combined z score
pval = 2*normcdf(-abs(z));

% get r back by averaging the z's and z-to-r transforming. Do weighted
% average by sample size
r = tanh(nansum(all_zs(:,1).*all_zs(:,5))./nansum(all_zs(:,5)));
fprintf('\nCombined r = %1.2f, p = %1.3f\n',r,pval);

%% Fisher's test to combine pvals
% This tests whether, in aggregate, the revision time has a
% disproportionately higher correlation relative to other times (regardless
% of whether the direction is different across patients)
X_2 = -2 * nansum(log(all_zs(:,4)));
sum_p = 1-chi2cdf(X_2,2*sum(~isnan(all_zs(:,4))));
fprintf('\nFisher combined p = %1.3f\n',sum_p);


%% Plot aggregate statistics
nexttile(tile_order(end),[2 1])
plot(all_zs(:,3),'o','markersize',20,'linewidth',2)
hold on
ylim([-1 1])
plot(xlim,[0 0],'k--')
yl = ylim;
xl = xlim;
text(xl(1)+0.8*(xl(2)-xl(1)),yl(1)+0.9*(yl(2)-yl(1)),...
    sprintf('Combined r = %1.2f\np = %1.3f',r,pval),'fontsize',15)
xticklabels([])
xlabel('Patient')
ylabel('Correlation coefficient')
set(gca,'fontsize',15)


%% Plot individual statistics
for i = 1:length(whichPts)
    axes(axesi(i))
    xl = xlim;
    yl = ylim;
    xpos = xl(1)+0.5*(xl(2)-xl(1));
    ypos = yl(1) + 0.9*(yl(2)-yl(1));
    text(xpos,ypos,sprintf('r = %1.2f\nMC p = %1.3f',all_zs(i,3),all_zs(i,4)),'fontsize',15)
end

%ylabel(ax,rtext,'fontsize',20)
%xlabel(ax,ptext,'fontsize',20)
print(gcf,[main_spike_results,which_resp,'_',which_pred],'-dpng')


end


function [true_rho,pval] = mc_corr(rate,ns,predictor,cblock,surround,nb,which_resp,corr_type,only_pre)

nblocks = size(rate,2);

%% Get true corr
% Identify pre and post times
%pre = cblock-surround:cblock - 1;
%post = cblock + 1: cblock+surround;
[pre,post] = get_surround_times(rate,cblock,surround);

% Calculate change in rate
switch which_resp
    case 'abs_rate'
        resp = nanmean(rate(:,post),2) - nanmean(rate(:,pre),2);
    case 'rel_rate'
        resp = (nanmean(rate(:,post),2) - nanmean(rate(:,pre),2))./nanmean(rate(:,pre),2);
    case 'ns'
        resp = nanmean(ns(:,post),2) - nanmean(ns(:,pre),2);
    case 'ns_rel'
        resp = (nanmean(ns(:,post),2) - nanmean(ns(:,pre),2))./nanmean(ns(:,pre),2);
end

true_rho = corr(resp,predictor,'Type',corr_type,'rows','pairwise');

mc_rho = nan(nb,1);
for ib = 1:nb
    
    % Make a fake change time
    if only_pre
        fchange = randi([surround+1,cblock-1]);
    else
        fchange = randi([surround+1,nblocks-surround]);
    end
    
    % recalculate change around this time
    %fpre = fchange-surround:fchange-1;
    %fpost = fchange+1:fchange+surround;
    [fpre,fpost] = get_surround_times(rate,fchange,surround);
    
    if length(fpre) == 1 || length(fpost) == 1
        mc_rho(ib) = nan;
        continue
    end
    
    switch which_resp
        case 'abs_rate'
            resp = nanmean(rate(:,fpost),2) - nanmean(rate(:,fpre),2);
        case 'rel_rate'
            resp = (nanmean(rate(:,fpost),2) - nanmean(rate(:,fpre),2))./nanmean(rate(:,fpre),2);
        case 'ns'
            resp = nanmean(ns(:,fpost),2) - nanmean(ns(:,fpre),2);
        case 'ns_rel'
            resp = (nanmean(ns(:,fpost),2) - nanmean(ns(:,fpre),2))./nanmean(ns(:,fpre),2);
    end
    
    mc_rho(ib) = corr(resp,predictor,'Type',corr_type,'rows','pairwise');
    
end


%% Determine number with as or more extreme a correlation
num_as_sig = sum(abs(mc_rho) >= abs(true_rho)); % those more positive or more negative

pval = (num_as_sig+1)/(nb+1);

if isnan(true_rho)
    pval = nan;
end

if 0
    figure
    plot(sort(mc_rho))
    hold on
    plot(xlim,[true_rho true_rho])
    title(sprintf('%1.3f',pval))
    pause
    close(gcf)
end


end