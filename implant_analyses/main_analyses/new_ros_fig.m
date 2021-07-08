function new_ros_fig(whichPts,saved_out)

%{
This is used to perform the analysis associated with Figure 3 of the
implant effect paper. It calculates spike and node strength stability from
pre-to-post revision and determines if the stability is lower than that
expected from random pseudo-revision times.
%}

%% User change parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
%all_surrounds = 12*[2];
main_surround = 3; %24 hour peri-revision surround
nb = 1e4; % number of monte carlo iterations 
ex_p = 5;
do_norm = 0; % doesn't seem to make much difference so I will keep 0 for simplicity (a node strength thign)
do_cat = 0;

%% Other info
n_surrounds = length(all_surrounds);
all_metrics = {'rate','ns'};
n_metrics = length(all_metrics);

% divide by 2 to get time in hours but then multiply by 2 to get the full
% pre+post surround
true_times = all_surrounds; 

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

%% Get names
names = cell(length(whichPts),1);
for i = 1:length(whichPts)
    names{i} = out(i).name;
end

%% Main rate correlation analysis
npts = length(whichPts);
all_all_true_r = nan(n_surrounds,length(whichPts),n_metrics);
all_all_mc_r = nan(n_surrounds,length(whichPts),nb,n_metrics);
all_all_p = nan(n_surrounds,n_metrics);
all_all_all_p = nan(n_surrounds,length(whichPts),n_metrics);

% Loop over metrics
for im = 1:n_metrics
    allcatlist = [];
    metric = all_metrics{im};
    
    % Loop over surround times
    for is = 1:n_surrounds

        surround = all_surrounds(is);

        % Initialize correlation stuff
        all_ps = nan(length(whichPts),1);

        for i = 1:length(whichPts)
            
            switch metric
                case 'rate'
                    rate = out(i).rate;
                case 'ns'
                    if do_norm
                        rate = out(i).metrics.ns_norm;
                    else
                        rate = out(i).metrics.ns;
                    end
            end
            
            
            cblock = out(i).change_block;

            % Remove EKG and scalp electrodes
            ekg = identify_ekg_scalp(out(i).unchanged_labels);
            rate(ekg,:) = [];

            % Get rho and statistics (in the math folder)
            [pval_curr,true_rho,mc_rho] = ...
                compare_rhos(rate,cblock,surround,nb,'Spearman',[],0,out(i).rate);

            all_all_mc_r(is,i,:,im) = mc_rho;
            all_all_true_r(is,i,im) = true_rho;
            all_all_all_p(is,i,im) = pval_curr;
            all_ps(i) = pval_curr;

             %% CAT analysis (I don't end up using this)
            if do_cat
            if is == main_surround
                % pre and post
                [pre,post] = get_surround_times(rate,cblock,surround);
                pre_mean = nanmean(rate(:,pre),2);
                post_mean = nanmean(rate(:,post),2);
                catlist = cat_plot(pre_mean,post_mean);
                length_diff = length(catlist) - size(allcatlist,1);

                if length_diff > 0 % new one longer
                    % pad old with nans
                    allcatlist = [allcatlist;nan(length_diff,size(allcatlist,2))];
                elseif length_diff < 0 % old one older
                    % pad new with nans
                    catlist = [catlist;nan(-length_diff,1)];
                end

         
                % add to the full patient list
                allcatlist = [allcatlist,catlist];
            end
            end

        end

        % Fisher test to combine pvalues
        if sum(isnan(all_ps)) > 0, error('why'); end
        X_2 = -2 * sum(log(all_ps));
        sum_p = 1-chi2cdf(X_2,2*length(all_ps));

        all_all_p(is,im) = sum_p;
        
       
        
    end
    
    if do_cat
        if im == 1
            all_cat_list = nan(size(allcatlist,1),size(allcatlist,2),2);
        end
        all_cat_list(:,:,im) = allcatlist;
    end
    
end


figure
set(gcf,'position',[100 87 733 710])
tiledlayout(4,2,'TileSpacing','tight','Padding','tight')
tileorder = [1,3,5,7,2,4,6,8];
count = 0;


%% plot stuff
% Loop over spike vs ns
for im = 1:n_metrics
    
    metric = all_metrics{im};
    
    if 1
    %% Raster example patient
    count = count+1;
    nexttile(tileorder(count))
    
    switch metric
        case 'rate'
            rate = out(ex_p).rate./out(ex_p).run_dur;
            rtext = 'Spikes';
        case 'ns'
            if do_norm
                rate = out(ex_p).metrics.ns_norm;
            else
                rate = out(ex_p).metrics.ns;
            end
            rtext = 'Node strength';
    end
    ekg = identify_ekg_scalp(out(ex_p).unchanged_labels);
    rate(ekg,:) = [];
    curr_times = (1:size(rate,2)) * out(ex_p).block_dur;
    curr_change = out(ex_p).change_block*out(ex_p).block_dur;
    h = turn_nans_white(rate);
    set(h,'XData',[0:curr_times(end)]);
    xlim([0 curr_times(end)])
    hold on
    cp = plot([curr_change curr_change],ylim,'r--','linewidth',3);
    %xticks([0 100 200 300])
    %xticklabels({'0 h','100 h','200 h','300 h'})
    xlabel('Hour')
    yticklabels([])
    if im == 1
        ylabel('Electrode')
    end
    set(gca,'fontsize',15)
    %c = colorbar('location','westoutside');
    %ylabel(c,rtext,'fontsize',15)
    title(rtext,'fontweight','normal')
    end
    
    %% Histogram
    count = count+1;
    nexttile(tileorder(count))
    all_rate_diff = [];
    surround = all_surrounds(1);
    for i = 1:npts
        switch metric
            case 'rate'
                rate = out(i).rate./out(i).run_dur;
                rtext = 'Rate change (spikes/min)';
            case 'ns'
                if do_norm
                    rate = out(i).metrics.ns_norm;
                else
                    rate = out(i).metrics.ns;
                end
                rtext = 'NS change';
        end
        cblock = out(i).change_block;
        ekg = identify_ekg_scalp(out(i).unchanged_labels);
        rate(ekg,:) = [];
        [pre,post] = get_surround_times(rate,cblock,surround);
        rate_pre = nanmean(rate(:,pre),2);
        rate_post = nanmean(rate(:,post),2);
        rate_diff = rate_post-rate_pre;
        all_rate_diff = [all_rate_diff;rate_diff];
    end
    histogram(all_rate_diff);
    %{
    if im == 1
        legend('# electrodes','fontsize',15,'location','northwest')
    else
        legend('# electrodes','fontsize',15,'location','northwest')
    end
    %}
    if im == 1
        ylabel('# Electrodes')
    end
    xlabel(rtext)
    set(gca,'fontsize',15)
    
    %% First surround corr
    count = count+1;
    nexttile(tileorder(count))
    is = main_surround;
    true_r = squeeze(all_all_true_r(is,:,im));
    mc_r = squeeze(all_all_mc_r(is,:,:,im));
    p = all_all_p(is,im);

    plot(1+0.05*rand(npts,1),true_r,'o','markersize',15,'linewidth',2)
    hold on
    plot(2+0.05*rand(npts,1),mean(mc_r,2),'o','markersize',15,'linewidth',2)
    yl = ylim;
    ylim([yl(1) yl(1)+1.3*(yl(2)-yl(1))])
    yl = ylim;
    plot([1 2],[yl(1)+0.83*(yl(2)-yl(1)) yl(1)+0.83*(yl(2)-yl(1))],'k','linewidth',2)
    text(1.5,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(p,1),...
        'horizontalalignment','center','fontsize',15)
    xlim([0.5 2.5])
    xticks([1 2])
    xticklabels({'True','Monte Carlo'})
    if im == 1
        ylabel({'Spike stability'})
    else
        ylabel({'NS stability'})
    end
    set(gca,'fontsize',15)

    if im == 1
        fprintf(['\n However, on a group level, the spike stability '...
            ' (mean %1.2f, std %1.2f) was significantly less than expected '...
            'by chance (Fisher’s method, p = %1.3f).\n'],...
            mean(true_r),std(true_r),p);
    elseif im == 2
        fprintf(['\n However, the group node strength stability '...
            ' (mean %1.2f, std %1.2f) was significantly less than expected '...
            'by chance (Fisher’s method, p = %1.3f).\n'],...
            mean(true_r),std(true_r),p);
    end

    %% All surrounds
    count = count+1;
    nexttile(tileorder(count))
    th = errorbar((1:n_surrounds)'-0.2,mean(all_all_true_r(:,:,im),2),std(all_all_true_r(:,:,im),[],2),...
        'o','linewidth',2,'markersize',10);
    hold on
    mch = errorbar((1:n_surrounds)'+0.2,mean(all_all_mc_r(:,:,:,im),[2 3]),std(all_all_mc_r(:,:,:,im),[],[2 3]),...
        'o','linewidth',2,'markersize',10);
    yl = ylim;
    ylim([yl(1) yl(1)+1.3*(yl(2)-yl(1))])
    yl = ylim;
    for is = 1:n_surrounds
        plot([is-0.2,is+0.2],[yl(1)+0.8*(yl(2)-yl(1)) yl(1)+0.8*(yl(2)-yl(1))],...
            'k','linewidth',2)
        text(is,yl(1)+0.9*(yl(2)-yl(1)),get_asterisks(all_all_p(is,im),1),...
        'horizontalalignment','center','fontsize',20)
    end
    
    xlabel('Peri-revision period (hours)')
    xticks(1:n_surrounds)
    xticklabels(true_times)
    xlim([0 n_surrounds+1])
    if im == 1
        ylabel({'Spike stability'})
        legend([th,mch],{'True','Monte Carlo'},'location','southeast','fontsize',15)
    else
        ylabel({'NS stability'})
    end
    set(gca,'fontsize',15)
    
    %% CAT analysis
    if do_cat
        count = count + 1;
        nexttile(tileorder(count));
        
        % Average over patients
        mean_cat = squeeze(nanmean(all_cat_list(:,:,im),2));
        std_cat = squeeze(nanstd(all_cat_list(:,:,im),[],2));
        
        % Plot
        shaded_error_bars(1:length(mean_cat),mean_cat,std_cat,[])
        xlabel('Number of electrodes')
        if im == 1
            ylabel('Spike stability')
        else
            ylabel('NS stability')
        end
        xlim([1 100])
        set(gca,'fontsize',15);
        
    end
end

%% Annotations
annotation('textbox',[0 0.90 0.1 0.1],'String','A','fontsize',20,'linestyle','none')
annotation('textbox',[0.5 0.90 0.1 0.1],'String','B','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.67 0.1 0.1],'String','C','fontsize',20,'linestyle','none')
annotation('textbox',[0.5 0.67 0.1 0.1],'String','D','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.42 0.1 0.1],'String','E','fontsize',20,'linestyle','none')
annotation('textbox',[0.5 0.42 0.1 0.1],'String','F','fontsize',20,'linestyle','none')
annotation('textbox',[0 0.18 0.1 0.1],'String','G','fontsize',20,'linestyle','none')
annotation('textbox',[0.5 0.18 0.1 0.1],'String','H','fontsize',20,'linestyle','none')

print(gcf,[main_spike_results,'new_ros'],'-dpng')

%% Save table of p-values
names = ['all';names];
all_spike = [all_all_p(:,1)';all_all_all_p(:,:,1)'];
all_ns = [all_all_p(:,2)';all_all_all_p(:,:,2)'];


spike_T = cell2table(arrayfun(@(x) sprintf('MC p = %1.3f',x),all_spike,...
    'UniformOutput',false),...
    'VariableNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'RowNames',names);
ns_T = cell2table(arrayfun(@(x) sprintf('MC p = %1.3f',x),all_ns,...
    'UniformOutput',false),...
    'VariableNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'RowNames',names);
%{
writetable(spike_T,[main_spike_results,'spike_ros.csv'],'WriteRowNames',true)  
writetable(ns_T,[main_spike_results,'ns_ros.csv'],'WriteRowNames',true)  
    %}

%% Also save these to a specific range in the Supplemental Table 1
agg_spike_T = cell2table(arrayfun(@(x) sprintf('MC %s',pretty_p_text(x)),all_all_p(:,1),...
    'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));
agg_ns_T = cell2table(arrayfun(@(x) sprintf('MC %s',pretty_p_text(x)),all_all_p(:,2),...
    'UniformOutput',false),...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false));
writetable(agg_spike_T,[main_spike_results,'spike_ros.csv'],'WriteRowNames',true)  
writetable(agg_ns_T,[main_spike_results,'ns_ros.csv'],'WriteRowNames',true)  

writetable(agg_spike_T,[main_spike_results,'Supplemental Table 1.xlsx'],'Range','C2:C12','WriteVariableNames',false)
writetable(agg_ns_T,[main_spike_results,'Supplemental Table 1.xlsx'],'Range','D2:D12','WriteVariableNames',false)

end